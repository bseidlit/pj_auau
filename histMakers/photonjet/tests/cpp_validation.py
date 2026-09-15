"""Backend-independent ROOT snapshots and full numerical-contract comparisons."""
from array import array
import math

CRITERIA = ("Exact inventories/types, axes/labels, integer counters and fill entries; "
            "weighted values, errors, Sumw2 and moments within 1e-12 relative/absolute; "
            "all flow bins, weighted-efficiency settings and RooUnfold components.")


def compare(reference, actual, path="", differences=None):
    differences = [] if differences is None else differences
    if isinstance(reference, dict):
        if not isinstance(actual, dict) or reference.keys() != actual.keys():
            differences.append((path, "different keys", list(reference), list(actual)))
        else:
            for key in reference:
                compare(reference[key], actual[key], path + "/" + key, differences)
    elif isinstance(reference, list):
        if not isinstance(actual, list) or len(reference) != len(actual):
            differences.append((path, "different lengths", len(reference), len(actual)))
        else:
            for i, (a, b) in enumerate(zip(reference, actual)):
                compare(a, b, path + "/" + str(i), differences)
    elif isinstance(reference, (float, int)) and not isinstance(reference, bool):
        exact = "/axes/" in path or path.endswith(("/entries", "/sumw2_size", "/statistic"))
        if any(name in path for name in ("h_pj_cutflow/", "h_pj_flagcheck/", "h_pj_input_parts/", "h_pj_auau_threshold_recipe/")):
            exact = True
        if not math.isfinite(actual) or (reference != actual if exact else abs(reference-actual) > 1e-12*max(1., abs(reference), abs(actual))):
            differences.append((path, reference, actual))
    elif reference != actual:
        differences.append((path, reference, actual))
    return differences


def histogram_snapshot(hist, titles=False):
    axes = []
    for get_axis in ("GetXaxis", "GetYaxis", "GetZaxis")[:hist.GetDimension()]:
        axis = getattr(hist, get_axis)()
        axes.append(dict(edges=[axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)],
                         labels=[str(axis.GetBinLabel(i)) for i in range(1, axis.GetNbins() + 1)]))
    stats = array("d", [0.] * 13)
    hist.GetStats(stats)
    result = dict(type=hist.ClassName(), axes=axes, entries=hist.GetEntries(), stats=list(stats),
                sumw2_size=hist.GetSumw2N(), sumw2=[hist.GetSumw2().At(i) for i in range(hist.GetSumw2N())],
                content=[hist.GetBinContent(i) for i in range(hist.GetNcells())],
                error=[hist.GetBinError(i) for i in range(hist.GetNcells())])

    if titles:
        result["title"] = str(hist.GetTitle())
        for axis, get_axis in zip(axes, ("GetXaxis", "GetYaxis", "GetZaxis")):
            axis["title"] = str(getattr(hist, get_axis)().GetTitle())
    return result


def root_snapshot(path, extra_metadata=(), titles=False):
    import ROOT
    ROOT.gSystem.Load("libRooUnfold")
    ROOT.gInterpreter.Declare("#include <RooUnfoldResponse.h>")
    file = ROOT.TFile.Open(str(path))
    assert file and not file.IsZombie()
    objects = {}
    metadata = {}
    for key in file.GetListOfKeys():
        name = key.GetName()
        obj = file.Get(name)
        if obj.InheritsFrom("TH1"):
            objects[name] = histogram_snapshot(obj, titles)
        elif obj.InheritsFrom("TEfficiency"):
            objects[name] = dict(type=obj.ClassName(), passed=histogram_snapshot(obj.GetPassedHistogram(), titles),
                total=histogram_snapshot(obj.GetTotalHistogram(), titles), weighted=bool(obj.UsesWeights()),
                statistic=int(obj.GetStatisticOption()), confidence=obj.GetConfidenceLevel(),
                beta_alpha=obj.GetBetaAlpha(), beta_beta=obj.GetBetaBeta(), weight=obj.GetWeight())
        elif obj.InheritsFrom("RooUnfoldResponse"):
            objects[name] = dict(type=obj.ClassName(), **{method: histogram_snapshot(getattr(obj, method)(), titles)
                for method in ("Hresponse", "Htruth", "Hmeasured", "Hfakes")})
        else:
            if name in extra_metadata and obj.InheritsFrom("TNamed"):
                metadata[name] = str(obj.GetTitle())
            else:
                assert obj.InheritsFrom("TObjString"), (name, obj.ClassName())
                metadata[name] = str(obj.GetString())
        if titles and name in objects:
            objects[name]["title"] = str(obj.GetTitle())
    file.Close()
    assert set(metadata) == {"config", "input_manifest", "photonjet_provenance", *extra_metadata}
    return dict(objects=objects, metadata=metadata)
