"""Original ROOT identities, bin layouts and sealed job-plan utilities."""
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import tempfile
from functools import lru_cache
from ruamel.yaml import YAML

PJ = Path(__file__).resolve().parents[1]
PRODUCTS = {"data": "data_files_list", "sim_signal": "sim_signal_files_list", "sim_inclusive": "sim_inclusive_files_list"}
SCHEMA = "PhotonJetTrees_v1"
BACKEND = "photonjet-cpp-v1"
PLAN_FORMAT_VERSION = 4
PROVENANCE_FORMAT_VERSION = 4
COMPLETION_FORMAT_VERSION = 3
INPUT_VERIFICATION = "worker-md5-stat-v1"
CODE_SOURCES = ("histmakers/PhotonJetHistMaker.C", "histmakers/PhotonJetReader.h",
                "histmakers/PhotonJetPhysics.h", "histmakers/PhotonJetHistograms.h",
                "histmakers/PhotonJetConfig.h", "histmakers/PhotonJetIO.h", "support/CrossSectionWeights.h")


def md5(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "md5").hexdigest()


def config(path):
    return YAML(typ="safe").load(Path(path).read_text())


def label(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", value) or value in (".", ".."):
        raise ValueError(f"invalid path label: {value!r}")
    return value


def product_name(cfg, product):
    system = cfg["photonjet"]["system"]
    if system not in ("pp", "auau") or product not in PRODUCTS:
        raise ValueError("unknown product/system")
    if product == "data":
        return "data/" + system
    return "simulation/" + ("pp_" if system == "pp" else "auau_embedded_") + ("photonjet" if product == "sim_signal" else "inclusive")


@lru_cache(maxsize=32)
def _master(path):
    """Freeze each master once per command; publication checks its digest again."""
    path = Path(path)
    parts = []
    for raw in path.read_text().splitlines():
        fields = raw.split("#", 1)[0].split()
        if not fields:
            continue
        if len(fields) != 1 or not Path(fields[0]).is_absolute():
            raise ValueError(f"master list needs one absolute path per row: {path}")
        parts.append(fields[0])
    if not parts or len(set(parts)) != len(parts):
        raise ValueError(f"empty or duplicate master input list: {path}")
    return tuple(parts), md5(path)


def master_parts(cfg, product):
    path = Path(cfg["photonjet"][PRODUCTS[product]]).resolve()
    return path, _master(str(path))[0]


def master_digest(cfg, product):
    path, _ = master_parts(cfg, product)
    return _master(str(path))[1]


def indices(text, size):
    """Explicit zero-based indices and half-open ranges, e.g. 0,4:7."""
    if text is None:
        return list(range(size))
    result = []
    for field in text.split(","):
        if ":" in field:
            first, last = map(int, field.split(":"))
            if first >= last:
                raise ValueError("empty/reversed part range")
            result.extend(range(first, last))
        else:
            result.append(int(field))
    if not result or len(set(result)) != len(result) or min(result) < 0 or max(result) >= size:
        raise ValueError("empty, duplicate or out-of-range original part selection")
    return result


def release_label(cfg, explicit=None):
    return label(explicit or os.environ.get("PJ_RELEASE_ID") or Path(cfg["photonjet"]["release"]).name)


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", dir=path.parent, prefix=path.name + ".tmp.", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(text)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def save(path, value):
    atomic_text(path, json.dumps(value, indent=2, allow_nan=False) + "\n")


def code_digests():
    return {name: md5(PJ / name) for name in CODE_SOURCES}


def bin_layout(cfg):
    """Mirror Config's finite ordered edges and versioned output layout."""
    import math
    system = cfg["photonjet"]["system"]
    if system not in ("pp", "auau"):
        raise ValueError("unknown collision system")
    eta = cfg["analysis"]["eta_bins"]
    centrality = cfg["photonjet"].get("centrality_bins", [])
    def edges(values, name):
        if len(values) < 2 or any(not math.isfinite(v) for v in values) or any(a >= b for a, b in zip(values, values[1:])):
            raise ValueError(f"{name} needs finite strictly increasing edges")
    edges(eta, "eta_bins")
    if system == "auau":
        edges(centrality, "centrality_bins")
        if centrality[0] < 0 or centrality[-1] > 100:
            raise ValueError("centrality edges must lie in [0,100]")
    elif centrality:
        raise ValueError("pp uses one inclusive centrality slot, without centrality edges")
    return dict(format_version=1, centrality_edges=centrality, eta_edges=eta,
                centrality_policy="inclusive-pp" if system == "pp" else "low-inclusive-high-exclusive",
                eta_policy="outer-exclusive-interior-to-higher",
                n_centrality=len(centrality)-1 if centrality else 1, n_eta=len(eta)-1)


def require_single_cell(layout, consumer):
    if layout["n_centrality"] != 1 or layout["n_eta"] != 1:
        raise ValueError(f"{consumer} does not support multiple centrality/eta cells; explicit cell selection and normalization are required")


def plan_source(path, index, cfg, product, release):
    """Record cheap source metadata; workers seal ROOT identity and contents."""
    path = Path(path)
    if not path.is_absolute() or index < 0:
        raise ValueError("original parts need an absolute path and nonnegative index")
    identity = path.resolve(strict=True)
    state = identity.stat()
    if not stat.S_ISREG(state.st_mode):
        raise ValueError(f"original input is not a regular file: {path}")
    return dict(original_part_index=index,
                original_part_id=f"{release}/{product_name(cfg, product)}/part/{index}",
                source_path=str(path), source_bytes=state.st_size,
                source_mtime_ns=state.st_mtime_ns)


def inspect_source(path, index, cfg, product, release):
    """Inspect ROOT headers using the same directory identity as TFile::GetUUID."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    record = plan_source(path, index, cfg, product, release)
    path = Path(record["source_path"])
    digest = md5(path)
    source = ROOT.TFile.Open(str(path))
    if not source or source.IsZombie() or source.TestBit(ROOT.TFile.kRecovered):
        if source:
            source.Close()
        raise ValueError(f"cannot open intact original ROOT file: {path}")
    try:
        def count(name, required=False):
            tree = source.Get(name)
            if not tree:
                if required:
                    raise ValueError(f"missing original tree {name}: {path}")
                return 0
            if not tree.InheritsFrom("TTree"):
                raise ValueError(f"wrong-type original tree {name}: {path}")
            return int(tree.GetEntries())
        record.update(source_uuid=str(source.GetUUID().AsString()), source_md5=digest,
                      events=count("events", True), photons=count("photons", True),
                      truth_photons=count("truthPhotons", product != "data"),
                      source_link_rows=count("recoTruthLinks"), truth_jet_rows=count("truthJets"),
                      links_available=bool(source.Get("recoTruthLinks")), truth_jets_available=bool(source.Get("truthJets")))
    finally:
        source.Close()
    state = path.stat()
    if state.st_size != record["source_bytes"] or state.st_mtime_ns != record["source_mtime_ns"]:
        raise ValueError(f"original file changed during inspection: {path}")
    return record


def dependencies(cfg, product):
    """Scientific file dependencies; inline response-prior settings live in the config."""
    settings = cfg["photonjet"]
    result = {}
    if product == "data":
        paths = []
    elif settings.get("weight_mode", "sample_map" if settings["system"] == "auau" else "stored") != "stored":
        sample_key = "sample_map_signal" if product == "sim_signal" else "sample_map_inclusive"
        paths = [settings.get(sample_key, settings.get("sample_map", "")), settings.get("vertex_weight_file", "")]
    else:
        paths = []
    paths.append(cfg["analysis"].get("run_list_file", ""))
    if product == "sim_signal":
        paths.append(settings.get("external_mbd_eff_file", ""))
    for path in filter(None, paths):
        result[str(Path(path).resolve())] = md5(path)
    return result
