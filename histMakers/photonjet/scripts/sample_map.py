#!/usr/bin/env python3
"""Part-index -> PYTHIA pT-hat sample map for the PhotonJetTrees simulation products.

The four simulation products are concatenations of contiguous blocks of PYTHIA
pT-hat samples and the trees carry no sample label:
  pp_photonjet             photon10 -> photon20 -> photon5
  pp_inclusive             jet12 -> jet20 -> jet30 -> jet40 -> jet8
  auau_embedded_photonjet  photon10 -> photon20
  auau_embedded_inclusive  jet12 -> jet20 -> jet30 -> jet40

Classifier (one cheap read per part):
  photon products: ``truthPhotons`` (prompt_class == 1), per-event leading
      truth_photon_pt, 5th percentile p5 (photon5 ~4.4, photon10 ~9-11,
      photon20 ~18.5 GeV).  Sample = nearest of {5,10,20} to p5 + 1.0.
  jet products:    ``truthJets`` (first --jet-entry-stop rows), per-event
      leading truth_jet_pt, 5th percentile p5 (jet8 ~7.5, jet12 ~11.7,
      jet20 ~19.9, jet30 ~30, jet40 ~40 GeV).  Sample = nearest of
      {8,12,20,30,40} to p5 + 0.5.

Scan strategy (blocks are contiguous, so all parts need not be opened): classify
every --stride-th part plus the last one, bisect every label transition down to
the exact boundary, then verify the 2 parts on each side of every boundary.

p+p cross-check: the stored ``event_weight`` (events tree) is 1.0 for most
events and (sigma_sample/sigma_ref) x 0.7088 otherwise; its maximum over
event_weight != 1 is compared with the value expected for the assigned sample
on the first and last part of every range.

Also writes, per range, the PPG12 sample name / truth window / weight from
efficiencytool/CrossSectionWeights.h.

Output: <out-dir>/sample_map_<product>.yaml (one per product).

Example:
  python3 efficiencytool/photonjet/scripts/sample_map.py                     # all four products
  python3 efficiencytool/photonjet/scripts/sample_map.py --products pp_inclusive --stride 20 -v
"""
import argparse
import datetime as _dt
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import uproot

DEFAULT_RELEASE = "/sphenix/tg/tg01/bulk/jbennett/PhotonJetTrees_v1"
PHOTONJET_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(os.path.dirname(PHOTONJET_DIR))
HEADER_PATH = os.path.join(REPO_ROOT, "efficiencytool", "CrossSectionWeights.h")

PRODUCTS = {
    "pp_photonjet": {"channel": "photon", "expected_order": ["photon10", "photon20", "photon5"], "weight_check": True},
    "pp_inclusive": {"channel": "jet", "expected_order": ["jet12", "jet20", "jet30", "jet40", "jet8"], "weight_check": True},
    "auau_embedded_photonjet": {"channel": "photon", "expected_order": ["photon10", "photon20"], "weight_check": False},
    "auau_embedded_inclusive": {"channel": "jet", "expected_order": ["jet12", "jet20", "jet30", "jet40"], "weight_check": False},
}
THRESHOLDS = {"photon": [5, 10, 20], "jet": [8, 12, 20, 30, 40]}
# leading truth jet p5 sits ~0.5 GeV below the pT-hat threshold, leading prompt photon p5 ~1 GeV below
P5_OFFSET = {"photon": 1.0, "jet": 0.5}
LOW_STAT_EVENTS = 200

# Tree convention for the p+p event_weight maximum (over event_weight != 1):
# (sigma_sample / sigma_ref) x 0.7088, sigma_ref = photon20 (photons) / jet50 (jets).
TREE_WEIGHT_FACTOR = 0.7088
TREE_EVENT_WEIGHT_MAX = {
    "photon5": 795.3, "photon10": 37.74, "photon20": 0.7088,
    "jet8": 1.115e6, "jet12": 1.4449e5, "jet20": 6071.0, "jet30": 245.3, "jet40": 13.14,
}
WEIGHT_TOL = 0.01  # relative tolerance for the event_weight cross-check


# ----------------------------------------------------------------------------
# YAML writer (ruamel -> PyYAML -> plain)
# ----------------------------------------------------------------------------
def dump_yaml(obj, path):
    try:
        from ruamel.yaml import YAML
        y = YAML()
        y.default_flow_style = False
        y.width = 4096
        with open(path, "w") as f:
            y.dump(obj, f)
        return "ruamel.yaml"
    except Exception:
        pass
    try:
        import yaml
        with open(path, "w") as f:
            yaml.safe_dump(obj, f, sort_keys=False, default_flow_style=False, width=4096)
        return "pyyaml"
    except Exception:
        pass
    import json
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)  # JSON is valid YAML
    return "json"


# ----------------------------------------------------------------------------
# CrossSectionWeights.h parsing
# ----------------------------------------------------------------------------
def parse_header(path=HEADER_PATH):
    """Return (constants{name: value}, samples{name: {...window/weight...}}, notes)."""
    txt = open(path).read()
    consts = {m.group(1): float(m.group(2))
              for m in re.finditer(r"constexpr\s+float\s+(\w+)\s*=\s*([0-9.eE+-]+)f?\s*;", txt)}
    samples = {}
    # blocks:  if (filetype == "photon5") { ... }   (first name of an || chain is taken)
    for m in re.finditer(r'filetype\s*==\s*"(\w+)"[^{]*\{(.*?)\}', txt, flags=re.S):
        name, body = m.group(1), m.group(2)
        if name in samples:
            continue
        d = {}
        for key in ("photon_pt_lower", "photon_pt_upper", "jet_pt_lower", "jet_pt_upper", "cluster_ET_upper"):
            mm = re.search(rf"c\.{key}\s*=\s*([0-9.]+)", body)
            if mm:
                d[key] = float(mm.group(1))
        mw = re.search(r"c\.weight\s*=\s*(\w+)\s*/\s*(\w+)", body)
        if mw:
            num, den = mw.group(1), mw.group(2)
            d["weight_expr"] = f"{num}/{den}"
            d["weight_num_name"], d["weight_den_name"] = num, den
            if num in consts and den in consts:
                d["weight_value"] = consts[num] / consts[den]
        elif re.search(r"c\.weight\s*=\s*1\.0f", body):
            d["weight_expr"] = "1.0"
            d["weight_value"] = 1.0
        d["isbackground"] = "isbackground = true" in body
        samples[name] = d
    return consts, samples


def header_vs_git_head(path=HEADER_PATH):
    """Constants that differ between the working-tree header and the committed (HEAD) version."""
    import subprocess
    rel = os.path.relpath(path, REPO_ROOT)
    try:
        head_txt = subprocess.run(["git", "-C", REPO_ROOT, "show", f"HEAD:{rel}"],
                                  capture_output=True, text=True, check=True).stdout
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}
    head_consts = {m.group(1): float(m.group(2))
                   for m in re.finditer(r"constexpr\s+float\s+(\w+)\s*=\s*([0-9.eE+-]+)f?\s*;", head_txt)}
    wt_consts = {m.group(1): float(m.group(2))
                 for m in re.finditer(r"constexpr\s+float\s+(\w+)\s*=\s*([0-9.eE+-]+)f?\s*;", open(path).read())}
    diff = {k: {"working_tree": wt_consts.get(k), "HEAD": head_consts.get(k)}
            for k in sorted(set(wt_consts) | set(head_consts)) if wt_consts.get(k) != head_consts.get(k)}
    return {"header": rel, "constants_differing_from_HEAD": diff,
            "note": ("values quoted in this file come from the working-tree header; "
                     "the p+p tree event_weight convention was cross-checked against both")}


def ppg12_info(sample, consts, samples):
    """PPG12 sample name, truth window and weight normalisation from the header."""
    s = samples.get(sample)
    if s is None:
        return {"ppg12_sample": None, "note": f"{sample} not found in CrossSectionWeights.h"}
    info = {"ppg12_sample": sample, "header": os.path.relpath(HEADER_PATH, REPO_ROOT)}
    if s.get("isbackground"):
        info["truth_window"] = {
            "variable": "leading truth jet pT (GeV)",
            "low": s.get("jet_pt_lower"), "high": s.get("jet_pt_upper"),
            "interval": f"[{s.get('jet_pt_lower'):g}, {s.get('jet_pt_upper'):g})",
            "cluster_ET_upper": s.get("cluster_ET_upper"),
        }
        ref_name = s.get("weight_den_name", "jet50cross")
    else:
        info["truth_window"] = {
            "variable": "max truth photon pT (GeV)",
            "low": s.get("photon_pt_lower"), "high": s.get("photon_pt_upper"),
            "interval": f"[{s.get('photon_pt_lower'):g}, {s.get('photon_pt_upper'):g})",
        }
        ref_name = s.get("weight_den_name", "photon20cross")
    xs_name = f"{sample}cross"
    info["weight"] = {
        "expr": s.get("weight_expr"),
        xs_name: consts.get(xs_name),
        ref_name: consts.get(ref_name),
        "value": s.get("weight_value"),
    }
    return info


# ----------------------------------------------------------------------------
# Per-part classifier
# ----------------------------------------------------------------------------
def _leading_per_event(arr, pt):
    """pt maximum per (source_file_index, event_id_hi, event_id_lo) group."""
    if len(pt) == 0:
        return np.zeros(0)
    key = np.stack([arr["source_file_index"].astype(np.int64),
                    arr["event_id_hi"].astype(np.int64),
                    arr["event_id_lo"].astype(np.int64)], axis=1)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    mx = np.full(inv.max() + 1, -np.inf)
    np.maximum.at(mx, inv, pt)
    return mx


def nearest_threshold(value, channel):
    th = THRESHOLDS[channel]
    d = sorted((abs(value - t), t) for t in th)
    best, second = d[0], d[1]
    return best[1], float(second[0] - best[0])


class PartClassifier:
    def __init__(self, release, product, channel, jet_entry_stop=200000, read_weights=True, verbose=False):
        self.product = product
        self.channel = channel
        self.jet_entry_stop = jet_entry_stop
        self.read_weights = read_weights
        self.verbose = verbose
        self.files_txt = os.path.join(release, "simulation", product, "files.txt")
        with open(self.files_txt) as f:
            self.files = [ln.strip() for ln in f if ln.strip()]
        self.cache = {}
        self.n_opened = 0

    def __len__(self):
        return len(self.files)

    def _classify(self, idx):
        path = self.files[idx]
        res = {"part": idx, "file": os.path.basename(path), "channel": self.channel}
        with uproot.open(path) as f:
            ev = f["events"]
            res["n_events"] = int(ev.num_entries)
            if self.read_weights and res["n_events"] > 0:
                w = ev["event_weight"].array(library="np")
                nonunit = w[w != 1.0]
                res["event_weight"] = {
                    "max_raw": float(w.max()),
                    "max_nonunit": float(nonunit.max()) if len(nonunit) else None,
                    "min_nonunit": float(nonunit.min()) if len(nonunit) else None,
                    "frac_unit": float((w == 1.0).mean()),
                    "n_unique": int(len(np.unique(w))),
                }
            if self.channel == "photon":
                t = f["truthPhotons"]
                a = t.arrays(["truth_photon_pt", "prompt_class", "source_file_index", "event_id_hi", "event_id_lo"],
                             library="np")
                sel = a["prompt_class"] == 1
                pt = a["truth_photon_pt"][sel]
                res["n_objects"] = int(len(pt))
                res["n_tree_entries"] = int(t.num_entries)
                sub = {k: v[sel] for k, v in a.items()}
                lead = _leading_per_event(sub, pt)
                if len(pt):
                    res["object_pt_p5"] = float(np.percentile(pt, 5))
                    res["object_pt_p50"] = float(np.percentile(pt, 50))
            else:
                t = f["truthJets"]
                a = t.arrays(["truth_jet_pt", "source_file_index", "event_id_hi", "event_id_lo"],
                             entry_stop=self.jet_entry_stop, library="np")
                pt = a["truth_jet_pt"]
                res["n_objects"] = int(len(pt))
                res["n_tree_entries"] = int(t.num_entries)
                lead = _leading_per_event(a, pt)
        res["n_events_grouped"] = int(len(lead))
        if len(lead) == 0:
            res["label"] = None
            res["reason"] = "empty (no truth objects)"
            return res
        p5 = float(np.percentile(lead, 5))
        p50 = float(np.percentile(lead, 50))
        res["leading_pt_p5"] = p5
        res["leading_pt_p50"] = p50
        lab_val, margin = nearest_threshold(p5 + P5_OFFSET[self.channel], self.channel)
        res["label"] = f"{self.channel}{lab_val}"
        res["margin_gev"] = margin
        res["low_stat"] = bool(len(lead) < LOW_STAT_EVENTS)
        return res

    def get(self, idx):
        if idx not in self.cache:
            r = self._classify(idx)
            self.cache[idx] = r
            self.n_opened += 1
            if self.verbose:
                print(f"    part {idx:6d}: {r.get('label')} p5={r.get('leading_pt_p5', float('nan')):.2f} "
                      f"p50={r.get('leading_pt_p50', float('nan')):.2f} nev={r['n_events']} "
                      f"ngrp={r['n_events_grouped']}", flush=True)
        return self.cache[idx]

    def label(self, idx):
        return self.get(idx)["label"]

    def prefetch(self, indices, workers=8):
        todo = [i for i in indices if i not in self.cache]
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(self._classify, i): i for i in todo}
            for fut in as_completed(futs):
                i = futs[fut]
                self.cache[i] = fut.result()
                self.n_opened += 1
                if self.verbose:
                    r = self.cache[i]
                    print(f"    part {i:6d}: {r.get('label')} p5={r.get('leading_pt_p5', float('nan')):.2f} "
                          f"p50={r.get('leading_pt_p50', float('nan')):.2f} nev={r['n_events']} "
                          f"ngrp={r['n_events_grouped']}", flush=True)


# ----------------------------------------------------------------------------
# Boundary search
# ----------------------------------------------------------------------------
def nearest_labelled(pc, m, lo, hi):
    """Index in (lo, hi) exclusive with a non-None label, nearest to m (or None)."""
    for d in range(0, hi - lo):
        for cand in (m - d, m + d):
            if lo < cand < hi and pc.label(cand) is not None:
                return cand
    return None


def find_boundaries(pc, a, b, ambiguities):
    """Boundaries (first index of a new block) in (a, b], given label(a) != label(b) (both non-None)."""
    la, lb = pc.label(a), pc.label(b)
    if la == lb:
        return []
    if b - a <= 1:
        return [b]
    m = nearest_labelled(pc, (a + b) // 2, a, b)
    if m is None:  # every part strictly between a and b is empty
        ambiguities.append({"between": [a, b], "note": "all intermediate parts empty; boundary placed at b"})
        return [b]
    lm = pc.label(m)
    out = []
    if lm != la:
        out += find_boundaries(pc, a, m, ambiguities)
    if lm != lb:
        out += find_boundaries(pc, m, b, ambiguities)
    return out


def map_product(release, product, cfg, stride, workers, jet_entry_stop, verbose, consts, samples):
    t0 = time.time()
    pc = PartClassifier(release, product, cfg["channel"], jet_entry_stop=jet_entry_stop,
                        read_weights=cfg["weight_check"], verbose=verbose)
    n = len(pc)
    print(f"\n=== {product}: {n} parts, channel={cfg['channel']}, stride={stride}", flush=True)

    coarse = sorted(set(range(0, n, stride)) | {n - 1})
    pc.prefetch(coarse, workers=workers)
    print(f"  coarse scan: {len(coarse)} parts classified ({pc.n_opened} opened, {time.time()-t0:.0f}s)", flush=True)

    # anchors = coarse parts with a label; empty coarse parts are skipped
    anchors = [i for i in coarse if pc.label(i) is not None]
    ambiguities = []
    boundaries = []
    for a, b in zip(anchors[:-1], anchors[1:]):
        if pc.label(a) != pc.label(b):
            boundaries += find_boundaries(pc, a, b, ambiguities)
    boundaries = sorted(set(boundaries))
    print(f"  boundaries at {boundaries} ({pc.n_opened} opened, {time.time()-t0:.0f}s)", flush=True)

    # verify 2 parts on each side of each boundary
    checks = []
    for bnd in boundaries:
        left = [i for i in (bnd - 2, bnd - 1) if i >= 0]
        right = [i for i in (bnd, bnd + 1) if i < n]
        pc.prefetch(left + right, workers=workers)
        ll = [pc.label(i) for i in left]
        rl = [pc.label(i) for i in right]
        l_lab = pc.label(bnd - 1)
        r_lab = pc.label(bnd)
        ok = all(x == l_lab for x in ll if x is not None) and all(x == r_lab for x in rl if x is not None) \
            and l_lab is not None and r_lab is not None and l_lab != r_lab
        checks.append({"boundary": bnd, "left_parts": left, "left_labels": ll,
                       "right_parts": right, "right_labels": rl, "consistent": bool(ok)})

    # ranges
    starts = [0] + boundaries
    ends = [b - 1 for b in boundaries] + [n - 1]
    ranges = []
    for s, e in zip(starts, ends):
        inside = sorted(i for i in pc.cache if s <= i <= e)
        labs = [pc.label(i) for i in inside]
        non_none = [x for x in labs if x is not None]
        lab = non_none[0] if non_none else None
        n_empty = sum(1 for x in labs if x is None)
        p5s = [pc.get(i)["leading_pt_p5"] for i in inside if pc.label(i) is not None]
        p50s = [pc.get(i)["leading_pt_p50"] for i in inside if pc.label(i) is not None]
        r = {
            "first": s, "last": e, "sample": lab, "n_parts": e - s + 1,
        }
        r.update(ppg12_info(lab, consts, samples) if lab else {"ppg12_sample": None})
        r["evidence"] = {
            "classifier": ("per-event leading prompt truth photon pT (truthPhotons, prompt_class==1)"
                           if cfg["channel"] == "photon" else
                           f"per-event leading truth jet pT (truthJets, first {jet_entry_stop} rows)"),
            "rule": f"nearest of {THRESHOLDS[cfg['channel']]} to p5 + {P5_OFFSET[cfg['channel']]}",
            "parts_classified_in_range": len(inside),
            "parts_empty_in_range": n_empty,
            "all_classified_parts_agree": bool(all(x == lab for x in non_none)),
            "leading_pt_p5_min": min(p5s) if p5s else None,
            "leading_pt_p5_max": max(p5s) if p5s else None,
            "leading_pt_p50_min": min(p50s) if p50s else None,
            "leading_pt_p50_max": max(p50s) if p50s else None,
            "parts": [
                {"part": i, "label": pc.label(i),
                 "p5": (round(pc.get(i)["leading_pt_p5"], 3) if pc.label(i) else None),
                 "p50": (round(pc.get(i)["leading_pt_p50"], 3) if pc.label(i) else None),
                 "n_events": pc.get(i)["n_events"], "n_events_grouped": pc.get(i)["n_events_grouped"]}
                for i in inside],
        }
        # event_weight cross-check (p+p): first and last part of the range
        if cfg["weight_check"] and lab:
            exp_tree = TREE_EVENT_WEIGHT_MAX.get(lab)
            hdr_w = samples.get(lab, {}).get("weight_value")
            exp_header = hdr_w * TREE_WEIGHT_FACTOR if hdr_w is not None else None
            hv = header_vs_git_head().get("constants_differing_from_HEAD", {})
            xs_name, ref_name = f"{lab}cross", samples.get(lab, {}).get("weight_den_name")
            exp_head = None
            if ref_name:
                xs_head = hv.get(xs_name, {}).get("HEAD", consts.get(xs_name))
                ref_head = hv.get(ref_name, {}).get("HEAD", consts.get(ref_name))
                if xs_head and ref_head:
                    exp_head = xs_head / ref_head * TREE_WEIGHT_FACTOR
            checked = []
            for i in sorted({s, e}):
                ew = pc.get(i).get("event_weight") or {}
                mx = ew.get("max_nonunit")
                ratio = (mx / exp_tree) if (mx is not None and exp_tree) else None
                # which sample does the observed maximum actually correspond to?
                best = None
                if mx is not None:
                    best = min(TREE_EVENT_WEIGHT_MAX.items(), key=lambda kv: abs(np.log(mx / kv[1])))[0]
                checked.append({"part": i, "max_nonunit": mx, "max_raw": ew.get("max_raw"),
                                "frac_unit": (round(ew["frac_unit"], 4) if "frac_unit" in ew else None),
                                "ratio_to_expected": (round(ratio, 5) if ratio is not None else None),
                                "closest_sample_by_weight": best,
                                "ok": bool(ratio is not None and abs(ratio - 1) < WEIGHT_TOL and best == lab)})
            r["event_weight_check"] = {
                "convention": "max(event_weight != 1) = (sigma_sample/sigma_ref) x 0.7088",
                "expected_tree": exp_tree,
                "expected_from_header_x0.7088": exp_header,
                "tree_over_header": (round(exp_tree / exp_header, 4) if exp_header else None),
                "expected_from_HEAD_header_x0.7088": exp_head,
                "tree_over_HEAD_header": (round(exp_tree / exp_head, 4) if exp_head else None),
                "parts": checked,
                "all_ok": bool(checked and all(c["ok"] for c in checked)),
            }
        ranges.append(r)

    found_order = [r["sample"] for r in ranges]
    result = {
        "product": product,
        "release": release,
        "header_vs_git_head": header_vs_git_head(),
        "files_txt": pc.files_txt,
        "generated": _dt.datetime.now().isoformat(timespec="seconds"),
        "n_parts": n,
        "parts_opened": pc.n_opened,
        "stride": stride,
        "expected_order": cfg["expected_order"],
        "found_order": found_order,
        "order_matches_expected": bool(found_order == cfg["expected_order"]),
        "boundaries": boundaries,
        "boundary_checks": checks,
        "all_boundaries_consistent": bool(all(c["consistent"] for c in checks)),
        "ambiguities": ambiguities,
        "ranges": ranges,
        "seconds": round(time.time() - t0, 1),
    }
    if cfg["weight_check"]:
        result["event_weight_check_all_ok"] = bool(all(r.get("event_weight_check", {}).get("all_ok", False)
                                                       for r in ranges))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", default=DEFAULT_RELEASE)
    ap.add_argument("--products", nargs="+", default=list(PRODUCTS), choices=list(PRODUCTS))
    ap.add_argument("--out-dir", default=os.path.join(PHOTONJET_DIR, "metadata", "sample_maps"))
    ap.add_argument("--stride", type=int, default=40)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--jet-entry-stop", type=int, default=200000)
    ap.add_argument("--header", default=HEADER_PATH, help="CrossSectionWeights.h to quote")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    consts, samples = parse_header(args.header)
    os.makedirs(args.out_dir, exist_ok=True)
    total_opened = 0
    for product in args.products:
        res = map_product(args.release, product, PRODUCTS[product], args.stride, args.workers,
                          args.jet_entry_stop, args.verbose, consts, samples)
        total_opened += res["parts_opened"]
        out = os.path.join(args.out_dir, f"sample_map_{product}.yaml")
        used = dump_yaml(res, out)
        print(f"  ranges:")
        for r in res["ranges"]:
            ev = r["evidence"]
            wc = r.get("event_weight_check")
            wtxt = ""
            if wc:
                wtxt = (f"  w_max={[c['max_nonunit'] for c in wc['parts']]} exp={wc['expected_tree']} "
                        f"ok={wc['all_ok']}")
            print(f"    [{r['first']:5d}, {r['last']:5d}] n={r['n_parts']:5d} {str(r['sample']):9s} "
                  f"p5 {ev['leading_pt_p5_min']:.2f}-{ev['leading_pt_p5_max']:.2f} "
                  f"({ev['parts_classified_in_range']} parts checked, agree={ev['all_classified_parts_agree']}){wtxt}")
        print(f"  order {res['found_order']} matches expected: {res['order_matches_expected']}; "
              f"boundaries consistent: {res['all_boundaries_consistent']}; "
              f"opened {res['parts_opened']}/{res['n_parts']} parts in {res['seconds']}s")
        print(f"  wrote {out} ({used})")
    print(f"\nTotal parts opened: {total_opened}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
