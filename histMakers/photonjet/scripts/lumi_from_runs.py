#!/usr/bin/env python3
"""Run census + integrated luminosity for a PhotonJetTrees data product.

Reads ONLY the ``run`` branch of the ``events`` tree of every part listed in the
product's ``files.txt`` (never lists the parts/ directory), builds a
run -> event-count census and, for p+p, sums Joey's per-run Bit30Corr /
Bit30UC luminosities (all-z list and 60 cm list) over the runs present in the
tree.  Also reports the per-period split at run 51274 (0 mrad: run < 51274,
1.5 mrad: run >= 51274).

Outputs (``<out-prefix>`` is a path prefix, e.g. ``efficiencytool/photonjet/metadata/luminosity/``):
  <out-prefix>run_list_<pp|auau>.txt   run  n_events        (one run per line)
  <out-prefix>lumi_<pp|auau>.yaml      census + lumi summary

Lumi list format (lumi/allzLumi_fromJoey.list, lumi/60cmLumi_fromJoey.list):
  header lines, then rows
  RN Bit10Corr Bit10UC Bit18Corr Bit18UC Bit22Corr Bit22UC Bit30Corr Bit30UC   (pb^-1)
  A file may contain several concatenated "Processing ..." blocks (the 60 cm
  list has two passes with slightly different values).  Blocks are parsed
  separately and the LAST one is used by default (--lumi-pass first|last).
  Short rows (e.g. "51953 0 0 0 0 0 0") are zero-lumi runs.

Example:
  python3 efficiencytool/photonjet/scripts/lumi_from_runs.py --product data/pp
  nohup python3 efficiencytool/photonjet/scripts/lumi_from_runs.py --product data/auau \
        > efficiencytool/photonjet/logs/lumi_auau.log 2>&1 &
"""
import argparse
import datetime as _dt
import os
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import uproot

DEFAULT_RELEASE = "/sphenix/tg/tg01/bulk/jbennett/PhotonJetTrees_v1"
PERIOD_SPLIT_RUN = 51274  # 0 mrad: run < 51274 ; 1.5 mrad: run >= 51274
LUMI_COLUMNS = ["Bit10Corr", "Bit10UC", "Bit18Corr", "Bit18UC",
                "Bit22Corr", "Bit22UC", "Bit30Corr", "Bit30UC"]

PHOTONJET_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(os.path.dirname(PHOTONJET_DIR))


def repo_path(p):
    """Absolute path; relative paths are taken relative to the repo root."""
    if p is None:
        return None
    return p if os.path.isabs(p) else os.path.join(REPO_ROOT, p)


# ----------------------------------------------------------------------------
# YAML writer (ruamel -> PyYAML -> minimal plain text)
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

    def _emit(o, ind, out):
        pad = "  " * ind
        if isinstance(o, dict):
            for k, v in o.items():
                if isinstance(v, (dict, list)) and v:
                    out.append(f"{pad}{k}:")
                    _emit(v, ind + 1, out)
                else:
                    out.append(f"{pad}{k}: {_scalar(v)}")
        elif isinstance(o, list):
            for v in o:
                if isinstance(v, (dict, list)) and v:
                    out.append(f"{pad}-")
                    _emit(v, ind + 1, out)
                else:
                    out.append(f"{pad}- {_scalar(v)}")

    def _scalar(v):
        if v is None:
            return "null"
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, (list, dict)) and not v:
            return "[]" if isinstance(v, list) else "{}"
        if isinstance(v, str):
            return repr(v) if any(c in v for c in ":#{}[],&*!|>'\"%@`") else v
        return str(v)

    lines = []
    _emit(obj, 0, lines)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return "plain"


# ----------------------------------------------------------------------------
# Lumi list parsing
# ----------------------------------------------------------------------------
def parse_lumi_list(path):
    """Return list of blocks; each block is {run: {column: value}} (pb^-1)."""
    blocks = []
    cur = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("Processing"):
                cur = {}
                blocks.append(cur)
                continue
            tok = s.split()
            if not tok[0].isdigit():
                continue  # header / Totals / "(int) 0"
            if cur is None:
                cur = {}
                blocks.append(cur)
            run = int(tok[0])
            vals = [float(x) for x in tok[1:]]
            row = {c: (vals[i] if i < len(vals) else 0.0) for i, c in enumerate(LUMI_COLUMNS)}
            row["_short_row"] = len(vals) < len(LUMI_COLUMNS)
            cur[run] = row
    return blocks


def block_sum(block, col, runs=None):
    if runs is None:
        return float(sum(r[col] for r in block.values()))
    return float(sum(block[r][col] for r in runs if r in block))


# ----------------------------------------------------------------------------
# Tree reading
# ----------------------------------------------------------------------------
def read_part_runs(idx, path, tree_name="events", branch="run"):
    """Return (idx, n_entries, Counter{run: n_events}). Only the run branch is read."""
    with uproot.open(path) as f:
        t = f[tree_name]
        n = int(t.num_entries)
        if n == 0:
            return idx, 0, Counter()
        runs = t[branch].array(library="np")
        u, c = np.unique(runs, return_counts=True)
        return idx, n, Counter({int(a): int(b) for a, b in zip(u, c)})


def read_files_txt(path):
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip()]


def census(files, workers=8, progress_every=500, max_parts=None):
    if max_parts is not None:
        files = files[:max_parts]
    n_parts = len(files)
    per_run = Counter()
    n_nonempty = 0
    n_events = 0
    failed = []
    part_nevents = np.zeros(n_parts, dtype=np.int64)
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(read_part_runs, i, p): (i, p) for i, p in enumerate(files)}
        for fut in as_completed(futs):
            i, p = futs[fut]
            done += 1
            try:
                idx, n, cnt = fut.result()
            except Exception as e:  # keep going, report at the end
                failed.append({"part": i, "path": p, "error": f"{type(e).__name__}: {e}"})
            else:
                part_nevents[idx] = n
                if n > 0:
                    n_nonempty += 1
                    n_events += n
                    per_run.update(cnt)
            if done % progress_every == 0 or done == n_parts:
                el = time.time() - t0
                print(f"  [{done}/{n_parts}] parts read, {n_events} events, {len(per_run)} runs, "
                      f"{el:.0f}s elapsed", flush=True)
    return {
        "n_parts": n_parts,
        "n_parts_nonempty": n_nonempty,
        "n_parts_failed": len(failed),
        "failed_parts": failed,
        "n_events": int(n_events),
        "per_run": dict(sorted(per_run.items())),
        "part_nevents": part_nevents,
        "seconds": time.time() - t0,
    }


# ----------------------------------------------------------------------------
def period_of(run):
    return "0mrad" if run < PERIOD_SPLIT_RUN else "1p5mrad"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--product", default="data/pp", choices=["data/pp", "data/auau"],
                    help="release product (default data/pp)")
    ap.add_argument("--release", default=DEFAULT_RELEASE, help="PhotonJetTrees release directory")
    ap.add_argument("--lumi-list", default="lumi/allzLumi_fromJoey.list",
                    help="all-z per-run lumi list (repo-relative or absolute); p+p only")
    ap.add_argument("--lumi-list-60cm", default="lumi/60cmLumi_fromJoey.list",
                    help="60 cm-fiducial per-run lumi list (repo-relative or absolute); p+p only")
    ap.add_argument("--lumi-pass", default="last", choices=["first", "last"],
                    help="which 'Processing' block of a multi-pass lumi list to use (default last)")
    ap.add_argument("--out-prefix", default=os.path.join(PHOTONJET_DIR, "metadata", "luminosity", ""),
                    help="output path prefix (repo-relative or absolute)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--progress-every", type=int, default=500)
    ap.add_argument("--max-parts", type=int, default=None, help="debug: only read the first N parts")
    args = ap.parse_args()

    tag = "pp" if args.product == "data/pp" else "auau"
    is_pp = tag == "pp"
    files_txt = os.path.join(args.release, args.product, "files.txt")
    files = read_files_txt(files_txt)
    out_prefix = repo_path(args.out_prefix)
    out_dir = os.path.dirname(out_prefix) if out_prefix.endswith("/") else os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    run_list_path = f"{out_prefix}run_list_{tag}.txt"
    yaml_path = f"{out_prefix}lumi_{tag}.yaml"

    print(f"[lumi_from_runs] product={args.product} parts={len(files)} workers={args.workers}", flush=True)
    c = census(files, workers=args.workers, progress_every=args.progress_every, max_parts=args.max_parts)
    per_run = c["per_run"]
    runs = sorted(per_run)
    n_events = c["n_events"]
    print(f"[lumi_from_runs] done in {c['seconds']:.0f}s: {c['n_parts_nonempty']}/{c['n_parts']} non-empty parts, "
          f"{n_events} events, {len(runs)} runs"
          + (f" [{runs[0]}..{runs[-1]}]" if runs else "") + f", {c['n_parts_failed']} failed parts", flush=True)

    # ---- run list --------------------------------------------------------
    with open(run_list_path, "w") as f:
        f.write(f"# run  n_events   ({args.product}, tree=events, branch=run, {files_txt})\n")
        for r in runs:
            f.write(f"{r} {per_run[r]}\n")

    # ---- summary ----------------------------------------------------------
    per_part = c["part_nevents"]
    summary = {
        "product": args.product,
        "release": args.release,
        "files_txt": files_txt,
        "tree": "events",
        "branch": "run",
        "generated": _dt.datetime.now().isoformat(timespec="seconds"),
        "n_parts": c["n_parts"],
        "n_parts_nonempty": c["n_parts_nonempty"],
        "n_parts_failed": c["n_parts_failed"],
        "n_events": n_events,
        "events_per_part": {
            "min": int(per_part.min()) if len(per_part) else 0,
            "max": int(per_part.max()) if len(per_part) else 0,
            "mean": float(per_part.mean()) if len(per_part) else 0.0,
        },
        "n_runs": len(runs),
        "run_min": runs[0] if runs else None,
        "run_max": runs[-1] if runs else None,
        "run_list_file": run_list_path,
    }
    if c["failed_parts"]:
        summary["failed_parts"] = c["failed_parts"]

    if is_pp:
        allz_path = repo_path(args.lumi_list)
        cm60_path = repo_path(args.lumi_list_60cm)
        allz_blocks = parse_lumi_list(allz_path)
        cm60_blocks = parse_lumi_list(cm60_path)
        pick = (lambda b: b[-1]) if args.lumi_pass == "last" else (lambda b: b[0])
        allz = pick(allz_blocks)
        cm60 = pick(cm60_blocks)

        present = [r for r in runs if r in allz]
        missing = [r for r in runs if r not in allz]
        missing_60 = [r for r in runs if r not in cm60]
        missing_ev = sum(per_run[r] for r in missing)
        not_in_tree = sorted(r for r in allz if r not in per_run)

        summary.update({
            "lumi_list_allz": allz_path,
            "lumi_list_60cm": cm60_path,
            "lumi_pass_used": args.lumi_pass,
            "lumi_list_allz_blocks": [
                {"block": i + 1, "n_runs": len(b), "bit30corr_sum_all_runs_pb": round(block_sum(b, "Bit30Corr"), 4)}
                for i, b in enumerate(allz_blocks)],
            "lumi_list_60cm_blocks": [
                {"block": i + 1, "n_runs": len(b), "bit30corr_sum_all_runs_pb": round(block_sum(b, "Bit30Corr"), 4)}
                for i, b in enumerate(cm60_blocks)],
            # sums over the runs actually present in the tree
            "lumi_bit30corr_allz_pb": round(block_sum(allz, "Bit30Corr", runs), 4),
            "lumi_bit30uc_allz_pb": round(block_sum(allz, "Bit30UC", runs), 4),
            "lumi_bit30corr_60cm_pb": round(block_sum(cm60, "Bit30Corr", runs), 4),
            "lumi_bit30uc_60cm_pb": round(block_sum(cm60, "Bit30UC", runs), 4),
            # whole-list totals for reference
            "lumi_list_total_bit30corr_allz_pb": round(block_sum(allz, "Bit30Corr"), 4),
            "lumi_list_total_bit30corr_60cm_pb": round(block_sum(cm60, "Bit30Corr"), 4),
            "runs_in_tree_and_lumi_list": len(present),
            "runs_missing_from_lumi_list": {
                "n": len(missing),
                "n_events": int(missing_ev),
                "event_share": (float(missing_ev) / n_events) if n_events else 0.0,
                "runs": missing,
            },
            "runs_missing_from_lumi_list_60cm": {
                "n": len(missing_60),
                "runs": missing_60,
            },
            "runs_in_lumi_list_not_in_tree": {
                "n": len(not_in_tree),
                "lumi_bit30corr_allz_pb": round(block_sum(allz, "Bit30Corr", not_in_tree), 4),
                "lumi_bit30corr_60cm_pb": round(block_sum(cm60, "Bit30Corr", [r for r in not_in_tree if r in cm60]), 4),
                "runs": not_in_tree,
            },
        })
        periods = {}
        for name, sel in (("0mrad", lambda r: r < PERIOD_SPLIT_RUN), ("1p5mrad", lambda r: r >= PERIOD_SPLIT_RUN)):
            pr = [r for r in runs if sel(r)]
            ev = sum(per_run[r] for r in pr)
            periods[name] = {
                "run_selection": (f"run < {PERIOD_SPLIT_RUN}" if name == "0mrad" else f"run >= {PERIOD_SPLIT_RUN}"),
                "n_runs": len(pr),
                "run_min": pr[0] if pr else None,
                "run_max": pr[-1] if pr else None,
                "n_events": int(ev),
                "event_share": (float(ev) / n_events) if n_events else 0.0,
                "lumi_bit30corr_allz_pb": round(block_sum(allz, "Bit30Corr", pr), 4),
                "lumi_bit30uc_allz_pb": round(block_sum(allz, "Bit30UC", pr), 4),
                "lumi_bit30corr_60cm_pb": round(block_sum(cm60, "Bit30Corr", pr), 4),
                "lumi_bit30uc_60cm_pb": round(block_sum(cm60, "Bit30UC", pr), 4),
                "runs_missing_from_lumi_list": len([r for r in pr if r not in allz]),
            }
        summary["periods"] = periods

        print(f"[lumi_from_runs] Bit30Corr allz = {summary['lumi_bit30corr_allz_pb']} pb^-1 "
              f"(UC {summary['lumi_bit30uc_allz_pb']}), 60cm = {summary['lumi_bit30corr_60cm_pb']} pb^-1 "
              f"(UC {summary['lumi_bit30uc_60cm_pb']}); list totals allz {summary['lumi_list_total_bit30corr_allz_pb']}, "
              f"60cm {summary['lumi_list_total_bit30corr_60cm_pb']}")
        for name, p in periods.items():
            print(f"  {name:8s}: {p['n_runs']} runs [{p['run_min']}..{p['run_max']}], {p['n_events']} events "
                  f"({100*p['event_share']:.1f}%), Bit30Corr allz {p['lumi_bit30corr_allz_pb']}, "
                  f"60cm {p['lumi_bit30corr_60cm_pb']} pb^-1")
        print(f"  runs in tree missing from lumi list: {len(missing)} ({missing_ev} events, "
              f"{100*summary['runs_missing_from_lumi_list']['event_share']:.3f}%): {missing}")
        print(f"  runs in lumi list not in tree: {len(not_in_tree)} "
              f"({summary['runs_in_lumi_list_not_in_tree']['lumi_bit30corr_allz_pb']} pb^-1 allz)")
    else:
        summary["note"] = "Au+Au: no per-run lumi list available; run census only."

    used = dump_yaml(summary, yaml_path)
    print(f"[lumi_from_runs] wrote {run_list_path}\n[lumi_from_runs] wrote {yaml_path} ({used})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
