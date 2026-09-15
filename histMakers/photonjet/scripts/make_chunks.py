#!/usr/bin/env python3
"""Assign original ROOT parts to jobs; seal exact identities and coverage in jobs.json."""
import argparse
import os
from pathlib import Path
from inputs import (PJ, PRODUCTS, BACKEND, SCHEMA, config, release_label, master_parts, indices,
                    plan_source, label, md5, master_digest, atomic_text, save, bin_layout,
                    code_digests, dependencies, PLAN_FORMAT_VERSION, INPUT_VERIFICATION)

DEFAULT_SIZES = {"data": 200, "sim_signal": 40, "sim_inclusive": 40}


def make(args):
    if getattr(args, "store", None) is not None:
        raise ValueError("prepared stores are retired; chunks now reads original ROOT inputs directly")
    cfg = config(args.config)
    variation = label(str(cfg["output"]["var_type"]))
    release = release_label(cfg, args.release_id)
    config_path = Path(args.config).resolve()
    directory = Path(args.outdir).resolve() / variation
    sizes = dict(DEFAULT_SIZES)
    for option in filter(None, args.chunk_size.split(",")):
        product, size = option.split("=")
        if product not in PRODUCTS or int(size) < 1:
            raise ValueError("chunk sizes require a known product and positive integer")
        sizes[product] = int(size)
    plan = dict(format_version=PLAN_FORMAT_VERSION, backend=BACKEND, schema_version=SCHEMA,
                input_verification=INPUT_VERIFICATION,
                config=str(config_path), config_md5=md5(config_path), layout=bin_layout(cfg),
                code=code_digests(), release_id=release, variation=variation, products={})
    writes = {}
    all_jobs = []
    products = args.products.split(",")
    if len(set(products)) != len(products) or any(p not in PRODUCTS for p in products):
        raise ValueError("unknown or duplicate product selection")
    for product in products:
        master, paths = master_parts(cfg, product)
        selected = indices(args.indices, len(paths))
        records = []
        source_paths = set()
        for index in selected:
            part = plan_source(paths[index], index, cfg, product, release)
            identity = Path(part["source_path"]).resolve(strict=True)
            if identity in source_paths:
                raise ValueError("duplicate original path in planned coverage")
            source_paths.add(identity)
            records.append(part)
        jobs = []
        for first in range(0, len(records), sizes[product]):
            tag = f"{product}_{len(jobs):03d}"
            parts = records[first:first+sizes[product]]
            chunk = directory / f"{tag}.list"
            text = "".join(f"{part['original_part_index']} {part['source_path']}\n" for part in parts)
            writes[chunk] = text
            job = dict(config=str(config_path), product=product, input=str(chunk), tag=tag, parts=parts)
            jobs.append(job)
            row = " ".join([str(config_path), product, str(chunk), tag])
            if any(any(ch.isspace() for ch in value) for value in (str(config_path), str(chunk))):
                raise ValueError("Condor job paths must contain no whitespace")
            all_jobs.append(row)
        writes[directory/f"jobs_{product}.list"] = "".join(" ".join([j["config"],product,j["input"],j["tag"]])+"\n" for j in jobs)
        plan["products"][product] = dict(master=str(master), master_md5=master_digest(cfg, product), indices=selected,
            scope="full" if sorted(selected)==list(range(len(paths))) else "explicit_subset",
            dependencies=dependencies(cfg, product), jobs=jobs)
    writes[directory/"jobs_all.list"] = "\n".join(all_jobs)+"\n"
    # Inputs have been validated for every product. The plan is the last write;
    # its list hashes make an interrupted update fail validation.
    for product, record in plan["products"].items():
        if md5(record["master"]) != record["master_md5"]:
            raise ValueError("master list changed during chunk planning")
        for source, digest in record["dependencies"].items():
            if md5(source) != digest:
                raise ValueError("dependency changed during chunk planning: " + source)
        for job in record["jobs"]:
            for part in job["parts"]:
                if plan_source(part["source_path"], part["original_part_index"], cfg, product, release) != part:
                    raise ValueError("original input metadata changed during chunk planning")
    if code_digests() != plan["code"]:
        raise ValueError("maker code changed during chunk planning")
    if md5(config_path) != plan["config_md5"]:
        raise ValueError("config changed during chunk planning")
    for path, text in writes.items():
        atomic_text(path, text)
    plan["list_md5"] = {str(path): md5(path) for path in writes}
    save(directory/"jobs.json", plan)
    print(f"[make_chunks] {len(all_jobs)} serial original-input jobs; exact coverage: {directory/'jobs.json'}", flush=True)
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config")
    parser.add_argument("--products", default=",".join(PRODUCTS))
    parser.add_argument("--chunk-size", default="")
    parser.add_argument("--store", help=argparse.SUPPRESS)
    parser.add_argument("--outdir", default=os.environ.get("PJ_CHUNK_DIR", str(PJ/"chunks")))
    parser.add_argument("--release-id")
    parser.add_argument("--indices", help="explicit original indices/half-open ranges, e.g. 0,4:7")
    make(parser.parse_args())


if __name__ == "__main__":
    main()
