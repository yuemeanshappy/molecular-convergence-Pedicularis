#!/usr/bin/env python3
"""Classify Trinity transcripts using existing GetOrganelle pt/mt paths."""
import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

VERSION = "1.0.0"
KINDS = ("pt_all", "mt_all", "mt_pt_intersect", "pt_only", "mt_only", "mt_pt_union", "nucl")
GENE = re.compile(r"(?P<gene>TRINITY_\S+_g[0-9]+)_i[0-9]+\Z")
PATH_HEADER = re.compile(r"(?:contig_[1-9][0-9]*--)?(?P<path>[1-9][0-9]*[+-](?:,[1-9][0-9]*[+-])*)(?:\(circular\))?\Z")


class InputError(ValueError):
    """An input is missing, ambiguous, or inconsistent."""


@dataclass(frozen=True)
class Sample:
    sample: str
    trinity: Path
    name_map: Path
    pt: Path
    mt: Path


def fasta_records(path):
    """Yield full header and sequence; line wrapping is not part of the sequence."""
    label = None
    parts = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if label is not None:
                    if not parts:
                        raise InputError(f"{path}: empty sequence for {label!r}")
                    yield label, "".join(parts)
                label = line[1:]
                if not label or label != label.strip() or "\t" in label:
                    raise InputError(f"{path}:{line_no}: malformed FASTA header")
                parts = []
            elif line.strip():
                if label is None or any(c.isspace() for c in line):
                    raise InputError(f"{path}:{line_no}: sequence before header or whitespace in sequence")
                parts.append(line)
        if label is not None:
            if not parts:
                raise InputError(f"{path}: empty sequence for {label!r}")
            yield label, "".join(parts)


def parse_path_header(label):
    """Parse numeric GFA paths; orientation affects membership, never raw sequence."""
    match = PATH_HEADER.fullmatch(label)
    if not match:
        raise InputError(f"unsupported GetOrganelle path header {label!r}; expected "
                         "[contig_N--]ID+/ID- comma-separated path, optionally ending (circular)")
    return [token[:-1] for token in match.group("path").split(",")]


def load_raw(path):
    records = []
    index = {}
    genes = defaultdict(set)
    for label, seq in fasta_records(path):
        sid = label.split(" ", 1)[0]
        match = GENE.fullmatch(sid)
        if not match:
            raise InputError(f"{path}: invalid Trinity isoform ID {sid!r}; expected TRINITY_..._gN_iN")
        if sid in index:
            raise InputError(f"{path}: duplicate Trinity ID {sid!r}")
        idx = len(records)
        index[sid] = idx
        records.append((label, seq, match.group("gene")))
        genes[match.group("gene")].add(idx)
    if not records:
        raise InputError(f"{path}: Trinity FASTA is empty")
    return records, index, genes


def load_map(path, records, index):
    ctgs = {}
    seen_ids = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not re.fullmatch(r"[1-9][0-9]*", fields[1]):
                raise InputError(f"{path}:{line_no}: expected FULL_FASTA_HEADER<TAB>positive_numeric_contig_ID")
            label, cid = fields
            sid = label.split(" ", 1)[0]
            if sid not in index:
                raise InputError(f"{path}:{line_no}: unknown Trinity ID {sid!r}")
            if cid in ctgs or sid in seen_ids:
                raise InputError(f"{path}:{line_no}: duplicate/conflicting mapping for {cid!r} or {sid!r}")
            idx = index[sid]
            if label != records[idx][0]:
                raise InputError(f"{path}:{line_no}: full header differs from Trinity input for {sid!r}")
            ctgs[cid] = idx
            seen_ids.add(sid)
    if len(seen_ids) != len(records):
        raise InputError(f"{path}: incomplete name_map: {len(seen_ids)} mapped / {len(records)} Trinity records")
    return ctgs


def load_seeds(path, ctgs, allow_empty):
    ids = set()
    labels = set()
    for label, _seq in fasta_records(path):
        if label in labels:
            raise InputError(f"{path}: duplicate scaffold header {label!r}")
        labels.add(label)
        for cid in parse_path_header(label):
            if cid not in ctgs:
                raise InputError(f"{path}: unknown contig ID {cid!r} in {label!r}")
            ids.add(cid)
    if not labels and not allow_empty:
        raise InputError(f"{path}: empty seed FASTA; use --allow-empty-seeds only for intentional no-hit inputs")
    return ids, len(labels)


def classify(spec, allow_empty=False):
    records, index, genes = load_raw(spec.trinity)
    ctgs = load_map(spec.name_map, records, index)
    pt, pt_headers = load_seeds(spec.pt, ctgs, allow_empty)
    mt, mt_headers = load_seeds(spec.mt, ctgs, allow_empty)

    def expand(seeds):
        result = set()
        for cid in seeds:
            result.update(genes[records[ctgs[cid]][2]])
        return result

    chosen = {"pt_all": expand(pt), "mt_all": expand(mt),
              "mt_pt_intersect": expand(pt & mt), "pt_only": expand(pt - mt),
              "mt_only": expand(mt - pt), "mt_pt_union": expand(pt | mt)}
    chosen["nucl"] = set(range(len(records))) - chosen["mt_pt_union"]
    overlap = chosen["pt_all"] & chosen["mt_all"]
    only_overlap = chosen["pt_only"] & chosen["mt_only"]
    report = {"tool_version": VERSION, "sample": spec.sample,
              "semantics": "contig sets first, then all isoforms of each Trinity gene",
              "inputs": {name: str(getattr(spec, name)) for name in ("trinity", "name_map", "pt", "mt")},
              "raw_records": len(records), "raw_genes": len(genes),
              "seed_counts": {"pt_scaffolds": pt_headers, "mt_scaffolds": mt_headers,
                              "pt_contigs": len(pt), "mt_contigs": len(mt), "shared_contigs": len(pt & mt)},
              "classes": {kind: {"records": len(chosen[kind]),
                                  "bases": sum(len(records[i][1]) for i in chosen[kind])} for kind in KINDS},
              "overlap": {"pt_all_and_mt_all": len(overlap), "pt_only_and_mt_only": len(only_overlap),
                          "shared_genes": sorted({records[i][2] for i in overlap})},
              "warnings": []}
    if overlap:
        report["warnings"].append(f"pt_all and mt_all share {len(overlap)} isoforms after gene expansion; "
                                  f"pt_only and mt_only share {len(only_overlap)}. These are historical, non-exclusive classes.")
    if not pt or not mt:
        report["warnings"].append("An explicitly allowed seed input is empty; verify that no hits was intentional.")
    return records, chosen, report


def output_paths(sample, out_dir):
    return [out_dir / f"{sample}.{kind}.fasta" for kind in KINDS] + [out_dir / f"{sample}.classification.json"]


def preflight(specs, out_dir):
    samples = set()
    inputs = set()
    for spec in specs:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", spec.sample):
            raise InputError(f"invalid sample name {spec.sample!r}; use letters, numbers, dot, underscore, hyphen")
        if spec.sample in samples:
            raise InputError(f"duplicate sample name {spec.sample!r}")
        samples.add(spec.sample)
        for key in ("trinity", "name_map", "pt", "mt"):
            path = getattr(spec, key)
            if not path.is_file():
                raise InputError(f"{spec.sample}: missing/non-file {key} input: {path}")
            inputs.add(path.resolve())
    for spec in specs:
        for path in output_paths(spec.sample, out_dir):
            if path.exists() or path.is_symlink() or path.resolve() in inputs:
                raise InputError(f"refusing to overwrite existing output or input: {path}")


def write_outputs(spec, out_dir, records, chosen, report):
    out_dir.mkdir(parents=True, exist_ok=True)
    made = []
    try:
        for kind, path in zip(KINDS, output_paths(spec.sample, out_dir)):
            with path.open("x", encoding="utf-8", newline="\n") as handle:
                made.append(path)
                for idx in sorted(chosen[kind]):
                    label, seq, _gene = records[idx]
                    handle.write(f">{label}\n{seq}\n")
        path = output_paths(spec.sample, out_dir)[-1]
        with path.open("x", encoding="utf-8", newline="\n") as handle:
            made.append(path)
            json.dump(report, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
    except BaseException:
        # Only remove files exclusively created by this sample's current attempt.
        for path in made:
            path.unlink()
        raise


def read_manifest(path):
    specs = []
    fields = ["sample", "trinity", "name_map", "pt", "mt"]
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != fields:
            raise InputError("manifest header must be exactly: " + "\t".join(fields))
        for line_no, row in enumerate(reader, 2):
            if None in row or any(not row.get(key) for key in fields):
                raise InputError(f"{path}:{line_no}: missing or extra manifest fields")
            paths = [Path(row[key]).expanduser() for key in fields[1:]]
            paths = [p if p.is_absolute() else path.resolve().parent / p for p in paths]
            specs.append(Sample(row["sample"], *paths))
    if not specs:
        raise InputError("manifest has no samples")
    return specs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    modes = parser.add_subparsers(dest="mode", required=True)
    single = modes.add_parser("single", help="process one sample with explicit input paths")
    single.add_argument("--sample", required=True)
    for key in ("trinity", "name-map", "pt", "mt"):
        single.add_argument("--" + key, required=True, type=Path)
    batch = modes.add_parser("batch", help="process a TSV manifest sequentially")
    batch.add_argument("--manifest", required=True, type=Path)
    for subparser in (single, batch):
        subparser.add_argument("--out-dir", required=True, type=Path)
        subparser.add_argument("--allow-empty-seeds", action="store_true",
                               help="explicitly permit existing, empty pt/mt FASTA inputs")
    args = parser.parse_args(argv)
    try:
        specs = (read_manifest(args.manifest) if args.mode == "batch" else
                 [Sample(args.sample, args.trinity, args.name_map, args.pt, args.mt)])
        preflight(specs, args.out_dir)
        for spec in specs:
            records, chosen, report = classify(spec, args.allow_empty_seeds)
            write_outputs(spec, args.out_dir, records, chosen, report)
            print(spec.sample + "\t" + "\t".join(f"{kind}={len(chosen[kind])}" for kind in KINDS))
            for message in report["warnings"]:
                print(f"WARNING [{spec.sample}]: {message}", file=sys.stderr)
            del records, chosen, report
        return 0
    except (OSError, UnicodeError, InputError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
