#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def safe_int(x: str) -> int:
    try:
        return int(float(x))
    except Exception:
        return 0


def safe_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return 0.0


def write_tsv(path: Path, header: List[str], rows: List[List[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv(path: Path):
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_genome_map(path: Path) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            otu, genome = line.split("\t", 1)
            rows.append((otu, genome))
    return rows


def read_target_taxids(path: Path) -> set[str]:
    keep = set()
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if line:
                keep.add(line)
    return keep


def read_otu_to_taxid_from_target_groups(path: Path) -> Dict[str, str]:
    otu_to_taxid: Dict[str, str] = {}
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("@") or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            _, otu, taxid, _ = parts[:4]
            otu_to_taxid.setdefault(otu, taxid)
    return otu_to_taxid


def parse_biobox_binning(binning_file: Path) -> Dict[str, List[str]]:
    bin_to_contigs: Dict[str, List[str]] = defaultdict(list)
    with open(binning_file, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("@") or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            contig_id, bin_id = parts[0], parts[1]
            bin_to_contigs[bin_id].append(contig_id)
    return dict(bin_to_contigs)


def iter_fasta_records(fasta_path: Path):
    header = None
    seq_chunks: List[str] = []
    with open(fasta_path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def build_fasta_index(fasta_path: Path) -> Dict[str, str]:
    return {header: seq for header, seq in iter_fasta_records(fasta_path)}


def write_fasta_record(handle, header: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i+width] + "\n")


def generate_mags_from_binning(
    binning_file: Path,
    assembly_fasta: Path,
    output_dir: Path,
    overwrite: bool = False,
    report_every: int = 200,
) -> Dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)

    if not binning_file.exists():
        raise FileNotFoundError(f"Binning file not found: {binning_file}")
    if not assembly_fasta.exists():
        raise FileNotFoundError(f"Assembly FASTA not found: {assembly_fasta}")

    bin_to_contigs = parse_biobox_binning(binning_file)
    fasta_index = build_fasta_index(assembly_fasta)

    written_bins = 0
    written_contigs = 0
    missing_contigs = 0
    missing_examples: List[str] = []

    for i, (bin_id, contigs) in enumerate(bin_to_contigs.items(), start=1):
        out_fasta = output_dir / f"{bin_id}.fasta"
        if out_fasta.exists() and not overwrite:
            continue

        wrote_any = False
        with open(out_fasta, "w") as out_handle:
            for contig_id in contigs:
                seq = fasta_index.get(contig_id)
                if seq is None:
                    missing_contigs += 1
                    if len(missing_examples) < 20:
                        missing_examples.append(contig_id)
                    continue
                write_fasta_record(out_handle, contig_id, seq)
                wrote_any = True
                written_contigs += 1

        if wrote_any:
            written_bins += 1

        if report_every and i % report_every == 0:
            print(f"[generate_mags] processed {i} bins from {binning_file.name}")

    return {
        "binning_file": str(binning_file),
        "assembly_fasta": str(assembly_fasta),
        "output_dir": str(output_dir),
        "bins_in_binning": len(bin_to_contigs),
        "written_bins": written_bins,
        "written_contigs": written_contigs,
        "missing_contigs": missing_contigs,
        "missing_examples": missing_examples,
    }


def fasta_lengths_quastlike(fasta_path: Path, min_len: int = 500) -> List[int]:
    lengths: List[int] = []
    current = 0
    with open(fasta_path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current >= min_len:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)
        if current >= min_len:
            lengths.append(current)
    return lengths


def n50_from_lengths(lengths: List[int]) -> int:
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    half = total / 2
    running = 0
    for x in lengths:
        running += x
        if running >= half:
            return x
    return 0


def fasta_stats_quastlike(fasta_path: Path, min_len: int = 500) -> Tuple[int, int, int]:
    lengths = fasta_lengths_quastlike(fasta_path, min_len=min_len)
    if not lengths:
        return (0, 0, 0)
    return (len(lengths), sum(lengths), n50_from_lengths(lengths))


def build_summary_full(
    genome_map: Path,
    hyb_dir: Path,
    short_dir: Path,
    out_path: Path,
) -> Dict[str, object]:
    header = [
        "OTU", "reference_genome", "hybrid_bin", "short_bin",
        "ref_contigs", "ref_total_bp", "ref_N50",
        "hyb_contigs", "hyb_total_bp", "hyb_N50",
        "short_contigs", "short_total_bp", "short_N50",
    ]

    rows: List[List[object]] = []
    missing: List[Tuple[str, bool, bool, bool]] = []

    for otu, genome in read_genome_map(genome_map):
        ref = Path(genome)
        hyb = hyb_dir / f"{otu}.fasta"
        short = short_dir / f"{otu}.fasta"

        if not ref.exists() or not hyb.exists() or not short.exists():
            missing.append((otu, ref.exists(), hyb.exists(), short.exists()))
            continue

        ref_stats = fasta_stats_quastlike(ref)
        hyb_stats = fasta_stats_quastlike(hyb)
        short_stats = fasta_stats_quastlike(short)

        rows.append([
            otu, str(ref), str(hyb), str(short),
            *ref_stats, *hyb_stats, *short_stats
        ])

    write_tsv(out_path, header, rows)
    return {
        "output": str(out_path),
        "rows": len(rows),
        "missing_triplets": len(missing),
        "missing_examples": missing[:20],
    }


def build_summary_full_with_taxid(
    summary_full: Path,
    otu_to_taxid: Dict[str, str],
    otu_to_taxid_out: Path,
    out_path: Path,
) -> Dict[str, object]:
    with open(otu_to_taxid_out, "w") as handle:
        for otu, taxid in sorted(otu_to_taxid.items()):
            handle.write(f"{otu}\t{taxid}\n")

    summary_rows = read_tsv(summary_full)
    header = [
        "OTU", "TAXID", "reference_genome", "hybrid_bin", "short_bin",
        "ref_contigs", "ref_total_bp", "ref_N50",
        "hyb_contigs", "hyb_total_bp", "hyb_N50",
        "short_contigs", "short_total_bp", "short_N50",
    ]

    rows: List[List[object]] = []
    for row in summary_rows:
        otu = row["OTU"]
        rows.append([
            otu,
            otu_to_taxid.get(otu, "NA"),
            row["reference_genome"],
            row["hybrid_bin"],
            row["short_bin"],
            row["ref_contigs"],
            row["ref_total_bp"],
            row["ref_N50"],
            row["hyb_contigs"],
            row["hyb_total_bp"],
            row["hyb_N50"],
            row["short_contigs"],
            row["short_total_bp"],
            row["short_N50"],
        ])

    write_tsv(out_path, header, rows)
    return {"output": str(out_path), "rows": len(rows)}


def build_summary_ranked(
    summary_full_with_taxid: Path,
    out_path: Path,
) -> Dict[str, object]:
    summary = read_tsv(summary_full_with_taxid)
    if not summary:
        raise ValueError("summary_full_with_taxid.tsv is empty")

    header = list(summary[0].keys()) + ["ref/hyb", "hyb/short"]
    rows: List[List[object]] = []

    for row in summary:
        ref_n50 = safe_float(row["ref_N50"])
        hyb_n50 = safe_float(row["hyb_N50"])
        short_n50 = safe_float(row["short_N50"])
        if hyb_n50 > 0 and short_n50 > 0 and ref_n50 > hyb_n50 and hyb_n50 > short_n50:
            rows.append([row[k] for k in summary[0].keys()] + [ref_n50 / hyb_n50, hyb_n50 / short_n50])

    rows.sort(key=lambda r: (float(r[-2]), float(r[-1])), reverse=True)
    write_tsv(out_path, header, rows)
    return {"output": str(out_path), "rows": len(rows)}


def build_summary_ranked_hq_bgcrich(
    summary_full_with_taxid: Path,
    target_taxids: Path,
    out_path: Path,
    max_ref_contigs: int = 5,
    min_hyb_total: int = 3_000_000,
    min_short_total: int = 2_000_000,
) -> Dict[str, object]:
    summary = read_tsv(summary_full_with_taxid)
    if not summary:
        raise ValueError("summary_full_with_taxid.tsv is empty")

    keep_taxids = read_target_taxids(target_taxids)
    header = list(summary[0].keys()) + ["ref/hyb", "hyb/short"]
    rows: List[List[object]] = []

    for row in summary:
        ref_contigs = safe_int(row["ref_contigs"])
        hyb_total = safe_int(row["hyb_total_bp"])
        short_total = safe_int(row["short_total_bp"])
        ref_n50 = safe_float(row["ref_N50"])
        hyb_n50 = safe_float(row["hyb_N50"])
        short_n50 = safe_float(row["short_N50"])
        taxid = row["TAXID"]

        if (
            taxid in keep_taxids and
            ref_contigs <= max_ref_contigs and
            hyb_total >= min_hyb_total and
            short_total >= min_short_total and
            hyb_n50 > 0 and short_n50 > 0 and
            ref_n50 > hyb_n50 and hyb_n50 > short_n50
        ):
            rows.append([row[k] for k in summary[0].keys()] + [ref_n50 / hyb_n50, hyb_n50 / short_n50])

    rows.sort(key=lambda r: (float(r[-2]), float(r[-1])), reverse=True)
    write_tsv(out_path, header, rows)
    return {"output": str(out_path), "rows": len(rows)}


def first_clusters_file(folder: Path) -> Optional[Path]:
    files = sorted(folder.glob("*.clusters.tsv"))
    return files[0] if files else None


def count_clusters(tsv_path: Path) -> int:
    with open(tsv_path, "r") as handle:
        n = sum(1 for _ in handle)
    return max(0, n - 1)


def summarize_gecco_dirs(gecco_root: Path, out_path: Optional[Path] = None):
    by_otu: Dict[str, Dict[str, Path]] = defaultdict(dict)

    for folder in gecco_root.iterdir():
        if not folder.is_dir():
            continue
        name = folder.name
        if name.startswith("gecco_run_on_") and name.endswith("_ref"):
            otu = name[len("gecco_run_on_"):-len("_ref")]
            by_otu[otu]["ref"] = folder
        elif name.startswith("gecco_run_on_") and name.endswith("_hybrid_MAG"):
            otu = name[len("gecco_run_on_"):-len("_hybrid_MAG")]
            by_otu[otu]["hyb"] = folder
        elif name.startswith("gecco_run_on_") and name.endswith("_short_MAG"):
            otu = name[len("gecco_run_on_"):-len("_short_MAG")]
            by_otu[otu]["short"] = folder

    rows = []
    for otu, parts in sorted(by_otu.items()):
        row: Dict[str, object] = {"OTU": otu}
        for label in ("ref", "hyb", "short"):
            folder = parts.get(label)
            if folder is None:
                row[f"{label}_clusters_file"] = "NA"
                row[f"{label}_bgcs"] = 0
                continue
            f = first_clusters_file(folder)
            row[f"{label}_clusters_file"] = str(f) if f else "NA"
            row[f"{label}_bgcs"] = count_clusters(f) if f else 0
        rows.append(row)

    if out_path:
        header = [
            "OTU",
            "ref_clusters_file", "ref_bgcs",
            "hyb_clusters_file", "hyb_bgcs",
            "short_clusters_file", "short_bgcs",
        ]
        tsv_rows = [[row.get(h, "") for h in header] for row in rows]
        write_tsv(out_path, header, tsv_rows)

    return rows


def make_compare_to_reference_commands(
    otu: str,
    reference_fasta: Path,
    ref_clusters: Path,
    hyb_dir: Path,
    short_dir: Path,
    hyb_gecco_dir: Path,
    short_gecco_dir: Path,
    quast_out_dir: Path,
    bgcquast_out_dir: Path,
) -> Dict[str, str]:
    safe = otu.replace(".", "_")
    hyb_safe_fasta = hyb_dir / f"hybrid_{safe}.fasta"
    short_safe_fasta = short_dir / f"short_{safe}.fasta"
    hyb_clusters = hyb_gecco_dir / f"hybrid_{safe}.clusters.tsv"
    short_clusters = short_gecco_dir / f"short_{safe}.clusters.tsv"

    quast_cmd = f'''quast.py \\
  {hyb_safe_fasta} \\
  {short_safe_fasta} \\
  -r {reference_fasta} \\
  -o {quast_out_dir} \\
  -t 16'''

    bgcquast_cmd = f'''python ../bgc-quast/bgc-quast.py \\
  {hyb_clusters} \\
  {short_clusters} \\
  --mode compare-to-reference \\
  --reference-mining-result {ref_clusters} \\
  --quast-output-dir {quast_out_dir} \\
  --output-dir {bgcquast_out_dir}'''

    return {
        "sanitized_name": safe,
        "quast_cmd": quast_cmd,
        "bgcquast_cmd": bgcquast_cmd,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plant CAMI II triplet workflow")
    p.add_argument("--hyb-binning", type=Path, required=True)
    p.add_argument("--short-binning", type=Path, required=True)
    p.add_argument("--hyb-gsa", type=Path, required=True)
    p.add_argument("--short-gsa", type=Path, required=True)
    p.add_argument("--genome-map", type=Path, required=True)
    p.add_argument("--target-groups", type=Path, required=True)
    p.add_argument("--target-taxids", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)

    p.add_argument("--generate-mags", action="store_true")
    p.add_argument("--build-summaries", action="store_true")
    p.add_argument("--summarize-gecco", action="store_true")
    p.add_argument("--print-commands-for-otu", type=str, default=None)
    p.add_argument("--overwrite-mags", action="store_true")
    p.add_argument("--report-every", type=int, default=200)
    p.add_argument("--max-ref-contigs", type=int, default=5)
    p.add_argument("--min-hyb-total", type=int, default=3_000_000)
    p.add_argument("--min-short-total", type=int, default=2_000_000)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    hyb_mag_out = outdir / "good_hybrid_assembly_based_MAGs_nanosim"
    short_mag_out = outdir / "poor_short_assembly_based_MAGs"

    summary_full = outdir / "summary_full.tsv"
    summary_full_with_taxid = outdir / "summary_full_with_taxid.tsv"
    otu_to_taxid_out = outdir / "otu_to_taxid.tsv"
    summary_ranked = outdir / "summary_ranked.tsv"
    summary_ranked_hq_bgcrich = outdir / "summary_ranked_hq_bgcrich.tsv"
    gecco_summary = outdir / "gecco_summary.tsv"

    run_summary: Dict[str, object] = {"outdir": str(outdir)}

    if args.generate_mags:
        print("[step] Generating hybrid MAGs")
        hyb_info = generate_mags_from_binning(
            binning_file=args.hyb_binning,
            assembly_fasta=args.hyb_gsa,
            output_dir=hyb_mag_out,
            overwrite=args.overwrite_mags,
            report_every=args.report_every,
        )
        print(json.dumps(hyb_info, indent=2))
        run_summary["hybrid_mag_generation"] = hyb_info

        print("[step] Generating short MAGs")
        short_info = generate_mags_from_binning(
            binning_file=args.short_binning,
            assembly_fasta=args.short_gsa,
            output_dir=short_mag_out,
            overwrite=args.overwrite_mags,
            report_every=args.report_every,
        )
        print(json.dumps(short_info, indent=2))
        run_summary["short_mag_generation"] = short_info

    if args.build_summaries:
        print("[step] Building summary_full.tsv")
        sf = build_summary_full(
            genome_map=args.genome_map,
            hyb_dir=hyb_mag_out,
            short_dir=short_mag_out,
            out_path=summary_full,
        )
        print(json.dumps(sf, indent=2))
        run_summary["summary_full"] = sf

        print("[step] Building otu_to_taxid.tsv and summary_full_with_taxid.tsv")
        otu_to_taxid = read_otu_to_taxid_from_target_groups(args.target_groups)
        sft = build_summary_full_with_taxid(
            summary_full=summary_full,
            otu_to_taxid=otu_to_taxid,
            otu_to_taxid_out=otu_to_taxid_out,
            out_path=summary_full_with_taxid,
        )
        print(json.dumps(sft, indent=2))
        run_summary["summary_full_with_taxid"] = sft

        print("[step] Building summary_ranked.tsv")
        sr = build_summary_ranked(summary_full_with_taxid=summary_full_with_taxid, out_path=summary_ranked)
        print(json.dumps(sr, indent=2))
        run_summary["summary_ranked"] = sr

        print("[step] Building summary_ranked_hq_bgcrich.tsv")
        srh = build_summary_ranked_hq_bgcrich(
            summary_full_with_taxid=summary_full_with_taxid,
            target_taxids=args.target_taxids,
            out_path=summary_ranked_hq_bgcrich,
            max_ref_contigs=args.max_ref_contigs,
            min_hyb_total=args.min_hyb_total,
            min_short_total=args.min_short_total,
        )
        print(json.dumps(srh, indent=2))
        run_summary["summary_ranked_hq_bgcrich"] = srh

    if args.summarize_gecco:
        print("[step] Summarizing GECCO directories")
        gecco_rows = summarize_gecco_dirs(outdir, gecco_summary)
        print(f"GECCO OTUs summarized: {len(gecco_rows)}")
        run_summary["gecco_summary"] = {"output": str(gecco_summary), "rows": len(gecco_rows)}

    if args.print_commands_for_otu:
        otu = args.print_commands_for_otu
        if not summary_full_with_taxid.exists():
            eprint("summary_full_with_taxid.tsv does not exist yet.")
        else:
            summary = read_tsv(summary_full_with_taxid)
            selected = next((row for row in summary if row["OTU"] == otu), None)
            if not selected:
                eprint(f"OTU not found: {otu}")
            else:
                otu_safe = otu.replace(".", "_")
                cmds = make_compare_to_reference_commands(
                    otu=otu,
                    reference_fasta=Path(selected["reference_genome"]),
                    ref_clusters=outdir / f"gecco_run_on_{otu}_ref" / f"{Path(selected['reference_genome']).stem}.clusters.tsv",
                    hyb_dir=hyb_mag_out,
                    short_dir=short_mag_out,
                    hyb_gecco_dir=outdir / f"gecco_run_on_{otu}_hybrid_MAG",
                    short_gecco_dir=outdir / f"gecco_run_on_{otu}_short_MAG",
                    quast_out_dir=outdir / f"assembly_to_ref_align_quast_{otu_safe}",
                    bgcquast_out_dir=outdir / f"bgcquast_{otu_safe}_compare_to_reference",
                )
                print("\n[QUAST command]\n")
                print(cmds["quast_cmd"])
                print("\n[BGC-QUAST command]\n")
                print(cmds["bgcquast_cmd"])
                run_summary["printed_commands_for_otu"] = otu

    run_summary_path = outdir / "run_summary.json"
    with open(run_summary_path, "w") as handle:
        json.dump(run_summary, handle, indent=2)
    print(f"\n[done] Wrote run summary to {run_summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
