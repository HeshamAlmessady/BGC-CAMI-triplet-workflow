# CAMI Triplet Workflow

A reproducible workflow for generating and ranking reference / hybrid MAG / short-read MAG triplets from CAMI benchmark datasets.

## What it does

This workflow:

- generates hybrid and short MAG FASTA files from CAMI ground-truth `.binning` files and GSA FASTA files
- computes assembly statistics:
  - contig count
  - total bp
  - N50
- builds ranked candidate triplet tables
- helps identify strong candidates for downstream GECCO, QUAST, and BGC-QUAST analysis

## Main script

```bash
triplet_workflow.py
```

## Example command

```bash
python triplet_workflow.py \
  --hyb-binning path to CAMI ground truth hybrid-reads binning file  \
  --short-binning path to CAMI ground truth short-reads binning file  \
  --hyb-gsa path to CAMI hybrid-reads gsa fasta file  \
  --short-gsa path to CAMI short-reads gsa fasta file  \
  --genome-map  path to genome_to_id.tsv \
  --target-groups path to target_groups.tsv \
  --target-taxids path to target_taxids.txt \
  --outdir path to output directory \
  --generate-mags \
  --build-summaries
```

## Outputs

The script generates files such as:

- `summary_full.tsv`
- `summary_full_with_taxid.tsv`
- `summary_ranked.tsv`
- `summary_ranked_hq_bgcrich.tsv`

and MAG directories such as:

- `good_hybrid_assembly_based_MAGs_nanosim/`
- `poor_short_assembly_based_MAGs/`

## Notes

The downstream GECCO, QUAST, and BGC-QUAST analysis can be run separately using a shell script, for example:

```bash
run_quast_bgcquast_top5.sh path to triplet workflow output directory
```

## Author

Hesham Almessady
