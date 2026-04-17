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
  --hyb-binning /net/projects/halmessa/BGC_project/genome_binning_cami2/plant_associated_dataset/ground_truth/rhizosphere_hybrid_pooled_gsa.binning \
  --short-binning /net/projects/halmessa/BGC_project/genome_binning_cami2/plant_associated_dataset/ground_truth/rhizosphere_short_read_pooled_gsa.binning \
  --hyb-gsa /net/projects/halmessa/BGC_project/plant_BGC_paper/rhimgCAMI2_hybrid_nanosim_pooled_gsa.fasta \
  --short-gsa /net/projects/halmessa/BGC_project/plant_BGC_paper/rhimgCAMI2_short_read_pooled_gsa.fasta \
  --genome-map /net/projects/halmessa/BGC_project/plant_BGC_paper/genome_to_id.tsv \
  --target-groups /net/projects/halmessa/BGC_project/plant_BGC_paper/target_groups.tsv \
  --target-taxids /net/projects/halmessa/BGC_project/plant_BGC_paper/target_taxids.txt \
  --outdir /net/projects/halmessa/BGC_project/plant_BGC_paper/run_from_script \
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
run_quast_bgcquast_top5.sh
```

## Author

Hesham Almessady
