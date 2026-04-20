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




# Triplet Selection Workflow for BGC-QUAST from CAMI2 plant dataset 

## Step 1 — get the ground truth hybrid-reads binning file for CAMI2 plant dataset 

```bash
cp /net/sgi/cami/data/CAMI2/second_challenge_evaluation/binning/genome_binning/plant_associated_dataset/data/ground_truth/rhizosphere_hybrid_pooled_gsa.binning.tar.gz .
tar -xzvf rhizosphere_hybrid_pooled_gsa.binning.tar.gz
```

## Step 2 — get the ground truth short-reads binning file for CAMI2 plant dataset 

```bash
cp /net/sgi/cami/data/CAMI2/second_challenge_evaluation/binning/genome_binning/plant_associated_dataset/data/ground_truth/rhizosphere_short_read_pooled_gsa.binning.tar.gz .
tar -xzvf rhizosphere_short_read_pooled_gsa.binning.tar.gz
``` 

## Step 3 — get target taxids for BGCs rich taxa (Actinomycetes, Cyanobacteria, Myxococcota)

```bash
awk -F'\t' '!/^@/ {print $3}' \
rhizosphere_hybrid_pooled_gsa.binning | \
sort -u | \
taxonkit lineage | \
taxonkit reformat | \
grep -Ei 'Actinomycet|Actinobacter|Cyanobacteri|Myxococc' | \
cut -f1 > target_taxids.txt
```

## Step 4 — filter the binning file using those taxids

```bash
awk -F'\t' 'NR==FNR {keep[$1]=1; next} /^@/ || ($3 in keep)' target_taxids.txt \
rhizosphere_hybrid_pooled_gsa.binning > target_groups.tsv
```

## Step 5 — get genome mapping file from the simulation directory

```bash
cp /net/sgi/cami/data/CAMI2/roter_data/lotus_rhizosphere/simulation_short_read/genome_to_id.tsv . 
```

## Step 6 — get short-reads gold standard assembly fasta file 

```bash
cp /net/sgi/cami/data/CAMI2/roter_data/lotus_rhizosphere/rhimgCAMI2/assembly/rhimgCAMI2_short_read_pooled_gsa.fasta.gz . 
gunzip rhimgCAMI2_short_read_pooled_gsa.fasta.gz 
```

## Step 7 — get hybrid-reads gold standard assembly fasta file 

```bash
cp /net/benchmarking/fmeyer/cami2/rhizosphere/assembly/rhimgCAMI2_hybrid_nanosim_pooled_gsa.fasta.gz . 
gunzip rhimgCAMI2_hybrid_nanosim_pooled_gsa.fasta.gz
``` 

## Step 8 — Run the triplet workflow python script 

```bash
python triplet_workflow.py \
  --hyb-binning rhizosphere_hybrid_pooled_gsa.binning  \
  --short-binning rhizosphere_short_read_pooled_gsa.binning  \
  --hyb-gsa rhimgCAMI2_hybrid_nanosim_pooled_gsa.fasta  \
  --short-gsa rhimgCAMI2_short_read_pooled_gsa.fasta  \
  --genome-map genome_to_id.tsv \
  --target-groups target_groups.tsv \
  --target-taxids target_taxids.txt \
  --outdir path to output directory \
  --generate-mags \
  --build-summaries
```

