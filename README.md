# TPanno

**A standardized and scalable workflow for comprehensive multi-omics analysis in TickPathDB**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Manuscript title:** TickPathDB: A Comprehensive Multi-Omics Platform for Ticks and Tick-Associated Pathogens
>
> **Author:** Yangfan Wu
>
> **Affiliation:** School of Clinical Medicine, Wannan Medical University, Wuhu, China
>
> **Email:** wuyf2025wmc@outlook.com
>
> **Laboratory:** http://www.wnmc-bioinfo.com/team/

---

## Overview

TPanno is the modular, command-line analysis workflow developed for **TickPathDB**, a multi-omics resource for ticks and tick-associated pathogens. It covers the complete path from raw genome assemblies and RNA-seq reads to publication-ready annotation and reproducible result tables:

1. **Genome quality assessment** with BUSCO/compleasm and CheckM2.
2. **Repeat annotation and genome masking** with EDTA and HiTE.
3. **Structural annotation** with EGAPx (eukaryotes), Bakta (prokaryotes), and geNomad (viruses).
4. **Functional annotation** with eggNOG-mapper, PANNZER, and Swiss-Prot BLASTP.
5. **Transcriptomic integration** with the SRA Toolkit and Trinity.
6. **Optional domain-based phylogeny and protein-family clustering** with the bundled [`domain_phylo_pipeline/`](domain_phylo_pipeline/README.md) module.

Every step is independent and can be run on its own. All commands are written explicitly so that the exact parameters used in a study can be recorded, reviewed, and reproduced.

### Workflow at a glance

| Step | Purpose | Main tools | Typical output |
|---|---|---|---|
| 1 | Genome quality assessment | BUSCO, compleasm, CheckM2 | Completeness and contamination reports |
| 2 | Repeat annotation and masking | EDTA, HiTE | TE annotation, soft-masked genome |
| 3 | Structural annotation | EGAPx, Bakta, geNomad | GFF3/GTF, CDS, proteins |
| 4 | Functional annotation | eggNOG-mapper, PANNZER, Swiss-Prot BLASTP | GO/KEGG/COG assignments, putative homolog descriptions |
| 5 | Transcriptomic integration | SRA Toolkit, Trinity | Assembled transcripts |
| 6 | Domain phylogeny (optional) | HMMER, MAFFT, trimAl, IQ-TREE | Domain trees, clusters, QC report |

---

## Repository layout

```text
TPanno/
├── README.md                     # This file
├── LICENSE                       # MIT License
├── VERSION                       # Release version
├── scripts/
│   ├── gca_fna.sh                # Batch-download genome FASTA by GCA accession
│   ├── gca_faa.sh                # Batch-download protein FASTA by GCA accession
│   ├── setup_local_cache.sh      # Download the EGAPx local cache with rsync
│   └── add_swissprot_annotation.py # Attach Swiss-Prot BLASTP hits to FAA headers
└── domain_phylo_pipeline/        # Domain-based phylogeny and clustering module
    ├── README.md
    ├── run_pipeline.sh
    ├── run_pfam_scan.sh
    ├── config/
    ├── bin/
    ├── docs/
    ├── env/
    └── test/
```

---

## 1. Prerequisites

### 1.1 Software

**Genome assessment**

- [BUSCO](https://busco.ezlab.org/) (v5.8.3 or later). The workflow also supports [compleasm](https://github.com/huangnengCSU/compleasm), a faster and more accurate reimplementation of BUSCO.
- [CheckM2](https://github.com/chklovski/checkm2)

**Repeat annotation and genome masking**

- [EDTA](https://github.com/oushujun/EDTA)
- [HiTE](https://github.com/CSU-KangHu/HiTE)

**Gene prediction and structural annotation**

- [EGAPx](https://github.com/ncbi/egapx) for eukaryotic genomes
- [Bakta](https://github.com/oschwengers/bakta) for prokaryotic genomes
- [geNomad](https://github.com/apcamargo/genomad) for viral genomes

**Functional annotation**

- [eggNOG-mapper](http://eggnog-mapper.embl.de/)
- [PANNZER](http://ekhidna2.biocenter.helsinki.fi/sanspanz/)

**Transcriptomic integration**

- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)

**Domain phylogeny (optional, bundled module)**

- [HMMER](http://hmmer.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [trimAl](http://trimal.cgenomics.org/)
- [IQ-TREE](http://www.iqtree.org/)
- [CD-HIT](https://github.com/weizhongli/cdhit)

**Utilities**

- [DIAMOND](https://github.com/bbuchfink/diamond)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- `wget`, `unzip`, and `rsync` for data download and cache synchronization
- Miniconda or Mamba for environment management (recommended)

**Tool note.** DIAMOND provides fast protein alignment for functional annotation (Buchfink et al., 2015), and BLAST+ provides the NCBI alignment suite (Camacho et al., 2009).

Install all tools in dedicated conda environments, as shown in the sections below. Record the exact version of every tool in the methods section and, where possible, use the version-pinned environments provided by the tool developers.

### 1.2 Databases

| Database | Used by | Source |
|---|---|---|
| `arachnida_odb12` | BUSCO / compleasm | https://busco-data.ezlab.org/v5/data/lineages/arachnida_odb12.2025-02-25.tar.gz |
| EGAPx `local_cache` | EGAPx | https://github.com/ncbi/egapx |
| `db-light` | Bakta | https://zenodo.org/record/7669534/files/db-light.tar.gz?download=1 |
| `genomad_db_v1.9` | geNomad | https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz?download=1 |
| `CheckM2_database` | CheckM2 | https://zenodo.org/api/files/fd3bc532-cd84-4907-b078-2e05a1e46803/checkm2_database.tar.gz |
| UniProt Swiss-Prot | Section 5.3 | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz |
| Pfam-A HMM library | `domain_phylo_pipeline` | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/ |
| `Pfam-A.clans.tsv` (optional) | `domain_phylo_pipeline` | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz |

### 1.3 Datasets

- Genome assemblies: [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/)
- RNA-seq reads: [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/)

The accession numbers analysed in the manuscript are listed in its supplementary tables. The scripts below retrieve the corresponding files directly from NCBI.

Paths such as `./tickdb/genome.fna` are placeholders; replace them with the actual locations on your system.

### 1.4 Batch download of genomes and proteins

The download scripts take a tab-separated sheet with two columns: a sample prefix and a GCA accession.

`assemblies_sheet.tsv`:

```tsv
prefix	assembly
sample01	GCA_000000000.0
sample02	GCA_000000001.0
```

Replace the placeholder accessions with real GCA accessions from NCBI.

- `prefix` is the sample name. It is also used as the prefix for output files and, during annotation, as the locus-tag prefix. Keep it at most 24 characters and use only letters, digits, underscores, dots, and hyphens.
- `assembly` is the GCA accession of the GenBank genome to download.

A header row such as `prefix<TAB>assembly` is optional and is skipped.

Download the genome (`.fna`) and protein (`.faa`) files:

```bash
bash scripts/gca_fna.sh assemblies_sheet.tsv   # writes genomes/<prefix>.fna
bash scripts/gca_faa.sh assemblies_sheet.tsv   # writes proteins/<prefix>.faa
```

Both scripts skip files that already exist, so they can be re-run safely after an interruption. They require `wget`, `unzip`, and network access to the NCBI Datasets API.

Genome and protein records are retrieved from NCBI Datasets (NCBI, 2026a).

---

## 2. Genome assessment

### 2.1 BUSCO

Inputs:

- `genome.fna`
- `arachnida_odb12` lineage directory

Create the environment:

```bash
conda create -n busco -c conda-forge -c bioconda busco -y
conda activate busco
```

Run BUSCO in genome mode:

```bash
busco -i ./tickdb/genome.fna \
      -l ./arachnida_odb12 \
      -o ./busco_output \
      -m geno \
      --cpu 32 \
      --offline
```

Inspect the short summary:

```bash
cat ./busco_output/short_summary.specific.arachnida_odb12.out.txt
```

**Tool note.** BUSCO estimates genome completeness by searching for lineage-specific single-copy orthologs (Simão et al., 2015); the v5 lineage sets and workflow are described by Manni et al. (2021).

### 2.2 compleasm

compleasm is a faster and more accurate reimplementation of BUSCO. The organization of its lineage files differs from that of BUSCO, so download the lineage in the compleasm format before running:

```bash
python compleasm.py run \
    -t 16 \
    -l arachnida \
    -L ./arachnida_odb12 \
    -a genome.fna \
    -o compleasm_output
```

**Tool note.** compleasm reimplements the BUSCO workflow with a protein-to-genome alignment step, which improves both speed and completeness estimates on large eukaryotic genomes (Huang & Li, 2023).

### 2.3 CheckM2

Inputs:

- `bacteriadb/` (one FASTA file per genome)
- `CheckM2_database/`

CheckM2 was designed to estimate the completeness and contamination of metagenome-assembled genomes (MAGs) and single-amplified genomes (SAGs), but it also works for isolate genomes (Chklovski et al., 2023).

```bash
checkm2 predict \
    --threads 16 \
    --input ./bacteriadb \
    --database_path ./CheckM2_database/uniref100.KO.1.dmnd \
    --output_directory ./checkm2_output
```

By default, the output directory contains `quality_report.tsv`, a tab-delimited table with the completeness and contamination of each genome. Add `--stdout` to print the table to standard output instead.

---

## 3. Repeat annotation and genome masking

EDTA and HiTE are complementary. EDTA builds a consensus TE library and produces a soft-masked genome for gene prediction; HiTE adds structure-aware TE discovery. At minimum, run EDTA and soft-mask the genome before structural annotation.

### 3.1 EDTA

Input: `genome.fna`

Run the following commands from the EDTA source directory (the directory containing `EDTA.pl`):

```bash
perl EDTA.pl \
    --genome ./tickdb/genome.fna \
    --overwrite 1 \
    -t 12 \
    --sensitive 1 \
    --anno 1 \
    --force 1 \
    --evaluate 1
```

Convert hard masking to soft masking:

```bash
perl ./util/make_masked.pl \
    -genome ./tickdb/genome.fna \
    -minlen 80 \
    -hardmask 0 \
    -t 10 \
    -rmout genome.fna.mod.EDTA.TEanno.out
```

The soft-masked genome is used as the input for EGAPx in Section 4.1.

**Tool note.** EDTA combines multiple TE detection methods and builds a non-redundant consensus library for whole-genome annotation and masking (Ou et al., 2019).

### 3.2 HiTE

Input: `genome.fna`

HiTE annotates the genome with its own TE library and writes `HiTE.out`, `HiTE.gff`, and `HiTE.tbl` when annotation is enabled (Hu et al., 2024). Run the following command from the HiTE source directory (the directory containing `main.py`):

```bash
python main.py \
    --genome ./tickdb/genome.fna \
    --outdir ../hite_output \
    --thread 12 \
    --annotate 1 \
    --plant 0
```

HiTE keeps many intermediate files so that an interrupted run can be resumed with `--recover 1`. A successful run produces a directory like:

```text
hite_output/
├── longest_repeats_*.fa
├── longest_repeats_*.flanked.fa
├── confident_tir_*.fa
├── confident_helitron_*.fa
├── confident_non_ltr_*.fa
├── confident_other_*.fa
├── confident_ltr_cut.fa.cons
├── confident_TE.cons.fa
├── HiTE.out      # requires --annotate 1
├── HiTE.gff      # requires --annotate 1
└── HiTE.tbl      # requires --annotate 1
```

---

## 4. Gene prediction and structural annotation

### 4.1 Eukaryotic structural annotation with EGAPx

Inputs:

- `genome.fna.masked` (the soft-masked genome from Section 3.1)
- `input.yaml`
- `egapx_0.4.0-alpha.sif`
- EGAPx `local_cache`
- NCBI TaxID of the target organism

> `taxid` is the NCBI Taxonomy identifier. If the species-level TaxID is unknown, a TaxID for a higher taxonomic level can be used. Look up the TaxID at https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi.

**Tool note.** EGAPx is NCBI's eukaryotic genome annotation pipeline (NCBI, 2026b). It runs as a Nextflow workflow (Di Tommaso et al., 2017) inside Singularity containers (Kurtzer et al., 2017), which is why both tools are installed in the environment below.

#### 4.1.1 Create the environment

```bash
conda create -n egapx -c conda-forge -c bioconda python=3.11 pyyaml nextflow singularity -y
conda activate egapx
```

Pull the EGAPx container:

```bash
singularity pull docker://ncbi/egapx:0.4.0-alpha
```

If the default registry is slow or unavailable, the same image can be pulled from a public mirror and tagged with the expected filename:

```bash
singularity pull --name egapx_0.4.0-alpha.sif docker://docker.1ms.run/ncbi/egapx:0.4.0-alpha
```

Clone the EGAPx repository:

```bash
git clone https://github.com/ncbi/egapx.git
cd egapx
```

#### 4.1.2 Populate the local cache

Run the commands in this subsection from the cloned `egapx/` directory.

Use the built-in downloader:

```bash
python3 ui/egapx.py -dl -lc ../local_cache
```

Alternatively, create the cache manually with rsync:

```bash
mkdir -p ../local_cache
cd ../local_cache
mkdir -p gnomon ortholog_references target_proteins taxonomy \
         reference_sets misc cmsearch busco_downloads/lineages

rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/gnomon/2 gnomon/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/ortholog_references/3 ortholog_references/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/target_proteins/2 target_proteins/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/taxonomy/2 taxonomy/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/reference_sets/3 reference_sets/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/misc/2 misc/
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/cmsearch/1 cmsearch/

cd ../egapx
```

The same steps are wrapped in a helper script. Run it from the repository root (not the `egapx/` directory); it creates the cache at the repository root, so `-lc ../local_cache` remains correct after you return to `egapx/`:

```bash
chmod +x scripts/setup_local_cache.sh
bash scripts/setup_local_cache.sh
```

#### 4.1.3 Write the input YAML

EGAPx reads a YAML file with the following required keys:

```yaml
genome: /path/to/genome.fna.masked
taxid: 6945
reads:
  - /path/to/reads_1.fastq.gz
  - /path/to/reads_2.fastq.gz
annotation_provider: TickPathDB
annotation_name_prefix: TickPathDB
locus_tag_prefix: tick
```

- `genome`: assembled genome in FASTA format (soft-masked is recommended).
- `taxid`: NCBI Taxonomy identifier of the target organism.
- `reads`: RNA-seq reads used as transcript evidence. Paired-end reads can be given as separate files; the pipeline also accepts SRA accessions.
- `annotation_provider`: name recorded in the GFF3 annotation.
- `annotation_name_prefix`: prefix used for annotation names.
- `locus_tag_prefix`: prefix used for locus tags. Keep it short and unique within the project.

#### 4.1.4 Run EGAPx

```bash
python3 ui/egapx.py input.yaml \
    -e singularity \
    -w anno \
    -o egapx_output \
    -lc ../local_cache
```

If the first run stops because no container is configured, point EGAPx to the SIF image and rerun:

```bash
echo "process.container = '$(pwd)/egapx_0.4.0-alpha.sif'" >> egapx_config/singularity.config
python3 ui/egapx.py input.yaml \
    -e singularity \
    -w anno \
    -o egapx_output \
    -lc ../local_cache
```

#### 4.1.5 Outputs

The main annotation file is `complete.genomic.gff` in the output directory:

```text
egapx_output/
├── annot_builder_output/
├── annotated_genome.asn
├── annotation_data.cmt
├── complete.cds.fna
├── complete.genomic.fna
├── complete.genomic.gff
├── complete.genomic.gtf
├── complete.proteins.faa
├── complete.transcripts.fna
├── nextflow.log
├── resume.sh
├── run.report.html
├── run.timeline.html
├── run.trace.txt
├── run_params.yaml
├── stats/
└── validated/
```

- `complete.genomic.gff`: final annotation in GFF3 format.
- `complete.genomic.gtf`: final annotation in GTF format.
- `complete.genomic.fna`: genome sequences in FASTA format.
- `complete.cds.fna`: annotated coding sequences (CDS) in FASTA format.
- `complete.transcripts.fna`: annotated transcripts in FASTA format, including UTRs.
- `complete.proteins.faa`: annotated protein products in FASTA format.
- `annotated_genome.asn`: final annotation in ASN.1 format.
- `run_params.yaml`: exact parameters used for the run, which should be archived with the results.

### 4.2 Prokaryotic structural annotation with Bakta

Inputs:

- `genome.fna`
- `db-light`

Optionally update the AMRFinderPlus database shipped with Bakta:

```bash
mkdir -p db-light/amrfinderplus-db
amrfinder_update --force_update --database db-light/amrfinderplus-db/
```

Run Bakta:

```bash
bakta \
    --db ./db-light \
    --threads 16 \
    --output ./bakta_output \
    ./bacteriadb/genome.fna
```

Bakta produces dbxref-rich, sORF-aware, taxon-independent annotations in JSON and standard bioinformatics formats (GFF3, GenBank, EMBL, TSV), which can be used directly in downstream analyses.

**Tool note.** Bakta performs alignment-free, taxon-independent bacterial genome annotation and is designed for rapid, standardized output (Schwengers et al., 2021).

### 4.3 Viral genome annotation with geNomad

Inputs:

- `genome.fna.gz`
- `genomad_db`

geNomad depends on a database of marker profiles, taxonomic information, and functional annotations. Download it into the current directory first:

```bash
genomad download-database .
```

The database is stored in `genomad_db/`. geNomad can run its modules one by one, but the end-to-end command runs the whole pipeline in a single call:

```bash
genomad end-to-end \
    --cleanup \
    --splits 8 \
    ./virusdb/genome.fna.gz \
    ./genomad_output \
    ./genomad_db
```

Results are written to `genomad_output/`. Use `--splits` to control the number of parallel chunks; lower it if memory is limited.

**Tool note.** geNomad combines gene content, marker profiles, and a deep neural network to identify plasmids and viruses and to annotate their genes (Camargo et al., 2024).

---

## 5. Functional annotation

### 5.1 eggNOG-mapper

Input: `pep.faa`

The web server is sufficient for small protein sets:

1. Open http://eggnog-mapper.embl.de/.
2. Upload `pep.faa`.
3. Select the appropriate taxonomic scope and annotation options.
4. Download the `.emapper.annotations` table.

For local, high-throughput runs, install eggNOG-mapper and its database, then run:

```bash
download_eggnog_data.py -y
emapper.py -i pep.faa --output eggnog -m diamond --cpu 16
```

The main output is `eggnog.emapper.annotations`, which contains GO terms, KEGG orthologs, COG categories, and free-text descriptions for each protein.

**Tool note.** eggNOG-mapper v2 assigns orthology and functional annotations from precomputed eggNOG orthologous groups, using DIAMOND for the sequence search (Cantalapiedra et al., 2021; Buchfink et al., 2015).

### 5.2 PANNZER

Input: `pep.faa`

1. Open http://ekhidna2.biocenter.helsinki.fi/sanspanz/.
2. Upload `pep.faa`.
3. Choose the appropriate annotation mode (PANNZER2 is recommended for general use).
4. Download the result table.

PANNZER returns GO terms, enzyme commission numbers, and descriptive annotations. Because it is a web service, record the submission date and the database version reported in the result file.

**Tool note.** PANNZER2 is a rapid functional annotation web server that predicts GO terms and enzyme commission numbers from protein sequence (Törönen et al., 2018).

### 5.3 Swiss-Prot homology-based annotation

Input: `protein.faa`

This step transfers the description of the best BLASTP hit (Camacho et al., 2009) against the manually reviewed Swiss-Prot section of UniProt (UniProt Consortium, 2025). It provides a useful **putative** annotation for proteins that are absent from, or poorly represented in, the other annotation sources.

Download the current Swiss-Prot FASTA and build a BLAST database:

```bash
mkdir -p db/swiss_prot
wget -P db/swiss_prot \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -f db/swiss_prot/uniprot_sprot.fasta.gz

makeblastdb \
  -in db/swiss_prot/uniprot_sprot.fasta \
  -dbtype prot \
  -title swiss_prot \
  -parse_seqids \
  -out db/swiss_prot/uniprot_sprot
```

If your BLAST+ build rejects `-parse_seqids` for the Swiss-Prot identifier format, omit that option; the search still works with the default database indexing.

Search the protein FASTA:

```bash
blastp \
  -query protein.faa \
  -db db/swiss_prot/uniprot_sprot \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -max_hsps 1 \
  -num_threads 50 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -out protein.blastp.tsv
```

Attach the best hit to each FASTA header and write a tabular summary:

```bash
python3 scripts/add_swissprot_annotation.py \
  --query protein.faa \
  --blast protein.blastp.tsv \
  --out protein.swissprot.faa \
  --summary protein.swissprot.tsv \
  --mark-missing
```

The annotated headers look like:

```text
>protein_1 swissprot=sp|P12345|EXAMPLE pident=62.5 length=200 evalue=1e-45 bitscore=180 Example protein
```

The helper keeps the lowest-E-value hit per query, breaks ties by bit score and percent identity, and can mark proteins without a retained hit as `swissprot=no_hit`. The summary table contains one row per query. Use `--min-pident` (for example, `--min-pident 30`) to discard low-identity transfers that pass the E-value threshold.

> **Interpretation caution:** a Swiss-Prot BLASTP hit is a homology-based annotation, not an orthology assignment. Check percent identity, alignment length, E-value, and the description before transferring a function; low-identity hits can be misleading.

---

## 6. Transcriptomic integration

### 6.1 SRA Toolkit

Input: SRA accession (for example, `SRRxxxxxxx`).

```bash
prefetch SRRxxxxxxx
fasterq-dump --split-3 --threads 8 SRRxxxxxxx
```

`--split-3` writes `SRRxxxxxxx_1.fastq` and `SRRxxxxxxx_2.fastq` for paired-end runs; for single-end runs it writes a single FASTQ file. Compress large files with `gzip` if disk space is limited.

**Tool note.** The SRA Toolkit retrieves runs from the NCBI Sequence Read Archive (Leinonen et al., 2011; NCBI, 2026c).

### 6.2 Trinity

Inputs:

- `SRRxxxxxxx_1.fastq`
- `SRRxxxxxxx_2.fastq`

Create the environment:

```bash
conda create -n trinity -c conda-forge -c bioconda python=3.8 trinity -y
conda activate trinity
```

Run Trinity:

```bash
Trinity \
    --seqType fq \
    --max_memory 200G \
    --left SRRxxxxxxx_1.fastq \
    --right SRRxxxxxxx_2.fastq \
    --CPU 28 \
    --trimmomatic \
    --full_cleanup
```

Two files are written to the working directory:

- `trinity_out_dir.Trinity.fasta`: the assembled transcripts.
- `trinity_out_dir.Trinity.fasta.gene_trans_map`: the transcript-to-gene map.

The FASTA file is the reference transcriptome used in downstream analyses.

**Tool note.** Trinity performs de novo transcriptome assembly from RNA-seq reads without a reference genome (Grabherr et al., 2011); its `--trimmomatic` option invokes Trimmomatic for read trimming (Bolger et al., 2014).

---

## 7. Domain-based phylogeny module

The [`domain_phylo_pipeline/`](domain_phylo_pipeline/README.md) directory contains an independent, reproducible pipeline for domain-based phylogeny and protein-family clustering. It takes one protein FASTA per sample (27 proteomes in the TickPathDB use case), identifies homologous domains with Pfam/HMMER, extracts only the HMM-defined domain region, and then performs MAFFT alignment, trimAl trimming, IQ-TREE maximum-likelihood inference, and patristic-distance clustering. It outputs per-family trees, subfamily assignments, support values, and a QC report.

**Tool note.** Domain identification uses HMMER (Eddy, 2011) with Pfam profiles (Mistry et al., 2021); alignment and trimming use MAFFT (Katoh & Standley, 2013) and trimAl (Capella-Gutiérrez et al., 2009); tree inference uses IQ-TREE 2 (Minh et al., 2020) with ModelFinder (Kalyaanamoorthy et al., 2017) and UFBoot2 (Hoang et al., 2018); dereplication uses CD-HIT (Fu et al., 2012), and the clustering step follows the cluster-diameter concept of TreeCluster (Balaban et al., 2019).

Quick start:

```bash
cd domain_phylo_pipeline
mamba env create -f env/environment.yml
conda activate domain-phylo
bash env/check_deps.sh

# Place one protein FASTA per sample in 00_raw/, then:
bash run_pipeline.sh --init
bash run_pipeline.sh --jobs 8 --threads 4
```

For a full Pfam-A scan with two-stage family screening:

```bash
cp /path/to/Pfam-A.hmm 01_hmms/
hmmpress -f 01_hmms/Pfam-A.hmm
MIN_SAMPLES=22 MIN_SEQS=22 MAX_FAMILIES=300 bash run_pfam_scan.sh
```

See [`domain_phylo_pipeline/README.md`](domain_phylo_pipeline/README.md) for the complete documentation, [`docs/methods.md`](domain_phylo_pipeline/docs/methods.md) for manuscript-ready methods text, and [`docs/interpretation.md`](domain_phylo_pipeline/docs/interpretation.md) for guidance on interpreting the results.

---

## 8. Reproducibility

The workflow is designed so that every analysis can be repeated from the archived code, the recorded parameters, and public data:

**Tool note.** The Python components run on Python 3.11 (Python Software Foundation, 2026) and use Biopython for sequence and tree handling (Cock et al., 2009), NumPy for numerical arrays (Harris et al., 2020), and SciPy for hierarchical clustering and distance calculations (Virtanen et al., 2020).

- **Explicit commands.** No hidden defaults: every tool invocation is written out in this README or in the step scripts.
- **Centralized parameters.** The domain pipeline keeps all thresholds in [`config/config.sh`](domain_phylo_pipeline/config/config.sh).
- **Recorded software versions.** Run `tool --version` for every tool and archive the output with the results. The domain pipeline writes `results/06_report/software_versions.txt` automatically.
- **Fixed random seeds.** IQ-TREE is run with `-seed 12345`; record any seed used by other stochastic tools.
- **Parameter snapshots.** EGAPx writes `run_params.yaml`; the domain pipeline writes `results/06_report/parameters.tsv`.
- **Checkpointing.** Long steps can be resumed: download scripts skip existing files, and the domain pipeline resumes completed families unless `--force` is given.
- **Versioned releases.** Tag each release that is used for an analysis, so that the exact code behind a published result can be retrieved later.

When reporting results, include the exact version numbers, the masking and annotation strategy, the database releases, and the random seeds. Small numerical differences between CPU counts or software builds are expected; fixing the thread count is recommended when bit-wise reproducibility is required.

---

## 9. Citation

If you use TPanno in your work, please cite the software release and the accompanying manuscript:

```text
Wu, Y. (2026). TPanno: A standardized and scalable workflow for comprehensive
multi-omics analysis in TickPathDB (Version 1.0.0) [Computer software].
https://github.com/wuyangfan641/TPanno
```


---

## 10. License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for the full text.

---

## 11. Contact

- **Author:** Yangfan Wu
- **Email:** wuyf2025wmc@outlook.com
- **Laboratory:** http://www.wnmc-bioinfo.com/team/
- **Issues:** please open an issue on the GitHub repository for bug reports and questions about the workflow.

---

## 12. References

References are formatted in a generic author-year style; adapt punctuation and journal abbreviations to the target journal's requirements.

1. Balaban M, Moshiri N, Mai U, Jia X, Mirarab S. TreeCluster: Clustering biological sequences using phylogenetic trees. *PLoS One*. 2019;14(8):e0221068. https://doi.org/10.1371/journal.pone.0221068
2. Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*. 2014;30(15):2114-2120. https://doi.org/10.1093/bioinformatics/btu170
3. Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. *Nat Methods*. 2015;12(1):59-60. https://doi.org/10.1038/nmeth.3176
4. Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421. https://doi.org/10.1186/1471-2105-10-421
5. Camargo AP, Roux S, Schulz F, et al. Identification of mobile genetic elements with geNomad. *Nat Biotechnol*. 2024;42(8):1303-1312. https://doi.org/10.1038/s41587-023-01953-y
6. Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. *Mol Biol Evol*. 2021;38(12):5825-5829. https://doi.org/10.1093/molbev/msab293
7. Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*. 2009;25(15):1972-1973. https://doi.org/10.1093/bioinformatics/btp348
8. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods*. 2023;20(8):1203-1212. https://doi.org/10.1038/s41592-023-01940-w
9. Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422-1423. https://doi.org/10.1093/bioinformatics/btp163
10. Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. *Nat Biotechnol*. 2017;35(4):316-319. https://doi.org/10.1038/nbt.3820
11. Eddy SR. Accelerated Profile HMM Searches. *PLoS Comput Biol*. 2011;7(10):e1002195. https://doi.org/10.1371/journal.pcbi.1002195
12. Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics*. 2012;28(23):3150-3152. https://doi.org/10.1093/bioinformatics/bts565
13. Grabherr MG, Haas BJ, Yassour M, et al. Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nat Biotechnol*. 2011;29(7):644-652. https://doi.org/10.1038/nbt.1883
14. Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357-362. https://doi.org/10.1038/s41586-020-2649-2
15. Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. UFBoot2: Improving the ultrafast bootstrap approximation. *Mol Biol Evol*. 2018;35(2):518-522. https://doi.org/10.1093/molbev/msx281
16. Hu K, Ni P, Xu M, et al. HiTE: a fast and accurate dynamic boundary adjustment approach for full-length transposable element detection and annotation. *Nat Commun*. 2024;15(1):5573. https://doi.org/10.1038/s41467-024-49912-8
17. Huang N, Li H. compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics*. 2023;39(10):btad595. https://doi.org/10.1093/bioinformatics/btad595
18. Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nat Methods*. 2017;14(6):587-589. https://doi.org/10.1038/nmeth.4285
19. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol*. 2013;30(4):772-780. https://doi.org/10.1093/molbev/mst010
20. Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. *PLoS One*. 2017;12(5):e0177459. https://doi.org/10.1371/journal.pone.0177459
21. Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. The Sequence Read Archive. *Nucleic Acids Res*. 2011;39(Database issue):D19-D21. https://doi.org/10.1093/nar/gkq1019
22. Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. *Mol Biol Evol*. 2021;38(10):4647-4654. https://doi.org/10.1093/molbev/msab199
23. Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol*. 2020;37(5):1530-1534. https://doi.org/10.1093/molbev/msaa015
24. Mistry J, Chuguransky S, Williams L, et al. Pfam: The protein families database in 2021. *Nucleic Acids Res*. 2021;49(D1):D412-D419. https://doi.org/10.1093/nar/gkaa913
25. NCBI. 2026a. NCBI Datasets. https://www.ncbi.nlm.nih.gov/datasets/ (accessed 20 November 2026).
26. NCBI. 2026b. EGAPx: Eukaryotic Genome Annotation Pipeline. https://github.com/ncbi/egapx (accessed 20 November 2026).
27. NCBI. 2026c. SRA Toolkit. https://github.com/ncbi/sra-tools (accessed 20 November 2026).
28. Ou S, Su W, Liao Y, et al. Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline. *Genome Biol*. 2019;20:275. https://doi.org/10.1186/s13059-019-1905-y
29. Python Software Foundation. 2026. Python 3.11. https://www.python.org (accessed 20 November 2026).
30. Schwengers O, Jelonek L, Dieckmann MA, et al. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom*. 2021;7(11):000685. https://doi.org/10.1099/mgen.0.000685
31. Simão FA, Waterhouse RM, Ioannidis P, Kriventseva EV, Zdobnov EM. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. *Bioinformatics*. 2015;31(19):3210-3212. https://doi.org/10.1093/bioinformatics/btv351
32. Törönen P, Medlar A, Holm L. PANNZER2: a rapid functional annotation web server. *Nucleic Acids Res*. 2018;46(W1):W84-W88. https://doi.org/10.1093/nar/gky350
33. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res*. 2025;53(D1):D609-D617. https://doi.org/10.1093/nar/gkae1010
34. Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods*. 2020;17(3):261-272. https://doi.org/10.1038/s41592-019-0686-2
