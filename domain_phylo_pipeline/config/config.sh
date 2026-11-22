#!/usr/bin/env bash
# =============================================================================
#  config.sh -- single parameter entry point for the whole pipeline
#
#  Every step (01_prepare -> 06_report) sources this file. Change parameters
#  here only; no step script needs to be edited. Some settings can also be
#  overridden on the command line.
# =============================================================================

# ---------- Paths ------------------------------------------------------------
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${CONFIG_DIR}/.." && pwd)"

RAW_DIR="${PROJECT_DIR}/00_raw"                 # input: 27 raw pep.fasta files
HMM_DIR="${PROJECT_DIR}/01_hmms"                # hmmsearch mode: selected .hmm files
SAMPLES_TSV="${PROJECT_DIR}/config/samples.tsv" # sample table (can be generated)
RESULTS="${PROJECT_DIR}/results"

PREP_DIR="${RESULTS}/00_prep"
SEARCH_DIR="${RESULTS}/02_search"
DOMAIN_DIR="${RESULTS}/03_domains"
PHYLO_DIR="${RESULTS}/04_phylo"
CLUSTER_DIR="${RESULTS}/05_clusters"
REPORT_DIR="${RESULTS}/06_report"
LOG_DIR="${RESULTS}/logs"

# hmmsearch mode: one run per .hmm file in HMM_DIR times one run per sample.
# hmmscan mode:   one run per sample against a pressed HMM database
#                 (suitable for the full Pfam-A library).
SEARCH_MODE="${SEARCH_MODE:-hmmsearch}"          # hmmsearch | hmmscan
HMM_SCAN_DB="${HMM_SCAN_DB:-${PROJECT_DIR}/01_hmms/Pfam-A.hmm}"  # hmmscan mode; must be hmmpress'ed

# ---------- Parallelism and resources ----------------------------------------
JOBS="${JOBS:-4}"                    # concurrent tasks (HMMER / MAFFT / IQ-TREE)
THREADS_PER_JOB="${THREADS_PER_JOB:-4}"  # threads per IQ-TREE task
HMMSEARCH_CPU="${HMMSEARCH_CPU:-1}"  # CPUs per hmmsearch task

# ---------- Input cleaning ---------------------------------------------------
MIN_PROTEIN_LEN=30        # sequences shorter than this are discarded (aa)
MAX_X_FRAC=0.50           # sequences with a higher X fraction are discarded
CLEAN_UO_TO_X=1           # convert U/O (selenocysteine/pyrrolysine) to X to avoid HMMER warnings
DROP_INTERNAL_STOP=0      # 0 = convert internal stops to X; 1 = discard the sequence

# ---------- HMMER search and domain filtering --------------------------------
HMM_CUT_MODE="ga"         # ga (Pfam gathering threshold) | nc | evalue
HMM_EVALUE_MAX="1e-3"     # fallback ceiling in ga/nc mode; 1e-5 is recommended in evalue mode
HMM_COV_MIN=0.70          # minimum domain coverage of the HMM model length
DOMAIN_COORD="ali"        # extract domains using ali (alignment) or env (envelope) coordinates
DOMAIN_FLANK=0            # extra residues retained on each side of the domain
DOMAIN_MIN_LEN=30         # minimum domain length (aa)
DOMAIN_MAX_X_FRAC=0.30    # maximum X fraction in the extracted domain
KEEP_ALL_DOMAINS=0        # 0 = keep only the best domain copy per protein per family
                          # 1 = keep all tandem copies that pass the filters (IDs get _d2/_d3)

# Family-level inclusion criteria
MIN_SEQS="${MIN_SEQS:-4}"          # fewer than 4 sequences cannot support a meaningful unrooted tree
MIN_SAMPLES="${MIN_SAMPLES:-2}"    # minimum number of samples in which the family occurs
FAMILY_LIST="${FAMILY_LIST:-}"     # optional: analyse only listed families, one per line
                                   # (file name or Pfam name)

# ---------- Two-stage screening for full Pfam scans (03b_screen_families.py) ----------
# With 27 species: most species >= 14; conserved families (>=80%) >= 22.
MAX_FAMILIES="${MAX_FAMILIES:-0}"              # 0 = no limit; 200-500 is recommended for the first round
MIN_MEDIAN_COV="${MIN_MEDIAN_COV:-0.75}"       # lower bound on median HMM coverage within a family
MIN_MEDIAN_DOMAIN_LEN="${MIN_MEDIAN_DOMAIN_LEN:-30}"  # lower bound on median domain length
MAX_MEDIAN_COPIES="${MAX_MEDIAN_COPIES:-0}"    # 0 = no limit; 1 = single-copy families only
MIN_SINGLE_COPY_FRAC="${MIN_SINGLE_COPY_FRAC:-0.0}"  # lower bound on single-copy sample fraction, 0-1
CLAN_MAP="${CLAN_MAP:-}"                       # Pfam-A.clans.tsv (optional; removes clan redundancy)
COLLAPSE_CLANS="${COLLAPSE_CLANS:-0}"          # 1 = keep only the highest-scoring model per clan

# ---------- Sequence dereplication (only for very large families) ------------
DEREPLICATE=1
DEREP_ID=0.99             # CD-HIT identity threshold
DEREP_MAX_SEQS=1500       # dereplicate when n exceeds this value; 0 = never dereplicate

# ---------- Multiple sequence alignment --------------------------------------
ALIGNER="mafft"           # mafft | hmmalign
MAFFT_MODE="auto"         # auto | linsi | einsi | ginsi
                          # auto: L-INS-i for n <= 200, otherwise --auto
MAFFT_MAXITERATE=1000
MAFFT_EXTRA=""            # for example "--anysymbol"

# hmmalign mode needs to fetch a single HMM by family name. HMM_SOURCE can be
# Pfam-A.hmm or any .hmm file that contains all models. Leave empty to fall
# back to MAFFT automatically.
HMM_SOURCE=""

# ---------- Alignment trimming -----------------------------------------------
TRIMAL_MODE="gappyout"    # gappyout | automated1 | strict | custom arguments (e.g. "-gt 0.8")
TRIMAL_MIN_COL=30         # discard the family if fewer columns remain after trimming
MIN_COL_FRACTION=0.30     # warn if retained columns / original columns fall below this value

# ---------- IQ-TREE ----------------------------------------------------------
IQTREE_BIN="iqtree2"      # change to iqtree if only the iqtree binary is available
IQTREE_MODEL="MFP"        # ModelFinder Plus; a fixed model such as "LG+F+G4" is also allowed
IQTREE_MSET=""            # restrict candidate models, e.g. "LG,WAG,JTT,JTTDCMut"
UFBOOT=1000
ALRT=1000
IQTREE_SEED=12345
OUTGROUP=""               # optional: root used for display; clustering itself uses unrooted trees
IQTREE_EXTRA=""           # for example "--runs 2"

# ---------- Tree-based clustering --------------------------------------------
CLUSTER_CUTOFF="${CLUSTER_CUTOFF:-0.50}"   # patristic distance threshold (substitutions per site)
# complete: guarantees that the maximum pairwise patristic distance within a
# cluster is <= t, close to the TreeCluster diameter definition. average tends
# to chain and is less controllable when comparing thousands of families.
CLUSTER_LINKAGE="${CLUSTER_LINKAGE:-complete}" # complete | average | single
MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-2}"
CLUSTER_SCAN="${CLUSTER_SCAN:-0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.75,1.00}"
MAX_REPORT_ROWS="${MAX_REPORT_ROWS:-200}"  # maximum number of rows listed in the Markdown report

# ---------- Run control ------------------------------------------------------
FORCE="${FORCE:-0}"             # 1 = overwrite existing results
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-1}"  # 0 = delete MAFFT/trimAl intermediate files (not recommended)

# ---------- Export to child processes ----------------------------------------
export PROJECT_DIR RAW_DIR HMM_DIR SAMPLES_TSV RESULTS PREP_DIR SEARCH_DIR \
       DOMAIN_DIR PHYLO_DIR CLUSTER_DIR REPORT_DIR LOG_DIR SEARCH_MODE HMM_SCAN_DB \
       JOBS THREADS_PER_JOB HMMSEARCH_CPU MIN_PROTEIN_LEN MAX_X_FRAC CLEAN_UO_TO_X \
       DROP_INTERNAL_STOP HMM_CUT_MODE HMM_EVALUE_MAX HMM_COV_MIN DOMAIN_COORD \
       DOMAIN_FLANK DOMAIN_MIN_LEN DOMAIN_MAX_X_FRAC KEEP_ALL_DOMAINS MIN_SEQS \
       MIN_SAMPLES FAMILY_LIST MAX_FAMILIES MIN_MEDIAN_COV MIN_MEDIAN_DOMAIN_LEN \
       MAX_MEDIAN_COPIES MIN_SINGLE_COPY_FRAC CLAN_MAP COLLAPSE_CLANS \
       DEREPLICATE DEREP_ID DEREP_MAX_SEQS \
       ALIGNER MAFFT_MODE MAFFT_MAXITERATE MAFFT_EXTRA HMM_SOURCE TRIMAL_MODE \
       TRIMAL_MIN_COL MIN_COL_FRACTION IQTREE_BIN IQTREE_MODEL IQTREE_MSET \
       UFBOOT ALRT IQTREE_SEED OUTGROUP IQTREE_EXTRA CLUSTER_CUTOFF \
       CLUSTER_LINKAGE MIN_CLUSTER_SIZE CLUSTER_SCAN FORCE KEEP_INTERMEDIATE \
       MAX_REPORT_ROWS
