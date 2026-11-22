# Interpreting the results and avoiding common pitfalls

## 1. Distinguish three kinds of "tree"

| Tree | Produced by this pipeline? | What it can answer |
|---|---|---|
| Domain tree | Yes, one per HMM family | Evolutionary relationships among copies of the same domain across proteins/species; subfamily delimitation |
| Gene tree | No (requires full-length genes/proteins and handling of alternative splicing) | Gene duplication, loss, and orthology/paralogy |
| Species tree | No | Species divergence order; requires concatenated single-copy orthologs or coalescent methods |

Treating a domain tree directly as a species tree is the most common mistake. A single protein can contain several domains, and a single domain can occur in several different protein families through domain shuffling, gene fusion, or horizontal transfer.

## 2. How to read branch support

- **UFBoot:** ultrafast bootstrap support, ranging from 0 to 100. It measures the frequency with which a branch appears after resampling; it is not a posterior probability. `-bnni` reduces overconfidence under model misspecification.
- **SH-aLRT:** SH-like approximate likelihood ratio test, ranging from 0 to 100.
- A commonly used threshold is UFBoot >= 95 together with SH-aLRT >= 80. When the two measures conflict, check whether the branch is driven by a few long branches or by an ambiguous region of the alignment; do not report only the more favourable number.
- A low-support "subfamily" can still be reported, but it should be described as an **unresolved lineage** with its support values, not as a confirmed subfamily.

## 3. Three possible origins of a cluster

1. **Single-copy ortholog group:** present across species with at most one copy per species, and with a sample composition consistent with the species tree.
2. **Paralog group:** several copies in the same species, or a sample composition that conflicts with the species tree; may result from gene duplication, loss, or horizontal transfer.
3. **Artifact:** caused by long-branch attraction, alignment error, inappropriate domain-boundary truncation, or contaminated sequences. Check:
   - whether `mean_within` is unusually large;
   - whether `min_between` is close to 0 (no separation between clusters);
   - whether the cluster is monophyletic and well supported;
   - whether members are near-identical sequences from the same species (sequencing or assembly redundancy).

## 4. Thresholds are not objective boundaries

`CLUSTER_CUTOFF=0.50` is only a default. Always inspect `<family>.cutoff_scan.tsv`:

- If the number of clusters is essentially unchanged between 0.3 and 0.6 (a plateau), the result is robust.
- If the number of clusters changes sharply with the threshold, the data do not contain a clear subfamily structure, and groups should not be forced.
- The absolute value of silhouette depends on cluster size and tree shape. Use it only for **relative comparisons**; do not apply a universal rule such as ">0.5 means good clustering".

## 5. Effects of domain boundaries

- The default `ali` coordinates are the regions that HMMER actually aligned; for modular proteins, this is usually the cleanest choice.
- If the domain has conserved N- or C-terminal extensions, `ali` will truncate them. In that case, set `DOMAIN_COORD="env"` and `DOMAIN_FLANK=10-20`, but confirm with a sensitivity analysis that the topology does not change.
- Hits covering less than 70% of the model (HMM coverage < 0.7) are filtered by default. For fragmented genomes or distant homologues, the threshold can be lowered to 0.5, but this must be stated in the results.

## 6. Sample/species imbalance

If some of the 27 samples contain far more proteins than others, the trees can be influenced by artificial differences in sampling density:

- paralogous copies from the same species tend to cluster together and inflate that species' internal branch lengths;
- bootstrap support can be dominated by the species with the largest sample.

Recommendations:

- report `n_samples` and `sample_counts` for every family;
- for a sensitivity analysis, retain only one copy per species per family (the default is "one domain per protein", not "one copy per species");
- if necessary, use `KEEP_ALL_DOMAINS=0` together with per-species dereplication.

## 7. Long branches and saturation

Very long branches at deep nodes may indicate:

- sequences that are too divergent for a reliable alignment;
- substitution saturation, especially when only a short domain is used;
- long-branch attraction (LBA).

Handling:

- check the ratio of alignment columns before and after trimAl;
- rerun with a site-heterogeneous model such as `LG+C20+F+G` (`IQTREE_MODEL`);
- remove long-branch sequences in a sensitivity analysis, but do not report only the "prettier" tree obtained after removal.

## 8. What to report

At a minimum:

- the number of sequences in the 27 samples and the number of domains entering each family;
- the HMM source and version, and whether `--cut_ga` or an E-value threshold was used;
- HMM coverage, domain length, and the rule for retaining copies per protein;
- the MAFFT algorithm, trimAl mode, and number of columns after trimming;
- the IQ-TREE version, the model selected by ModelFinder, UFBoot/SH-aLRT values, and the random seed;
- the clustering method, distance definition, threshold, and threshold-scan results;
- the sample composition and support values of each cluster;
- the availability of data and code (repository URL and data accessions).

## 9. Special notes for full-Pfam mode

The output of a full-Pfam `hmmscan` is not "all protein families of the 27 species". It is the set of hits to **known Pfam domains** in each proteome. When interpreting the results:

1. **Families absent from Pfam do not appear.** This is an HMM-guided survey of known domains, not unsupervised clustering from sequence data; do not claim that "all protein families were discovered".
2. **A protein can hit several domain families**, and a domain family can occur in several different protein architectures. A domain-family tree describes the evolution of "domain copies", not a protein-family tree, and certainly not a species tree.
3. **Clan redundancy:** several HMMs within the same Pfam clan may cover the same domain. `03b --collapse-clans` retains only the highest-scoring model. If this is not done, merge the results by clan downstream; do not treat multiple trees for the same domain as independent evidence.
4. **Copy number:** a low `single_copy_fraction` indicates that the family has multiple copies in most samples and that the tree mixes paralogs. Such families can be used for domain-family clustering but not directly for orthology or species-tree analyses.
5. **Sample-coverage thresholds:** `MIN_SAMPLES=14` means "most species", whereas `22` means ">=80% conserved". Do not conflate the two. Absence may also reflect incomplete assembly or annotation rather than true loss.
6. **Computational scale:** `MAX_FAMILIES` is a first-round cap, not a biological threshold. Before increasing it, confirm that memory, disk, and IQ-TREE runtime are sufficient.
