# MIDAS2_scripts
Scripts for MIDAS2 (Metagenomic Intra-Species Diversity Analysis System 2)
- call_consensus_v2.py (need to install python-lz4 by `pip install lz4`)
Performs the same calling as MIDAS. Classes are pasted and works as a single file.

- export_strainfacts.py
Convert MIDAS2 output to the input for StrainFacts. Discussion [here](https://github.com/bsmith89/StrainFacts/issues/5). This script output two files, `[CANDIDATE_SPECIES].metagenotype.tsv` containing ref and alt count for each position and `[CANDIDATE_SPECIES]_sample.pickle` containing dictionary of sample name and integer number.

Usage:
```bash
python export_strainfacts.py \
[SAMPLE_DIR] \
[CANDIDATE_SPECIES] \
[DB_DIR]
```
- `[SAMPLE_DIR]`: directory containing results of `midas2 run_snps` for each sample
- `[MERGE_DIR]`: output directory of `midas2 merge`
- `[CANDIDATE_SPECIES]` species name e.g., 100053
- `[DB_DIR]` used database directory e.g., MIDAS2/my_midasdb_gtdb
