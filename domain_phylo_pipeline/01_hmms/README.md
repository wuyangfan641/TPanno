# HMM directory

When `SEARCH_MODE=hmmsearch`, place the selected `.hmm` files here.

Extract models from Pfam-A with:

```bash
bash bin/00_prepare_hmms.sh --pfam /path/to/Pfam-A.hmm \
    --names pfam_ids.txt --out 01_hmms/selected.hmm
```
