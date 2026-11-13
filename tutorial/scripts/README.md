# Instructor notes: rebuilding the workshop files

Students don't need anything here. Everything below is already built: the starting BAM, the "Lost?"
catch-up archives and the *SNRPN* track hub are live in the public bucket
(`https://s3.kopah.uw.edu/fiberseq-workshop`). Rebuild only if a step's output changes (e.g. a new
`ft` or FiberHMM version) or for a new workshop.

## Rebuild and upload

```bash
conda activate fiberseq-workshop && cd tutorial
HUB_EMAIL=<contact> bash scripts/make_backup.sh   # starting BAM, all steps, SNRPN hub → backup/
bash scripts/make_bucket.sh                       # → .bucket/ (starting BAM, chr15, stepNN.tar.gz, hub/)
rclone copy .bucket/ <remote>:<public bucket>/ --s3-acl public-read -P   # then set BUCKET= in the tutorial README
```

`.bucket/` and `backup/` step outputs are not tracked in git (see the top-level `.gitignore`).

## Scripts

| Script | Who | What it does |
|---|---|---|
| `check_setup.sh` | students | checks every tool on 20 fibers before the workshop (PASS / FAIL) |
| `cpg_by_hap.py` | tutorial step 10 | % methylated per CpG, all fibers and per haplotype |
| `fisher_hap_test.py` | tutorial step 9 | Fisher's exact test per peak, hap1 vs hap2 (re-creates the FIRE pipeline's test) |
| `make_input_bam.sh` | instructors | builds the starting BAM (*SNRPN* + *NAPA*) from the public GM12878 FIRE CRAM |
| `run_all.sh` | instructors | runs tutorial steps 1–11 end to end |
| `make_backup.sh` | instructors | starting BAM + `run_all.sh` + hub → `backup/` |
| `make_hub.sh` | instructors | builds the step 12a UCSC hub (see `../hub/README.md`) |
| `make_bucket.sh` | instructors | packages everything for the public bucket → `.bucket/` |

Each script starts with a comment block explaining its inputs, outputs and usage.
