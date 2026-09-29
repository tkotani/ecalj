# Job scripts for the ISSP supercomputers

`job_kugui.sh` (system C, a GPU node, PBS) and `job_ohtaka.sh` (system B, CPU, Slurm) run one QSGW iteration of
`../inas2gasb2` (8 atoms) or `../inas4gasb4` (16 atoms). The installation is described in ecaljdoc
[installISSP](https://ecalj.github.io/ecaljdoc/install/installISSP).

```bash
cp -r ../inas2gasb2 <your work directory>/        # the work files are large
cp job_kugui.sh <your work directory>/inas2gasb2/
cd <your work directory>/inas2gasb2 && qsub job_kugui.sh      # ohtaka: sbatch job_ohtaka.sh
```

The run has ended well when `lgwsc` ends with `OK! ==== All calclation finished for  gwsc ====`; compare `QPU.1run`
with the one of the sample.

The scripts were brought to the present input (`ctrlg.<sname>.toml`) and options on 2026-09-30 and were not run on the
ISSP machines since.
