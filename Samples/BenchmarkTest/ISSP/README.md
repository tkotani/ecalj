# Job scripts for the ISSP supercomputers

`job_kugui.sh` (system C, a GPU node, PBS) and `job_ohtaka.sh` (system B, CPU, Slurm) run one QSGW iteration of
`../inas2gasb2` (8 atoms) or `../inas4gasb4` (16 atoms). The installation is described in ecaljdoc
[installISSP](https://ecalj.github.io/ecaljdoc/install/installISSP).

```bash
cp -r ../inas2gasb2 <your work directory>/        # the work files are large
cp job_kugui.sh <your work directory>/inas2gasb2/
cd <your work directory>/inas2gasb2 && qsub job_kugui.sh      # ohtaka: sbatch job_ohtaka.sh
```

The run has ended well when `lgwsc` ends with `OK! ==== All calclation finished for  gwsc ====`. The sample contains
`QPU.1run`, so the iteration is counted as 2: compare the new `QPU.2run` with the sample's `QPU.1run`. The scripts pass
`-c kugui_gpu` / `-c ohtaka` (launcher settings in `SRC/exec/clusters.toml`; 2026-11-14).

The scripts were brought to the present input (`ctrlg.<sname>.toml`) and options on 2026-09-30 and were not run on the
ISSP machines since.
