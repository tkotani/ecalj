# Sample files

Temporary file sizes will be large; therefore, it is recommended to perform calculations in the work directory.
Use a copied Sample directory for calculations to keep the original files.

Note inas2gasb2 means inas2gasb2_ohtaka or inas2gasb2_kugui.
Take diff inas2gasb2_ohtaka or inas2gasb2_kugui. It shows only the difference is the memory contral option given in GWinput.

```
cp -r ./inas2gasb2 YOUR WORK DIRECOTRY
```

### (InAs)_2/(GaSb)_2
- ./inas2gasb2
- InAs 2layers and GaSb 2layers repeated system
- 8 atoms/cell
- job submission command on ISSP system B
```bash
sbatch job_ohtaka.sh
```
- job submission command on ISSP system C
```bash
qsub job_kugui.sh
```
The GPU node will be used on `job_kugui.sh`

### (InAs)_4/(GaSb)_4
- ./inas4gasb4
- InAs 4layers and GaSb 4layers repeated system
- 16 atoms/cell
> [!NOTE]
> This calculation will not finish within 30 minutes, but debug queue (maximum running time:30 min) is set in`job_ohtaka.sh` and `job_kugui.sh`.
 
- job submission command on ISSP system B
```bash
sbatch job_ohtaka.sh
```
- job submission command on ISSP system C
```bash
qsub job_kugui.sh
```
