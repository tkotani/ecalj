# Using the GPU Version

The GPU can be utilized in QSGW calculations:
- The GPU is used for the calculation of self-energy and screened Coulomb interaction.
- `lmf` is executed on the CPU.

## Requirements:
1. NVIDIA GPU
2. Compiler: NVHPC with GPU-accelerated libraries: cuBLAS, cuSolver

## Compilation
```bash
FC=nvfortran ./InstallAll --gpu
```

## Execution
- The GPU version can be used by specifying the `--gpu` option with `gwsc`.
- The number of MPI processes for the GPU code can be specified with `-np2`. This should typically match the number of GPUs available.
Example:
```
gwsc -np 64 -np2 4 --gpu 1 $id > lgwsc
```

## Tips

1. The GPU version uses fewer MPI processes compared to the CPU version, so the available memory per MPI process is larger
   This allows for larger batch sizes in the calculation, potentially improving computation speed.
   The variables controlling the batch size are `NEMnmbatch` and `zeml_max_size`, which are specified in `GWinput`.
   For GPU calculations, set these values to around 4 (representing 4GB).
   ```text GWinput
   zmel_max_size 4
   MEMnmbatch 4
   ```
   > [!TIP]
   > If not set, the default values for CPU calculations will be used.

2. When handling large systems, add the following to `GWinput` to prevent memory exhaustion:
   ```
   keepEigen .false.
   ```

## Notes
- It is recommended to run GPU calculations on a single node.
- When using multiple nodes, ensure that MPI processes are correctly distributed across nodes.
