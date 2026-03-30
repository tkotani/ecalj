#!/bin/bash
# ecalj build script for nvfortran (gaugm.f90 ICE workaround)
set -e

cd ~/ecaljdeveloper
./InstallAll.py --gpu --fc nvfortran 2>&1 || true

cd ~/ecaljdeveloper/SRC/exec/build

# gaugm.f90 を各バリアント用に手動コンパイル (-O1, -acc/-Mvectなし)
echo "=== Workaround: compiling gaugm.f90 manually ==="

mpifort -O1 -fpic -Mbackslash -cpp \
  -module nvfortran/mod \
  -c ~/ecaljdeveloper/SRC/subroutines/gaugm.f90 \
  -o CMakeFiles/ecaljF.dir/home/takao/ecaljdeveloper/SRC/subroutines/gaugm.f90.o

mpifort -O1 -fpic -Mbackslash -cpp -D__MP \
  -module nvfortran/mod_mp \
  -c ~/ecaljdeveloper/SRC/subroutines/gaugm.f90 \
  -o CMakeFiles/ecaljF_mp.dir/home/takao/ecaljdeveloper/SRC/subroutines/gaugm.f90.o

mpifort -O1 -fpic -Mbackslash -cpp -D__GPU \
  -module nvfortran/mod_gpu \
  -c ~/ecaljdeveloper/SRC/subroutines/gaugm.f90 \
  -o CMakeFiles/ecaljF_gpu.dir/home/takao/ecaljdeveloper/SRC/subroutines/gaugm.f90.o

mpifort -O1 -fpic -Mbackslash -cpp -D__GPU -D__MP \
  -module nvfortran/mod_mp_gpu \
  -c ~/ecaljdeveloper/SRC/subroutines/gaugm.f90 \
  -o CMakeFiles/ecaljF_mp_gpu.dir/home/takao/ecaljdeveloper/SRC/subroutines/gaugm.f90.o

echo "=== Manual compile done. Freezing gaugm.f90 timestamp ==="
touch -t 202001010000 ~/ecaljdeveloper/SRC/subroutines/gaugm.f90

echo "=== Resuming build ==="
gmake -j32

echo "=== Restoring gaugm.f90 timestamp ==="
touch ~/ecaljdeveloper/SRC/subroutines/gaugm.f90

echo "=== Build complete ==="
