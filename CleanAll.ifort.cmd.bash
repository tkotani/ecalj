#!/bin/bash

cd lm7K/
make PLATFORM=ifort_mpik.cmd cleanall
make PLATFORM=ifort.cmd cleanall
cd ../fpgw/exec/
make PLATFORM=ifort.cmd cleanall
cd ../Wannier/
make PLATFORM=ifort.cmd cleanall



