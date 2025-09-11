#!/usr/bin/env python3
import argparse,os, subprocess
parser = argparse.ArgumentParser(prog='lmfham')
parser.add_argument('--ctrl', required=True, help='foobar for ctrl.foobar')
args = parser.parse_args()
ctrl = args.ctrl 

def run(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

run(f"lmfa {ctrl} > llmf")
run(f"mpirun -np 8 lmf {ctrl} >> llmf")
run(f"job_band {ctrl} -np 8")
os.system(f"job_ham  {ctrl} -np 8")
run(f"lmfham2 {ctrl} --job=0")
run(f"lmfham2 {ctrl} --job=1")
run("gnuplot -p bandplot_MLO.isp1.glt")
