#!/bin/bash
# pops materials from queue.txt (flock) and runs them with np=$1
cd $(dirname $0)
while true; do
  m=$(flock queue.lock bash -c 'm=$(head -1 queue.txt); [ -n "$m" ] && sed -i 1d queue.txt; echo $m')
  [ -z "$m" ] && break
  ./run_one.sh $m $1
done
echo "$(date '+%F %T') worker np=$1 finished" >> status.log
