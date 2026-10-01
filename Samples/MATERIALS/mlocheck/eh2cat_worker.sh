#!/bin/bash
# takes materials from queue.txt once their baseline is DONE in mlocheck_auto_20261001/status_redo.log
cd $(dirname $0)
while true; do
  m=$(flock queue.lock bash -c 'm=$(head -1 queue.txt); [ -n "$m" ] && sed -i 1d queue.txt; echo $m')
  [ -z "$m" ] && break
  until grep -qE " $m DONE" /home/takao/work/mlocheck_auto_20261001/status_redo.log; do sleep 20; done
  ./eh2cat_one.sh $m $1
done
echo "$(date '+%F %T') eh2cat worker finished" >> status.log
