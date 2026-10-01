#!/bin/bash
cd $(dirname $0)
while true; do
  m=$(flock redo_queue.lock bash -c 'm=$(head -1 redo_queue.txt); [ -n "$m" ] && sed -i 1d redo_queue.txt; echo $m')
  [ -z "$m" ] && break
  ./redo_one.sh $m $1
done
echo "$(date '+%F %T') redo worker finished" >> status_redo.log
