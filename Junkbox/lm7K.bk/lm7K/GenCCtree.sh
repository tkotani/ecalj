#!/bin/bash
# You need graphviz
sort callcaller.dotdata |uniq >temp
echo "digraph CCMap {" >ccmap.dot
echo 'graph [ rankdir = LR];' >>ccmap.dot
cat temp             >>ccmap.dot
echo "}"             >>ccmap.dot
dot -Tps ccmap.dot -o ccmap.ps
echo "CallCaller tree map (ccmap.ps) is generated."
