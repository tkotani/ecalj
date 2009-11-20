#!/usr/bin/env python
from api import parse
#fff='./suham.F'
fff='./xxx.F'
#fff='./bndfp.F'
#fff='./x.F'
#/home/takao/ecal//lm-7.0betaK001/fp/bndfp.F'
#print fff
tree = parse(fff,isfree=False,isstrict=False,ignore_comments=False)
print tree.content
#print tree
