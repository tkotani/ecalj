#!/usr/bin/env python3
#edit glt 
#add newline to get graph

import sys,os

def edit(ma,mid):

    mate = ma
    mpid = mid

    a1 = "set terminal postscript enhanced color eps"
    a2 = "set output \"bandplot.isp1.band.eps\""


    matepath = mate + "/bandplot.isp1.glt"
    newmatepath = mate + "/new.glt"
    idpath = mpid + "/bandplot.isp1.glt"
    newidpath = mpid + "/new.glt"

    #print(os.path.isfile(matepath))
    #sys.exit()

    if os.path.isfile("bandplot.isp1.glt") == True:
        print("1")
        with open("bandplot.isp1.glt","r") as f:
            line = f.read()

        newline = a1 + "\n"+ a2 + "\n" + line

        with open("new.glt","w") as l:
            l.write(newline)     
     
    elif os.path.isfile(matepath) == True:
        print("2")
        with open(matepath,"r") as f:
            line = f.read()

        newline = a1 + "\n"+ a2 + "\n" + line

        with open(newmatepath,"w") as l:
            l.write(newline)

    elif os.path.isfile(idpath) == True:
        print("3")
        with open(idpath,"r") as f:
            line = f.read()

        newline = a1 + "\n"+ a2 + "\n" + line

        with open(newidpath,"w") as l:
            l.write(newline)

    else:
        print("4")
        print("where is bandplot.isp1.glt")
        print("")
