#!/usr/bin/python

f = open("SEComg_int.UP",'r')
lines = f.readlines()
f.close()
f.close

ib = lines[0].split()[1]
ikp = lines[0].split()[2]
contents = ""

for line in lines:
      if line.strip() == "":
            continue
      line_sp = line.split()
      if ib == line_sp[1] and ikp == line_sp[2]:
            contents += line
      else :
            f = open("SEComg_int_n"+ib+"k"+ikp+".UP",'w')
            f.write(contents)
            f.close()
            contents = "" + line
            ib = line_sp[1]
            ikp = line_sp[2]
f = open("SEComg_int_n"+ib+"k"+ikp+".UP",'w')
f.write(contents)
f.close()
