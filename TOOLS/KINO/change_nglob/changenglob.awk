BEGIN{ IGNORECASE=1; program="Changenglob" ; gvarname="globalvariables"}
/^ .*nglob *\(/{
  str=$0
  while ( match(str,"nglob") ) {
  str=change_nglob(str)
  }
  print program $0
  print str
  next
}
/^ .*dglob *\(/{
  str=$0
#  while ( match(str,"dglob") ) {
  str=change_dglob(str)
#  }
  print program $0
  print str
  next
}
{ print }
function change_dglob(original,   n,s0,s1,s2,s2a,s2b,s2c,s2c1,s2c2,s2c3,key1,
  s2c3a,s2c2b,s2c3c,s0a,s0b,s0c ,ret){
#  print original
  n=match(original,"dglob\\(")
#  print original,RSTART,RLENGTH
  s0=substr(original,0,RSTART-1)
  s1=substr(original,RSTART,RLENGTH)
  s2=substr(original,RSTART+RLENGTH)
#  print "0=",s0
#  print "1=",s1
#  print "2=",s2
#  print s0 s1 s2
  n=match(s2,",")
#  print n
  s2a=substr(s2,0,RSTART-1)
  s2b=substr(s2,RSTART,RLENGTH)
  s2c=substr(s2,RSTART+RLENGTH)
#  print "2a=",s2a   # key1
#  print "2b=",s2b   # , 
#  print "2c=",s2c
#  print s2a s2b s2c
  n=match(s2c,",")
  s2c1=substr(s2c,0,RSTART-1)
  s2c2=substr(s2c,RSTART,RLENGTH)
  s2c3=substr(s2c,RSTART+RLENGTH)
#  print "2c1=",s2c1  # key2
#  print "2c2=",s2c2  # ,
#  print "2c3=",s2c3  
#  print s2c1,s2c2,s2c3
  n=match(s2c3,"\\)")
  s2c3a=substr(s2c3,0,RSTART-1)
  s2c3b=substr(s2c3,RSTART,RLENGTH)
  s2c3c=substr(s2c3,RSTART+RLENGTH)
#  print "2c3a=",s2c3a # key3
#  print "2c3b=",s2c3b # )
#  print "2c3c=",s2c3c # other
#  print s2c3a s2c3b s2c3c

  key1=s2a
  gsub("'","",key1)

  n=match(s0,"[0-9a-zA-Z]")
  s0a=substr(s0,0,RSTART-1)
  s0b=substr(s0,RSTART,RLENGTH)
  s0c=substr(s0,RSTART+RLENGTH)

  ret= s0a gvarname"%" key1 " = " s2c1 s2c3c  "; "
  ret = ret  gvarname "%l_" key1 " = " gvarname"%l_" key1 " +1; "
  ret=ret  s0b s0c s2c1
  return ret
}
function change_nglob(original,   n,s0,s1,s2,s2a,s2b,s2c) {
  n=match(original,"nglob\\(")
#  print original,RSTART,RLENGTH
  s0=substr(original,0,RSTART-1)
  s1=substr(original,RSTART,RLENGTH)
  s2=substr(original,RSTART+RLENGTH)
#  print "0=",s0
#  print "1=",s1
#  print "2=",s2
#  print s0 s1 s2
  n=match(s2,"\\)")
#  print n
  s2a=substr(s2,0,RSTART-1)
  s2b=substr(s2,RSTART,RLENGTH)
  s2c=substr(s2,RSTART+RLENGTH)
#  print "2a=",s2a
#  print "2b=",s2b
#  print "2c=",s2c
#  print s2a s2b s2c
  gsub("'","",s2a)
#  print "=",s2a
  return  s0 gvarname "%" s2a s2c
   
}
