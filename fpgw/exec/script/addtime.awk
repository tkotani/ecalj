BEGIN {
   start=START
   poplist=0
   thisprogram="addtime.awk"
   nokey="__NOKEY__"
}
/^\!TIME0/{ 
  original=$0
  key=$1
  if (length(key)>6) {
    gsub("[!]TIME0","",key);
  }
  else {
   key=nokey
  } 
#  print "debug: key=",key
  printf("       call realtimediff(%i,'')\n",start); 
  push(start,key); start+=2;next 
}
/^\!TIME1/{ 
  filename=FILENAME
  original=$0
  key=$1
  if (length(key)>6) {
    gsub("[!]TIME1","",key);
  }
  else {
   key=nokey
  } 
  msg=sprintf("'%s_%s'",filename,key);
  if (NF>=2) {
     msg=$2
  }
#  print "debug: key=",key, "msg=",msg
  i=pop(key); 
  printf("       call realtimediff(%i,%s)\n",i,msg);next 
}
/^\!TIMESHOW/ { 
            print  "       call print_realtimediff()" ; next}
{ print }
func push(x,key) {
   list[poplist,1]=x 
   list[poplist,2]= key
   poplist++
}
func pop(key) {
   poplist--;
   i= list[poplist,1]
   k= list[poplist,2]
   if ( k!=key) { 
    print " ERROR message of ",thisprogram
    print " ERROR: inconsistent key, key=",key
    print " ERROR: but the last key in the stack is ",k
    print " ERROR: original line is"
    print original
    exit(1000)
   }
   return i
}
END{
 if (poplist!=0) {
  print ""
  print thisprogram,": possible error because poplist!=0"
  print "poplist="
  for (i=0;i<poplist;i++) {
    print i,list[i,1],list[i,2]
  }
 }
}
