BEGIN {
   start=START

   poplist=0
}
/^\!TIME0/{ 
  original=$0
  key=$1; 
  gsub("[!]TIME0","",key);
  printf("       call realtimediff(%i,'')\n",start); 
  push(start,key); start+=2;next 
}
/^\!TIME1/{ 
  original=$0
  delchar=$1
  key=$1; gsub("[!]TIME1","",key);
  gsub(delchar,"");
  i=pop(key); 
  printf("       call realtimediff(%i,%s)\n",i,$0);next 
}
/^\!TIMESHOW/ { 
            print  "       call print_realtimediff()" ; next}
#/^\!KINO/{ gsub("\!KINO","     "); print; next }
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
   if (k!=key) { 
    print " ERROR message of addtime.awk"
    print " ERROR: inconsistent key, key=",key
    print " ERROR: but the last key in the stack is ",k
    print " ERROR: original line is"
    print original
    exit(1000)
   }
   return i
}
