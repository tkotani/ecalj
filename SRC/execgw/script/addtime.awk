BEGIN {
   thisprogram="addtime.awk"
   start=1
   if (length(START)>0) {
     start=START
   }
   poplist=0
}
#
# format :
#  !TIME0_key  anycomment
#  !TIME1_key  "comment"
#  !TIMESHOW 
#
#  The nesting is allowed. Spaghettie jump code isn't allowed. 
#
#  You can change key of !TIME0_key and !TIME1_key. 
#
#  Don't put spaces in "comment".
#  e.g., this program doesn't accept "this is a comment". It accepts "this_is_a_comment."
#
/^\!TIME0/{ 
  original=$0
# extract key of !TIME0key 
  key=$1
  if (length(key)>7) {
    gsub("[!]TIME0_","",key);
  }
  else {
    print " ERROR message of ",thisprogram
    print " ERROR: no key for TIME0"
    print " ERROR: the original line is"
    print original
    exit(100)
  } 
#  print "debug: key=",key
  printf("       call realtimediff(%i,'')\n",start); 
  push(start,key); start+=2;next 
}
/^\!TIME1/{ 
  filename=FILENAME
  original=$0
# extract key of !TIME1key 
  key=$1
  if (length(key)>7) {
    gsub("[!]TIME1_","",key);
  }
  else {
    print " ERROR message of ",thisprogram
    print " ERROR: No key for TIME1"
    print " ERROR: The original line is"
    print original
    exit(200)
  } 
# get msg of !TIME1 msg
  if (NF>=2) {
     msg=$2
  }
  else {
    print " ERROR message of ",thisprogram
    print " ERROR: No message for TIME1"
    print " ERROR: The original line is"
    print original
    exit(210)
  }
# msg must be "..." or '...'
  flag=0
  if ( msg ~/^"[^"]*"$/ ) { flag=1 }
  if ( msg ~/^'[^']*'$/ ) { flag=1 }
  if (flag==0) {
    print " ERROR message of ",thisprogram
    print " ERROR: message for TIME1 must be \"...\"."
    print " ERROR: The original line is"
    print original
    exit(220)
  }
#  print "debug: key=",key, "msg=",msg
  i=pop(key); 
  printf("       call realtimediff(%i,%s)\n",i,msg);next 
}
/^\!TIMESHOW/ { print  "       call print_realtimediff()" ; next}
{ print }
func push(x,key) {
   list[poplist,1]=x 
   list[poplist,2]= key
   poplist++
}
func pop(key, i,k) {
   poplist--;
   i= list[poplist,1]
   k= list[poplist,2]
   if ( k!=key) { 
    print " ERROR message of ",thisprogram
    print " ERROR: Inconsistent key between TIME0 and TIME1, key=",key
    print " ERROR: The last key in the stack is ",k
    print " ERROR: The original line is"
    print original
    exit(1000)
   }
   return i
}
#
# This function isn't used. But I add it.
# It returns  "... ..." 
func get_string(s,start,  n,i,msg){
 n=split(s,a)
 is=start; ie=0
 for  (i=start;i<=n;i++) {
   print "i=",i,a[i]
   if (a[i]~/^"/) { is=i }
   if (a[i]~/"$/) { ie=i;break; }
 }
 if (ie==0) {
    return -1
 }
 msg=""
 for (i=is;i<=ie;i++) {
   msg = sprintf("%s %s",msg,a[i])
 }
 return msg
}
END{
 if (poplist!=0) {
  print ""
  print " ",thisprogram," shows an error at the last of the process."
  print " ",thisprogram,": A possible error in TIME because poplist!=0"
  print " The stack list="
  for (i=0;i<poplist;i++) {
    print i,list[i,1],list[i,2]
  }
 }
}

