#
# gawk -f thisfunc < foo.F
# read 
#      type(...):: name 
# and 
#      call mpibc1(name2,...)
#
#  if name ==name2, then print it as ERROR
#
BEGIN {
  IGNORECASE=1
  g_n_typelist=0
}
/^[Cc#]/ { next; }
/^ *subroutine/{
  routine_name=$2
  g_n_typelist=0
  print $0
}
/^ *type/{
  line=$0
  gsub("::"," :: ",line)
  gsub("\\("," ( ",line)
  gsub("\\)"," ) ",line)
  n=split(line,a)
  print line,n
  for (i=1;i<=n;i++) {
     if (a[i]=="::") {
       ipos_sep=i
       break
     }
  }
  typelist_push( a[ipos_sep+1] )
}
/^ *call *mpibc1/{
  print $0
  line=$0
  gsub(","," , ",line)
  gsub("\\("," ( ",line)
  gsub("\\)"," ) ",line)
  n=split(line,a)
  if (n<4) {
     print "error n<4",$0
     exit(10)
  }
  if (a[1]=="call"  && a[2]=="mpibc1" && a[3]=="(") {
    name=a[4]
    typelist_check(name)
  } 
}

func typelist_push(name   )
{
  g_n_typelist++
  g_typelist[g_n_typelist]=name
}

func typelist_check(name )
{
  for (i=1;i<=g_n_typelist;i++) {
     if (g_typelist[i]==name) {
        print "ERROR, found", name, ", that is type, but in mpibc1 (",routine_name,")"
        print "ERROR, found", name, ", that is type, but in mpibc1 (",routine_name,")" > "/dev/stderr"
     }
  }
}


