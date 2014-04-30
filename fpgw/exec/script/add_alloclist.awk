BEGIN{ iomp=0 }
/^ .*allocate *[(]/{
  if (npop()) { popall() }
  ialloc=1
#  print 
  push($0)
  next
}
/^\!\$OMP *parallel/{
  if (npop()) { popall() }
  iomp++
  print "C omp level=",iomp
}
/^\!\$OMP *end *parallel/{
  if (npop()) { popall() }
  iomp--;
  print "C omp level=",iomp
}
/^ .*deallocate *[(]/{
  if (npop()) { popall() }
  ialloc=1
#  print
  push($0)
  next
}
/^     [^ ]/{
  if (ialloc) {
  #print 
  push($0)
  next
  }
}
{
  if (npop()) { popall() }
  ialloc=0
  print 
}
func push(str) {
  iline++
  line[iline]=str
}
func npop() { 
  return iline
}

func popall(   all,i,s) {
  all=""
  for (i=1;i<=iline;i++) {
    s=line[i]
    gsub("\!.*$","",s);
    all = all substr(s,7)
  }
  #print all
  process_allocate(all) 
  for (i=1;i<=iline;i++) {
    print line[i]
  }
  iline=0
}

func process_allocate(str ) {
  if (str~/deallocate *[(]/) { dealloc=1 }
  else { dealloc=0 }
  gsub("^.*allocate *[(]","",str);
  gsub("[)][ ;]*$","",str);
  #print "P>",str
  level0comma(str) 
}

func level0comma(str, name,istart, n,level,i,s,result,pos,dim,ndim ) {
 result[1]=0; result[2]=0
 n=length(str)
 level=0
 idx=0
# print
# print n
 pos[idx++]=istart
 for (i=istart;i<=n;i++) {
    s=substr(str,i,1)
    if (s=="(") { if(level==0) {result[1]=i}; level++ }
    if (s==")") { level--; if (level==0) {result[2]=i;} }
    if (s=="," && level==0) { #print ">",i; 
        pos[idx++]= i+1 }
#    print i,level,s
 }
 pos[idx++]=n+2;
 for (i=0;i<idx-1;i++) {
 #    print "i",substr(str,pos[i],  pos[i+1]-pos[i]-1)
     dim[i]=substr(str,pos[i],  pos[i+1]-pos[i]-1)
     sub(",$","",dim[i])
     gsub(" ","",dim[i])
 }
 ndim=idx-1
 for (i=0;i<ndim;i++) {
    if (dim[i]~/stat=/) { continue; }
    print "C--------------------------"
    print "C var",i,dim[i]#,dimsize1(dim[i])
    name=dim[i]
    sub("[(].*$","",name);
    print "C name=",name
    g_name=name
   if (dealloc) {
    print "C--------------------------"
    s=g_name
    if (iomp>0) { s= s"/omp"iomp }
      if (iomp>0) { print "!$OMP critical"}
      printf( "      call del_alloclist( \"%s\") !omplevel%i\n",s,iomp )
      if (iomp>0) { print "!$OMP end critical"}
   }
   else {
    s=dim[i]
    sub("^[^(]*[(]","",s)
    sub("[)] *$","",s)
    level0comma_size(s)
   }
 }
}

func level0comma_size(str, name,istart, n,level,i,s,result,pos,dim,ndim ) {

 result[1]=0; result[2]=0
 n=length(str)
 level=0
 idx=0
 pos[idx++]=istart
 for (i=istart;i<=n;i++) {
    s=substr(str,i,1)
    if (s=="(") { if(level==0) {result[1]=i}; level++ }
    if (s==")") { level--; if (level==0) {result[2]=i;} }
    if (s=="," && level==0) { #print ">",i;
        pos[idx++]= i+1 }
#    print i,level,s
 }
 pos[idx++]=n+2;
 for (i=0;i<idx-1;i++) {
 #    print "i",substr(str,pos[i],  pos[i+1]-pos[i]-1)
     dim[i]=substr(str,pos[i],  pos[i+1]-pos[i]-1)
     sub(",$","",dim[i])
 }
 ndim=idx-1
 for (i=0;i<ndim;i++) {
    print "C dim",i,dim[i],dimsize1(dim[i])
 }
    print "C--------------------------"
 s=g_name
 if (iomp>0) { s= s"/omp"iomp }
 if (iomp>0) { print "!$OMP critical" }
 printf("      call add_alloclist(\"%s\",storage_size(%s)/8,\n",s,g_name) 
 printf("     & int(") 
 for (i=0;i<ndim;i++) {
   printf("(%s)",dimsize1(dim[i]))
   if (i!=ndim-1) { printf("*"); }
 }
 printf(",kind=8))  !omplevel%i\n",iomp);
 if (iomp>0) { print "!$OMP end critical" }
}


func dimsize1(str, n) {
 if (str~/:/) {
    n=split(str,a,":")
    return "(("a[2]")-("a[1]")+1)"
 }
 return str
}


#func process_allocate0(str,   n,i,a,name) {
#  
#  gsub("^     [^ ]","",str)
#  gsub("\!.*$","",str);
#  gsub("^ .*allocate[(]","",str);
#  gsub("[)][ ;]*$","",str);
#  gsub("[(][^(]*[)]"," ",str)
#  gsub("," , " ",str);
#  n=split(str,a)
#  for (i=1;i<=n;i++) {
#    name=a[i]
#    printf("%scall add_patable(\"%s\",int(storage_size(%s)/8),int8(dim(%s)))\n", spc,name,name,name)
#  }
#}

