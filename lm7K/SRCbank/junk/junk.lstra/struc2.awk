#
# usage
# gawk -f cont.awk lstra.F | gawk -f thisfile
#
# internal temporary files
#   fname=sprintf("___list___%s",subname);
#  fname2=sprintf("___list2___%s",subname);
#  They are removed automatically if this script is successfully done. 
#
BEGIN{
  read_uspec();
  IGNORECASE=1;
  n_list=0;
  n_names=0;
  startname=0;
}
/^ +subroutine/{
 gsub("\\("," ");
 n=split($0,a);
 subname=a[2];
 print;
}
/^Cr +off/{
  print"Cr --- explanation of ",subname, "---";
  print;
  startname=1;
  n_names=0;
  n_list=0;
  while (1) {
    getline ;
    print
    pos= match($0,"----------");
    if (pos>0) { start=0; break; }
    parse_remark($0);
  }

}
/^C ------------/{
  if (startname) {
  n_names=n_list;  
  for (i=0;i<n_list;i++) {
    names[i]=list[i];
  }
  startname=0;
  }
}


/^ +data +ilists/ {
  str=$0
 gsub("ilists","",str);
   parse_lists(str);

   n_ilists=n_list;
   for(i=0;i<n_list;i++) {
      ilists[i]=list[i]; 
   }

}

/^ +data +casts/ {
  str=$0
 gsub("casts","",str);
   parse_lists(str)

   n_casts=n_list;
   for(i=0;i<n_list;i++) {
    casts[i]=list[i];
   }   

   write_type()
   n_names=0;
   n_ilists=0;
   n_casts=0;
}


function parse_lists(str,n,i,a) {
 print str;
 gsub("data","",str);
 gsub("/"," ",str);
 gsub(" ","",str);
 n=split(str,a,",");
 n_list=0
 for (i=1;i<=n;i++) {
    if (match(a[i],"\*")) {
    parse_ast(a[i]);
    }
    else {
      list[n_list]=a[i]; 
      n_list++;
    }
 }

}

function parse_ast(str, n,a,i) {
  # str= "10*2" and so on
  n=split(str,a,"\*");
  for (i=1;i<=a[1];i++) {
     list[n_list]=a[2];
     n_list++;
  }
}




function write_type()
{
  if (n_ilists==0 || n_casts==0  || n_names==0 ) {
    print "error: n=0",n_ilists,n_casts,n_names
    return;
  }
  if ( n_ilists!=n_casts+1 ) {
    print "error: n_ilists!=n_casts+1 ",n_ilists,n_casts
  }
  if ( n_names!=n_casts) {
    print "error: n_names!=n_casts",n_names,n_casts
  }
#  print n_ilists, n_casts;
# ilists, casts
  sort_list();
  print "type=",subname, n_casts;

  for (i=0;i<n_casts;i++) {
    com="";
    if (match(names[i],"^o")>0) { com="possible_structure";}
    if ( subname=="usite") {
       usite_nelt(names[i],ilists[i],ilists[i+1]);
    printf("%2i %3i %2i %10s %2i %s\n", i+1,ilists[i], casts[i],names[i],ilists[i+1]-ilists[i],com);
    } else if ( subname=="uoptic") {
       uoptic_nelt(names[i],ilists[i],ilists[i+1]);
    printf("%2i %3i %2i %10s %2i %s\n", i+1,ilists[i], casts[i],names[i],ilists[i+1]-ilists[i],com);
    } else  if ( subname=="uspec" ) {
       nelt=uspec_nelt(names[i],ilists[i],ilists[i+1]);
    printf("%2i %3i %2i %10s %2i %s\n", i+1,ilists[i], casts[i],names[i],nelt,com);
    } else {
    printf("%2i %3i %2i %10s %2i %s\n", i+1,ilists[i], casts[i],names[i],ilists[i+1]-ilists[i],com);
    }
 } 

 print "end_type"

}

function uoptic_nelt(name, i1,i2)
{
   if (name=="axes") {
      printf("#note: size_must_be_changed, %s , size=3*max(nchi2,0)\n",name);
      printf("#note: if (mode>6) print 'uoptic: increase size of axes'\n");
   }

}

function usite_nelt(name, i1,i2)
{
   if (name=="delta") {
      printf("#note: %s, size = nint( size, is_th )\n",name);
   }
}

function uspec_nelt(name,ilists ,ilists2 ,nelt0,nelt,inew)
{
      n0=10;nkap0=3;

      ochfa=4;op=48;opz=51;oq=52;ocxbs=5
      oidx=22;orbp=46;ongcut=40;ocoreq=7
      oidu=21;ojh=25;ouh=71;oiq1=23;oivso=24;ovso=73
      oorcfa=57;

   nelt0=ilists2-ilists 
   nelt=nelt0

   inew=0;

  il=search_uspec(name);
  if (ilists>70) { nelt=n0; inew=1 }
  if (ilists == oorcfa) {nelt = 2; inew=1 }
  if (il==ocxbs) { nelt = 3; inew=1 }
  if (il==op || il==opz || il==oq || il==ochfa) { nelt = 2*n0; inew=1 }
  if (il==oidx || il==ongcut) { nelt = nkap0*n0; inew=1 }
  if (il==orbp) { nelt = nkap0*n0*2; inew=1 }
  if (il==oidu || il==ojh || il==ouh || il==oiq1 || il==oivso || il==ovso) { nelt = 4; inew=1 }
  if (il==ocoreq) { nelt = 2; inew=1 }

  if (inew) {

  if (nelt!=nelt0) {
     printf ("#name= %s , size changed from %d to %d\n",name,nelt0,nelt)
   }
   else {
     printf ("#name= %s , size unchanged from %d to %d\n",name,nelt0,nelt)
   }

   }

  return nelt
}


function sort_list(  fname,fname2,i,ret,cmd){
  print "sort list"
  fname=sprintf("___list___%s",subname);
  fname2=sprintf("___list2___%s",subname);
  for (i=0;i<n_casts;i++) {
    print i+1,ilists[i], casts[i],names[i]>fname;
  }
  close(fname);
  cmd =sprintf("sort -k 2 -n %s >%s",fname,fname2);
   system(cmd);
  i=0;
  while (1) { 
    ret=getline < fname2 
    if (ret<=0) {break;}
    ilists[i] =$2
    casts[i]=$3;
    names[i]=$4;
    i++;
  }
  close(fname2)
  cmd=sprintf("rm -f %s %s",fname,fname2); 
  system(cmd);
}

function parse_remark(str,   str1,str2,n,str3,a)
{
  str1=substr(str,3);
  str2=substr(str1,1,12);
  n=split(str2,a);
  if (n==0) { 
    return;
  }
  str3=substr(str1,13);
  n=split(str3,a);
  list[n_list]=a[1]; 
#  print n_list,a[1];
  n_list++;
}


func read_uspec( file,i,ret,id,name)
{
  file = "uspec"
  i=0
  while (1) {
    ret=getline < file;
    if (ret<=0) { break; }
    i++
    id=$3;
    name=$4;
    list_uspec[i,1]=id;
    list_uspec[i,2]=name;
  }
  n_list_uspec=i;
  close(file)
}


func search_uspec(name ,i)
{
  for (i=1;i<=n_list_uspec;i++) {
     if (name==list_uspec[i,2]) { return list_uspec[i,1]; }
  }
  print "error: search_uspec: failed to find",name;
  exit(10);
}



