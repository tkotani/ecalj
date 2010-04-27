#input m_struc_def.F
#output struc_sub.F
# $ gawk -f make_sub.awk m_struc_def.F  > struc_sub.F
#discard foo.F
BEGIN{
   IGNORECASE=1

   spc="       "
   flag_expr=0
   flag_struc=0

   n_subnamelist=0

   prefix="struc_"

   strucfile="foo.F"

   sizefile="foo.F"

   print spc, "module m_struc_def" > strucfile
   print spc > strucfile
   print spc,"type s_ri" > strucfile
   print spc," integer:: t" >strucfile
   print spc," integer(8):: i" >strucfile
   print spc," real(8):: r" >strucfile
   print spc,"end type" > strucfile
   print spc > strucfile
}
/^ +module/{ next; }
/^ +end +module/{ next; }
/^ +type s_ri/{ next; }
/^ +type .*$/{
#   print;
   subname=$2;
   sub("s_","u",subname);
#   n_name=$3;
   print "C",subname
   n_name=0
   flag_struc=1
   idx=0

   n_subnamelist++;
   subnamelist[n_subnamelist]=subname
   next;
}
/^#/{ next; }
/^C/{ next; }
/^ +end type .*$/{

   if (name[1]!="size") {
       print "Error: the first item must be 'size'"
       abortthisprogram(10);
   }
   if (kind[1]!=4) {
       print "Error: the first item must be 'real(8)'"
       abortthisprogram(10);
   }

   totalsize=1
   strucname=subname
   sub("^u","s_",strucname);
   print spc,"type",strucname > strucfile
   print spc,"real(8):: size" > strucfile
   for (i=1;i<=n_name;i++) {
      comment=""
      if ( name[i]=="size") { continue; }
      type="none"
      if (kind[i]==-1) { type= "integer(8)"; }
      if (kind[i]==2) { type= "integer(8)"}
      if (kind[i]==4) { type= "real(8)"}
      if (kind[i]==1) { type= "real(8)"; comment="!string"}

      str=name[i];
      gsub("\\.","",str)
      if (size[i]==1) {
      print spc,type,"::", str,comment > strucfile
      totalsize++
      }
      else {
      print spc,type,"::", str"("size[i]")",comment > strucfile
      totalsize+=size[i];
      }
   }   
   print spc,"end type",strucname > strucfile
   print spc > strucfile

    varname=strucname;
    sub("s_","v_",varname)

   newsubname=subname

   print spc
   print "C offe=-1 is only supported"
   if ( newsubname=="usite" || newsubname=="uspec" ) {
   print spc,"subroutine "newsubname"(struc,offe,nli,is,off,cast,nelt)"
   } else  {
   print spc,"subroutine "newsubname"(struc,offe,nli,off,cast,nelt)"
   }
   print spc, "use m_struc_def"
   print spc,"implicit none"
   print spc,"type("strucname"):: struc"
   print spc,"integer:: offe(1),nli,is,off,cast,nelt"
   print spc,"integer::val"
   print spc,"if (offe(1).eq.-1) then"
   print spc,"  call "newsubname"_init(struc)"
   print spc,"else"
   print spc,"  write(*,*) '"newsubname" offe=',offe(1),' not suppoted'"
   print spc,"  stop"
   print spc,"endif"
   print spc,"end subroutine "newsubname

#make subroutine
   print spc
    newsubname=subname"_size"
   print spc,"integer function "newsubname"()"
   print spc,"integer:: n"
   print spc,"n=",totalsize
       print "C  +2 is margin"
   print spc,"n=n+2"
   print spc,"if (mod(n,2).eq.1) n=n+1"
   print spc,newsubname"=n"
   print spc,"end function "newsubname

   print spc
    newsubname=subname"_show"
   print spc,"subroutine "newsubname"(struc)"
    print spc,"use m_struc_def"
   print spc,"implicit none"
    print spc,"type("strucname"):: struc"
   print spc,"write(*,*) struc"
   print spc,"end subroutine "newsubname
#

   print spc
    newsubname=subname"_init"
   print spc,"subroutine "newsubname"(struc)"
    print spc,"use m_struc_def"
   print spc,"implicit none"
    print spc,"type("strucname"):: struc"
    print spc,"integer:: "subname"_size"
   for (i=1;i<=n_name;i++) {
      if (name[i]=="size") { break; }
   }
      membername="struc%"name[i]
   if (kind[i]==1 || kind[i]==4 || kind[i]==-1 ) {  val ="0.0d0"}
   else if (kind[i]==2 ) { val=0}
   else {
     print "error kind[i]",kind[i],i
     abortthisprogram(10);
   }
   if (size[i]==1) { brace=""; } else  { brace="(:)"; }
   print spc,membername brace"="subname"_size()"


   for (i=1;i<=n_name;i++) {
      if (name[i]=="size") { continue; }
      membername="struc%"name[i]

   if (kind[i]==1 || kind[i]==4 ) {  val ="0.0d0"}
   else if (kind[i]==2) { val=0}
   else {
     print "error kind[i]",kind[i],i
     abortthisprogram(10);
   }
   if (size[i]==1) { brace=""; } else  { brace="(:)"; }
   print spc,membername brace"="val
   }
   print spc,"end subroutine "newsubname

#

   print spc
    newsubname=subname"_io"
    print spc, "subroutine "prefix newsubname"("varname",name,arr,rw,range1,range2,val,v_ri)"
    print "C use val or v_ri"
    print "C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*"
    print "C arr='s': use v_ri and return type and size of varname"
    print "C arr='a': use val and pack/unpack"
    print spc,"use m_struc_def"
    print spc,"use m_struc_func"
    print spc,"implicit none"
    print spc,"type("strucname"):: "varname
    print spc,"character(*),intent(in):: name"
    print spc,"character,intent(in):: rw,arr"
    print spc,"integer:: range1,range2"
    print spc,"integer:: val"
    print spc,"type(s_ri):: v_ri"
    print spc,"integer:: n,v_iv(10)"

   for (i=1;i<=n_name;i++) {
      if (name[i]=="size") { break; }
   }
      print spc, "select case (trim(adjustl(name)))"
      print spc, "case ('"name[i]"')"
      membername=strucname"%"name[i]
      sub("s_","v_",membername);
      print spc, "  if (rw.eq.'p') then"
      print spc, "  v_iv(1)=-1"
      print spc, "  call "subname"("varname",v_iv,0,0,0,0)"
      print spc, "  else"
      print spc, "  call "prefix"eval_io(",membername",rw,1,range1,range2,val)"
      print spc, "  endif"

   for (i=1;i<=n_name;i++) {
      if (name[i]=="size") { continue; }
      print spc, "case ('"name[i]"')"
      membername=strucname"%"name[i]
      sub("s_","v_",membername);

      sizename=strucname
      sub("s_","",sizename);
      varname=strucname
      sub("s_","v_",varname);
      sizename=prefix sizename"_"name[i]"_size"
      print spc, "  if (arr.eq.'1') then"
      print spc, "    if ( range2.eq.-1 ) range2 = range1"
      if (kind[i]==4 || kind[i]==1) {
      print spc, "    call "prefix"eval_io_r8_realbody(",membername",rw,1,range1,range2,v_ri%r)"
      print spc, "    v_ri%t=",kind[i]
      print spc, "    v_ri%i= nint(v_ri%r)"
      } else {
      print spc, "    call "prefix"eval_io_i8_realbody(",membername",rw,1,range1,range2,v_ri%i)"
      print spc, "    v_ri%t="kind[i]
      print spc, "    v_ri%r= v_ri%i"
      }
      print spc,"   else if (arr.eq.'s') then"
      print spc,"     v_ri%t=",kind[i]
      if (size[i]==1) {
      print spc,"     v_ri%i= 1"
      } else {
      print spc,"     v_ri%i= size("membername")"
      }
      print spc,"     v_ri%r= v_ri%i"
      print spc, "  else"
      if (size[i]==1) {
      print spc,"     n= 1"
      } else {
      print spc,"     n= size("membername")"
      }
      print spc, "    call "prefix"eval_io(",membername",rw,n,range1,range2,val)"
      print spc, "  endif"
   }
   print spc,"case default"
   print spc,"  write(*,*) '"prefix newsubname": failed to find ',name"
   print spc,"  stop"
   print spc, "end select"
   print spc,"end subroutine "prefix newsubname
   print spc


   
#end make subroutine

   flag_struc=0
   idx=0
   next;
}
{
   if (flag_struc==1) {
     gsub("::", " :: ");
     gsub("!", " ! ");
     gsub(" *[(]", "(");
     gsub(" *[)]", ")");
     print "C",$0;
     n=split($0,a);
     idx++;

     kind[idx]=0
     if (a[1]=="integer(8)") { kind[idx] = 2; }
     if (a[1]=="real(8)") { kind[idx] = 4; }
     if (a[5]=="string") { kind[idx] = 1; }
     if (kind[idx]==0) {
       print "error: kind must be integer(8) or real(8), or must specify 'string'"
       abortthisprogram(0);
     }

     sname=a[3]
     sub("[(][0-9]+[)]","",sname);
     name[idx]=sname

     sname=a[3]
     if ( match(sname,"[(]") ) {
        gsub("[(]", " ( ",sname);
        gsub("[)]", " ) ",sname);
        n=split(sname,b);
        if (n==4) {
           size[idx]=b[3];
        }
        else {
           print "error split ", sname
           print b[1],",",b[2],",",b[3],",",b[4]
           abortthisprogram(10);
        }
     }
     else {
        size[idx]=1;
     } 

     n_name++

#     print idx,kind[idx],name[idx],size[idx]
   }
}


END{
     print spc,"end module m_struc_def"> strucfile
     close(strucfile)

     print spc,"subroutine struc_packupack_val1(subname,name,struc,arr,rw,range1,range2,value,v_ri)"
     print spc,"use m_struc_def"
     print spc,"implicit none"
     print spc,"real(8):: struc"
     print spc,"integer:: range1,range2"
     print spc,"integer:: value"
     print spc,"character(*):: subname,name"
     print spc,"character:: rw,arr"
     print spc,"type(s_ri):: v_ri"
     print spc,"select case (trim(adjustl(subname)))"
  for (i=1;i<=n_subnamelist;i++) {
     str=subnamelist[i];
     sub("^u","",str);
     print spc,"  case ('"str"')"
     print spc,"    call "prefix subnamelist[i]"_io(struc,name,arr,rw,range1,range2,value,v_ri)"
  }
  print spc,"case default"
  print spc,"  write(*,*)'struc_packupack_val1: failed to find str=',subname"
  print spc,"  stop"
  print spc,"end select"
  print spc,"end subroutine struc_packupack_val1"
  print spc

}

