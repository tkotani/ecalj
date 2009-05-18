#input: lstra.F
#output: m_struc_def.F
# $ gawk -f cont.awk lstra.F | gawk -f struc2.awk | gawk -f make_def.awk
#it outputs m_struc_def.F : The next program uses this.
#discard the stdout
BEGIN{
   spc="       "
   flag_expr=0
   flag_struc=0

   n_subnamelist=0

   prefix="struc_"

   strucfile="m_struc_def.F"

   sizefile="m_struc_size.F"

   print spc, "module m_struc_def" > strucfile
   print spc > strucfile
   print spc,"type s_ri" > strucfile
   print spc," integer:: t" >strucfile
   print spc," integer(8):: i" >strucfile
   print spc," real(8):: r" >strucfile
   print spc,"end type" > strucfile
   print spc > strucfile
}
/explanation/ {
   flag_expr=1
}
/^ +data/{
   flag_expr=0
   print "C",$0 > strucfile 
   next;
}
/^ *subroutine/{
   subdef=$0
}
/sort list/{ next; }
/^type=/{
#   print;
   subname=$2;
   n_name=$3;
   print "C",subname,n_name
   flag_struc=1
   idx=0

   n_subnamelist++;
   subnamelist[n_subnamelist]=subname
   next;
}
/^#/{ next; }
/^end_type *$/{
   if (n_name!=idx) {
      print "error size different", subname,n_name,idx
      exit(10);
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
     exit(10);
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
     exit(10);
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
   if (flag_expr==1) {
     print > strucfile
   }
   if (flag_struc==1) {
     idx++;
     kind[idx]=$3; 
     str=$4;
     gsub("\\.","",str);
     name[idx]=str;
     size[idx]=$5;
     print "C",$0 > strucfile 
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

