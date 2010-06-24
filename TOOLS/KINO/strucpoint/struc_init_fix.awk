BEGIN{
 IGNORECASE=1
 in_subroutine=0
 in_show_subroutine=0
 comment="CINITFIX"
 g_n_struc_name=0
}
/^ +subroutine/{
  in_subroutine=1
  if (match($0,"_show")) {
    in_show_subroutine=1
    print comment,$0
    next
  }
  print 
  next
}
/^ +end *subroutine/{
  in_subroutine=0
  if (in_show_subroutine) {
    print comment,$0
  }
  else {
    print
  }
  in_show_subroutine=0
  next
}
/^ +type\(/{
  if (in_show_subroutine) {
    print comment,$0
    next
  }
  print
  line=$0
  if (match(line,"::")){

    gsub("::"," :: ",line)
    gsub("\\("," ( ",line)
    gsub("\\)"," ) ",line)
    n=split(line,a)
    strucname=a[3]

  }
  next;
}
/^ +type/{
 if (in_show_subroutine) {
    print comment,$0
    next
 }
 if (NF==2 && match($0,"\\(")==0) { g_in_struc=1 }
 if (g_in_struc==1) {
   g_strucname=$2
   print
   next
 }
}
/^ +end *type/{
 if (in_show_subroutine) {
    print comment,$0
    next
 }
 print
 g_in_struc=0
 next
} 
{
  if (in_show_subroutine) {
    print comment,$0
    next
  }
  else if (in_subroutine) {
     if (match($0,"^ +.*%.*=")) {
        line=$0
        gsub("%"," % ",line)
        gsub("="," = ",line)
        n=split(line,a)
        name=a[3]
        gsub("\\(.*\\)","",name)
        id=struc_element_find(strucname,name)
        if (id>0) {
           print
        }
        else {
           print comment,$0
        }
     }
     else {
        print
     }
  }
  else if (g_in_struc) {
     print
     gsub("::"," :: ")
     n=split($0,a)
     if (n==3) {
        g_n_struc_name++
        g_struclist[g_n_struc_name,1]=g_strucname
        name=a[3]
        gsub("\\(.*\\)","",name)
        g_struclist[g_n_struc_name,2]=name
     }
  }
  else {
     print
  }
}

func struc_element_find(struc,elem)
{
   for (i=1;i<=g_n_struc_name;i++) {
     if (g_struclist[i,1]=struc &&  g_struclist[i,2]==elem) {
         return i
     }
   }
   return -1
}
