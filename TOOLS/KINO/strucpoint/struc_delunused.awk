BEGIN { 
 IGNORECASE=1 
 g_in_struc=0

 spc6="      "
 comment="CCHSTRUC"
 structype_read()
}
/^ +type/{
 if (NF==2 && match($0,"\\(")==0) { g_in_struc=1 }
 if (g_in_struc==1) {
   g_strucname=$2
   print
   next
 }
}
/^ +end *type/{
 print
 g_in_struc=0
 next
} 

{
  if (g_in_struc) {
     gsub("::"," :: ")
     n=split($0,a)
     if (n>=3) {
        #print  a[1],"|",a[2],"|",a[3]
        id=structype_find(g_strucname,a[3])
        if (id>=1 && g_structype[id,3]=="undef" ) {
           print  comment,$0,"! not used because of undef"
           next
        }
#        if (id>=1 && g_structype[id,4]=="change") {
#            if ( g_structype[id,3]!="undef"){
#              print comment,$0
#              #print  g_structype[id,1],g_structype[id,2],g_structype[id,3]
#              prefix= get_prefix(g_structype[id,3])
#              if ( length(prefix)>0) {
#                print spc6,g_structype[id,3],", pointer :: ", prefix a[3],"(:)"
#              }
#            }
#            else {
#              print $0,"! not used"
#            }
#        }
#        else {
            print
#        }
     }
     else {
        print
     }
  }
  else {
     print
  }
}

func structype_read(  itotal,i)
{
  itotal=0
  file="struc_types"
  while (1) {
    ret=getline < file
    if (ret<=0) { break }
    if (match($0,"^#")==0) {
      itotal++
      for (i=1;i<=3;i++) {
        g_structype[itotal,i]=$i
      }
      if (NF>=5 && $5=="unchange") {
        g_structype[itotal,4]="unchage"
      }
      else {
        g_structype[itotal,4]="change"
      }
    }
  }
  close(file)
  g_n_structype=itotal
}

func structype_find(struc,elem,    i)
{
  for (i=1;i<=g_n_structype;i++) {
     if  (g_structype[i,1]==struc && g_structype[i,2]==elem) {
         return i
     }
  }
#  print "failed to find struc,elem=",struc,elem
#  exit(10)
  return -1
}



func get_prefix(type) {
 if (type=="integer") {
    return "iv_p_"
 }
 else if (type=="real(8)") {
    return "rv_p_"
 }
 else if (type=="complex(8)") {
    return "zv_p_"
 }
 else if (type=="undef") {
    return ""
 }
 else {
    print "unknown type:",type
    exit(10)
 }
}

