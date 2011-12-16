#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_fp_lmfp_3
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_opos/{
  if ( match($0,"real\\(8\\),pointer") ) { 
     str=$0
     gsub("real.*$","",str)
     print comment, "add allocatable"
     print str "real(8),allocatable:: rv_a_opos(:)"
  }
  if ( (FNR>=1019 && FNR<=1259) || (FNR>=1264 && FNR<=1297)  || FNR==1328) {
  print comment,$0
#  if (match($0,"pointer")) { next;}
#  if (match($0,"rv_p_opos *=> *slat%rv_p_opos")) { next;}
#  if (match($0,"slat%rv_p_opos *=> *rv_p_opos")) { next;}
  gsub("rv_p_opos","rv_a_opos")
  gsub("associated","allocated")
  if (match($0,"[^e]allocate[(]")) {
    str=$0
    gsub("allocate.*$","",str)
    print str"if (allocated(rv_a_opos)) deallocate(rv_a_opos)"
  }
#  if (match($0,"slat%slat%")) {next;}
  }
}
{print $0}
