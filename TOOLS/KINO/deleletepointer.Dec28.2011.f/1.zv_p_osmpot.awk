BEGIN {
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "

}
/^ .*zv_p_osmpot/{
  po="zv_p_osmpot"
  al="zv_a_osmpot"
  print comment,$0

  gsub(po,al)
  sub("associated","allocated")
  if (match($0,"allocate\\(") ) {
    s=$0
    sub("^.*allocate\\(","",s)
    sub("\\(.*$","",s)
    spc=$0
    sub("allocate.*$","",spc)
    print spc "if (allocated(" s ")) deallocate(" s ")"
  }
  
  if (match($0,"complex\\(8\\) *, *pointer")) {
    sub("pointer","allocatable")
    sub("=>NULL\\(\\)","")
  }
}
{print }
