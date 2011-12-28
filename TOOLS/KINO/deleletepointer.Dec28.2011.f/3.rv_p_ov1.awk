BEGIN {
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "

}


/^ .*rv_p_ov1/{
  print comment,$0
  po="rv_p_ov1"
  al="rv_a_ov1"
  gsub(po,al)
  if (match($0,"real\\(8\\) *, *pointer")) {
     sub("pointer","allocatable")
     sub("=> *NULL *\\(\\)","")
  }

  sub("associated","allocated")
  if (match($0,"[ )]allocate\\(")) {
    s=$0
    sub("^.*allocate\\(","",s)
    sub(al ".*$","",s)
    s= s al
    spc=$0
    sub("allocate.*$","",spc)
    print spc "if (allocated(" s ")) deallocate(" s ")"
  }
}
/^ .*rv_p_ov0/{
  print comment,$0
  po="rv_p_ov0"
  al="rv_a_ov0"
  gsub(po,al)
  if (match($0,"real\\(8\\) *, *pointer")) {
     sub("pointer","allocatable")
     sub("=> *NULL *\\(\\)","")
  }

  sub("associated","allocated")
  if (match($0,"[ )]allocate\\(")) {
    s=$0
    sub("^.*allocate\\(","",s)
    sub(al ".*$","",s)
    s= s al
    spc=$0
    sub("allocate.*$","",spc)
    print spc "if (allocated(" s ")) deallocate(" s ")"
  }
}

{print}
