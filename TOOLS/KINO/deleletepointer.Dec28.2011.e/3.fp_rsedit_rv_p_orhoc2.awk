BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "

}
/^ .*rv_p_orhoc2/{
  print comment,$0
  if (match($0,"real\\(8\\),pointer")){ next }
  if (match($0,"=>")) { next }
  gsub("rv_p_orhoc2","spec2(js)%rv_p_orhoc")
}
/^ .*rv_p_orhoc1/{
  print comment,$0
  if (match($0,"real\\(8\\),pointer")){ next }
  if (match($0,"=>")) { next }
  gsub("rv_p_orhoc1","spec1(is)%rv_p_orhoc")
}
/^ .*rv_p_ov01/{
  print comment,$0
  if (match($0,"real\\(8\\),pointer")){ next }
  if (match($0,"=>")) { next }
  gsub("rv_p_ov01","site1(ib)%rv_p_ov0")
}
/^ .*rv_p_ov11/{
  print comment,$0
  if (match($0,"real\\(8\\),pointer")){ next }
  if (match($0,"=>")) { next }
  gsub("rv_p_ov11","site1(ib)%rv_p_ov1")
}
/^ .*rv_p_ov02/{
  print comment,$0
  if (match($0,"real\\(8\\),pointer")){ next }
  if (match($0,"=>")) { next }
  gsub("rv_p_ov02", "site2(jb)%rv_p_ov0")
}
/^ .*rv_p_ov12/{
  print comment,$0
  if (match($0,"real\\(8\\),pointer")){ next }
  if (match($0,"=>")) { next }
  gsub("rv_p_ov12","site2(jb)%rv_p_ov1")
}
{print}



