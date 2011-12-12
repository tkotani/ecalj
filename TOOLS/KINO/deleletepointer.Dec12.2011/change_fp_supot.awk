BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *slat%rv_p_opos")) { next;}
  gsub("rv_p_opos","slat%rv_p_opos")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_odlv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_odlv *=> *slat%rv_p_odlv")) { next;}
  gsub("rv_p_odlv","slat%rv_p_odlv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oag/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oag *=> *slat%rv_p_oag")) { next;}
  gsub("rv_p_oag","slat%rv_p_oag")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oqlv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqlv *=> *slat%rv_p_oqlv")) { next;}
  gsub("rv_p_oqlv","slat%rv_p_oqlv")
  if (match($0,"slat%slat%")) {next;}
}
{print}
