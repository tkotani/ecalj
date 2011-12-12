BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_odlv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_odlv *=> *slat%rv_p_odlv")) { next;}
  if (match($0,"slat%rv_p_odlv *=> *rv_p_odlv")) { next;}
  gsub("rv_p_odlv","slat%rv_p_odlv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oqlv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqlv *=> *slat%rv_p_oqlv")) { next;}
  if (match($0,"slat%rv_p_oqlv *=> *rv_p_oqlv")) { next;}
  gsub("rv_p_oqlv","slat%rv_p_oqlv")
  if (match($0,"slat%slat%")) {next;}
}
{print}
