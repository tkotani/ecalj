BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *slat%rv_p_opos")) { next;}
  if (match($0,"slat%rv_p_opos *=> *rv_p_opos")) { next;}
  gsub("rv_p_opos","slat%rv_p_opos")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_okv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_okv *=> *slat%iv_p_okv")) { next;}
  if (match($0,"slat%iv_p_okv *=> *iv_p_okv")) { next;}
  gsub("iv_p_okv","slat%iv_p_okv")
  if (match($0,"slat%slat%")) {next;}
}
{print}
