BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_okv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_okv *=> *slat%iv_p_okv")) { next;}
  if (match($0,"slat%iv_p_okv *=> *iv_p_okv")) { next;}
  gsub("iv_p_okv","slat%iv_p_okv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_ogv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ogv *=> *slat%rv_p_ogv")) { next;}
  if (match($0,"slat%rv_p_ogv *=> *rv_p_ogv")) { next;}
  gsub("rv_p_ogv","slat%rv_p_ogv")
  if (match($0,"slat%slat%")) {next;}
}
{print}
