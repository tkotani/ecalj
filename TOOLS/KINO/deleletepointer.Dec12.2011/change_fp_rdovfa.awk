BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_okv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_okv *=> *slat%iv_p_okv")) { next;}
  gsub("iv_p_okv","slat%iv_p_okv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_ogv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ogv *=> *slat%rv_p_ogv")) { next;}
  gsub("rv_p_ogv","slat%rv_p_ogv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*zv_p_osmrho/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"zv_p_osmrho *=> *spot%zv_p_osmrho")) { next;}
  gsub("zv_p_osmrho","spot%zv_p_osmrho")
  if (match($0,"spot%spot%")) {next;}
}
{print}
