BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_ocy/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ocy *=> *slat%rv_p_ocy")) { next;}
  if (match($0,"slat%rv_p_ocy *=> *rv_p_ocy")) { next;}
  gsub("rv_p_ocy","slat%rv_p_ocy")
  if (match($0,"slat%slat%")) {next;}
}
{print}
