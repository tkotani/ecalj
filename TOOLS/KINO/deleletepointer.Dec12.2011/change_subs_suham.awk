BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_oag/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oag *=> *slat%rv_p_oag")) { next;}
  gsub("rv_p_oag","slat%rv_p_oag")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *slat%rv_p_opos")) { next;}
  gsub("rv_p_opos","slat%rv_p_opos")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oipc/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oipc *=> *sarray%iv_p_oipc")) { next;}
  gsub("iv_p_oipc","sarray%iv_p_oipc")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*iv_p_oips/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oips *=> *sarray%iv_p_oips")) { next;}
  gsub("iv_p_oips","sarray%iv_p_oips")
  if (match($0,"sarray%sarray%")) {next;}
}
{print}
