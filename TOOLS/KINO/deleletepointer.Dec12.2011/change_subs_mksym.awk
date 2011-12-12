BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *sarray%rv_p_opos")) { next;}
  gsub("rv_p_opos","sarray%rv_p_opos")
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
