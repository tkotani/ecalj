BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_ov0/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ov0 *=> *ssite\\(ib\\)%rv_p_ov0")) { next;}
  if (match($0,"ssite\\(ib\\)%rv_p_ov0 *=> *rv_p_ov0")) { next;}
  gsub("rv_p_ov0","ssite(ib)%rv_p_ov0")
  if (match($0,"ssite\\(ib\\)%ssite\\(ib\\)%")) {next;}
}
{print}
