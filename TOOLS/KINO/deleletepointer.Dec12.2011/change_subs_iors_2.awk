#please set LANG=C before invoking this script to show date in English
#input=  replace_subs_iors_2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_ov1/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ov1 *=> *ssite\\(ib\\)%rv_p_ov1")) { next;}
  if (match($0,"ssite\\(ib\\)%rv_p_ov1 *=> *rv_p_ov1")) { next;}
  gsub("rv_p_ov1","ssite(ib)%rv_p_ov1")
  if (match($0,"ssite\\(ib\\)%ssite\\(ib\\)%")) {next;}
}
/^ .*rv_p_ov0/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ov0 *=> *ssite\\(ib\\)%rv_p_ov0")) { next;}
  if (match($0,"ssite\\(ib\\)%rv_p_ov0 *=> *rv_p_ov0")) { next;}
  gsub("rv_p_ov0","ssite(ib)%rv_p_ov0")
  if (match($0,"ssite\\(ib\\)%ssite\\(ib\\)%")) {next;}
}
/^ .*rv_p_oves/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oves *=> *spot%rv_p_oves")) { next;}
  if (match($0,"spot%rv_p_oves *=> *rv_p_oves")) { next;}
  gsub("rv_p_oves","spot%rv_p_oves")
  if (match($0,"spot%spot%")) {next;}
}
{print}
