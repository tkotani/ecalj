BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_ojcg/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_ojcg *=> *slat%iv_p_ojcg")) { next;}
  if (match($0,"slat%iv_p_ojcg *=> *iv_p_ojcg")) { next;}
  gsub("iv_p_ojcg","slat%iv_p_ojcg")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oidxcg/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oidxcg *=> *slat%iv_p_oidxcg")) { next;}
  if (match($0,"slat%iv_p_oidxcg *=> *iv_p_oidxcg")) { next;}
  gsub("iv_p_oidxcg","slat%iv_p_oidxcg")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_ocy/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ocy *=> *slat%rv_p_ocy")) { next;}
  if (match($0,"slat%rv_p_ocy *=> *rv_p_ocy")) { next;}
  gsub("rv_p_ocy","slat%rv_p_ocy")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_ocg/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ocg *=> *slat%rv_p_ocg")) { next;}
  if (match($0,"slat%rv_p_ocg *=> *rv_p_ocg")) { next;}
  gsub("rv_p_ocg","slat%rv_p_ocg")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oag/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oag *=> *slat%rv_p_oag")) { next;}
  if (match($0,"slat%rv_p_oag *=> *rv_p_oag")) { next;}
  gsub("rv_p_oag","slat%rv_p_oag")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oistab/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oistab *=> *slat%iv_p_oistab")) { next;}
  if (match($0,"slat%iv_p_oistab *=> *iv_p_oistab")) { next;}
  gsub("iv_p_oistab","slat%iv_p_oistab")
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
/^ .*rv_p_ogv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ogv *=> *slat%rv_p_ogv")) { next;}
  if (match($0,"slat%rv_p_ogv *=> *rv_p_ogv")) { next;}
  gsub("rv_p_ogv","slat%rv_p_ogv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*zv_p_obgv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"zv_p_obgv *=> *slat%zv_p_obgv")) { next;}
  if (match($0,"slat%zv_p_obgv *=> *zv_p_obgv")) { next;}
  gsub("zv_p_obgv","slat%zv_p_obgv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oips0/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oips0 *=> *slat%iv_p_oips0")) { next;}
  if (match($0,"slat%iv_p_oips0 *=> *iv_p_oips0")) { next;}
  gsub("iv_p_oips0","slat%iv_p_oips0")
  if (match($0,"slat%slat%")) {next;}
}
{print}
