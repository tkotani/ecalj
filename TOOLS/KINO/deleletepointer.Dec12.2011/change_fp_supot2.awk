BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_omad/{
  print comment,$0
  if (match($0,"spot%rv_p_omad => rv_p_omad") ) {next}
  if (match($0,"pointer *::")){next}
  gsub("rv_p_omad","spot%rv_p_omad")
}
/^ .*rv_p_ogv/{
  print comment,$0
  if (match($0,"slat%rv_p_ogv => rv_p_ogv")) {next}
  if (match($0,"pointer *::")){next}
  gsub("rv_p_ogv","slat%rv_p_ogv")
}
/^ .*iv_p_okv/{
  print comment,$0
  if (match($0,"slat%iv_p_okv => iv_p_okv")) {next}
  if (match($0,"pointer *::")){next}
  gsub("iv_p_okv","slat%iv_p_okv")
}
/^ .*zv_p_osmrho/{
  print comment,$0
  if (match($0,"spot%zv_p_osmrho => zv_p_osmrho")){next}
  if (match($0,"pointer *::")){next}
  gsub("zv_p_osmrho","spot%zv_p_osmrho")
}
/^ .*zv_p_osmpot/{
  print comment,$0
  if (match($0,"spot%zv_p_osmpot => zv_p_osmpot")){next}
  if (match($0,"pointer *::")){next}
  gsub("zv_p_osmpot","spot%zv_p_osmpot")
}
/^ .*sv_p_oorhat/{
  print comment,$0
  if (match($0,"spot%sv_p_oorhat => sv_p_oorhat")){next}
  if (match($0,"pointer *::")){next}
  gsub("sv_p_oorhat","spot%sv_p_oorhat")
}
/^ .*iv_p_oips0/{
  print comment,$0
  if (match($0,"slat%iv_p_oips0 => iv_p_oips0")){next}
  if (match($0,"pointer *::")){next}
  gsub("iv_p_oips0","slat%iv_p_oips0")
}
/^ .*zv_p_obgv/{
  print comment,$0
  if (match($0,"slat%zv_p_obgv => zv_p_obgv")){next}
  if (match($0,"pointer *::")){next}
  gsub("zv_p_obgv","slat%zv_p_obgv")
}

{print}
