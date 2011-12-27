#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_lmv7util_2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_ovrmax/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_ovrmax *=> *spot%rv_p_ovrmax")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_ovrmax *=> *rv_p_ovrmax")) { print comment,$0;next;}
  gsub("rv_p_ovrmax","spot%rv_p_ovrmax")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*iv_p_oics/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"iv_p_oics *=> *sarray%iv_p_oics")) { print comment,$0;next;}
  if (match($0,"sarray%iv_p_oics *=> *iv_p_oics")) { print comment,$0;next;}
  gsub("iv_p_oics","sarray%iv_p_oics")
  gsub("sarray%sarray%","sarray%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_oqnu/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_oqnu *=> *spot%rv_p_oqnu")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_oqnu *=> *rv_p_oqnu")) { print comment,$0;next;}
  gsub("rv_p_oqnu","spot%rv_p_oqnu")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_oqpp/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_oqpp *=> *spot%rv_p_oqpp")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_oqpp *=> *rv_p_oqpp")) { print comment,$0;next;}
  gsub("rv_p_oqpp","spot%rv_p_oqpp")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_oqt/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_oqt *=> *spot%rv_p_oqt")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_oqt *=> *rv_p_oqt")) { print comment,$0;next;}
  gsub("rv_p_oqt","spot%rv_p_oqt")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_obxc/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_obxc *=> *spot%rv_p_obxc")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_obxc *=> *rv_p_obxc")) { print comment,$0;next;}
  gsub("rv_p_obxc","spot%rv_p_obxc")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_oves/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_oves *=> *spot%rv_p_oves")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_oves *=> *rv_p_oves")) { print comment,$0;next;}
  gsub("rv_p_oves","spot%rv_p_oves")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_oqc/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_oqc *=> *spot%rv_p_oqc")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_oqc *=> *rv_p_oqc")) { print comment,$0;next;}
  gsub("rv_p_oqc","spot%rv_p_oqc")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_ova/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_ova *=> *spot%rv_p_ovintr")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_ovintr *=> *rv_p_ova")) { print comment,$0;next;}
  gsub("rv_p_ova","spot%rv_p_ovintr")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
{print}
