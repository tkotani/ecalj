#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_lmv7util
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_ormax/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_ormax *=> *sarray%rv_p_ormax")) { print comment,$0;next;}
  if (match($0,"sarray%rv_p_ormax *=> *rv_p_ormax")) { print comment,$0;next;}
  gsub("rv_p_ormax","sarray%rv_p_ormax")
  gsub("sarray%sarray%","sarray%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_osop/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_osop *=> *spot%rv_p_osop")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_osop *=> *rv_p_osop")) { print comment,$0;next;}
  gsub("rv_p_osop","spot%rv_p_osop")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_opprel/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_opprel *=> *spot%rv_p_opprel")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_opprel *=> *rv_p_opprel")) { print comment,$0;next;}
  gsub("rv_p_opprel","spot%rv_p_opprel")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_orhrmx/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_orhrmx *=> *spot%rv_p_orhrmx")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_orhrmx *=> *rv_p_orhrmx")) { print comment,$0;next;}
  gsub("rv_p_orhrmx","spot%rv_p_orhrmx")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_orhos/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_orhos *=> *spot%rv_p_orhos")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_orhos *=> *rv_p_orhos")) { print comment,$0;next;}
  gsub("rv_p_orhos","spot%rv_p_orhos")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_ogrrme/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_ogrrme *=> *spot%rv_p_ogrrme")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_ogrrme *=> *rv_p_ogrrme")) { print comment,$0;next;}
  gsub("rv_p_ogrrme","spot%rv_p_ogrrme")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_oclabl/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_oclabl *=> *sarray%rv_p_oclabl")) { print comment,$0;next;}
  if (match($0,"sarray%rv_p_oclabl *=> *rv_p_oclabl")) { print comment,$0;next;}
  gsub("rv_p_oclabl","sarray%rv_p_oclabl")
  gsub("sarray%sarray%","sarray%")
  if (old != $0 ) {print comment,$0}
}
/^ .*rv_p_ovdif/{
  old=$0
  if (match($0,"pointer")) { print comment,$0;next;}
  if (match($0,"rv_p_ovdif *=> *spot%rv_p_ovdif")) { print comment,$0;next;}
  if (match($0,"spot%rv_p_ovdif *=> *rv_p_ovdif")) { print comment,$0;next;}
  gsub("rv_p_ovdif","spot%rv_p_ovdif")
  gsub("spot%spot%","spot%")
  if (old != $0 ) {print comment,$0}
}
{print}
