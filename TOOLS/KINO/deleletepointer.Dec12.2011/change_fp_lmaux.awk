BEGIN {
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";

}

/^ .*rv_p_oqnu/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqnu *=> *spot%rv_p_oqnu")) { next;}
  gsub("rv_p_oqnu","spot%rv_p_oqnu")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*rv_p_oqc/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqc *=> *spot%rv_p_oqc")) { next;}
  gsub("rv_p_oqc","spot%rv_p_oqc")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*rv_p_opp/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opp *=> *spot%rv_p_opp")) { next;}
  gsub("rv_p_opp","spot%rv_p_opp")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*rv_p_opnu/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opnu *=> *spot%rv_p_opnu")) { next;}
  gsub("rv_p_opnu","spot%rv_p_opnu")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*rv_p_oqc/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqc *=> *spot%rv_p_oqc")) { next;}
  gsub("rv_p_oqc","spot%rv_p_oqc")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*iv_p_oinitc/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oinitc *=> *sarray%iv_p_ohave")) { next;}
  gsub("iv_p_oinitc","sarray%iv_p_ohave")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*rv_p_ovrmax/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_ovrmax *=> *spot%rv_p_ovrmax")) { next;}
  gsub("rv_p_ovrmax","spot%rv_p_ovrmax")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*rv_p_oves/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oves *=> *spot%rv_p_oves")) { next;}
  gsub("rv_p_oves","spot%rv_p_oves")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*rv_p_oqlv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqlv *=> *slat%rv_p_oqlv")) { next;}
  gsub("rv_p_oqlv","slat%rv_p_oqlv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_oics/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oics *=> *sarray%iv_p_oics")) { next;}
  gsub("iv_p_oics","sarray%iv_p_oics")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*iv_p_oics/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oics *=> *sarray%iv_p_oics")) { next;}
  gsub("iv_p_oics","sarray%iv_p_oics")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*rv_p_oclabl/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oclabl *=> *sarray%rv_p_oclabl")) { next;}
  gsub("rv_p_oclabl","sarray%rv_p_oclabl")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*rv_p_opos/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_opos *=> *sarray%rv_p_opos")) { next;}
  gsub("rv_p_opos","sarray%rv_p_opos")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*iv_p_oipc/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oipc *=> *sarray%iv_p_oipc")) { next;}
  gsub("iv_p_oipc","sarray%iv_p_oipc")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*rv_p_odlv/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_odlv *=> *slat%rv_p_odlv")) { next;}
  gsub("rv_p_odlv","slat%rv_p_odlv")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*rv_p_oag/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oag *=> *slat%rv_p_oag")) { next;}
  gsub("rv_p_oag","slat%rv_p_oag")
  if (match($0,"slat%slat%")) {next;}
}
/^ .*iv_p_onrc/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_onrc *=> *sarray%iv_p_onrc")) { next;}
  gsub("iv_p_onrc","sarray%iv_p_onrc")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*iv_p_oips/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_oips *=> *sarray%iv_p_oips")) { next;}
  gsub("iv_p_oips","sarray%iv_p_oips")
  if (match($0,"sarray%sarray%")) {next;}
}
/^ .*rv_p_orhrmx/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_orhrmx *=> *spot%rv_p_orhrmx")) { next;}
  gsub("rv_p_orhrmx","spot%rv_p_orhrmx")
  if (match($0,"spot%spot%")) {next;}
}
/^ .*rv_p_oqt/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"rv_p_oqt *=> *spot%rv_p_oqt")) { next;}
  gsub("rv_p_oqt","spot%rv_p_oqt")
  if (match($0,"spot%spot%")) {next;}
}
{print}
