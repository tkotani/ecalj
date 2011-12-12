BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_oqp/{
  print comment,$0
  if (match($0,"sbz%rv_p_oqp => rv_p_oqp")){next}
  if (match($0,"pointer *::")){next}
  gsub("rv_p_oqp","sbz%rv_p_oqp")
}
/^ .*rv_p_owtkp/{
  print comment,$0
  if (match($0,"pointer *::")){next}
  if (match($0,"sbz%rv_p_owtkp => rv_p_owtkp")){next}
  gsub("rv_p_owtkp","sbz%rv_p_owtkp")
}
/^ .*iv_p_oidtet/{
  print comment,$0
  if (match($0,"pointer *::")){next}
  if (match($0,"nullify")){next}
  if (match($0,"sbz%iv_p_oidtet => iv_p_oidtet")){next}
  gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*iv_p_ogstar/{
  print comment,$0
  if (match($0,"pointer *::")){next}
  if (match($0,"sbz%iv_p_ostar => iv_p_ogstar")){next}
  gsub("iv_p_ogstar","sbz%iv_p_ostar")
}

{print}
