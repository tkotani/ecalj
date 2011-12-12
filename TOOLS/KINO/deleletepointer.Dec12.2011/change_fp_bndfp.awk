BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_oqp/{
  print comment,$0
  if (match($0,"rv_p_oqp => sbz%rv_p_oqp")){next}
  if (match($0,"pointer *::")){next}
  gsub("rv_p_oqp","sbz%rv_p_oqp")
}
/^ .*rv_p_oswtk/{
  print comment,$0
  gsub("rv_p_oswtk","rv_a_oswtk")
}
/^ .*rv_p_owtkb/{
  print comment,$0
  gsub("rv_p_owtkb","rv_a_owtkb")
}


{print}
