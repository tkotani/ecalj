BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_ogv/{
  print comment,$0
  gsub("rv_p_ogv","rv_a_ogv")
}
/^ .*rv_p_okv/{
  print comment,$0
  gsub("rv_p_okv","rv_a_okv")
}

{print}
