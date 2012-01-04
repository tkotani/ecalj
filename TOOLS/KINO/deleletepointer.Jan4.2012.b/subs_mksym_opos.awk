BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}

/^ .*sarray%rv_p_opos/{
  print comment,$0
  gsub("sarray%rv_p_opos","slat%rv_p_opos")
}
{print}
