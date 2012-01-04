BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}

/^ .*%rv_p_opos/{
 print comment,$0
 gsub("%rv_p_opos","%rv_a_opos")
}
{print}
