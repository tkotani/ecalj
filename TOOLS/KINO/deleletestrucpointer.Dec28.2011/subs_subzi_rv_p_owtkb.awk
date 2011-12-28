BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/rv_p_owtkb/{
 print comment,$0
 gsub("rv_p_owtkb","rv_a_owtkb")
}
/rv_p_oswtk/{
 print comment,$0
 gsub("rv_p_oswtk","rv_a_oswtk")
}
{print}
