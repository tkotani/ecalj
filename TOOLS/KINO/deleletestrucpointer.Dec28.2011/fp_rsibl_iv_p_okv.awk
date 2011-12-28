BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/iv_p_okv/{
 print comment,$0
 gsub("iv_p_okv","iv_a_okv")
}
{print}
