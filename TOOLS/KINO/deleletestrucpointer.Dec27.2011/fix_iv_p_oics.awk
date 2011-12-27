
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*iv_p_oics/{
 print comment,$0
 gsub("integer , pointer ::","integer , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub("iv_p_oics","iv_a_oics")
}
{print}
