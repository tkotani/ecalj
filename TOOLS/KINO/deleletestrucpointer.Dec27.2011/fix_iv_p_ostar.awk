
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*iv_p_ostar/{
 print comment,$0
 gsub("integer , pointer ::","integer , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub("iv_p_ostar","iv_a_ostar")
}
{print}
