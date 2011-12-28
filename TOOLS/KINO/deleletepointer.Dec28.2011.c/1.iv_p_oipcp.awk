BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*iv_p_oipcp/{
 print comment,$0
 if (match($0,"integer, *pointer")) { next }
 if (match($0,"=>")) { next }
 gsub("iv_p_oipcp","iv_p_oipc")
}
{print}
