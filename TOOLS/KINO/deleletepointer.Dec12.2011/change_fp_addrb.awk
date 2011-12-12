BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_oidxsh/{
 print comment,$0
 if (match($0,"pointer")) {next}
 if (match($0,"iv_p_oidxsh => sham%iv_p_oindxo")){next}
 gsub("iv_p_oidxsh","sham%iv_p_oindxo")
}
{
 print 
}
