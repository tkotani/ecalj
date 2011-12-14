BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}

/^ .*sham%iv_p_oiaxs/{
 print comment,$0
 gsub("sham%iv_p_oiaxs","sham%iv_a_oiaxs")
}
/^ .*sham%iv_p_ontabs/{
 print comment,$0
 gsub("sham%iv_p_ontabs","sham%iv_a_ontabs")
}
{print }
