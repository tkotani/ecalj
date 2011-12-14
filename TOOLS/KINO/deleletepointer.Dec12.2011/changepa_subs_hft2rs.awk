BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_oiax/{
print comment,$0
sub("integer,pointer","integer,allocatable")
gsub("iv_p_oiax","iv_a_oiax")
}
/^ .*iv_p_ontab/{
print comment,$0
sub("integer,pointer","integer,allocatable")
gsub("iv_p_ontab","iv_a_ontab")
}
{print }
