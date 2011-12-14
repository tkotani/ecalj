BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_oiaxo/{
print comment,$0
sub("integer,pointer","integer,allocatable")
gsub("iv_p_oiaxo","iv_a_oiaxo")
sub("=> *NULL\\( *\\)","")
}
/^ .*iv_p_ontabo/{
print comment,$0
sub("integer,pointer","integer,allocatable")
gsub("iv_p_ontabo","iv_a_ontabo")
sub("=> *NULL\\( *\\)","")
}
{print }
