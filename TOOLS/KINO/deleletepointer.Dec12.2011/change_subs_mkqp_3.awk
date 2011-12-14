BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/iv_p_oipq/{
print comment,$0
if (match($0,"pointer")) { next }
if (match($0,"sbz%iv_p_oipq => iv_p_oipq")) {next}
gsub("iv_p_oipq","sbz%iv_p_oipq")
}
{print}
