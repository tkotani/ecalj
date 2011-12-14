BEGIN {
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}

/^ .*integer,pointer :: iv_p_owk/{
print  comment, $0
next
}
{print }
