BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*iv_p_oipc/{
 print comment,$0
 gsub("iv_p_oipc","iv_a_oipc")
 if (match($0,"integer, *pointer")) {
   gsub("integer, *pointer","integer, allocatable")
   gsub("=> *NULL\\(\\)","")
 }
 gsub("associated[(]sarray%iv_a_oipc","allocated(sarray%iv_a_oipc")
}
{print}
