
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*iv_p_oipq/{
 print comment,$0
 gsub("integer , pointer ::","integer , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub("iv_p_oipq","iv_a_oipq")
 if (match($0,"nullify")) { 
  print "        if (allocated(sbz%iv_a_oipq)) deallocate(sbz%iv_a_oipq)"
  next; 
 }
 if (match($0,"allocate\\(sbz%iv_a_oipq")) {
  print "        if (allocated(sbz%iv_a_oipq)) deallocate(sbz%iv_a_oipq)"
 }

 gsub("associated\\(","allocated(")
}
{print}
