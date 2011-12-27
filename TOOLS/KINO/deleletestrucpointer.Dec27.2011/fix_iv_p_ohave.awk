
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*iv_p_ohave/{
 print comment,$0
 gsub("integer , pointer ::","integer , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub("iv_p_ohave","iv_a_oipq")
 if (match($0,"nullify")) { 
  print "        if (allocated(sbz%iv_a_oipq)) deallocate(sbz%iv_a_oipq)"
  next; 
 }
 gsub("associated\\(","allocated(")
}
{print}
