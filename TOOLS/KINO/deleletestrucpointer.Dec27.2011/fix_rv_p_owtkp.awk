
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*rv_p_owtkp/{
 print comment,$0
 gsub("real\\(8\\) , pointer ::","real(8) , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub("rv_p_owtkp","rv_a_owtkp")
 gsub("associated\\(sbz%rv_a_owtkp\\)","allocated(sbz%rv_a_owtkp)")
}
{print}
