
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*rv_p_oclabl/{
 print comment,$0
 if (match($0,"nullify\\(sarray%rv_p_oclabl")) { 
    print "         if (allocated(sarray%rv_a_oclabl)) deallocate(sarray%rv_a_oclabl)"
    next
 }
 gsub("real\\(8\\) , pointer ::","real(8) , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub("rv_p_oclabl","rv_a_oclabl")
}
{print}
