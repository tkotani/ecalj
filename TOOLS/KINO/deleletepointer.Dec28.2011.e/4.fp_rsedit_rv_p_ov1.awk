BEGIN{
start=0
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "


}
/^ .*subroutine prsed4/{
  start=1
}
/^ .*rv_p_ov[01]/{
 if (start==1) {
   print comment,$0
   if (match($0,"real *\\(8\\) *, *pointer")) {next}
   if (match($0,"=>")){next}

   gsub("rv_p_ov1","ssite(ib)%rv_p_ov1")
   gsub("rv_p_ov0","ssite(ib)%rv_p_ov0")
 }
}
{print}
