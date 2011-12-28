BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*zv_p_osmpot/{
 if (FNR<=180) {
   print comment,$0
   next
 } 
}
{print}
