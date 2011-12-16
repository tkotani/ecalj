#please set LANG=C before invoking this script to show date in English
#input=  tmp/replace_subs_asados_2
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*rv_p_oqp/{
  if ( FNR<=447 ) {
    print comment,$0
    if (match($0,"pointer")) { 
       str=$0
       sub("real.*$","",str)
       print str"real(8),allocatable :: rv_a_oqp(:)"
       next
    }
    gsub("rv_p_oqp","rv_a_oqp")
  }
}
{print}
