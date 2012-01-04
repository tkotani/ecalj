BEGIN{
str=strftime("%b.%d.%Y");
comment="ckino " str ": "
}

/^ .*v_sspec/{
 print comment,$0
 sub("pointer","allocatable")
 sub("=> *NULL *\\( *\\)","")
}
/^ .*v_ssite/{
 print comment,$0
 sub("pointer","allocatable")
 sub("=> *NULL *\\( *\\)","")
}
{print}
