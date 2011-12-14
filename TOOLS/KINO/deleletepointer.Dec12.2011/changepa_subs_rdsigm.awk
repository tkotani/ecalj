BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_oiaxs/{
print comment,$0
sub("integer,pointer","integer,allocatable")
gsub("iv_p_oiaxs","iv_a_oiaxs")
sub("=> *NULL\\( *\\)","")
}
/^ .*iv_p_ontabs/{
print comment,$0
sub("integer,pointer","integer,allocatable")
gsub("iv_p_ontabs","iv_a_ontabs")
sub("=> *NULL\\( *\\)","")
}
/^ *sham%iv_a_ontabs *=> *iv_a_ontabs/{
  print comment,$0
  gsub("=>"," => ")
  n=split($0,a)
  print "      call move_alloc(from=",a[3],",to=",a[1],")"
  next
}
/^ *sham%iv_a_oiaxs *=> *iv_a_oiaxs/{
  print comment,$0
  gsub("=>"," => ")
  n=split($0,a)
  print "      call move_alloc(from=",a[3],",to=",a[1],")"
  next
} 
{print }
