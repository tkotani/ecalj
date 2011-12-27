BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*real\(8\),pointer :: rv_p_oqp\(:\) =>NULL\(\)/{
print 
print comment
gsub("pointer","allocatable,target")
gsub("=>NULL\\(\\)","")
gsub("rv_p_oqp","rv_a_oqp")
print
next
}
/^ .*allocate\(rv_p_oqp\(abs/{
str=$0
sub("allocate.*$","",str)
print comment,$0
print str "if (allocated(rv_a_oqp)) deallocate(rv_a_oqp)"
gsub("rv_p_oqp","rv_a_oqp")
print
print str "rv_p_oqp=>rv_a_oqp"
next
}
/^ .*if \(associated\(rv_p_oqp\)\) deallocate\(rv_p_oqp/{
print comment,$0
str=$0
sub("if.*","",str)
print str "nullify(rv_p_oqp)"
next 
}
{print}
