BEGIN{
spc0="       "
}
/^ *if.*.*deallocate *[(]/{
 print "C",$0
 i=index($0,"deallocate")
 print substr($0,1,i-1)" then"
 print spc0,substr($0,i)
 print spc0,"endif"
 next
}
/^ *if.*allocate *[(]/{
 print "C",$0
 i=match($0,"allocate *[(]")
 print substr($0,1,i-1)" then"
 print spc0,substr($0,i)
 print spc0,"endif"
 next
}

{print }
