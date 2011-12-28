
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*iv_p_oidtet/{

 st="sbz%"
 po="iv_p_oidtet"
 al="iv_a_oidtet"
 ty="integer"

 print comment,$0
 ty_sub=ty
 sub("\\(","\\(",ty_sub)
 sub("\\)","\\)",ty_sub)
 gsub(ty_sub " *, *pointer *::",ty " , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub(po,al)
 if (match($0,"nullify *\\( *" st al)) { 
  print "        if (allocated(" st al ")) deallocate(" st al ")"
  next; 
 }
 if (match($0,"allocate *\\( *" st po)) {
  print "        if (allocated("st al")) deallocate(" st al ")"
 }

 gsub("associated *\\( *" st al,"allocated(" st al)
}
{print}
