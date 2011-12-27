
BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/^ .*sv_p_oorhat/{

 st="spot%"
 po="sv_p_oorhat"
 al="sv_a_oorhat"
 ty="type(s_rv1)"

 print comment,$0
 ty_sub=ty
 sub("\\(","\\(",ty_sub)
 sub("\\)","\\)",ty_sub)
 gsub(ty_sub " *, *pointer *::",ty " , allocatable ::")
 gsub("=>NULL\\(\\)","")
 gsub(po,al)
 if (match($0,"nullify")) { 
  print "        if (allocated(" st al ")) deallocate(" st al ")"
  next; 
 }
 if (match($0,"allocate\\(" st po)) {
  print "        if (allocated("st al")) deallocate(" st al ")"
 }

 gsub("associated\\(","allocated(")
}
{print}
