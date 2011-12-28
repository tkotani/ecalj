BEGIN{
  s="zv_p_osmrh2"
  t="spot2%zv_p_osmrho"

 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "


}
/^ .*zv_p_osmrh2/{
 print comment,$0
 if (match($0,"complex\\(8\\),pointer *:: *" s)) { next }
 if (match($0,"=>")) { next; }

 if (match($0,"allocate *\\( *"s)) { 
    spc=$0
    sub("allocate.*$","",spc)
    print spc "if (associated(" t "))  deallocate(" t ")"
 }

 gsub("zv_p_osmrh2","spot2%zv_p_osmrho")
}
{print }
