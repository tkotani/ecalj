BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "

}
/^ .*zv_p_osmrho/{

 if (match($0,"complex\\(8\\) *, *pointer *:: *zv_p_osmrho") ) { print comment,$0;next }
 if (match($0,"zv_p_osmrho *=> *spot%zv_p_osmrho")) { print comment,$0;next }
 if (FNR==437 || FNR==462||FNR==541||FNR==542||FNR==610||FNR==739||FNR==757) {
 print comment,$0
  gsub("zv_p_osmrho","spot%zv_p_osmrho")
 }

}
{print}
