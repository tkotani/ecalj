BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*iv_p_owk/{
  print comment,$0
  if (match($0,"pointer")) { next;}
  if (match($0,"iv_p_owk *=> *none")) { next;}
  if (match($0,"none *=> *iv_p_owk")) { next;}
  gsub("iv_p_owk","none")
  if (match($0,"none%none%")) {next;}
}
{print}
