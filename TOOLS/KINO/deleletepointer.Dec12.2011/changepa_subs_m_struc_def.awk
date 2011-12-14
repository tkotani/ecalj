BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
 sham=0
}
/^ *type /{
 if ($2=="s_ham") { sham=1 }
 else { sham=0}
}
/^ *integer, pointer :: iv_p_oiaxs/{
 if (sham==1) {
  print comment,$0
  sub("pointer","allocatable")
  sub("iv_p_oiaxs","iv_a_oiaxs")
  sub("=>NULL\\(\\)","")
 }
}
/^ *integer, pointer :: iv_p_ontabs/{
 if (sham==1) {
  print comment,$0
  sub("pointer","allocatable")
  sub("iv_p_ontabs","iv_a_ontabs")
  sub("=>NULL\\(\\)","")
 }
}
{print}
