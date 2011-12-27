BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}
/^ .*real\(8\),pointer :: rv_p_opp\(:\) =>NULL/{
 print comment,$0
next
}
/^ .*allocate\(rv_p_opp\(abs/{
  print comment,$0
 gsub("rv_p_opp","spot%rv_p_opp") 
}
/^ .*if \(-6\*nlspc<0\) rv_p_opp\(:\)=0.0d0/{
  print comment,$0
  gsub("rv_p_opp","spot%rv_p_opp") 
}
/^ .*spot%rv_p_opp => rv_p_opp/{
  print comment,$0
  next
}
/^ .*real\(8\),pointer :: rv_p_opp\(:\) =>NULL/{
  print comment,$0
  next
}
/^ .*rv_p_opp => spot%rv_p_opp/{
  print comment,$0
  next
}
/^ .*call dpscop \( rv_p_opp , olpp/{
  print comment,$0
 gsub("rv_p_opp","spot%rv_p_opp") 
}
/^ .*call dpscop \( olpp , rv_p_opp , k/{
  print comment,$0
 gsub("rv_p_opp","spot%rv_p_opp") 
}
{print}


