BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*integer,pointer :: iv_p_oidtet/{
print comment,$0
next
}
/^ .*iv_p_oidtet => sbz%iv_p_oidtet/{
print comment,$0
next
}
/^ .*\.     , 2 , ntet , iv_p_oidtet , sev , dum/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*eferm , 2 , ntet , iv_p_oidtet , sev , dum/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*.         \( 2 \) , nkabc \( 3 \) , nkp , ntet , iv_p_oidtet , qval - qbg/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*if \(.not.associated\(iv_p_oidtet\)\)then/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*allocate\(iv_p_oidtet\(1\)\)/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*         \( 2 \) , nkabc \( 3 \) , nkp , ntet , iv_p_oidtet , qval - qbg/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*deallocate\(iv_p_oidtet\); iv_p_oidtet=>NULL\(\)/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}
/^ .*.       2 \) , dos_rv , ndos , eferm , 1 , ntet , iv_p_oidtet , dum/{
print comment,$0
gsub("iv_p_oidtet","sbz%iv_p_oidtet")
}

{print}
