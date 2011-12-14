BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}


/real\(8\),pointer :: rv_p_ov0/{
 print comment,$0
 next
}
/rv_p_ov0 => ssite\(ib\)%rv_p_ov0/{
 print comment,$0
 next
}
/call makusp \( n0 , z , nsp , rmt , lmxa , rv_p_ov0/{
 print comment,$0
 sub("rv_p_ov0","ssite(ib)%rv_p_ov0")
}
{print }

