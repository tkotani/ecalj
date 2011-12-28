BEGIN{
 str=strftime("%b.%d.%Y")
 comment="ckino " str ": "
}

/complex\(8\),pointer :: zv_p_osmpt2\(:\) =>NULL/{
print comment,$0
next
}
/allocate\(zv_p_osmpt2\(abs/{
 print comment,$0
 gsub("zv_p_osmpt2","spot2%zv_p_osmpot")
}
/if \(-k1\*k2\*k3\*nsp<0) zv_p_osmpt2\(:\)=0.0d0/{
 print comment,$0
 gsub("zv_p_osmpt2","spot2%zv_p_osmpot")
}
/spot2%zv_p_osmpot => zv_p_osmpt2/{
 print comment,$0
 next
}
{print}
