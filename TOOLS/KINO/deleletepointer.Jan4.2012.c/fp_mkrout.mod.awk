BEGIN{
str=strftime("%b.%d.%Y");
comment="ckino " str ": "
}
/^ .*rv_p_oqhh *=> *sv_p_oqkkl\(3,ib\)%v/{
 print comment,$0
 print "#define OQHH sv_p_oqkkl\(3,ib\)%v"
 next
}
/^ .*rv_p_oqhh *=> *sv_p_oeqkkl\(3,ib\)%v/{
 print comment,$0
 print "#define OQHH sv_p_oeqkkl\(3,ib\)%v"
 next
}
/^ .*rv_p_oqhh/{
 print comment,$0
 gsub("rv_p_oqhh","OQHH")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oqhp *=> *sv_p_oqkkl\(2,ib\)%v/{
 print comment,$0
 print "#define OQHP sv_p_oqkkl\(2,ib\)%v"
 next
}
/^ .*rv_p_oqhp *=> *sv_p_oeqkkl\(2,ib\)%v/{
 print comment,$0
 print "#define OQHP sv_p_oeqkkl\(2,ib\)%v"
 next
}
/^ .*rv_p_oqhp/{
 print comment,$0
 gsub("rv_p_oqhp","OQHP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oqpp *=> *sv_p_oqkkl\(1,ib\)%v/{
 print comment,$0
 print "#define OQPP sv_p_oqkkl\(1,ib\)%v"
 next
}
/^ .*rv_p_oqpp *=> *sv_p_oeqkkl\(1,ib\)%v/{
 print comment,$0
 print "#define OQPP sv_p_oeqkkl\(1,ib\)%v"
 next
}
/^ .*rv_p_oqpp/{
 print comment,$0
 gsub("rv_p_oqpp","OQPP")
 if (match($0,"pointer")) {next}
}
{print}
