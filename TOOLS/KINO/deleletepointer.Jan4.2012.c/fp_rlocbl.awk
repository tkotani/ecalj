BEGIN{
str=strftime("%b.%d.%Y");
comment="ckino " str ": "
}
/^ .*rv_p_oqpp *=> *sv_p_oqkkl\(1,ia\)%v/{
 print comment,$0
 print "#define OQPP sv_p_oqkkl\(1,ia\)%v"
 next
}
/^ .*rv_p_oqpp/{
 print comment,$0
 gsub("rv_p_oqpp","OQPP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oqhp *=> *sv_p_oqkkl\(2,ia\)%v/{
 print comment,$0
 print "#define OQHP sv_p_oqkkl\(2,ia\)%v"
 next
}
/^ .*rv_p_oqhp/{
 print comment,$0
 gsub("rv_p_oqhp","OQHP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oqhh *=> *sv_p_oqkkl\(3,ia\)%v/{
 print comment,$0
 print "#define OQHH sv_p_oqkkl\(3,ia\)%v"
 next
}
/^ .*rv_p_oqhh/{
 print comment,$0
 gsub("rv_p_oqhh","OQHH")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oeqpp *=> *sv_p_oeqkkl\(1,ia\)%v/{
 print comment,$0
 print "#define OEQPP sv_p_oeqkkl\(1,ia\)%v"
 next
}
/^ .*rv_p_oeqpp/{
 print comment,$0
 gsub("rv_p_oeqpp","OEQPP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oeqhp *=> *sv_p_oeqkkl\(2,ia\)%v/{
 print comment,$0
 print "#define OEQHP sv_p_oeqkkl\(2,ia\)%v"
 next
}
/^ .*rv_p_oeqhp/{
 print comment,$0
 gsub("rv_p_oeqhp","OEQHP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oeqhh *=> *sv_p_oeqkkl\(3,ia\)%v/{
 print comment,$0
 print "#define OEQHH sv_p_oeqkkl\(3,ia\)%v"
 next
}
/^ .*rv_p_oeqhh/{
 print comment,$0
 gsub("rv_p_oeqhh","OEQHH")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oppipp *=> *sv_p_oppi\(1,ia\)%v/{
 print comment,$0
 print "#define OPPIPP sv_p_oppi\(1,ia\)%v"
 next
}
/^ .*rv_p_oppipp/{
 print comment,$0
 gsub("rv_p_oppipp","OPPIPP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_osigpp *=> *sv_p_osig\(1,ia\)%v/{
 print comment,$0
 print "#define OSIGPP sv_p_osig\(1,ia\)%v"
 next
}
/^ .*rv_p_osigpp/{
 print comment,$0
 gsub("rv_p_osigpp","OSIGPP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oppihp *=> *sv_p_oppi\(2,ia\)%v/{
 print comment,$0
 print "#define OPPIHP sv_p_oppi\(2,ia\)%v"
 next
}
/^ .*rv_p_oppihp/{
 print comment,$0
 gsub("rv_p_oppihp","OPPIHP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_osighp *=> *sv_p_osig\(2,ia\)%v/{
 print comment,$0
 print "#define OSIGHP sv_p_osig\(2,ia\)%v"
 next
}
/^ .*rv_p_osighp/{
 print comment,$0
 gsub("rv_p_osighp","OSIGHP")
 if (match($0,"pointer")) {next}
}
{print}
