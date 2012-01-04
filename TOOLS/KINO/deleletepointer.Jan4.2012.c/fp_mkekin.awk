BEGIN{
str=strftime("%b.%d.%Y");
comment="ckino " str ": "
}
/^ .*rv_p_oqhh *=> *sv_p_oqkkl\(3,ib\)%v/{
 print comment,$0
 print "#define OQHH sv_p_oqkkl\(3,ib\)%v"
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
/^ .*rv_p_oqpp/{
 print comment,$0
 gsub("rv_p_oqpp","OQPP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_otauhh *=> *sv_p_otau\(3,ib\)%v/{
 print comment,$0
 print "#define OTAUHH sv_p_otau\(3,ib\)%v"
 next
}
/^ .*rv_p_otauhh/{
 print comment,$0
 gsub("rv_p_otauhh","OTAUHH")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_otauhp *=> *sv_p_otau\(2,ib\)%v/{
 print comment,$0
 print "#define OTAUHP sv_p_otau\(2,ib\)%v"
 next
}
/^ .*rv_p_otauhp/{
 print comment,$0
 gsub("rv_p_otauhp","OTAUHP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_otaupp *=> *sv_p_otau\(1,ib\)%v/{
 print comment,$0
 print "#define OTAUPP sv_p_otau\(1,ib\)%v"
 next
}
/^ .*rv_p_otaupp/{
 print comment,$0
 gsub("rv_p_otaupp","OTAUPP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_osighh *=> *sv_p_osig\(3,ib\)%v/{
 print comment,$0
 print "#define OSIGHH sv_p_osig\(3,ib\)%v"
 next
}
/^ .*rv_p_osighh/{
 print comment,$0
 gsub("rv_p_osighh","OSIGHH")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_osighp *=> *sv_p_osig\(2,ib\)%v/{
 print comment,$0
 print "#define OSIGHP sv_p_osig\(2,ib\)%v"
 next
}
/^ .*rv_p_osighp/{
 print comment,$0
 gsub("rv_p_osighp","OSIGHP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_osigpp *=> *sv_p_osig\(1,ib\)%v/{
 print comment,$0
 print "#define OSIGPP sv_p_osig\(1,ib\)%v"
 next
}
/^ .*rv_p_osigpp/{
 print comment,$0
 gsub("rv_p_osigpp","OSIGPP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oppihh *=> *sv_p_oppi\(3,ib\)%v/{
 print comment,$0
 print "#define OPPIHH sv_p_oppi\(3,ib\)%v"
 next
}
/^ .*rv_p_oppihh/{
 print comment,$0
 gsub("rv_p_oppihh","OPPIHH")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oppihp *=> *sv_p_oppi\(2,ib\)%v/{
 print comment,$0
 print "#define OPPIHP sv_p_oppi\(2,ib\)%v"
 next
}
/^ .*rv_p_oppihp/{
 print comment,$0
 gsub("rv_p_oppihp","OPPIHP")
 if (match($0,"pointer")) {next}
}
/^ .*rv_p_oppipp *=> *sv_p_oppi\(1,ib\)%v/{
 print comment,$0
 print "#define OPPIPP sv_p_oppi\(1,ib\)%v"
 next
}
/^ .*rv_p_oppipp/{
 print comment,$0
 gsub("rv_p_oppipp","OPPIPP")
 if (match($0,"pointer")) {next}
}
{print}
