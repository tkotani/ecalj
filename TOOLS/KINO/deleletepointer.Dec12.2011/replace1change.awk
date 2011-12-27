BEGIN{
  print "#please set LANG=C before invoking this script to show date in English"
  print "#input= ",ARGV[1]
  print "BEGIN{"
  print " str=strftime(\"%b.%d.%Y\")"
  print " comment=\"ckino \" str \": \"";
  print "}"
}
{
 if (NF==0) {next}
  str=$0
  gsub("^ *","",str)
  gsub(" *$","",str)
  strmatch=str
  gsub("\\(","\\\\(",strmatch)
  gsub("\\)","\\\\)",strmatch)

  if (match($0,"=>")) {
  sub("=>","")
  n=split($0,a)
  }
  if (match($0,"gsub")){
  gsub("gsub"," ")
  gsub("\""," ")
  gsub(","," ")
  gsub("\\("," ")
  gsub("\\)"," ")
  n=split($0,a)
  }
  name1=a[1]
  name2=a[2]
  n=split(name2,a,"%")

  name1match=name1
  gsub("\\(","\\\\(",name1match)
  gsub("\\)","\\\\)",name1match)
  name2match=name2
  gsub("\\(","\\\\(",name2match)
  gsub("\\)","\\\\)",name2match)
  struc= a[1] "%"
  strucmatch=struc
  gsub("\\(","\\\\(",strucmatch)
  gsub("\\)","\\\\)",strucmatch)

  print "/^ .*" name1 "/{"
  print "  old=$0"
#  print "  if (match($0,\"" strmatch "\")) {next}"
  print "  if (match($0,\"pointer\")) { print comment,$0;next;}"
  print "  if (match($0,\"" name1match " *=> *" name2match "\")) { print comment,$0;next;}"
  print "  if (match($0,\"" name2match " *=> *" name1match "\")) { print comment,$0;next;}"
  print "  gsub(\"" name1 "\",\"" name2 "\")"
#  print "  gsubexact0($0,\"" name1 "\",\"" name2 "\")"
#  print "  if (match($0,\""  strucmatch strucmatch "\")) {next;}"
  print "  gsub(\"" strucmatch strucmatch "\",\"" strucmatch "\")"
  print "  if (old != $0 ) {print comment,$0}"
  print "}"
}
END{ 
print "{print}"
}
