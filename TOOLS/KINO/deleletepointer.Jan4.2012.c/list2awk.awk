BEGIN{
 print "BEGIN{"
 print "str=strftime(\"%b.%d.%Y\");"
 print "comment=\"ckino \" str \": \"";
 print "}"

}
/=>/{
  gsub("=>"," => ")
  n=split($0,a)
  file1=a[1]
  file2=a[3]

  file1u=file1
  sub("rv_p_","",file1u)
  file1u=toupper(file1u)

  gsub("\\(","\\(",file2)
  gsub("\\)","\\)",file2)

  print "/^ .*"file1" *=> *"file2"/{"
  print " print comment,$0"
  print " print \"#define "file1u " "file2"\""
  print " next"
  print "}"
  print "/^ .*"file1"/{"
  print " print comment,$0"
  print " gsub(\""file1"\",\""file1u"\")"
  print " if (match($0,\"pointer\")) {next}"
  print "}"
}
END{
  print "{print}"
}
