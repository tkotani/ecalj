BEGIN{

 dir="tmp/"
 file="x"
 cmd="ls " dir "fix*.awk > " file
 system(cmd)
 i=0
 while (1) {
   ret=getline <file
   if (ret<=0) { break; }
   i++
   sub(dir,"",$1)
   script[i]=$1
 }
 close(file)

 nscript=i
 for (i=1;i<=nscript;i++) {
   print "make ", script[i]
   do_script(script[i])
 }
}

func do_script(scriptfile,  var,file,cmd,i,ret,nffile,runfile) {

  var=scriptfile
  sub(dir,"",var)
  sub("fix_","",var)
  sub("\\.awk","",var)
  file="x"
  cmd="grep -i \"^ .*" var "\" */*.F *.F -l > " file
  system(cmd)
  i=0
  while (1) {
    ret=getline <file
    if (ret<=0) {break}
    i++
    ffile[i]=$1
  }
  close(file)

  nffile=i
  runfile=scriptfile
  sub("\\.awk","",runfile)
  runfile = dir "run_" runfile
  
  print dir
  print scriptfile
  print runfile
  print "script=" dir scriptfile >runfile 
  printf "for n in " >runfile
  for (i=1;i<=nffile;i++){
   printf ffile[i] " " >runfile
  }
  print "" >runfile
  print "do gawk -f $script $n >x;mv x $n; done" > runfile
  close(runfile)
}



