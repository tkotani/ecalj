# usage
# gawk -f thisscript gwsc > gwsc.ref 
#
BEGIN{
COMMENT="###make_ref.awk"
}
/intelMPI.csh)/ {
print COMMENT,$0
print "%HEADER"
next
}
/^[^#].*\$nfpgw\//{ 
 print COMMENT,$0
 echoinput=""
 output=""
 program=""
 ntarget=0
 mpi=0
if (match($0,"echo.*\\|")>0 ) {
# print "match",$0
 i = match($0,"\\|")
 s1=substr($0,1,i-1)
 echoinput=extract_echo(s1)
 s2=substr($0,i+1)
}
else{
 s2=$0
}

 output=extract_output(s2)
 program=extract_program(s2)
 extract_target(s2)
 if (match(s2,"mpirun")) { mpi=1 }


 print_cmd()
 next
}
{
 print
}

func print_cmd() {

 printf("%CMD " )
 if (mpi>0) { printf("mpi=1 ") }
 if (length(program)>0) { printf("program=%s ",program) }
 if (length(echoinput)>0) { printf("echoinput=%s ",echoinput) }
 if (ntarget>0) { 
   printf("target=%s ",target[1]) 
   for (i=2;i<=ntarget;i++) {
     printf("target%i=%s ",i,target[i])
   }
  }
 if (length(output)>0) { printf("output=%s ",output) }
 print ""

}

func extract_echo(s1){
  sub("echo *","",s1)
  sub("^ *","",s1)
  return s1
}

func extract_output(s2) {
  sub("\.*>","",s2)
  sub("^ *","",s2)
  return s2
}

func extract_program(s2,  n,a) {
  sub("mpirun *-np *\\$MPI_SIZE","",s2)
  sub("\>.*$","",s2)
  sub("\\$TARGET","",s2)
  sub("^ *","",s2)
  sub("\\$nfpgw\/","",s2)
  n=split(s2,a)
  return a[1]
}


func extract_target(s2) {
  sub("^.*\\$nfpgw\/[^ ]+","",s2)
  sub("\>.*$","",s2)
  sub("^ *","",s2)
  ntarget=split(s2,target)
  return ntarget
}
