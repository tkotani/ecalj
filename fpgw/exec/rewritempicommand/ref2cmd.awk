# usage
# gawk -f thisscript gwsc.ref > gwsc.new
#
BEGIN {
 COMMENT="###ref2cmd.awk"

 tempfile="_IN_"    # temporary file name 
 mpirun="mpirun"    # MPI command 
 programdir="$nfpgw/"   # programpath 
# header="if(-e /home/etc/intelMPI.csh) source /home/etc/intelMPI.csh"   # custom extra script for some system, e.g. Kyushu-U, tatara
 header=""   # if uncessary blank it.

# and also change  make_cmd()

}
/%HEADER/{
print COMMENT,$0
 if (length(header)>0) {
   print header
 }
next
}
/%CMD/{
 print COMMENT,$0
 reset_item()
 sub("#.*$","")
 for (i=1;i<=NF;i++) {
   extract_item($i)
 }
 make_cmd() 
  next
}
{
print
}

func make_cmd( ){
# print "->",mpi,program,echoinput,target,output
  print "#>>>"
  if (length(echoinput)>0) { printf("echo %s > %s\n", echoinput,tempfile)}
  if (mpi>0) { printf("%s -np $MPI_SIZE ",mpirun) } 

  if (program=="hx0fp0_sc") { runoption="-nq 2 -nm 2" } # for -np 4 , q=2,m=2 parallel  
  else { runoption="" }

  if (length(program)>0) { printf("%s %s ",programdir program,runoption); }
  if (length(target)>0) { printf("%s ",target); }
  if (length(target2)>0) { printf("%s ",target2); }
  if (length(target3)>0) { printf("%s ",target3); }
  if (length(echoinput)>0) { printf("< %s ",tempfile) }
  if (length(output)>0) { printf("> %s ",output); }
  printf("\n")

 # postprocessing
 print "  if ( $? != 0 )  then"
 printf( "    echo Error in \'%s \' ",program);
 if (length(echoinput)>0) {
 printf( "input=\'%s\' ",echoinput);
 }
 print "output=\'"output"\'"
 print "    exit 10"
 print "  endif"
  print "#<<<"


}

func reset_item(){
 program=""
 echoinput=""
 output=""
 target=""
 target2=""
 target3=""
 mpi=0
}

func extract_item(s, a,n) {
 if (s=="%CMD") {return }
 n=split(s,a,"=")
 if (n!=2) { print "ERROR",s,n; exit(10)}
# print a[1],a[2]
 if (a[1]=="program"){ program=a[2] }
 if (a[1]=="echoinput"){ echoinput=a[2] }
 if (a[1]=="output"){ output=a[2] }
 if (a[1]=="target"){ target=a[2] }
 if (a[1]=="target2"){ target2=a[2] }
 if (a[1]=="target3"){ target3=a[2] }
 if (a[1]=="mpi"){ mpi=a[2] }
}

