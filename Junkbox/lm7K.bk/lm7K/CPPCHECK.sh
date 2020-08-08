#!/bin/bash
#first version: 03/04/09, Hiori Kino


function test_FC {
afile=$1
if [ -e $afile ]; then
 rm -f $afile
fi

echo \
"       program test
        integer :: iarg
        write(*,*) iarg
        end" > $afile
$FC -o $exe  $afile 2> /dev/null
ret=$?
if [ -e $afile ]; then
rm -f $afile 
fi
if [ $ret != 0 ]; then
 echo "Error: Is \$FC correct? I can't compile/link a simple program."
 exit
fi
if [ -e $exe ]; then
rm -f $exe
fi 
}


 function find_argc { 

argc="NONE"

afile=$1 
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
 echo "check nargc"
 fi 
echo \
"       program test
        integer :: iarg,nargc
        iarg=nargc()
        end"   > $afile
$FC -o $exe  $afile  2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASNARGC "
   rm -f $afile
   return
fi 

if [ -e $afile ]; then
rm -f $afile
fi 
if [ $verb -gt 0 ]; then
echo "check iargc"
fi 
echo \
"       program test
        integer :: iarg,iargc
        iarg=iargc()
        end"   > $afile
$FC -o $exe  $afile  2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
  cppcmd=$cppcmd"-DHASIARGC "
  rm -f $afile
  return
fi

if [ -e $afile ]; then
rm -f $afile
fi
if [ $verb -gt 0 ]; then
echo "check command_argument_count"
fi
echo \
"       program test
        integer(4) :: iarg,command_argument_count
        iarg=command_argument_count()+1
        end"   > $afile
$FC -o $exe  $afile  2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
  cppcmd=$cppcmd"-DHASCOMMANDARGUMENTCOUNT "
  rm -f $afile
  return
fi

errmsg=$errmsg"Error: You must have subroutine/function to get number of command line strings.\n"

}



#-----------------------------------------------------------------

function find_argv {
argv="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check gtargc"
fi 
echo \
"       program test
        integer :: iarg
        character(80):: str
        iarg=1
        call gtargc(iarg,str)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASGTARGC "
   rm -f $afile
   return
fi

if [ -e $afile ]; then
rm -f $afile
fi
if [ $verb -gt 0 ]; then
echo "check getarg"
fi 
echo \
"       program test
        integer :: iarg
        character(80):: str
        iarg=1
        call getarg(iarg,str)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASGETARG "
   rm -f $afile
   return
fi

if [ -e $afile ]; then
rm -f $afile
fi
if [ $verb -gt 0 ]; then
echo "check get_command_argument"
fi
echo \
"       program test
        integer :: iarg
        character(80):: str
        iarg=1
        call get_command_argument(iarg,str)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASGETCOMMANDARGUMENT "
   rm -f $afile
   return
fi


errmsg=$errmsg"Error: You must have subroutine/function to get command line strings.\n"

}

#-----------------------------------------------------------------

function find_fdate {

fdate="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check fctime"
fi 
echo \
"       program test
        character(80):: str
        call fctime(str)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DFCTIME "
   rm -f $afile
   return
fi


if [ -e $afile ]; then
 rm -f $afile
fi
if [ $verb -gt 0 ]; then
echo "check fdate"
fi 
echo \
"       program test
        character(80):: str
        call fdate(str)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DFDATE "
   rm -f $afile
   return
fi

if [ -e $afile ]; then
 rm -f $afile
fi
if [ $verb -gt 0 ]; then
echo "check \$fdate"
fi 
echo \
"       program test
        character(80):: str
        call \$fdate(str)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DDOL_FDATE "
   rm -f $afile
   return
fi

errmsg=$errmsg"Warning: better to have subroutine/function to get date.\n"

}

function find_time {

time="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check gettimeofday"
fi
echo \
"       program test
        integer:: ierr
        integer(4):: v(2)
        call gettimeofday(v,ierr)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASGETTIMEOFDAY "
   rm -f $afile
   return
fi


if [ -e $afile ]; then
 rm -f $afile
fi
if [ $verb -gt 0 ]; then
echo "check cpu_time"
fi
echo \
"       program test
        real:: x
        call cpu_time(x)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASCPUTIME "
   rm -f $afile
   return
fi

errmsg=$errmsg"Warning: better to have subroutine/function to get time.\n"

}


function find_getenv {

getenv="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check getenvqq"
fi
echo \
"       program test
        character(40)::pnam,pval
        call getenvqq(pnam,pval)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASGETENVQQ "
   rm -f $afile
   return
fi

if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check get_environment_variable"
fi
echo \
"       program test
        character(40)::pnam,pval
        call get_environment_variable(pnam,pval)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASGETENVIRONMENTVARIABLE "
   rm -f $afile
   return
fi

if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check getenv"
fi
echo \
"       program test
        character(40)::pnam,pval
        call getenv(pnam,pval)
        end"   > $afile
#$FC -o $exe  $afile -lf90c 2>/dev/null
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASGETENV "
   rm -f $afile
   return
fi



errmsg=$errmsg"Error: You must have subroutine/function to get environment.\n"

}


function find_setenv {

setenv="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check setenvqq"
fi
echo \
"       program test
        character(40)::pnam,pval
        call setenvqq(pnam,pval)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASSETENVQQ "
   rm -f $afile
   return
fi

}


function find_ifport {

ifport="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check ifport"
fi
echo \
"       program test
        use ifport
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASIFPORT "
   rm -f $afile
   return
fi

}


function find_etime {


afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check ifport"
fi
echo \
"       program test
        real t(2)
        call etime(t)
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASETIME "
   rm -f $afile
   return
fi


if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check ifport"
fi
echo \
"       program test
        real t(2)
        t(1)= second()
        end"   > $afile
$FC -o $exe  $afile -lf90c 2>/dev/null
ret=$?
#echo $ret
if [ $ret  == 0 ]; then
   cppcmd=$cppcmd"-DHASSECOND "
   libcmd=$libcmd"-lf90c "

   rm -f $afile
   return
fi


}



function find_noquad {

noquad="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check quad"
fi
echo \
"       program test
        real*16 x
        end"   > $afile
$FC -o $exe  $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  != 0 ]; then
   cppcmd=$cppcmd"-DNOQUAD "
   rm -f $afile
   return
fi

}

function find_ovar {

ovar="NONE"

afile=$1
exefile=$2
if [ -e $afile ]; then
 rm -f $afile
fi

if [ $verb -gt 0 ]; then
echo "check overlap varriable"
fi
echo \
"       REAL FUNCTION R1MACH(I)
        INTEGER SMALL(2)
       INTEGER LARGE(2)
       INTEGER RIGHT(2)
       INTEGER DIVER(2)
       INTEGER LOG10(2)
       REAL RMACH(5)
       EQUIVALENCE (RMACH(1),SMALL(1))
       EQUIVALENCE (RMACH(2),LARGE(1))
       EQUIVALENCE (RMACH(3),RIGHT(1))
       EQUIVALENCE (RMACH(4),DIVER(1))
       EQUIVALENCE (RMACH(5),LOG10(1))
       DATA SMALL(1) /     8388608 /
       DATA LARGE(1) /  2139095039 /
       DATA RIGHT(1) /   864026624 /
       DATA DIVER(1) /   872415232 /
       DATA LOG10(1) /  1050288283 /
       r1mach=l
        end"   > $afile
$FC -o $exe -c $afile 2>/dev/null
ret=$?
#echo $ret
if [ $ret  != 0 ]; then
   cppcmd=$cppcmd"-DNOT_OVERLAP_VAR "
   rm -f $afile
   return
fi

}



#function printcpp {
## show result
#while [ "$#" -gt 0 ]
#do
#if [ $1 != "NONE" ]; then
#echo -n "-D"$1 " "
#fi
#shift
#done
#}


#-----------------------------------------------------------------
#main start
if [ $# -ne 1 ]; then
  echo " usage: CPPCHECK.sh {name of your fortran comilar}"
  echo "       e.g,  CPPCHECK.sh ifort"
  echo "   This script judges CPP_SW (non MPI mode)."
  echo "   You can use these for CPP_SW in Make.inc."
  echo "   Let us know if this did not work in your machine."
  exit 1
fi
FC=$1
#echo $FC

# write debug info if verb=1 
# if not, verb=0
verb=0

#FC is necessary
if [ -z $FC ]; then
echo You must set FC to run this script
echo "E.g., export FC=ifort; bash this_script "
exit
fi 

#output file=$afile
afile="__ff__.f"
exe="a___.out"
cppcmd=""
errmsg=""
libcmd=""

#function call

#try compiling, to check FC is correct or not
test_FC $afile


find_argc  $afile $exe
find_argv $afile $exe
find_fdate $afile $exe
find_time $afile $exe
find_getenv $afile $exe
find_setenv $afile $exe
find_ifport $afile $exe
find_noquad $afile $exe
find_ovar $afile $exe
find_etime $afile $exe


#cleanup $exe
if [ -e $exe ]; then
rm -f $exe
fi 
if [ -e $afile ]; then
rm -f $afile
fi 

echo -n '       :' CPP_SW=
echo -n $cppcmd
echo ""

echo -n ForMPI : CPP_SW=
echo -n $cppcmd
echo -n " -UMPE -UMPIK -DMPI"
echo ""

echo -n ForMPIK: CPP_SW=
echo -n $cppcmd
echo -n " -UMPE -DMPIK -UMPI" 
echo ""

# error 
if [ -n "$errmsg" ]; then
echo -e $errmsg
fi

# lib
if [ -n "$libcmd" ]; then
echo  "Add library: "$libcmd
fi 

