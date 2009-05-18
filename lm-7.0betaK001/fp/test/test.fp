#!/bin/csh -f

# A shell script testing operation of fp suite
# set verbose

alias call 'set retcall = \!\!:2 ; set callarg = \!\!:3 ; goto \!\!:1'
alias runjob 'set retcall = \!\!:1; set outfile = \!\!:2 ; set callarg = \!\!:3 ; goto runjob'
alias runrdcmd 'set retcall = \!\!:1; set rdcmdfmt = \!\!:2 ; set outfile = \!\!:3 ; set callarg = \!\!:4 ; goto runrdcmd'
alias findcmd  'set retcall = \!\!:1 ; set prog_cmd = \!\!:2 ; set path_name = \!\!:3 ; set make_path = \!\!:4 ; goto findcmd'
alias extract_res_n 'set retcall = \!\!:1; set testvar = \!\!:2 ; set refvar = \!\!:3 ; set keyword = \!\!:4  ; set arg_number = \!\!:5 ; set occur_number = \!\!:6 ; set sed_strn = \!\!:7 ; goto extract_res_n'
alias compare_res 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set refvar = \!\!:4 ; set tol = \!\!:5 ; set passvar = \!\!:6 ; goto compare_res'
alias compare_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set tol = \!\!:4 ; set passvar = \!\!:5 ; goto compare_res_0'
alias compare_resf 'set retcall = \!\!:1; set testvar = \!\!:2 ; set refvar = \!\!:3 ; set keyword = \!\!:4  ; set arg_number = \!\!:5 ; set occur_number = \!\!:6 ; set sed_strn = \!\!:7 ; goto compare_resf'
#alias zcmpmfiles_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set tol = \!\!:4 ; set passvar = \!\!:5 ; set ndig = \!\!:6 ; set srcfile = \!\!:7 ; set reffile = \!\!:8 ; goto zcmpmfiles_res_0 '
alias zcmpmfiles_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set tol = \!\!:3 ; set passvar = \!\!:4 ; set ndig = \!\!:5 ; set srcfile = \!\!:6 ; set reffile = \!\!:7 ; goto zcmpmfiles_res_0 '
alias cnvt_d_fmt  'set retcall = \!\!:1; set testvar = \!\!:2 ; set testval = \!\!:3 ; goto cnvt_d_fmt'
alias query 'set retcall = \!\!:1 ; set retcall2 = \!\!:2 ; set callarg = \!\!:3 ; goto query'

set allargs = ($argv)
set a
set slow
set testfile = $0
set testdir = $testfile:h
#set topdir  = `$testdir/../../startup/absolute-path $testdir/../..`
set topdir  = `$testdir/../../TOOLS/absolute-path  $testdir/../..`
set maindir = $topdir/main
set space = '        '
set failed = 0
alias zcat 'gunzip -c'

# Prepend current working-directory, top-level dir and maindir to path
set path = ($cwd $topdir $maindir $path)

set plot = `which fplot`
if (-x "$plot") then
  if `$plot --h | head -1 | awk '{print ($3 == "fplot")}'` set have_fplot
endif
set mc = `which mc`
if (-x "$mc") then
  if `$mc --h |& head -1 | awk '{print ($7 == "(vsn" && ($8 * 1 >= 1.04))}'` set have_mc
endif
set pldos = `which pldos`
if (-x "$pldos") then
  if `$pldos --h | head -1 | awk '{print ($2 == "pldos")}'` set have_pldos
endif
# see if ghostscript is available
set gs = `which gs`
if (-x "$gs") then
  if `$gs --help | head -1 | awk '{print ($2 == "Ghostscript")}'` set have_ghostscript
endif
# see if gnu grep is available
echo X | grep -A 1 X > & /dev/null
set retval = $status
if ($retval == 0) set gnu_grep

# --- Pick off switches ---
while (`echo $1 | sed -e 's/\(.\).*/\1/' `  ==  "-")

  set arg1 = $1; shift
  if ($?verb) echo test.fp: parsing switch $arg1
  switch ($arg1)
    case "--quiet":
      set quiet
      unset slow
      breaksw
    case "--add0":
      set ladd0
      breaksw
    case "--clean":
      set clean
      breaksw
    case "--veryclean":
      set clean
      set veryclean
      breaksw
    case "--no-iact*":
      unset slow
      breaksw
    case "--verb*":
      set verb = 1
      breaksw

    case "--all":
      set mater_lst = (copt te zrt co cr3si6 fe cu srtio3 felz gasls gdn eras c crn)
      set joblist
      while (`echo $1 | sed -e 's/\([0-9][0-9]*\)/-/'`  ==  "-")
        set joblist = ($joblist $1)
        shift
      end
      set pass
      set failed
      foreach i ($mater_lst)
        $testfile `echo $allargs | sed s/--all//g | sed -e 's/\([0-9][0-9]*\)//g' | sed -e 's/-add/-add0/g'` $i $joblist
        set retval = $status
        if ($retval != 0) then
          unset pass
          set failed = ($failed $i)
#  	  echo " $testfile : failed test $i ... aborting"
#            exit -1
        endif
      end
      if ($?pass) then
        echo "$space all tests PASSED ($mater_lst)"
        exit
      else
        echo "$space checks FAILED for the following materials:  $failed"
        exit -1
      endif

    default:
      echo unrecognized switch $arg1
      goto usage
  endsw

end

echo ' '
echo "         ---- test.fp: test FP program lmf ---"

# --- use copt as default in the absence of specific choice ---
if ($#argv == 0) then
  set ext = copt
  echo "$space .... no file extension specified; use input file ctrl.copt"
else
  set ext = $argv[1]
  shift
endif

if (! -e $testdir/ctrl.$ext) then
   echo ' '
   echo " test.fp aborting ... missing file $testdir/ctrl.$ext"
   goto usage
endif

if ($ext == "copt") then
  echo '         Case copt: a distorted L12 environment with four atoms.'
#    echo '         Other checks:'
#    echo '         spin-pol, tetrahedron+metal=3, forces(mode 12),'
#    echo '         2-kappa basis, inequiv kmxa,lmax,rmt'
  set cplst = ($testdir/{ctrl.copt,spec.prop})
# set fitbas  = ' Co RSMH= 2.425 2.717 1.017 EH= -0.360 -0.200 -0.222'
# set fitbas2 = ' Pt RSMH= 2.461 3.042 1.085 EH= -0.441 -0.200 -0.200'
  set dfmax1tol1 = 0.1
  set dfmaxntol1 = 0.1
else if ($ext == "fe") then
  echo '         Case Fe: spin-polarized Fe in bcc structure'
  set cplst = ($testdir/{ctrl.fe})
  set dosmulltol = 7e-3
  set dosclstol = .0015
else if ($ext == "felz") then
  echo '         Case felz: spin-polarized Fe spin-orbit coupling'
  set cplst = ($testdir/lzsz/{ctrl,rsta}.felz)
  set dorbmtol = 0.00001
else if ($ext == "gasls") then
  echo '         Case GaAs: GaAs with spin-orbit coupling'
  set cplst = ($testdir/ls/{ctrl,rsta,syml}.gasls)
  set gmtol = 0.0001
# lineeval specifies which line evals of interest are in AFTER bndfp: tag
# eval1 eval2 specifies ranges of columns evals are in
  set lineeval=3  eval1=4  eval2=9  evalso=5
else if ($ext == "gaslc") then
  echo '         Case GaAs: GaAs with spin-orbit coupling, local orbitals and scGW self-energy'
  set cplst = ($testdir/ls/{ctrl.gaslc,rst.gaslc,syml.gaslc,sigm.gaslc,semi.mater})
  set gmtol = 0.0001
# lineeval specifies which line evals of interest are in AFTER bndfp: tag
# eval1 eval2 specifies ranges of columns evals are in
  set lineeval=4  eval1=5  eval2=9  evalso=6
else if ($ext == "gdn") then
  echo '         Case GdN: Test of LDA+U, and also LDA with spin polarized 4f core'
  set cplst = ($testdir/{ctrl.gdn,occnum.gdn})
  set gmtol = 0.0001
  set dosmulltol = 2e-4
else if ($ext == "cdte") then
  echo '         Case CdTe: Test of LDA+U in different modes'
  set cplst = ($testdir/{ctrl.cdte,occnum.cdte,semi.mater})
else if ($ext == "eras") then
  echo '         Case ErAs: Test of LDA+U'
  set cplst = ($testdir/{ctrl.eras,occnum.eras})
  set gmtol = 0.0001
else if ($ext == "er") then
  echo '         Case Er: Test of LDA+U'
  set cplst = ($testdir/{ctrl.er,site.er,occnum.er,syml.er,specialspec1,atparms})
  set gmtol = 0.0001
else if ($ext == "zrt") then
  echo '         Case zrt: ZrO_2 fluorite in tetragonal setting, with tetragonal distortion'
  set cplst = ($testdir/{ctrl.zrt})
  set dfmax1tol1 = 0.01
  set dehf1toln = 6e-6
else if ($ext == "co") then
  echo '         Case co: a hexagonal environment with two equivalent atoms.'
#    echo '         Other checks:'
#    echo '         spin-polarized, tetrahedron+metal=4'
  set cplst = ($testdir/{ctrl.co,syml.co})
  set lmfargs1 = "-vmet=4 -vlmf=1 -vnk=8 -vnit=10 --pr31"
# set fitbas = ' A RSMH= 2.375 2.722 1.047 EH= -0.342 -0.200 -0.236'
  set drmsqtol1  = 1e-6
  set pdostol    = 0.01
else if ($ext == "cr3si6") then
  echo '         Case cr3si6: a hexagonal environment with several atoms, two classes'
  echo '         Other checks:'
  echo '         verbose output, insulator, forces(mode 1)'
  set cplst = ($testdir/{ctrl.cr3si6})
  set lmfargs1 = "--pr51 -vnit=2 --no-iactiv --time=6"
# set fitbas  = ' Cr RSMH= 2.959 3.303 1.275 EH= -0.307 -0.200 -0.200'
# set fitbas2 = ' Si RSMH= 1.729 1.727 -1.000 EH= -0.609 -0.200 0.000'
else if ($ext == "te") then
  echo '         Case te: molecular statics in an open structure'
#    echo '         There is only one symmetry-allowed degree of freedom.'
  set cplst = ($testdir/{ctrl.te})
# set fitbas = ' X1 RSMH= 1.410 1.370 -1.000 EH= -1.028 -0.277 0.000'
  set lmfargs2 = '-vminx=t --rs=1,1 -vnk=3 -vnit=20 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1d-4 -verefc=0 -verefa=0'
  set dfmaxntol1 = 0.002
  set drmsqtol1  = 1e-6
else if ($ext == "srtio3") then
  echo '         Case srtio3: an oxide with local orbitals.'
  set cplst = ($testdir/{ctrl.srtio3})
  set lmfargs1 = " "
  set dfmax1tol1 = 1e-5
  set drmsqtol1 = 2e-6
else if ($ext == "tio2") then
  echo '         Case tio2: example of relaxation'
  set cplst = ($testdir/{ctrl,site}.tio2)
  set lmfargs1 = " "
  set dfmax1tol1 = 1e-5
  set drmsqtol1 = 2e-6
else if ($ext == "cu") then
  echo '         Case cu: illustration of high-lying local orbitals'
  echo '                  and bands of Cu up to ~50 eV.'
  set cplst = ($testdir/{ctrl.cu,syml.cu})
  set lmfargs1 = " "
  set drmsqtol1 = 1e-6
else if ($ext == "na") then
  echo '         Case na: illustration of low- and high-lying local orbitals'
  set cplst = ($testdir/{ctrl.na})
  set lmfargs1 = " "
  set drmsqtol1 = 1e-6
else if ($ext == "c") then
  echo '         Case C: test of homogeneous background'
  set cplst = ($testdir/{ctrl.c})
  set lmfargs1 = " "
  set drmsqtol3 = 5e-6
else if ($ext == "crn") then
  echo '         Case CrN: test of CLS with core hole'
  set cplst = ($testdir/{ctrl.crn})
  set lmfargs1 = " "
  set drmsqtol1 = 1e-6
else if ($ext == "gas") then
  echo '         Case GaAs: GaAs with local d orbital.  Overlapping,'
  echo '                    space-filling spheres to compare partial and total dos'
  set cplst = ($testdir/{ctrl.gas})
else
  echo test.fp: No test case for $ext
  exit -1
endif
endif

if ( $?joblist == 0 ) then
set joblist = ($argv)
if ( $#joblist == 0 ) set joblist = (1 2 3 4 5)
endif

echo $joblist | grep 1 >/dev/null
if ($status) goto chk1e
cat <<EOF

         --- Test 1.  Basic check of programs lmfa,lmf ---
         Checks that program lmfa produces a sensible atm file
         and that program lmf iterates to the proper energy.

EOF
endif
if ($?quiet) then
else if ($ext == "gdn") then
cat <<EOF
         The gdn test illustrates the LDA+U implementation, and an LDA calculation with:

         the majority 4f levels in the core
         the minority 4f levels in the core with no charge (tests partial core occupancy)

EOF
else if ($ext == "cdte") then
cat <<EOF
         The cdte test illustrates the LDA+U implementation in various modes:

         It begins with a self-consistent LDA+U calculation with U=4 eV on Cd d, FLL limit.
         The main effects are:
           * to push down the Cd d levels by 2 eV relative to the LDA
           * reduce total energy is about 20mRy less binding relative to LDA (-0.389970 Ry)
             (there is presumably a corresponding shift in the free atom energy, but
             no attempt was made to calculate it)
           * increase the bandgap by 0.1 eV relative to LDA (0.52 eV) to 0.62 eV

         Starting from this density and density-matrix, three one-shot calculations are performed.

         1. A potential shift -2 eV on Cd d (IDU=4) is used in place of U=4 eV.
            The test verifies that both the total energy and the bandgap are unchanged,
            and that the density is almost self-consistent.

         2. An additional potential shift +1 eV on Cd s (IDU=4) is included.
            It has the effect of increasing the gap to 0.99 eV, and reducing the energy by 70 mRy.
            Ueff is determined to be 0.337 Ry to generate this potential shift.
            (Note that this pass somewhat reduces the Cd s density-matrix)

         3. A normal LDA+U is included on the Cd s orbital, U=0.337 in the FLL.
            The test verifies that the gap, small change in output density, and EHK
            are essentially identical to the results of test 2.
         
EOF
else if ($ext == "eras") then
cat <<EOF
         The ErAs test also illustrates the LDA+U implementation for:

         1. a case when the LDA puts f orbitals at the Fermi energy.

         2. demonstration of U on both d and f orbitals

         3. convergence to a metastable solution with a reasonable spin moment
            but wrong orbital moment.  

            For a much more stable solution, continue with
               cp $testdir/occnum2.eras occnum.eras
               rm mixm.eras wkp.eras dmats.eras
               lmf eras -vnit=30
            This aligns the orbital moment antiparallel to the spin moment.

            For a still more stable solution, continue with
               cp $testdir/occnum3.eras occnum.eras
               rm mixm.eras wkp.eras dmats.eras
               lmf eras -vnit=30
            This aligns the orbital moment parallel to the spin moment.

EOF
else if ($ext == "er") then
cat <<EOF
         The Er test also illustrates the LDA+U implementation with:

           1. spin-orbit coupling

           2. the initial occupation number matrix given in spherical harmonics

           3. use of local orbitals in LDA+U

EOF
else if ($ext == "cu") then
cat <<EOF
         The cu test also illustrates the following:

         1.  high-lying local orbitals (Cu 5s,5p,4d are included as local orbitals)

         2.  METAL=3 for BZ integration

         3.  bands mode (see command-line argument --band in last lmf invocation)

EOF
else if ($ext == "na") then
cat <<EOF
         The na test also illustrates the following:

         1.  compare the total energy for conventional and extended Na 2p orbitals
             After the test finishes, compare the three energies with
               grep ^c out.lmf.na
EOF
else if ($ext == "srtio3") then
cat <<EOF
         The srtio3 test also illustrates the following:

         1.  how to freeze augmented wave functions for a particular
             species (cf token FRZWF=t for species Sr).

         2.  local orbitals (Sr 4p and Ti 4p are included as local orbitals)
             
         3.  extended local orbitals
             (last step recalculated with extended Sr 4p and Ti 4p local orbitals)
             
         4.  APW addition to basis

         5.  Renormalized free O atom for better starting density

         6.  use of restart file with shifted atom positions

EOF
else if ($ext == "copt") then
cat <<EOF
         The copt test also checks and illustrates the following:

         1.  a spin-polarized case

         2.  Sampling with Fermi function at fairly high temperature (5mRy)

         3.  forces with correction mode 12 (FORCES=12)

         4.  two-kappa basis, with second kappa consisting of Co p only

	 5.  inequivalent kmxa,lmxa,rmt

EOF
else if ($ext == "te") then
cat <<EOF
         The te test also checks and illustrates the following:

         1.  a simple molecular relaxation with one symmetry-allowed degree of freedom

         2.  an spd-sp basis augmented by floating orbitals, and also by q-dependent APWs

         3.  Comparison of forces and total energy with and without floating orbitals, 
             and with and without APWs.

         lmf will first relax the atoms with the basis including floating orbitals.
	 
         After relaxation, the a new calculation is performed that remove floating orbitals
         but adding plane waves to the basis, so the total energy and forces may be compared. 
         The basis size is variable, but averages ~80 orbitals, a little more than the floating
         orbitals case (~70 orbitals). About 3 mRy is gained relative to the floating orbitals case.

         Note that KMXA=5 is used with the PW calculation.  It isn't necessary in this case; still the
         user is cautioned to monitor this parameter when using energy cutoffs higher than 3 Ry or so.

         As a last step the calculation is repeated at the relaxed position with only atom-centered
         MTO's (neither floating orbitals nor plane waves).  Elimination of floating orbitals reduces 
         the basis from 75 to 39 orbitals, and reduces the LDA total energy by about 4 mRy.
	 
         The forces are not strongly affected by the local orbitals (or APWs), as can be seen by
         looking at the maximum force after the last step.

EOF

else if ($ext == "zrt") then
cat <<EOF
         The zrt test also checks and illustrates the following:

         1.  the use of restart files in both binary and ascii forms

         2.  two-kappa basis

EOF

else if ($ext == "co") then
set strn = "Install the fplot plotting package to have this script create"
if ($?have_pldos && $?have_fplot) set strn = "This script will prompt you to see whether it should create "

cat <<EOF
         The co test also checks and illustrates the following:

         1.  a spin-polarized case

         2.  mixed tetrahedron/sampling method (METAL=4)

         3.  Constrained mixing (first spin is kept frozen, then charge, then both are allowed to change)

         4.  Broyden mixing

         5.  bands mode (see command-line argument --band in last lmf invocation)
             SO coupling and color weights for projection onto Co d channels are included.
             Two separate weights are made (one each for d majority and minority bands.)
             $strn
             a figure from the bands file.

EOF

else if ($ext == "cr3si6") then
cat <<EOF
         The cr3si6 test also checks and illustrates the following:

         1.  insulator

         2.  forces with correction mode 1 (FORCES=1)

         3.  verbose output

         4.  A case with low kmxa (kmxa=2)
EOF

endif

set refout=$testdir/out.lmf.$ext.gz testout=out.lmf.$ext
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk1e
endif
set pass
# goto chk12
query chk11 chk1e 'run this test'
chk11:
# ... Look for executables
findcmd chk11a rdcmd "$path" "optional"
chk11a:
findcmd chk11b lmf "$path" "$topdir"
chk11b:
findcmd chk11c lmfa "$path" "optional"
chk11c:

# goto chk1ch

# ... Setup: remove existing files and copy new ones
echo "$space rm -f {atm,eula,moms,mixm,rst,rsta,save,log,hssn,wkp,cv,bsmv,bnds,dmats,dmats-save,sigm,qpp}.$ext $testout"
             rm -f {atm,eula,moms,mixm,rst,rsta,save,log,hssn,wkp,cv,bsmv,bnds,dmats,dmats-save,sigm,qpp}.$ext $testout
if (! $?clean) then
echo "$space cp $cplst ."
             cp $cplst .
endif
# ... Run lmf program
if (! $?clean) then
  runrdcmd chk12 %11f $testout "-cat:TESTLMF --noerr ctrl.$ext"
else
  if (-e ctrl.$ext) then
    runrdcmd chk1e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk1e
endif
chk12:

# ... Extract total energies, forces, magnetic moments 1st and last iter
extract_res_n chk12a efa erfa "etot=" 2 0 etot=
chk12a:
set ehf1  =  `cat $testout | grep ehf= | egrep -v '^   it' | head -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1  =  `cat $testout | grep ehk= | egrep -v '^   it' | head -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehf1r =  `zcat $refout | grep ehf= | egrep -v '^   it' | head -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1r =  `zcat $refout | grep ehk= | egrep -v '^   it' | head -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set dq1   =  `cat $testout | grep 'RMS DQ=' | head -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dq1r   = `zcat $refout | grep 'RMS DQ=' | head -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqn   =  `cat $testout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`

egrep ' pwmode=[^0]' $testout >/dev/null
if (! $status) then
  set epw  =  `cat $testout | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set epwr =  `zcat $refout | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
endif

grep 'Etot(LDA+U)' $testout >/dev/null
if (! $status) then
  set eldau  =  `cat $testout | grep 'Etot(LDA+U)' | tail -1 | awk '{print $NF}'`
  set eldaur =  `zcat $refout | grep 'Etot(LDA+U)' | tail -1 | awk '{print $NF}'`
endif

grep 'Maximum Harris force' $testout >/dev/null
if (! $status) then
  set fmax1  = `cat $testout | grep 'Maximum Harris force' | head -1 | awk '{print $5}'`
  set fmax1r = `zcat $refout | grep 'Maximum Harris force' | head -1 | awk '{print $5}'`
  set fmaxn  = `cat $testout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
  set fmaxnr = `zcat $refout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
endif

grep mmom= $testout >/dev/null
if (! $status) then
set mmom1  =  `cat $testout      | grep mmom= | head -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmom1r =  `zcat $refout | grep mmom= | head -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomn  =  `cat $testout      | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomnr =  `zcat $refout | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
endif

if ($ext == "te") then
set ebig  = `cat $testout | grep ^C | head -1 | awk '{sub(".*ehf=","");sub("ehk=.*",""); print $0}'`
set ebigr = `zcat $refout | grep ^C | head -1 | awk '{sub(".*ehf=","");sub("ehk=.*",""); print $0}'`
endif

set ediff = `echo $efa $erfa  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} printf "%10.2E", k}'`
if (! $?quiet) then
  echo " "
  echo "$space Total energy last free atom      = $efa"
  echo "$space Total energy of reference        = $erfa"
  echo "$space                    difference    =  $ediff"
  echo ' '

  echo "$space first iteration Harris energy    = $ehf1"
  echo "$space first iteration reference energy = $ehf1r"
  set ediff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo ' '

  echo "$space first iteration K-Sham energy    = $eks1"
  echo "$space first iteration reference energy = $eks1r"
  set ediff = `echo $eks1 $eks1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehf1 $eks1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

  echo ' '
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmax1) then
  echo "$space first iteration maximum force    = $fmax1"
  echo "$space first iteration reference force  = $fmax1r"
  endif
  if ($?mmom1) then
    echo "$space first iteration magnetic moment  = $mmom1"
    echo "$space first iteration reference moment = $mmom1r"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif
  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo " "

  if ($?ebig) then
  echo "$space K-Sham energy, big basis         = $ebig"
  echo "$space corresponding reference energy   = $ebigr"
  set ediff = `echo $ebig $ebigr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo " "
  endif

  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  echo " "
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?ebig) then
  set ediff = `echo $ebig $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference large->small basis    = $ediff"
  echo " "
  endif

  if ($?epw) then
    echo "$space last iteration E(MTO + APW)      = $epw"
    echo "$space last iteration ref E(MTO + APW)  = $epwr"
    set ediff = `echo $epw $epwr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $ediff"
    echo " "
  endif

  if ($?eldau) then
    echo "$space last iteration Etot(LDA+U)       = $eldau"
    echo "$space last iteration ref E(LDA+U)      = $eldaur"
    set ediff = `echo $eldau $eldaur  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $ediff"
  endif

  if ($?fmaxn) then
  echo "$space last iteration maximum force     = $fmaxn"
  echo "$space last iteration reference force   = $fmaxnr"
  endif
  if ($?mmom1) then
  echo "$space last iteration magnetic moment   = $mmomn"
  echo "$space last iteration reference moment  = $mmomnr"
  set mdiff = `echo $mmomn $mmomnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"
  endif
  echo "$space last iter RMS input-output drho  = $dqn"
  echo "$space last iter reference RMS drho     = $dqnr"
  echo " "

  zcat $refout | grep RELAX >/dev/null
  if ($status == 0) then
    call showout chk13 RELAX
chk13:
    echo ' '
  endif

  call zdiffiles chk14 "CPU -1 $testout $refout"
chk14:
endif

if ($?have_pldos && $?have_fplot && ! $?quiet) then

  egrep '^PLOTBND' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk14c
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk14a chk14c "generate a picture for energy bands"
  chk14a:
  runrdcmd chk14b %11f . "-cat:PLOTBND --noerr ctrl.$ext"
  chk14b:
   echo "$space $plot -disp -pr10 -f plot.plbnds"
                $plot -disp -pr10 -f plot.plbnds
  chk14c:
endif


# ... Check that FA fit basis set is within tol of reference
#  echo "$fitbas" > tmp.$ext.ref
#  if ($?fitbas2) then
#    echo "$fitbas2" >> tmp.$ext.ref
#  endif
#  grep RSMH atm.$ext  > tmp.$ext.lmfa
#  cmp tmp.$ext.lmfa tmp.$ext.ref >/dev/null
#  set retval = $status
#  echo -n "$space lmfa fit basis identical to reference ? ..."
#  if ($retval == 0) then
#   echo yes
#  else
#    echo no
#    unset pass
#  endif

if ($?defatol1 == 0) set defatol1 = 2e-6
if ($?dehf1tol1 == 0) set dehf1tol1 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol1 == 0) set dmom1tol1 = 1e-4
if ($?dmomntol1 == 0) set dmomntol1 = 1e-4
if ($?dfmax1tol1 == 0) set dfmax1tol1 = 0.1
if ($?dfmaxntol1 == 0) set dfmaxntol1 = 0.1
if ($?drmsqtol1 == 0) set drmsqtol1 = 1e-4

# pass checks
chk1c:

# ... Check that FA total energy is within tol of reference
compare_res chk1ca "FA etot (last species)" $efa $erfa $defatol1  pass
chk1ca:

compare_res chk1cb "1st  iter ehf" $ehf1 $ehf1r $dehf1tol1 pass
chk1cb:

if (! $?fmax1) goto chk1cc
compare_res chk1cc "1st  iter max force" $fmax1 $fmax1r $dfmax1tol1 pass
chk1cc:

if (! $?mmom1) goto chk1cd
compare_res chk1cd "1st  iter mmom" $mmom1 $mmom1r $dmom1tol1 pass
chk1cd:

compare_res chk1ce "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chk1ce:

if ($?eldau) then
compare_res chk1ci "last iter E(LDA+U)" $eldau $eldaur $dehf1toln pass
chk1ci:
endif

if ($?epw) then
compare_res chk1cj "last iter E(MTO+PW)" $epw $epwr $dehf1toln pass
chk1cj:
endif

if ($?fmaxn) then
compare_res chk1cf "last iter max force" $fmaxn $fmaxnr $dfmaxntol1 pass
chk1cf:
endif

if ($?mmomn) then
compare_res chk1cg "last iter mmom" $mmomn $mmomnr $dmomntol1 pass
chk1cg:
endif

compare_res chk1ch "last iter RMS dq" $dqn $dqnr $drmsqtol1 pass
chk1ch:

# compare bnds to reference
if (-e bnds.$ext) then
if ! ($?bndstol) set bndstol = 1e-4
zcmpmfiles_res_0 chk1ck "Max deviation in bnds.$ext from reference" $bndstol pass 4 bnds.$ext $testdir/bnds.$ext.gz
chk1ck:
#  chk25a:
#  echo -n "$space ... files bnds.$ext and $testdir/bnds.$ext.gz equivalent to $ndig digits? ... "
#  if ($retval == 0) then
#    echo  yes
#  else
#  #    set ndig = 4
#  #    call zcmpnfiles chk1ck "$ndig bnds.$ext $testdir/bnds.$ext.gz"
#  #    chk1ck:
#  #    echo -n "no ... to $ndig digits? ... "
#    if ($retval == 0) then
#      echo yes
#    else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
#      echo ok "($retval difference(s) of $ncharfile)"
#    else
#      echo no "($retval difference(s) remaining of $ncharfile)"
#      unset pass
#    endif
#  endif
endif

if ($?pass) then
    echo "$space test 1 PASSED ($ext)"
else
    echo "$space test 1 FAILED ($ext)"
    set failed = ($failed 1)
endif

chk1e:


echo $joblist | grep 2 >/dev/null
if ($status) goto chk2e
cat <<EOF

         --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---

EOF
if ($?quiet) then
else if ($ext == "crn") then
cat <<EOF
         The CrN test case generates core-level spectroscopy for the
         1s state in N.  The self-consistent calculation proceeds with
         an electron missing from the N 1s core, which corresponds to
         the 'sudden approximation' (system relaxes instantanously
         from electron exited out of hole).

EOF
else if ($ext == "co") then
cat <<EOF
         The Co test case illustrates partial dos resolved by both l and m.
         Note: because the routine generating partial DOS within MT spheres
         does not properly symmetrize it, symmetry operations must be
         suppressed when resolving DOS by m.

EOF
else if ($ext == "gdn") then
set strn = "Install the fplot plotting package to have this script create"
if ($?have_pldos && $?have_fplot) set strn = "After the calculation completes, this script will ask whether you want to make"
cat <<EOF
         The GdN test case illustrates use of Mulliken analysis in 
         the spin-coupled case.  Mulliken analysis is a useful way
         to resolve total DOS into spin components in the spin-coupled case.
         This test will generate files tdos.gdn and dos.gdn, which you can 
         use with a graphics packages to emphasize orbital or spin character of the DOS.

         $strn
         a figure of total DOS with majority DOS blue, minority DOS red.

EOF
else if ($ext == "fe") then
set strn = "Install the fplot plotting package to have this script create"
if ($?have_pldos && $?have_fplot) set strn = "After the calculation completes, this script will ask whether you want to make"
cat <<EOF
         The fe test case illustrates both core-level spectroscopy and
         Mulliken analysis resolved by l and m.
         Note: lmf does not properly symmetrize the output in either calculation,
         so symmetry operations must be suppressed.

         From the Mulliken DOS, you can resolve the total DOS into orbital contributions.
         $strn
         a picture of the DOS resolved into spin and three groups of orbitals:
         (s+p+f), (d states of t_2 symmetry), (d states of e_g symmetry).

         Mulliken analysis is also useful for the spin-coupled case, enabling
         the resolution of total DOS into spin components.  The following will generate 
         a picture of the total DOS, where the e2 part of the d channel is colored in red.

         lmf --rs=0 -vso=t --mull:mode=2 -vnk=6 -vnit=1 fe
         mv dos.fe tdos.fe
         lmdos --nosym -vso=t --mull:mode=2 --dos:npts=1001:window=-.7,.8 -vnk=6 fe
         mv dos.fe dos-mull.fe
         echo 40 7 -9 10 | pldos -ef=0 -escl=13.6 -fplot '-lst=13,17' -ref:fn=tdos.fe:chan=1:scale dos-mull.fe

EOF
endif
set refout=$testdir/out.lmf-dos.$ext.gz testout=out.lmf-dos.$ext
#set pass; echo 'TEST'; goto chk256e
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk2e
endif
if ($ext == "gdn" && ! $?clean) then
if (! -e dmats-save.gdn || ! -e rst-save.gdn) then
  echo "$space ... skipping test : missing one or both of files dmats-save.gdn rst-save.gdn "
  echo "$space     Run test 1 for GdN to generate them automatically."
  goto chk2e
endif
endif

set pass
query chk21 chk2e 'run this test'
chk21:
# ... Look for executables
findcmd chk21a rdcmd "$path" "optional"
chk21a:
findcmd chk21b lmf "$path" "$topdir"
chk21b:
findcmd chk21c lmfa "$path" "optional"
chk21c:
findcmd chk21d lmdos "$path" "optional"
chk21d:

echo "$space cp $cplst ."
             cp $cplst .
echo "$space rm -f {atm,mixm,rst,save,log,hssn,wkp,dos,tdos,pdos,dos-mull,qpp}.$ext $testout"
             rm -f {atm,mixm,rst,save,log,hssn,wkp,dos,tdos,pdos,dos-mull,qpp}.$ext $testout
if (! $?clean) then
  egrep '^TESTCLS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk25
  runrdcmd chk22 %11f $testout "-cat:TESTCLS --noerr ctrl.$ext"
else
  if (-e ctrl.$ext) then
  runrdcmd chk2e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk2e
endif

chk22:
echo ' '
call zdiffiles chk23 "CPU -1 $testout $refout"
chk23:

if ($?have_pldos && $?have_fplot && ! $?quiet) then
  egrep '^PLOTCLS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk24c
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk24a chk24c "generate a picture for CLS"
  chk24a:
  runrdcmd chk24b %11f . "-cat:PLOTCLS --noerr ctrl.$ext"
  chk24b:
   echo "$space $plot -disp -pr10 -f plot.dos"
                $plot -disp -pr10 -f plot.dos
  chk24c:

  egrep '^PLOTMUL' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk24f
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk24d chk24f "generate a picture for Mulliken"
  chk24d:
  runrdcmd chk24e %11f . "-cat:PLOTMUL --noerr ctrl.$ext"
  chk24e:
   echo "$space $plot -disp -pr10 -f plot.dos"
                $plot -disp -pr10 -f plot.dos
  chk24f:
endif

# compare dos-cls to reference
if (! -e "$testdir/dos-cls.$ext.gz") goto chk256e
if ! ($?dosclstol) set dosclstol = 1e-4
call zcmpmfiles chk256 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# call zcmpnfiles chk256 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
chk256:
compare_res_0 chk25a "Max deviation in dos-cls from reference" $retval $dosclstol pass
chk25a:
# echo -n "$space ... files dos-cls.$ext and $testdir/dos-cls.$ext.gz equivalent to 6 digits? ... "
# if ($retval == 0) then
#   echo  yes
# else
#   call zcmpnfiles chk253 "4 dos-cls.$ext $testdir/dos-cls.$ext.gz"
#   chk253:
#   echo -n "no ... to 4 digits? ... "
#   if ($retval == 0) then
#     echo yes
#   else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
#     echo ok "($retval difference(s) of $ncharfile)"
#   else
#     echo no "($retval difference(s) remaining of $ncharfile)"
#     unset pass
#   endif
# endif
chk256e:

# Compare dos-mull to reference
if (! -e "$testdir/dos-mull.$ext.gz") goto chk25
if ! ($?dosmulltol) set dosmulltol = 1e-4
# call zcmpmfiles chk266 "6 dos-mull.$ext $testdir/dos-mull.$ext.gz"
# chk266:
# compare_res_0 chk26a "Max deviation in dos-mull from reference" $retval $dosmulltol pass
# chk26a:
call zcmpnfiles chk266 "6 dos-mull.$ext $testdir/dos-mull.$ext.gz"
chk266:
echo -n "$space ... files dos-mull.$ext and $testdir/dos-mull.$ext.gz equivalent to 6 digits? ... "
if ($retval == 0) then
  echo  yes
else
  call zcmpnfiles chk263 "3 dos-mull.$ext $testdir/dos-mull.$ext.gz"
  chk263:
  echo -n "no ... to 3 digits? ... "
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif

# --- check partial dos ---
chk25:
egrep '^TSTPDOS' ctrl.$ext >/dev/null
set retval = $status
if ($retval != 0) goto chk29
if (! $?clean) then
  runrdcmd chk26 %11f $testout "-cat:TSTPDOS --noerr ctrl.$ext"
else
  if (-e ctrl.$ext) then
  runrdcmd chk2e %11f $testout "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk2e
endif
chk26:

echo ' '
call zdiffiles chk27 "CPU -1 $testout $refout"
chk27:
if ($?have_pldos && $?have_fplot && ! $?quiet) then
  egrep '^PLOTDOS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk27c
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk27a chk27c "generate a picture for DOS"
  chk27a:
  runrdcmd chk27b %11f . "-cat:PLOTDOS --noerr ctrl.$ext"
  chk27b:
   echo "$space $plot -disp -pr10 -f plot.dos"
                $plot -disp -pr10 -f plot.dos
  chk27c:
endif

# Compare pdos to reference
if ! ($?pdostol) set pdostol = 1e-4
call zcmpmfiles chk274 "4 pdos.$ext $testdir/pdos.$ext.gz"
call zcmpnfiles chk274 "4 pdos.$ext $testdir/pdos.$ext.gz"
chk274:
compare_res_0 chk274a "Max deviation in pdos from reference" $retval $pdostol pass
chk274a:
#  if ($?quiet) then
#  else if ($ext == "co") then
#  cat <<EOF
#           Note: generation of m-resolved dos is numerically a rather delicate
#           procedure.  The check to determine whether the partial DOS agrees
#           with the reference file only compares a few digits.
#  EOF
#  endif
#  if ($?ndosdig == 0) set ndosdig = 3
#  call zcmpnfiles chk274 "4 pdos.$ext $testdir/pdos.$ext.gz"
#  chk274:
#  echo -n "$space ... files pdos.$ext and $testdir/pdos.$ext.gz equivalent to 4 digits? ... "
#  if ($retval == 0) then
#    echo  yes
#  else
#    call zcmpnfiles chk272 "$ndosdig pdos.$ext $testdir/pdos.$ext.gz"
#    chk272:
#    echo -n "no ... to $ndosdig digits? ... "
#    if ($retval == 0) then
#      echo yes
#    else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
#      echo ok "($retval difference(s) of $ncharfile)"
#    else
#      echo no "($retval difference(s) remaining of $ncharfile)"
#      unset pass
#    endif
#  endif

chk29:
if ($?pass) then
    echo "$space test 2 PASSED ($ext)"
else
    echo "$space test 2 FAILED ($ext)"
    set failed = ($failed 2)
endif

chk2e:

echo $joblist | grep 3 >/dev/null
if ($status) goto chk3e
cat <<EOF

         --- Test 3.  Check of miscellaneous special features, programs lmfa,lmf ---

EOF
endif
if ($?quiet) then
else if ($ext == "c") then
cat <<EOF

         The C test tests the code's implementation of homogeneous background mode
         It checks that program lmf generates the correct total energy and 
         ionization potential for a single C atom.

         Note that:

	*The total energy of the neutral atom computed by lmf (-74.996 Ry) is very close
         to the free-atom energy computed by lmfa (-74.995 Ry), and that the free-atom 
         potential is approximately self-consistent; 

	*Also the SUM of total energy of the ionized system (-74.443) and the
         interaction of a point charge with a homogeneous background E^I,
         estimated from E^I=9/5R, with 4*pi*R^3/3=vol => R=6.20 and E^I = 0.290 Ry
         evaluates to -74.443+0.290 = -74.153 Ry, is close to the total energy of
         the positively charged free ion as computed by lmfa (-74.171 Ry).

EOF
else if ($ext == "felz") then
cat <<EOF

         The Fe test tests the code's implementation of fixed-spin moment
         with and without spin-orbit coupling

EOF
endif
set refout=$testdir/out.lmf.neutral.$ext.gz testout=out.lmf.$ext
if ($ext == "c") then
set testhom
endif
if ($ext == "felz") then
set testfsm
set refout=$testdir/out.lmf.fsmom.$ext.gz testout=out.lmf.$ext
endif
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk3e
endif
set pass
query chk31 chk3e 'run this test'
chk31:
# ... Look for executables
findcmd chk31a rdcmd "$path" "optional"
chk31a:
findcmd chk31b lmf "$path" "$topdir"
chk31b:
findcmd chk31c lmfa "$path" "optional"
chk31c:

#  goto chk32

# ... Setup: remove existing files and copy new ones
echo "$space rm -f {mixm,rst,save,log,hssn,wkp,bsmv,bnds}.$ext"
             rm -f {mixm,rst,save,log,hssn,wkp,bsmv,bnds}.$ext
echo "$space cp $cplst ."
             cp $cplst .

#  goto xxxx

# ... Run lmf program
if (! $?clean  && $?testfsm) then
  runrdcmd chk32 %11f $testout "-cat:TESTFSM --noerr ctrl.$ext"
else if (! $?clean) then
  runrdcmd chk32 %11f $testout "-cat:TESTLMF --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk3e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk3e
endif
chk32:

# ... Extract total energies, forces, magnetic moments 1st and last iter
extract_res_n chk32a efa erfa "etot=" 2 0 etot=
chk32a:
set ehf1  =  `cat $testout | grep ehf= | egrep -v '^   it' | head -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1  =  `cat $testout | grep ehk= | egrep -v '^   it' | head -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehf1r =  `zcat $refout | grep ehf= | egrep -v '^   it' | head -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1r =  `zcat $refout | grep ehk= | egrep -v '^   it' | head -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set dq1   =  `cat $testout | grep 'RMS DQ=' | head -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dq1r   = `zcat $refout | grep 'RMS DQ=' | head -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqn   =  `cat $testout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`

grep 'Maximum Harris force' $testout >/dev/null
if (! $status) then
  set fmax1  = `cat $testout | grep 'Maximum Harris force' | head -1 | awk '{print $5}'`
  set fmax1r = `zcat $refout | grep 'Maximum Harris force' | head -1 | awk '{print $5}'`
  set fmaxn  = `cat $testout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
  set fmaxnr = `zcat $refout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
endif

grep mmom= $testout >/dev/null
if (! $status) then
set mmom1  =  `cat $testout      | grep mmom= | head -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmom1r =  `zcat $refout | grep mmom= | head -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomn  =  `cat $testout      | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomnr =  `zcat $refout | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
if ($?testfsm) then
compare_resf chk32b mmom1u mmom1ur 'Mag. moment:' 3 1 zzz
chk32b:
compare_resf chk32c mmom1 mmom1r 'Mag. moment:' 3 2 zzz
chk32c:
endif
endif

set ediff = `echo $efa $erfa  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} printf "%10.2E", k}'`
if (! $?quiet) then
  echo " "
  echo "$space Total energy last free atom      = $efa"
  echo "$space Total energy of reference        = $erfa"
  echo "$space                    difference    =  $ediff"
  echo ' '

  echo "$space first iteration Harris energy    = $ehf1"
  echo "$space first iteration reference energy = $ehf1r"
  set ediff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space first iteration K-Sham energy    = $eks1"
  echo "$space first iteration reference energy = $eks1r"
  set ediff = `echo $eks1 $eks1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehf1 $eks1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmax1) then
  echo "$space first iteration maximum force    = $fmax1"
  echo "$space first iteration reference force  = $fmax1r"
  endif
  if ($?mmom1) then
    echo "$space first iteration magnetic moment  = $mmom1"
    echo "$space first iteration reference moment = $mmom1r"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif
  if ($?testfsm) then
    echo "$space first iter unconstr. moment      = $mmom1u"
    echo "$space first iter unconstr. ref moment  = $mmom1ur"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmaxn) then
  echo "$space last iteration maximum force     = $fmaxn"
  echo "$space last iteration reference force   = $fmaxnr"
  endif
  if ($?mmom1) then
  echo "$space last iteration magnetic moment   = $mmomn"
  echo "$space last iteration reference moment  = $mmomnr"
  set mdiff = `echo $mmomn $mmomnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"
  endif

  echo "$space last iter RMS input-output drho  = $dqn"
  echo "$space last iter reference RMS drho     = $dqnr"
  echo " "

  zcat $refout | grep RELAX >/dev/null
  if ($status == 0) then
    call showout chk33 RELAX
chk33:
    echo ' '
  endif

  call zdiffiles chk34 "CPU -1 $testout $refout"
chk34:
endif

# ... Check that FA fit basis set is within tol of reference
#  echo "$fitbas" > tmp.$ext.ref
#  if ($?fitbas2) then
#    echo "$fitbas2" >> tmp.$ext.ref
#  endif
#  grep RSMH atm.$ext  > tmp.$ext.lmfa
#  cmp tmp.$ext.lmfa tmp.$ext.ref >/dev/null
#  set retval = $status
#  echo -n "$space lmfa fit basis identical to reference ? ..."
#  if ($retval == 0) then
#   echo yes
#  else
#    echo no
#    unset pass
#  endif

if ($?defatol3 == 0) set defatol3 = 2e-6
if ($?dehf1tol3 == 0) set dehf1tol3 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol3 == 0) set dmom1tol3 = 1e-4
if ($?dmomntol3 == 0) set dmomntol3 = 1e-4
if ($?dfmax1tol3 == 0) set dfmax1tol3 = 0.1
if ($?dfmaxntol3 == 0) set dfmaxntol3 = 0.1
if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-4

# pass checks
chk3c:

# ... Check that FA total energy is within tol of reference
compare_res chk3ca "FA etot (last species)" $efa $erfa $defatol3  pass
chk3ca:

compare_res chk3cb "1st  iter ehf" $ehf1 $ehf1r $dehf1tol3 pass
chk3cb:

if (! $?fmax1) goto chk3cc
compare_res chk3cc "1st  iter max force" $fmax1 $fmax1r $dfmax1tol3 pass
chk3cc:

if ($?testfsm) then
compare_res chk3ccc "1st  iter unconst. mmom" $mmom1u $mmom1ur $dmom1tol3 pass
chk3ccc:
endif

if (! $?mmom1) goto chk3cd
compare_res chk3cd "1st  iter mmom" $mmom1 $mmom1r $dmom1tol3 pass
chk3cd:

compare_res chk3ce "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chk3ce:

if ($?fmaxn) then
compare_res chk3cf "last iter max force" $fmaxn $fmaxnr $dfmaxntol3 pass
chk3cf:
endif

if ($?mmomn) then
compare_res chk3cg "last iter mmom" $mmomn $mmomnr $dmomntol3 pass
chk3cg:
endif

compare_res chk3ch "last iter RMS dq" $dqn $dqnr $drmsqtol3 pass
chk3ch:

# compare bnds to reference
if (-e bnds.$ext) then
set ndig = 4
call zcmpnfiles chk3ci "$ndig bnds.$ext $testdir/bnds.$ext.gz"
chk3ci:
echo -n "$space ... files bnds.$ext and $testdir/bnds.$ext.gz equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
#    set ndig = 4
#    call zcmpnfiles chk3cj "$ndig bnds.$ext $testdir/bnds.$ext.gz"
#    chk3cj:
#    echo -n "no ... to $ndig digits? ... "
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif

set saveehfn = $ehfn

#  xxxx:
#  echo 'FIX'
#  
#  set saveehfn = -74.997344
#  set ehfn = -74.424166
#  set refout=$testdir/out.lmf.ionized.$ext.gz testout=out.lmf.$ext
#  zcat $testdir/out.lmf.ionized.$ext.gz >out.lmf.$ext
#  goto chkx32

if ($?testfsm) goto chk3e
echo " "
echo "$space ... repeat for ionized case"
echo " "
set refout=$testdir/out.lmf.ionized.$ext.gz testout=out.lmf.$ext

#goto chkx32

# ... Setup: remove existing files and copy new ones
echo "$space rm -f {mixm,rst,save,log,hssn,wkp,bsmv,bnds}.$ext"
             rm -f {mixm,rst,save,log,hssn,wkp,bsmv,bnds}.$ext
echo "$space cp $cplst ."
             cp $cplst .

# ... Run lmf program
if (! $?clean) then
  runrdcmd chkx32 %11f $testout "-cat:TESTION --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk3e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk3e
endif
chkx32:

# ... Extract total energies, forces, magnetic moments 1st and last iter
extract_res_n chkx32a efa erfa "etot=" 2 0 etot=
chkx32a:
set ehf1  =  `cat $testout | grep ehf= | egrep -v '^   it' | head -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1  =  `cat $testout | grep ehk= | egrep -v '^   it' | head -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehf1r =  `zcat $refout | grep ehf= | egrep -v '^   it' | head -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1r =  `zcat $refout | grep ehk= | egrep -v '^   it' | head -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set dq1   =  `cat $testout | grep 'RMS DQ=' | head -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dq1r   = `zcat $refout | grep 'RMS DQ=' | head -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqn   =  `cat $testout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`

grep 'Maximum Harris force' $testout >/dev/null
if (! $status) then
  set fmax1  = `cat $testout | grep 'Maximum Harris force' | head -1 | awk '{print $5}'`
  set fmax1r = `zcat $refout | grep 'Maximum Harris force' | head -1 | awk '{print $5}'`
  set fmaxn  = `cat $testout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
  set fmaxnr = `zcat $refout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
endif

grep mmom= $testout >/dev/null
if (! $status) then
set mmom1  =  `cat $testout      | grep mmom= | head -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmom1r =  `zcat $refout | grep mmom= | head -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomn  =  `cat $testout      | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomnr =  `zcat $refout | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
endif

set ediff = `echo $efa $erfa  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} printf "%10.2E", k}'`

set ewald = `grep 'Energy for background' out.lmf.$ext | tail -1 | awk '{print $NF}'`

if (! $?quiet) then
  echo " "
  echo "$space Total energy last free atom      = $efa"
  echo "$space Total energy of reference        = $erfa"
  echo "$space                    difference    =  $ediff"
  echo ' '

  echo "$space first iteration Harris energy    = $ehf1"
  echo "$space first iteration reference energy = $ehf1r"
  set ediff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space first iteration K-Sham energy    = $eks1"
  echo "$space first iteration reference energy = $eks1r"
  set ediff = `echo $eks1 $eks1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehf1 $eks1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmax1) then
  echo "$space first iteration maximum force    = $fmax1"
  echo "$space first iteration reference force  = $fmax1r"
  endif
  if ($?mmom1) then
    echo "$space first iteration magnetic moment  = $mmom1"
    echo "$space first iteration reference moment = $mmom1r"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"

  if ($?fmaxn) then
  echo "$space last iteration maximum force     = $fmaxn"
  echo "$space last iteration reference force   = $fmaxnr"
  endif
  if ($?mmom1) then
  echo "$space last iteration magnetic moment   = $mmomn"
  echo "$space last iteration reference moment  = $mmomnr"
  set mdiff = `echo $mmomn $mmomnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"
  endif
  echo "$space last iter RMS input-output drho  = $dqn"
  echo "$space last iter reference RMS drho     = $dqnr"
  echo " "

  echo "$space Energy of charged system         = $ehfn"
  echo "$space Estat energy q*q/9/5/<r>         = $ewald"
  set ediff = `echo $ehfn $ewald | awk '{k=$1+$2; print k}'`
  echo "$space Corrected charged system energy  = $ediff"
  echo "$space Energy of neutral system         = $saveehfn"
  set ediff = `echo $ediff $saveehfn | awk '{k=$1-$2; print k}'`
  echo "$space difference                       = $ediff"

  zcat $refout | grep RELAX >/dev/null
  if ($status == 0) then
    call showout chkx33 RELAX
chkx33:
    echo ' '
  endif

  call zdiffiles chkx34 "CPU -1 $testout $refout"
chkx34:
endif

# ... Check that FA fit basis set is within tol of reference
#  echo "$fitbas" > tmp.$ext.ref
#  if ($?fitbas2) then
#    echo "$fitbas2" >> tmp.$ext.ref
#  endif
#  grep RSMH atm.$ext  > tmp.$ext.lmfa
#  cmp tmp.$ext.lmfa tmp.$ext.ref >/dev/null
#  set retval = $status
#  echo -n "$space lmfa fit basis identical to reference ? ..."
#  if ($retval == 0) then
#   echo yes
#  else
#    echo no
#    unset pass
#  endif

if ($?defatol3 == 0) set defatol3 = 2e-6
if ($?dehf1tol3 == 0) set dehf1tol3 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol3 == 0) set dmom1tol3 = 1e-4
if ($?dmomntol3 == 0) set dmomntol3 = 1e-4
if ($?dfmax1tol3 == 0) set dfmax1tol3 = 0.1
if ($?dfmaxntol3 == 0) set dfmaxntol3 = 0.1
if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-4

# pass checks
chkx3c:

# ... Check that FA total energy is within tol of reference
compare_res chkx3ca "FA etot (last species)" $efa $erfa $defatol3  pass
chkx3ca:

compare_res chkx3cb "1st  iter ehf" $ehf1 $ehf1r $dehf1tol3 pass
chkx3cb:

if (! $?fmax1) goto chkx3cc
compare_res chkx3cc "1st  iter max force" $fmax1 $fmax1r $dfmax1tol3 pass
chkx3cc:

if (! $?mmom1) goto chkx3cd
compare_res chkx3cd "1st  iter mmom" $mmom1 $mmom1r $dmom1tol3 pass
chkx3cd:

compare_res chkx3ce "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chkx3ce:

if ($?fmaxn) then
compare_res chkx3cf "last iter max force" $fmaxn $fmaxnr $dfmaxntol3 pass
chkx3cf:
endif

if ($?mmomn) then
compare_res chkx3cg "last iter mmom" $mmomn $mmomnr $dmomntol3 pass
chkx3cg:
endif

compare_res chkx3ch "last iter RMS dq" $dqn $dqnr $drmsqtol3 pass
chkx3ch:

# compare bnds to reference
if (-e bnds.$ext) then
set ndig = 4
call zcmpnfiles chkx3ci "$ndig bnds.$ext $testdir/bnds.$ext.gz"
chkx3ci:
echo -n "$space ... files bnds.$ext and $testdir/bnds.$ext.gz equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
#    set ndig = 4
#    call zcmpnfiles chkx3cj "$ndig bnds.$ext $testdir/bnds.$ext.gz"
#    chkx3cj:
#    echo -n "no ... to $ndig digits? ... "
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif


if ($?pass) then
    echo "$space test 3 PASSED ($ext)"
else
    echo "$space test 3 FAILED ($ext)"
    set failed = ($failed 3)
endif

chk3e:

echo $joblist | grep 4 >/dev/null
if ($status) goto chk4e

cat <<EOF

         --- Test case 4:  Spin-orbit coupling ---

EOF

if ($?quiet) then
else if ($ext == "felz") then
cat <<EOF
         The felz test computes the orbital moment in Fe.
         lmf calculates the orbital magnetic moment.

           * In the first part of this test only LzSz is used.
        
           * The APW basis with LzSz is also checked.

	   * In the second part the FULL SPIN ORBIT is used.

           * Symmetry operations must be suppressed at present.

           * Only 4x4x4 k points are used in this test.

EOF
else if ($ext == "gasls") then
cat <<EOF
         The GaAs test computes the energy bands at (0,0,0) (Gamma point),
         (1/4,1/4,1/4) and (1/2,1/2,1/2) (L point).
	 The spin-orbit splitting of the valence states is tested.

         This test checks SO coupling in conjunction with conventional local orbitals.

EOF
else if ($ext == "gaslc") then
cat <<EOF
	 The GaAs test with local orbitals and GW self-energy computes 
	 the energy bands at (0,0,0) (Gamma point),(1/4,1/4,1/4) 
	 and (1/2,1/2,1/2) (L point).
	 The spin-orbit splitting of the valence states is tested.

         This test also checks SO coupling in conjunction with 
         extended local orbitals and floating orbitals.

EOF
#  goto chk4c
endif
set pass
set refout=$testdir/out.lmf.lzsz.$ext.gz testout=out.lmf.lzsz.$ext
if (! -e $refout) then
  echo "$space ... skipping orbital moments test : missing reference file $refout"
  goto chk4c
endif
echo ' '
query chk41 chk4e 'run this test'
chk41:
set pass
if ($a == "s") goto chk4e
# ... Look for executables
findcmd chk41a rdcmd "$path" "$maindir"
chk41a:
findcmd chk41b lmf "$path" "$topdir"
chk41b:
findcmd chk41c lmfa "$path" "$topdir"
chk41c:

#goto chk42

# ... remove related files
echo "$space rm -f {atm,ctrl,fs,moms,mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds}.$ext"
             rm -f {atm,ctrl,fs,moms,mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds}.$ext
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  set testout=out.lmf.ls.$ext
  echo "$space rm -f $testout"
               rm -f $testout
  set testout=out.lmf.$ext
  echo "$space rm -f $testout"
               rm -f $testout
endif
# ... copy required files
echo "$space cp $cplst ."
             cp $cplst .
# ... Run lmf program
if (! $?clean) then
  runrdcmd chk42 %11f $testout "-cat:TESTSZ --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk4e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk4e
endif
chk42:

  set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`

egrep ' pwmode=[^0]' $testout >/dev/null
if (! $status) then
  set epw  =  `cat $testout | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set epwr =  `zcat $refout | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
endif

  set orbm  = `cat $testout      |egrep "total moment" | tail -1 | awk '{print $6}'`
  set orbmr = `gunzip -c $refout |egrep "total moment" | tail -1 | awk '{print $6}'`

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  echo " "
  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

  if ($?epw) then
    echo " "
    echo "$space last iteration E(MTO + APW)      = $epw"
    echo "$space last iteration ref E(MTO + APW)  = $epwr"
    set ediff = `echo $epw $epwr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $ediff"
  endif

  echo " "
  echo "$space Orbital magnetic moment          = $orbm"
  echo "$space        reference moment          = $orbmr"
  set ediff = `echo $orbm $orbmr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

echo ' '
call zdiffiles chk43 "CPU -1 $testout $refout"
chk43:

# pass checks
chk4achk:
if ($?dehf1toln == 0) set dehf1toln = 2e-6

compare_res chk4achka "Orbital moment" $orbm $orbmr $dorbmtol  pass
chk4achka:
compare_res chk4achkb "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chk4achkb:
compare_res chk4achkc "last iter ehk" $eksn $eksnr $dehf1toln pass
chk4achkc:
if ($?epw) then
compare_res chk4achkd "last iter E(MTO+PW)" $epw $epwr $dehf1toln pass
chk4achkd:
endif

if ($?pass) then
    echo "$space test 4a PASSED"
else
    echo "$space test 4a FAILED"
    set failed = ($failed 4)
endif

echo ' '

echo "$space Calculate orbital moment with FULL SPIN ORBIT"

set refout=$testdir/out.lmf.ls.$ext.gz testout=out.lmf.ls.$ext

# ... Run lmf program for full L\dotS
 
  runrdcmd chk44 %11f $testout "-cat:TESTSO --noerr ctrl.$ext" 

chk44:
  set orbm1  = `cat $testout      |egrep "total moment" | awk '{print $6}'`
  set orbm1r = `gunzip -c $refout |egrep "total moment" | awk '{print $6}'`
  echo " "
  echo "$space Orbital magnetic moment WITH FULL SPIN ORBIT= $orbm1"
  echo "$space         reference moment                    = $orbm1r"

echo ' '
call zdiffiles chk45 "CPU -1 $testout $refout"
chk45:

# ... Check that orbital moment is within tol of reference
compare_res chk4bchk "Orbital moment" $orbm $orbmr $dorbmtol  pass
chk4bchk:

if ($?pass) then
    echo "$space test 4b PASSED"
else
    echo "$space test 4b FAILED"
    set failed = ($failed 4)
endif

#  third SO test
chk4c:

set refout=$testdir/out.lmf.ls-bands.$ext.gz testout=out.lmf.$ext
if (! -e $refout) then
  echo "$space ... skipping band splitting test :  missing reference file $refout"
  goto chk4e
endif
echo ' '
query chk4c1 chk4e 'run this test'
chk4c1:
set pass
if ($a == "s") goto chk4e
# ... Look for executables
findcmd chk4c1a rdcmd "$path" "$maindir"
chk4c1a:
findcmd chk4c1b lmf "$path" "$topdir"
chk4c1b:
findcmd chk4c1c lmfa "$path" "$topdir"
chk4c1c:

# ... remove related files
echo "$space rm -f {ctrl,rst,syml,wkp,bnds}.$ext"
             rm -f {ctrl,rst,syml,wkp,bnds}.$ext
# ... copy required files
echo "$space cp $cplst ."
             cp $cplst .

# ... Run lmf program
if (! $?clean) then
  runrdcmd chk4c2 %11f $testout "-cat:TESTSO --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk4e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  goto chk4e
endif

chk4c2:
#  echo 'uncomment this line'
set refout=$testdir/out.lmf.ls-bands.$ext.gz testout=out.lmf.$ext

#  set strn    = `./extract-lines --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $testout | tail -1`
#  set statesG = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
#  set delta   = `echo $strn | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
#  set strn    = `./extract-lines --gzip --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $refout | tail -1`
#  set deltar  = `echo $strn | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
#  set strn    = `./extract-lines --n=$lineeval  'k=  0.50000  0.50000  0.50000' 1 $testout | tail -1`
#  set statesL = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
#  set strn    = `./extract-lines --n=$lineeval  'k=  0.25000  0.25000  0.25000' 1 $testout | tail -1`
#  set statesq = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `

#    echo " "
#    echo "$space The top valence states at:"
#    echo "$space Gamma point        : $statesG Ry"
#    echo "$space (1/4,1/4,1/4) point: $statesq Ry"
#    echo "$space L point            : $statesL Ry" 
#    echo "$space Spin-Orbit splitting at the Gamma point = $delta eV"
#    echo "$space        reference splitting              = $deltar eV"

if (! $?quiet) then
echo " "
echo "$space The top valence states at the calculated k-points"
set i = 0
set n = `grep -c ' bndfp:  kpt' $testout`
while ($i < $n)
 @ i = $i + 1
set strn = `./extract-lines --n=$lineeval  'bndfp:  kpt' $i $testout | tail -1`
set states = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
./extract-lines --n=1 'bndfp:  kpt' $i $testout | awk '{printf "     k = %s %s %s", $7,$8,$9}'
  echo " : $states Ry"
end
endif

set strn    = `./extract-lines --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $testout | tail -1`
set delta   = `echo $strn | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
set strnr   = `./extract-lines --gzip --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $refout | tail -1`
set deltar  = `echo $strnr | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
  echo "$space Spin-Orbit splitting at the Gamma point = $delta eV ( states" `echo $strn | awk -vn=$evalso '{printf "%8.4f and %8.4f \n", $(n),$(n+1)}'` ")"
  echo "$space        reference splitting              = $deltar eV"

echo ' '
call zdiffiles chk4c3 "CPU -1 $testout $refout"
chk4c3:

echo -n "$space Levels at the Gamma point match reference? ... " 
if ("$strn" == "$strnr") then
  echo  yes
else
  echo no
  unset pass
endif
compare_res chk4chkb "Splitting at the Gamma point" $delta $deltar $gmtol  pass
chk4chkb:

if ($?pass) then
    echo "$space test 4c PASSED"
else
    echo "$space test 4c FAILED"
    set failed = ($failed 4c)
endif


chk4e:

# --- Summary ---
echo ' '
if ($#failed <= 1) then
    echo "$space $testfile : all tests PASSED ($ext)"
    echo " "
    exit 0
else
    shift failed
    echo "$space $testfile : These tests FAILED:" $failed
    echo " "
    exit -1
endif

# ---------------- runjob --------------
exit
runjob:
  set quitjob=$retcall
  if ($outfile == ".") then
    echo "$space $callarg"
    echo " "
    $callarg
    set retval = $status
    if ($retval != 0) goto cleanup
    goto $quitjob
  endif

  if (`echo $outfile | awk '{print substr($1,1,2)}'` == '>>') then
    set appfile = `echo $outfile | awk '{print substr($1,3)}'`
    echo "$space $callarg  >> $appfile"
    $callarg >> $appfile
    set retval = $status
  else
    echo "$space $callarg  > $outfile"
    $callarg > $outfile
    set retval = $status
  endif
  if ($retval != 0) goto cleanup
  goto $quitjob

# ---------------- compare_resf --------------
# Extracts one element of a line in files $testout and $refout containing a keyword.
# Variables testout and refout point to file names and must be set beforehand ($refout is gzipped file)
# usage: compare_resf retcall testvar refvar keyword arg_number occur_number sed_strn
#   Variables testout and refout referring to file names must be set
#   testvar      : put result from file $testout into this variable
#   refvar       : put result from file $refout (compressed) into this variable
#   keyword    	 : string line must contain
#   arg_number 	 : extracts $arg_number'th entry in line, as defined by awk
#   occur_number : argument from $occur_number'th line; if zero, use last line
#   sed_strn     : purge this string from result before assigning
exit
compare_resf:
  set quitjob=$retcall
# echo $retcall $testvar $refvar $keyword $arg_number $occur_number $sed_strn
  set $testvar = `grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  set $refvar = `zcat $refout | grep "$keyword" | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  goto $quitjob

# ---------------- extract_res_n --------------
# Extracts nth token in a line containing a keyword
# usage: extract_res_n retcall testvar refvar keyword arg_number occur_number sed_strn
#   Variables testout and refout referring to file names must be set ($refout is gzipped file)
#   keyword      : string line must contain
#   testvar      : put result from file $testout into this variable
#   refvar       : put result from file $refout (compressed) into this variable
#   arg_number   : extracts $arg_number'th entry in line, as defined by awk
#   occur_number : argument from $occur_number'th line; if zero, use last line
#   sed_strn     : delete this string with from result before assigning
exit
extract_res_n:
  set quitjob=$retcall
#  echo $retcall $testvar $refvar $keyword $arg_number $occur_number $sed_strn
  set $testvar = `grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  set $refvar = `gunzip -c $refout | grep "$keyword" | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  goto $quitjob

# ---------------- runrdcmd --------------
exit
runrdcmd:
  set quitjob=$retcall
  if ($outfile == ".") then
    $rdcmd -f:$rdcmdfmt $callarg
    set retval = $status
    echo ' '
    if ($retval == 0) then
      echo "$space Job(s) completed successfully"
      goto $quitjob
    endif
  else
    if (`echo $outfile | awk '{print substr($1,1,2)}'` == '>>') then
      set appfile = `echo $outfile | awk '{print substr($1,3)}'`
      echo "$space $callarg  >> $appfile"
      exit
#      $callarg >> $appfile
      set retval = $status
    else
      echo "$space ... the following job(s) will be executed by invoking "\""rdcmd $callarg"\"
      $rdcmd -f:$rdcmdfmt --n $callarg
      echo "$space ... starting invocation of rdcmd:"
      echo "$space $rdcmd '-f:#rdcmd:%2f' $callarg  >& $outfile"
      $rdcmd '-f:rdcmd:%2f' $callarg >& $outfile
      set retval = $status
    endif
  endif

  if ($retval == 0) then
    echo "$space Job(s) completed successfully; output in $outfile"
    if ($?ladd0) then
      echo -n "         ..." ; $testdir/add0 $testout
    endif
    goto $quitjob
  else
    echo "$space ...oops... the following command returned with nonzero exit status:"
    echo -n "$space   "
    grep rdcmd: $outfile | tail -1 | sed 's/rdcmd:  //'
    goto cleanup
  endif

# ---------------- cleanup --------------
exit
cleanup:
  if ($retval != 0) echo "$space job returned with error status $retval"
  if ($retval != 0) echo "$space ... $testfile aborting"
  exit $retval

# ---------------- diffiles --------------
# calling argument should consist of four strings:
# 1st string = string that terminates diff
# 2nd string = integer that counts how many times terminator should occur before terminating
# 3nd string = first file name
# 4th string = second file name
# example: call diffiles chk69 "CPU 3 $testout $refout"
exit
diffiles:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  set files = ($callarg)
  set endstr = $files[1]
  shift files
  set nend = $files[1]
  shift files
  if ($nend == "-1") then
    set nend = `grep "$endstr" $files[1] | wc | awk '{print $1}'`
  endif

#    echo difffiles : $quitjob $nend
#    grep $endstr $files[1]

  query diff11 $quitjob "compare $files"
diff11:
  diff $files | awk -v endstr=$endstr -v nend=$nend -v endl=0 -v endr=0 '{if ($1 == "<" && endl < nend) print ; if ($1 == ">" && endr < nend) print ; if ($1 == ">" || $1 == "<" || endl >= nend && endr >= nend) ; else {print} ; if ($1 == "<" && $2 == endstr) {endl+=1}; if ($1 == ">" && $2 == endstr) {endr+=1};}' | head -50
  goto $quitjob

# ---------------- zdiffiles --------------
# calling argument should consist of four strings:
# 1st string = string that terminates zdiff
# 2nd string = integer that counts how many times terminator should occur before terminating
#              -1 -> last occurence
# 3nd string = first file name
# 4th string = second file name
# example: call zdiffiles chk69 "CPU 3 $testout $refout"
exit
zdiffiles:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  set files = ($callarg)
  set endstr = $files[1]
  shift files
  set nend = $files[1]
  shift files
  if ($nend == "-1") then
    set nend = `grep "$endstr" $files[1] | wc | awk '{print $1}'`
  endif

#    echo zdiffiles : $quitjob $nend
#    grep $endstr $files[1]

  query zdiff11 $quitjob "compare $files"
zdiff11:
  $testdir/zdiff $files | awk -v endstr="$endstr" -v nend=$nend -v endl=0 -v endr=0 '{if ($1 == "<" && endl < nend) print ; if ($1 == ">" && endr < nend) print ; if ($1 == ">" || $1 == "<" || endl >= nend && endr >= nend) ; else {print} ; if ($1 == "<" && $2 == endstr) {endl+=1}; if ($1 == ">" && $2 == endstr) {endr+=1};}' | head -50
  echo " "
  goto $quitjob

# ---------------- compare_res --------------
# Compares two numbers $testvar-$refvar and unsets $passvar if |testvar-refvar|<tol
# usage: compares_res retcall keyword testvar refvar tol passvar
#   keyword      : label (for printout)
#   testvar      : first number
#   refvar       : second number
#   tol          : tolerance
#   passvar      : $passvar is unset if |testvar-refvar|<tol
exit
compare_res:
  set quitjob=$retcall
# echo $retcall $keyword $testvar $refvar $tol $passvar
  echo -n "$space $keyword ($testvar) within tol ($tol) of reference ($refvar)? ... "
  if (`echo $testvar $refvar | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjob

# ---------------- compare_res_0 --------------
# Compares a number $testvar and unsets $passvar if |testvar|<tol
# usage: compares_res_0 retcall keyword testvar tol passvar
# Example:
# compare_res_0 chk274a "Max deviation in pdos from reference" $retval $pdostol pass
#   keyword      : label (for printout)
#   testvar      : first number
#   tol          : tolerance
#   passvar      : $passvar is unset if |testvar|<tol
exit
compare_res_0:
  set quitjob=$retcall
#  echo $retcall $keyword $testvar $tol $passvar
  echo -n "$space $keyword ($testvar) within tol ($tol)? ... "
  if (`echo $testvar 0 | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjob

# ---------------- zcmpmfiles_res_0 --------------
# Compares two files, stripping all but numerical fields.
# Checks for max absolute difference and unsets $passvar if difference<$tol
# Files with .gz or .Z extensions are assumed to be gzipped.
# usage: zcmpmfiles_res_0 retcall keyword testvar tol passvar ndig srcfile reffile
#   retcall      : return to this point in script on exit
#   keyword      : label (for printout)
#   tol          : tolerance in maximum allowed deviation
#   passvar      : $passvar is unset if |testvar|<tol
#   ndig         : number of digits numbers in file are stripped to 
#   srcfile      : first file to compare
#   reffile      : second file to compare
# Example:
# zcmpmfiles_res_0 chk1ck "Max deviation in bnds.$ext from reference" $bndstol pass 4 bnds.$ext $testdir/bnds.$ext.gz
exit
zcmpmfiles_res_0:
  set quitjobl=$retcall
# echo $retcall $keyword $tol $?passvar $ndig $srcfile $reffile

  unset retval
  call zcmpmfiles zcmpmfilesx "$ndig $srcfile $reffile"
zcmpmfilesx:
  echo -n "$space $keyword ($retval) within tol ($tol)? ... "
  if (`echo $retval 0 | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjobl

# ---------------- zcmpnfiles --------------
# Compares two files, treating each field as a number.
# call arguments should contain 3 strings: no-digits test-file reference-file
# Files with .gz or .Z extensions are assumed to be gzipped.
# Returns with retval = number of differences in reduced files
# Example :  call zcmpnfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $testdir/tmp1 $testdir/tmp2
exit
zcmpnfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)
  set digits = $zcmpnargs[1]
# set a = ' { for (i = NF; i > 0; --i) printf " %.'$digits'f", $i; printf "\n" }'
  set a = ' { for (i = 1; i <= NF; i++) { k = sprintf("%.'$digits'f",$i); if (k+k == 0) k = 0 ; printf "%s ", k}; printf "\n" }'

  set fn1 = $testdir/tmp_compnfile_1
  set fn2 = $testdir/tmp_compnfile_2
  if ("$zcmpnargs[2]:e" == 'gz' || "$zcmpnargs[2]:e" == 'Z') then
    set cat1 = 'gunzip -c'
  else    
    set cat1 = cat
  endif
  if ("$zcmpnargs[3]:e" == 'gz' || "$zcmpnargs[3]:e" == 'Z') then
    set cat2 = 'gunzip -c'
  else    
    set cat2 = cat
  endif

  $cat1  $zcmpnargs[2] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
  $cat2  $zcmpnargs[3] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2
  set ncharfile = `wc $fn1 | awk '{print $3}'`
  cmp $fn1 $fn2 >/dev/null
  set retval = $status

  if ($retval == 0) rm -f $fn1 $fn2 
  if ($retval == 0) goto $quitjob

  set retval = `cmp -l $fn1 $fn2 |& grep -v EOF | wc | awk '{printf "%d", $1}'`
  if ($retval == 0) set retval = '-1'
  rm -f $fn1 $fn2 
  goto $quitjob

# ---------------- zcmpmfiles --------------
# Compares two files, treating each field as a number.
# Call arguments should contain 3 strings: no-digits test-file reference-file
# files with .gz or .Z extensions are assumed to be gzipped.
# Returns with retval = max numerical difference
# Example :  call zcmpmfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $testdir/tmp1 $testdir/tmp2
exit
zcmpmfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)
  set digits = $zcmpnargs[1]
# set a = ' { for (i = NF; i > 0; --i) printf " %.'$digits'f", $i; printf "\n" }'
  set a = ' { for (i = 1; i <= NF; i++) { k = sprintf("%.'$digits'f",$i); if (k+k == 0) k = 0 ; printf "%s ", k}; printf "\n" }'

  set fn1 = $testdir/tmp_compnfile_1
  set fn2 = $testdir/tmp_compnfile_2
  if ("$zcmpnargs[2]:e" == 'gz' || "$zcmpnargs[2]:e" == 'Z') then
    set cat1 = 'gunzip -c'
  else    
    set cat1 = cat
  endif
  if ("$zcmpnargs[3]:e" == 'gz' || "$zcmpnargs[3]:e" == 'Z') then
    set cat2 = 'gunzip -c'
  else    
    set cat2 = cat
  endif

  $cat1  $zcmpnargs[2] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
  $cat2  $zcmpnargs[3] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2

  set retval = `diff -y --width=300 $fn1 $fn2 | grep '|' | awk -v top=0 '{n=split($0,a,"|"); n1=split(a[1],b1); n2=split(a[2],b2); { j=0; while (j++ < n1) if (j <= n1 && j<=n2) {x = (b1[j]-b2[j])>0?(b1[j]-b2[j]):(b2[j]-b1[j]); top = (top-x)>0?top:x; }}} END {printf "%12.4e\n", top}'`
  rm -f $fn1 $fn2 
  goto $quitjob

# ---------------- qprint (print only quiet not set) --------------
exit
qprint:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  echo "$callarg"
  goto $quitjob

# ---------------- showout --------------
exit
showout:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  echo ' '
  echo "$space ... Compare $callarg to line(s) in file $refout":
  grep "$callarg" $testout
  if (`cat $testout | grep "$callarg" | wc | awk '{print $1}'` > 1) echo ' ---'
  zcat $refout | grep "$callarg"
  goto $quitjob

# ---------------- findcmd --------------
# Finds an executable program within the supplied path
# Usage: findcmd return_label executable_command path_name make_path
# If $executable_command is not found, findcmd does one of the following:
# If make_path = 'no' : returns silently.
# Otherwise findcmd aborts with a message, which assumes
# $make_path is the path where $executable_command is made.
exit
findcmd:
set found = 'no'
foreach ac_dir ($path_name)
 if (-x $ac_dir/$prog_cmd) then
   set $prog_cmd = $ac_dir/$prog_cmd
   set found = 'yes'
   break
 endif
end
if (! $?quiet) then
  if ($found == 'yes') echo "$space ... using executable $ac_dir/$prog_cmd"
  if ($found == 'no')  echo "$space ... no executable $prog_cmd found in path"
endif
if ($found == 'no' && $make_path != "no") then
  echo "  "
  echo "  Sorry, $testfile cannot find program '"$prog_cmd"' it needs to execute."
  echo "  '"$prog_cmd"' was not found in supplied path, or in the following:"
  echo "        $topdir $maindir"
# echo "  ... This script ($testfile) requires binary "'"rdcmd"'" to run."
  echo "  You must create or put '"$prog_cmd"' in your path before invoking this script."
  echo "  Normally '"$prog_cmd"' is created as part of the installation process."
  echo "  Invoking '"make $prog_cmd"' in $make_path should create it."
  echo "  $testfile aborting ..."
  exit -1
endif
goto $retcall

# ---------------- query --------------
exit
query:
  unset skip
  if ($?slow != 0) then
    echo "$space"'*'"hit <return> to $callarg, s <return> to skip it."
    set a = ($<)
    if ($a == "") goto $retcall
    switch ($a)
      case "quit":
      case "q":
      case "a":
        exit
      case "i":
        unset slow
        breaksw
      case "s":
        set skip
        breaksw
      case "t":
        time
        goto query
      default:
        echo 'q to quit; i unsets slow; s skips this job, t shows time'
        goto query
    endsw
  endif
  if ($?skip) goto $retcall2
  goto $retcall

# ---------------- usage: --------------
usage:
cat <<EOF
 usage: test.fp [switches] [file-extension] [testcase-list]
        e.g., "test.fp copt 1"
        If file-extension is missing, test.fp uses copt
        Switches:
        --all        run through all the test cases set up
        --clean      clean up files generated by this script
        --add0       add suppressed zeros to fortran output
        --no-iactive runs tests without prompting user
        --quiet runs tests with minimal output and without prompting user
#       --verbose    script prints out extra information
EOF
exit -1
