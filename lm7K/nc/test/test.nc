#!/bin/csh -f

# This file is a shell script testing operation of noncollinear ASA package

alias call 'set retcall = \!\!:2 ; set callarg = \!\!:3 ; goto \!\!:1'
alias runjob 'set retcall = \!\!:1; set outfile = \!\!:2 ; set callarg = \!\!:3 ; goto runjob'
alias runrdcmd 'set retcall = \!\!:1; set rdcmdfmt = \!\!:2 ; set outfile = \!\!:3 ; set callarg = \!\!:4 ; goto runrdcmd'
alias findcmd  'set retcall = \!\!:1 ; set prog_cmd = \!\!:2 ; set path_name = \!\!:3 ; set make_path = \!\!:4 ; goto findcmd'
alias query 'set retcall = \!\!:1 ; set retcall2 = \!\!:2 ; set callarg = \!\!:3 ; goto query'
alias compare_res 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set refvar = \!\!:4 ; set tol = \!\!:5 ; set passvar = \!\!:6 ; goto compare_res'
alias compare_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set tol = \!\!:4 ; set passvar = \!\!:5 ; goto compare_res_0'
alias compare_resf 'set retcall = \!\!:1; set testvar = \!\!:2 ; set refvar = \!\!:3 ; set keyword = \!\!:4  ; set arg_number = \!\!:5 ; set occur_number = \!\!:6 ; set sed_strn = \!\!:7 ; goto compare_resf'
alias cnvt_d_fmt  'set retcall = \!\!:1; set testvar = \!\!:2 ; set testval = \!\!:3 ; goto cnvt_d_fmt'

set allargs = ($argv)
set a
set slow
set testfile = $0
set testdir = $testfile:h
#set topdir  = `cd $testdir/../..; pwd`
set topdir  = `$testdir/../../startup/absolute-path $testdir/../..`
set maindir = $topdir/main
set space = '        '
set failed = 0
alias zcat 'gunzip -c'
set eojob

# Prepend current working-directory, top-level dir and maindir to path
set path = ($cwd $topdir $maindir $path)

# --- Pick off switches ---
while (`echo $1 | sed -e 's/\(.\).*/\1/' `  ==  "-")

  set arg1 = $1; shift
  if ($?verb) echo test.pgf: parsing switch $arg1
  switch ($arg1)
    case "--quiet":
      set quiet
      unset slow
      breaksw
    case "--clean":
      set clean
      breaksw
    case "--add0":
      set ladd0
      breaksw
    case "--no-iact*":
      unset slow
      breaksw
    case "--verb*":
      set verb = 1
      breaksw

    case "--all":
      echo "    ... invoke $testdir/test.so for spin-orbit coupling checks ..."
      $testdir/test.so `echo $allargs | sed s/--all//g`
      set retval = $status
      echo "    ... completed $testdir/test.so (spin-orbit coupling checks)"
      if ($retval != 0) then
        echo " $testfile : failed so tests ... aborting"
        exit -1
      endif

      set joblist
      while (`echo $1 | sed -e 's/\([0-9][0-9]*\)/-/'`  ==  "-")
        set joblist = ($joblist $1)
        shift
      end

      set pass
      set failed
      $testfile `echo $allargs | sed s/--all//g | sed -e 's/\([0-9][0-9]*\)//g'` $joblist
      set retval = $status
      if ($retval != 0) then
        unset pass
#       set failed = ($failed $i)
        echo " $testfile : failed nc tests ... aborting"
        exit -1
      endif
      echo "     ... $testfile : passed all checks (so, nc)"
      exit

    default:
      echo unrecognized switch $arg1
      goto usage
  endsw

end

set joblist = (`echo $argv | sed s/fe2//`)
if ($#joblist == 0 ) set joblist = (1 2 3 4 5 6 7)

echo $joblist | grep 1 >/dev/null
if ($status) goto chk1e

cat <<EOF

         --- Test case 1: 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests that equilibrium 3k structure is self-consistent.
         Test DESTROYS or overwrites *.fccfe
EOF
if ($?quiet) then
else
cat <<EOF
         Inspect input file for other magnetic structures to test.
EOF
endif
set pass
set refout=$testdir/out.fccfe.3k.equil.gz testout=out.fccfe
echo ' '
query chk11 chk1e 'run this test'
chk11:
# ... Look for executables
findcmd chk11a rdcmd "$path" "$maindir"
chk11a:
findcmd chk11b lm "$path" "$topdir"
chk11b:
findcmd chk11c lmstr "$path" "$maindir"
chk11c:

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f *.fccfe"
             rm -f *.fccfe
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk1e
endif
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .
runjob chk12 $testout "lmstr fccfe"
chk12:
runjob chk13 $testout "lm -vnit=0 -vsdyn=t --keepsignm fccfe --no-iactive"
chk13:
runjob chk14 $testout "lm -vnit=1 -vsdyn=t --keepsignm fccfe --no-iactive"
chk14:
echo "$space Program lm returned successfully."
if ($?ladd0) then
 echo -n "         ..." ; $testdir/add0 $testout
endif

set dq = `cat $testout | grep SV: | tail -1 | awk '{print $3}' | sed 's/D/e/'`
awk '{if ($1 == "MAGTRQ:") {getline;{getline;if ($3*$3 > 1e-12) {exit -1}};{getline;if ($3*$3 > 1e-12) {exit -1}};{getline;if ($3*$3 > 1e-12) {exit -1}};{getline;if ($3*$3 > 1e-12) {exit -1}}}}' $testout
set forcesmall = (! $status)

if ($?quiet) then
  else
  echo ' '
  echo "$space ...show magnetic forces"
  awk '{if ($1 == "MAGTRQ:") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' $testout

  echo ' '
  call zdiffiles chk15 "CPU 1 $testout $refout"
chk15:

call showout chk15a SV
chk15a:

endif

set refout=$testdir/out.fccfe.3k.ss.equil.gz testout=out.fccfe
echo ' '
echo "$space ... repeat, for 3k+SS"
echo "$space rm -f *.fccfe"
             rm -f *.fccfe
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .
runjob chk16a $testout "lmstr fccfe"
chk16a:
runjob chk16b $testout "lm -vnit=0 fccfe -vqss=.5 --keepsignm -vsdyn=f fccfe --no-iactive"
chk16b:
runjob chk16c $testout "lm -vnit=1 fccfe -vqss=.5 --keepsignm -vsdyn=f fccfe --no-iactive"
# runjob chk16c $testout "mpix -np=8 lm-MPIK -vmet=2 -vnit=1 fccfe -vqss=.5 --keepsignm -vsdyn=f fccfe --no-iactive"
chk16c:
echo "$space Program lm returned successfully."
if ($?ladd0) then
 echo -n "         ..." ; $testdir/add0 $testout
endif

set dqss = `cat $testout | grep SV: | tail -1 | awk '{print $3}' | sed 's/D/e/'`

if ($?quiet) then
else

  call showout chk18a SV
chk18a:

  echo ' '
  call zdiffiles chk18 "CPU 1 $testout $refout"
chk18:

endif

set refout=$testdir/out.fccfe.2ss.dnf.gz testout=out.fccfe
echo ' '
echo "$space ... Check ++-- spin structure using SS, downfolding Fe p orbitals"
echo "$space rm -f *.fccfe"
             rm -f *.fccfe
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .

runjob chk19a $testout "lmstr fccfe -cstrx=2ss"
chk19a:
runjob chk19b $testout "lm -vnit=0 fccfe -vqss=0.5 -cstrx=2ss --keepsignm -vsdyn=f --iactiv=no -vnk=32 -vidxp=2"
chk19b:
#  echo "$space cp a.fccfe a1.fccfe"
#               cp a.fccfe a1.fccfe
#  echo "$space cp a.fccfe a2.fccfe"
#               cp a.fccfe a2.fccfe
runjob chk19c $testout "lm -vnit=1 fccfe -vqss=0.5 -cstrx=2ss --keepsignm -vsdyn=f --iactiv=no -vnk=32 -vidxp=2"
#  asa-soft-link-str fccfe
# runjob chk19c $testout "mpix -np=8 lm-MPIK -vmet=2 -vnit=1 fccfe -vqss=0.5 -cstrx=2ss --keepsignm -vsdyn=f --iactiv=no -vnk=32 -vidxp=2"
chk19c:
echo "$space Program lm returned successfully."
if ($?ladd0) then
 echo -n "         ..." ; $testdir/add0 $testout
endif

set ehf2ss = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehfref = `zcat $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
#  set ehf2ss = `cat $testout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`
#  set ehfref = `zcat $refout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`

set ehfref48 = -0.8884846
set ediff = `echo $ehf2ss $ehfref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediff48 = `echo $ehf2ss $ehfref48  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`


if ($?quiet) then
else

  echo " "
  echo "$space calculated harris energy =          $ehf2ss"
  echo "$space compare to reference (32 divisions) $ehfref"
  echo "$space difference                           $ediff"
  echo "$space compare to reference (48 divisions) $ehfref48"
  echo "$space difference                           $ediff48"
  echo "$space NB: compare to total energy of ++-- spin configuration calculated from 4-atom structure (next test)"

  echo ' '
  call zdiffiles chk19d "CPU 1 $testout $refout"
chk19d:

endif

set refout=$testdir/out.fccfe.4ss.dnf.gz testout=out.fccfe
echo ' '
echo "$space ... Check ++-- spin structure using 4-atom cell, downfolding Fe p orbitals."
echo "$space     Starting moments are identical to SS. Confirm that the 4-atom energy exactly doubles that"
echo "$space     of the 2-atom SS structure, and that the two are self-consistent in the same potential."
echo "$space rm -f *.fccfe"
             rm -f *.fccfe
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .

runjob chk19e $testout "lmstr fccfe -cstrx=4ss"
chk19e:
runjob chk19f $testout "lm -vnit=0 fccfe -vqss=0.0 -cstrx=4ss --keepsignm -vsdyn=f --iactiv=no -vnk=32 -vidxp=2"
chk19f:
runjob chk19g $testout "lm -vnit=1 fccfe -vqss=0.0 -cstrx=4ss --keepsignm -vsdyn=f --iactiv=no -vnk=32 -vidxp=2"
chk19g:
echo "$space Program lm returned successfully."
if ($?ladd0) then
 echo -n "         ..." ; $testdir/add0 $testout
endif

set ehf4ss = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehfref = `zcat $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
#  set ehf4ss = `cat $testout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`
#  set ehfref = `zcat $refout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`

set ehfref48 = -1.7769874
set ediff4 = `echo $ehf4ss $ehfref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediff48 = `echo $ehf4ss $ehfref48  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`


if ($?quiet) then
else

  echo " "
  echo "$space calculated harris energy =          $ehf4ss"
  echo "$space compare to reference (32 divisions) $ehfref"
  echo "$space difference                           $ediff4"
  echo "$space compare to reference (48 divisions) $ehfref48"
  echo "$space difference                           $ediff48"
  echo "$space NB: last energy matches ++-- spin configuration calculated from 4ss structure"

  echo ' '
  call zdiffiles chk19h "CPU 1 $testout $refout"
chk19h:

endif

# ... automatic pass checks
chk1p:
echo ' '
call qprint chk1pa "$space ... automatic pass checks :"
chk1pa:

echo -n "$space rms dq (=$dq) < 1e-6 in 3k structure ? ... "
if (`echo $dq | awk '{print ($1 < 1e-6)}'`) then
  echo  yes
else
  echo no
  unset pass
endif

echo -n "$space magnetic forces < 1e-6 in 3k structure ? ... "
if ($forcesmall) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space rms dq (=$dqss) < 1e-6 in 3k+SS structure ? ... "
if (`echo $dqss | awk '{print ($1 < 1e-6)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space difference in Harris energy for 2SS (=$ediff) < 2e-6 ? ... "
if (`echo $ediff | awk '{print ($1 < 2e-6)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
if ($?pass) then
    echo "$space test 1  PASSED"
else
    echo "$space test 1 FAILED"
    set failed = ($failed 1)
endif

chk1e:

echo $joblist | grep 2 >/dev/null
if ($status) goto chk2e
cat <<EOF

         --- Test case 2: SD in 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests spin statics using sdmod=1 for 4 atoms/cell Fe.
         Test DESTROYS or overwrites *.fccfe.
         NB: the HF and HK energy functionals don't exactly agree because
         the atom potentials are artifically averaged (GRP2=1).

EOF
set eojob=chk2p  sdmod=1  eulafile=eula.ran
goto sdjob
chk2p:
if ($?clean) then
  goto chk2e
endif
echo ' '
if ($?quiet) then
else
echo "$space ... automatic pass checks :"
endif

if ($?drmsqtol2 == 0) set drmsqtol2 = 1e-5
compare_res chk2cb "last iter RMS dq" $dqn $dqnr $drmsqtol2 pass
chk2cb:

if ($?drmsetol2 == 0) set detol2 = 2e-6
compare_res chk2cc "last iter ehf" $etest $eref $detol2 pass
chk2cc:

echo ' '
if ($?pass) then
    echo "$space test 2 PASSED"
else
    echo "$space test 2 FAILED"
    set failed = ($failed 2)
endif

chk2e:
if ($eojob == "chk3p") goto chk3e
if ($eojob == "chk4p") goto chk4e

echo $joblist | grep 3 >/dev/null
if ($status) goto chk3e
cat <<EOF

         --- Test case 3: SD in 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests spin statics using sdmod=11 for 4 atoms/cell Fe.
         Test DESTROYS or overwrites *.fccfe
         NB: the HF and HK energy functionals don't exactly agree because
         the atom potentials are artifically averaged (GRP2=1).

EOF
set eojob=chk3p sdmod=11  eulafile=eula.ran
goto sdjob
chk3p:
if ($?clean) then
  goto chk3e
endif
echo ' '
if ($?quiet) then
else
echo "$space ... automatic pass checks :"
endif

if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-5
compare_res chk3cb "last iter RMS dq" $dqn $dqnr $drmsqtol3 pass
chk3cb:

if ($?drmsetol3 == 0) set detol3 = 1e-5
compare_res chk3cc "last iter ehf" $etest $eref $detol3 pass
chk3cc:


echo ' '
if ($?pass) then
    echo "$space test 3 PASSED"
else
    echo "$space test 3 FAILED"
    set failed = ($failed 3)
endif

chk3e:

echo $joblist | grep 4 >/dev/null
if ($status) goto chk4e
cat <<EOF

         --- Test case 4: SD in 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests spin statics using sdmod=1011 for 4 atoms/cell Fe.
         Test also uses l-dependent Euler angles.
         Test DESTROYS or overwrites *.fccfe
EOF
set eojob=chk4p sdmod=1011 eulafile=eula-l.ran
goto sdjob
chk4p:
if ($?clean) then
  goto chk4e
endif
echo ' '
if ($?quiet) then
else
echo "$space ... automatic pass checks :"
endif

if ($?drmsmtol4 == 0) set drmsmtol4 = 1e-5
compare_res chk4cb "last iter RMS change in euler angles" $dmn $dmnr $drmsmtol4 pass
chk4cb:

if ($?drmsetol4 == 0) set detol4 = 1e-5
compare_res chk4cc "last iter ehf" $etest $eref $detol4 pass
chk4cc:

echo ' '
if ($?pass) then
    echo "$space test 4 PASSED"
else
    echo "$space test 4 FAILED"
    set failed = ($failed 4)
endif

chk4e:

echo $joblist | grep 5 >/dev/null
if ($status) goto chk5e
cat <<EOF

         --- Test case 5: Applied magnetic field ($testdir/ctrl.fe)  ---
         Tests the application of an external magnetic field 4 atoms/cell Fe
         Test DESTROYS or overwrites *.fe
         See $testdir/Notes-testing-bfield for further checks.
EOF
set pass
set refout=$testdir/out.fe.matchcoll.gz testout=out.fe
echo ' '
query chk51 chk5e 'run this test'
chk51:
# ... Look for executables
findcmd chk51a rdcmd "$path" "$maindir"
chk51a:
findcmd chk51b lm "$path" "$topdir"
chk51b:
findcmd chk51c lmstr "$path" "$maindir"
chk51c:

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f {ctrl,moms,wkp,site2,log,fe,fe2,eula,bfield,save,sdot,str,sv}.fe"
             rm -f {ctrl,moms,wkp,site2,log,fe,fe2,eula,bfield,save,sdot,str,sv}.fe
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk5e
endif
echo "$space cp $testdir/{ctrl,site2}.fe ."
             cp $testdir/{ctrl,site2}.fe .
runjob chk52a $testout "lmstr -vnit=0 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1"
chk52a:
runjob chk52b $testout "lm -vnit=0 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1"
chk52b:

echo ' '
echo "$space ... Verify that collinear and noncollinear branches produce same result:"
cat >eula.fe <<in
.1 .2 .3
.1 .2 .3
in
echo "$space using Euler angles file:"
cat eula.fe
runjob chk53a "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=0 -vbf=0 --quit=band"
chk53a:
runjob chk53b ">>$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=0 --quit=band"
chk53b:

set sumevc =      `cat $testout | grep 'sumev=' | head -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set sumev0 =      `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set sumevr = `gunzip -c $refout | grep 'sumev=' | head -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf0   =      `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk0   =      `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehfr   = `gunzip -c $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `gunzip -c $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif
set mz0 = `$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $5}'`

echo "$space    collinear band structure energy = $sumevc"
echo "$space noncollinear band structure energy = $sumev0"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumevc $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf0"
echo "$space          reference   Harris energy = $ehfr"
set ehfdiff = `echo $ehf0 $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk0"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk0 $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf0 $ehk0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"

echo ' '
call qprint chk54a "$space ... automatic pass checks :"
chk54a:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev0 $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf0 $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

call zdiffiles chk55a "CPU 1 $testout $refout"
chk55a:


# ... Stoner susceptibility 
set refout=$testdir/out.fe.stonerI.gz testout=out.fe
echo ' '
echo "$space ... Check longitudinal (Stoner susceptibility)"
cat >eula.fe <<in
0 0 0
0 0 0
in
echo "$space using Euler angles file:"
cat eula.fe
cat >bfield.fe <<in
0 0 .001
0 0 0
in
echo "$space using B-field file:"
cat bfield.fe

gunzip -c $refout >$testout
set mzref  = `$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $5}'`
set sumevr = `cat $testout | grep 'sumev=' | head -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehfr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

runjob chk56a "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
chk56a:

set mz     = `$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $5}'`
set sumev  = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif

echo "$space noncollinear band structure energy = $sumev"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumev $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf"
echo "$space         reference    Harris energy = $ehfr"
set ehfdiff = `echo $ehf $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf $ehk | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"
echo "$space output moment on atom 1            = $mz"
set mzdiff = `echo $mz $mz0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space induced moment = M(B)-M(B=0)       = $mzdiff"
echo "$space reference output moment            = $mzref"
set mzdiff = `echo $mz $mzref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $mzdiff"

echo ' '
call qprint chk57a "$space ... automatic pass checks :"
chk57a:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
call zdiffiles chk58a "CPU 1 $testout $refout"
chk58a:

# ... transverse susceptibility
set refout=$testdir/out.fe.bfieldx.gz testout=out.fe
echo ' '
echo "$space ... Check transverse susceptibility"
cat >eula.fe <<in
0 0 0
0 0 0
in
echo "$space using Euler angles file:"
cat eula.fe
cat >bfield.fe <<in
.001 0 0
0 0 0
in
echo "$space using B-field file:"
cat bfield.fe

gunzip -c $refout >$testout
set mzref  = `$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $3}'`
set sumevr = `cat $testout | grep 'sumev=' | head -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehfr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

runjob chk59a "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
chk59a:

set mz     = `$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $3}'`
set sumev  = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
$testdir/../../testing/extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif

echo "$space noncollinear band structure energy = $sumev"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumev $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf"
echo "$space         reference    Harris energy = $ehfr"
set ehfdiff = `echo $ehf $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf $ehk | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"
echo "$space output moment on atom 1            = $mz"
echo "$space output (induced) moment on atom 1  = $mz"
echo "$space reference output moment            = $mzref"
set mzdiff = `echo $mz $mzref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $mzdiff"

echo ' '
call qprint chk5aa "$space ... automatic pass checks :"
chk5aa:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
call zdiffiles chk5ba "CPU 1 $testout $refout"
chk5ba:


echo ' '
if ($?pass) then
    echo "$space test 5 PASSED"
else
    echo "$space test 5 FAILED"
    set failed = ($failed 5)
endif

chk5e:

echo $joblist | grep 6 >/dev/null
if ($status) goto chk6e

set ext = fe2

cat <<EOF

         --- Test case 6: Exchange interactions  ---
         Makes a small rotation of a spin in the presence of an external
         constraining field to zero out the d component of the spin density
         matrix.  Calculation compared to rotation in the absence of
         constraining field.  NB: a careful calculation would need to
         pay attention to k-point convergence.

EOF

if ($ext == "fe2") then
cat <<EOF
         Test case Fe2 consists of 2 Fe atoms (simple cubic unit cell). 
         One atom is rotated by -theta/2 about z; the other is rotated by theta/2.
         Constraining field is taken to be along x axis.

EOF
endif

set pass

set refout=$testdir/out.$ext.coll.gz testout=out.$ext
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk6e
endif
echo ' '
query chk61 chk6e 'run this test'
chk61:
# ... Look for executables
findcmd chk61a rdcmd "$path" "$maindir"
chk61a:
findcmd chk61b lm "$path" "$topdir"
chk61b:
findcmd chk61c lmstr "$path" "$maindir"
chk61c:

findcmd chk61d rdfile "$path" "no"
chk61d:
if ("$found" == "no") then
  echo "$space ... no rdfile in path ... skipping test"
  goto chk6e
endif

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f {a,a2,eula0,b0,ctrl,site2,mixm,moms,log,wkp,bfield,eula,save,sdot,str,sv}.$ext"
             rm -f {a,a2,eula0,b0,ctrl,site2,mixm,moms,log,wkp,bfield,eula,save,sdot,str,sv}.$ext
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk6e
endif
echo "$space cp $testdir/{eula0,b0,ctrl,site2}.$ext ."
             cp $testdir/{eula0,b0,ctrl,site2}.$ext .

echo "$space ... generate total energy for collinear case"
runrdcmd chk62a %11f $testout "-cat:JOBCOLL --noerr ctrl.$ext"
chk62a:
echo -n "$space shortening output ... "
echo "$topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
      $topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext
call zdiffiles chk62b "CPU 2 $testout $refout"
chk62b:

set ehf0ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk0ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehf0    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk0    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ediffhk = `echo $ehk0 $ehk0ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediffhf = `echo $ehf0 $ehf0ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`


if (! $?quiet) then
echo ' '
echo "$space collinear ehk = $ehk0"
echo "$space reference     = $ehk0ref"
echo "$space difference    = $ediffhk"

echo "$space collinear ehf = $ehf0"
echo "$space     reference = $ehf0ref"
echo "$space difference    = $ediffhf"
echo ' '
endif

compare_res chk62c "ehk" $ehk0 $ehk0ref 1e-6 pass
chk62c:
compare_res chk62d "ehf" $ehf0 $ehf0ref 1e-6 pass
chk62d:

echo " "
echo "$space ... rotation in the absence of constraining field"
set refout=$testdir/out.$ext.nc1.gz
echo "$space rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log}.fe2"
             rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log}.fe2
echo "$space cp $testdir/{eula0,b0,ctrl,site2}.fe2 ."
             cp $testdir/{eula0,b0,ctrl,site2}.fe2 .
runrdcmd chk63a %11f $testout "-cat:JOBNC1 --noerr ctrl.$ext"
chk63a:
echo -n "$space shortening output ... "
echo "$topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
      $topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext

call zdiffiles chk63b "CPU 1 $testout $refout"
chk63b:

set ehf1ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk1ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehf1    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk1    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ediffhk = `echo $ehk1 $ehk1ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediffhf = `echo $ehf1 $ehf1ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

set mx1in    = `$topdir/startup/catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
set mx1dout  = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $3}'`
set mx1doutl = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $6}'`
set mx1out   = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
zcat $refout >$testout{}~
set rmx1in    = `$topdir/startup/catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
set rmx1dout  = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $3}'`
set rmx1doutl = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $6}'`
set rmx1out   = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
rm $testout{}~

if (! $?quiet) then

echo "$space input  Mx     = $mx1in"
echo "$space output Mx     = $mx1out"
echo "$space output Mx(d)  = $mx1dout"
echo "$space in loc. coord = $mx1doutl"
echo "$space reference     = $rmx1doutl"

echo ' '
echo "$space ehk           = $ehk1"
echo "$space reference     = $ehk1ref"
echo "$space difference    = $ediffhk"
echo "$space ehk-ehk(coll) = "`echo $ehk1 $ehk0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

echo ' '
echo "$space ehf           = $ehf1"
echo "$space reference     = $ehf1ref"
echo "$space difference    = $ediffhf"
echo "$space ehf-ehf(coll) = "`echo $ehf1 $ehf0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo ' '
endif

compare_res chk63c "ehk" $ehk1 $ehk1ref 1e-6 pass
chk63c:
compare_res chk63d "ehf" $ehf1 $ehf1ref 1e-6 pass
chk63d:
compare_res chk63e "Mx" $mx1out $rmx1out 1e-5 pass
chk63e:
compare_res chk63f "d part of Mx, local coordinates" $mx1doutl $rmx1doutl 1e-5 pass
chk63f:

echo " "
echo "$space ... rotation in the presence of constraining field"
set refout=$testdir/out.$ext.nc3.gz
echo "$space rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log}.fe2"
             rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log}.fe2
echo "$space cp $testdir/{eula0,b0,ctrl,site2}.fe2 ."
             cp $testdir/{eula0,b0,ctrl,site2}.fe2 .
runrdcmd chk64a %11f $testout "-cat:JOBNC3 --noerr ctrl.$ext"
chk64a:
echo -n "$space shortening output ... "
echo "$topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
      $topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext

call zdiffiles chk64b "CPU 1 $testout $refout"
chk64b:

set ehf3ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk3ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehf3    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk3    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ediffhk = `echo $ehk3 $ehk3ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediffhf = `echo $ehf3 $ehf3ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

set mx3in    = `$topdir/startup/catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
set mx3dout  = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $3}'`
set mx3doutl = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $6}'`
set mx3out   = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
zcat $refout >$testout{}~
set rmx3in    = `$topdir/startup/catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
set rmx3dout  = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $3}'`
set rmx3doutl = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $6}'`
set rmx3out   = `$topdir/startup/catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
rm $testout{}~

if (! $?quiet) then

echo "$space input  Mx     = $mx3in"
echo "$space output Mx     = $mx3out"
echo "$space output Mx(d)  = $mx3dout"
echo "$space in loc. coord = $mx3doutl"
echo "$space reference     = $rmx3doutl"

echo ' '
echo "$space ehk           = $ehk3"
echo "$space reference     = $ehk3ref"
echo "$space difference    = $ediffhk"
echo "$space ehk-ehk(coll) = "`echo $ehk3 $ehk0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space ehk-ehk(noB)  = "`echo $ehk3 $ehk1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

echo ' '
echo "$space ehf           = $ehf3"
echo "$space reference     = $ehf3ref"
echo "$space difference    = $ediffhf"
echo "$space ehf-ehf(coll) = "`echo $ehf3 $ehf0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space ehf-ehf(noB)  = "`echo $ehf3 $ehf1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo ' '
endif

compare_res chk64c "ehk" $ehk3 $ehk3ref 1e-6 pass
chk64c:
compare_res chk64d "ehf" $ehf3 $ehf3ref 1e-6 pass
chk64d:
compare_res chk64e "Mx" $mx3out $rmx3out 1e-5 pass
chk64e:
compare_res chk64f "d part of Mx, local coordinates" $mx3doutl $rmx3doutl 1e-5 pass
chk64f:

if ($?pass) then
    echo "$space test 6 PASSED ($ext)"
else
    echo "$space test 6 FAILED ($ext)"
    set failed = ($failed 6)
endif
chk6e:



echo $joblist | grep 7 >/dev/null
if ($status) goto chk7e

set ext = er

cat <<EOF

         --- Test case 7: noncollinear LDA+U hamiltonian ---

EOF

if ($ext == "er") then
cat <<EOF
         Test Er compares a collinear, antiferromagnetic ASA LDA+U
         calculation in hcp Er (one spin up, the other down) to a
         two kinds of equivalent noncollinear calculations:

         1.  Moments are equal and opposite in sign.
             Spins are aligned parallel, but not along z.

         2.  Moments are equal and identical in sign.
             The spinor of the second atom is rotated 180 degrees about y
             to recover the AFM result.
EOF
endif

set pass

set testout=out.lm.{$ext}.coll refout=$testdir/out.lm.$ext.coll.gz 
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk7e
endif
echo ' '
query chk71 chk7e 'run this test'
chk71:
# ... Look for executables
findcmd chk71a rdcmd "$path" "$maindir"
chk71a:
findcmd chk71b lm "$path" "$topdir"
chk71b:
findcmd chk71c lmstr "$path" "$maindir"
chk71c:

#  findcmd chk71d rdfile "$path" "no"
#  chk71d:
#  if ("$found" == "no") then
#    echo "$space ... no rdfile in path ... skipping test"
#    goto chk7e
#  endif

echo "$space ... set up ASA strux and starting potential"
touch ctrl.$ext
echo "$space rm -f *.$ext"
             rm -f *.$ext
if ($?clean) then
  echo "$space rm -f out.lm.{$ext}.coll out.lm.{$ext}.nc1 out.lm.{$ext}.nc2 specialspec1 atparms"
               rm -f out.lm.{$ext}.coll out.lm.{$ext}.nc1 out.lm.{$ext}.nc2 specialspec1 atparms
  goto chk7e
endif
echo "$space cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 $testdir/atparms ."
             cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 $testdir/atparms .
echo "$space cp $testdir/dmats.$ext.afm dmats.$ext"
             cp $testdir/dmats.$ext.afm dmats.$ext
echo "$space cp $testdir/rsta.$ext.afm rsta.$ext"
             cp $testdir/rsta.$ext.afm rsta.$ext


echo " "
echo "$space ... AFM calculation, collinear case"
runrdcmd chk72a %11f $testout "-cat:JOBCOLL --noerr ctrl.$ext"
chk72a:
#  echo -n "$space shortening output ... "
#  echo "$topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
#        $topdir/startup/catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext
#  call zdiffiles chk72b "CPU 2 $testout $refout"
#  chk72b:
call zdiffiles chk7c9 "CPU -1 $testout $refout"
chk7c9:

compare_resf chk7c1 efc efcref "Fermi" 4 0 ";"
chk7c1:
compare_resf chk7c2 sevc sevcref "sumev=" 4 0 "sumev="
chk7c2:
compare_resf chk7c3 ehkc ehkcref "LM:" 3 0 "ehk="
chk7c3:
compare_resf chk7c4 amomc amomcref "ATOM=" 6 0 mom=
chk7c4:

compare_res chk7c5 "Fermi level" $efc $efcref 1e-6 pass
chk7c5:
compare_res chk7c6 "Mag. moment" $amomc $amomcref 1e-5 pass
chk7c6:
compare_res chk7c7 "sum of evals" $sevc $sevcref 1e-6 pass
chk7c7:
compare_res chk7c8 "ehk" $ehkc $ehkcref 1e-6 pass
chk7c8:


echo " "
echo "$space ... Equivalent AFM calculation, both spins rotated by a common (random) angle (see site.$ext)"
set testout=out.lm.{$ext}.nc1 refout=$testdir/out.lm.$ext.nc1.gz
echo "$space rm -f *.$ext"
             rm -f *.$ext
echo "$space cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 ."
             cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 .
echo "$space cp $testdir/dmats.$ext.afm dmats.$ext"
             cp $testdir/dmats.$ext.afm dmats.$ext
echo "$space cp $testdir/rsta.$ext.afm rsta.$ext"
             cp $testdir/rsta.$ext.afm rsta.$ext

runrdcmd chk7r0 %11f $testout "-cat:JOBNC1 --noerr ctrl.$ext"
chk7r0:
call zdiffiles chk7r9 "CPU -1 $testout $refout"
chk7r9:

compare_resf chk7r1 efr efrref "Fermi" 4 0 ";"
chk7r1:
compare_resf chk7r2 sevr sevrref "sumev=" 4 0 "sumev="
chk7r2:
compare_resf chk7r3 ehkr ehkrref "LM:" 3 0 "ehk="
chk7r3:
compare_resf chk7r4 amomr amomrref "ATOM=" 6 0 mom=
chk7r4:
compare_resf chk7r4a Myr Myrref "<My>=" 2 0 "<My>="
chk7r4a:

compare_res chk7r5 "Fermi level" $efr $efrref 1e-6 pass
chk7r5:
compare_res chk7r6 "Mag. moment" $amomr $amomrref 1e-5 pass
chk7r6:
compare_res chk7r7 "sum of evals" $sevr $sevrref 1e-6 pass
chk7r7:
compare_res chk7r8 "ehk" $ehkr $ehkrref 1e-6 pass
chk7r8:

echo " "
echo "$space ... AFM calculation by noncollinear rotation second spin through 180 degrees"
set testout=out.lm.{$ext}.nc2 refout=$testdir/out.lm.$ext.nc2.gz
echo "$space rm -f *.$ext"
             rm -f *.$ext
echo "$space cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 ."
             cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 .
echo "$space cp $testdir/dmats.$ext.nc2 dmats.$ext"
             cp $testdir/dmats.$ext.nc2 dmats.$ext
echo "$space cp $testdir/rsta.$ext.nc2 rsta.$ext"
             cp $testdir/rsta.$ext.nc2 rsta.$ext
echo "$space cp $testdir/site.$ext.nc2 site.$ext"
             cp $testdir/site.$ext.nc2 site.$ext
runrdcmd chk7n0 %11f $testout "-cat:JOBNC1 --noerr ctrl.$ext"
chk7n0:
call zdiffiles chk7n9 "CPU -1 $testout $refout"
chk7n9:

compare_resf chk7n1 efn efnref "Fermi" 4 0 ";"
chk7n1:
compare_resf chk7n2 sevn sevnref "sumev=" 4 0 "sumev="
chk7n2:
compare_resf chk7n3 ehkn ehknref "LM:" 3 0 "ehk="
chk7n3:
compare_resf chk7n4 amomn amomnref "ATOM=" 6 0 mom=
chk7n4:
compare_resf chk7n4a Myn Mynref "<My>=" 2 0 "<My>="
chk7n4a:

if (! $?quiet) then
echo ' '
echo "$space ... The following numbers are extracted from the last iteration:"
echo "$space collinear Fermi level            = $efc"
echo "$space rotated collinear Fermi level    = $efr"
echo "$space noncollinear Fermi level         = $efn"
echo "$space reference                        = $efnref"
set diff = `echo $efn $efnref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
echo "$space collinear Er mag. moment         = $amomc"
echo "$space rotated collinear Er mag. moment = $amomr"
echo "$space noncollinear Er mag. moment      = $amomn"
echo "$space reference                        = $amomnref"
set diff = `echo $amomn $amomnref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
echo "$space rotated collinear Er <My>        = $Myr"
echo "$space noncollinear Er <My>             = $Myn"

echo ' '
echo "$space collinear sum of evals           = $sevc"
echo "$space rotated collinear sum of evals   = $sevr"
echo "$space noncollinear sum of evals        = $sevn"
echo "$space reference                        = $sevnref"
set diff = `echo $sevn $sevnref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
echo "$space collinear ehk                    = $ehkc"
echo "$space rotated collinear ehk            = $ehkr"
echo "$space noncollinear ehk                 = $ehkn"
echo "$space reference                        = $ehknref"
set diff = `echo $ehkn $ehknref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
endif

compare_res chk7n5 "Fermi level" $efn $efnref 1e-6 pass
chk7n5:
compare_res chk7n6 "Mag. moment" $amomn $amomnref 1e-5 pass
chk7n6:
compare_res chk7n7 "sum of evals" $sevn $sevnref 1e-6 pass
chk7n7:
compare_res chk7n8 "ehk" $ehkn $ehknref 1e-6 pass
chk7n8:

if ($?pass) then
    echo "$space test 7 PASSED ($ext)"
else
    echo "$space test 7 FAILED ($ext)"
    set failed = ($failed 7)
endif
chk7e:

echo ' '
if ($#failed <= 1) then
    echo "$space $testfile : all tests PASSED"
    echo " "
    exit 0
else
    shift failed
    echo "$space $testfile : These tests FAILED:" $failed
    echo " "
    exit -1
endif

# ---------------- sdjob --------------
# runs sd test for input $sdmod; returns to label eojob
exit
sdjob:
if ($?quiet) then
else
endif
set pass
set refout=$testdir/out.fccfe.3k.sdmod=$sdmod.gz testout=out.fccfe
echo ' '
query chk21 chk2e 'run this test'
chk21:
# ... Look for executables
findcmd chk21a rdcmd "$path" "$maindir"
chk21a:
findcmd chk21b lm "$path" "$topdir"
chk21b:
findcmd chk21c lmstr "$path" "$maindir"
chk21c:

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f *.fccfe"
             rm -f *.fccfe
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto $eojob
endif
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .
echo "$space cp $testdir/$eulafile eula.fccfe"
             cp $testdir/$eulafile eula.fccfe
runjob chk22 $testout "lmstr fccfe"
chk22:
runjob chk23 $testout "lm -vnit=0 -vsdyn=t --keepsignm fccfe --no-iactive"
chk23:
runjob chk24 $testout "lm -vmet=1 -vnk=5 -vnit=15 -vsdmod=$sdmod -vsdyn=t --keepsignm fccfe --no-iactive"
#asa-soft-link-str fccfe
#runjob chk24 $testout "mpix -np=9 lm-MPIK -vmet=3 -vnk=5 -vnit=15 -vsdmod=$sdmod -vsdyn=t --keepsignm fccfe --no-iactive"
chk24:

set dqn   =  `cat $testout | grep 'RMS DQ=' | grep 'file mq' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dmn   =  `cat $testout | grep 'RMS DQ=' | grep 'file ma' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | grep 'file mq' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dmnr  =  `zcat $refout | grep 'RMS DQ=' | grep 'file ma' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`

#set etest = `grep "LM: it" $testout | awk '{print substr($6,5)}' | tail -1`
#set eref = `gunzip -c $refout | grep "LM: it" | awk '{print substr($6,5)}' | tail -1`
set etest = `grep ' ehf= '  $testout | grep -v last | awk '{print $6}' | tail -1`
set eref = `gunzip -c $refout | grep ' ehf= '  | grep -v last | awk '{print $6}' | tail -1`
#  if ($sdmod >= 1000) then
#    set dq = `grep 'file ma' $testout | tail -1 | awk '{print substr($9,4)}'`
#    set etest = `grep 'ehf=' $testout | tail -1 | awk '{print substr($2,5)}'`
#    set eref = `gunzip -c $refout | grep ehf= | tail -1 | awk '{print substr($2,5)}'`
#  endif


if ($?quiet) then
else
  echo ' '
  echo "$space rms difference output-input charge dq, last iteration = $dqn"
  echo "$space rms change in euler angles dm, last iteration =         $dmn"
  echo "$space ehf,    last iteration =                               $etest"

  echo ' '
  echo "$space Compare mag. forces 1st iteration to $refout"
  cat $testout | awk '{if ($1 == "MAGTRQ:") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | head -6
  echo  ---
  gunzip -c $refout | awk '{if ($1 == "MAGTRQ:") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | head -6

  echo ' '
  echo "$space Compare last iterations of SV line to file $refout"
  cat $testout | grep SV: | tail -5
  echo "---"
  gunzip -c $refout | grep SV: | tail -5
chk28a:

  echo ' '
  echo "$space Compare last iterations of av. mag. to file $refout"
  cat $testout | grep amagnc | tail -5
  echo "---"
  gunzip -c $refout | grep amagnc | tail -5
chk28b:

  call showout chk28c CPU
chk28c:

  echo ' '
  call zdiffiles chk28d "CPU 1 $testout $refout"
chk28d:

endif
goto $eojob


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
      echo "$space Invoking rdcmd will execute the following job(s):"
      $rdcmd -f:$rdcmdfmt --n $callarg
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
  if (`echo $testvar $refvar | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjob

# ---------------- compare_res_0 --------------
# Compares a number $testvar and unsets $passvar if |testvar|<tol
# usage: compares_res_0 retcall keyword testvar tol passvar
#   keyword      : label (for printout)
#   testvar      : first number
#   tol          : tolerance
#   passvar      : $passvar is unset if |testvar|<tol
exit
compare_res_0:
  set quitjob=$retcall
#  echo $retcall $keyword $testvar $tol $passvar
 echo -n "$space $keyword ($testvar) smaller than tol ($tol)? ... "
  if (`echo $testvar 0 | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
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
# grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//"
  set $testvar = `grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  set $refvar = `zcat $refout | grep "$keyword" | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  goto $quitjob

# ---------------- cnvt_d_fmt --------------
# converts exponential format #.##D## or #.##d## to #.##E##
# usage: cnvt_d_fmt retcall testvar testval
exit
cnvt_d_fmt:
  set quitjob = $retcall
  set $testvar = `echo $testval | sed s/D/E/ | sed s/d/E/`
  goto $quitjob

# ---------------- zcmpnfiles --------------
# Compares two files, treating each field as a number.
# call arguments should contain 3 strings: no-digits test-file reference-file
# Files with .gz or .Z extensions are assumed to be gzipped.
# Example :  call zcmpnfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $testdir/tmp1 $testdir/tmp2
exit
zcmpnfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)
  set digits = $zcmpnargs[1]
  set a = ' { for (i = NF; i > 0; --i) printf " %.'$digits'f", $i; printf "\n" }'

  set fn1 = $testdir/tmp1
  set fn2 = $testdir/tmp2
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

  $cat1  $zcmpnargs[2] | awk "$a" > $fn1
  $cat2  $zcmpnargs[3] | awk "$a" > $fn2
  set ncharfile = `wc $fn1 | awk '{print $3}'`
  cmp $fn1 $fn2 >/dev/null
  set retval = $status

  if ($retval == 0) then
    rm -f $fn1 $fn2
    goto $quitjob
  endif

  set retval = `cmp -l $fn1 $fn2 |& grep -v EOF | wc | awk '{printf "%d", $1}'`
  if ($retval == 0) set retval = '-1'
  rm -f $fn1 $fn2
  goto $quitjob

# ---------------- diffiles --------------
exit
diffiles:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  set files = "$callarg"
  query diff11 $quitjob "compare $files"
diff11:
  diff $files | head -100
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
    set nend = `grep $endstr $files[1] | wc | awk '{print $1}'`
  endif

#    echo zdifffiles : $quitjob
#    grep $endstr $files[1]
  query zdiff11 $quitjob "compare $files"
zdiff11:
#    $testdir/zdiff $files | awk -v endstr="$endstr" -v nend=$nend -v endl=0 -v endr=0 '{if ($1 == "<" && endl < nend) print ; if ($1 == ">" && endr < nend) print ; if ($1 == ">" || $1 == "<" || endl >= nend && endr >= nend) ; else {print} ; if ($1 == "<" && $2 == endstr) {endl+=1}; if ($1 == ">" && $2 == endstr) {endr+=1};}' | head -50
  $testdir/zdiff $files | awk -v endstr="$endstr" -v nend=$nend -v endl=0 -v endr=0 '{if (endl < nend) print; if ($1 == ">" && $2 == endstr) {endl+=1}}' | head -50
  echo " "
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
  grep $callarg $testout
  if (`cat $testout | grep $callarg | wc | awk '{print $1}'` > 1) echo ' ---'
  gunzip -c $refout | grep $callarg
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
 usage: test.nc [switches] [file-extension] [testcase-list]
        e.g., "test.gf co 1 2"
        If file-extension is missing, test.gf uses co
        Switches:
        --no-iactive runs tests without prompting user
        --quiet runs tests with minimal output and without prompting user
        --clean runs cleanup
#       --verbose    script prints out extra information
EOF
exit -1
