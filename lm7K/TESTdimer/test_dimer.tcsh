#!/bin/tcsh

foreach ext ( fe2 ge2 li2 mn2 mg2 co2 sc2 zn2 kr2 ne2 ni2 cr2 na2 ti2 ar2 kr2 ne2 ni2 cr2 na2 ti2 ar2)
set da=0

if( $ext == 'al2' ) then
  set elt=Al
  set mmom='MMOM=0 1 0 0'
  set nspin=2
  set disang=1.246237
else if( $ext == 'fe2') then
  set elt=Fe
  set mmom='MMOM=0 0 4 0'
  set nspin=2
  set disang=1.005793
else if( $ext == 'ge2') then
  set elt=Ge
  set mmom='MMOM=0 2 0 0'
  set nspin=2
  set disang=1.211273
else if( $ext == 'li2') then
  set elt=Li
  set mmom='MMOM=0 0 0 0'
  set nspin=1
  set disang=1.367254
else if( $ext == 'mn2') then
  set elt=Mn
  set mmom='MMOM=0 0 2 0'
  set nspin=2
  set disang=0.824049
else if( $ext == 'mg2') then
  set elt=Mg
  set mmom='MMOM=1 0 0 0'
  set nspin=2
  set disang=1.751421
else if( $ext == 'co2') then
  set elt=Co
  set mmom='MMOM=0 0 4 0'
  set nspin=2
  set disang=1.216687
else if( $ext == 'sc2') then
  set elt=Sc
  set mmom='MMOM=0 0 2 0'
  set nspin=2
  set disang=1.156531
else if( $ext == 'zn2') then
  set elt=Zn
  set mmom='MMOM=0 0 0 0'
  set nspin=1
  set disang=1.604071
else if( $ext == 'kr2' ) then
  set elt=Kr
  set mmom='MMOM=0 0 0 0'
  set nspin=1
  set disang=1.878011
else if( $ext == 'ne2') then
  set elt=Ne
  set mmom='MMOM=0 0 0 0'
  set nspin=1
  set disang=1.275151
else if( $ext == 'ni2') then
  set elt=Ni
  set mmom='MMOM=0 0 2 0'
  set nspin=2
  set disang=1.063304
else if( $ext == 'cr2') then
  set elt=Cr
  set mmom='MMOM=0 0 3 0'
  set nspin=2
  set disang=0.797734
else if( $ext == 'na2') then
  set elt=Na
  set mmom='MMOM=0 0 0 0'
  set nspin=1
  set disang=1.545153
else if( $ext == 'ti2') then
  set elt=Ti
  set mmom='MMOM=0 0 2 0'
  set nspin=2
  set disang=0.948721
else if( $ext == 'ar2') then
  set elt=Ar
  set mmom='MMOM=0 0 0 0'
  set nspin=2
  set disang=1.719724
endif

pushd .
mkdir $ext
cd $ext
cat >ctrls.$ext <<EOF
% const len_au=0.529177
% const disexp=$disang*2/len_au a=10.0/len_au dis=0.0  dd=(disexp+dis)/a
% const da=$da 
STRUC   ALAT={a} DALAT={da} PLAT= 1 0 0  0 1 0   0 0 1
SITE    ATOM=$elt POS= 0 0 {dd}*.5
        ATOM=$elt POS= 0 0 -{dd}*.5
EOF
ctrlgen.py $ext |tee out lctrlgen
sed -e "s/MMOM=0 0 0 0/$mmom/g" \
    -e 's/\#TETRA=0/TETRA=0/g'  \
    -e 's/\#N=-1/N=-1/g'  \
    -e 's/\#W=0.01/W=0.01/g'  \
    -e "s/nspin=1/nspin=$nspin/g"  \
    -e 's/nk=2/nk=1/g'  \
    -e 's/nit=30/nit=40/g'  \
    ctrlgen.ctrl.$ext  > ctrl.$ext
lmfa $ext >llmfa 
lmf  $ext >llmf
popd

end 
