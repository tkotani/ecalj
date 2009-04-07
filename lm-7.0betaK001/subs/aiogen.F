      logical function aiogen(alabl,z,rmax,lmx,nsp,lrel,nr,a,
     . qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
C- File I/O for atomic general data.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Cio  alabl :class label
Cio  z     :nuclear charge
Cio  rmax  :augmentation radius, in a.u.,
Cio  lmx   :lmx(j) = maximum l for atom j
Cio  nsp   :2 for spin-polarized case, otherwise 1
Cio  lrel  :0 for non-relativistic, 1 for scalar relativistic
Cio  nr    :number of radial mesh points
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio  qc    :sphere core charge
Cio  dq    :dq(ic) = charge on atom of class ic
Cio  vrmax :average potential at wsr
Cio  sumec :sum of core eigenvalues
Cio  sumev :sum of valence eigenvalues, calculated from energy moments
Cio  thrpv :3PV
Cio  ekin  :sphere kinetic energy
Cio  utot  :sphere electrostatic energy
Cio  rhoeps:sphere XC energy
Cio  etot  :sphere total energy, calculated from energy moments
Cio  ifi   :file logical unit, but >0 for read, <0 for write
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   14 Mar 03 Save rhoeps to 7 decimal places
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lmx,nsp,nr,ifi
      double precision z,rmax,a,
     . qc,dq,vrmax(2),sumec,sumev,thrpv,ekin,utot,rhoeps,etot
      character*8 alabl
      logical lrel
C Local parameters
      logical scat

      aiogen = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'GEN:',':',.true.)) return
        read(ifi,200) alabl,lmx,nsp,rmax,lrel,nr,a
        read(ifi,201) z,qc,dq,vrmax,rhoeps
        read(ifi,202) sumec,utot,ekin,sumev,etot,thrpv
        aiogen = .true.
      else
        write(-ifi,'(''GEN:'')')
        write(-ifi,100) alabl,lmx,nsp,rmax,lrel,nr,a
        write(-ifi,101) int(z),int(qc),dq,vrmax,rhoeps
        write(-ifi,102) sumec,utot,ekin,sumev,etot,thrpv
      endif
  100 format(3X,A4,'  LMX=',I1,'  NSPIN=',I1,'  RMAX=',F9.6,'  REL=',L1,
     .       '  NR=',I4,'  A=',F5.3)
  101 format(3X,'Z=',I2,  '  QC=',I2,'  QTOT=',F9.6,'  VRMAX=',2F9.6,
     .          '  RHOEPS=',F13.7)
  102 format('   SUMEC=',F15.7,'  UTOT=',F15.7,'  EKIN= ',F15.7/
     .       '   SUMEV=',F15.7,'  ETOT=',F15.7,'  THRPV=',F15.7)

  200 format(3X,A4,6X,I1,8X,I1,7X,F9.6,6X,L1,5X,I4,4X,F5.3)
  201 format(3X,2X,F2.0,5X,F2.0,7X,F9.6,8X,2F9.6,9X,F12.6)
  202 format(2(9X,F15.7,7X,F15.7,8X,F15.7/))

      end

