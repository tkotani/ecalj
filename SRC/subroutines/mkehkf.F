      subroutine mkehkf(mode,sev,valmom,sumtv,etot)!sham,

      use m_struc_def  !Cgetarg
      use m_globalvariables
      use m_lmfinit,only: ham_eterms
C- Make Harris energy or Kohn-Sham energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 Harris Foulkes energy
Ci         :2 Kohn-Sham energy
Ci   sham  :struct containing information about the hamiltonian
Ci   sev   :band sum
Cio Inputs/Outputs
Cio  sumtv :valence kinetic energy :
Cio        :output for mode=1
Cio        :input for mode=2
Co   sham->eterms
Ci         :... the following quantities are input
Ci         :(3)  utot   = total electrostatic energy
Ci         :(4)  valves = valence rho * estat potential (not used)
Ci         :(5)  cpnves = core+nuc * estat potential (not used)
Ci         :(6)  rhoexc = rho * exc
Ci         :(7)  rhovxc = rho * vxc (not used)
Ci         :(8)  sumec  = sum of foca=0 core eigenvalues
Ci         :(9)  sumtc  = sum of core kinetic energies
Ci         :(10) xcore  = rhoc * total potential
Ci         :(11) valvef = smrhov * vsm + sum_ib valvef_ib
Ci                        valvef_ib = rhov * vtrue - smrhov * vsm)_ib
Ci         :(12) sumt0  = sum of frozen core kinetic energies
Ci         :
Co         :(15) amom   = system magnetic moment
Ci         :(19) rhosig = d.c. from self-energy sigma-vxc (eks only)
Ci                        Not used if rhosig=-99
Co         :... the following quantities are stored:
Co         :(1 )  eh    = etot = Harris energy (mode=1)
Co         :(2 )  eks   = etot = Hohnberg-Kohn-Sham energy (mode=2)
Co         :(13)  sumev = sum of eigenvalues
Co         :(15)  amom  = magnetic moment from bands (mode=1)
Co Outputs
Co   etot  :Harris energy (mode = 1)
Co         :Hohnberg-Kohn-Sham energy (mode = 2)
Cr Remarks
Cr   Information related to total energy is printed out, depending on mode
Cu Updates
Cu   02 Jan 06 stores or prints magnetic moments
Cu   11 Jan 05 double-counting term rho*sig subtracted from ehks.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode
      real(8):: sev , etot , sumtv , valmom
c      type(s_ham)::sham

C ... Local parameters
      integer stdo,stdl,ipr,ipl,nsp,nglob
      double precision eh,eks,rhoexc,sumec,utot,valvef,xcore,
     .ttcor,sumt0,sumev,sumtc,ekin,rhosig,amom
C     double precision cpnves,rhovxc,sumtc,valves
C     This ordering must match sham->eterms; see uham
      double precision eterms(20)
      equivalence (eterms(1),eh)
      equivalence (eterms(2),eks)
      equivalence (eterms(3),utot)
C     equivalence (eterms(4),valves)
C     equivalence (eterms(5),cpnves)
      equivalence (eterms(6),rhoexc)
C     equivalence (eterms(7),rhovxc)
      equivalence (eterms(8),sumec)
      equivalence (eterms(9),sumtc)
      equivalence (eterms(10),xcore)
      equivalence (eterms(11),valvef)
      equivalence (eterms(12),sumt0)
      equivalence (eterms(13),sumev)
      equivalence (eterms(15),amom)
      equivalence (eterms(19),rhosig)

C     stdo = lgunit(1)
C     stdl = lgunit(2)
Changenglob      nsp  = nglob('nsp')
      nsp  = globalvariables%nsp
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
Changenglob      stdl = nglob('stdl')
      stdl = globalvariables%stdl
      call getpr(ipr)
      ipl = 1

c      i_copy_size=size(ham_eterms)
c      call dcopy(i_copy_size,ham_eterms,1,eterms,1)
      eterms=ham_eterms
      sumev = sev

      if (mode .eq. 1) then
        sumtv = sev - valvef
        ttcor = sumec - xcore + sumt0
        ekin  = sumtv + ttcor
        etot  = ekin + utot + rhoexc
        eh    = etot
        amom  = valmom
ccccccccccccccc
c      print *,'vvv ekin=',ekin,sumtv,ttcor
ccccccccccccccc
      elseif (mode .eq. 2) then
C   ... double counting from sigma
        if (rhosig .ne. -99) then
          sumtv = sumtv - rhosig
        endif
        ekin = sumtv + sumtc
        etot = ekin + utot + rhoexc
        eks  = etot
ccccccccccccccc
c      print *,'vvv ekin=',ekin,sumtv,sumtc
ccccccccccccccc
      endif

      i_copy_size=size(ham_eterms)
      call dcopy(i_copy_size,eterms,1,ham_eterms,1)

      if (mode .eq. 0) return

C      print *, sumev(1,1)-sev
C      print *, 'xxsev=',sev
C      print *, 'xxvalvef*-1=',-valvef
C      print *, 'xxttcor=',ttcor
C      print *, 'xxutot=',utot
C      print *, 'xxrhoexc=',rhoexc

C --- Printout ---
      if (mode .eq. 1) then
        if (ipr .ge. 30) then
          write(stdo,660) sev,valvef,sumtv,sumec,xcore,ttcor,rhoexc,utot,
     .    eh
  660     format(/' Harris energy:'
     .    /' sumev= ',f15.6,'  val*vef=',f15.6,'   sumtv=',f15.6
     .    /' sumec= ',f15.6,'  cor*vef=',f15.6,'   ttcor=',f15.6
     .    /' rhoeps=',f15.6,'     utot=',f15.6,'    ehar=',f15.6)
C 660 format(/' sumev= ',f17.8,'  val*vef=',f17.8,'   sumtv=',f17.8
C    .       /' sumec= ',f17.8,'  cor*vef=',f17.8,'   ttcor=',f17.8
C    .       /' rhoeps=',f17.8,'     utot=',f17.8,'    etot=',f17.8)
        elseif (ipr .ge. 10) then
          call awrit4('%N ekin=%,6;6d  rho*v=%,6;6d  sumev=%,6;6d'//
     .    ' ehf =%,6;6d',' ',80,stdo,ekin,utot+rhoexc,sev,eh)
        endif
        if (ipr.gt.0 .and. ipl.gt.0) write (stdl,720) utot,ekin,sev,etot
  720   format('fp Har: U',f15.6,'  T',f15.6,'  sev',f14.6,'  EH ',f14.6)

      elseif (mode .eq. 2) then
        if (ipr .ge. 30) then
          if (rhosig .ne. -99 .and. rhosig .ne. 0) then
            write (stdo,411) sumtv,sumtc,ekin,rhoexc,utot,rhosig,eks
  411       format(/' Kohn-Sham energy:'
     .      /' sumtv= ',f15.6,'  sumtc=',f17.6,'   ekin=',f16.6
     .      /' rhoep= ',f15.6,'   utot=',f17.6,' rhosig=',f16.6
     .      /'  ehks=', f16.6)
          else
            write (stdo,410) sumtv,sumtc,ekin,rhoexc,utot,eks
  410       format(/' Kohn-Sham energy:'
     .      /' sumtv= ',f15.6,'  sumtc=',f17.6,'   ekin=',f16.6
     .      /' rhoep= ',f15.6,'   utot=',f17.6,'   ehks=',f16.6)
          endif
          if (nsp .eq. 2) write (stdo,412) valmom,amom
  412     format(' mag. mom=',f13.6,'  (bands)',f16.6,'  (output rho)')

        elseif (ipr .ge. 10) then
          call awrit4('%N ekin=%,6;6d  rho*v=%,6;6d ehf =%,6;6d'//
     .    '  ehks =%,6;6d',' ',80,stdo,ekin,utot+rhoexc,eh,eks)
        endif
        if (ipr.gt.0 .and. ipl.gt.0) write (stdl,721) utot,ekin,rhoexc,eks
  721   format('fp KS: U',f16.7,' T',f16.7,' Exc',f15.7,' EKS',f15.7)

      endif

C      ifet=1456
C      open(ifet,file='ETOTLDA')
C      write(ifet,"(d23.16,a)")  eks,    ' ! EKS  (Ry)'
C      write(ifet,"(d23.16,a)")  rhoexc, ' ! \int rho exc '
C      write(ifet,"(d23.16,a)")  utot,   ' ! U '
C      write(ifet,"(d23.16,a)")  ekin,   ' ! T '
C      close(ifet)

      end subroutine mkehkf


