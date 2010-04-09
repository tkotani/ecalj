Cgetarg...info...           structure ['asetot', 'ham', 'sham']
Cgetarg...info...           use_to_add ['asetot', 'm_struc_def']
      subroutine asetot(mode,sham,sev,etot)
       
       use m_struc_def  !Cgetarg

C- Make ASA Harris energy or Kohn-Sham energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 Harris Foulkes energy
Ci         :2 Kohn-Sham energy
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: eterms seref
Ci     Stored:    eterms
Ci   sham->eterms various integrals for the total energy are used:
Ci         :(3)  utot   = total electrostatic energy
Ci         :(6)  rhoexc = rho * exc
Ci         :(8)  sumec  = sum-of-core eigenvalues
Ci         :(10) xcore  = rhoc * total potential
Ci         :(11) valvef = rhov * total potential
Ci         :(17) rinvxt = double counting terms, input rho * appl field
Ci         :(18) rouvxt = double counting terms, output rho * appl field
Ci   sev   :sum of eigenvalues
Ci   eref  :reference energy to be subtracted from total
Co Outputs
Co   etot  :Harris energy (mode = 1)
Co         :Hohnberg-Kohn-Sham energy (mode = 2)
Co   sham->eterms various integrals for the total energy are stored:
Co         :(1 )  eh    = etot = Harris energy (mode=1)
Co         :(2 )  eks   = etot = Hohnberg-Kohn-Sham energy (mode=2)
Co         :(13)  sumev = sum of eigenvalues
Cl Local variables
Cl         :
Cr Remarks
Cr
Cr   Total energy is sum of K.E., Hartree energy, XC energy:
Cr      etot = ekin + utot + rhoexc
Cr   The kinetic energy is computed via double-counting terms
Cr     ekin = sev + sumec - rhov
Cu Updates
Cu   08 Mar 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode
Cgetarg       double precision sham(1),sev,etot
       real(8):: sev , etot 
       type(s_ham)::sham

C ... Local parameters
      integer stdo,ipr,lgunit,i,j
      double precision eterms(20),sumec,rhov,ekin,utot,rhoexc,eref,
     .  edcvxt
      character nam1*9,nam2*5

      call isanrg(mode,1,2,' asetot','mode',.true.)
      stdo = lgunit(1)
      call getpr(ipr)
C     ipr = 51

Cgetarg...info...           real(8) :: eterms(20) , real(8) :: sham%eterms 20
Cgetarg...info...           real(8) :: eref , real(8) :: sham%seref 1
Cgetarg       call upack('ham eterms seref',sham,eterms,eref,0,0,0)
       
       i_copy_size=size(sham%eterms) 
       call dcopy(i_copy_size,sham%eterms,1,eterms,1) 
       eref=sham%seref


      eterms(16) = sev

      if (mode .eq. 1) edcvxt = eterms(17)
      if (mode .eq. 2) edcvxt = eterms(18)
      sumec = eterms(8)
C     sumtc = eterms(9)
      rhov = eterms(10) + eterms(11) + edcvxt
      ekin = sev + sumec - rhov

      utot = eterms(3)
      rhoexc = eterms(6)
      etot = ekin + utot + rhoexc
      etot = etot - eref

      if (mode .eq. 1) then
        eterms(1) = etot
        nam1 = 'Harris'
        nam2 = 'ehar='
      else
        eterms(2) = etot
        nam1 = 'Kohn-Sham'
        nam2 = 'ehks='
      endif

Cgetarg...info...           real(8) :: eterms(20) , real(8) :: sham%eterms 20
Cgetarg       call pack1('ham eterms',sham,eterms)
       
       i_copy_size=size(sham%eterms) 
       call dcopy(i_copy_size,eterms,1,sham%eterms,1) 


C --- Printout ---
      if (ipr .gt. 40) then
        call word(nam1,1,i,j)
        write(stdo,310) nam1(i:j),sev
        write(stdo,311) rhov,ekin,eref,rhoexc,utot,nam2,etot
  310   format(/1x,a,' energy:':'  sumev=',f10.6)
  311   format(' rhov=',  f17.6,'   ekin=',f17.6,'   eref=',f16.6
     .        /' rhoep= ',f15.6,'   utot=',f17.6,3x,a5,     f16.6)
      endif

Cgetarg       end
       end subroutine asetot 


