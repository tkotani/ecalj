      subroutine getcor(mode,z,a,pnu,pnz,nr,lmax,rofi,v0,kcor,lcor,qcor,
     .sumec,sumtc,rhoc,ncore,ecore,gcore,nmcore)
      use m_globalvariables
C- Generate the core density for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do not return ecore,gcore
Ci         :1 return core eigenvalues ecore and wave functions gcore
Ci         :2 return core charge only
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnz   :boundary conditions for (optional) second p.q.n.
Ci   nr    :number of radial mesh points
Ci   lmax  :augmentation l-cutoff
Ci   rofi  :radial mesh points
Ci   v0    :spherical potential, excluding nuclear part
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci          if mode=2 on input, returns core charge and exits
Co Outputs
Co   sumec  :sum of single-particle energies
Co   sumtc  :core kinetic energy
Co   rhoc   :core density
Co   ncore  :number of core states
Co   ecore  :(mode=1) core eigenvalues
Co   gcore  :(mode=1) core wave functions
Cr Remarks
Cu Updates
Cu   27 Jan 07 extend to qcor<0 on input (see above)
Cu   28 Aug 01 Extended to local orbitals.  Altered argument list.
Cu   19 Apr 01 core evals and wave functions may be saved
Cu   19 Jun 00 spin polarized
Cu   26 May 00 adapted from nfp getcor.f
C ----------------------------------------------------------------------
      use m_rhocor
C     implicit none
C ... Passed parameters
      integer mode,nr,nrmx,n0,kcor,lcor,lmax,ncore
      parameter (nrmx=1501,n0=10)
      double precision a,qcor(2),sumec,sumtc,z,ecore(*),gcore(nr,2,*)
      double precision v0(nr),rhoc(nr),rofi(nr),pnu(n0,2),pnz(n0,2)
C ... Local parameters
      integer stdo,nglob,ipr,isp,jpr,k,kkk,l,lgunit,nsp,nsc,
     .konf(n0),konfsc(n0),isc(n0)
      double precision deg,qcore,qsc,tol,g(2*nrmx),ec(100),b,smec(2),
     .smtc(2)
      character ch*2
      data ch/' *'/

      stdo = lgunit(1)
      call getpr(ipr)
Changenglob      nsp = nglob('nsp')
      nsp = globalvariables%nsp
      qcore = 0d0
      ncore = 0
      qsc = 0d0
      nsc = 0
      do  l = 0, lmax
        deg = dble((2*(2*l+1))/nsp)
        do  isp = 1, nsp
          konfsc(l+1) = pnu(l+1,isp)
          konf(l+1)   = mod(pnz(l+1,isp),10d0)
          if (konf(l+1).eq.0 .or. konf(l+1).gt.konfsc(l+1))
     .    konf(l+1)= konfsc(l+1)
          do  kkk = l+1, konf(l+1)-1
            ncore = ncore+1
            if (mode .ne. 2) then
              ec(ncore) = -5d0
            endif
            qcore = qcore+deg
          enddo
          do  kkk = l+1, konfsc(l+1)-1
            nsc = nsc+1
            qsc = qsc+deg
          enddo
        enddo
        isc(l+1) = 1 + konfsc(l+1)-konf(l+1)
      enddo

      if (ipr .ge. 40) write(stdo,850) qcore,qsc,
     .(konf(k),ch(isc(k):isc(k)), k=1,lmax+1)
  850 format(' getcor:  qcore=',f6.2,'  qsc=',f6.2,'  konf =',10(i2,a1))

      if (mode .eq. 2) then
        qcor(1) = qcore
        return
      endif

      jpr = 0
      if (ipr .ge. 30) jpr = 1
      if (ipr .ge. 45) jpr = 2
      tol = 1d-8
      call dpzero(rhoc,   nr*nsp)
      b = rofi(nr)/(dexp(a*(nr-1))-1)
      call rhocor(mode,z,lmax,nsp,konf,a,b,nr,rofi,v0,g,kcor,lcor,
     .qcor,tol,ec,smec,smtc,rhoc,gcore,nmcore,jpr)
      if (mode .eq. 1) call dcopy(ncore,ec,1,ecore,1)

      sumec = smec(1)
      sumtc = smtc(1)
      if (nsp .eq. 2) then
        sumec = smec(1)+smec(2)
        sumtc = smtc(1)+smtc(2)
      endif

      end

