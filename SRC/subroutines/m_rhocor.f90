!>core density
module m_rhocor 
  use m_lgunit,only:stdo
  use m_rseq,only: rseq
  public:: getcor,rhocor
  private
contains
  subroutine getcor(mode,z,a,pnu,pnz,nr,lmax,rofi,v0,kcor,lcor,qcor, & !- Generate the core density for one site
       sumec,sumtc,rhoc,ncore,ecore,gcore,nmcore)
    use m_lmfinit,only: nsp
    use m_lgunit,only:stdo
    !i Inputs
    !i   mode  :0 do not return ecore,gcore
    !i         :1 return core eigenvalues ecore and wave functions gcore
    !i         :2 return core charge only
    !i   z     :nuclear charge
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   pnz   :boundary conditions for (optional) second p.q.n.
    !i   nr    :number of radial mesh points
    !i   lmax  :augmentation l-cutoff
    !i   rofi  :radial mesh points
    !i   v0    :spherical potential, excluding nuclear part
    !i   kcor  :(partial core occupation) p.q.n for occupation
    !i   lcor  :(partial core occupation) l quantum for occupation
    !i   qcor  :(partial core occupation) core charge and moment
    !i          if mode=2 on input, returns core charge and exits
    !o Outputs
    !o   sumec  :sum of single-particle energies
    !o   sumtc  :core kinetic energy
    !o   rhoc   :core density
    !o   ncore  :number of core states
    !o   ecore  :(mode=1) core eigenvalues
    !o   gcore  :(mode=1) core wave functions
    !r Remarks
    !u Updates
    !u   27 Jan 07 extend to qcor<0 on input (see above)
    !u   28 Aug 01 Extended to local orbitals.  Altered argument list.
    !u   19 Apr 01 core evals and wave functions may be saved
    !u   19 Jun 00 spin polarized
    !u   26 May 00 adapted from nfp getcor.f
    ! ----------------------------------------------------------------------
    implicit none
    integer :: mode,nr,nrmx,n0,kcor,lcor,lmax,ncore,nmcore
    parameter (nrmx=1501,n0=10)
    double precision :: a,qcor(2),sumec,sumtc,z,ecore(*),gcore(nr,2,*)
    double precision :: v0(nr,nsp),rhoc(nr,nsp),rofi(nr),pnu(n0,2),pnz(n0,2)
    integer :: ipr,isp,jpr,k,kkk,l,nsc,konf(n0),konfsc(n0),isc(n0)
    double precision :: deg,qcore,qsc,tol,g(2*nrmx),ec(100),b,smec(2), smtc(2)
    character ch*2
    data ch/' *'/
    call getpr(ipr)
    qcore = 0d0
    ncore = 0
    qsc = 0d0
    nsc = 0
    do  l = 0, lmax
       deg = dble((2*(2*l+1))/nsp)
       do  isp = 1, nsp
          konfsc(l+1) = pnu(l+1,isp)
          konf(l+1)   = mod(pnz(l+1,isp),10d0)
          if (konf(l+1) == 0 .OR. konf(l+1) > konfsc(l+1)) &
               konf(l+1)= konfsc(l+1)
          do  kkk = l+1, konf(l+1)-1
             ncore = ncore+1
             if (mode /= 2) then
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
    if (ipr >= 40) write(stdo,850) qcore,qsc, (konf(k),ch(isc(k):isc(k)), k=1,lmax+1)
850 format(' getcor:  qcore=',f6.2,'  qsc=',f6.2,'  konf =',10(i2,a1))
    if (mode == 2) then
       qcor(1) = qcore
       return
    endif
    jpr = 0
    if (ipr >= 30) jpr = 1
    tol = 1d-8
    rhoc=0d0
    b = rofi(nr)/(dexp(a*(nr-1))-1)
    call rhocor(mode,z,lmax,nsp,konf,a,b,nr,rofi,v0,g,kcor,lcor,qcor,tol,ec,smec,smtc,rhoc,gcore,nmcore,jpr)
    if (mode == 1) call dcopy(ncore,ec,1,ecore,1)
    sumec = sum(smec(1:nsp))
    sumtc = sum(smtc(1:nsp))
  end subroutine getcor
  subroutine rhocor(isw,z,lmax,nsp,konfig,a,b,nr,rofi,v,g,kcor,lcor, & !Generates the (spherical) charge density from the core states
       qcor,tol,ec,sumec,sumtc,rho,gcore,nmcore,ipr) !nmcore jun2012
    !i Inputs:
    !i   isw   :1 return the core wave functions in gcore
    !i   z     :nuclear charge
    !i   lmax  :maximum l
    !i   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
    !i   konfig:core configuration. Core orbitals are specified by:
    !i         :  1, 2, ..., konf(0)-1 for s
    !i         :  2, 3, ..., konf(1)-1 for p
    !i         :  3, 4, ..., konf(2)-1 for d, and so on.
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   v     :spherical potential (electronic contribution)
    !i   g     :work array holding normalized wave function times r
    !i   kcor  :(partial core occupation) p.q.n for occupation
    !i   lcor  :(partial core occupation) l quantum for occupation
    !i   qcor  :(partial core occupation) core charge and moment
    !i   tol   :wave function tolerance
    !i   ipr   :0 no printout
    !i         :1 summary printout
    !i         :2 detailed printout
    ! o Inputs/Outputs:
    ! o  ec    :guessed core eigenvalues (input)
    ! o         core eigenvalues (output)
    !o Outputs:
    !o   rho   :spherical charge density times 4*pi*r*r
    !o         :the spherical charge density of the core states is added
    !o   sumec :sum of core eigenvalues
    !o   sumtc :core kinetic energy
    !o   gcore :(isw=0) not used
    !o         :(isw=1) core wave functions
    !l Local variables
    !l   deg   :orbital occupation number
    !r Remarks:
    !u Updates
    !u   19 Apr 01 core wave functions may be saved in gcore
    !  ---------------------------------------------------------------------
    implicit none
    integer:: isw,ipr,lmax,konfig(0:lmax),nr,nsp,kcor,lcor
    double precision :: a,b,z,g(nr,2),rho(nr,nsp),rofi(nr),tol,qcor(2), &
         sumec(nsp),sumtc(nsp),v(nr,nsp),ec(*)
    real(8),optional:: gcore(nr,2,*)
    integer,optional:: nmcore
    character(1) :: pqn(9),ang(6)
    integer :: icore,isp,konf,l,ll,nodes,nre,iprint,ir
    double precision :: deg,dlml,e1,e2,ecor0,ecore,qcore,rhorim,rmax, &
         rorim,slo,sumx,tcore,val,vrho,kappa2
    ! ... Heap
    data pqn /'1','2','3','4','5','6','7','8','9'/
    data ang /'s','p','d','f','g','h'/
    real(8),allocatable:: vaa(:,:)
    logical :: oncewrite
    !! --- spin-averaged potential ---
    !      print *,'rhocor: nmcore=',nmcore
    allocate(vaa(nr,nsp))
    if(nmcore==1 .AND. iprint()>10) then
       if(oncewrite(6))print *,'NOTE: rhocor: core density is spin-independent now for any MMOM, june2012.'
       if(oncewrite(7)) print *,'      We use spin-avaraged potential to calculate core density rhoc.'
       if(oncewrite(8)) print *,'      Thus diff of core eigen is pot. due to valence electron. I/O:iofa.F'
       do ir=1,nr
          vaa(ir,1:nsp)=sum(v(ir,1:nsp))/dble(nsp)
       enddo
    else
       vaa=v
    endif
    !! -------------------------------
    rmax  = rofi(nr)
    b     = rmax/(dexp(a*nr-a)-1.d0)
    e1    = -2.5d0*z*z - 5d0
    e2    = 20d0
    icore = 0
    do  80  isp = 1, nsp
       sumec(isp) = 0d0
       sumtc(isp) = 0d0
       qcore = 0d0
       rhorim = 0d0
       if (ipr >= 2) write(stdo,757)
       if(ipr >= 2 .AND. nsp == 2) &
            write(stdo,'('' spin'',i2,'':'')') isp
757    format(/' state  chg          ecor0',10x,'ecore',10x,'tcore', &
            4x,'nre',2x,'rho(rmax)')
       do  101  l = 0, lmax
          do  10  konf = l, konfig(l)-2
             deg = (2*(2*l+1))/nsp
             if (konf+1 == kcor .AND. l == lcor) then
                deg = deg + (qcor(1)+(3-2*isp)*qcor(2))/nsp
             endif
             icore = icore+1
             nodes = konf - l
             ecor0 = ec(icore)
             val = 1d-30
             slo = -val
             call rseq(e1,e2,ecor0,tol,z,l,nodes,val,slo,vaa(1,isp), &
                  g,sumx,a,b,rofi,nr,nre)
             ecore = ecor0
             !     ... Correct core energy by using hankel bc's
             kappa2 = ecor0 - vaa(nr,isp) + 2*z/rmax
             if (nre == nr .AND. kappa2 < 0d0) then
                dlml = -1d0 - dsqrt(-kappa2)*rmax
                do  31  ll = 1, l
                   dlml = -kappa2*rmax*rmax/dlml - (2*ll+1)
31              enddo
                slo = val*(dlml+l+1)/rmax
                call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,vaa(1,isp), &
                     g,sumx,a,b,rofi,nr,nre)
             endif
             ec(icore) = ecore
             if (isw == 1) call dcopy(2*nr,g,1,gcore(1,1,icore),1)
             !     --- Add to rho, make integral v*rho ---
             call xyrhsr(ecore,l,z,a,b,nr,nre,g,rofi,vaa(1,isp),rho(1,isp), &
                  deg,vrho,rorim)
             rhorim = rhorim + rorim
             tcore = ecore - vrho
             qcore = qcore + deg
             sumec(isp) = sumec(isp) + deg*ecore
             sumtc(isp) = sumtc(isp) + deg*tcore
             if (ipr >= 2) write(stdo,758) &
                  pqn(konf+1),ang(l+1),deg,ecor0,ecore,tcore,nre,rorim
758          format(1x,2a1,f8.2,3f15.6,i7,f9.5)
10        enddo
101    enddo
       if (ipr > 0) write(stdo,230) qcore,sumec(isp),sumtc(isp),rhorim
230    format(' sum q=',f5.2,'  sum ec=',f12.5,'  sum tc=',f12.5, &
            '  rho(rmax)',f8.5)
80  enddo
  end subroutine rhocor
  subroutine xyrhsr(ecore,l,z,a,b,nr,nre,g,rofi,v,rho,deg,vrho, rhormx)!Make density and integrate potential*density for one core state
    use m_lmfinit,only: c=>cc
    !i Inputs
    !i   ecore :core eigenvalue
    !i   l     :l quantum number
    !i   z     :nuclear charge
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   nre   :Make density to the smaller of nr and nre
    !i   g     :normalized wave function times r
    !i   rofi  :radial mesh points
    !i   v     :electronic part of spherical potential
    !i   deg   :occupation number
    !o Outputs
    !o   rho   :contribution core density from this state from 1..min(nr,nre)
    !o   vrho  :integral
    !o   rhormx:contribution to true density at nr from this state
    !r Remarks
    !r   xyrhsr makes the density at points 1..min(nr,nre), and the integral
    !r   of rho*density from 1..nre.  Thus:
    !r     if nre .eq. nr, the density and integral are the same
    !r     if nre .lt. nr  g**2 is negligible at nr (deep core state)
    !r     if nre .gt. nr  (large sphere) contribution to rho*v is included
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    integer :: nr,nre,l
    double precision :: a,b,z,ecore,g(nr,2),rofi(*),v(nr),rho(nr),rhormx
    integer :: nrmx,ir
    double precision :: fpi,fllp1,vrho,r,wgt,tmc,gfac,rhoir,deg,rmax
    fpi = 16d0*datan(1d0)
    fllp1 = l*(l+1)
    vrho = 0
    nrmx = min0(nr,nre)
    ! ... Make rho, and integrate vrho for points 1..nrmx
    rhoir=-1d99
    do  11  ir = 2, nrmx
       r = rofi(ir)
       wgt = 2*(mod(ir+1,2)+1)
       if (ir == nre) wgt = 1
       tmc = c - (v(ir) - 2*z/r - ecore)/c
       gfac = 1 + fllp1/(tmc*r)**2
       rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
       vrho = vrho + wgt*rhoir * (v(ir)-2*z/r) * (r+b)
       rho(ir) = rho(ir) + deg*rhoir
11  enddo
    rmax = rofi(nr)
    rhormx = 0
    !     nfp has the following line:
    !     if (nre .ge. nr) rhormx = deg*g(nr,1)**2/(rmax*rmax*fpi)
    if (nre >= nr) rhormx = deg*rhoir/(fpi*rmax**2)
    ! ... Integrate rho*v from nrmx+1 .. nre
    do  12  ir = nrmx+1, nre
       r = rofi(ir)
       wgt = 2*(mod(ir+1,2)+1)
       if (ir == nre) wgt = 1
       tmc = c - (v(ir) - 2*z/r - ecore)/c
       gfac = 1 + fllp1/(tmc*r)**2
       rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
       vrho = vrho + wgt*rhoir * (v(ir)-2*z/r) * (r+b)
12  enddo
    vrho = vrho*a/3
  end subroutine xyrhsr
end module m_rhocor
