subroutine rhogkl ( ib1 , ib2 , nsp , mode , sspec , sv_p_orhoat , kmax , qkl )
  use m_lgunit,only:stdo
  use m_struc_def  !Cgetarg
  use m_lmfinit,only: ispec
  !- G_kL expansion of valence sphere densities
  ! ----------------------------------------------------------------------
  !i Inputs
  !i  ib1,ib2: compute expansion coefficents for sites ib1..ib2
  !i   nsp   :1 make qkl for first spin (possibly the only one)
  !i         :2 make qkl combining spins 1 and 2
  !i   mode  : a compound of digits specifying what is to be included
  !i         : in the expansion coefficients
  !i         : 1s   digit = 1 include local density rho1-rho2
  !i         :              2 include local density rho1
  !i         :              3 include local density rho2
  !i         : 10s  digit = 1 include core density rhoc
  !i                        2 include -1 * core density from sm-hankel
  !i                        3 combination 1+2
  !i         : 100s digit = 1 add -1 * nuclear density Y0 Z delta(r)
  !i   kmax  :make expansion coffs to polynomial cutoff kmax
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxl z qc a nr rmt rg
  !i     Stored:    *
  !i     Passed to: corprm
  !i   orhoat:vector of offsets containing site density
  !o Outputs
  !o   qkl  :Expansion coefficients, stored as a single long vector.
  !o        := integral pkl Y_L integrand
  !o        :where integrand is according to mode
  !r Remarks
  !r    In the spin-polarized case, up- and down- spins are combined.
  !u Updates
  !u   19 Oct 01 Adapted from rhomom.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer:: ib1 , ib2 , nsp , mode , kmax
  type(s_rv1) :: sv_p_orhoat(3,1)

  real(8):: qkl(0:kmax,1)
!  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)

  ! ... Local parameters
  integer:: ipr , iprint, nrmx , j1 , ib , is &
       , igetss , lmxl , nr , nlml , ilm , j , lfoc , k , l , m
  real(8) ,allocatable :: rofi_rv(:)
  real(8) ,allocatable :: rwgt_rv(:)
  real(8) ,allocatable :: h_rv(:)

  parameter( nrmx=1501)
  double precision :: z,qc,a,rmt,qcorg,qcorh,qsc,cofg,cofh,rg, &
       ceh,rfoc,df(0:20)
  ! ... Heap

  ! --- Setup ---
  ipr  = iprint()
  !      stdo = lgunit(1)
  allocate(rofi_rv(nrmx))

  allocate(rwgt_rv(nrmx))

  call stdfac(20,df)
  if (ipr >= 40) write(stdo,221)

  ! --- Loop over sites ---
  j1 = 1
  do  ib = ib1, ib2
     is = ispec(ib)! int(ssite(ib)%spec)
     lmxl=sspec(is)%lmxl
     z=sspec(is)%z
     qc=sspec(is)%qc
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     rg=sspec(is)%rg

     if (lmxl == -1) goto 10
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     qc = qcorg+qcorh
     nlml = (lmxl+1)**2
     call radmsh ( rmt , a , nr , rofi_rv )
     call radwgt ( rmt , a , nr , rwgt_rv )
     allocate(h_rv(nr*(kmax+1)*(lmxl+1)))
     call pvrgkl ( mode , kmax , nlml , nr , nsp , rofi_rv , rwgt_rv &
          , sv_p_orhoat( 1 , ib )%v , sv_p_orhoat( 2 , ib )%v , sv_p_orhoat( 3 , ib )%v &
          , h_rv , cofh , rg , ceh , rfoc , z , qkl ( 0 , j1 ) )
     deallocate(h_rv)

     if (ipr >= 40) then
        write(stdo,222) ib,0,1,(qkl(k,j1), k=0,kmax)
        ilm = 1
        do  l = 1, lmxl
           do  m = -l, l
              ilm = ilm+1
              j = j1+ilm-1
              if (dabs(qkl(0,j))*df(2*l+1) > 1d-6) write(stdo,220) 0,ilm, &
                   (qkl(k,j)*df(2*l+1),k=0,kmax)

           enddo
        enddo
     endif
222  format(2x,'ib=',i3,i5,i6,10f12.6)
220  format(9x,i4,i6,f12.6,10f12.6)
221  format(/' rhogkl:    k   ilm      qkl (2l+1)!! ...')
     !$$$#if DEBUG
     !$$$        if (ib .eq. 1) then
     !$$$          print *, 'ib=',ib
     !$$$          call prtrkl ( mode , kmax , rg , nr , nlml , nsp , rofi_rv ,
     !$$$     .     sv_p_orhoat( 1 , ib )%v , sv_p_orhoat( 2 , ib )%v , sv_p_orhoat( 3 , ib )%v
     !$$$     .     , qkl ( 0 , j1 ) )
     !$$$
     !$$$
     !$$$        endif
     !$$$#endif
     j1 = j1+nlml
10   continue
  enddo

  if (allocated(rwgt_rv)) deallocate(rwgt_rv)
  if (allocated(rofi_rv)) deallocate(rofi_rv)

end subroutine rhogkl


subroutine pvrgkl(mode,kmax,nlml,nr,nsp,rofi,rwgt,rho1,rho2,rhoc, pkl,cofh,rg,ceh,rfoc,z,qkl)
  use m_hansr,only:hansmr
  !- Multipole moments for one site
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  : a compound of digits specifying what is to be included
  !i         : in the expansion coefficients
  !i         : 1s   digit = 1 include local density rho1-rho2
  !i         :              2 include local density rho1
  !i         :              3 include local density rho2
  !i         : 10s  digit = 1 include core density rhoc
  !i                        2 include -1 * core density from sm-hankel
  !i                        3 combination 1+2
  !i         : 100s digit = 1 add -1 * nuclear density Y0 Z delta(r)
  !i   kmax  :k-cutoff for polynomial expansion of radial part
  !i   nlml  :L-cutoff for charge
  !i   nr    :number of radial mesh points
  !i   nsp   :number of spins
  !i   rofi  :radial mesh points
  !i   rwgt  :radial integration weights
  !i   rho1  :local true density*r**2, tabulated on a radial mesh
  !i   rho2  :local smoothed density*r**2, tabulated on a radial mesh
  !i   rhoc  :core density
  !i   cofh  :coefficient to Hankel part of pseudocore density (corprm)
  !i   rg    :smoothing radius for compensating gaussians
  !i   ceh   :energy of hankel function to fit core tail
  !i   rfoc  :smoothing radius for hankel head fitted to core tail
  !i   z     :nuclear charge
  !o Outputs
  !o   qkl  :expansion coefficients for rho
  !w Workarea:
  !w   pkl:
  !r Remarks
  !r   Q_kL = integral p_kl (rho1-rho2) + l=0 contr. from core spillout
  !r   The core spillout term is:
  !r      qcore(rhoc)-z  - sm_qcore-sm_qnuc
  !r   pvrgkl makes this Q_kL when mode=131; partial contr for other modes
  !r   NB: p0l = a**l and scaling factor for k=0 is 4*pi/(a**l * (2l+1)!!)
  !r       => q0l = 4*pi/(2l+1)!! q_l, where q_l is the multipole moment
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,kmax,nlml,nr,nsp
  double precision :: ceh,cofh,rfoc,rg,z
  double precision :: rofi(1),rwgt(1),qkl(0:kmax,nlml), &
       rhoc(nr,nsp),rho1(nr,nlml,nsp),rho2(nr,nlml,nsp), &
       pkl(nr,0:kmax,0:*)
  ! ... Local parameters
  integer :: n0,i,ilm,l,m,lmxl,ll,isp,k
  parameter (n0=10)
  double precision :: ag,fac,y0,xi(0:n0),fpi,factk,dfact, &
       df(0:20),wk(nr),smrch,f1,f2
  !     double precision gkl(0:kmax,0:nlml),wk2(0:kmax,0:nlml),rl
  !     double precision sumh,samh,y0

  !     call prrmsh('rho1',rofi,rho1,nr,nr,nlml*nsp)
  !     call prrmsh('rho2',rofi,rho2,nr,nr,nlml*nsp)

  fpi  = 16d0*datan(1d0)
  y0   = 1d0/dsqrt(fpi)
  lmxl = ll(nlml)
  call stdfac(20,df)
  call vecpkl(rofi,rg,nr,kmax,lmxl,nr,kmax,wk,1,pkl,pkl)
  !     Non-vectorized form ... should be able to integrate with
  !     either pkl -> gkl exp, or with gkl -> pkl exp, but
  !     something is wrong ...  doesn't work
  !      do  i = 1, nr
  !        call radgkl(rofi(i),rg,kmax,lmxl,kmax,wk2)
  !C       call radpkl(rofi(i),rg,kmax,lmxl,kmax,wk2)
  !        rl = 1
  !        do  l = 0, lmxl
  !          pkl(i,0:kmax,l) = wk2(0:kmax,l)*rl
  !          rl = rl * rofi(i)
  !        enddo
  !      enddo
  call dpzero(qkl,nlml*(kmax+1))

  ! ... rho1-rho2 contribution (or rho1, or rho2, depending on mode)
  if (mod(mode,10) > 0) then
     if (mod(mode,10) == 1) then
        f1 = 1
        f2 = -1
     elseif (mod(mode,10) == 2) then
        f1 = 1
        f2 = 0
     elseif (mod(mode,10) == 3) then
        f1 = 0
        f2 = 1
     else
        call rx('rhogkl: bad mode')
     endif
     ilm = 0
     do  l = 0, lmxl
        do  m = -l, l
           ilm = ilm+1
           do  k = 0, kmax
              do  i = 1, nr
                 !             If rg is small enough, these should all integrate to 1
                 !             call radgkl(rofi(i),rg,kmax,lmxl,kmax,gkl)
                 !              if (m .eq. -l) qkl(k,ilm) = qkl(k,ilm) +
                 !     .          rwgt(i)*rofi(i)**(2+l) * gkl(k,l) * pkl(i,k,l)
                 do  isp = 1, nsp
                    qkl(k,ilm) = qkl(k,ilm) + rwgt(i) * pkl(i,k,l) * &
                         (f1*rho1(i,ilm,isp) + f2*rho2(i,ilm,isp))
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  ! ... Core part (spec'd by 10s digit mode)
  if (mod(mode/10,10) > 0) then
     do  k = 0, kmax
        do  isp = 1, nsp
           !         Case 1 or 3: add core density
           if (mod(mod(mode/10,10),2) /= 0) then
              do  i = 1, nr
                 qkl(k,1) = qkl(k,1) + y0*rwgt(i)*rhoc(i,isp)*pkl(i,k,0)
              enddo
           endif
           !         Case 2 or 3: subtract core density from sm. Hankel
           if (mod(mode/10,10) >= 2) then
              do  i = 1, nr
                 call hansmr(rofi(i),ceh,1/rfoc,xi,1)
                 smrch = cofh*xi(0)*rofi(i)**2
                 qkl(k,1) = qkl(k,1) - rwgt(i)*smrch*pkl(i,k,0)
              enddo
           endif
        enddo
     enddo
  endif

  ! ... Nuclear part (spec'd by 100s digit mode)
  if (mod(mode/100,10) == 1) then
     do  k = 0, kmax
        qkl(k,1) = qkl(k,1) - y0*z*pkl(1,k,0)
     enddo
  endif

  ! ... Scale to get coefficients of the G_kL; see radpkl
  ag = 1/rg
  ilm = 0
  dfact = 1
  do  l = 0, lmxl
     do  m = -l, l
        ilm = ilm+1
        factk = 1d0
        do  k = 0, kmax
           fac = fpi / ((4*ag*ag)**k * ag**l * factk * dfact)
           qkl(k,ilm) = qkl(k,ilm) * fac
           factk = factk*(k+1)
        enddo
     enddo
     dfact = dfact*(2*l+3)
  enddo

end subroutine pvrgkl

!$$$#if DEBUG
!$$$      subroutine prtrkl(opt,kmax,rg,nr,nlml,nsp,rofi,rho1,rho2,rhoc,qkl)
!$$$
!$$$C- Printout pkl expansion of rho, for debugging
!$$$C  Example for integrating tabulated moment:
!$$$C  mc out.te -e4 x1 x2 'x5*x1' 'x7*x1*x1' -int 0 2.113465
!$$$      implicit none
!$$$      integer opt,kmax,nr,nlml,nsp
!$$$      double precision rg,rofi(nr),rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),
!$$$     .qkl(0:kmax,nlml),rhoc(nr,nsp)
!$$$      integer lmxx,n0
!$$$      parameter(lmxx=6)
!$$$      double precision g(0:kmax,0:lmxx)
!$$$      double precision ,allocatable:: p(:,:,:), rhop(:,:,:)
!$$$      integer lmxl,ll,isp,ilm,ir,k,l,m
!$$$      double precision wk(nr)
!$$$
!$$$C     call prrmsh('rhoc',rofi,rhoc,nr,nr,1)
!$$$      if (mod(opt,10) .eq. 1) then
!$$$        call daxpy(nr*nlml*nsp,-1d0,rho2,1,rho1,1)
!$$$      endif
!$$$      call prrmsh('given rho1-rho2',rofi,rho1,nr,nr,nlml)
!$$$      lmxl = ll(nlml)
!$$$
!$$$      allocate (p(nr,0:kmax,0:lmxl),rhop(nr,nlml,nsp))
!$$$      rhop = 0
!$$$      isp = 1
!$$$      do  ir = 1, nr
!$$$C       call radpkl(rofi(ir),rg,kmax,lmxl,kmax,g)
!$$$        call radgkl(rofi(ir),rg,kmax,lmxl,kmax,g)
!$$$        ilm = 0
!$$$        do  l = 0, lmxl
!$$$          do  m = -l, l
!$$$            ilm = ilm+1
!$$$            do  k = 0, kmax
!$$$              rhop(ir,ilm,1) = rhop(ir,ilm,1) +
!$$$     .        qkl(k,ilm) * g(k,l) * rofi(ir)**(2+l)
!$$$            enddo
!$$$          enddo
!$$$        enddo
!$$$      enddo
!$$$      call prrmsh('fit rho',rofi,rhop,nr,nr,nlml)
!$$$      deallocate (p,rhop)
!$$$      if (mod(opt,10) .eq. 1) then
!$$$        call daxpy(nr*nlml*nsp,1d0,rho2,1,rho1,1)
!$$$      endif
!$$$      end subroutine prtrkl
!$$$
!$$$#endif

