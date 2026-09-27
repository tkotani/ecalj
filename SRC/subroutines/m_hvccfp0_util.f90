module m_hvccfp0_util
  public mkb0,strxq,strxq_all
  private
contains
subroutine mkb0( q, lxx,lx,nxx,nx, aa,bb, nrr,nrx,rprodx, alat,bas,nbas,nbloch, b0mat)
  !--make the matrix elementes < B_q | exp(iq r)>
  use m_ll,only: ll
  implicit none
  integer :: nlx,l,n,m,nr,ir,lm,ibl1,ibas,nrx,nbloch, nbas,lxx, lx(nbas), nxx, nx(0:lxx,nbas),nrr(nbas)
  real(8)    :: rprodx(nrx,nxx,0:lxx,nbas),aa(nbas),bb(nbas), phi(0:lxx),psi(0:lxx), bas(3,nbas), &
       alat, pi,fpi,tpiba,qg1(3),q(3),absqg,r2s,a,b
  complex(8) :: b0mat(nbloch),img=(0d0,1d0) ,phase
  integer(4),allocatable:: ibasbl(:), nbl(:), lbl(:), lmbl(:)
  real(8),allocatable :: ajr(:,:),rofi(:),rob0(:,:,:),cy(:),yl(:)
  complex(8),allocatable :: pjyl(:,:)
  write(6,*)'mkb0:'
  pi   = 4d0*datan(1d0)
  fpi  = 4*pi
  nlx  = (lxx+1)**2
  tpiba = 2*pi/alat
  qg1(1:3) = tpiba * q(1:3)
  absqg    = sqrt(sum(qg1(1:3)**2))
  allocate(ajr(1:nrx,0:lxx), pjyl(nlx,nbas),rofi(nrx), &
       ibasbl(nbloch), nbl(nbloch), lbl(nbloch), lmbl(nbloch), &
       cy(nlx),yl(nlx),rob0(nxx,0:lxx,nbas))
  call sylmnc(cy,lxx)
  call sylm( qg1/absqg,yl,lxx,r2s) !spherical factor Y( q+G )
  do ibas = 1,nbas
     a = aa(ibas)
     b = bb(ibas)
     nr= nrr(ibas)
     rofi(1)    = 0d0
     do ir      = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
        call bessl(absqg**2*rofi(ir)**2,lx(ibas),phi,psi)
        do l  = 0,lx(ibas)
           ajr(ir,l) = phi(l)* rofi(ir) **(l +1 ) !  ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l , where j_l is the Bessel function
        enddo
     enddo
     ! ... Coefficients for j_l yl  on MT  in the expantion of of exp(i q r).
     phase = exp( img*sum(qg1(1:3)*bas(1:3,ibas))*alat  )
     do lm = 1,(lx(ibas)+1)**2
        l = ll(lm)
        pjyl(lm,ibas) = fpi *img**l *cy(lm)*yl(lm) *phase  *absqg**l
     enddo
     ! ... rob0
     do l = 0,lx(ibas)
        do n = 1,nx(l,ibas)
           call gintxx( ajr(1,l), rprodx(1,n,l,ibas), a,b,nr, rob0(n,l,ibas) )
        enddo
     enddo
  enddo
  ! ... index (mx,nx,lx,ibas) order.
  ibl1 = 0
  do ibas= 1, nbas
     do l   = 0, lx(ibas) ! write(6,'(" l ibas nx =",3i5)') l,nx(l,ibas),ibas
        do n   = 1, nx(l,ibas)
           do m   = -l, l
              ibl1  = ibl1 + 1
              ibasbl(ibl1) = ibas
              nbl   (ibl1) = n
              lbl   (ibl1) = l
              lmbl  (ibl1) = l**2 + l+1 +m ! write(6,*)ibl1,n,l,m,lmbl(ibl1)
           enddo
        enddo
     enddo
  enddo
  ! ... pjyl * rob0
  do ibl1= 1, nbloch
     ibas= ibasbl(ibl1)
     n   = nbl  (ibl1)
     l   = lbl  (ibl1)
     lm  = lmbl (ibl1)
     b0mat(ibl1) = pjyl(lm,ibas) * rob0(n,l,ibas)
  enddo
  deallocate(ajr, pjyl,rofi, ibasbl, nbl, lbl, lmbl, cy,yl,rob0)
end subroutine mkb0
subroutine strxq(mode,e,q,p,nlma,nlmh,ndim,alat,vol,awald,nkd,nkq,dlv,qlv,cg,indxcg,jcg, s,sd)
  use m_ll,only: ll
  use m_hamindex,only:   plat,qlat
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  use m_hsmq,only: hsmq,hsmqe0
  !- One-center expansion coefficents to j of Bloch summed h (strux)
  ! ----------------------------------------------------------------
  !r  Onsite contribution is not contained in the bloch sum in the case of p=0. See job=1 for hsmq.
  !i Inputs:
  !i   mode  :1's digit (not implemented)
  !i         :1: calculate s only
  !i         :2: calculate sd only
  !i         :any other number: calculate both s and sdot
  !i   e     :energy of Hankel function.  e must be <=0
  !i   q     :Bloch wave number
  !i   p     :position of Hankel function center;
  !i         :structure constants are for expansion about the origin
  !i   nlma  :Generate coefficients S_R'L',RL for L' < nlma
  !i   nlmh  :Generate coefficients S_R'L',RL for L  < nlmh
  !i   ndim  :leading dimension of s,sdot
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   vol   :cell volume
  !i   awald :Ewald smoothing parameter
  !i   nkq   :number of direct-space lattice vectors
  !i   nkq   :number of reciprocal-space lattice vectors
  !i   dlv   :direct-space lattice vectors, units of alat
  !i   qlv   :reciprocal lattice vectors, units of 2pi/alat
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L quantum number for the C.G. coefficients (scg.f)
  !o Outputs
  !o   s      :structure constant matrix S_R'L',RL
  !o   sd     :Energy derivative of s
  ! ----------------------------------------------------------------
  !1.  Bloch phase.  For translation vectors T, it's sum_T exp(+i q T)

  !2.  Methfessel's definitions of Hankels and Bessel functions:

  !    h_0 = Re e^(ikr)/r and j = sin(kr)/kr, k = sqrt(e), Im k >=0.
  !    H_L = Y_L(-grad) h(r);   J_L = E^(-l) Y_L (-grad) j(r)

  !    They are related to the usual n_l and j_l by factors (I think)
  !       H_L  =  (i k)^(l+1) n_l (kr) Y_L   (E < 0)
  !       J_L  =  (i k)^(-l)  j_l (kr) Y_L   (E < 0)

  !   which explains how the energy-dependence is extracted out.
  !   Also cases e .ne. 0 and e .eq. 0 have the same equations.

  !r Expansion Theorem: H_{RL}(r) = H_L(r-R)
  !r   H_{RL}(E,r) = J_{R'L'}(E,r) * S_{R'L',RL}
  !r   S_R'L',RL = 4 pi Sum_l" C_{LL'L"} (-1)^l (-E)^(l+l'-l")/2 H_L"(E,R-R')
  ! ---
  implicit none
  integer :: mode,ndim,nlma,nlmh
  integer :: indxcg(*),jcg(*),nkd,nkq
  double precision :: p(3),q(3),alat,awald,vol,e,cg(*), dlv(*),qlv(*)
  double complex s(ndim,nlmh),sd(ndim,nlmh)
  integer :: lmxx,nrxmx,nlm0
  double precision :: fpi,p1(3),sp
  real(8),allocatable :: yl(:),efac(:)
  complex(8),allocatable :: dl(:),dlp(:)
  integer(4),allocatable :: sig(:)
  double complex phase,sumx,sud !dl(nlm0),dlp(nlm0)
  integer :: icg,icg1,icg2,ii,indx,ipow,l,lmax,nrx,nlm, ilm,ilma,la,ilmb,lh !sig(0:lmxx),
  logical :: ldot
  integer(4) :: job
  integer ::lmax_(1)
  real(8):: e_(1),rsm_(1),pp(3)
  ldot = .false.
  lmax = ll(nlma)+ll(nlmh)
  nlm = (lmax+1)**2
  nrx  = max(nkd,nkq)
  fpi  = 16d0*datan(1d0)
  if (nlma > ndim) call rxi('strxq: increase ndim: need',nlma)
  lmxx = lmax
  nlm0 =(lmxx+1)**2
  nrxmx= nrx
  allocate( yl(nrxmx*(lmxx+1)**2), efac(0:lmxx),sig(0:lmxx),dl(nlm0),dlp(nlm0))
  pp= matmul(transpose(qlat),p)
  call shortn3_plat(pp)
  p1 = matmul(plat,pp+nlatout(:,1))
  sp = fpi/2*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  job = 0
  if( sum(abs(p))<1d-10 ) job = 1
  lmax_(1)=lmax
  e_(1)=e
  rsm_(1)=0d0
  if (e < 0) then
     call hsmq(1,0,lmax_,e_,rsm_,job,q,p1,nrx,nlm0,yl,awald,alat,qlv,nkq,dlv,nkd,vol,dl,dlp)
  else
     call hsmqe0(lmax,0d0,job,q,p1,nrx,nlm0,yl, awald,alat,qlv,nkq,dlv,nkd,vol,dl)
     ldot = .false.
  endif
  if (sp /= 0d0) then
     do  20  ilm = 1, nlm
        dl(ilm) = phase*dl(ilm) ! ... Put in phase to undo shortening
        if (ldot) dlp(ilm) = phase*dl(ilm)
20   enddo
  endif
  ! --- Combine with Clebsch-Gordan coefficients ---
  ! ... efac(l)=(-e)**l; sig(l)=(-)**l
  efac(0) = 1
  sig(0) = 1
  do  l = 1, lmax
     efac(l) = -e*efac(l-1)
     sig(l) = -sig(l-1)
  enddo
  do  11  ilma = 1, nlma
     la = ll(ilma)
     do  14  ilmb = 1, nlmh
        lh = ll(ilmb)
        ii = max0(ilma,ilmb)
        indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1)-1
        sumx = 0d0
        sud = 0d0
        if (ldot) then
           do  16  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sumx = sumx + cg(icg)*efac(ipow)*dl(ilm)
              sud = sud + cg(icg)*efac(ipow)*(dlp(ilm)+ipow*dl(ilm)/e)
16         enddo
        else
           do  15  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sumx  = sumx + cg(icg)*efac(ipow)*dl(ilm)
15         enddo
        endif
        s(ilma,ilmb) = fpi*sig(lh)*dconjg(sumx)
        if (ldot) sd(ilma,ilmb) = fpi*dconjg(sud)*sig(lh)
14   enddo
11 enddo
  if (allocated(yl))deallocate(yl)
  if (allocated(efac))deallocate(efac)
  if (allocated(sig))deallocate(sig)
  if (allocated(dl))deallocate(dl)
  if (allocated(dlp))deallocate(dlp)
end subroutine strxq
subroutine strxq_all(e,q,nbas,bas,lx,lxx,alat,vol,awald,nkd,nkq,dlv,qlv,cg,indxcg,jcg, strx)
  ! strx(L1,ibas1,L2,ibas2) = 4pi*s of strxq for all pairs of atoms at once (e<0), on the device when strx is there
  ! (hvccfp0 keeps it there for vcoulq_4).  The Ewald sums of hsmq (rsm=0, only the value, not the e-derivative):
  !   Q-space part: sum_G Y_L(q+G) w(|q+G|) exp(i(q+G).p1) for all pairs p1 as one matrix product,
  !   real-space part: Y_L(p1-T) chi_l(|p1-T|) (hansr4) for each lattice vector T and pair, summed with exp(iqT),
  ! then the Clebsch-Gordan sums of strxq.  The pairs are p=0 (once, for the blocks ibas1=ibas2) and ibas1<ibas2;
  ! strx(L2,ibas2,L1,ibas1) = conj(strx(L1,ibas1,L2,ibas2)).  Y_L is ropyln's (r^l times real harmonics).
  use m_ll,only: ll
  use m_hamindex,only: plat,qlat
  use m_shortn3_plat,only: shortn3_plat,nlatout
#ifdef __GPU
  use m_blas, only: zmm => zmm_d, dmm => dmm_d, m_op_T
#else
  use m_blas, only: zmm => zmm_h, dmm => dmm_h, m_op_T
#endif
  implicit none
  integer,intent(in):: nbas,lxx,nkd,nkq,lx(nbas),indxcg(*),jcg(*)
  real(8),intent(in):: e,q(3),bas(3,nbas),alat,vol,awald,dlv(3,nkd),qlv(3,nkq),cg(*)
  complex(8):: strx((lxx+1)**2,nbas,(lxx+1)**2,nbas)
  integer,parameter:: lmxs=32
  real(8),parameter:: pi=4d0*datan(1d0), fpi=4d0*pi, y0=1d0/dsqrt(4d0*pi)
  integer:: lmax,nlm,nlxx,np,ip,ipc,ib,ib1,ib2,l,m,ilm,it,ig,istat,npc,ip0,ilma,ilmb,la,lh,ii,indx,icg,kk, &
       lncg,lnxcg,nla,nlb,ncgx
  integer,allocatable:: ipair(:,:),llx(:),nlat(:),jcgd(:),indxcgd(:)
  real(8):: tpiba,gam,a,a2,ah,rsm,akap,arsm,earsm,erfcarsm,pp(3),p(3),sp,cx0,f2m,pf,efac(0:lmxs),sgn(0:lmxs), &
       x,y,z,r2,r,ra,h0,wk,xx,xa,um,up,wk2,w,qq,q1,q2,cm,sm,cmx,chi(-1:lmxs),chi0m1,chi00
  real(8),allocatable:: p1(:,:),cx1(:,:),ca(:,:),cb(:,:),c0(:),c1(:),phr(:,:),mt(:,:,:),rr(:,:),cgd(:)
  complex(8),allocatable:: phs(:),aq(:,:),bq(:,:),dl(:,:)
  complex(8):: cof0,val,sumx
  complex(8),parameter:: img=(0d0,1d0)
  lmax = 2*lxx
  if(lmax>lmxs) call rx('strxq_all: 2*lxx > lmxs')
  if(e>=0d0)    call rx('strxq_all: e<0 only')
  nlm  = (lmax+1)**2
  nlxx = (lxx+1)**2
  ! --- pairs: p=0, then ibas1<ibas2 with p shortened as in strxq ---
  np = 1 + (nbas*(nbas-1))/2
  allocate(ipair(2,np), p1(3,np), phs(np))
  ipair(:,1) = 0; p1(:,1) = 0d0; phs(1) = 1d0
  ip = 1
  do ib1 = 1, nbas
    do ib2 = ib1+1, nbas
      ip = ip+1
      ipair(:,ip) = [ib1,ib2]
      p = bas(:,ib2)-bas(:,ib1)
      pp = matmul(transpose(qlat),p)
      call shortn3_plat(pp)
      p1(:,ip) = matmul(plat,pp+nlatout(:,1))
      sp = 2d0*pi*sum(q*(p-p1(:,ip)))
      phs(ip) = dcmplx(dcos(sp),dsin(sp))
    enddo
  enddo
  ! --- coefficients of the recursion for Y_L in ropyln ---
  allocate(cx1(0:lmax+1,0:lmax), ca(0:lmax,0:lmax), cb(0:lmax,0:lmax), c0(0:lmax), c1(0:lmax), llx(nlm), nlat(nbas))
  ca = 0d0; cb = 0d0
  f2m = 1d0   ! (2m)!
  pf  = 1d0   ! (2m-1)!!
  do m = 0, lmax
    if(m>0) f2m = f2m*(2*m-1)*(2*m)
    cx0 = dsqrt(1/fpi)
    if (m >0) cx0 = dsqrt((2*m+1)*2/fpi/f2m)
    cx1(m,m) = cx0
    do l = m, lmax
      cx1(l+1,m) = cx1(l,m)*dsqrt(dble((l+1-m)*(2*l+3))/dble((l+1+m)*(2*l+1)))
    enddo
    c0(m) = pf*cx1(m,m)
    pf = pf*(2*m+1)
    c1(m) = pf*cx1(m+1,m)
    do l = m+2, lmax
      ca(l,m) = -(l+m-1d0)/(l-m)*cx1(l,m)/cx1(l-2,m)
      cb(l,m) = (2*l-1d0)/(l-m)*cx1(l,m)/cx1(l-1,m)
    enddo
  enddo
  llx = [(ll(ilm),ilm=1,nlm)]
  nlat = (lx+1)**2
  efac(0) = 1d0; sgn(0) = 1d0
  do l = 1, lmax
    efac(l) = -e*efac(l-1)
    sgn(l) = -sgn(l-1)
  enddo
  ! --- constants of hsmq and hansr4 (a=awald, rsm=1/a) ---
  a = awald; a2 = a*a; gam = 0.25d0/a2; tpiba = 2d0*pi/alat
  rsm = 1d0/a; ah = 1d0/rsm
  akap = dsqrt(-e); arsm = akap*rsm/2; earsm = dexp(-arsm**2)/2
  erfcarsm = erfc(arsm)
  chi0m1 = -erfcarsm/akap                                        ! r=0: -h^s_-1 and -h^s_0 (hansr4)
  chi00  = akap*erfcarsm - 4d0*ah*earsm/dsqrt(4d0*datan(1d0))
  cof0 = fpi*dexp(gam*e)/vol
  ncgx = maxval(indxcg(1:(nlxx*(nlxx+1))/2+1))
  npc = max(1, min(np, int(2.5d8/(8d0*nkd*nlm))))                 ! pairs per chunk of the real-space table mt
  allocate(aq(nkq,nlm), bq(nkq,np), dl(nlm,np), phr(nkd,2), mt(nkd,nlm,npc), rr(nlm*npc,2))
  allocate(cgd(ncgx), jcgd(ncgx), indxcgd((nlxx*(nlxx+1))/2+1))
  cgd = cg(1:ncgx); jcgd = jcg(1:ncgx); indxcgd = indxcg(1:(nlxx*(nlxx+1))/2+1)
  phr(:,1) = [(dcos(2d0*pi*sum(q*dlv(:,it))), it=1,nkd)]
  phr(:,2) = [(dsin(2d0*pi*sum(q*dlv(:,it))), it=1,nkd)]
  !$acc data create(aq,bq,dl,mt,rr) copyin(ipair,p1,phs,cx1,ca,cb,c0,c1,llx,nlat,efac,sgn,phr,qlv,dlv,cgd,jcgd,indxcgd) present(strx)
  ! --- Q-space part: aq(G,L) = Y_L(q+G) (-exp(-gam |q+G|^2)/(e-|q+G|^2)) (-i)^l cof0 (pvhsmq), bq(G,p) = exp(i(q+G).p) ---
  !$acc parallel loop gang vector private(cm,sm,cmx,q1,q2,qq,w,x,y,z,r2,kk,l,m,ilm)
  do ig = 1, nkq
    x = tpiba*(q(1)+qlv(1,ig)); y = tpiba*(q(2)+qlv(2,ig)); z = tpiba*(q(3)+qlv(3,ig))
    r2 = x**2+y**2+z**2
    w = -dexp(-gam*r2)/(e-r2)
    cm = 1d0; sm = 0d0
    do m = 0, lmax
      if(m>0) then
        cmx = x*cm - y*sm
        sm  = y*cm + x*sm
        cm  = cmx
      endif
      q1 = 0d0; q2 = 0d0
      do l = m, lmax
        kk = l-m
        if(kk==0) then;     qq = c0(m)
        elseif(kk==1) then; qq = c1(m)*z
        else;               qq = ca(l,m)*r2*q2 + cb(l,m)*z*q1
        endif
        q2 = q1; q1 = qq
        aq(ig,l*(l+1)+1+m) = cm*qq*w*(-img)**l*cof0
        if(m/=0) aq(ig,l*(l+1)+1-m) = sm*qq*w*(-img)**l*cof0
      enddo
    enddo
  enddo
  !$acc parallel loop gang vector collapse(2)
  do ip = 1, np
    do ig = 1, nkq
      bq(ig,ip) = exp(img*alat*tpiba*sum((q(:)+qlv(:,ig))*p1(:,ip)))
    enddo
  enddo
  !$acc host_data use_device(aq,bq,dl)
  istat = zmm(aq, bq, dl, m=nlm, n=np, k=nkq, opA=m_op_T)
  !$acc end host_data
  ! --- real-space part, chunks of npc pairs: mt(T,L,p) = Y_L(p-T) chi_l(|p-T|), then sum_T mt exp(iqT) ---
  do ip0 = 1, np, npc
    npc = min(npc, np-ip0+1)
    !$acc parallel loop gang vector collapse(2) private(chi,cm,sm,cmx,q1,q2,qq,x,y,z,r2,r,ra,h0,wk,xx,xa,um,up,wk2,kk,l,m,ip)
    do ipc = 1, npc
      do it = 1, nkd
        ip = ip0+ipc-1
        x = alat*(p1(1,ip)-dlv(1,it)); y = alat*(p1(2,ip)-dlv(2,it)); z = alat*(p1(3,ip)-dlv(3,it))
        r2 = x**2+y**2+z**2
        if(r2 < 1d-12) then     ! the unsmoothed Hankel at T=p1 is left out (hansr4)
          chi(-1) = chi0m1
          chi(0) = chi00
          chi(1:lmax) = 0d0
        else
          r = dsqrt(r2); ra = r*ah
          h0 = dexp(-akap*r)/r
          wk = y0*dexp(-r2*a2)
          xx = earsm*wk/r
          xa = ra - arsm
          if(xa>0d0) then; um = h0-xx*erfcee_d(xa)
          else;            um = xx*erfcee_d(xa)
          endif
          up = xx*erfcee_d(ra + arsm)
          chi(-1) = (h0 - um - up)*r/akap
          chi(0) = h0 - um + up
          wk2 = 8*ah*earsm*wk
          do l = 1, lmax
            chi(l) = ((2*l-1)*chi(l-1) - e*chi(l-2) + wk2)/r2
            wk2 = 2d0*ah**2*wk2
          enddo
        endif
        cm = 1d0; sm = 0d0
        do m = 0, lmax
          if(m>0) then
            cmx = x*cm - y*sm
            sm  = y*cm + x*sm
            cm  = cmx
          endif
          q1 = 0d0; q2 = 0d0
          do l = m, lmax
            kk = l-m
            if(kk==0) then;     qq = c0(m)
            elseif(kk==1) then; qq = c1(m)*z
            else;               qq = ca(l,m)*r2*q2 + cb(l,m)*z*q1
            endif
            q2 = q1; q1 = qq
            mt(it,l*(l+1)+1+m,ipc) = cm*qq*chi(l)
            if(m/=0) mt(it,l*(l+1)+1-m,ipc) = sm*qq*chi(l)
          enddo
        enddo
      enddo
    enddo
    !$acc host_data use_device(mt,phr,rr)
    istat = dmm(mt, phr, rr, m=nlm*npc, n=2, k=nkd, opA=m_op_T)
    !$acc end host_data
    !$acc parallel loop gang vector collapse(2)
    do ipc = 1, npc
      do ilm = 1, nlm
        dl(ilm,ip0+ipc-1) = (dl(ilm,ip0+ipc-1) + dcmplx(rr(ilm+nlm*(ipc-1),1),rr(ilm+nlm*(ipc-1),2)))*phs(ip0+ipc-1)
      enddo
    enddo
  enddo
  ! --- Clebsch-Gordan sums of strxq for each pair; strx = 4pi*s ---
  !$acc kernels
  strx = 0d0
  !$acc end kernels
  !$acc parallel loop gang vector collapse(3) private(la,lh,ii,indx,sumx,val,ib1,ib2,nla,nlb)
  do ip = 1, np
    do ilmb = 1, nlxx
      do ilma = 1, nlxx
        la = llx(ilma); lh = llx(ilmb)
        ii = max(ilma,ilmb)
        indx = (ii*(ii-1))/2 + min(ilma,ilmb)
        sumx = 0d0
        do icg = indxcgd(indx), indxcgd(indx+1)-1
          sumx = sumx + cgd(icg)*efac((la+lh-llx(jcgd(icg)))/2)*dl(jcgd(icg),ip)
        enddo
        val = fpi*(fpi*sgn(lh)*dconjg(sumx))
        if(ip==1) then
          do ib = 1, nbas
            if(ilma<=nlat(ib) .and. ilmb<=nlat(ib)) strx(ilma,ib,ilmb,ib) = val
          enddo
        else
          ib1 = ipair(1,ip); ib2 = ipair(2,ip)
          if(ilma<=nlat(ib1) .and. ilmb<=nlat(ib2)) then
            strx(ilma,ib1,ilmb,ib2) = val
            strx(ilmb,ib2,ilma,ib1) = dconjg(val)
          endif
        endif
      enddo
    enddo
  enddo
  !$acc end data
  deallocate(ipair,p1,phs,cx1,ca,cb,c0,c1,llx,nlat,aq,bq,dl,phr,mt,rr,cgd,jcgd,indxcgd)
end subroutine strxq_all
pure real(8) function erfcee_d(ra) ! erfcee of util.f90 for the device: erfc(|x|)/y0/exp(-x*x)
  !$acc routine seq
  implicit none
  real(8),intent(in):: ra
  real(8):: w
  real(8),parameter:: &
       t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0, &
       t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0, &
       t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0, &
       t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0, &
       t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0, &
       t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0, &
       t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0, &
       t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0, &
       b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0, &
       b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0, &
       b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0, &
       b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0, &
       b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0, &
       b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0, &
       b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0, &
       b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0
  if (abs(ra) > 1.3d0) then   ! y0*dexp(-x*x)*f2(w=x-2) is erfc(x) for x>1.3
     w = abs(ra) - 2d0
     erfcee_d = (((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)*w+t21)*w+t20) &
          /  ((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)*w+b23)*w+b22)*w+b21)*w+1)
  else                        ! y0*dexp(-x*x)*f1(w=x-1/2) is erfc(x) for x<1.3
     w = abs(ra) - .5d0
     erfcee_d = (((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)*w+t11)*w+t10) &
          /  ((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)*w+b13)*w+b12)*w+b11)*w+1)
  endif
end function erfcee_d
endmodule m_hvccfp0_util
