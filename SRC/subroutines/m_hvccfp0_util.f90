module m_hvccfp0_util
  public mkradmatch,pmatorth,zgesvdnn2,mkb0,strxq
contains
subroutine mkradmatch( p, nxdim, rdmatch)
  !i  p(1,i): phi     at mt for i-th basis
  !i  p(2,i): dphi/dr at mt for i-th basis
  !o rdmatch(nxdim,nxdim)
  !-------
  !r    phinew_j(r) =sum_i phi_i(r)* rdmatch (i,j)
  !r     phinew_1(rmt)    =1      phinew_2(rmt)   =0
  !r   d phinew_1(rmt)/dr =0    d phinew_2(rmt)/dr=1
  !r for k >=3
  !r     phinew_k(rmt)    =0
  !r   d phinew_k(rmt)/dr =0
  !----------------------------------------------------
  implicit none
  integer(4):: nxdim,lbas,i,i1,i2,ix
  real(8):: p(1:2, 1:nxdim), rdmatch(1:nxdim,1:nxdim)
  real(8):: pd,p1,p1d,p2,p2d,s,t, eps=1d-3,delta
  !r                                       old     new
  !      write(6,"('mkradmatch: nxdim=',i4)") nxdim
  if(nxdim <=0) return
  ! top2rx 2013.08.09 kino      if(nxdim ==1) stop 'mkradmatch err nxdim==1'
  if(nxdim ==1) call rx( 'mkradmatch err nxdim==1')
  rdmatch=0d0
  !... pivot--- get better set of phi for augmentation
  do
     i1= nxdim
     i2= nxdim-1
     p1 = p(1, i1)
     p2 = p(1, i2)
     p1d= p(2, i1)
     p2d= p(2, i2)
     write(6,"('mkradmatch: i1 p1 p1d=',i3,2d13.6)") i1,p1,p1d
     write(6,"('mkradmatch: i2 p2 p2d=',i3,2d13.6)") i2,p2,p2d
     delta = p1*p2d-p2*p1d
     if(abs(delta) <eps*p1*p2) then
        if(i2==1) then
           write(6,"(' i1 i2=',2i5,2d13.6)") i1,i2,p1d/p1,p2d/p2
           ! top2rx 2013.08.09 kino            stop'mkradmatch: err poor linear dep'
           call rx( 'mkradmatch: err poor linear dep')
        endif
        i2=i2-1
     endif
     exit
  enddo
  !...
  call phimatch(1d0,0d0,  p1,p1d,p2,p2d, s,t)
  rdmatch(i1, 1)=  s
  rdmatch(i2, 1)=  t
  write(6,"('mkradmatch: 1 0    st=',2d13.5)") s,t
  call phimatch(0d0,1d0,  p1,p1d,p2,p2d, s,t)
  rdmatch(i1, 2)=  s
  rdmatch(i2, 2)=  t
  write(6,"('mkradmatch: 0 1    st=',2d13.5)") s,t

  ix=2
  do i= 1,nxdim
     if(i==i1 .OR. i==i2) cycle
     ix=ix+1
     !        write(6,"('mkradmatch: i p pd=',i3,2d13.5)") i,p(1,i),p(2,i)
     call phimatch(p(1,i),p(2,i),  p1,p1d,p2,p2d, s,t)
     rdmatch(i,  ix)=  1d0
     rdmatch(i1, ix)=  -s
     rdmatch(i2, ix)=  -t
     write(6,"('mkradmatch: ix st=',i3,2d13.5)") ix,s,t
  enddo
end subroutine mkradmatch

subroutine phimatch(p,pd, p1,p1d,p2,p2d, s,t)
  ! --- match for given p and pd
  !   phi = s phi1 + t phi2 !slope and value are at MT
  !     p  = s p1  + t p2
  !     pd = s pd1 + t pd2
  implicit none
  real(8):: matinv(2,2),p,pd,p1,p1d,p2,p2d,s,t,delta,ddd1,ddd2
  delta = p1*p2d-p2*p1d
  matinv(1,1) = 1/delta *  p2d
  matinv(1,2) = 1/delta * (-p2)
  matinv(2,1) = 1/delta * (-p1d)
  matinv(2,2) = 1/delta *  p1
  s = matinv(1,1) *p  + matinv(1,2) *pd
  t = matinv(2,1) *p  + matinv(2,2) *pd
  !... check
  ddd1 = abs(s*p1  + t*p2   -  p )
  ! top2rx 2013.08.09 kino      if(  ddd1 >1d-8 ) stop 'phimatch: ddd1 err'
  if(  ddd1 >1d-8 ) call rx( 'phimatch: ddd1 err')
  ddd2 = abs(s*p1d + t*p2d  -  pd)
  ! top2rx 2013.08.09 kino      if(  ddd2 >1d-8 ) stop 'phimatch: ddd2 err'
  if(  ddd2 >1d-8 ) call rx( 'phimatch: ddd2 err')
end subroutine phimatch

subroutine pmatorth(oo,oon,pmat,no,nn, pomat)
  ! get conversion matrix from old mixed basis(no) to augmented mixed basis(nn).
  ! pmatorth contains
  !   oo^{-1}_IJ
  implicit none
  integer(4):: no,nn,io,in,i
  complex(8):: pmat(no,nn),pomat(nn,no),oo(no,no),oon(nn,nn)
  complex(8),allocatable:: ooninv(:,:)
  real(8),allocatable:: eb(:)
  allocate(ooninv(nn,nn))
  ooninv = oon
  call matcinv(nn,ooninv) !generate ooninv
  !      pomat = matmul(ooninv, matmul(dconjg(transpose(pmat)),oo))
  pomat = transpose (matmul( oo, matmul(pmat,ooninv)))
  deallocate(ooninv)
end subroutine pmatorth
!      allocate(pp(nn,nn),ppin(nn,nn),eb(nn),zz(nn,nn),zze(nn,nn))
!      ppin = pp
!      call diagcvh(ppin,nn,eb,zz)
!      do i=1,nn
!        zze(:,i) =  zz(:,i)* sqrt(eb(i))
!      enddo
!      pomat = matmul(pmat, matmul(zze,dconjg(transpose(zz))))

! subroutine diagcvh(hh,ngb,eb,zz)
!   implicit none
!   integer(4):: nmx,nev,i,ngb
!   complex(8):: hh(ngb,ngb),oo(ngb,ngb),zz(ngb,ngb)
!   real(8):: eb(ngb)
!   nmx=ngb
!   oo = 0d0
!   do i=1,ngb
!      oo(i,i) = 1d0
!   enddo
!   call diagcv(oo,hh,zz,ngb, eb,nmx,1d99,nev)
!   write(6,*)' diagcvv: ngb,nev=',ngb,nev
!   do i=1,nev
!      write(6,'(i4,d23.16)')i, eb(i)
!   enddo
! end subroutine diagcvh
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine zgesvdnn2(no,nn, nnmx,epsmx, &
     pmat, &
     nnn)
  ! pmat(no,nn) ---> pmat(no,nnn)
  ! o input          pmat(no,nn)
  ! o output reduced pmat(no,nnn)
  implicit none
  integer(4):: lwork,info,nn,no,nnn,nnmx,i
  complex(8)::  pmat(no,nn),uu(no,no),vt(nn,nn)
  real(8):: ss(nn),epsmx
  real(8),allocatable:: rwork(:)
  complex(8),allocatable:: work(:),vtt(:,:),pmatx(:,:)
  !      write(6,*)' sumchk pmat=',sum(abs(pmat(1:no,1:nn)))
  lwork=4*no
  allocate(work(LWORK),rwork(5*no),pmatx(no,nn))
  pmatx =pmat
  call zgesvd('A','A',no,nn,pmat,no,SS,UU,no,VT,nn,work,lwork,rwork,info)
  nnn=-999
  do i=1,nn
     write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
     !         vtt(i,:)=ss(i)*vt(i,:)
     if(nnn==-999 .AND. ss(i)<epsmx) nnn = i-1
  enddo
  !      write(6,*) 'nnn=',nnn
  ! top2rx 2013.08.09 kino      if(nnn==0) stop 'strange: nnn=0'
  if(nnn==0) call rx( 'strange: nnn=0')
  if(nnn>nnmx) nnn=nnmx
  pmat=pmatx
  !      pmat(:,1:nnn) = uu(:,1:nnn)
  !      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
  !     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
  !      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk)))
  !     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
  !      deallocate(vtt)
end subroutine zgesvdnn2
subroutine mkb0( q, lxx,lx,nxx,nx, aa,bb, nrr,nrx,rprodx, &
     alat,bas,nbas,nbloch, &
     b0mat)
  !--make the matrix elementes < B_q | exp(iq r)>
  use m_ll,only: ll
  implicit none
  integer(4) :: nlx,l,n,m,nr,ir,lm,ibl1,ibas,nrx,nbloch
  integer(4) :: nbas,lxx, lx(nbas), nxx, nx(0:lxx,nbas),nrr(nbas)
  real(8)    :: rprodx(nrx,nxx,0:lxx,nbas),aa(nbas),bb(nbas), &
       phi(0:lxx),psi(0:lxx), bas(3,nbas), &
       alat, &
       pi,fpi,tpiba,qg1(3),q(3),absqg,r2s,a,b
  complex(8) :: b0mat(nbloch),img=(0d0,1d0) ,phase
  integer(4),allocatable:: ibasbl(:), nbl(:), lbl(:), lmbl(:)
  real(8),allocatable :: ajr(:,:),rofi(:),rob0(:,:,:)
  real(8),allocatable::cy(:),yl(:)
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
           ! ... bessel function
           ajr(ir,l) = phi(l)* rofi(ir) **(l +1 )
           ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
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
           call gintxx( ajr(1,l), rprodx(1,n,l,ibas), a,b,nr, &
                rob0(n,l,ibas) )
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
  deallocate(ajr, pjyl,rofi, &
       ibasbl, nbl, lbl, lmbl, &
       cy,yl,rob0)
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
!  if (allocated(wk))deallocate(wk)
  if (allocated(yl))deallocate(yl)
  if (allocated(efac))deallocate(efac)
  if (allocated(sig))deallocate(sig)
  if (allocated(dl))deallocate(dl)
  if (allocated(dlp))deallocate(dlp)
end subroutine strxq
endmodule m_hvccfp0_util
