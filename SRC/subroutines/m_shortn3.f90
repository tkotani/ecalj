!> subroutine pack related to shortest vectors. 
module m_shortvec ! core of m_shortn3: Find shortest vectors in modulo of rlat 
  public:: shortvec,shortvecinitialize
contains
  subroutine shortvec(pin,rlatp,xmx2,noutmx,nout,nlatout)
    !i pin is the fractional coodinate on rlat.
    !i rlatp,xmx2 are passed from shortvecinitialize
    !   rlatp(i,j)= sum( rlat(:,i)*rlat(:,j) )
    !   rlat(3,i) i-th vertor for modulo
    !i noutmax: upper limit of nlatout
    !o nout
    !o nlatout
    ! Shortest vectors are
    !    pin+nlatout(:,ix), where ix=1:nout, is the shortest vectors.
    !    We may have multiple nlatout (# is nout).
    implicit none
    integer:: noutmx
    integer:: nlatout(3,noutmx)
    integer:: nmax(3),nknknk,ik1,ik2,ik3,nout,nk,ik,i,j
    real(8):: rmax2,pin(3),eps=1d-8,rlat(3,3),xmax2(3),rr(3),rmin,nrmax(3)
    integer,allocatable:: nlat0(:,:)
    real(8),allocatable:: rnorm(:)
    real(8):: rlatp(3,3),xmx2(3)
    rmax2 = sum(pin*matmul(rlatp,pin)) + eps  ! eps is to make degeneracy safe.
    nrmax(:) =  sqrt(rmax2*xmx2(:))+abs(pin(:)) ! range of ix
    nmax =  nrmax
    ! we are looking for shortest vectors
    ik=0
    nknknk= (2*nmax(1)+1)*(2*nmax(2)+1)*(2*nmax(3)+1)
    allocate( nlat0(3, nknknk), rnorm(nknknk) )
    do ik1=-nmax(1),nmax(1)
       do ik2=-nmax(2),nmax(2)
          do ik3=-nmax(3),nmax(3)
             ik=ik+1
             nlat0(:,ik) = [ik1,ik2,ik3]
             rr= pin + nlat0(:,ik)
             rnorm(ik) = sum(rr*matmul(rlatp,rr))
          enddo
       enddo
    enddo
    nk=ik
    rmin=minval(rnorm(1:nk))
    nout=0
    do ik=1,nk
       rr= pin + nlat0(:,ik)
       if(rnorm(ik)<rmin+eps) then
          nout=nout+1
          if(nout>noutmx) call rx('shortn3: enlarge noutmx')
          nlatout(:,nout)=nlat0(:,ik)
       endif
    enddo
    deallocate(rnorm,nlat0)
    return
  end subroutine shortvec
  subroutine shortvecinitialize(rlat,rlatp,xmx2)    !!== Set translation vactors rlat(:,i),i=1,3 ==
    ! i rlat
    ! o rlatp,xmx2: these are passed to shortn3
    integer:: i,j
    real(8):: rlat(3,3)
    real(8):: rlatp(3,3),xmx2(3)
    do i=1,3
       do j=1,3
          rlatp(i,j) = sum(rlat(:,i)*rlat(:,j))
       enddo
    enddo
    call ellipsoidxmax(rlatp,xmx2)
  end subroutine shortvecinitialize
endmodule m_shortvec
module m_shortn3 
  ! Set rlat at first, then call shortn3(pin)
  ! Shortest p= pin + matmul(rlat(:,:),nlatout(:,i)) for i=1,nout
  ! Only when pin is on the Volonoi boundaries, nout>1.
  implicit none
  public:: gennlat,mshsiz,gvlst2,gvctof !shortn3_initialize, shortn3
  integer,parameter,private:: noutmx=48
  integer,private:: nout,nlatout(3,noutmx)
  real(8),private:: rlatp(3,3),xmx2(3)
contains
  subroutine gvctof(iopt,alat,plat,pos,n1,n2,n3,gmax,ng)! Makes k-space cutoff associated with mesh orders n1,n2,n3
    use m_ftox
    use m_lgunit,only:stdo
!    use m_shortn3,only: shortn3_initialize,shortn3,nout,nlatout
    !i Inputs
    !i   iopt       0 Use default (smaller of iopt=1,2)
    !i              1 use Nyquist cutoff
    !i              2 use cutoff for largest sphere in BZ
    !i   alat,plat  Real-space lattice vectors
    !i   n1,n2,n3   no. divisions along the three lattice vectors
    !o Outputs
    !o   gmax       Energy cutoff for these n1,n2,n3
    !o   ng         Number of lattice vectors
    !r Remarks
    !r   Adapted from nfp gvcutoff.f
    !u Updates
    !u   07 Feb 01 changed gmax tolerance to be consistent with gvlst2
    implicit none
    integer :: n1,n2,n3,ng,iopt
    double precision :: gmax,alat,plat(3,3),pos(3)
    ! ... Local parameters
    integer :: ipr,iprint,m,k1,k2,k3,ipv,nv1,nv2,ipr0,ipr1
    parameter (ipr0=30,ipr1=40)
    double precision :: pi,vol,tpiba,tol,volg,volgs,gvol,h1, &
         h2,h3,gg,h0,voln,gbot,g1,g2,g3,gmax1,vol1,qlat(3,3), &
         g(3),gs(3),plat1(3,3),qlat1(3,3),gmax2
    call pshpr(iprint()-20)
    ipr = iprint()
    pi = 4*datan(1d0)
    call dinv33(plat,1,qlat,vol)
    vol = dabs(alat**3*vol)
    tpiba = 2*pi/alat
    tol = 1d-8
    ! --- Basis vectors for real-space mesh and Q-space supercell ---
    plat1(:,1) = plat(:,1)/n1
    plat1(:,2) = plat(:,2)/n2
    plat1(:,3) = plat(:,3)/n3
    qlat1(:,1) = qlat(:,1)*n1
    qlat1(:,2) = qlat(:,2)*n2
    qlat1(:,3) = qlat(:,3)*n3
    ! --- Cutoff corresponding to recip supercell volume ---
    volg = (2*pi)**3 / vol
    volgs = n1*n2*n3* volg
    gvol = (volgs*3/(4*pi))**(1d0/3d0)
    ! --- Get shortest mesh pts where ki .ne. 0, gmax is pi/(longest) ---
    h1 = 1d10
    h2 = 1d10
    h3 = 1d10
    do   k1 = -5, 5
       do   k2 = -5, 5
          do   k3 = -5, 5
             g(:) = k1*plat1(:,1)+k2*plat1(:,2)+k3*plat1(:,3)
             gg = alat*dsqrt(g(1)**2+g(2)**2+g(3)**2)
             if (k1 /= 0) h1 = dmin1(h1,gg)
             if (k2 /= 0) h2 = dmin1(h2,gg)
             if (k3 /= 0) h3 = dmin1(h3,gg)
          enddo
       enddo
    enddo
    h0 = dmax1(h1,h2,h3)
    gmax = pi/h0
    if(ipr>=ipr0) write(stdo,ftox)' GVCTOF: mesh division',n1,n2,n3,'length',h1,h2,h3
    voln = (4*pi/3)*gmax**3
    ipv = int(100*voln/volgs)
    if (ipr >= ipr1) write(stdo,"('   Nyquist cutoff pi/h    ',f10.3,'    (volume', f10.2,',',i4,'%)')") gmax,voln,ipv
    gbot = 1d10
    do   k1 = -5, 5 !! --- Alternative: non-overlapping spheres on recip superlattice ---
       do   k2 = -5, 5
          do   k3 = -5, 5
             g1 = k1*qlat1(1,1)+k2*qlat1(1,2)+k3*qlat1(1,3) - pos(1)
             g2 = k1*qlat1(2,1)+k2*qlat1(2,2)+k3*qlat1(2,3) - pos(2)
             g3 = k1*qlat1(3,1)+k2*qlat1(3,2)+k3*qlat1(3,3) - pos(3)
             if (k1 /= 0 .OR. k2 /= 0 .OR. k3 /= 0) gbot = dmin1(gbot,g1*g1+g2*g2+g3*g3)
          enddo
       enddo
    enddo
    gmax1 = 0.5d0 * tpiba*dsqrt(gbot)
    vol1  = (4*pi/3)*gmax1**3
    ipv   = int(100*vol1/volgs)
    if(ipr >= ipr1) write(stdo,"('   largest sphere within BZ ',f8.3,'    (volume', f10.2,',',i4,'%)')") gmax1,vol1,ipv
    if(ipr >= ipr1) write(stdo,"('   cutoff for largest vector',f8.3,'    (volume', f10.2,',',i4,'%)')") gvol,volgs,100
    if(gmax1 < gmax .AND. iopt == 0) gmax = gmax1
    if(iopt == 2) gmax = gmax1 ! --- Count g-vectors within sphere of radius gmax ---
    gmax2 = (gmax-tol)**2 !
    nv1 = 0
    nv2 = 0
    s3: block
      real(8):: gg(3)
      call shortn3_initialize(qlat1) !initialization for shoten3
      do   k1 = 0, n1-1
         do   k2 = 0, n2-1
            do   k3 = 0, n3-1
               g = matmul(qlat,[k1,k2,k3]) - pos
               gg= matmul(g,plat1)
               call shortn3(gg)!,noutmx,nout,nlatout) !gg+nlatout are shortests of module qlat1.
               gs = matmul(qlat1,gg+nlatout(:,1))
               if ((tpiba*tpiba)*sum(gs**2) <= gmax2) nv2 = nv2+1
               nv1 = nv1+1
            enddo
         enddo
      enddo
    endblock s3
    if(ipr>=ipr0) write(stdo,ftox)' Reciprocal lattice: use sphere of radius',ftof(gmax),'keeping',nv2,'vectors of',nv1
    ng = nv2
    call poppr
  end subroutine gvctof
  subroutine mshsiz(alat,plat,gmax,ngabc,ng) !Finds dimensions for a mesh of G vectors that satisfy a cutoff
    use m_lgunit,only:stdo
    ! io  gmax  :On input, cutoff in reciprocal lattice vectors
    ! o        :Energy cutoff is gmax**2.
    ! o        :If input gmax is zero, it will be generated it from ngabc
    ! o        :(It is an error for both gmax and any of ngabc to be zero.)
    ! o        :On output, gmax is computed if input n1..n3 are nonzero,
    ! o        :or if job = 1
    ! io  ngabc :On input, max # divisions along the three recip lattice vecs.
    ! o        :(It is an error for both gmax and any of ngabc to be zero.)
    ! o        :Otherwise, those input ngabc initially zero are found that
    ! o        :satisfy cutoff criterion.
    ! o   ng    :number of G-vectors
    implicit none
    integer :: ngabc(3),ng,job,j
    double precision :: alat,plat(3,3),gmax
    logical :: change
    integer :: npfac,pfac(10),fmax,i,indx,iprint,i1,i2,i3,nn,nmin,nmx(3),nmxn(3),k,ngabcn(3),ngn,nginit,PRTG,kxx(1,1)
    double precision :: gmaxn,xx,q(3),gmax0,qlat(3,3),tpiba,facg,tolg,s1,s2,s3,ddot,gxx(1,1)
    character(256) :: outs
    parameter (fmax=600,tolg=1d0,PRTG=30)
    integer :: mshlst(0:3*fmax)
    logical:: fullmesh
    call tcn('mshsiz')
    q=0d0
    if (gmax == 0) then
       call gvctof(0,alat,plat,q,ngabc(1),ngabc(2),ngabc(3),gmax,ng)
       goto 999
    endif
    tpiba = 2*4d0*datan(1d0)/alat
    gmax0  = gmax/tpiba ! ... gmax0 = dimensionless gmax
    call dinv33(plat,1,qlat,xx)
    mshlst(0) = 0
    call ppfac(fmax,1,mshlst(1),mshlst(0))  !list of all allowed values of n1..n3
    if (mshlst(1) == 0) call rx('mshsiz: null set of allowed mesh points')
    ! --- Upper bound for n1..n3 : guaranteed to hold all G<gmax ---
    do  i = 1, 3
       i2 = mod(i,3)+1
       i3 = mod(i2,3)+1
       call gvlstn(qlat(1,i),qlat(1,i2),qlat(1,i3),q,mshlst,gmax0,nn)
       nmx(i)=nn
       if (ngabc(i) == 0) ngabc(i) = nn
    enddo
    gmaxn = gmax0
    facg = .995d0
    call pshpr(iprint()-30)
    call gvlst2(alat,plat,q,ngabc(1),ngabc(2),ngabc(3), 0d0,gmax,[0],000, 0,ng,kxx,gxx,kxx)
    ! ... Count the number of G vectors ng for initial n1..n3
    nginit = ng
    Getngloop: do !Reduce ng slowly until at least one of n1..n3 changes
       do  
          gmaxn = gmaxn * facg
          change = .false.
          do  i = 1, 3
             i2 = mod(i,3)+1
             i3 = mod(i2,3)+1
             call gvlstn(qlat(1,i),qlat(1,i2),qlat(1,i3),q,mshlst,gmaxn,nn)
             nmxn(i)=nn
             if (nmxn(i) /= nmx(i)) change = .TRUE. 
             !       The granularity of gvlstn may be too coarse.
             !       Don't assign, ngabcn(i) = nn but find next one smaller in mshlst
             indx = 1
             !call hunti(mshlst(1),mshlst,ngabc(i),0,indx)
             indx= findloc([(mshlst(j)<ngabc(i).and.ngabc(i)<=mshlst(j+1),j=1,mshlst(0)-1)],value=.true.,dim=1)
             indx = max(indx,1)
             ngabcn(i) = mshlst(indx)
          enddo
          if (change) exit
       enddo
       if(fullmesh()) then
          ngabc=ngabcn
          gmax=1d10
          ng=ngabc(1)*ngabc(2)*ngabc(3)
          exit
       endif   ! ... Count the number of G vectors for (smaller) trial ngabcn
       call gvlst2(alat,plat,q,ngabcn(1),ngabcn(2),ngabcn(3), 0d0,gmax,[0],000, 0,ngn,kxx,gxx,kxx)
       if (dble(ngn) >= nginit*tolg) then
          ng = ngn
          ngabc=ngabcn
       else
          exit
       endif
    enddo Getngloop
    call poppr
    i1 = ngabc(1)
    i2 = ngabc(2)
    i3 = ngabc(3)
    if (max0(i1,i2,i3) == fmax .AND. iprint() >= 10) then
       write(stdo,301)
       write(stdo,300) fmax
       write(stdo,301)
    endif
300 format(' WARNING!'/' At least one of the mesh divisions reached its maximal value fmax = ',i4/ &
         ' You might need to increase parameter fmax in mshsiz')
301 format(/1x,79('*')/)
    if (iprint() >= PRTG) then
       s1 = alat*sqrt(ddot(3,plat(1,1),1,plat(1,1),1))/i1
       s2 = alat*sqrt(ddot(3,plat(1,2),1,plat(1,2),1))/i2
       s3 = alat*sqrt(ddot(3,plat(1,3),1,plat(1,3),1))/i3
       write(stdo,"('MSHSIZ: mesh has ',i0,' x ',i0,' x ',i0,' divisions; length =',3f10.3)")i1,i2,i3,s1,s2,s3
       write(stdo,"('      generated from gmax (a.u.)=',f12.4,': ',i0,' vectors of ', &
            i0,' (',i0,'%)')")  gmax,ng,i1*i2*i3,(ng*100)/(i1*i2*i3)
    endif
999 continue
    call tcx('mshsiz')
  end subroutine mshsiz
  subroutine gvlst2(alat,plat,q,n1,n2,n3,gmin,gmax,mshlst,job,ngmx, ng,kv,gv,igv)!Setup a list of recip vectors within cutoff |q+G| < gmax
    use m_ftox
!    use m_shortn3,only: shortn3_initialize,shortn3,nout,nlatout
    use m_lgunit,only:stdo
    use m_mksym,only:  ngrp,gsym=>symops
    implicit none
    intent(in)::    alat,plat,q,n1,n2,n3,gmin,gmax,mshlst,job,ngmx
    intent(out)::                                                      kv,gv,igv 
    intent(inout)::                                                 ng
    !i   alat     Lattice constant
    !i   plat     Real-space primitive lattice vectors
    !i   q        vectors |q+G|<gmax are included in list
    !i   mshlst   if first entry is nonzero, a list of allowed values
    !i            n1,n2,n3 may take.  First entry is the size of
    !i            the list; then follows the list itself.
    !i   job      1s digit
    !i              0 return ng only
    !i              1 return kv and igv
    !i              8 return kv and gv, and then sorted
    !i                any combination of 1,8 is allowed
    !i          +500: gv without q 
    !i          +1000: symcheck
    !i   ngmx     Leading dimension of kv,gv,igv
    !i   gmin     Lower bound for reciprocal lattice vectors, in a.u.
    !i  n1..3    On input, max # divisions along the three lattice vectors.
    !i           (It is an error for both gmax and n1..n3 to be zero.)
    !i           Otherwise, input n1..n3 additionally constrain which
    !i           vectors are added to list; see Remarks.
    !i
    !i   gmax    On input, cutoff for reciprocal lattice vectors, in a.u.
    !i           Energy cutoff is gmax**2.
    !i           If input gmax is zero, gvlst2 will generate it from n1..n3
    !i           (It is an error for both gmax and n1..n3 to be zero.)
    !io   ng       Number of lattice vectors
    !o   kv       indices for gather/scatter operations.
    !o            kv(ig,i=1,2,3) for vector ig point to which entry
    !o            (i1,i2,i3) vector ig belongs
    !o   gv       list of reciprocal lattice vectors G
    !o   igv      list of reciprocal lattice vectors G, represented as
    !o            three integers (the multiples of qlat)
    !o            gv and igv are related by:
    !o              gv(1:3,1:ng) = 2*pi/alat * (qlat * igv(1:ng))
    !r Remarks
    !r   Collects a list of q + reciprocal lattice vectors (G+q) that lie
    !r   within a cutoff gmax.  List is optionally sorted by length (logical: lsort)
    !r   Vectors G are integer multiples of primitive r.l.v. qlat.
    !r
    !r   Additional constraints may be imposed, namely that the number of
    !r   multiples in each of axes 1..3 not exceed input values n1..n3.
    !r
    !r   For copying the list of points to/from a regular mesh, use:
    !r     call gvgetf(ng,1,kv,k1,k2,k3,c,c0)
    !r       to collect elements into list c0 from 3D array c
    !r     call gvputf(ng,1,kv,k1,k2,k3,c0,c)
    !r       to poke elements from list c0 into 3D array c
    !r     call gvaddf(ng,kv,k1,k2,k3,c0,c)
    !r       to add elements from list c0 into 3D array c
    !r
    !u Updates See git log after 2009
    !u   11 Jul 08 New argument gmin
    !u   01 Jun 01 revised gvlst2 together with gvlist.  They form similar
    !u             operations but with different functions; see gvlist.f
    !u   26 Mar 01 Another bug fix for input n1..n3 and gmax nonzero
    !u   06 Mar 01 Bug fix for input n1..n3 and gmax nonzero
    !u   Routine was adapted from T. Kotani, routine getgv2
    ! ----------------------------------------------------------------------
    integer :: n1,n2,n3,ng,ngmx,kv(ngmx,3),igv(ngmx,3),mshlst(0:*) !,igv2(3,*)
    double precision :: alat,gmin,gmax,gv(ngmx,3),plat(3,3),q(3)
    integer :: ig,n1max,n1min,n2max,n2min,n3max,n3min,i1,i2,i3,nn,i
    integer :: n1l,n2l,n3l
    integer :: PRTG,PRTG2,iset(3),ipr,job,job0,job1,job8,k1,k2,k3
    double precision :: qlat(3,3),vol,pi,tpiba,qpg(3),q2
    double precision :: gmin0,gmax0,gmin2,gmax2,h1,h2,h3,ddot,tol
    character(256) :: outs
    parameter (PRTG=30,PRTG2=100,tol=1d-8)
    real(8),parameter:: tolg2=1d-2
    logical::  lgpq,lgv,lsort,ligv !,ligv2
    real(8):: plat1(3,3),qlat1(3,3),gg,gs(3)
    integer:: j1,j2,j3,m,jj1,jj2,jj3,nn1,nn2,nn3,i123(3),jjj(3),jg,igrp,jjg
    real(8):: rlatp(3,3),xmx2(3),gvv(3),diffmin
    integer :: nginit,kv_tmp(ngmx,3),igv_tmp(ngmx,3),ips(ngmx),jx,nxx,itemp(ngmx),ix,iprint
    real(8):: gv_tmp(ngmx,3)
    call tcn('gvlst2')
    call pshpr(iprint()-30)
    call getpr(ipr)
    call dinv33(plat,1,qlat,vol)
    pi = 4d0*datan(1d0)
    tpiba = 2*pi/alat
    gmin0  = gmin/tpiba
    gmax0  = gmax/tpiba
    if (gmin < 0)  call rx('gvlst2: input gmin <= 0')
    if (gmax <= 0) call rx('gvlst2: input gmax <= 0')
    job0 = mod(job,100)
    ligv  = mod(job0,2)/=0    !1,9 !Get igv
    lgv   = mod(job0/4,4)/=0  !8,9 !Get gv
    lsort = mod(job0/8,2)/=0  !8,9 sorted
    lgpq  = mod(job/100,10)>4 !+500 or not

    !! ... Basis vectors for real-space mesh and recip-space supercell
    nn1=n1
    nn2=n2
    nn3=n3
    n1min=0; n1max=nn1-1
    n2min=0; n2max=nn2-1
    n3min=0; n3max=nn3-1
    if(nn1*nn2*nn3==0) then
       print *,'nn1,nn2,nn3=',nn1,nn2,nn3
       call rx('gvlst2: require n1*n2*n3/=0')
    endif
    !  if(nn1==0) call gvlstn(qlat(1,1),qlat(1,2),qlat(1,3),q,mshlst,gmax0,nn1) 
    !  if(nn2==0) call gvlstn(qlat(1,2),qlat(1,3),qlat(1,1),q,mshlst,gmax0,nn2) 
    !  if(nn3==0) call gvlstn(qlat(1,3),qlat(1,1),qlat(1,2),q,mshlst,gmax0,nn3) 
    do m = 1, 3
       plat1(m,:) = plat(m,:)*[1d0/nn1,1d0/nn2,1d0/nn3]
       qlat1(m,:) = qlat(m,:)*[nn1,nn2,nn3]
    enddo
    n1l=nn1
    n2l=nn2
    n3l=nn3
    gmax2 = (gmax0-tol)**2 !this tol is needed for gmax given by gvctof
    gmin2 = gmin0**2
    call shortn3_initialize(qlat1) !initialization for m_shoten3 for qlat1
    ig=0
    j1loop: do j1 = 0,nn1-1 !! --- g vectors, shorten, count and keep if within gmax ---
       j2loop: do j2 = 0,nn2-1 
          j3loop: do j3 = 0,nn3-1
             jjj=[j1,j2,j3] 
             qpg= [1d0/nn1, 1d0/nn2, 1d0/nn3] * (jjj+ matmul(q,plat(:,:))) !frac corrdinate on qlat1.
             call shortn3(qpg) ! return nout,nlatout
             gs = matmul(qlat1(:,:), (qpg+nlatout(:,1)))
             gg = sum(gs**2)
             if(gmin2<= gg .AND. gg < gmax2) then
                ig = ig+1
                if (job0==0) cycle
                if (ig > ngmx) write(stdo,*) 'ERROR: ig ngmx=',ig,ngmx
                if (ig > ngmx) call rx('gvlist2: ng exceeds ngmx')
                kv_tmp(ig,:) = jjj+1
                if(ligv) igv_tmp(ig,:) = jjj + [nn1,nn2,nn3]*nlatout(:,1) !integer index of Gvec
                if(lgv)  gv_tmp(ig,:)= merge(gs,gs-q,mask=lgpq)
             endif
          enddo j3loop
       enddo j2loop
    enddo j1loop
    ng = ig
    if(job0/=0) then !symmetry checker
       nginit = ng
       ng=0
       ips=0
       do ig = 1,nginit
          if(job>999.and.ips(ig)==0) then 
             itemp=0
             ix=0
             do igrp = 1, ngrp
                gvv = matmul(gsym(:,:,igrp),gv_tmp(ig,:)) !  ... gvv = g(k) gv
                !write(6,ftox) 'iii ig igrp gvv=',ig,igrp,ftof(gvv)
                jg = findloc([(sum(abs(gvv-gv_tmp(jg,:)))<tolg2,jg=1,nginit)],value=.true.,dim=1)
                if(jg/=0) then
                   ix=ix+1
                   itemp(ix)=jg
                else
                   goto 70
                endif
             enddo
             ips(itemp(1:ix))=1 !jg=itemp(1:ix) need to be included.
          endif
          ng = ng+1
          kv(ng,:) = kv_tmp(ig,:)
          if(ligv) igv(ng,:)=igv_tmp(ig,:)
          if(lgv )  gv(ng,:)= gv_tmp(ig,:) 
70        continue
       enddo !   write(stdo,ftox)'gmax ng nginit ngmx=',ftof(gmax),ng,nginit,ngmx,ftof(q)
    endif
    if(lsort) then
       gvsort:block !Sort vectors -- !call gvlsts(ng,gv(1:ng,1:3),kv(1:ng,1:3),igv(1:ng,1:3),ligv) 
         integer:: iprm(ng)
         call dvshel(1,ng, sum(gv(1:ng,1:3)**2,dim=2)*[((1d0 + 1d-15*ig),ig=1,ng)], iprm,1)
         gv(1:ng,1:3) = gv(iprm+1,1:3)
         kv(1:ng,1:3) = kv(iprm+1,1:3)
         if(ligv) igv(1:ng,1:3) =igv(iprm+1,1:3)
       endblock gvsort !write(stdo,ftox)'gmax gv=',ftof(gv(1,1:3))
    endif
    if(ipr>=PRTG)write(stdo,ftox)'gvlst2: gmax=',ftof(gmax,3),'a.u. created',ng,'vectors of',n1l*n2l*n3l,&
         '(',(ng*100)/(n1l*n2l*n3l),'%)'
    if(ipr >= PRTG2 .AND. ng > 0 .AND. ligv) then
       write(stdo,"(' G vectors (multiples of reciprocal lattice vectors)'/ '   ig    G1   G2   G3     E')")
       do  ig = 1, ng
          i123 = igv(ig,:)
          qpg= q + matmul(qlat,i123) 
          write(stdo,"(i5,1x,3i5,2x,f8.4)")  ig,i123,sum(qpg**2) *tpiba**2
          write(stdo,"('q  qpg=',3d13.5,3x,3d13.5)")  q,qpg
       enddo
    endif
    call poppr
    call tcx('gvlst2')
  end subroutine gvlst2
  subroutine gennlat(pos,nbas,plat,nk1,nk2,nk3,npairmx,npair,nlat,nqwgt,nlatS,nlatE,ok)
!    use m_shortn3,only: shortn3_initialize,shortn3,nout,nlatout
    implicit none
    intent(in)::       pos,nbas,plat,nk1,nk2,nk3,npairmx
    intent(out)::                                        ok!,npair,  t,nqwgt
    integer:: nbas,nk1,nk2,nk3,npairmx ,npair(nbas,nbas)
    integer:: ib1,ib2,nnn,ik,nmax(3),ix,i,j,ik1,ik2,ik3,ni
    integer:: nd(3),nlat(3,npairmx,nbas,nbas)
    integer:: debug=0,nqwgt(npairmx,nbas,nbas)
    real(8):: pos(3,nbas),plat(3,3), pi,qlat(3,3),q(3),pin(3),rmax2 !,nrmax(3)
    real(8):: eps=1d-8
    real(8):: posp(3),rxlat(3,3),rxprod(3,3),dummy,rrr,rmax
    logical :: ok
    real(8):: rlatp(3,3),xmx2(3)
    integer:: nlatS(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas)
    integer:: nlatE(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas)
    call tcn('gennlat')
    !      print *,'gennlat:'
    ok=.true.
    pi=4d0*datan(1d0)
    call dinv33(plat,1,qlat,dummy)
    ! periodicity for R mesh for nk1 nk2 nk3 division.
    rxlat(:,1)=plat(:,1)*nk1
    rxlat(:,2)=plat(:,2)*nk2
    rxlat(:,3)=plat(:,3)*nk3
    call shortn3_initialize(rxlat)
    !      do i=1,3
    !      do j=1,3
    !        rxprod(i,j) = sum(rxlat(:,i)*rxlat(:,j))
    !      enddo
    !      enddo
    !      print *, 'goto ellips',rxprod
    !      call ellipsoidxmax(rxprod,xmx2)
    !      print *, 'end of ellips xmx2=',xmx2
    ! et Maxmum x_i**2 for ellipsoid for 1d0 = sum_{i,j} x_i rxprod(i,j) x_j
    ! Maximum value for x_i for ellipsoid
    ! Ellipsoid is given as 1d0 = sum x_i nn(i,j) x_j.

    !      allocate(nlatS(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas))
    !      allocate(nlatE(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas))
    npair=0
    do ib1=1,nbas
       do ib2=1,nbas
          if(debug>0) print *
          if(debug>0) print *,' ---- ib1 ib2=',ib1,ib2
          posp = matmul( pos(:,ib1)-pos(:,ib2), qlat(:,:)) !posp is on plat-coodinate
          nnn=0
          do ik1=0,nk1-1
             do ik2=0,nk2-1
                do ik3=0,nk3-1
                   ! rint *,' ik1,ik2,ik3=',ik1,ik2,ik3
                   nd=[ik1,ik2,ik3]
                   if( ik1 >nk1/2) nd= nd-[nk1, 0,  0]
                   if( ik2 >nk2/2) nd= nd-[0, nk2,  0]
                   if( ik3 >nk3/2) nd= nd-[0, 0,  nk3]
                   ! above procedures are not necessary, but these may accelarate shortn3 a little.
                   pin(1)= (posp(1)+nd(1))/nk1
                   pin(2)= (posp(2)+nd(2))/nk2
                   pin(3)= (posp(3)+nd(3))/nk3
                   call shortn3(pin) ! return nout,nlatout
                   if(nnn+nout>npairmx) then
                      ok=.false.
                      return
                   endif
                   nlatS(ik1,ik2,ik3,ib1,ib2)= nnn+1
                   nlatE(ik1,ik2,ik3,ib1,ib2)= nnn+nout
                   do ix=1,nout
                      nlat(:,nnn+ix,ib1,ib2)= &
                           nd + [nk1*nlatout(1,ix), nk2*nlatout(2,ix), nk3*nlatout(3,ix)]
                      nqwgt(nnn+ix,ib1,ib2)= nout
                   enddo
                   do ix=1,nout
                      nqwgt(nnn+ix,ib1,ib2)= nout
                   enddo
                   nnn = nnn+nout
                enddo
             enddo
          enddo
          npair(ib1,ib2)= nnn
          !$$$            write(6,"(a,2i4,2x,i6)") ' ib1 ib2 npair=',ib1,ib2,npair(ib1,ib2)
          !$$$            do ni = 1,npair(ib1,ib2)
          !$$$              posp =  pos(:,ib1)-pos(:,ib2) + matmul(plat,nlat(:,ni,ib1,ib2))
          !$$$              rrr = sqrt(sum(posp**2))
          !$$$              write(6,"(i6,3x,3i3,f8.3,i5)") ni,nlat(1:3,ni,ib1,ib2),rrr,nqwgt(ni,ib1,ib2)
          !$$$            enddo
       enddo
    enddo
    do ib1=1,nbas
       do ib2=1,nbas
          if(abs(sum(1d0/nqwgt(1:npair(ib1,ib2),ib1,ib2))-nk1*nk2*nk3)>1d-8) call rx('bug:nqwgt sum is not unity')
       enddo
    enddo
    call tcx('gennlat')
  end subroutine gennlat
  subroutine gvlstn(q0,q1,q2,qp,mshlst,gmax0,nn)! Multiples of r.l.v. that bound cutoff gmax0  (r.l.v. reciplocal lattice vector)
    implicit none
    intent(in) ::   q0,q1,q2,qp,mshlst,gmax0
    intent(out)::                            nn
    !i   q0    :r.l.v. for which bounds nmin and nmax are to be computed
    !i   q1    :first  r.l.v. different from q0
    !i   q2    :second r.l.v. different from q0
    !i   qp    :k-point added to G vectors: sphere is centered G=qp.
    !i   mshlst:An ordered list of integers used to restrict the assignment
    !i         :of nn to one of an allowed list of points, should a value
    !i         :of nn be assigned (see description for nn below).
    !i         :If first entry is nonzero, mshlst(1..) = list of allowed
    !i         :values nn may take.  mshlst(0) is the size of mshlst.
    !i         :Certain fourier transform routines have restrictions
    !i         :on the allowed mesh sizes; this constraint is designed
    !i         :to handle that restriction.
    !i  gmax0  :cutoff G
    !o Outputs
    !o     nn :   meshsize
    !r Remarks
    !r   q0,q1,q2,qp and G are all dimensionless (units of 2*pi/a)
    ! ----------------------------------------------------------------------
    integer :: nmin,nmax,nn,mshlst(0:*), indx,i
    double precision :: q0(3),q1(3),q2(3),qp(3),gmax0,qperp(3),ddot,qqperp
    qperp=[q1(2)*q2(3)-q1(3)*q2(2),  q1(3)*q2(1)-q1(1)*q2(3), q1(1)*q2(2)-q1(2)*q2(1)] 
    qperp  = qperp/sum(qperp**2)**.5 ! = q1 x q2 / |q1 x q2| 
    qqperp = sum(q0*qperp)           
    nmax =  gmax0/abs(Qqperp) - sum(qp*qperp)/Qqperp + 1
    nmin = -gmax0/abs(Qqperp) - sum(qp*qperp)/Qqperp - 1
    nn = 2*max(iabs(nmin),nmax)+1
    if(mshlst(0) /= 0) then
       indx = 1
       indx= findloc([(mshlst(i)<nn.and.nn<=mshlst(i+1),i=1,mshlst(0)-1)],value=.true.,dim=1)
       nn = mshlst(min(indx+1,mshlst(0)))
    endif
  end subroutine gvlstn
  subroutine shortn3_initialize(rlat)
    use m_shortvec,only: shortvecinitialize
    real(8):: rlat(3,3)
    call shortvecinitialize(rlat,rlatp,xmx2)
  end subroutine shortn3_initialize
  subroutine shortn3(pin) !pin is in fractional coordinate in rlat(:,i)
    use m_shortvec,only: shortvec
    real(8):: pin(3)
    call shortvec(pin,rlatp,xmx2,noutmx,nout,nlatout)
  end subroutine shortn3
  subroutine ppfac(fmax,job,fac,nfac)   !- Find all products of prime factors within some maximum
    !i   npfac :number of prime factors
    !i   pfac  :list of prime factors
    !i   fmax  :maximum value of product
    !i   job   :0 list is returned unsorted
    !i         :1 list is returned sorted
    !o Outputs
    !o   fac   :all products of pfac <= fmax
    !o         :List is returned sorted.
    !o         :fac must be dimensioned at least nfac
    !o         :or 3*nfac if fac is returned sorted.
    !o   nfac  :length of fac
    ! ----------------------------------------------------------------------
    implicit none
    integer,parameter:: npfac=5,pfac(1:npfac)=[2,3,5,7,11]
    integer :: fac(*),fmax,job,nfac
    integer :: i,m,k,fack,nfaci
    nfac = 0
    nfaci = 0
    do  i = 1, npfac
       nfac = nfac+1
       nfaci = nfac
       fac(nfac) = pfac(i)
       do  k  = 1, nfaci
          fack = fac(k)
          do  m = 1, fmax
             fack = fack*pfac(i)
             if (fack <= fmax) then
                nfac = nfac+1
                fac(nfac) = fack
             else
                goto 10
             endif
          enddo
10        continue
       enddo
    enddo
    if (job /= 0) call ivheap(1,nfac,fac,fac(1+2*nfac),0)
  end subroutine ppfac
endmodule m_shortn3
subroutine ellipsoidxmax(nn, xmx2)!== Maximum value for x_i for ellipsoid ==
  implicit none
  !!  Ellipsoid is given as 1d0 = sum x_i nn(i,j) x_j.
  !i nn(3,3)
  !o xmx2(i)  Maximum of x_i**2
  real(8):: nn(3,3),v2(2),ainv(2,2), rmax2, xmx2(3),det,fac1,fac2,fac3,nv2(2)
  associate(&
       n11=>nn(1,1),n12=>nn(1,2),n13=>nn(1,3), &
       n21=>nn(2,1),n22=>nn(2,2),n23=>nn(2,3), &
       n31=>nn(3,1),n32=>nn(3,2),n33=>nn(3,3) )
    ainv(1,1:2)=  [ n33,-n23]
    ainv(2,1:2)=  [-n32, n22]
    fac1 = n11-sum([n12,n13]*matmul(ainv,[n12,n13]))/(n22*n33-n23*n32)
    ainv(1,1:2)=  [ n11,-n31]
    ainv(2,1:2)=  [-n13, n33]
    fac2 = n22-sum([n23,n21]*matmul(ainv,[n23,n21]))/(n33*n11-n31*n13)
    ainv(1,1:2)=  [ n22,-n12]
    ainv(2,1:2)=  [-n21, n11]
    fac3 = n33-sum([n31,n32] *matmul(ainv,[n31,n32]))/(n11*n22-n12*n21)
    xmx2 = [1d0/fac1,1d0/fac2,1d0/fac3]
  endassociate
end subroutine ellipsoidxmax
