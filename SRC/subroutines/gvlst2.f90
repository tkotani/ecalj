!module m_gvlst2
!  public gvlst2,mshsiz
!  private
!  contains
subroutine gvlst2(alat,plat,q,n1,n2,n3,gmin,gmax,mshlst,job,ngmx, ng,kv,gv,igv,igv2)
  use m_ftox
  use m_shortn3,only: shortn3_initialize,shortn3,nout,nlatout
  use m_lgunit,only:stdo
  implicit none
  intent(in)::    alat,plat,q,n1,n2,n3,gmin,gmax,mshlst,job,ngmx  
  intent(out)::                                                   ng,kv,gv,igv,igv2 
  !- Set up a list of recip vectors within cutoff |q+G| < gmax
  ! ----------------------------------------------------------------------
  !i   alat     Lattice constant
  !i   plat     Real-space primitive lattice vectors
  !i   q        vectors |q+G|<gmax are included in list
  !i   mshlst   if first entry is nonzero, a list of allowed values
  !i            n1,n2,n3 may take.  First entry is the size of
  !i            the list; then follows the list itself.
  !i   job      1s digit
  !i              0 return ng only
  !i              1 return kv and igv
  !i              2 return kv and igv2
  !i              4 return kv and gv
  !i              8 return kv and gv, and sort list
  !i                any combination of 1,2,4 is allowed
  !i          +500: gv without q (sorted when lsort=T)
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
  !o   ng       Number of lattice vectors
  !o   kv       indices for gather/scatter operations.
  !o            kv(ig,i=1,2,3) for vector ig point to which entry
  !o            (i1,i2,i3) vector ig belongs
  !o   gv       list of reciprocal lattice vectors G
  !o   igv      list of reciprocal lattice vectors G, represented as
  !o            three integers (the multiples of qlat)
  !o            gv and igv are related by:
  !o              gv(1:3,1:ng) = 2*pi/alat * (qlat * igv(1:ng))
  !o   igv2     same as igv except first and second columns are permuted
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
  integer :: n1,n2,n3,job,ng,ngmx,kv(ngmx,3),igv(ngmx,3),igv2(3,*),mshlst(0:*)
  double precision :: alat,gmin,gmax,gv(ngmx,3),plat(3,3),q(3)
  integer :: ig,n1max,n1min,n2max,n2min,n3max,n3min,i1,i2,i3,nn,i
  integer :: n1l,n2l,n3l
  integer :: PRTG,PRTG2,iset(3),ipr,job0,job1,job2,job8,k1,k2,k3
  double precision :: qlat(3,3),vol,pi,tpiba,qpg(3),q2
  double precision :: gmin0,gmax0,gmin2,gmax2,h1,h2,h3,ddot,tol
  character(256) :: outs
  parameter (PRTG=30,PRTG2=100,tol=1d-8)
  logical::  lgpq,lgv,lsort,ligv,ligv2
  real(8):: plat1(3,3),qlat1(3,3),gg,gs(3)
  integer:: j1,j2,j3,m,jj1,jj2,jj3,nn1,nn2,nn3
  real(8):: rlatp(3,3),xmx2(3)
  call getpr(ipr)
  call dinv33(plat,1,qlat,vol)
  pi = 4d0*datan(1d0)
  tpiba = 2*pi/alat
  gmin0  = gmin/tpiba
  gmax0  = gmax/tpiba
  if (gmin < 0)  call rx('gvlst2: input gmin <= 0')
  if (gmax <= 0) call rx('gvlst2: input gmax <= 0')
  job0 = mod(job,100)
  job1 = mod(job0,2)
  job2 = mod(job0/2,2)
  ligv  = mod(job0,2)/=0
  ligv2 = mod(job0/2,2)/=0
  lgv   = mod(job0/4,4)/=0
  lsort = mod(job0/8,2)/=0
  lgpq  = mod(job/100,10)>4
  !! ... Basis vectors for real-space mesh and recip-space supercell
  nn1=n1
  nn2=n2
  nn3=n3
  n1min=0; n1max=nn1-1
  n2min=0; n2max=nn2-1
  n3min=0; n3max=nn3-1
  if(nn1==0) call gvlstn(qlat(1,1),qlat(1,2),qlat(1,3),q,mshlst,gmax0,nn1) 
  if(nn2==0) call gvlstn(qlat(1,2),qlat(1,3),qlat(1,1),q,mshlst,gmax0,nn2) 
  if(nn3==0) call gvlstn(qlat(1,3),qlat(1,1),qlat(1,2),q,mshlst,gmax0,nn3) 
  do m = 1, 3
     plat1(m,1) = plat(m,1)/nn1
     plat1(m,2) = plat(m,2)/nn2
     plat1(m,3) = plat(m,3)/nn3
     qlat1(m,1) = qlat(m,1)*nn1
     qlat1(m,2) = qlat(m,2)*nn2
     qlat1(m,3) = qlat(m,3)*nn3
  enddo
  n1l=nn1
  n2l=nn2
  n3l=nn3
  !! --- Loop through g vectors, shorten, count and keep if within gmax ---
  ig = 0
  gmax2 = (gmax0-tol)**2
  gmin2 = gmin0**2
  call shortn3_initialize(qlat1) !initialization for m_shoten3
  do  212  j1 = 0,nn1-1 
     do  211  j2 = 0,nn2-1 
        do  21  j3 = 0,nn3-1 
           qpg= [j1, j2, j3] + matmul(q,plat(:,:))
           qpg= [qpg(1)/dble(nn1), qpg(2)/dble(nn2), qpg(3)/dble(nn3)]
           ! qpg is in the fractional corrdinate on qlat1.
           call shortn3(qpg) ! return nout,nlatout
           gs= matmul(qlat1(:,:), (qpg+nlatout(:,1)))
           gg = (gs(1)**2+gs(2)**2+gs(3)**2)
           k1 = j1+1 
           k2 = j2+1 
           k3 = j3+1 
           if(gmin2<= gg .AND. gg < gmax2) then
              ig = ig+1
              if (job0 /= 0) then
                 if (ig > ngmx) then
                    print *,' ig ngmx=',ig,ngmx
                    call rx('gvlist: ng exceeds ngmx')
                 endif
                 kv(ig,1) = k1!j1+1
                 kv(ig,2) = k2!j2+1
                 kv(ig,3) = k3!j3+1
                 if (ligv .OR. ligv2) then
                    jj1= j1 + nn1*nlatout(1,1)
                    jj2= j2 + nn2*nlatout(2,1)
                    jj3= j3 + nn3*nlatout(3,1)
                    if (ligv) then
                       igv(ig,1) = jj1
                       igv(ig,2) = jj2
                       igv(ig,3) = jj3
                    endif
                    if (ligv2) then
                       igv2(1,ig) = jj1
                       igv2(2,ig) = jj2
                       igv2(3,ig) = jj3
                    endif
                 endif
                 if (lgv) then
                    if (lgpq) then
                       gv(ig,1) = gs(1)
                       gv(ig,2) = gs(2)
                       gv(ig,3) = gs(3)
                    else
                       gv(ig,1) = gs(1) - q(1)
                       gv(ig,2) = gs(2) - q(2)
                       gv(ig,3) = gs(3) - q(3)
                    endif
                 endif
              endif
           endif
21      enddo
211  enddo
212 enddo
  ng = ig
  if (ipr >= PRTG .AND. n1l*n2l*n3l == 0) then
     write(stdo,ftox)' GVLST2: gmax=',ftof(gmax0*tpiba,3),'created',ng,' recip. lattice vectors'
  elseif (ipr >= PRTG) then
     h1 = alat*sqrt(ddot(3,plat(1,1),1,plat(1,1),1))/n1l
     h2 = alat*sqrt(ddot(3,plat(1,2),1,plat(1,2),1))/n2l
     h3 = alat*sqrt(ddot(3,plat(1,3),1,plat(1,3),1))/n3l
     write(stdo,ftox)'gvlst2: gmax=',ftof(gmax,3),'a.u. created',ng,'vectors of',n1l*n2l*n3l,&
          '(',(ng*100)/(n1l*n2l*n3l),'%)'
  endif
  if(lsort) call gvlsts(ngmx,ng,gv,kv,igv,igv2,job1,job2)! --- Sort the list of vectors --
  if (ipr >= PRTG2 .AND. ng > 0 .AND. job1+job2 /= 0) then
     write(stdo,333)
333  format(' G vectors (multiples of reciprocal lattice vectors)'/ '   ig    G1   G2   G3     E')
     do  ig = 1, ng
        if (job1 /= 0) then
           i1 = igv(ig,1)
           i2 = igv(ig,2)
           i3 = igv(ig,3)
        endif
        if (job2 /= 0) then
           i1 = igv2(1,ig)
           i2 = igv2(2,ig)
           i3 = igv2(3,ig)
        endif
        do  i = 1, 3
           qpg(i)= q(i) + qlat(i,1)*i1 + qlat(i,2)*i2 + qlat(i,3)*i3
        enddo
        q2 = (qpg(1)**2+qpg(2)**2+qpg(3)**2) *tpiba**2
        write(stdo,334)  ig,i1,i2,i3,q2
        write(stdo,"('q  qpg=',3d13.5,3x,3d13.5)")  q,qpg
334     format(i5,1x,3i5,2x,f8.4)
     enddo
  endif
end subroutine gvlst2
subroutine gvlstn(q0,q1,q2,qp,mshlst,gmax0,nn)
  implicit none
  intent(in) ::   q0,q1,q2,qp,mshlst,gmax0
  intent(out)::                            nn
  !- Multiples of r.l.v. that bound cutoff gmax0  (r.l.v. reciplocal lattice vector)
  ! ----------------------------------------------------------------------
  !i Inputs
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
  !u Updates
  !u   15 Apr 05 Bug fix when nn > max mshlst
  ! ----------------------------------------------------------------------
  integer :: nmin,nmax,nn,mshlst(0:*), indx
  double precision :: q0(3),q1(3),q2(3),qp(3),gmax0,qperp(3),ddot,qqperp
  ! ... qperp = q1 x q2 / |q1 x q2| ; qqperp = q . qperp
  qperp(1)  = q1(2)*q2(3) - q1(3)*q2(2)
  qperp(2)  = q1(3)*q2(1) - q1(1)*q2(3)
  qperp(3)  = q1(1)*q2(2) - q1(2)*q2(1)
  call dscal(3,1/sqrt(ddot(3,qperp,1,qperp,1)),qperp,1)
  qqperp = ddot(3,q0,1,qperp,1)
  nmax =  gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp + 1
  nmin = -gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp - 1
  nn = 2*max(iabs(nmin),nmax)+1
  if (mshlst(0) /= 0) then
     indx = 1
     call hunti(mshlst(1),mshlst,nn,0,indx)
     nn = mshlst(min(indx+1,mshlst(0)))
  endif
end subroutine gvlstn
subroutine gvlsts(ngxx,ng,gv,kv,igv,igv2,job1,job2)
  implicit none
  integer:: ng,kv(ng,3),igv(ng,3),igv2(3,ng),job1,job2,ig,m,jg,iprm(ng),ngxx
  real(8):: gv(ng,3)
  call dvshel(1,ng, sum(gv**2,dim=2)*[((1d0 + 1d-15*ig),ig=1,ng)], iprm,1)
  gv(:,:) = gv(iprm+1,:)
  kv(:,:) = kv(iprm+1,:)
  if(job1/=0) igv(:,:) =igv(iprm+1,:) 
  if(job2/=0) igv2(:,:)=igv2(iprm+1,:)
end subroutine gvlsts
subroutine gvgetf(ng,n,kv,k1,k2,k3,c,c0)!- Gathers Fourier coefficients from 3D array c into list c0.
  implicit none
  integer :: ng,n,k1,k2,k3,kv(ng,3)
  complex(8):: c0(ng,n),c(k1,k2,k3,n)
  integer :: ig,i,j1,j2,j3
  do i=1,n
     c0(:,i) = [(c(kv(ig,1),kv(ig,2),kv(ig,3),i), ig=1,ng)]
  enddo   
end subroutine gvgetf
subroutine gvputf(ng,n,kv,k1,k2,k3,c0,c)!- Pokes Fourier coefficients from list c0 into 3D array c.
  implicit none
  integer :: ng,n,k1,k2,k3,kv(ng,3)
  complex(8):: c0(ng,n),c(k1,k2,k3,n)
  integer :: ig,i,j1,j2,j3
  c=0d0
  do ig=1,ng
     c(kv(ig,1),kv(ig,2),kv(ig,3),:) = c0(ig,:)
  enddo   
end subroutine gvputf
subroutine gvaddf(ng,kv,k1,k2,k3,c0,c)! Adds Fourier coefficients from list c0 into 3D array c.
  implicit none
  integer :: ng,k1,k2,k3,kv(ng,3)
  complex(8):: c0(ng),c(k1,k2,k3)
  integer :: ig,j1,j2,j3
  do  10  ig = 1, ng
     j1 = kv(ig,1)
     j2 = kv(ig,2)
     j3 = kv(ig,3)
     c(j1,j2,j3) = c(j1,j2,j3) + c0(ig)
10 enddo
end subroutine gvaddf
subroutine mshsiz(alat,plat,job,gmax,ngabc,ng)
  use m_lgunit,only:stdo
  !- Finds dimensions for a mesh of G vectors that satisfy a cutoff
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   plat  :primitive lattice vectors, in units of alat
  !i   job   :0 only change input gmax if input is zero.
  !i         :1 set output gmax to that generated from ngabc
  !i         :not implemented ... probably doesn't make sense
  ! o Inputs/Outputs
  ! o  gmax  :On input, cutoff in reciprocal lattice vectors
  ! o        :Energy cutoff is gmax**2.
  ! o        :If input gmax is zero, it will be generated it from ngabc
  ! o        :(It is an error for both gmax and any of ngabc to be zero.)
  ! o        :On output, gmax is computed if input n1..n3 are nonzero,
  ! o        :or if job = 1
  ! o  ngabc :On input, max # divisions along the three recip lattice vecs.
  ! o        :(It is an error for both gmax and any of ngabc to be zero.)
  ! o        :Otherwise, those input ngabc initially zero are found that
  ! o        :satisfy cutoff criterion.
  !o Outputs
  !o   ng    :number of G-vectors
  !r Remarks
  !b Bugs
  !b   job not implemented.
  !u Updates
  !u   31 Jul 06 fmax increased 160 -> 600, a warning added
  !u   15 Apr 05 Bug fix when ngabc(i)=2
  !u   01 Jun 01 Redesigned with new gvlist
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ngabc(3),ng,job
  double precision :: alat,plat(3,3),gmax
  ! ... Local parameters
  logical :: change
  integer :: npfac,pfac(10),fmax,i,indx,iprint,i1,i2,i3,nn,nmin,nmx(3),nmxn(3),k,ngabcn(3),ngn,nginit,PRTG,kxx(1,1)
  double precision :: gmaxn,xx,q(3),gmax0,qlat(3,3),tpiba,facg,tolg,s1,s2,s3,ddot,gxx(1,1)
  character(256) :: outs
  parameter (fmax=600,tolg=1d0,PRTG=30)
  integer :: mshlst(0:3*fmax)
  logical:: fullmesh
  !     print *,'mshsiz: gmax=',gmax
  !      stdo = lgunit(1)
  call dpzero(q,3)
  ! ... If input gmax is zero assume n1..n3 available, and return w/ gmax
  if (gmax == 0) then
     call pshpr(iprint()-20)
     call gvctof(0,alat,plat,q,ngabc(1),ngabc(2),ngabc(3),gmax,ng)
     call poppr
     return
  endif
  !      print *,'mmm2 mshsiz:'
  ! ... gmax0 = dimensionless gmax
  tpiba = 2*4d0*datan(1d0)/alat
  gmax0  = gmax/tpiba
  call dinv33(plat,1,qlat,xx)
  ! ... list of all allowed values of n1..n3
  mshlst(0) = 0
  call gtpfac(npfac,pfac)
  call ppfac(npfac,pfac,fmax,1,mshlst(1),mshlst)
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
  ! ... Count the number of G vectors for initial n1..n3
  call pshpr(iprint()-30)
  !      call gvlist(alat,plat,q,ngabc(1),ngabc(2),ngabc(3),gmax,000,0,ng,xx,xx,xx,xx)
  call gvlst2(alat,plat,q,ngabc(1),ngabc(2),ngabc(3), 0d0,gmax,[0],000, 0,ng,kxx,gxx,kxx,kxx)
  !     Hold on to the upper bound to know when we fall below tolerance
  nginit = ng
  ! akao
  !      if(fullmesh()) then
  !        print *, ' Full mesh mode: taka all G for charge density. n1*n2*n3'
  !        ng = ngabc(1)*ngabc(2)*ngabc(3)
  !      endif

  ! ... Reduce gmax slowly until at least one of n1..n3 changes
10 continue
  do  k = 1, 99
     gmaxn = gmaxn * facg
     change = .false.
     do  i = 1, 3
        i2 = mod(i,3)+1
        i3 = mod(i2,3)+1
        call gvlstn(qlat(1,i),qlat(1,i2),qlat(1,i3),q,mshlst,gmaxn,nn) !,    nmin,nmxn(i))
        nmxn(i)=nn
        if (nmxn(i) /= nmx(i)) change = .TRUE. 
        !       The granularity of gvlstn may be too coarse.
        !       Don't assign, ngabcn(i) = nn but find next one smaller in mshlst
        indx = 1
        call hunti(mshlst(1),mshlst,ngabc(i),0,indx)
        indx = max(indx,1)
        ngabcn(i) = mshlst(indx)
     enddo
     if (change) goto 12
  enddo
12 continue
  ! akao
  ! 18   continue
  if(fullmesh()) then
     ngabc=ngabcn
     gmax=1d10
     ng=ngabc(1)*ngabc(2)*ngabc(3)
     !       print *,' zzz mshsiz ngn=',ngn
     goto 21
  endif
  ! ... Count the number of G vectors for (smaller) trial n1..n3
  ! cccccccccccccccccccccc
  !      print *,' uuuuuuuu mshsiz ngabcn=',ngabcn(1),ngabcn(2),ngabcn(3)
  !      call gvlist(alat,plat,q,ngabcn(1),ngabcn(2),ngabcn(3),gmax,000,0,ngn,xx,xx,xx,xx)
  call gvlst2(alat,plat,q,ngabcn(1),ngabcn(2),ngabcn(3), 0d0,gmax,[0],000, 0,ngn,kxx,gxx,kxx,kxx)
  !      print *,' uuuuuuuu mshsiz ngn=',ngn
  ! ccccccccccccccccccccccc
  !     If the G vector count doesn't fall by more than tolg, use trial
  if (dble(ngn) >= nginit*tolg) then
     ng = ngn
     call icopy(3,ngabcn,1,ngabc,1)
     goto 10
  endif

  ! ... Check if fmax is too small
21 continue !takao
  call poppr
  i1 = ngabc(1)
  i2 = ngabc(2)
  i3 = ngabc(3)
  if (max0(i1,i2,i3) == fmax .AND. iprint() >= 10) then
     write(stdo,301)
     write(stdo,300) fmax
     write(stdo,301)
  endif
300 format(' WARNING!'/' At least one of the mesh divisions ', &
       'reached its maximal value fmax = ',i4/ &
       ' You might need to increase parameter fmax in mshsiz')
301 format(/1x,79('*')/)
  ! ... Printout
  if (iprint() >= PRTG) then
     s1 = alat*sqrt(ddot(3,plat(1,1),1,plat(1,1),1))/i1
     s2 = alat*sqrt(ddot(3,plat(1,2),1,plat(1,2),1))/i2
     s3 = alat*sqrt(ddot(3,plat(1,3),1,plat(1,3),1))/i3
     write(stdo,"('MSHSIZ: mesh has ',i0,' x ',i0,' x ',i0, &
          ' divisions; length =',3f10.3)")i1,i2,i3,s1,s2,s3
     write(stdo,"('      generated from gmax (a.u.)=',f12.4,': ',i0,' vectors of ', &
          i0,' (',i0,'%)')")  gmax,ng,i1*i2*i3,(ng*100)/(i1*i2*i3)
  endif
end subroutine mshsiz

subroutine gtpfac(npfac,pfac)
  !- Returns allowed prime factors in integers for uniform mesh
  !     implicit none
  integer :: npfac,pfac(5)
  npfac = 5
  pfac(1) = 2
  pfac(2) = 3
  pfac(3) = 5
  pfac(4) = 7
  pfac(5) = 11
end subroutine gtpfac

subroutine ppfac(npfac,pfac,fmax,job,fac,nfac)
  !- Find all products of prime factors within some maximum
  ! ----------------------------------------------------------------------
  !i Inputs
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
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  integer :: npfac,pfac(npfac),fac(*),fmax,job,nfac
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
10      continue
     enddo
  enddo
  if (job /= 0) call ivheap(1,nfac,fac,fac(1+2*nfac),0)
end subroutine ppfac

!end module m_gvlst2
