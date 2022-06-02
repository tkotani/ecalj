! FCPP#define F90 1
subroutine gvlst2(alat,plat,q,n1,n2,n3,gmin,gmax,mshlst,job,ngmx, &
     ng,kv,gv,igv,igv2)
  use m_ftox
  use m_shortn3,only: shortn3_initialize,shortn3
  use m_lgunit,only:stdo
  !- Set up a list of recip vectors within cutoff |q+G| < gmax
  ! ----------------------------------------------------------------------
  !i Inputs
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
  !i            100s digit
  !i              1 to return internally generated values for n1,n2,n3;
  !i                see description of n1,n2,n3 below.
  !i            1000s digit
  !i              2  not used ... OLD Do not change input gmax if it is nonzero.
  !i   ngmx     Leading dimension of kv,gv,igv
  !i   gmin     Lower bound for reciprocal lattice vectors, in a.u.
  ! o Inputs/Outputs
  ! o   gmax    On input, cutoff for reciprocal lattice vectors, in a.u.
  ! o           Energy cutoff is gmax**2.
  ! o           If input gmax is zero, gvlst2 will generate it from n1..n3
  ! o           (It is an error for both gmax and n1..n3 to be zero.)
  ! o           On output, gmax may be altered; see Remarks.
  ! o  n1..3    On input, max # divisions along the three lattice vectors.
  ! o           (It is an error for both gmax and n1..n3 to be zero.)
  ! o           Otherwise, input n1..n3 additionally constrain which
  ! o           vectors are added to list; see Remarks.
  ! o           On output, any n1..n3 initially zero are found from gmax
  !o Outputs
  !o   ng       Number of lattice vectors
  !o   gv       list of reciprocal lattice vectors G
  !o   igv      list of reciprocal lattice vectors G, represented as
  !o            three integers (the multiples of qlat)
  !o            gv and igv are related by:
  !o              gv(1:3,1:ng) = 2*pi/alat * (qlat * igv(1:ng))
  !o   igv2     same as igv except first and second columns are permuted
  !o   kv       indices for gather/scatter operations.
  !o            kv(ig,i=1,2,3) for vector ig point to which entry
  !o            (i1,i2,i3) vector ig belongs
  !r Remarks
  !r   Collects a list of q + reciprocal lattice vectors (G+q) that lie
  !r   within a cutoff gmax.  List is optionally sorted by length.
  !r   Vectors G are integer multiples of primitive r.l.v. qlat.
  !r
  !r   Cutoff gmax may be input (preferred mode), or generated from
  !r   input values n1..n3.
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
  !u Updates
  !u   11 Jul 08 New argument gmin
  !u   01 Jun 01 revised gvlst2 together with gvlist.  They form similar
  !u             operations but with different functions; see gvlist.f
  !u   26 Mar 01 Another bug fix for input n1..n3 and gmax nonzero
  !u   06 Mar 01 Bug fix for input n1..n3 and gmax nonzero
  !u   Routine was adapted from T. Kotani, routine getgv2
  ! ----------------------------------------------------------------------
  implicit none
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
  integer:: j1,j2,j3,noutmx,m,nout,jj1,jj2,jj3,nn1,nn2,nn3
  integer,allocatable:: nlatout(:,:)
  real(8):: rlatp(3,3),xmx2(3)
  !      stdo = lgunit(1)
  call getpr(ipr)
  call dinv33(plat,1,qlat,vol)
  pi = 4d0*datan(1d0)
  tpiba = 2*pi/alat
  gmin0  = gmin/tpiba
  gmax0  = gmax/tpiba
  if (gmin < 0) call rx('gvlst2: input gmin <= 0')
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
  if(nn1==0) call gvlstn(qlat(1,1),qlat(1,2),qlat(1,3),q,mshlst,gmax0,nn1) !,n1min,n1max)
  if(nn2==0) call gvlstn(qlat(1,2),qlat(1,3),qlat(1,1),q,mshlst,gmax0,nn2) !,n2min,n2max)
  if(nn3==0) call gvlstn(qlat(1,3),qlat(1,1),qlat(1,2),q,mshlst,gmax0,nn3) !,n3min,n3max)
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
  !      print *,' jjjjjjj qlat qlat1=',qlat,qlat1
  !! --- Loop through g vectors, shorten, count and keep if within gmax ---
  ig = 0
  gmax2 = (gmax0-tol)**2
  gmin2 = gmin0**2
  noutmx=48
  call shortn3_initialize(qlat1) !initialization for shoten3
  allocate(nlatout(3,noutmx))
  !      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  n1l =',n1l,n1min,n1max
  !      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  n2l =',n2l,n2min,n2max
  !      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  n3l =',n3l,n3min,n3max
  !      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  gmax2 =',gmax2
  do  212  j1 = 0,nn1-1 !n1min, n1max
     do  211  j2 = 0,nn2-1 !n2min, n2max
        do  21  j3 = 0,nn3-1 !n3min, n3max
           qpg= (/j1, j2, j3/) + matmul(q,plat(:,:))
           qpg= (/ qpg(1)/dble(nn1), qpg(2)/dble(nn2), qpg(3)/dble(nn3) /)
           ! pg on the supercell corrdinate of qlat1.
           call shortn3(qpg, noutmx, nout,nlatout)
           gs= matmul(qlat1(:,:), (qpg+nlatout(:,1)))
           gg = (gs(1)**2+gs(2)**2+gs(3)**2)
           ! cccccccccccccccccccccccccccccccccccccc
           !          write(*,"(a,2d13.5,3i5,3d13.5)") 'gggggggggggg gg=',gg,gmax2,nlatout(:,1),sum(abs(qlat1)),sum(abs(qlat1)),sum(abs(qpg))
           ! cccccccccccccccccccccccccccccccccccccc
           ! qpg+nlatout(:,ix) are shortest on the corrdinate of
           ! qlat1 (number of vectors are nlatout).
           ! We use nlatout(1:3,1).
           !$$$ccccccccccccccccccccccccccccccccccccccccccccccccccc
           !$$$          print *,'jjjjj nlatout',j1,j2,j3,nout,nlatout(:,1)
           !$$$            do  i = 1, 3
           !$$$              qpg(i)= q(i) + qlat(i,1)*j1 + qlat(i,2)*j2 + qlat(i,3)*j3
           !$$$            enddo
           !$$$            gg = qpg(1)**2+qpg(2)**2+qpg(3)**2
           ! cccccccccccccccccccccccccccccccccccccccccccccc
           k1 = j1+1 !mod(j1+n1l,n1l)+1
           k2 = j2+1 !mod(j2+n2l,n2l)+1
           k3 = j3+1 !mod(j3+n3l,n3l)+1
           !          write(*,"('wwwww1111 j1j2j3=',7i4,f13.5)")ig,j1,j2,j3,k1,k2,k3,gg
           if(gmin2<= gg .AND. gg < gmax2) then
              ig = ig+1
              !            print *,' zzz ig=',ig,ngmx
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
  deallocate(nlatout)
  ng = ig
  ! --- Printout ---
  if (ipr >= PRTG .AND. n1l*n2l*n3l == 0) then
     write(stdo,ftox)' GVLST2: gmax=',ftof(gmax0*tpiba,3),'created',ng,' recip. lattice vectors'
  elseif (ipr >= PRTG) then
     h1 = alat*sqrt(ddot(3,plat(1,1),1,plat(1,1),1))/n1l
     h2 = alat*sqrt(ddot(3,plat(1,2),1,plat(1,2),1))/n2l
     h3 = alat*sqrt(ddot(3,plat(1,3),1,plat(1,3),1))/n3l
     write(stdo,ftox)' GVLST2: gmax = ',ftof(gmax,3), &
          ' a.u. created ',ng,' vectors of ',n1l*n2l*n3l,'(',(ng*100)/(n1l*n2l*n3l),'%)'
     !$$$        iset(1) = n1
     !$$$        iset(2) = n2
     !$$$        iset(3) = n3
     !$$$        i = iset(1)*iset(2)*iset(3)
     !$$$        call awrit7('%a%N%9f%?#n#(input) ##mesh has %i x %i x %i'//
     !$$$     .  ' divisions; length %,3;3d, %,3;3d, %,3;3d',outs,
     !$$$     .  len(outs),0,i,n1l,n2l,n3l,h1,h2,h3)
     !$$$        if (i .eq. 0 .and. iset(1)**2+iset(2)**2+iset(3)**2 .ne. 0) then
     !$$$          call awrit3('%a%N%9fgenerated from input mesh with ('//
     !$$$     .    '%?#n#%-1j%i#*#,%?#n#%-1j%i#*#,'//
     !$$$     .    '%?#n#%-1j%i#*#) divisions',
     !$$$     .    outs,len(outs),0,iset(1),iset(2),iset(3))
     !$$$        endif
     !$$$        call awrit0('%a',outs,len(outs),-stdo)
  endif
  !! --- Sort the list of vectors --
  if(lsort) then
     call gvlsts(ngmx,ng,gv,kv,igv,igv2,job1,job2)
  endif
  if (ipr >= PRTG2 .AND. ng > 0 .AND. job1+job2 /= 0) then
     write(stdo,333)
333  format(' G vectors (multiples of reciprocal lattice vectors)'/ &
          '   ig    G1   G2   G3     E')
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

subroutine gvlstn(q0,q1,q2,qp,mshlst,gmax0,nn) !,nmin,nmax)
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
  !i Inputs/Outputs
  !i   nn    :On input, if nonzero, maximum number of mesh points
  !i         :allowed.  Thus upper limit of nmax is nn/2; lower limit of nmin is -nn/2.
  !i         :On output, if input is zero, nn = 2*max(|nmin|,nmax)+1
  ! o  gmax0 :On input, cutoff G
  ! o        :On output, gmax0 may be reduced because constraints
  ! o        :on nmin,nmax cause lattice vectors to be neglected
  ! o        :that are smaller than input gmax0 (if input nn nonzero).
  !o Outputs
  !o     nn :   meshsize
  ! xo   nmin  :search for lattice vectors limited to (nmin..nmax)*q0
  ! xo   nmax  :search for lattice vectors limited to (nmin..nmax)*q0
  !r Remarks
  !r   q0,q1,q2,qp and G are all dimensionless (units of 2*pi/a)
  !u Updates
  !u   15 Apr 05 Bug fix when nn > max mshlst
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nmin,nmax,nn,mshlst(0:*)
  double precision :: q0(3),q1(3),q2(3),qp(3),gmax0
  ! ... Local parameters
  integer :: indx
  double precision :: qperp(3),ddot,qqperp
  ! ... qperp = q1 x q2 / |q1 x q2| ; qqperp = q . qperp
  qperp(1)  = q1(2)*q2(3) - q1(3)*q2(2)
  qperp(2)  = q1(3)*q2(1) - q1(1)*q2(3)
  qperp(3)  = q1(1)*q2(2) - q1(2)*q2(1)
  !     anorm = sqrt(ddot(3,q1,1,q1,1))
  !     bnorm = sqrt(ddot(3,q2,1,q2,1))
  !     call dscal(3,1/anorm/bnorm,qperp,1)
  call dscal(3,1/sqrt(ddot(3,qperp,1,qperp,1)),qperp,1)
  qqperp = ddot(3,q0,1,qperp,1)
  !   10 continue
  !     print *, gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp + 1
  !     print *, -gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp - 1
  nmax =  gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp + 1
  nmin = -gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp - 1
  ! ... Assign nn, if input value is zero
  if (nn == 0) then
     nn = 2*max(iabs(nmin),nmax)+1
     if (mshlst(0) /= 0) then
        indx = 1
        call hunti(mshlst(1),mshlst,nn,0,indx)
        nn = mshlst(min(indx+1,mshlst(0)))
     endif
     !      else !bugfix added oct2013
     !         print *,' nnnnn1111 nmax nmin=',nmax,nmin
     !         nmax= min(nmax,  nn/2)
     !         nmin= nmax-nn+1
     !         print *,' nnnnn2222 nmax nmin=',nmax,nmin
  endif
end subroutine gvlstn

subroutine gvlsts(ngmx,ng,gv,kv,igv,igv2,job1,job2)
  !- Kernel called by gvlst2 to sort gv and kv
  !     implicit none
  integer :: ngmx,ng,kv(ngmx,3),igv(ngmx,3),igv2(3,ngmx),job1,job2
  double precision :: gv(ngmx,3)
  ! Local variables
  integer :: ig,m,jg
  ! FCPP#if F90 | AUTO_ARRAY
  integer :: kk(ngmx),iprm(ngmx)
  double precision :: gg(ngmx)
  ! FCPP#else
  ! FCPP      integer ngmxx
  ! FCPP      parameter (ngmxx=20000)
  ! FCPP      integer kk(ngmxx),iprm(ngmxx)
  ! FCPP      double precision gg(ngmxx)
  ! FCPP      if (ng .gt. ngmxx) call rxi('gvlst2: increase ngmx, need',ng)
  ! FCPP#endif
  do  ig = 1, ng
     gg(ig) = gv(ig,1)**2 + gv(ig,2)**2 + gv(ig,3)**2
     ! akao test case1 Apr2009
     gg(ig)= gg(ig) *(1d0 + 1d-15*ig)
  enddo
  call dvshel(1,ng,gg,iprm,1)
  !     call dvheap(1,ng,gg,iprm,0d0,11)
  ! ... Rearrange gv,kv
  do  20  m = 1, 3
     do  22  ig = 1, ng
        jg = iprm(ig)+1
        gg(ig) = gv(jg,m)
        kk(ig) = kv(jg,m)
22   enddo
     do  24  ig = 1, ng
        gv(ig,m) = gg(ig)
        kv(ig,m) = kk(ig)
24   enddo
20 enddo
  ! ... Rearrange igv
  if (job1 /= 0) then
     do  30  m = 1, 3
        do  32  ig = 1, ng
           jg = iprm(ig)+1
           kk(ig) = igv(jg,m)
32      enddo
        do  34  ig = 1, ng
           igv(ig,m) = kk(ig)
34      enddo
30   enddo
  endif
  ! ... Rearrange igv2
  if (job2 /= 0) then
     do  40  m = 1, 3
        do  42  ig = 1, ng
           jg = iprm(ig)+1
           kk(ig) = igv2(m,jg)
42      enddo
        do  44  ig = 1, ng
           igv2(m,ig) = kk(ig)
44      enddo
40   enddo
  endif
  !      do   ig = 1, ng
  !        gg(ig) = gv(ig,1)**2 + gv(ig,2)**2 + gv(ig,3)**2
  !        print 550, ig,gv(ig,1),gv(ig,2),gv(ig,3),
  !     .     kv(ig,1),kv(ig,2),kv(ig,3),sqrt(gg(ig))
  !  550   format(i5,3f11.5,3i6,f11.5)
  !      enddo
  !      pause
end subroutine gvlsts

