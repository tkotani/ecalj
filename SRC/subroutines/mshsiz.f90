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
  integer :: npfac,pfac(10),fmax,i,indx,iprint, &
       i1,i2,i3,nn,nmin,nmx(3),nmxn(3),k,ngabcn(3),ngn,nginit, &
       PRTG
  double precision :: gmaxn,xx,q(3),gmax0,qlat(3,3),tpiba,facg,tolg,s1, &
       s2,s3,ddot
  character(256) :: outs
  !     parameter (fmax=360,tolg=1d0,PRTG=30)
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
  call gvlst2(alat,plat,q,ngabc(1),ngabc(2),ngabc(3), 0d0,gmax,0,000, 0,ng,xx,xx,xx,xx)
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
  call gvlst2(alat,plat,q,ngabcn(1),ngabcn(2),ngabcn(3), 0d0,gmax,0,000, 0,ngn,xx,xx,xx,xx)
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

