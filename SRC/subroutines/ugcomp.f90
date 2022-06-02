subroutine ugcomp(nbas,ssite,sspec,qmom,gpot0,hpot0,ugg,f) !slat,
  use m_struc_def           !Cgetarg
  use m_lgunit,only:stml
  !- Part of the smooth estatic energy from compensating G's alone.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct containing site-specific information
  !i   sspec :struct containing species-specific information
  !i   slat  :struct containing information about the lattice
  !i   qmom  :multipole moments of on-site densities (rhomom.f)
  ! o Inputs/Outputs
  ! o  Let n0  = smooth potential without compensating gaussians
  ! o      n0~ = smooth potential with compensating gaussians
  ! o    phi0  = ves[n0]
  ! o    phi0~ = ves[n0~]
  ! o    g_RL  = gaussian in RL channel
  ! o    h_R   = l=0 sm hankel in RL channel, for core density
  ! o  Then:
  ! o  gpot0 :On input, integrals g_RL * phi0
  ! o        :On output, integrals g_RL * phi0~
  ! o  hpot0 :On input, integrals h_R * phi0
  ! o        :On output, integrals h_R * phi0~
  !o Outputs
  !o   ugg   :electrostatic energy integral [n0~-n0]*[phi0~-phi0]
  !i   f     :contribution to forces is added
  !r Remarks
  !u Updates
  !u   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   22 Apr 00 Adapted from nfp ugcomp
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
! #if MPI | MPIK
!   include "mpif.h"
!   integer :: procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
!   integer :: MAX_PROCS
!   parameter (MAX_PROCS = 100)
!   integer :: resultlen
!   character*(MPI_MAX_PROCESSOR_NAME) name
!   character(10) :: shortname(0:MAX_PROCS-1)
!   character(20) :: ext
!   character(26) :: datim
!   integer :: namelen(0:MAX_PROCS-1)
!   double precision :: starttime, endtime
!   logical :: mlog,cmdopt
!   character(120) :: strn
! #endif
  integer :: nbas
  real(8):: qmom(*) , gpot0(*) , f(3,nbas) , hpot0(nbas) , ugg
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  !      type(s_lat)::slat

  ! ... Local parameters
  integer :: ndim,ndim0,i,ib,ilm1,ilm2,is,iv0,jb,js,jv0,nvl,l1,l2, &
       lfoc1,lfoc2,ll,lmax1,lmax2,m,nlm1,nlm2
  parameter (ndim=49, ndim0=2)
  double precision :: ceh1,ceh2,cof1,cof2,cofg1,cofg2,cofh1,cofh2,fpi, &
       pi,qcorg1,qcorg2,qcorh1,qcorh2,qsc1,qsc2,qm1,qm2,rg1,rg2,rh1, &
       rh2,srfpi,y0,z1,z2
  double precision :: df(0:20),ff(3),tau1(3),tau2(3)
  double complex s(ndim,ndim),ds(ndim,ndim,3),s0(ndim0,ndim0), &
       ds0(ndim0,ndim0,3),wk(ndim0,ndim0),dwk(ndim0,ndim0,3)
  ! ... For parallel threads
  integer :: nlmx,npmx,ip,mp,nbmx


  ! #ifndef SGI_PARALLEL
  parameter (nlmx=64, npmx=1, nbmx=256)
  ! FCPP#if F90 | AUTO_ARRAY
  double precision :: xf(3,nbas,npmx),xhpot0(nbas,npmx), &
       xgpot0(nlmx*nbas,npmx),xugg(npmx)
  ! FCPP#else
  ! FCPP      double precision xf(3,nbmx,npmx),xhpot0(nbmx,npmx),
  ! FCPP     .xgpot0(nlmx*nbmx,npmx),xugg(npmx)
  ! FCPP#endif
  ! #else
  !      parameter (nlmx=64, npmx=32)
  !      double precision xf(3,nbas,npmx),xhpot0(nbas,npmx),
  !     .  xgpot0(nlmx*nbas,npmx),xugg(npmx)
  ! #endif



! #if MPI | MPIK
!   integer :: , dimension(:), allocatable :: bproc
!   double precision :: , dimension(:), allocatable :: buffer
!   integer :: nvl0,iiv0(nbas)
! #endif

  ! ... Heap

  integer:: ibini,ibend
  call tcn('ugcomp')

! #if MPI | MPIK
!   call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
!   call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
!   call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
!   call strcop(shortname(procid),name,10,'.',i)
!   namelen(procid) = i-1
!   master = 0
!   mlog = cmdopt('--mlog',6,0,strn)
! #endif

  call stdfac(20,df)
  pi = 4d0*datan(1d0)
  fpi = 4d0*pi
  srfpi = dsqrt(fpi)
  y0 = 1d0/srfpi
  ! #if ! (F90 | AUTO_ARRAY | SGI_PARALLEL)
  ! FCPP#if ! (F90 | AUTO_ARRAY)
  ! FCPP      if (nbas .gt. nbmx) call rx('ugcomp: increase nbkmx')
  ! FCPP#endif

  ! ... Setup array iiv0 = (vector of iv0 for parallel); allocate work arrays
  mp = 1

! #if MPI | MPIK
!   call setofl(0,ssite,sspec,nbas,nvl0,iiv0)
!   if (nlmx*nbas < nvl0) call rx('ugcomp: increase nlmx')
! #endif

  if (npmx < mp) call rxi('ugcomp: increase npmx, needed',mp)

  ! --- Loop over sites where charge lump making pot is centered ---
  ugg = 0d0
  iv0 = 0
  ip = 1
  call dpzero(xugg, mp)
  call dpzero(xgpot0, nlmx*nbas*mp)
  call dpzero(xf, 3*nbas*mp)
  call dpzero(xhpot0, nbas*mp)
! #if MPI | MPIK
!   allocate (bproc(0:numprocs), stat=ierr)
!   call dstrbp(nbas,numprocs,1,bproc(0))
!   !      do  ib = bproc(procid), bproc(procid+1)-1
!   ibini= bproc(procid)
!   ibend= bproc(procid+1)-1
! #else
  !      do  ib = 1, nbas
  ibini=1
  ibend=nbas
!#endif

  do ib=ibini,ibend
! #if MPI |MPIK
!      if (mlog .AND. ib == bproc(procid)) then
!         call gettime(datim)
!         write(stml,"(a,i5,' of ',i5,a,i5,' to ',i5)")' ugcomp '//datim//' Process ',procid,numprocs, &
!              ' on '//shortname(procid)(1:namelen(procid))//' starting atoms ',bproc(procid),bproc(procid+1)-1
!      endif
!      iv0 = iiv0(ib)
! #endif

     is=ssite(ib)%spec
     !        i_copy_size=size(ssite(ib)%pos)
     !        call dcopy(i_copy_size,ssite(ib)%pos,1,tau1,1)
     tau1=ssite(ib)%pos
     lmax1=sspec(is)%lmxl
     rg1=sspec(is)%rg

     call corprm(sspec,is,qcorg1,qcorh1,qsc1,cofg1,cofh1,ceh1,lfoc1, &
          rh1,z1)
     nlm1 = (lmax1+1)**2

     !   ... Loop over sites where charge lump sees the potential
     if (lmax1 > -1) then
        jv0 = 0
        do  jb = 1, nbas
           js=ssite(jb)%spec
           !            i_copy_size=size(ssite(jb)%pos)
           !            call dcopy(i_copy_size,ssite(jb)%pos,1,tau2,1)
           tau2=ssite(jb)%pos
           lmax2=sspec(js)%lmxl
           rg2=sspec(js)%rg

           if (lmax2 > -1) then
              call corprm(sspec,js,qcorg2,qcorh2,qsc2,cofg2,cofh2,ceh2, &
                   lfoc2,rh2,z2)
              nlm2 = (lmax2+1)**2
! #if MPI | MPIK
!               jv0 = iiv0(jb)
! #endif
              if (nlm1 > ndim) call rxi('ugcomp: ndim < nlm1=',nlm1)
              if (nlm2 > ndim) call rxi('ugcomp: ndim < nlm2=',nlm2)
              call ggugbl(tau1,tau2,rg1,rg2,nlm1,nlm2,ndim,ndim,s,ds) !,slat

              ff(1) = 0d0
              ff(2) = 0d0
              ff(3) = 0d0
              do  ilm1 = 1, nlm1
                 l1 = ll(ilm1)
                 qm1 = qmom(iv0+ilm1)
                 if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
                 cof1 = qm1*fpi/df(2*l1+1)
                 do  ilm2 = 1, nlm2
                    l2 = ll(ilm2)
                    qm2 = qmom(jv0+ilm2)
                    if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
                    cof2 = qm2*fpi/df(2*l2+1)
                    xugg(ip) = xugg(ip) + cof1*cof2*s(ilm1,ilm2)
                    xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip) &
                         + s(ilm1,ilm2)*cof1*fpi/df(2*l2+1)
                    !         ... Forces
                    ff(1) = ff(1) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,1)
                    ff(2) = ff(2) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,2)
                    ff(3) = ff(3) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,3)
                 enddo
              enddo

              !     --- Additional h*h, h*g, g*h terms for foca ---
              if (lfoc1 > 0 .OR. lfoc2 > 0) then
                 call hhugbl(0,tau1,tau2,rh1,rh2,ceh1,ceh2,1,1,ndim0,ndim0, &
                      wk,dwk,s0,ds0) !slat,
                 xugg(ip) = xugg(ip) + cofh1*s0(1,1)*cofh2
                 xhpot0(jb,ip) = xhpot0(jb,ip) + cofh1*s0(1,1)
                 ff(1) = ff(1) + 0.5d0*cofh1*cofh2*ds0(1,1,1)
                 ff(2) = ff(2) + 0.5d0*cofh1*cofh2*ds0(1,1,2)
                 ff(3) = ff(3) + 0.5d0*cofh1*cofh2*ds0(1,1,3)

                 call hgugbl(tau1,tau2,rh1,rg2,ceh1,1,nlm2,ndim,ndim, &
                      s,ds) !slat,
                 do  ilm2 = 1, nlm2
                    l2 = ll(ilm2)
                    qm2 = qmom(jv0+ilm2)
                    if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
                    cof2 = qm2*fpi/df(2*l2+1)
                    xugg(ip) = xugg(ip) + cofh1*s(1,ilm2)*cof2
                    ff(1) = ff(1) + 0.5d0*cofh1*cof2*ds(1,ilm2,1)
                    ff(2) = ff(2) + 0.5d0*cofh1*cof2*ds(1,ilm2,2)
                    ff(3) = ff(3) + 0.5d0*cofh1*cof2*ds(1,ilm2,3)
                    xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip) &
                         + s(1,ilm2)*cofh1*fpi/df(2*l2+1)
                 enddo

                 call hgugbl(tau2,tau1,rh2,rg1,ceh2,1,nlm1,ndim,ndim, &
                      s,ds) !slat,
                 do  ilm1 = 1, nlm1
                    l1 = ll(ilm1)
                    qm1 = qmom(iv0+ilm1)
                    if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
                    cof1 = qm1*fpi/df(2*l1+1)
                    xugg(ip) = xugg(ip) + cof1*s(1,ilm1)*cofh2
                    ff(1) = ff(1) - 0.5d0*cof1*cofh2*ds(1,ilm1,1)
                    ff(2) = ff(2) - 0.5d0*cof1*cofh2*ds(1,ilm1,2)
                    ff(3) = ff(3) - 0.5d0*cof1*cofh2*ds(1,ilm1,3)
                    xhpot0(jb,ip) = xhpot0(jb,ip) + cof1*s(1,ilm1)
                 enddo
              endif

              if (jb /= ib) then
                 do  m = 1, 3
                    xf(m,ib,ip) = xf(m,ib,ip) - ff(m)
                    xf(m,jb,ip) = xf(m,jb,ip) + ff(m)
                 enddo
              endif

              jv0 = jv0+nlm2
           endif
        enddo
        iv0 = iv0+nlm1
     endif
  enddo
! #if MPI | MPIK
!   nvl = nvl0
! #else
  nvl = iv0
!#endif

  ! ... Assemble data from separate threads
! #if MPI | MPIK
!   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   allocate(buffer(1:nvl), stat=ierr)
!   call MPI_ALLREDUCE(xgpot0,buffer,nvl, &
!        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   if (mlog) then
!      call gettime(datim)
!      write(stml,"(a,i5,' of ',i5,a,i5)") ' ugcomp '//datim//' Process ',procid,numprocs, &
!           ' on '//shortname(procid)(1:namelen(procid))//' allreduce gpot0 nvl=',nvl
!   endif
!   call daxpy(nvl,1d0,buffer,1,gpot0,1)
!   deallocate(buffer, stat=ierr)

!   allocate(buffer(1:nbas), stat=ierr)
!   call MPI_ALLREDUCE(xhpot0,buffer,nbas, &
!        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   if (mlog) then
!      call gettime(datim)
!      write(stml,"(a,i5,' of ',i5,a,i5)") &
!           ' ugcomp '//datim//' Process ',procid,numprocs, &
!           ' on '//shortname(procid)(1:namelen(procid))//' allreduce hpot0 nbas=',nbas
!   endif
!   call daxpy(nbas,1d0,buffer,1,hpot0,1)
!   deallocate(buffer, stat=ierr)

!   allocate(buffer(1:3*nbas), stat=ierr)
!   call MPI_ALLREDUCE(xf,buffer,3*nbas, &
!        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   if (mlog) then
!      call gettime(datim)
!      write(stml,"(a,i5,' of ',i5,a,i5)") &
!           ' ugcomp '//datim//' Process ',procid,numprocs, &
!           ' on '//shortname(procid)(1:namelen(procid))//' allreduce f 3nbas=',3*nbas
!   endif
!   call daxpy(3*nbas,1d0,buffer,1,f,1)
!   deallocate(buffer, stat=ierr)

!   call MPI_ALLREDUCE(xugg,ugg,1, &
!        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   if (mlog) then
!      call gettime(datim)
!      write(stml,"(a,i5,' of ',i5,a)") &
!           ' ugcomp '//datim//' Process ',procid,numprocs, &
!           ' on '//shortname(procid)(1:namelen(procid))//' allreduce ugg'
!   endif
!   deallocate(bproc, stat=ierr)
! #else
  do  80  ip = 1, mp
     do  82  ib = 1, nbas
        f(1,ib) = f(1,ib) + xf(1,ib,ip)
        f(2,ib) = f(2,ib) + xf(2,ib,ip)
        f(3,ib) = f(3,ib) + xf(3,ib,ip)
        hpot0(ib) = hpot0(ib) + xhpot0(ib,ip)
82   enddo
     do  84  i = 1, nvl
        gpot0(i) = gpot0(i) + xgpot0(i,ip)
84   enddo
     ugg = ugg + xugg(ip)
80 enddo
!#endif
  call tcx('ugcomp')
end subroutine ugcomp


