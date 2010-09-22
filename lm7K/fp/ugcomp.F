CSFCPP#define F90 1
      subroutine ugcomp(nbas,ssite,sspec,slat,qmom,gpot0,hpot0,ugg,f)
      use m_struc_def  !Cgetarg
#if MPE
      use m_events !Event numbers: include "events.ins"
#endif
C- Part of the smooth estatic energy from compensating G's alone.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Ci   slat  :struct containing information about the lattice
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Cio Inputs/Outputs
Cio  Let n0  = smooth potential without compensating gaussians
Cio      n0~ = smooth potential with compensating gaussians
Cio    phi0  = ves[n0]
Cio    phi0~ = ves[n0~]
Cio    g_RL  = gaussian in RL channel
Cio    h_R   = l=0 sm hankel in RL channel, for core density
Cio  Then:
Cio  gpot0 :On input, integrals g_RL * phi0
Cio        :On output, integrals g_RL * phi0~
Cio  hpot0 :On input, integrals h_R * phi0
Cio        :On output, integrals h_R * phi0~
Co Outputs
Co   ugg   :electrostatic energy integral [n0~-n0]*[phi0~-phi0]
Ci   f     :contribution to forces is added
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   22 Apr 00 Adapted from nfp ugcomp
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI | MPIK
#if MPE
      include "mpef.h"
#endif
      include "mpif.h"
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1),lgunit
      double precision starttime, endtime
      logical mlog,cmdopt
      character*120 strn
#endif
c#if MPE
cC Event numbers:
c      include "events.ins"
c#endif
      integer nbas
      real(8):: qmom(*) , gpot0(*) , f(3,nbas) , hpot0(nbas) , ugg
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

C ... Local parameters
      integer ndim,ndim0,i,ib,ilm1,ilm2,is,iv0,jb,js,jv0,nvl,l1,l2,
     .lfoc1,lfoc2,ll,lmax1,lmax2,m,nlm1,nlm2
      parameter (ndim=49, ndim0=2)
      double precision ceh1,ceh2,cof1,cof2,cofg1,cofg2,cofh1,cofh2,fpi,
     .pi,qcorg1,qcorg2,qcorh1,qcorh2,qsc1,qsc2,qm1,qm2,rg1,rg2,rh1,
     .rh2,srfpi,y0,z1,z2
      double precision df(0:20),ff(3),tau1(3),tau2(3)
      double complex s(ndim,ndim),ds(ndim,ndim,3),s0(ndim0,ndim0),
     .ds0(ndim0,ndim0,3),wk(ndim0,ndim0),dwk(ndim0,ndim0,3)
C ... For parallel threads
      integer nlmx,npmx,ip,mp,nbmx


c #ifndef SGI_PARALLEL
      parameter (nlmx=64, npmx=1, nbmx=256)
CSFCPP#if F90 | AUTO_ARRAY
      double precision xf(3,nbas,npmx),xhpot0(nbas,npmx),
     .xgpot0(nlmx*nbas,npmx),xugg(npmx)
CSFCPP#else
CSFCPP      double precision xf(3,nbmx,npmx),xhpot0(nbmx,npmx),
CSFCPP     .xgpot0(nlmx*nbmx,npmx),xugg(npmx)
CSFCPP#endif
c #else
c      parameter (nlmx=64, npmx=32)
c      double precision xf(3,nbas,npmx),xhpot0(nbas,npmx),
c     .  xgpot0(nlmx*nbas,npmx),xugg(npmx)
c #endif



#if MPI | MPIK
      integer , dimension(:), allocatable :: bproc
      double precision , dimension(:), allocatable :: buffer
      integer nvl0,iiv0(nbas)
#endif

C ... Heap

      integer:: ibini,ibend
      call tcn('ugcomp')

#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#endif

      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      srfpi = dsqrt(fpi)
      y0 = 1d0/srfpi
c #if ! (F90 | AUTO_ARRAY | SGI_PARALLEL)
CSFCPP#if ! (F90 | AUTO_ARRAY)
CSFCPP      if (nbas .gt. nbmx) call rx('ugcomp: increase nbkmx')
CSFCPP#endif

C ... Setup array iiv0 = (vector of iv0 for parallel); allocate work arrays
      mp = 1

#if MPI | MPIK
      call setofl(0,ssite,sspec,nbas,nvl0,iiv0)
      if (nlmx*nbas .lt. nvl0) call rx('ugcomp: increase nlmx')
#endif

      if (npmx .lt. mp) call rxi('ugcomp: increase npmx, needed',mp)

C --- Loop over sites where charge lump making pot is centered ---
      ugg = 0d0
      iv0 = 0
      ip = 1
      call dpzero(xugg, mp)
      call dpzero(xgpot0, nlmx*nbas*mp)
      call dpzero(xf, 3*nbas*mp)
      call dpzero(xhpot0, nbas*mp)
#if MPI | MPIK
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_START_UGCOMP,procid,"ugcomp")
#endif
      allocate (bproc(0:numprocs), stat=ierr)
      call dstrbp(nbas,numprocs,1,bproc(0))
c      do  ib = bproc(procid), bproc(procid+1)-1
      ibini= bproc(procid)
      ibend= bproc(procid+1)-1
#else
c      do  ib = 1, nbas
      ibini=1
      ibend=nbas
#endif

      do ib=ibini,ibend
#if MPI |MPIK
        if (mlog .and. ib .eq. bproc(procid)) then
          call gettime(datim)
          call awrit4(' ugcomp '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' starting atoms %i to %i',' ',256,lgunit(3),
     .    procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        iv0 = iiv0(ib)
#endif

        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pos)
        call dcopy(i_copy_size,ssite(ib)%pos,1,tau1,1)


        lmax1=sspec(is)%lmxl
        rg1=sspec(is)%rg

        call corprm(sspec,is,qcorg1,qcorh1,qsc1,cofg1,cofh1,ceh1,lfoc1,
     .  rh1,z1)
        nlm1 = (lmax1+1)**2

C   ... Loop over sites where charge lump sees the potential
        if (lmax1 .gt. -1) then
          jv0 = 0
          do  jb = 1, nbas

            js=ssite(jb)%spec
            i_copy_size=size(ssite(jb)%pos)
            call dcopy(i_copy_size,ssite(jb)%pos,1,tau2,1)


            lmax2=sspec(js)%lmxl
            rg2=sspec(js)%rg

            if (lmax2 .gt. -1) then
              call corprm(sspec,js,qcorg2,qcorh2,qsc2,cofg2,cofh2,ceh2,
     .        lfoc2,rh2,z2)
              nlm2 = (lmax2+1)**2
#if MPI | MPIK
              jv0 = iiv0(jb)
#endif
              if (nlm1 .gt. ndim) call rxi('ugcomp: ndim.lt.nlm1=',nlm1)
              if (nlm2 .gt. ndim) call rxi('ugcomp: ndim.lt.nlm2=',nlm2)
              call ggugbl(tau1,tau2,rg1,rg2,nlm1,nlm2,ndim,ndim,slat,s,ds)

              ff(1) = 0d0
              ff(2) = 0d0
              ff(3) = 0d0
              do  ilm1 = 1, nlm1
                l1 = ll(ilm1)
                qm1 = qmom(iv0+ilm1)
                if (ilm1 .eq. 1) qm1 = qm1 + y0*(qcorg1-z1)
                cof1 = qm1*fpi/df(2*l1+1)
                do  ilm2 = 1, nlm2
                  l2 = ll(ilm2)
                  qm2 = qmom(jv0+ilm2)
                  if (ilm2 .eq. 1) qm2 = qm2 + y0*(qcorg2-z2)
                  cof2 = qm2*fpi/df(2*l2+1)
                  xugg(ip) = xugg(ip) + cof1*cof2*s(ilm1,ilm2)
                  xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip)
     .            + s(ilm1,ilm2)*cof1*fpi/df(2*l2+1)
C         ... Forces
                  ff(1) = ff(1) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,1)
                  ff(2) = ff(2) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,2)
                  ff(3) = ff(3) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,3)
                enddo
              enddo

C     --- Additional h*h, h*g, g*h terms for foca ---
              if (lfoc1.gt.0 .or. lfoc2.gt.0) then
                call hhugbl(0,tau1,tau2,rh1,rh2,ceh1,ceh2,1,1,ndim0,ndim0,
     .          slat,wk,dwk,s0,ds0)
                xugg(ip) = xugg(ip) + cofh1*s0(1,1)*cofh2
                xhpot0(jb,ip) = xhpot0(jb,ip) + cofh1*s0(1,1)
                ff(1) = ff(1) + 0.5d0*cofh1*cofh2*ds0(1,1,1)
                ff(2) = ff(2) + 0.5d0*cofh1*cofh2*ds0(1,1,2)
                ff(3) = ff(3) + 0.5d0*cofh1*cofh2*ds0(1,1,3)

                call hgugbl(tau1,tau2,rh1,rg2,ceh1,1,nlm2,ndim,ndim,
     .          slat,s,ds)
                do  ilm2 = 1, nlm2
                  l2 = ll(ilm2)
                  qm2 = qmom(jv0+ilm2)
                  if (ilm2 .eq. 1) qm2 = qm2 + y0*(qcorg2-z2)
                  cof2 = qm2*fpi/df(2*l2+1)
                  xugg(ip) = xugg(ip) + cofh1*s(1,ilm2)*cof2
                  ff(1) = ff(1) + 0.5d0*cofh1*cof2*ds(1,ilm2,1)
                  ff(2) = ff(2) + 0.5d0*cofh1*cof2*ds(1,ilm2,2)
                  ff(3) = ff(3) + 0.5d0*cofh1*cof2*ds(1,ilm2,3)
                  xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip)
     .            + s(1,ilm2)*cofh1*fpi/df(2*l2+1)
                enddo

                call hgugbl(tau2,tau1,rh2,rg1,ceh2,1,nlm1,ndim,ndim,
     .          slat,s,ds)
                do  ilm1 = 1, nlm1
                  l1 = ll(ilm1)
                  qm1 = qmom(iv0+ilm1)
                  if (ilm1 .eq. 1) qm1 = qm1 + y0*(qcorg1-z1)
                  cof1 = qm1*fpi/df(2*l1+1)
                  xugg(ip) = xugg(ip) + cof1*s(1,ilm1)*cofh2
                  ff(1) = ff(1) - 0.5d0*cof1*cofh2*ds(1,ilm1,1)
                  ff(2) = ff(2) - 0.5d0*cof1*cofh2*ds(1,ilm1,2)
                  ff(3) = ff(3) - 0.5d0*cof1*cofh2*ds(1,ilm1,3)
                  xhpot0(jb,ip) = xhpot0(jb,ip) + cof1*s(1,ilm1)
                enddo
              endif

              if (jb .ne. ib) then
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
#if MPI | MPIK
      nvl = nvl0
#else
      nvl = iv0
#endif

C ... Assemble data from separate threads
#if MPI | MPIK
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_UGCOMP,procid,"ugcomp")
      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
#endif
      allocate(buffer(1:nvl), stat=ierr)
      call MPI_ALLREDUCE(xgpot0,buffer,nvl,
     .MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' ugcomp '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce gpot0 nvl=%i',' ',256,lgunit(3),
     .  procid,numprocs,nvl)
      endif
      call daxpy(nvl,1d0,buffer,1,gpot0,1)
      deallocate(buffer, stat=ierr)

      allocate(buffer(1:nbas), stat=ierr)
      call MPI_ALLREDUCE(xhpot0,buffer,nbas,
     .MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' ugcomp '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce hpot0 nbas=%i',' ',256,lgunit(3),
     .  procid,numprocs,nbas)
      endif
      call daxpy(nbas,1d0,buffer,1,hpot0,1)
      deallocate(buffer, stat=ierr)

      allocate(buffer(1:3*nbas), stat=ierr)
      call MPI_ALLREDUCE(xf,buffer,3*nbas,
     .MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' ugcomp '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce f 3nbas=%i',' ',256,lgunit(3),
     .  procid,numprocs,3*nbas)
      endif
      call daxpy(3*nbas,1d0,buffer,1,f,1)
      deallocate(buffer, stat=ierr)

      call MPI_ALLREDUCE(xugg,ugg,1,
     .MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit2(' ugcomp '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce ugg',' ',256,lgunit(3),
     .  procid,numprocs)
      endif
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
#endif
      deallocate(bproc, stat=ierr)
#else
      do  80  ip = 1, mp
        do  82  ib = 1, nbas
          f(1,ib) = f(1,ib) + xf(1,ib,ip)
          f(2,ib) = f(2,ib) + xf(2,ib,ip)
          f(3,ib) = f(3,ib) + xf(3,ib,ip)
          hpot0(ib) = hpot0(ib) + xhpot0(ib,ip)
   82   continue
        do  84  i = 1, nvl
          gpot0(i) = gpot0(i) + xgpot0(i,ip)
   84   continue
        ugg = ugg + xugg(ip)
   80 continue
#endif

      call tcx('ugcomp')
C      print *, ugg
C      call prmx('qmom',qmom,nvl,nvl,1)
C      call prm('gpot0',gpot0,nvl,nvl,1)
      end subroutine ugcomp


