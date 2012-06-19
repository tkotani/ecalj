      subroutine fsmbl(nbas,ssite,sspec,slat,vavg,q,ndimh,nlmto,iprmb,
     .numq,nevec,evl,evec,ewgt,f)

      use m_struc_def  !Cgetarg
      use m_globalvariables
#if MPE
      use m_events !Event numbers: include "events.ins"
#endif
C- Force from smoothed hamiltonian (constant potential) and overlap
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: uspecb
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy
Ci     Stored:    *
Ci     Passed to: hhigbl
Ci   vavg  :constant potential (MT zero) to be added to h
Ci   q     :Bloch wave vector
Ci   ndimh :dimension of hamiltonian
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   numq  :number of trial Fermi levels
Ci   nevec :number of occupied eigenvectors
Ci   evl   :eigenvalues
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Co Outputs
Co   f
Cr Remarks
Cu Updates
Cu   05 Jul 08 Decouple ndimh from nlmto, for PW basis
Cu   10 Apr 02 Redimensioned eh,rsmh to accommodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   10 Sep 01 Extended to local orbitals.
Cu   23 May 00 Adapted from nfp fsm_q.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI
#if MPE
      include "mpef.h"
#endif
      include "mpif.h"
      integer pid, procid, master, numprocs, ierr,
     .status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen,i,lgunit
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog,cmdopt
      character*120 strn
#endif
c#if MPE
cC Event numbers:
c      include "events.ins"
c#endif
      integer nbas,ndimh,nlmto,nevec,numq,iprmb(nlmto)
      real(8):: q(3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double precision evl(ndimh),f(3,nbas,numq),ewgt(numq,nevec),vavg
      double complex evec(ndimh,ndimh)
C ... Local parameters
      integer nlms,k0,n0,nkap0
      parameter (nlms=25, k0=1, n0=10, nkap0=3)
      integer:: i1 , i2 , ib1 , ib2 , ilm1 , ilm2 , io1 , io2 , iq 
     .    , is1 , is2 , l1 , l2 , ik1 , ik2 , ivec , m , nglob , nlm1 , 
     .    nlm2
ckino Dec.9.2011:            integer,pointer :: iv_p_ojcg(:) =>NULL()

ckino Dec.9.2011:           integer,pointer :: iv_p_oidxcg(:) =>NULL()

ckino Dec.9.2011:          real(8),pointer :: rv_p_ocy(:) =>NULL()

ckino Dec.9.2011:         real(8),pointer :: rv_p_ocg(:) =>NULL()

      integer lh1(nkap0),lh2(nkap0),nkap1,nkap2,nlm21,nlm22,nlm11,nlm12
      integer norb1,ltab1(n0*nkap0),ktab1(n0*nkap0),offl1(n0*nkap0),
     .blks1(n0*nkap0),ntab1(n0*nkap0)
      integer norb2,ltab2(n0*nkap0),ktab2(n0*nkap0),offl2(n0*nkap0),
     .blks2(n0*nkap0),ntab2(n0*nkap0)
      double precision e1(n0,nkap0),rsm1(n0,nkap0),p1(3),wt,
     .e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx(n0)
      double complex  s(nlms,nlms,0:k0,nkap0,nkap0)
      double complex ds(nlms,nlms,0:k0,3,nkap0,nkap0),ccc,sum
#if MPI
      double precision, dimension(:,:,:), allocatable :: xf
      double precision, dimension(:),     allocatable :: buffer
      integer, dimension(:,:), allocatable :: index
      integer iloop,ib
#endif
C ... Heap
      integer:: iloopend

C --- Setup ---
      if (nevec .le. 0) return
      call tcn ('fsmbl')

#if MPI
      allocate(xf(1:3,1:nbas,1:numq), stat=ierr)
      call dcopy(3*nbas*numq,0d0,0,xf,1)
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
      if (mlog) then
        do  pid = 0, numprocs-1
          call MPI_BCAST(shortname(pid),10,MPI_CHARACTER,pid,
     .    MPI_COMM_WORLD,ierr)
          call MPI_BCAST(namelen(pid),1,MPI_INTEGER,pid,
     .    MPI_COMM_WORLD,ierr)
        enddo
      endif
#endif

Changenglob      nbas  = nglob('nbas')
      nbas  = globalvariables%nbas

ckino Dec.9.2011:         rv_p_ocg => slat%rv_p_ocg

ckino Dec.9.2011:         iv_p_ojcg => slat%iv_p_ojcg

ckino Dec.9.2011:         iv_p_oidxcg => slat%iv_p_oidxcg

ckino Dec.9.2011:         rv_p_ocy => slat%rv_p_ocy



C --- Loop over first and second site indices ---
#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_START_FSMBL,procid,"fsmbl")
#endif
      allocate (index(0:numprocs-1,0:nbas-1), stat=ierr)
      call dstrbp(nbas,numprocs,-1,index(0,0))
c      do  iloop = 1, index(procid,0)
      iloopend=index(procid,0)
#else
      iloopend=nbas
#endif


      do  iloop = 1, iloopend
#if MPI
        ib1 = index(procid,iloop)
        if (mlog) then
          call gettime(datim)
          call awrit4(' fsmbl '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' starting atom %i of %i',' ',256,lgunit(3),
     .    procid,numprocs,ib1,index(procid,0))
        endif
#else
        ib1=iloop
#endif

        is1=ssite(ib1)%spec
        i_copy_size=size(ssite(ib1)%pos)
        call dcopy(i_copy_size,ssite(ib1)%pos,1,p1,1)

        call uspecb(0,2,sspec,is1,is1,lh1,rsm1,e1,nkap1)
C       Row info telling fsmbl where to poke s0 made by hhibl
        call orbl(ib1,0,nlmto,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
        call gtbsl1(4+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)

        do  ib2 = ib1+1, nbas

          is2=ssite(ib2)%spec
          i_copy_size=size(ssite(ib2)%pos)
          call dcopy(i_copy_size,ssite(ib2)%pos,1,p2,1)

          call uspecb(0,2,sspec,is2,is2,lh2,rsm2,e2,nkap2)
C         Column info telling fsmbl where to poke s0 made by hhibl
          call orbl(ib2,0,nlmto,iprmb,norb2,ltab2,ktab2,xx,offl2,xx)
          call gtbsl1(4+16,norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)

C     ... M.E. <1> and <KE> between all envelopes connecting ib1 and ib2
          do  i1 = 1, nkap1
            do  i2 = 1, nkap2
              nlm1 = (lh1(i1)+1)**2
              nlm2 = (lh2(i2)+1)**2
              if (nlm1 .gt. nlms .or. nlm2 .gt. nlms)
     .        call rx('fsmbl: increase nlms')
              call hhigbl ( 11 , p1 , p2 , q , rsm1 ( 1 , i1 ) , rsm2 ( 1 , 
     .            i2 ) , e1 ( 1 , i1 ) , e2 ( 1 , i2 ) , nlm1 , nlm2 , 1 , nlms
ckino Dec.9.2011:       .            , nlms , k0 , rv_p_ocg , iv_p_oidxcg , iv_p_ojcg , rv_p_ocy ,
ckino Dec.9.2011:       .            , nlms , k0 , rv_p_ocg , iv_p_oidxcg , slat%iv_p_ojcg , rv_p_ocy ,
ckino Dec.9.2011:       .            , nlms , k0 , rv_p_ocg , slat%iv_p_oidxcg , slat%iv_p_ojcg , rv_p_ocy ,
ckino Dec.9.2011:       .            , nlms , k0 , rv_p_ocg , slat%iv_p_oidxcg , slat%iv_p_ojcg , slat%rv_p_ocy ,
ckino Dec.27.2011:       .            , nlms , k0 , slat%rv_p_ocg , slat%iv_p_oidxcg , slat%iv_p_ojcg , slat%rv_p_ocy ,
ckino Dec.27.2011:       .            , nlms , k0 , slat%rv_p_ocg , slat%iv_a_oidxcg , slat%iv_p_ojcg , slat%rv_p_ocy ,
ckino Dec.27.2011:       .            , nlms , k0 , slat%rv_a_ocg , slat%iv_a_oidxcg , slat%iv_p_ojcg , slat%rv_p_ocy ,
ckino Dec.28.2011:       .            , nlms , k0 , slat%rv_a_ocg , slat%iv_a_oidxcg , slat%iv_p_ojcg , slat%rv_a_ocy , 
     .            , nlms , k0 , slat%rv_a_ocg , slat%iv_a_oidxcg , slat%iv_a_ojcg , slat%rv_a_ocy , 
     .            slat , s ( 1 , 1 , 0 , i1 , i2 ) , ds ( 1 , 1 , 0 , 1 , i1 , 
     .            i2 ) )




            enddo
          enddo

C     ... Loop over pairs of orbital groups, multiply bloc of gradients
          do  io2 = 1, norb2
            if (blks2(io2) .ne. 0) then
C           l2,ik2 = l and kaph indices, needed to locate block in s0
              l2  = ltab2(io2)
              ik2 = ktab2(io2)
              nlm21 = l2**2+1
              nlm22 = nlm21 + blks2(io2)-1
              do  io1 = 1, norb1
                if (blks1(io1) .ne. 0) then
C             l1,ik1 = l and kaph indices, needed to locate block in s0
                  l1  = ltab1(io1)
                  ik1 = ktab1(io1)
                  nlm11 = l1**2+1
                  nlm12 = nlm11 + blks1(io1)-1

C         ... Loop over occupied eigenstates and x,y,z
                  do  ivec = 1, nevec
                    do  m = 1, 3

C           ... Loop over orbital pairs within the groups
                      sum = 0d0
C               i2 = hamiltonian offset
                      i2 = offl2(io2)
                      do  ilm2 = nlm21, nlm22
                        i2 = i2+1
C                 i1 = orbital index in iprmb order
                        i1 = offl1(io1)
                        do  ilm1 = nlm11, nlm12
                          i1 = i1+1
                          ccc = vavg*ds(ilm1,ilm2,0,m,ik1,ik2)
     .                    -      ds(ilm1,ilm2,1,m,ik1,ik2)
     .                    - evl(ivec)*ds(ilm1,ilm2,0,m,ik1,ik2)
                          sum = sum + dconjg(evec(i1,ivec))*ccc*evec(i2,ivec)
                        enddo
                      enddo
                      do  iq = 1, numq
                        wt = ewgt(iq,ivec)
#if MPI
                        xf(m,ib1,iq) = xf(m,ib1,iq) - 2*wt*sum
                        xf(m,ib2,iq) = xf(m,ib2,iq) + 2*wt*sum
#else
                        f(m,ib1,iq) = f(m,ib1,iq) - 2*wt*sum
                        f(m,ib2,iq) = f(m,ib2,iq) + 2*wt*sum
#endif
                      enddo
                    enddo
                  enddo

                endif
              enddo
            endif
          enddo

        enddo
      enddo

#if MPI
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_FSMBL,procid,"fsmbl")
      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
#endif
      allocate(buffer(1:3*nbas*numq), stat=ierr)
      call MPI_ALLREDUCE(xf,buffer,3*nbas*numq,
     .MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' fsmbl '//datim//' Process %i of %i on '
     .  //shortname(procid)(1:namelen(procid))//
     .  ' allreduce f 3*nbas*numq=%i',' ',256,lgunit(3),
     .  procid,numprocs,3*nbas*numq)
      endif
      call daxpy(3*nbas*numq,1d0,buffer,1,f,1)
      deallocate(buffer, stat=ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
#endif
      deallocate(index, stat=ierr)
#endif

C      write(6,*) '---- END OF FSMBL ---'
C      do  ib = 1, nbas
C        do  iq = 1, numq
C          write(stdo,220) ib,iq,(f(m,ib,iq),m=1,3)
C  220     format(2i5,3f12.6)
C        enddo
C      enddo

      call tcx ('fsmbl')

      end subroutine fsmbl


