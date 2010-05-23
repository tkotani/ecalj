#define F90 1
      subroutine smhsbl(mode,ssite,sspec,slat,vavg,q,ndimh,iprmb,
     . napw,igapw,h,s)
      use m_struc_def  !Cgetarg
#if MPE
      use m_events !Event numbers: include "events.ins"
#endif
C- Smoothed Bloch Hamiltonian (constant potential) and overlap matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 compute both hamiltonian and overlap
Ci         :  otherwise, compute overlap only.
Ci         :  In this case, vavg is not used
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
Ci     Passed to: hhibl
Ci   vavg  :constant potential (MT zero) to be added to h
Ci   q     :Bloch wave vector
Ci   ndimh :dimension of hamiltonian
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Co Outputs
Co   h     :smooth Bloch hamiltonian added to h (= s * vavg)
Co   s     :smooth Bloch overlap added to s
Cb Bugs
Cb   Use of qpgv(1,1:napw) is poor design ... rewrite
Cr Remarks
Cr  *How orbital information is extracted and deployed.
Cr   Orbital specification requires the following information:
Cr     1. orbital type, specified by the triplet (l,rsmh,eh,pz)
Cr        (pz is relevant only for local orbitals)
Cr        Information is stored in spec->orbp, unpacked by uspecb
Cr        For local orbitals spec->pz is also required.
Cr        Note that each orbital type has 2*l+1 orbitals.
Cr     3. Location of orbitals (offsets) in hamiltonian
Cr        Information is stored in iprmb (passed array)
Cr
Cr   Routines needing this information unpack the data as follows:
Cr
Cr      call    Purpose, Data and format
Cr     (Input)
Cr     -------  ---------------------------------
Cr
Cr     uspecb    Extract orbital parameters for all orbital types.
Cr    (sspec)    rsmh(l,ik),eh(l,ik), l=0,lh(ik), ik=1..nkapi.
Cr               Entries for which rsmh(l,ik)>0 have envelopes.
Cr               Entries for which ik=nkapi and pz(l)>0 are local orbitals.
Cr               Note that these possiblities can simultaneously occur.
Cr
Cr     orbl      norb : (total number of orbital types)
Cr     (iprmb)   ltab(io),ktab(io) : specifes l and index ik to
Cr               rsmh(l,ik),eh(l,ik) for orbital type io (io=1..norb).
Cr               Thus, the parameters
Cr                 l=ltab(io),ik=ktab(io),rsmh(l,ik),eh(l,ik)
Cr               completely specify envelope information for type io
Cr               offl(io): hamiltonian offset for orbital type io
Cr
Cr     gtbsl1    Serves two purposes.
Cr     (norb,    1. Blocks orbitals with common (rsmh,eh) and
Cr     ltab,     consecutive l so that mesh tabulations are more
Cr     ktab,     efficiently generated.
Cr     rsmh,eh)  2. Marks each orbital type as being :
Cr                  a. valence function with envelope;
Cr                  b. local orbital with envelope;
Cr                  c. local orbital without envelope.
Cr               blks(io) = size of contiguous block of orbitals.  Note
Cr               that if io+1, io+2, ... are contiguous to io,
Cr               blks(io+1)=blks(io+2)=...=0
Cr
Cr     A typical loop that does NOT block orbitals therefore looks like
Cr
Cr     do  io= 1, norb
Cr       l  = ltab(io)
Cr       ik = ktab(io)
Cr   C   off = orbital index in iprmb order
Cr       off = offl(io)
Cr       do  ilm1 = l**2+1, (l+1)**2
Cr         off = off+1
Cr         ...
Cr       enddo
Cr     enddo
Cr
Cr     A typical loop that DOES block orbitals therefore looks like
Cr
Cr     do  io= 1, norb
Cr     if (blks(io) .ne. 0) then
Cr       l   = ltab(io)
Cr       ik  = ktab(io)
Cr       nlm1 = l**2+1
Cr       nlm2 = nlm1 + blks1(io1)-1
Cr   C   off = orbital index in iprmb order
Cr       off = offl(io)
Cr       do  ilm1 = nlm1, nlm2
Cr         off = off+1
Cr         ...
Cr       enddo
Cr     endif
Cr     enddo
Cr
Cr  *The Bloch sum of s(q) is computed as s = sum_T s_T exp(iq.T)
Cr   where T are the set of lattice vectors.
Cr
Cr  *NB: The matrix subblock in s for sites (i,j) is computed with
Cr       the connecting vector pi-pj.  This convention is different
Cr       from the tight-binding convention, where s is computed for
Cr       pj-pi.  Thus this (i,j) matrix subblock is the hermitian
Cr       conjugate of the same subblock in the tight-binding case.
Cr       tight-binding case.  Note, however, that the entire s
Cr       would not correspond to the hermitian conjugate of the tight
Cr       binding case.
Cr
Cr  *MPI
Cr   See remarks to hsibl. Buffers for h and s are F90 ALLOCATEd as
Cr   they need to be used locally as dimensioned. A buffer is taken
Cr   from the heap for ALLREDUCE.
Cu Updates
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   14 Aug 02 Added overlap-only option
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu   02 Mar 01 Bug fix for multiple-kappa case
Cu   19 May 00 Adapted from nfp smhs_q.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI
#if MPE
      include "mpef.h"
#endif
      include "mpif.h"
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog,cmdopt
      integer i,lgunit
      character*120 strn
#endif
      integer procid,master
      integer mode,ndimh,iprmb(1),napw,igapw(3,napw)
      real(8):: q(3) , vavg
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex h(ndimh,ndimh),s(ndimh,ndimh)
C ... Local parameters
      integer nlms,kdim,n0,nkap0
      parameter (nlms=25, kdim=1, n0=10, nkap0=3)
      integer nbas,nglob,nlmto
      integer ocg,ocy,oidxcg,ojcg
      integer i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,is1,is2,nlm1,nlm2,l1,l2,ig
      integer lh1(nkap0),lh2(nkap0),nkap1,nkap2
      integer ltab1(n0*nkap0),ktab1(n0*nkap0),offl1(n0*nkap0),norb1,ik1,
     .        blks1(n0*nkap0),ntab1(n0*nkap0)
      integer ltab2(n0*nkap0),ktab2(n0*nkap0),offl2(n0*nkap0),norb2,ik2,
     .        blks2(n0*nkap0),ntab2(n0*nkap0)
      double precision e1(n0,nkap0),rsm1(n0,nkap0),p1(3),
     .                 e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx
      double complex s0(nlms,nlms,0:kdim,nkap0,nkap0)
C     Debugging
C     double complex ds(nlms,nlms,0:6)
C ... for PW part
      integer lmxax,lmxa,nlmax
      double precision qpg2,alat,plat(3,3),qlat(3,3),vol,srvol,tpiba,pi,
     .  denom,gam,fpi,ddot
#if F90
      real(8),allocatable:: yl(:),ylv(:,:),qpgv(:,:),qpg2v(:)
      complex(8),allocatable:: srm1l(:)
      complex(8):: ovl,srm1,phase,fach
      parameter (srm1=(0d0,1d0))
#endif
#if MPI
      double complex, dimension(:,:), allocatable :: hbuf, sbuf
      integer, dimension(:,:), allocatable :: index
      integer obuff
#endif
c#if MPE
cC Event numbers:
c      include "events.ins"
c#endif
C ... Heap
      integer w(1)
      common /w/ w

      integer:: iloop,iloopmx

      call tcn('smhsbl')

#if MPI
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#else
      procid = 0
      master = 0
#endif

C --- Setup ---
      nlmto = ndimh-napw
      nbas  = nglob('nbas')

      ocg=slat%ocg
      ojcg=slat%ojcg
      oidxcg=slat%oidxcg
      ocy=slat%ocy


C ... Setup below needed for PWs
      if (napw .gt. 0) then

        alat=slat%alat
        i_copy_size=size(slat%plat)
        call dcopy(i_copy_size,slat%plat,1,plat,1)
        vol=slat%vol

        pi = 4d0*datan(1d0)
        tpiba = 2d0*pi/alat
        srvol = dsqrt(vol)
        fpi = 4*pi
        call dinv33(plat,1,qlat,vol)
        vol = dabs(vol)*(alat**3)
C       Find largest lmxa ... should be made elsewhere
        lmxax = -1
        do  ib1 = 1, nbas

          is1=ssite(ib1)%spec


          lmxa=sspec(is1)%lmxa

          lmxax = max(lmxax,lmxa)
        enddo
        nlmax=(lmxax+1)**2
#if F90
        allocate(ylv(napw,nlmax),yl(nlmax),qpgv(3,napw),qpg2v(napw))
C       Note: rearrange order
        do  ig = 1, napw
          qpgv(:,ig) = tpiba * (q + matmul(qlat, igapw(:,ig)))
        enddo
        call ropyln(napw,qpgv(1,1:napw),qpgv(2,1:napw),qpgv(3,1:napw),
     .    lmxax,napw,ylv,qpg2v)
        allocate(srm1l(0:lmxax))
        srm1l(0) = 1d0
        do  l1 = 1, lmxax
          srm1l(l1) = (srm1)**l1
        enddo
#else
        call rx('PW branch requires F90')
#endif
      endif

C --- Loop over first and second site indices ---
C$DOACROSS LOCAL(i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,is1,is2,nlm1,nlm2,
C$&              lh1,lh2,nkap1,nkap2,
C$&              ltab1,ktab1,offl1,norb1,ik1,
C$&              ltab2,ktab2,offl2,norb2,ik2,l1,l2,
C$&              p1,p2,s0,
C$&              e1,rsm1,
C$&              e2,rsm2,xx)
C$&        SHARED(ssite,sspec,slat,vavg,q,ndimh,iprmb,h,s,nbas,nsp)
C$&        MP_SCHEDTYPE=RUNTIME
      if (nlmto .ne. 0) then
#if MPI
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_START_SMHSBL,procid,"smhsbl")
#endif
        allocate (hbuf(ndimh,ndimh), stat=ierr)
        allocate (sbuf(ndimh,ndimh), stat=ierr)
        call dcopy(2*ndimh*ndimh,0d0,0,hbuf,1)
        call dcopy(2*ndimh*ndimh,0d0,0,sbuf,1)
        allocate (index(0:numprocs-1,0:nbas-1), stat=ierr)
        call dstrbp(nbas,numprocs,-1,index(0,0))
cccccccccccccccccccc
c      do  iloop = 1, index(procid,0)
c        ib1 = index(procid,iloop)
c        print *,'iloop, procid ib1=',iloop,procid,ib1
c      enddo
c      call MPI_abort(MPI_COMM_WORLD,ierr)
cccccccccccccccccccc
ccccccccccccccccccccccc
c      print *,'aaaaaaaaaaaa111111111'
c      call uspecb(0,2,sspec,1,1,lh1,rsm1,e1,nkap1)
c      print *,'aaaaaaaaaaaa222222222222'
c      call MPI_abort(MPI_COMM_WORLD,ierr)
cccccccccccccccccccccc
c      print *,'procid=',procid,index(procid,0)
ccccccccccccccccccccc
        iloopmx= index(procid,0)
#else
        iloopmx= nbas
#endif

        do iloop = 1,iloopmx
#if MPI
          ib1 = index(procid,iloop)
          if (mlog) then
            call gettime(datim)
            call awrit4(' smhsbl '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting atom %i of %i',' ',256,lgunit(3),
     .        procid,numprocs,ib1,index(procid,0))
          endif
#else
          ib1=iloop
#endif


          is1=ssite(ib1)%spec
          i_copy_size=size(ssite(ib1)%pos)
          call dcopy(i_copy_size,ssite(ib1)%pos,1,p1,1)

          call uspecb(0,2,sspec,is1,is1,lh1,rsm1,e1,nkap1)
C       Row info telling smhsbl where to poke s0 made by hhibl
          call orbl(ib1,0,nlmto,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
          call gtbsl1(8+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)

          do  ib2 = ib1, nbas

            is2=ssite(ib2)%spec
            i_copy_size=size(ssite(ib2)%pos)
            call dcopy(i_copy_size,ssite(ib2)%pos,1,p2,1)

            call uspecb(0,2,sspec,is2,is2,lh2,rsm2,e2,nkap2)
C         Column info telling smhsbl where to poke s0 made by hhibl
            call orbl(ib2,0,nlmto,iprmb,norb2,ltab2,ktab2,xx,offl2,xx)
            call gtbsl1(8+16,norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)

C     ... M.E. <1> and <T> between all envelopes connecting ib1 and ib2
            do  i1 = 1, nkap1
              do  i2 = 1, nkap2
                nlm1 = (lh1(i1)+1)**2
                nlm2 = (lh2(i2)+1)**2
                if (nlm1 .gt. nlms .or. nlm2 .gt. nlms)
     .          call rx('smhsbl: increase nlms')
                call hhibl(11,p1,p2,q,rsm1(1,i1),rsm2(1,i2),e1(1,i1),
     .          e2(1,i2),nlm1,nlm2,1,nlms,nlms,w(ocg),w(oidxcg),w(ojcg),
     .          w(ocy),slat,s0(1,1,0,i1,i2))

C              print *, '!!'; k0=1
C              call hhigbl(11,p1,p2,q,rsm1(1,i1),rsm2(1,i2),e1(1,i1),
C     .          e2(1,i2),nlm1,nlm2,1,nlms,nlms,k0,w(ocg),w(oidxcg),
C     .          w(ojcg),w(ocy),slat,s0(1,1,0,i1,i2),ds(1,1,0))
C              s0(:,:,0,i1,i2) = ds(:,:,0)

              enddo
            enddo

C     ... Loop over orbital indices, poke block of integrals into s,h
            do  io2 = 1, norb2
              if (blks2(io2) .ne. 0) then
C           l2,ik2 = l and kaph indices, needed to locate block in s0
                l2  = ltab2(io2)
                ik2 = ktab2(io2)
C           i2 = orbital index in iprmb order
                i2 = offl2(io2)

                do  ilm2 = l2**2+1, (l2+1)**2
                  i2 = i2+1

                  if (mode .eq. 0) then
                    do  io1 = 1, norb1
                      if (blks1(io1) .ne. 0) then
C               l1,ik1 = l and kaph indices, needed to locate block in s0
                        l1  = ltab1(io1)
                        ik1 = ktab1(io1)
C               i1 = orbital index in iprmb order
                        i1 = offl1(io1)
                        do  ilm1 = l1**2+1, (l1+1)**2
                          i1 = i1+1
#if MPI
                          sbuf(i1,i2) = sbuf(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
                          hbuf(i1,i2) = hbuf(i1,i2) - s0(ilm1,ilm2,1,ik1,ik2)
     .              + vavg*s0(ilm1,ilm2,0,ik1,ik2)
#else
                          s(i1,i2) = s(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
                          h(i1,i2) = h(i1,i2) - s0(ilm1,ilm2,1,ik1,ik2)
     .              + vavg*s0(ilm1,ilm2,0,ik1,ik2)
#endif
                        enddo
                      endif
                    enddo
                  else
                    do  io1 = 1, norb1
                      if (blks1(io1) .ne. 0) then
C               l1,ik1 = l and kaph indices, needed to locate block in s0
                        l1  = ltab1(io1)
                        ik1 = ktab1(io1)
C               i1 = orbital index in iprmb order
                        i1 = offl1(io1)
                        do  ilm1 = l1**2+1, (l1+1)**2
                          i1 = i1+1
#if MPI
                          sbuf(i1,i2) = sbuf(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
#else
                          s(i1,i2) = s(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
#endif
                        enddo
                      endif
                    enddo
                  endif

                enddo
              endif
            enddo
C ... end loop over ib2
          enddo
C   ... Hsm (i1) \times PW (i2)  Takao. Similar logic in fsmbl
#if F90
          do  ig = 1, napw
            i2 = ig + nlmto
            qpg2 = qpg2v(ig)
            phase = exp(srm1*alat*ddot(3,qpgv(1,ig),1,p1,1))
C         phase = exp( srm1 * sum(qpgv(:,ig)*p1)*alat  )
            do  io1 = 1, norb1
              if (blks1(io1) == 0) cycle
              l1  = ltab1(io1)
              ik1 = ktab1(io1)
              i1  = offl1(io1)
              denom = e1(l1+1,ik1) - qpg2
              gam   = 1d0/4d0*rsm1(l1+1,ik1)**2
C           Note: fach depends on l.
              fach  = -fpi/denom * phase * srm1l(l1) * exp(gam*denom)
              do  ilm1 = l1**2+1, (l1+1)**2
                i1 = i1+1
C             <Hsm^bloch | exp(i q+G r)/srvol > and
C             <Hsm^bloch | nabla + vavg |  exp(i q+G r)/srvol > in a cell.
                ovl = fach * ylv(ig,ilm1)/srvol ! JMP Eq.(9.4)
C             gradient of ovl
C             ovl = srm1*qpgv(1,ig) * ovl
#if MPI
                sbuf(i1,i2) = sbuf(i1,i2) + ovl
                hbuf(i1,i2) = hbuf(i1,i2) + qpg2*ovl + vavg*ovl
#else
                s(i1,i2) = s(i1,i2) + ovl
                h(i1,i2) = h(i1,i2) + qpg2*ovl + vavg*ovl
#endif
              enddo
            enddo
          enddo
#endif

C ... end loop over ib1
        enddo

#if MPI
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_END_SMHSBL,procid,"smhsbl")
        ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
        ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
#endif
        if (mode .eq. 0) then
          call defcc(obuff, ndimh*ndimh)

          call MPI_ALLREDUCE(hbuf,w(obuff),2*ndimh*ndimh,
     .     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          if (mlog) then
            call gettime(datim)
            call awrit3(' smhsbl '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce h ndimh=%i',' ',256,lgunit(3),
     .        procid,numprocs,ndimh)
          endif
          call daxpy(2*ndimh*ndimh,1d0,w(obuff),1,h,1)
        endif
        call MPI_ALLREDUCE(sbuf,w(obuff),2*ndimh*ndimh,
     .     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit3(' smhsbl '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce s ndimh=%i',' ',256,lgunit(3),
     .        procid,numprocs,ndimh)
        endif
        call daxpy(2*ndimh*ndimh,1d0,w(obuff),1,s,1)
        call rlse(obuff)
        deallocate(index, stat=ierr)
        deallocate(hbuf, stat=ierr)
        deallocate(sbuf, stat=ierr)
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
#endif
#endif
      endif



C ... PW x PW part (diagonal matrix)
      do  ig = 1, napw
        i2 = ig + nlmto
        s(i2,i2) = s(i2,i2) + 1
        h(i2,i2) = h(i2,i2) + qpg2v(ig) + vavg
      enddo

#if F90
      if (napw .gt. 0) then
        deallocate(yl,ylv,qpgv,qpg2v,srm1l)
      endif
#endif

C ... Occupy second half of matrix
      if (mode .eq. 0) then
        do  i1 = 1, ndimh
          do  i2 = i1, ndimh
            h(i2,i1) = dconjg(h(i1,i2))
            s(i2,i1) = dconjg(s(i1,i2))
          enddo
        enddo
      else
        do  i1 = 1, ndimh
          do  i2 = i1, ndimh
            s(i2,i1) = dconjg(s(i1,i2))
          enddo
        enddo
      endif

      call tcx('smhsbl')

C      call zprm('h-warp',2,h,ndimh,ndimh,ndimh)
C      call zprm('s-warp',2,s,ndimh,ndimh,ndimh)
      end subroutine smhsbl


