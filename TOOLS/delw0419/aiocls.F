Cgetarg...info...           structure ['aiocls', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['aiocls', 'pot', 'spot']
Cgetarg...info...           structure ['aiocls', 'array', 'sarray']
Cgetarg...info...           structure ['aiocls', 'spec', 'sspec']
Cgetarg...info...           structure ['aiocls', 'lat', 'slat']
Cgetarg...info...           use_to_add ['aiocls', 'm_struc_def']
#define ONE_ATOM_UNIT 1
      subroutine aiocls(lio,mode,sarray,sctrl,sham,spot,sspec,slat,ic1,
     .  ic2)
       
       use m_struc_def  !Cgetarg

C- File I/O atomic data for classes ic1..ic2
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio    F for read, T for write
Ci   mode   ones digit
Ci          1 use data from first class corresponding to same species,
Ci            if data from own class is missing
Ci          2 like 1, but attempt to read data from disk anyway
Ci          3 make no attempt to read from disk but copy info from
Ci            first class corresponding to same species, if data missing
Ci          4 Add 4 if to use default P,Q when not otherwise supplied
Ci          10s digit
Ci          1 assemble background rho
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp oclabl ohave oics onrc
Ci     Stored:
Ci     Passed to:
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nl nspin zbak nbas nclass
Ci     Stored:    zbak
Ci     Passed to: lgors
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci          (not used now)
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opp osop ogrrme opnu oqnu orhrmx ovrmax ovintr opmpol
Ci                oves
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idmod p q
Ci     Stored:
Ci     Passed to:
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: vol
Ci     Stored:
Ci     Passed to:
Ci   ic1,ic2: range of classes to read data
Cr Remarks
Cr   Right now, read always takes data from file if available
Cu Updates
Cu   09 Nov 07 Corrected sign of default moment (paioc2)
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   26 Apr 03 Added MPI calls
Cu   07 Feb 03 adjusted for redimensioned sop
Cu   30 May 02 Assign better default P
Cu   28 Apr 98 I/O of radial matrix elements of grad
Cu   28 Sep 00 Added setting default P,Q
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lio
      integer mode,ic1,ic2
Cgetarg       double precision sarray(1),sctrl(1),sspec(1),slat(1),sham(1),
Cgetarg      .  spot(1)
       type(s_array)::sarray
       type(s_ctrl)::sctrl
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_ham)::sham
       type(s_pot)::spot

C heap
      integer w(1)
      common /w/ w
C Local variables
      character*8 clabl,alabel, outs1*20, outs2*20, outs3*20, outs*80
      logical sw,aiomom,aiopar,aiopot,aiova,lpot,lcor,scat,lgors,
     .  aiogen,aiosop,aiorme,aiocor,aiomp,lrell,lgen
      logical havepq,havepp,haveso,haveop,haveva,havemp,
     .        readpq,readpp,readso,readop,readva,readmp
      integer opnu,oqnu,orhrmx,ovrmax,opp,opprel,oves,ova,omp,osop,olpp,
     .  olppr,olsop,olgrme,olva,olmp,opot,ocor,oics,ohave,oclabl,onrc,
     .  ogrrme,n0,nglob
      parameter(n0=10)
      integer ic,is,lmx,k,nl,nsp,ifi,jfi,fopn,lmxx,nspx,nrx,nr,
     .  idmod(n0),isw,bitand,i2,nclasp,iclbsj,icmap,jc,nbas,nclass,
     .  nclspp,iprint,lgunit,mode0,mode00,lrel
      integer mpipid,procid
      double precision rhrmx,vrmax(2),ves,z,rmxx,ax,qc,dq,vrmxx(2),
     .  sumec,sumtc,sumev,thrpv,ekin,utot,rhoeps,etot,a,rmax,dval,
     .  zbak(2),dgets,pdf(n0,2),qdf(n0,2),pnuloc(100),qnuloc(100)

Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg       call upack('ctrl nl nspin nclass',sctrl,nl,nsp,nclass,0,0)
       
       nl=sctrl%nl
       nsp=sctrl%nspin
       nclass=sctrl%nclass

Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg...info...           integer :: opprel , integer(8) :: spot%opprel 1
Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
Cgetarg...info...           integer :: ogrrme , integer(8) :: spot%ogrrme 1
Cgetarg       call upack('pot opp opprel osop ogrrme',spot,opp,opprel,osop,
Cgetarg      .  ogrrme,0)
       
       opp=spot%opp
       opprel=spot%opprel
       osop=spot%osop
       ogrrme=spot%ogrrme

Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg       call upack('pot opnu oqnu orhrmx ovrmax',spot,opnu,oqnu,orhrmx,
Cgetarg      .  ovrmax,0)
       
       opnu=spot%opnu
       oqnu=spot%oqnu
       orhrmx=spot%orhrmx
       ovrmax=spot%ovrmax

Cgetarg...info...           integer :: ova , integer(8) :: spot%ovintr 1
Cgetarg...info...           integer :: omp , integer(8) :: spot%opmpol 1
Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg       call upack('pot ovintr opmpol oves',spot,ova,omp,oves,0,0)
       
       ova=spot%ovintr
       omp=spot%opmpol
       oves=spot%oves

Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: ohave , integer(8) :: sarray%ohave 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg       call upack('array nclasp oclabl ohave oics',sarray,nclasp,oclabl,
Cgetarg      .  ohave,oics,0)
       
       nclasp=sarray%nclasp
       oclabl=sarray%oclabl
       ohave=sarray%ohave
       oics=sarray%oics

      lrel = nglob('lrel')

      mode0 = mod(mode,10)
      mode00 = mod(mode0,4)

C ... MPI: only master does sphere program
      procid = mpipid(1)
      if (procid .eq. 0) then

      olppr = 1
      call defdr(olpp,100)
      call defdr(olppr,3000)
      call defdr(olsop,500)
      call defdr(olgrme,500)
      call defdr(olva,100)
      call defdr(olmp,3000)
      call defdr(opot,3000)
      call defdr(ocor,3000)
      i2 = ic2
      if (i2 .eq. 0) i2 = nclasp
c     call awrit2('%n:1i',' ',100,6,nclasp,w(oics))
      do  10  ic = ic1, i2
        is = w(oics+ic-1)
        icmap = iclbsj(is,w(oics),-nclasp,1)
        if (icmap .eq. ic .or. mode00 .eq. 0) icmap = 0
Cgetarg...info...           integer :: lmx , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idmod(n0) , integer(8) :: sspec%idmod 10
Cgetarg...info...           real(8) :: pdf(n0,2) , real(8) :: sspec%p 20
Cgetarg...info...           real(8) :: qdf(n0,2) , real(8) :: sspec%q 20
Cgetarg         call upack('spec lmxa idmod p q',sspec,is,lmx,idmod,pdf,qdf)
         
         lmx=sspec(is)%lmxa
         i_copy_size=size(sspec(is)%idmod) 
         call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 
         i_copy_size=size(sspec(is)%p) 
         call dcopy(i_copy_size,sspec(is)%p,1,pdf,1) 
         i_copy_size=size(sspec(is)%q) 
         call dcopy(i_copy_size,sspec(is)%q,1,qdf,1) 


        call dpzero(pnuloc,100)
        call dpzero(qnuloc,100)

C   --- Open the atom file ---
        outs1 = ' '
        outs2 = ' '
        outs3 = ' '
        call r8tos8(dval(w(oclabl),ic),clabl)
        if (mode00 .ne. 3) then
#if ONE_ATOM_UNIT
          ifi = fopn(clabl)
#else
         ifi = fopna(clabl,30+ic,0)
#endif
        endif

C  --- Copy what is passed through to holding arrays ---
        k = nl*nsp
        call pvaioc(w(ohave),1,ic,icmap,havepq,jc)
        readpq = .not. havepq .or. mode00 .eq. 2 .and. jc .ne. ic
        readpq = readpq .and. mode00 .lt. 3
        if (havepq) then
          call dpscop(w(opnu),pnuloc,k,1+(jc-1)*k,1,1d0)
          call dpscop(w(oqnu),qnuloc,3*k,1+(jc-1)*3*k,1,1d0)
          call dpscop(w(orhrmx),rhrmx,1,jc,1,1d0)
          call dpscop(w(ovrmax),vrmax,2,2*jc-1,1,1d0)
          if (ic.ne.jc) call awrit0('%a pq,',outs1,len(outs1),0)
        endif
        k = 6*nl*nsp
        call pvaioc(w(ohave),2,ic,icmap,havepp,jc)
        readpp = .not. havepp .or. mode00 .eq. 2 .and. jc .ne. ic
        readpp = readpp .and. mode00 .lt. 3
        if (havepp) then
          call dpscop(w(opp),w(olpp),k,1+(jc-1)*k,1,1d0)
          if (lrel .eq. 2) then
            k = 4*nl*2*nl*2*2
            call dpscop(w(opprel),w(olppr),k,1+(jc-1)*k,1,1d0)
          endif
          call dpscop(w(oves),ves,1,jc,1,1d0)
          if (ic.ne.jc) call awrit0('%a pp,',outs1,len(outs1),0)
        endif
        k = (nl*nsp)**2
        call pvaioc(w(ohave),8,ic,icmap,haveva,jc)
        haveva = haveva .and. ova .gt. 0
        readva = .not. haveva .or. mode00 .eq. 2 .and. jc .ne. ic
        readva = readva .and. ova .gt. 0 .and. mode00 .lt. 3
        if (haveva) then
          call dpscop(w(ova),w(olva),k,1+(jc-1)*k,1,1d0)
          if (ic.ne.jc) call awrit0('%a va,',outs1,len(outs1),0)
        endif
        k = nl**2*(2*nl-1)*3*nsp
        call pvaioc(w(ohave),16,ic,icmap,havemp,jc)
        havemp = havemp .and. omp .gt. 0
        readmp = .not. havemp .or. mode00 .eq. 2 .and. jc .ne. ic
        readmp = readmp .and. omp .gt. 0 .and. mode00 .lt. 3
        if (havemp) then
          call dpscop(w(omp),w(olmp),k,1+(jc-1)*k,1,1d0)
          if (ic.ne.jc) call awrit0('%a mp,',outs1,len(outs1),0)
        endif
        k = nl*nsp*nsp*9
        call pvaioc(w(ohave),4,ic,icmap,haveso,jc)
        haveso = haveso .and. osop.gt. 0
        readso = .not. haveso .or. mode00 .eq. 2 .and. jc .ne. ic
        readso = readso .and. osop .gt. 0 .and. mode00 .lt. 3
        if (haveso) then
          call dpscop(w(osop),w(olsop),k,1+(jc-1)*k,1,1d0)
          if (ic.ne.jc) call awrit0('%a sop,',outs1,len(outs1),0)
        endif
        k = 16*nl*nsp
        call pvaioc(w(ohave),32,ic,icmap,haveop,jc)
        haveop = haveop .and. ogrrme .gt. 0
        readop = .not. haveop .or. mode00 .eq. 2 .and. jc .ne. ic
        readop = readop .and. ogrrme .gt. 0 .and. mode00 .lt. 3
        if (haveop) then
          call dpscop(w(ogrrme),w(olgrme),k,1+(jc-1)*k,1,1d0)
          if (ic.ne.jc) call awrit0('%a opp,',outs1,len(outs1),0)
        endif

C   --- File WRITE ---
        if (lio) then
          lgen = .false.
          lpot = .false.
          lcor = .false.
C     ... Pick up GEN and POT, if available, to save again
          if (scat(iabs(ifi),'GEN:',':',.true.)) then
            lgen = aiogen(alabel,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,dq,
     .        vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          endif
          if (scat(iabs(ifi),'POT:',':',.true.)) then
            read(ifi,102) nr,nsp,a,rmax
  102       format(2i5,2f12.5)
            lpot = aiopot(nr,nsp,a,rmax,-99d0,w(opot),ifi)
          endif
          lcor = aiocor(nr,nsp,a,rmxx,w(ocor),sumec,sumtc,ifi)

          rewind ifi
          jfi = -ifi
          if (lgen) sw = aiogen(clabl,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,
     .      dq,vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,jfi)
          if (havepq) sw = aiomom(clabl,pnuloc,qnuloc,idmod,
     .      nl,lmx,nsp,rhrmx,vrmax,jfi)
          if (havepp) sw = aiopar(clabl,lrel,w(olpp),w(olppr),ves,nl,
     .      lmx,nsp,jfi)
          if (haveva) sw = aiova(clabl,w(olva),nl,lmx,nsp,jfi)
          if (havemp) sw = aiomp(clabl,w(olmp),nl,2*nl-2,nsp,jfi)
          if (haveso) sw = aiosop(clabl,w(olsop),nl,lmx,nsp,jfi)
          if (haveop) sw = aiorme(clabl,w(olgrme),nl,nsp,jfi)
          if (lpot)   sw = aiopot(nr,nsp,a,rmax,-99d0,w(opot),jfi)
          if (lcor) lcor = aiocor(nr,nsp,a,rmxx,w(ocor),sumec,sumtc,jfi)

C   --- File READ ---
        else

C     ... Copy whatever is available on disk to holding arrays
          if (readpq .or. mode0 .ge. 4) then
            if (readpq) rewind ifi
            if (readpq) readpq = aiomom(clabl,pnuloc,qnuloc,idmod,
     .        nl,lmx,nsp,rhrmx,vrmax,ifi)
C           Couldn't read from atom file ; take default values
            if (readpq) call awrit0('%a pq,',outs2,len(outs2),0)
            if (mode0 .ge. 4 .and. .not. (readpq .or. havepq)) then
C             call dmcpy(pdf,n0,1,pnuloc,nl,1,nl,nsp)
              call paioc2(nsp,nl,n0,pdf,qdf,pnuloc,qnuloc)
              call awrit0('%a pq,',outs3,len(outs2),0)
              call dvset(vrmax,1,2,-.7d0)
              rhrmx = .1d0
              readpq = .true.
            endif
          endif
          if (readpp) then
            readpp = aiopar(clabl,lrel,w(olpp),w(olppr),ves,nl,lmx,nsp,
     .        ifi)
            if (readpp) call awrit0('%a pp,',outs2,len(outs2),0)
          endif
          if (readso .and. osop .gt. 0) then
            readso = aiosop(clabl,w(olsop),nl,lmx,nsp,ifi)
            if (readso) call awrit0('%a so,',outs2,len(outs2),0)
          endif
          if (readop .and. ogrrme .gt. 0) then
            readop = aiorme(clabl,w(olgrme),nl,nsp,ifi)
            if (readop) call awrit0('%a op,',outs2,len(outs2),0)
          endif
          if (readva .and. ova .gt. 0) then
            readva = aiova(clabl,w(olva),nl,lmx,nsp,ifi)
            if (readva) call awrit0('%a va,',outs2,len(outs2),0)
          endif
          if (readmp .and. omp .gt. 0) then
            readmp = aiomp(clabl,w(olmp),nl,2*nl-2,nsp,ifi)
            if (readmp) call awrit0('%a mp,',outs2,len(outs2),0)
          endif

C     ... Update what parameters are available
          w(ohave+ic-1) = isw(haveop.or.readop)*32+
     .                    isw(havemp.or.readmp)*16+
     .                    isw(haveva.or.readva)*8 +
     .                    isw(haveso.or.readso)*4 +
     .                    isw(havepp.or.readpp)*2 +
     .                    isw(havepq.or.readpq)*1 +
     .                    w(ohave+ic-1) - bitand(w(ohave+ic-1),63)

          k = nl*nsp
          if (havepq .or. readpq) then
            call dpscop(pnuloc,w(opnu),k,1,1+(ic-1)*k,1d0)
            call dpscop(qnuloc,w(oqnu),3*k,1,1+(ic-1)*3*k,1d0)
            call dpscop(rhrmx,w(orhrmx),1,1,ic,1d0)
            call dpscop(vrmax,w(ovrmax),2,1,2*ic-1,1d0)
          endif
          k = 6*nl*nsp
          if (havepp .or. readpp) then
            call dpscop(w(olpp),w(opp),k,1,1+(ic-1)*k,1d0)
            if (lrel .eq. 2) then
              k = 4*nl*2*nl*2*2
              call dpscop(w(olppr),w(opprel),k,1,1+(ic-1)*k,1d0)
            endif
            call dpscop(ves,w(oves),1,1,ic,1d0)
          endif
          k = (nl*nsp)**2
          if ((haveva .or. readva) .and. ova .gt. 0)
     .      call dpscop(w(olva),w(ova),k,1,1+(ic-1)*k,1d0)
          k = nl**2*(2*nl-1)*3*nsp
          if ((havemp .or. readmp) .and. omp .gt. 0)
     .      call dpscop(w(olmp),w(omp),k,1,1+(ic-1)*k,1d0)
          k = nl*nsp*nsp*9
          if ((haveso .or. readso) .and. osop .gt. 0)
     .      call dpscop(w(olsop),w(osop),k,1,1+(ic-1)*k,1d0)
          k = 16*nl*nsp
          if ((haveop .or. readop) .and. ogrrme .gt. 0)
     .      call dpscop(w(olgrme),w(ogrrme),k,1,1+(ic-1)*k,1d0)
        endif

        if (mode00 .lt. 3) call fclr(clabl,ifi)
        if (iprint() .gt. 40) then
        outs = ' '
        if (outs1 .ne. ' ') then
          call awrit1('%x '//clabl//'%a: copied '//outs1//'%a%b from '//
     .      'class %i',outs,len(outs),0,jc)
        endif
        if (outs2 .ne. ' ') then
          call awrit0('%x '//clabl//'%a: read '//outs2//'%a%b from '//
     .      'disk',outs,len(outs),0)
        endif
        if (outs3 .ne. ' ') then
          if (outs2 .eq.  ' ') then
            call awrit0('%x '//clabl//'%a: use defaults for: '//outs3//
     .        '%a%b',outs,len(outs),0)
          else
            call awrit0('%a; use defaults for: '//outs3//'%a%b',
     .      outs,len(outs),0)
          endif
        endif
        if (outs .eq. ' ')
     .    call awrit0(' '//clabl//'%a: nothing read',outs,len(outs),0)
        call awrit0(' aiocls class'//outs,' ',-len(outs),lgunit(1))
       endif
   10 continue
C     call awrit2('%n:1i',' ',100,6,nclasp,w(ohave))
      call rlse(olpp)
C     End of MPI master-only branch
      endif

C --- MPI broadcast everything passed out of aiocls ---
C     sarray and contents of sarray pointers
      call mpibc1(w(ohave),nclasp,2,.false.,'aiocls','initc')
C     sctrl and contents
      call mpibc1(zbak,2,4,.false.,'aiocls','zbak')
C     spot and and contents of spot pointers
      call mpibc1(w(opp),6*nl*nsp*nclasp,4,.false.,'aiocls','pp')
      if (lrel .eq. 2) then
        k = 4*nl*2*nl*2*2
        call mpibc1(w(opprel),k*nclasp,4,.false.,'aiocls','pp')
      endif
C     if (igets('ctrl lncol',sctrl) .gt. 0) then
      if (osop .gt. 0) then
       call mpibc1(w(osop),9*nl*nsp*nsp*nclasp,4,.false.,'aiocls','sop')
      endif
      if (ogrrme .gt. 0) then
       call mpibc1(w(ogrrme),8*nl*nsp*nclasp,4,.false.,'aiocls','grrme')
      endif
      call mpibc1(w(opnu),nl*nsp*nclasp,4,.false.,'aiocls','pnu')
      call mpibc1(w(oqnu),3*nl*nsp*nclasp,4,.false.,'aiocls','qnu')
      call mpibc1(w(orhrmx),nclasp,4,.false.,'aiocls','rhrmx')
      nclspp = 2*nclasp-nclass
      call mpibc1(w(ovrmax),2*nclspp,4,.false.,'aiocls','vrmax')
      if (ova .gt. 0) then
        k = nclasp*(nl*nsp)**2
        call mpibc1(w(ova),k,4,.false.,'aiocls','vintr')
      endif
C     if (lgors('ctrl lasa,32',sctrl)) then
      if (omp .gt. 0) then
        k = (2*nl-1)*nl**2*3*nsp*nclasp
        call mpibc1(w(omp),k,4,.false.,'aiocls','pmpol')
      endif
      call mpibc1(w(oves),nclspp,4,.false.,'aiocls','ves')

C --- Assemble background rho (rmt) ---
Cgetarg       if (mod(mode/10,10) .ne. 0 .and. lgors('ctrl lasa,64',sctrl)) then
       if ( mod ( mode/10 , 10 ) .ne. 0 .and. iand(64,int(sctrl%lasa)) .ne.0 
     . ) then 

Cgetarg...info...           real(8) :: zbak(2) , real(8) :: sctrl%zbak 2
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg         call upack('ctrl zbak nbas nclass',sctrl,zbak,nbas,nclass,0,0)
         
         i_copy_size=size(sctrl%zbak) 
         call dcopy(i_copy_size,sctrl%zbak,1,zbak,1) 
         nbas=sctrl%nbas
         nclass=sctrl%nclass

Cgetarg...info...           integer :: onrc , integer(8) :: sarray%onrc 1
Cgetarg         call upack1('array onrc',sarray,onrc)
         
         onrc=sarray%onrc

C        call awrit2('%n:1d',' ',100,6,nclass,w(orhrmx))
C        call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
        if (zbak(2) .eq. 0) then
          do  15  ic = 1, nclass
            z = dval(w(oclabl),ic)
            call r8tos8(z,clabl)
            z = dval(w(orhrmx),ic)
            if (iprint() .ge. 30 .and. z .eq. 0) print *,
     .        ' (AIOCLS, mtcor): rho(rmax) for class ',clabl,' is 0'
            zbak(2) = zbak(2) + (w(onrc+ic-1)*z)/nbas
   15     continue
Cgetarg           zbak(2) = zbak(2) * dgets('lat vol',slat)
           zbak ( 2 ) = zbak ( 2 ) * slat%vol 

        endif
Cgetarg...info...           real(8) :: zbak(2) , real(8) :: sctrl%zbak 2
Cgetarg         call pack1('ctrl zbak',sctrl,zbak)
         
         i_copy_size=size(sctrl%zbak) 
         call dcopy(i_copy_size,zbak,1,sctrl%zbak,1) 

      endif

Cgetarg       end
       end subroutine aiocls 

      subroutine pvaioc(initc,mask,ic0,icmap,lhave,ic)

C- Find whether data avaiable either in class or mapped class
C     implicit none
      logical lhave
      integer initc(1),mask,ic0,icmap
      integer ic
      ic = ic0
      lhave = mod(initc(ic)/mask,2) .eq. 1
      if (.not. lhave .and. icmap .ne. 0) then
        lhave = mod(initc(icmap)/mask,2) .eq. 1
        ic = icmap
      endif
Cgetarg       end
       end subroutine pvaioc 

      subroutine paioc2(nsp,nl,n0,pat,qat,pnu,qnu)

C- Widget to copy pat,qat to pnu,qnu
C     implicit none
      integer n0,nl,nsp
      double precision pat(n0,2),qat(n0,2),pnu(nl,nsp),qnu(3,nl,nsp)
      integer i,il

      do  i = 1, nsp
      do  il = 1, nl
C       pnu(il,i) = int(pat(il,i)) + .5d0
        pnu(il,i) = pat(il,i)
        qnu(1,il,i) = qat(il,1)/nsp
        if (nsp .eq. 2) then
          if (pat(il,i) .eq. 0) pnu(il,i) = pnu(il,1)
          qnu(1,il,i) = qat(il,1)/nsp + qat(il,2)/2*dble(3-2*i)
        endif
        qnu(2,il,i) = 0d0
        qnu(3,il,i) = 0d0
      enddo
      enddo
Cgetarg       end
       end subroutine paioc2 


