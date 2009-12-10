      subroutine iomomx(ifi,mode,nl,nsp,nspc,nkp,ldim,nfstg,iq,
     .  nband,nchan,nchan2,nchds,nevmx,nlistc,listc,eb,accwt,doswt,
     .  dosw2,efermi,vmtz)
C- Read selected channels of data from moments file
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ifi   :file handle
Ci   mode  : 1's digit 0 read nl,nsp,nspc,nkp,ldim,nfstg
Ci         :           1 require a match in nl,nsp; read nspc,ldim,nkp
Ci         :           2 require a match in nl,nsp,nspc,nkp,ldim
Ci         :10's digit 0 exit after reading header info
Ci         :       1,2,3 read number of iq before EOF encountered
Ci         :         2,3 read all information sought by nfstg
Ci         :           3 read efermi,vmtz if available
Ci   nband :leading dimension of eb,dosw2,doswt,accts
Ci   nchds :number of DOS channels to read
Ci   nlistc:number of classes, species or sites (see Remarks)
Ci   listc :class, species or site list
Ci   nchan :number of channels (l+m+class) for accwt,doswt
Ci   nchan2:number of channels for dosw2
Cio Inputs/Outputs
Cio    ... The following are read for 1s digit mode=0; else they are input
Cio  nl    :(global maximum l) + 1
Cio  nsp   :2 for spin-polarized case, otherwise 1
Cio  nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Cio  nkp   :number of irreducible k-points (bzmesh.f)
Cio  ldim  :dimension of hamiltonian matrix (makidx.f)
Cio  nfstg :describes information contained in moments file (see iomoms).
Cio        :iomomq passes nfstg to iomoms; moments file (ifi) may contain
Cio        :more information than sought by nfstg; but if some information
Cio        :is missing and 10's digit of mode is 2 or 3, iomomq aborts.
Cio        :For 1s digit mode > 0, nfstg is an input.
Co Outputs:
Co   iq     :(1s digit mode>0) number of qp available in moments file.
Co   nevmx  :largest # evals encountered
Cf  Files:
Cf    Moments file has the following records
Cf    1.   nl nsp nspc nkp ldim nfstg
Cf    ... For each qpt (and spin), the following records:
Cf    2.   nchan  nev (if nfstg nonzero)
Cf         eband      (if 1s   digit of nfstg nonzero)
Cf         accwt      (if 10s  digit of nfstg 1)
Cf         doswt      (if 10s  digit of nfstg 2)
Cf         dosw2      (if 100s digit of nfstg 1)
Cf    3.   efermi, vmtz
Cb Bugs
Cb   dosw2 shouldn't be dimensioned nchds!
Cr Remarks
Cr   Same as iomomq, but reads a subset of DOS channels.
Cr   The dos channels are assumed be ordered by class, species, or site
Cr   (it doesn't matter for this routine's purposes) with nl channels
Cr   per class, species, or site.
Cu Updates
Cu    7 Apr 04 dosw2 can have different nchan than doswt
Cu   27 Apr 02 Bug fix
Cu   18 Jan 02 Cleaner conventions involving 10s digit nfstg.
Cu             Old accwt and doswt merged into one array.
C ----------------------------------------------------------------------
C     implicit none
      integer ifi,mode,nl,nsp,nspc,nkp,ldim,nfstg,iq,nband,nchan,nchan2,
     .  nchds,nevmx,nlistc,listc(nlistc)
      double precision eb(nband*nsp,1),doswt(nchds,nband*nspc,nspc,1),
     .  accwt(nchds,nband*nspc,nspc,3,1),dosw2(2*nchds*nband*nspc*3,1),
     .  efermi,vmtz
C Local variables
      double precision xx
      integer nlf,nspf,nspcf,nkpf,ldimf,nfstgf,iprint,i1mach,jq,iomoms,
     .  nev,oeb,oaccm,odosw,nschan
      character outs*80
C ... Heap
      integer w(1)
      common /w/ w

C --- Read header, checking req'd matches and copying the rest ---
      rewind ifi
      if (mod(mode,10) .eq. 0) then
        read (ifi,err=999,end=999) nl, nsp, nspc, nkp, ldim, nfstg
        if (mode .eq. 0) return
      else
        read (ifi,err=999,end=999) nlf, nspf, nspcf, nkpf, ldimf, nfstgf
        if (mod(mode,10) .eq. 1) then
          nkp  = nkpf
          nspc = nspcf
          ldim = ldimf
        endif
      endif
      iq = 0
      rewind ifi
      if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,nfstg,1,iq,1,nevmx,nevmx,
     .  nchan,nchan2,nevmx,eb,accwt,doswt,dosw2,efermi,vmtz) .lt. 0)
     .  goto 999
      if (mod(mode/10,10) .eq. 0) return

C --- Determine number of qp available ---
      iq = 0
      do  12  jq = 1, nkp*(nsp/nspc)
        if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,0,1,jq,1,nband*nspc,
     .    nband*nspc,nchan,nchan2,nev,eb,accwt,doswt,dosw2,efermi,vmtz)
     .    .lt. 0) goto 14
        iq = jq-1
   12 continue
C ... If read all qp, and another record also present, let iq=nkp
      read (ifi,err=14,end=14) xx
      iq = iq+1
   14 continue
      iq = (nspc*iq)/nsp
      call awrit1(' IOMOMX: read %i qp',outs,80,0,iq)
      if (mod(mode/10,10) .eq. 1) goto 99

C --- Read info spec'd by nfstg for each qp until error ---
      nfstgf = nfstg
      nschan = mod(nfstg/10,10)
      rewind ifi
      read (ifi,err=999,end=999) nlf
      nevmx = 0
      call defrr(oeb,  nband)
      call defrr(oaccm,nchan*nband*nfstg)
      call defrr(odosw,nchan*nband*nfstg)
      if (nfstg .ge. 100) call rxi('iomomx not ready for nfstg=',nfstg)
      do  10  jq = 1, iq*(nsp/nspc)
        if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,nfstg,nschan,1,1,nband*nspc,
     .    nband*nspc,nchan,nchan2,nev,w(oeb),w(oaccm),w(odosw),dosw2,
     .    efermi,vmtz).lt. 0)goto 999
        call pviomx(nchan,nchds,nband*nspc,nl,nspc,listc,nlistc,nfstg,
     .    nschan,jq,w(oeb),w(odosw),eb,doswt)
        nevmx = max(nev,nevmx)
   10 continue
      call rlse(oeb)

      if (mod(mode/10,10) .eq. 2 .or. iq .ne. nkp) goto 99

C --- Read efermi, vmtz if sought and info available ---
      read (ifi,err=999,end=999) efermi,vmtz
      call awrit2('%a  efermi=%,6d  vmtz=%,6d',outs,80,0,efermi,vmtz)

C --- Exit ---
   99 continue
      if (iprint() .lt. 30) return
      call awrit0('%a from moms file',outs,80,-i1mach(2))
      return

C --- Error exit ---
  999 if (iprint() .gt. 10)
     .  print *, 'IOMOMX (warning): empty or bad moments file'
      iq = -1
      end

      subroutine pviomx(nchan,nchds,ndim,nl,nspc,listc,nlistc,nfstg,
     .  ldwt3,iq,eb0,dosw0,eb,doswt)
C- Copy a subset of channels
C     implicit none
      integer nchan,nchds,ndim,nspc,nl,ldwt3,nfstg,iq,nlistc,
     .  listc(nlistc)
      double precision eb0(ndim),dosw0(nchan,ndim,nspc)
      double precision eb(ndim,iq),doswt(nchds,ndim,ldwt3,iq)
      integer il,ic,ichan,ichds,imax,i
      if (nspc .ne. 1) call rx('pviomx implemented only for nspc=1')

      il = 1
      ichan = -nl
      ichds = -nl
      imax = listc(nlistc)
      if (mod(nfstg,10) .eq. 1) call dcopy(ndim,eb0,1,eb(1,iq),1)
      do  10  ic = 1, imax
        ichan = ichan + nl
        if (ic .lt. listc(il)) goto 10
        ichds = ichds + nl
        if (mod(nfstg/10,10) .ge. 1) then
          do  22  i = 1, nl
            call dcopy(ndim*nspc,dosw0(ichan+i,1,1),nchan,
     .        doswt(ichds+i,1,1,iq),nchds)
   22     continue
        endif

        il = il+1
   10 continue
      end

