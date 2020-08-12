      subroutine iomomq(ifi,mode,nl,nsp,nspc,nkp,ldim,nfstg,iq,
     .nband,nchan,nchan2,nevmx,eb,accwt,doswt,dosw2,efermi,vmtz)
C- Read data from moments file, depending on mode and nfstg
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   mode  :1's digit  0 read nl,nsp,nspc,nkp,ldim,nfstg
Ci         :           1 require a match in nl,nsp; read nspc,ldim,nkp
Ci         :           2 require a match in nl,nsp,nspc,nkp,ldim
Ci         :10's digit 0 exit after reading header info
Ci         :      1,2,3 read number of iq before EOF encountered
Ci         :        2,3 read all information sought by nfstg
Ci         :          3 read efermi,vmtz if available
Ci   nband :leading dimension of eb,dosw2,doswt,accts
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
Cu Updates
Cu    7 Apr 04 dosw2 can have different nchan than doswt
Cu   18 Jan 02 Cleaner conventions involving 10s digit nfstg.
Cu             Old accwt and doswt merged into one array.
Cu   21 Nov 95 Made iomomq return w/ iq=-1 if error
Cu    5 Dec 95 Adapted changed argument list in iomoms.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,mode,nl,nsp,nspc,nkp,ldim,nfstg,iq,nband,nchan,nchan2,
     .nevmx
      double precision eb(nband*nsp,1),dosw2(2*nchan2*nband*nspc*3,1),
     .doswt(nchan,nband*nspc,nspc,1),efermi,vmtz,
     .accwt(nchan,nband*nspc,nspc,2+nspc,1)
C ... Local parameters
      double precision xx
      integer nlf,nspf,nspcf,nkpf,ldimf,nfstgf,jq,iomoms,nev,nschan
      character outs*80

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
     .nchan,nchan2,nevmx,eb,accwt,doswt,dosw2,efermi,vmtz) .lt. 0)
     .goto 999
      if (mod(mode/10,10) .eq. 0) return

C --- Determine number of qp available ---
      iq = 0
      do  12  jq = 1, nkp*(nsp/nspc)
        if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,0,1,jq,1,nband*nspc,
     .  nband*nspc,nchan,nchan2,nev,eb,accwt,doswt,dosw2,efermi,vmtz)
     .  .lt.0) goto 14
        iq = jq-1
   12 continue
C ... If read all qp, and another record also present, let iq=nkp
      read (ifi,err=14,end=14) xx
      iq = iq+1
   14 continue
      iq = (nspc*iq)/nsp
      call awrit1(' IOMOMQ: read %i qp',outs,80,0,iq)
      if (mod(mode/10,10) .eq. 1) goto 99

C --- Read info spec'd by nfstg for each qp until error ---
      nfstgf = nfstg
      nschan = mod(nfstg/10,10)
      rewind ifi
      read (ifi,err=999,end=999) nlf
      nevmx = 0
      do  10  jq = 1, iq*(nsp/nspc)
        if (iomoms(ifi,nl,nsp,nspc,nkp,ldim,nfstg,nschan,jq,1,
     .  nband*nspc,nband*nspc,nchan,nchan2,nev,eb,accwt,doswt,dosw2,
     .  efermi,vmtz).lt. 0) goto 999
        nevmx = max(nev,nevmx)
   10 continue

      if (mod(mode/10,10) .eq. 2 .or. iq .ne. nkp) goto 99

C --- Read efermi, vmtz if sought and info available ---
      read (ifi,err=999,end=999) efermi,vmtz
      call awrit2('%a  efermi=%,6d  vmtz=%,6d',outs,80,0,efermi,vmtz)

C --- Exit ---
   99 continue
      call info0(30,0,0, outs//'%a')
      return

C --- Error exit ---
  999 call info0(10,0,0,' IOMOMQ (warning): empty or bad moments file')
      iq = -1
      end


