      integer function iomoms(ifi,nl,nsp,nspc,nkp,ldim,nfstg,ldwt3,iq,
     .isp,ndev,nevmx,nchan,nchan2,nev,eb,accwt,doswt,dosw2,efermi,vmtz)
C- File I/O for binary moments file.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi   :logical unit: positive for read, negative for write
Ci   nl    :(global maximum l) + 1 (not used if iq>0)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Co   nfstg: a compound of digits specifying what is to be I/O
Co        : in moments file; See Remarks
Ci         :1s digit =>  eb I/O
Ci         :10s digit => doswt I/O; digit's number specifies number of
Ci                       'spin' channels to I/O.
Ci                       See Remarks in iomomn for how this number
Ci                       is set, and the meaning of each channel.
Ci         :100s digit => dosw2 I/O
Ci         :1000s digit => I/O in parallel mode
Ci   ldwt3 :3rd dimension of doswt
Ci   iq    :index to current qp and spin
Ci   isp   :index to current spin (usind
Ci   ndev  :dimensions eb
Ci   nevmx :dimensions accwt,doswt
Ci   nchan :number of channels (l+m+class) for accwt,doswt
Ci   nchan2:number of channels for dosw2
Ci   nspc:  1 for uncoupled spins, 2 for coupled spins
Cio Inputs/Outputs.
Cio   iomoms reads or writes one of the these groups, depending on iq
Cio     iq=0:          nl,nsp,nspc,nkp,ldim,nfstg
Cio     0<iq<nspc*nkp: one or more of eb,accwt,doswt,dosw2
Cio                    belonging to this qp and spin (see Remarks)
Cio     iq=nspc*nkp+1: vmtz,efermi
Cio Outputs
Co    iomoms: 0 if read successful, -1 if not
Co    nev:    number of evals read at this k-point:
Co            required that nev <= nevmx
Cr Remarks
Cr   iomoms reads in one set of data (some scalars, and optionally evals,
Cr   moments, spin density matrix...) from the secular matrix for one
Cr   q-point and spin channel.  For the coupled spin case, the two spin
Cr   channels are read as one, since they cannot be uncoupled.
Cr   The moments file has the following records, depending on nfstg:
Cr     1.   nl nsp nspc nkp ldim nfstg
Cr     ...  For each qpt and spin:
Cr     2.   nchan  nev (if nfstg nonzero)
Cr          eb         (if 1s    digit of nfstg nonzero)
Cr          doswt      (if 10s   digit of nfstg nonzero)
Cr                     (this digit specifies the number of 'spin')
Cr                     (channels to read or write; see iomomn below
Cr                     (for how they are specified).
Cr          dosw2      (if 100s  digit of nfstg 1)
Cr     3.   efermi, vmtz
Cr   Caller may read less information than is available.
Cr   iomoms fails if caller tries to read more than is available.
Cr   On first call (iq=0), iomoms checks to ensure match between
Cr   passed nl nsp nspc nkp ldim nfstg and file values.
Cu Updates
Cu   07 Jul 08 Extra dimensioning parameter ndev to decouple dimension of
Cu             eb from that of accwt,dowst
Cu    9 Jan 06 file write, 'parallel mode'
Cu    7 Apr 04 dosw2 can have different nchan than doswt
Cu   18 Jan 02 Cleaner conventions involving 10s digit nfstg.
Cu             Old accwt and doswt merged into one array.
Cu    7 Dec 00 Added space for orbital moments
Cu    5 Dec 95 Use nevmx for dimensioning eb etc, rather than old nband
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,nl,nsp,nspc,nkp,ldim,nfstg,iq,ndev,nevmx,nchan,nchan2,
     .nev,ldwt3
      double precision eb(ndev,iq),efermi,vmtz,
     .accwt(nchan,nevmx,nspc,2+nspc,iq),doswt(nchan,nevmx,ldwt3,iq),
     .dosw2(2*nchan2*nevmx*3,iq)
C ... Local parameters
      integer nlf,nspf,nspcf,nkpf,ldimf,nfstgf,iprint,nchani,j,
     .nschan,nschnf,isp
      logical lddump
      double precision xx
      character outs*80
      save nfstgf

C ------------------ File read -----------------
      if (ifi .gt. 0) then
        iomoms = -1
C   --- Read dimensions ---
        if (iq .eq. 0) then
          read (ifi,err=999,end=999) nlf,nspf,nspcf,nkpf,ldimf,nfstgf
          if (nl .ne. nlf .or. nsp .ne. nspf .or. nspc .ne. nspcf .or.
     .    nkp .ne. nkpf .or. ldim .ne. ldimf) then
            print '(a)', '   IOMOMS:  nl nsp nspc  nkp ldim'
            print 333, ' expected ', nl, nsp, nspc, nkp, ldim
            print 333, 'but found ',nlf,nspf,nspcf,nkpf,ldimf
  333       format(1x,a,i3,i4,i5,i5,i5,i4)
            call rx('IOMOMS: bad input')
          endif
          if (iprint() .gt. 40) then
            outs = ' IOMOMS: file contains:'
            if (nfstgf .eq. 0) call awrit0('%a  nothing',outs,80,0)
            if (mod(nfstgf,10) .eq. 1)
     .      call awrit0('%a  bands',outs,80,0)
            if (mod(nfstgf/10,10) .ge. 1)
     .      call awrit0('%a  doswts',outs,80,0)
            if (mod(nfstgf/10,10) .ge. 3)
     .      call awrit0('%a  moments',outs,80,0)
            if (mod(nfstgf/100,10) .ge. 1)
     .      call awrit0('%a  spin density matrix',outs,80,0)
            call info0(30,0,0, outs//'%a')
            call info5(30,0,0,'         nl=%i  nsp=%i  nspc=%i  '//
     .      'nkp=%i  ldim=%i',nlf,nspf,nspcf,nkpf,ldimf)
          endif

C   --- Read data for this qp ---
        elseif (iq .le. nkp*nsp .and. nfstgf .eq. 0) then
          call rx('IOMOMS: file is missing bands')
        elseif (iq .le. nkp*nsp .and. nfstg .ge. 0) then
          read (ifi,err=999,end=999) nchani,nev
          if (nchani .ne. nchan .and. nfstg .gt. 1)
     .    call rx('IOMOMS: nchan mismatch')
          if (nev .gt. nevmx .and. nfstg .gt. 0)
     .    call rx('IOMOMS: nev gt nevmx')
C     ... Read eb
          if (nfstg .eq. 0) then
            read (ifi,err=999,end=999) xx
          else
            if (.not. lddump(eb(1,iq),nev,ifi)) goto 999
          endif
C     ... Read, or read over doswt
          nschan = mod(nfstg/10,10)
          nschnf = mod(nfstgf/10,10)
          if (nschan .gt. nschnf)
     .    call rx('IOMOMS: file has fewer doswts than required')
          if (nschnf .gt. 0) then
            do  j = 1, nschnf
              if (j .le. nschan) then
                if(.not.lddump(doswt(1,1,j,iq),nchani*nev,ifi)) goto 999
              else
                read (ifi,err=999,end=999) xx
              endif
            enddo
          endif
C     ... Read, or read over rhos
          if (mod(nfstg/100,10) .eq. 0) then
            if (mod(nfstgf/100,10) .eq. 1) read(ifi,err=999,end=999) xx
          elseif (mod(nfstg/100,10) .eq. 1) then
            if (mod(nfstgf/100,10) .eq. 1) then
              if (.not. lddump(dosw2(1,iq),2*nchan2*nev*3,ifi)) goto 999
C              if (nev .ne. nevmx) call rx('iomoms can''t read dosw2')
            else
              call rx('IOMOMS: file missing spin-density-matrix')
            endif
          endif
        elseif (iq .eq. nkp*nsp+1) then
          read (ifi,err=999,end=999) efermi, vmtz
        endif
        iomoms = 0
  999   return

C ------------------ File write -----------------
      else
        iomoms = 1
        if (iq .eq. 0) then
          write (-ifi) nl, nsp, nspc, nkp, ldim, nfstg
C       write, Parallel mode
        elseif (nfstg .gt. 1000) then
          write (-ifi) nchan,nev,iq,isp
          call dpdump(eb,nev,ifi)
          nschan = mod(nfstg/10,10)
          do  j = 1, nschan
            call dpdump(doswt(1,1,j,1),nchan*nev,ifi)
          enddo
          if (mod(nfstg/100,10) .eq. 1) then
            call dpdump(dosw2(1,1),2*nchan2*nev*3,ifi)
          endif
C       write, Serial mode
        elseif (iq .lt. nkp*nsp .and. nfstg .gt. 0) then
          write (-ifi) nchan,nev
          call dpdump(eb(1,iq),nev,ifi)
          nschan = mod(nfstg/10,10)
          do  j = 1, nschan
            call dpdump(doswt(1,1,j,iq),nchan*nev,ifi)
          enddo
          if (mod(nfstg/100,10) .eq. 1) then
            call dpdump(dosw2(1,iq),2*nchan2*nev*3,ifi)
          endif
        elseif (iq .eq. nkp*nsp+1) then
          write (-ifi) efermi, vmtz
        endif
      endif
      end
      subroutine iomomn(metal,ldos,lrhos,nevmx,nspc,nmom,lwtkb,nfstg)
C- Returns nfstg appropriate to conditions specified by caller
C ----------------------------------------------------------------------
Ci Inputs
Ci   metal :T => metal, F => nonmetal
Ci   ldos =a combination of the following integers
Ci         1 make dos
Ci         2 generate weights for partial dos
Ci         4 generate weights for m-decompos'n of pdos
Ci   lrhos :make off-diagonal (in spin space) spin-density matrix
Ci   nevmx :largest number of eigenvectors to find
Ci         :nevmx<0 => no eigenvectors are made
Ci   lwtkb :0 weights are neither required a priori nor available
Ci         :1 weights are required a priori, and are available
Ci         :  (having read from disk or from a prior band pass)
Ci         :-1 weights are required a priori, but were not read
Ci   nmom  :number of energy moments + orbital moments
Ci         :typical numbers:
Ci         :  ASA collinear case, 3 (energy moments)
Ci         :  ASA noncollear case, 2 or 4 (energy moms+orbital moment)
Co Outputs
Co   nfstg: a compound of digits specifying what is to be I/O
Co        : in moments file
Co        :1s digit   energy bands I/O
Co        :10s digit  doswt I/O; digit's number specifies number of
Co                    (moments + spin channels) to I/O (nmom*nspc)
Co                    Only set if metal=T and ldos>1.
Co        :100s digit => dosw2 I/O
Cr Remarks
Cr
Cu Updates
Cu   08 Jun 06 separated nmom to be computed externally
Cu   16 Jan 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical metal,lrhos
      integer nfstg,ldos,lwtkb,nevmx,nspc,nmom
C ... Local parameters

C     Save bands only
      nfstg = 1
C     10s digit for the number of different moments * spin channels
      if (metal .or. ldos .gt. 1) then
        nfstg = 1 + 10*nmom*nspc
      endif
      if (lrhos) nfstg = nfstg + 100
C     No eigenvectors; save bands only
      if (nevmx .lt. 0 .or. lwtkb .eq. -1) nfstg = 1
      end

