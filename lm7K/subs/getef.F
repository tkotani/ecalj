      integer function getef(ifi,isw,efermi)
C- Widget to get Fermi energy from end of BAND file
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit
Ci   isw   :compound switch
Ci         :1s digit
Ci         :  0 Return getef nonzero if fermi level is not read
Ci         :  1 Abort if fermi level is not read
Ci         :10s digit
Ci         :  0 read from dos weights style file
Ci         :  1 read from end of moments-style file
Ci         :100s digit
Ci         :  1 look for command-line override
Ci   nsp (nspc) : number of (coupled) spins
Co Outputs
Co   efermi : Fermi energy
Cr Remarks
Cr
Cu Updates
Cu   18 Apr 01 Added switches; getef is now a function
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,isw
      double precision efermi
C ... Local parameters
      logical cmdopt,a2bin
      character*80 strn
      integer iq,nl,nspx,nspcx,nkp,ldim,nfstg,nevmx,iobzwt,i
      double precision xx

C --- File read ---
      if (mod(isw,100) .le. 10) then
        i = 1
        iq = iobzwt(1,i,i,i,efermi,xx,ifi)
      else
        iq = 0
        call pshpr(0)
C       Extract the Fermi level from moments file
        call iomomq(ifi,0,nl,nspx,nspcx,nkp,ldim,nfstg,iq,1,1,1,nevmx,
     .    xx,xx,xx,xx,efermi,xx)
        call iomomq(ifi,32,nl,nspx,nspcx,nkp,ldim,0,iq,ldim,1,1,nevmx,
     .    xx,xx,xx,xx,efermi,xx)
        call poppr
        call rxx(nevmx .gt. ldim,'GETEF: file mismatch')
        rewind ifi
      endif

C --- Override Fermi level with command-line value ---
      if (isw .ge. 100) then
        i = 4
        if (cmdopt('-ef=',i,0,strn)) then
          if (.not. a2bin(strn,efermi,4,0,' ',i,-1)) call
     .      rxs2('GETEF: failed to parse "',strn(1:30),'%a"')
          iq = 1
        endif
      endif

C ... Exit if missing efermi when it is required
      if (iq .lt. 0) then
        call rxx(mod(isw,10),'GETEF: Missing Fermi energy')
        efermi = 999d0
      endif

      getef = iq
      end

