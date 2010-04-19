      subroutine iores(nRLc,nsp,nspc,ifi,p)
C- Input/output of response matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nRLc  :dimension of response matrix
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ifi   :file handle
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Co Outputs
Co   p     :is read from disk
C ----------------------------------------------------------------------
C     implicit none
      integer nRLc,nsp,nspc,ifi
      double precision p(nRLc,nspc,nRLc,nsp)
C Local variable
      integer rdm,nspf,nspcf,nr,nrf,ncf,nRLcc,i,j,lgunit
      logical pars1v,parstr
      character*100 s,outs
      data nspf /0/


C --- Setup ---
      call dpzero(p,nRLc*nspc*nRLc*nsp)
      rewind ifi
C ... Pick up file values of nsp, nspc
      call rdfile(ifi,'#{}%',s,1,outs,len(s),nr)
C     strip '"' from label
      i = 0
      if (parstr(s,'l="',len(s),3,'"',i,j)) then
        s(j:j) = ' '
        if (parstr(s,'"',len(s),1,'"',j,i)) then
          s(i:i) = ' '
        endif
      endif
      if (.not. pars1v(s,len(s),'nsp=','=',2,nspf)) nspf = 1
      if (.not. pars1v(s,len(s),'nspc=','=',2,nspcf)) nspcf = 1
C     Require internal dimensions of matrix match nRLc
      nrf = 0
      ncf = 0
      rewind ifi
      if (rdm(ifi,0,0,' ',p,nrf,ncf) .ne. 1) call rx('IORES:  bad file')
      if (nRLc*nspcf .ne. nrf .or. nrf .ne. ncf)
     .  call rx('IORES:  file mismatch')
      rewind ifi

      nRLcc = nRLc*nspc

C ... Case file nspcf matches input nspc
      if (nspc .eq. nspcf) then
        if (rdm(ifi,0,nRLcc*nRLcc,' ',p,nRLcc,nRLcc) .ne. 1) call
     .    rx('IORES: failed to read response matrix')
C       The next two cases can only occur when nspc=1
        if (nsp .eq. 2 .and. nspf .eq. 1) then
          call dscal(nRLc**2,.5d0,p,1)
          call dpcopy(p,p(1,1,1,2),1,nRLc**2,1d0)
        elseif (nsp .eq. 2 .and. nspc .eq. 1) then
          if (rdm(ifi,0,nRLc*nRLc,' ',p(1,1,1,2),nRLc,nRLc) .ne. 1)
     .      call rx('IORES: failed to read response matrix')
        endif

C ... Case file nspcf is less than nspc : read and copy
      elseif (nspc .gt. nspcf) then
C       Read 1st spin channel and copy to (1,2) block
        if (rdm(ifi,0,nRLc*nRLc,' ',p,nRLc,nRLc) .ne. 1)
     .    call rx('IORES: failed to read response matrix')
        if (nsp .eq. 2 .and. nspf .eq. 1) call dscal(nRLc**2,.5d0,p,1)
        call dmscop(p(1,1,1,2),nRLcc,p,nRLc,1,nRLc,1,nRLc,1,1,1d0)
C       Read 2nd spin channel and copy to (2,2) block
        if (.not. (nsp .eq. 2 .and. nspf .eq. 1)) then
        if (rdm(ifi,0,nRLc*nRLc,' ',p,nRLc,nRLc) .ne. 1)
     .    call rx('IORES: failed to read response matrix')
        endif
        call dmscop(p(1,2,1,2),nRLcc,p,nRLc,1,nRLc,1,nRLc,1,1,1d0)
C       Copy to (1,2) block to (1,1) block
        call dmscop(p,nRLcc,p(1,1,1,2),nRLcc,1,nRLc,1,nRLc,1,1,1d0)
C       Zero out(1,2) block
        do  20  i = 1, nRLc
        do  20  j = 1, nRLc
        p(i,2,j,1) = 0
   20   p(i,1,j,2) = 0
      elseif (nspc .lt. nspcf) then
        call rx('IORES:  file mismatch : file contains nspc=2')
      endif

C ... debugging
C      call prm('p, spin %i',1,p,nRLcc,nRLcc,nRLcc)
C      if (nsp .eq. 2 .and. nspc .eq. 1)
C     .  call prm('spin %i',2,p(1,1,1,2),nRLc,nRLc,nRLc)

      call awrit2(' iores: read psta, file nsp=%i nspc=%i',' ',80,
     .  lgunit(1),nspf,nspcf)

      end

