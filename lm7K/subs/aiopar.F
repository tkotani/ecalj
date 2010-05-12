      logical function aiopar(alabl,lrel,pp,pprel,ves,nl,lmax,nsp,ifi)
C- File I/O for potential parameters.
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl :class label
Ci   nl    :(global maximum l) + 1
Ci   lmax  :maximum l for a given site
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio   pp   :potential parameters (atomsr.f)
Cio  ves   :estat potential at rmax
Cio  aiopar:true unless read error or category not found
Cr Remarks
Cu Updates
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   26 Apr 03 Added MPI calls
Cr   11 Apr 94 Added convention: reading 1/p=0 => set p to 0
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character alabl*8
      integer nl,lmax,nsp,ifi,lrel
      double precision ves,pp(6,nl,nsp),pprel(4,0:nl-1,2*nl,2,2)
C Local parameters
      double precision xx,x2(0:7),mu
      equivalence (xx,x2)
      integer i,l,ip,ll,k,iprint,i1mach,a2vec,ix2(8),nl2,ipr,nsp2,n,lmx
      integer mpipid,procid,imu,m,mumax
      logical scat,sw,rdstrn,pars1v
      character s*100

      call getpr(ipr)
      aiopar = .false.
      sw = .false.
      procid = mpipid(1)

C --- File READ ---
      if (ifi .gt. 0) then
        if (procid .eq. 0) then

C   ... return unless file has PPAR: category
          if (.not. scat(ifi,'PPAR:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
          backspace ifi
          read(ifi,'(a72)') s
          if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
          if (lmax .ne. nl2-1 .and. ipr .ge. 10)
     .    print *, 'aiopar (warning) mismatch in nl, class '//alabl
          if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
          if (nsp .ne. nsp2 .and. ipr .ge. 10)
     .    print *, 'aiopar (warning) mismatch in nsp, class '//alabl
          if (.not. pars1v(s,len(s),'ves=','=',4,ves)) 
     .    print *, 'aiopar (warning) failed to read ves, class '//alabl
          call dpzero(pp,6*nl*nsp)
          lmx = min(nl-1,nl2-1,lmax)
C   ... read ppar
          n = min(nsp,nsp2)
          read(ifi,'(a72)') s
          do  10  i = 1, n
            do  16  l = 0, nl2-1
              if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
              if (l .gt. lmx) goto 16
              ip = 0
              k = a2vec(s,len(s),ip,4,' ',1,-2,-7,ix2,x2)
C    ...  Abort if failed to read 7 numbers
              if (k .ne. 7)
     .      call rx('AIOPAR, class'//alabl//'%a: failed to parse '//s)
              ll = x2(0)
C     ... no error if file has extra l's
              if (ll .gt. l) goto 16
C     ... but pot pars must be available up to lmax
              if (ll .ne. l) call rx('AIOPAR: bad l quantum number')
C Map delta into sqrt(delta), preserving sign and 1/sqrt(p) into p
              x2(3) = dsign(1.d0,x2(3))*dsqrt(dabs(x2(3)))
              if (x2(4) .eq. 0) then
                if (iprint() .gt. 40)
     .        print *, 'aiopar: encountered 1/p=0 ... set p to 0'
              else
                x2(4) = 1/x2(4)**2
              endif
              do  15  k = 1, 6
                pp(k,l+1,i) = x2(k)
                pp(k,l+1,nsp) = x2(k)
   15         continue
   16       continue
   10     continue
        endif
C     Read relativistic ppar's
        if (procid .eq. 0 .and. lrel .eq. 2) then
C   ... return unless file has PPREL: category
          if (.not. scat(ifi,'PPREL:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
          backspace ifi
          read(ifi,'(a72)') s
          if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
          if (lmax .ne. nl2-1 .and. ipr .ge. 10)
     .    print *, 'aiopar (warning) mismatch in nl, class '//alabl
          if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
          if (nsp .ne. nsp2 .and. ipr .ge. 10)
     .    print *, 'aiopar (warning) mismatch in nsp, class '//alabl
          if (.not. pars1v(s,len(s),'ves=','=',4,ves)) 
     .    print *, 'aiopar (warning) failed to read ves, class '//alabl
          call dpzero(pprel,4*nl*2*nl*4)
          lmx = min(nl-1,nl2-1,lmax)
C   ... read ppar
          n = min(nsp,nsp2)
          read(ifi,'(a72)') s

          do  l = 0, lmx
            mumax = 2*(l+1)
            do  imu = 1, mumax
              do  i = 1, 2
                do  k = 1, 2
                  if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
                  if (l .gt. lmx) goto 116
                  ip = 0
                  m = a2vec(s,len(s),ip,4,' ',1,-2,-8,ix2,x2)
C               Abort if failed to read 8 numbers
                  if (m .ne. 8) call rx('AIOPAR (rel) class'//alabl//
     .            '%a: failed to parse '//s)
                  mu = dble(imu-l) - 1.5d0
                  if (dabs(mu-x2(1)) .gt. 1d-12) goto 18
                  ll = x2(0)
C           ... no error if file has extra l's
                  if (ll .gt. l) goto 116
C           ... but pot pars must be available up to lmax
                  if (ll .ne. l) call rx('AIOPAR: bad l quantum number')
                  do  m = 1, 4
                    pprel(m,l,imu,i,k) = x2(m+3)
                  enddo
                enddo
              enddo
            enddo
          enddo
  116     continue

        endif

C     call mpibc1(pp,6*nl*nsp,4,.false.,'aiopar','pp')
C     call mpibc1(ves,1,4,.false.,'aiopar','ves')

      else
        write(-ifi,21) alabl, lmax+1, nsp, ves
        do  11  i = 1, nsp
        do  11  l = 0, lmax
          xx = pp(3,l+1,i)
          if (pp(4,l+1,i) .eq. 0) then
            write(-ifi,20) l, (pp(k,l+1,i), k=1,2), xx**2*dsign(1d0,xx),
     .        0d0, (pp(k,l+1,i), k=5,6)
            sw = .true.
          else
            write(-ifi,20) l, (pp(k,l+1,i), k=1,2), xx**2*dsign(1d0,xx),
     .        1/dsqrt(pp(4,l+1,i)), (pp(k,l+1,i), k=5,6)
          endif
   11   continue
        if (iprint() .ge. 30  .and. sw)
     .    write(i1mach(2),*) 'aiopar: encountered p=0 ... wrote 1/p=0'
        aiopar = .true.

C       Write relativistic ppar's
        if (lrel .eq. 2) then
          write(-ifi,22) alabl, lmax+1, nsp, ves
          do  l = 0, lmax
            mumax = 2*(l+1)
            do  imu = 1, mumax
              mu = dble(imu-l) - 1.5d0
              do  i = 1, 2
                do  k = 1, 2
                  write(-ifi,23) l,mu,i,k, (pprel(m,l,imu,i,k), m=1,4)
                enddo
              enddo
            enddo
          enddo
        endif

      endif

      aiopar = .true.
      return

   18 continue
      if (iprint() .gt. 0) then
        print *, 'aiopar: (input skipped) bad syntax, class '//alabl
      endif
      return

   20 format(i2,3f12.8,f12.7,3f12.8)
   23 format(i3,f5.1,2i4,4f13.8)
   21 format('PPAR:  ',a,'  nl=',i1,'  nsp=',i1,'  ves=',f12.8/' l',5x,
     .  'e_nu',10x,'C',8x,'+/-del',5x,'1/sqrt(p)',6x,'gam',9x,'alp')
   22 format('PPREL:  ',a,'  nl=',i1,'  nsp=',i1,'  ves=',f12.8/
     .  '  l  mu   ms1 ms2',5x,'C',12x,'W',12x,'gam',10x,'p')

      end

