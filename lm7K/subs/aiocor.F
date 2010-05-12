      logical function aiocor(nr,nsp,a,rmax,rhoc,sumec,sumtc,ifi)
C- File I/O for core charge density.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci   nr,nsp,a,rmax
Ci   rhoc, if file write
Co Outputs
Co   rhoc, if file read
Cr Remarks
Cr    Format for core in atomic file begins with
Cr    category 'CORE:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, aiocor returns the value true and rhoc if:
Cr       the category is found, and
Cr       the file's value of a and nr match input and rmax is close to
Cr         file's value, and
Cr       the density is read without error.
Cu Updates
Cu   26 Apr 03 Added MPI calls.  Does not broadcast file read
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,rhoc(nr,nsp),sumec,sumtc
C Local parameters
      integer i,isp,nr2,nsp2
      integer mpipid,procid
      double precision a2,rmax2
      logical scat

      procid = mpipid(1)

      if (ifi .gt. 0) then
        aiocor = .false.
        if (procid .eq. 0) then
          if (.not. scat(ifi,'CORE:',':',.true.)) goto 15
          read(ifi,102,err=15) nr2,nsp2,a2,rmax2,sumec,sumtc
          if (nr .eq. 0) nr = nr2
          if (nsp .eq. 0) nsp = nsp2
          if (a .eq. 0) a = a2
          if (rmax .eq. 0) rmax = rmax2
          if (a2 .ne. a .or. nr .ne. nr2 .or.
     .      dabs(rmax2-rmax) .gt. 1d-5) goto 15
          do  10  isp = 1, min0(nsp2,nsp)
            read(ifi,101) (rhoc(i,isp), i=1,nr)
   10     continue
          if (nsp .lt. nsp2) call dscal(nr,2d0,rhoc,1)
          if (nsp .gt. nsp2) then
            call dscal(nr,.5d0,rhoc,1)
            call dcopy(nr,rhoc,1,rhoc(1,2),1)
          endif
          aiocor = .true.
        endif
   15   continue
      else
        if (procid .eq. 0) then
          write(-ifi,'(''CORE:'')')
          write(-ifi,102) nr,nsp,a,rmax,sumec,sumtc
          do  11  isp = 1, nsp
            write(-ifi,101) (rhoc(i,isp),i = 1,nr)
   11     continue
        endif
        aiocor = .true.
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f12.5,2f18.8)
      end

