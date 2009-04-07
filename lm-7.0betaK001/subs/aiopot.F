      logical function aiopot(nr,nsp,a,rmax,bhat,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci   nr,nsp,a,rmax
Ci   bhat: orientation of magnetic field (optional)
Ci         If file write AND bhat(1) <> -99, write bhat in POT: line
Ci   v, if file write
Co Outputs
Co   v   :potential read from idst if ifi>0 (file read)
Co   bhat:orientation of magnetic field (optional)
Co        If file read AND bhat(1)<>-99, read bhat from POT: line
Cr Remarks
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, aiopot returns true only if the category is found,
Cr    the file's value of a and nr match input and rmax is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
Cu Updates
Cu    4 Apr 04 Optional I/O of bfield orientation
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,bhat(3),v(nr,nsp)
C Local parameters
      integer i,j,isp,nr2,nsp2,ix2(3),a2vec
      double precision a2,rmax2
      logical scat
      character strn*72

      aiopot = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'POT:',':',.true.)) return
        if (bhat(1) .ne. -99d0) then
          bhat(1) = 0
          bhat(2) = 0
          bhat(3) = 0
          backspace ifi
          strn = ' '
          read(ifi,'(a72)') strn
          call words(strn,i)
          if (i .ge. 5) then
            call word(strn,3,i,j)
            i = i-1
            j = a2vec(strn,len(strn),i,4,' ',1,-2,-3,ix2,bhat)
            if (j .ne. 3) call rxs(
     .        'aiopot: failed to parse Bfield in line: ',strn)
C           print *, bhat
          endif
        endif
        read(ifi,102,err=15) nr2,nsp2,a2,rmax2
        if (nr .eq. 0) nr=nr2
        if (nsp .eq. 0) nsp=nsp2
        if (a .eq. 0) a=a2
        if (rmax .eq. 0) rmax=rmax2
        if (a2 .ne. a .or. nr .ne. nr2 .or.
     .      dabs(rmax2-rmax) .gt. .03) goto 15
        do  10  isp = 1, min0(nsp2,nsp)
c          read(ifi,101,err=15) (v(i,isp),i = 1,nr)
          read(ifi,101) (v(i,isp),i = 1,nr)
          do  10  i = 1, nr
          v(i,nsp) = v(i,isp)
   10   continue
        aiopot = .true.
   15   continue
      else
        if (bhat(1) .eq. -99d0) then
          write(-ifi,'(''POT:'')')
        else
          write(-ifi,'(''POT:   bhat='',3f12.7)') bhat
        endif
        write(-ifi,102) nr,nsp,a,rmax
        do  11  isp = 1, nsp
          write(-ifi,101) (v(i,isp),i = 1,nr)
   11   continue
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f12.5)
      end

