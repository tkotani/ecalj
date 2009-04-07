      subroutine ywrm(opt,filel,icast,ifi,fmt,s,ofi,ns,nr,nc)
C- Writes complex matrix to file ifi.  File is not rewound before writing
C ----------------------------------------------------------------
Ci Inputs:
Ci   opt:   1s digit:
Ci           0 writes in ascii mode
Ci           1 writes in binary mode
Ci  filel:   reference string
Ci  icast:   0 integer
Ci           1 double precision
Ci           2 double complex with imaginary following real
Ci           3 double complex
Ci           4 double complex with imaginary following real in columns
Ci      Add 10 to indicate symmetric or hermitian
Ci    ifi:   file logical unit
Ci    fmt:   fortran format to write ascii string, e.g. '(5f12.6)'
Ci           You can use awrite format, e.g. '(%5,6g)'
Ci      s:   matrix to be printed out
Ci  ns,nr:   leading dimension of s and number of rows to write
Ci nsc,nc:   column dimension of s and number of columns to write
Cr Remarks
Cr   Binary write first record is
Cr     nr  nc  cast  optional-string-length
Cr   If a optional-string-length > 0, second record contains string
Cr   Next record is entire array, written as nr*nc elements
Cu Updates
Cu   3 Sep 99  argument list changed: offset to Im s is passed.
Cu  17 Mar 99  ywrm can write integer array
Cu  19 Jul 99  implement icast=4 format for ascii write
C ----------------------------------------------------------------
C     implicit none
      integer opt,icast,ofi,ns,nr,nc,ifi
      character*(*) filel,fmt
      double precision s(ns,nc)
      character outs*80,lfmt*40,sout*500
      integer i,j,nw,iw(100),lf,ic,jc,k,ip,ir,imm,lbin
      logical a2bin,ltmp,gettok

      lbin = mod(opt,10)

C --- Binary write ---
      if (lbin .gt. 0) then
        if (filel .eq. ' ') then
          write(ifi) nr,nc,icast,0
        else
          write(ifi) nr,nc,icast,len(filel)
          write(ifi) filel(1:len(filel))
        endif
        if (nr .ne. ns) call rx('ywrm not ready for nr.ne.ns')
        if (mod(icast,10) .eq. 0)  then
          call ipdump(s,nr*nc,-ifi)
        else
          k = 1
          if (icast .gt. 1) k = 2
          if (icast .eq. 2 .and. ofi .ne. nr*nc .or.
     .        icast .eq. 4 .and. ofi .ne. nr)
     .      call rx('ywrm not ready for this ofi')
          call dpdump(s,nr*nc*k,-ifi)
        endif
        return
      endif
      if (mod(icast,10) .eq. 0)
     .  call rx('ywrm not ready for ascii write of integer')
      outs = ' '
      if (icast .eq. 1)  outs = ' real l="' //filel
      if (icast .eq. 11) outs = ' symm l="' //filel
      if (icast .ge. 2 .and. icast .le. 4) outs = ' complex l="' //filel
      if (icast .ge.12 .and. icast .le.14) outs = ' herm l="' //filel
      call word(outs,2,i,j)
      if (outs(i:) .eq. 'l="') then
        outs(i:) = ' '
      else
        call skpblb(outs,len(outs),i)
        outs(i+2:i+2) = '"'
      endif
      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)

C --- Case awrite format ---
      if (fmt(2:2) .eq. '%') then
        if (icast .eq. 3 .or. icast .eq. 13 .or.
     .      icast .eq. 4 .or. icast .eq. 14)
     .   call rx('ywrm cannot handle %% fmt for complex*16')

C   ... copy to lfmt to get around SGI compiler bug
        lfmt = fmt
        lf = len(lfmt)
C   ... Find end of first expression
        j = 2
   30   ltmp = gettok(lfmt,j,' ',lf)
        if (ltmp) goto 30
        nw = 1
        if (j .gt. 2) then
          i = 2
          call rxx(.not. a2bin(lfmt,nw,2,0,' ',i,j-1),
     .      'ywrm failed to parse format %'//lfmt(2:j+2))
        endif
C   ... Loop through all rows to find col-dependent size
        nw = min(nw,100)
        call iinit(iw,nw)
        call skpblb(lfmt,lf,i)
        outs = '%'//lfmt(j+1:i)
        do  32  i = 1, nr
          call prmwid(s(i,1),nw,outs,iw,ns,nc)
          if (mod(icast,10) .gt. 1)
     .    call prmwid(s(i+ofi,1),nw,outs,iw,ns,nc)
   32   continue
*       call awrit2('%n:1i',' ',80,6,nw,iw)

C   ... Write out matrix
        imm = 0
   44   continue
        do  40  ir = 1, nr
        do  40  ic = 1, nc, nw
          jc = min(ic+nw-1,nc)
          k = 0
          sout = ' '
          ip = 0
          do  42  i = ic, jc
            k = k+1
            call awrit1(outs,lfmt,len(lfmt),0,s(ir+imm*ofi,i))
            call skpblb(lfmt,iw(k),j)
            call awrit2('%np'//outs,sout,len(sout),0,
     .        ip+iw(k)-j,s(ir+imm*ofi,i))
            ip = ip + iw(k)
   42     continue
          call awrit0('%a',sout,-len(sout),-ifi)
   40   continue
        if (imm .eq. 0 .and. mod(icast,10) .gt. 1) then
          imm = 1
          write(ifi,'(1x)')
          goto 44
        endif
        return
      endif

      if (mod(icast,10) .eq. 3) then
        do  50  i = 1, nr
   50   write(ifi,fmt) (s(2*i-1,2*j-1), j=1,nc)
        write(ifi,'(1x)')
        do  52  i = 1, nr
   52   write(ifi,fmt) (s(2*i,2*j-1), j=1,nc)

      elseif (mod(icast,10) .eq. 4) then
        do  60  i = 1, nr
   60   write(ifi,fmt) (s(i,2*j-1), j=1,nc)
        write(ifi,'(1x)')
        do  62  i = 1, nr
   62   write(ifi,fmt) (s(i,2*j), j=1,nc)

      else
        do  10  i = 1, nr
   10   write(ifi,fmt) (s(i,j), j=1,nc)
        if (mod(icast,10) .gt. 1) then
          write(ifi,'(1x)')
          do  20  i = 1, nr
   20     write(ifi,fmt) (s(i+ofi,j), j=1,nc)
        endif
      endif
      end

      subroutine prmwid(s,nw,fmt,iw,ns,nc)
C- Determine the widths needed to aligning in columns
C     implicit none
      integer nw,ns,nc,iw(nw)
      character*(*) fmt
      double precision s(ns,nc)
      integer ic,jc,i,k,j1,j2
      character*50 sout

      do  10  ic = 1, nc, nw
        jc = min(ic+nw-1,nc)
        k = 0
        do  20  i = ic, jc
          sout = ' '
          call awrit1(fmt,sout,len(sout),0,s(1,i))
          k = k+1
          call word(sout,1,j1,j2)
          iw(k) = max(iw(k),j2-j1+2)
*         print *, k,j1,j2,iw(k),' ',sout(j1:j2)
   20   continue
   10 continue
      end

