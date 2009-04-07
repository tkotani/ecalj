      subroutine crdbf(l1,l2,s,array,nrow,ncol,cast,ifi)
C- Conditional read of an array in a binary format, in standard form
C ----------------------------------------------------------------------
Ci Inputs
Ci   l1    :F array is not sought; read past record
Ci         :T array is sought
Ci   l2    :F array is not on disk
Ci         :T array is on disk
Ci   s     :name, for used in case of error
Ci   nrow  :number of rows to read
Ci   ncol  :number of columns to read
Ci   cast  :0 integer array
Ci         :1 real array
Ci         :2 complex array
Ci   ifi   :file handle
Co Outputs
Ci   array :array is read in, provided the following criteria are met:
Ci          1. input nrow = file nrow
Ci          2. input ncol = file ncol
Ci          3. input cast matches file value
Ci          4. file is on disk
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical l1,l2,lddump,lidump
      character *(*) s
      integer cast,nrow,ncol,ifi,array(nrow,ncol)
C ... Local parameters
      integer n,nrowr,ncolr,castr
      character*20 outs
C ... Heap
      integer w(1)
      common /w/ w

      if (l1 .and. .not. l2) goto 99
      if (.not. l2) return
      if (.not. l1) then
        read(ifi)
        read(ifi)
        return
      endif
      read(ifi) nrowr,ncolr,castr
      outs = ', array '//s
      call skpblb(outs,len(outs),n)
      n = n+1
      call isanrg(nrow,nrowr,nrowr,'crdfb: ','nrow'//outs(1:n),.true.)
      call isanrg(ncol,ncolr,ncolr,'crdfb: ','ncol'//outs(1:n),.true.)
      call isanrg(cast,castr,castr,'crdfb: ','cast'//outs(1:n),.true.)

      n = nrowr*ncolr

      if (mod(cast,10) .eq. 2) then
        if (.not. lddump(array,2*n,ifi)) goto 99
      elseif (mod(cast,10) .eq. 1) then
        if (.not. lddump(array,n,ifi)) goto 99
      else
        if (.not. lidump(array,n,ifi)) goto 99
      endif
      return

   99 call rxs('CRDFB: file mismatch, missing array ',s)

      end


