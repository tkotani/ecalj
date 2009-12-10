      integer function mkdlst(strn,fuzz,npmx,xp)
C- Resolve list (ascii string) into a vector of double-precision numbers
C ----------------------------------------------------------------
Ci Inputs
Ci   strn  :string holding list of numbers see Remarks for syntax
Ci Inputs
Ci   fuzz  :uncertainty in the precision.  In a number sequence,
Ci         :the upper bound is known only to precision fuzz.
Ci         :Upper limit is replaced by upper limit + fuzz.
Ci         :fuzz<0 => use |fuzz|*|max element from strn|
Ci   npmx  :npmx>0  => maximum number of points allowed.
Ci          npmx<=0 => mkdlst returns np without filling xp
Co Outputs:
Co   xp(1..) the vector of numbers (read only if npmx>0)
Co   mkdlst: number of points read, or would try to read if npmx>0.
Co           If failed to parse ascii input, returns -1.
Co           If np>npmx and npmx>0, returns -npmx
Cr Remarks
Cr   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
Cr   low:high:step
Cr   low, high, and step are real expressions specifying the sequence
Cr     low, low+step, low+2*step, ... high.
Cr   If :step is missing, the step size defaults to 1.  If also :high
Cr   is missing,  the sequence reduces to a single integer. Thus,
Cr     '5+1'       becomes a single number, 6.
Cr     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
Cr     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
Cr   Sequences may be strung together separated by commas, eg
Cr     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
Cu Updates
Cu   25 Aug 04 Added fuzz.  New argument list
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npmx
      character*(*) strn
      double precision fuzz,xp(*)
C ... Local parameters
      double precision dv(1000),d1mach,dx,tfuzz
      integer it(1000),a2vec,ip,i,j,k,n,iprint,i1mach,np

      ip = 0
      call skipbl(strn,len(strn),ip)
      k = a2vec(strn,len(strn),ip,4,',: ',3,3,30,it,dv)
      tfuzz = fuzz
      if (tfuzz .lt. 0) then
        call idmxmn(k,dv,1,i,j)
        tfuzz = max(dabs(dv(i)),dabs(dv(j))) * abs(fuzz)
      endif
      mkdlst = -1
      if (k .lt. 1) return
      mkdlst = -npmx
C ... loop over all iv
      np = 0
      i = 0
   14 i = i+1
C ... Case dv => a single number
      if (it(i) .ne. 2) then
        np = np+1
        if (npmx .gt. 0) then
          if (np .gt. npmx) return
          xp(np) = dv(i)
        endif
C ... Case dv => n1:n2[:n3]
      else
        dx = 1
        if (it(i+1) .eq. 2 .and. dv(i+2) .ne. 0) dx = dv(i+2)
        n = int((dv(i+1)+tfuzz-dv(i))/dx)
        do  12  j = 0, n
          np = np+1
          if (npmx .gt. 0) then
            if (np .gt. npmx) return
            xp(np) = dv(i) + j*dx
          endif
   12   continue
        i = i+1
        if (it(i) .eq. 2) i = i+1
      endif
      if (i .lt. k) goto 14


C --- Entry for np = npmax
C  20 continue
      mkdlst = np
      if (iprint() .gt. 100 .and. npmx .gt. 0) call
     .  awrit3(' MKDLST: %i points:%n:1g',' ',100,i1mach(2),np,np,xp)
        
      end
C --- Test of mkdlst ---
C      subroutine fmain
CC      implicit none
C      integer np,npmx,mkdlst
C      parameter (npmx=20)
C      double precision xp(npmx)
C
C      call pshpr(101)
C      np = mkdlst('9,3:5:.2 ',npmx,xp)
C      print *, np
C      end

