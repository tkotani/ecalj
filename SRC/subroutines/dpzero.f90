      subroutine cinit(array,leng)
C- Initializes complex array to zero
      integer leng
      real array(2*leng)
      call ainit(array,leng+leng)
      end
      subroutine zinit(array,leng)
C- Initializes complex array to zero
      integer leng
      double precision array(2*leng)
      call dpzero(array,leng+leng)
      end
      subroutine iinit(array,leng)
C- Initializes integer array to zero
      integer leng,i
      integer array(leng)
      do   i=1, leng
         array(i) = 0
      enddo
      end
      subroutine ainit(array,leng)
C- Initializes real array to zero
      integer leng,i
      real array(leng)
      do   i=1, leng
         array(i) = 0
      enddo   
      end
      subroutine dpzero(array,leng)
C- Initializes double precision array to zero
      integer leng,i
      double precision array(leng)
      do   i=1, leng
         array(i) = 0
      enddo   
      end
      real function rval(array,index)
      integer index
C- Returns the real value of ARRAY(INDEX)
      real array(index)
      rval = array(index)
      end
      double precision function dval(array,index)
      integer index
C- Returns the double precision value of ARRAY(INDEX)
      double precision array(index)
      dval = array(index)
      end
      integer function ival(array,index)
C- Returns the integer value of ARRAY(INDEX)
      integer index
      integer array(index)
      ival = array(index)
      end
      integer function ival2(array,nda,i1,i2)
C- Returns the integer value of ARRAY(i1,i2)
C     implicit none
      integer nda,i1,i2,array(nda,1)
      ival2 = array(i1,i2)
      end
      logical function logval(array,index)
C- Returns the integer value of ARRAY(INDEX)
      integer index
      logical array(index)
      logval = array(index)
      end
      complex function cval(array,index)
C- Returns the complex value of ARRAY(INDEX)
      integer index
      complex array(index)
      cval = array(index)
      end
      subroutine dvset(array,i1,i2,val)
C- Sets some elements of double precision array to value
      integer i1,i2
      double precision array(i2),val
      integer i
      do   i = i1, i2
         array(i) = val
      enddo   
      end
      subroutine ivset(array,i1,i2,val)
C- Sets some elements of integer array to value
      integer i1,i2,array(1),val,i
      do   i = i1, i2
         array(i) = val
      enddo   
      end
      subroutine lvset(array,i1,i2,val)
C- Sets some elements of logical array to value
      integer i1,i2,i
      logical array(1),val
      do  i = i1, i2
         array(i) = val
      enddo   
      end
C$$$      subroutine redfrr(oname,leng)
C$$$C- Release to pointer oname, reallocate double oname of length leng
C$$$C     implicit none
C$$$      integer oname,leng
C$$$      call rlse(oname)
C$$$      call defrr(oname,leng)
C$$$      end
C$$$      subroutine redfi(oname,leng)
C$$$C- Release to pointer oname, reallocate double oname of length leng
C$$$C     implicit none
C$$$      integer oname,leng
C$$$      call rlse(oname)
C$$$      call defi(oname,leng)
C$$$      end

