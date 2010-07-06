      integer function ipsw(defsw,sw,ndef,isw)
C- Unpack a set of integers switches from a string
C ----------------------------------------------------------------------
Ci Inputs
Ci   sw   : string containing 1-digit integers
Ci   defsw: string of defaults.  For sw which takes value ndef, defsw is
Ci          used instead of sw.  sw and defsw may be identical.
Co Outputs
Co   isw: vector of integers
Cr Remarks:
Cr   The length of defsw determines the number of integers converted.
C ----------------------------------------------------------------------
      implicit none
      integer ndef,isw(7)
      character*(*) defsw,sw
      integer id,jd,is,js,i,n1,n2

      call word(defsw,1,id,jd)
      n1 = jd-id+1
      call word(sw,1,is,js)
      n2 = js-is+1
      read(defsw(id:jd),333) (isw(i), i=1,n1)
      read(sw(is:js),333)    (isw(i), i=1,n2)
      do  10  i = 1, min(n1,n2)
   10 if (isw(i) .eq. ndef) read(defsw(id+i-1:id+i-1),333) isw(i)
  345 format(999i2)
  333 format(999i1)
      ipsw = max(n1,n2)
      end
C ... Test ipsw
C      subroutine fmain
C      integer upsw(30),ipsw,n,i
C
C      call iinit(upsw,30)
C      n = ipsw('3194821','123',3,upsw)
C      print 345, (upsw(i), i=1,n)
C
C
C      call iinit(upsw,30)
C      n = ipsw('3194821','1234567890',3,upsw)
C      print 345, (upsw(i), i=1,n)
C
C
C  345 format(999i2)
C      end

