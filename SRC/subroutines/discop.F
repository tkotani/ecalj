      subroutine discop(afrom,ito,nel,nf,nt,mode)
C- Copies double precision array to integer array
Ci   nf,nt    offsets for afrom and ito
Ci   nel      number of elements to copy
Ci   mode     1,  use ito = nint(afrom)
Co   ito      copy ifrom into this array
C     implicit none
      integer nf,nt,i,iadd,ntop,nel,ito(1),mode
      double precision afrom(1)

      iadd = nt-nf
      ntop = nf+nel-1

      if (mode .eq. 1) then
        do    i = nf, ntop
           ito(i+iadd) = nint(afrom(i))
        enddo   
        return
      else
        do   i = nf, ntop
           ito(i+iadd) = afrom(i)
        enddo   
        return
      endif
      end

      subroutine idscop(ifrom,to,nel,nf,nt)
C- Copies integer array to double precision array
Ci   ifrom
Ci   nf,nt    offsets for ifrom and to (1 for 1st array elt)
Ci   nel      number of elements to copy
Co   to       copy ifrom into this array
C     implicit none
      integer nf,nt,i,iadd,ntop,nel,ifrom(1)
      double precision to(1)

      iadd = nt-nf
      ntop = nf+nel-1
      do  i = nf, ntop
         to(i+iadd) = ifrom(i)
      enddo   
      end

C      subroutine fmain
CC      implicit none
C      integer ito(10)
C      double precision from(10)
C      integer i
C
C      do  10  i = 1, 10
C      ito(i) = 0
C   10 from(i) = i + 2d0/i
C
C      call discop(from,ito,4,3,5,0)
C      print 332, from
C      print 333, ito
C      call discop(from,ito,4,3,5,1)
C      print 333, ito
C  332 format(10f8.3)
C  333 format(10i4)
C
C      call idscop(ito,from,4,3,5)
C      print 332, from
C
C      end

