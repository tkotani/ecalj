      subroutine fmain
      implicit none
      double precision res,rev(10)
      integer ival,i1mach,a2vec,ip,nvec,nsep,iterm,n,i
      character*1 nam*20

      do  10  i = 1, 10
   10 rev(i) = i

      call addsvv('abc',10,ival)
      call addsvv('three',6,ival)
      call addsvv('five',5,ival)
      call lodsvv('abc',ival,0,1,10,rev)
      call lodsvv('three',ival,0,3,6,rev)
c      call lodsvv('five',ival,0,1,5,rev(5))
      call lodsvv(' ',3,1,1,5,rev(5))

C      print *,ival
      call numsvv(i)
      call shosvv(1,i,i1mach(2))


      do  110  i = 1, 10
  110 rev(i) = 99
      print 333, rev

      call getsvv('three',ival,0,2,6,rev)
      print 333, rev
  333 format(10f8.4)

      call getsvv('xxx',ival,0,2,7,rev)
      print *, ival

      call getsvv('abc',ival,0,2,7,rev)
      print *, ival

      print 333, rev

      do  210  i = 1, 10
  210 rev(i) = 99
      call getsvv(' ',1,1,2,7,rev)
      print 333, rev



      end
