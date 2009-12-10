C ... Test ipsw
      subroutine fmain
      integer upsw(30),ipsw,n,i

      call iinit(upsw,30)
      n = ipsw('3194821','123',3,upsw)
      print 345, (upsw(i), i=1,n)

      call iinit(upsw,30)
      n = ipsw('3194821','1234567890',3,upsw)
      print 345, (upsw(i), i=1,n)

  345 format(999i2)
      end
