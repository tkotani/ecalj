      subroutine fmain
      implicit none
      character *40 sout
      integer nout
C ... global change 10*0 -> 10*1 to check suppress leading zero ..

  333 format(3i3,1x,a,a,a)
C --- test 'f' representation ---
C ... Standard case, test absolute precision with rounding
      call pretty('+9.123456',0,0,4,0+10*0,sout,nout)
      print 333, 0,4,0, '+9.123456 ->|', sout(1:nout), '|'
C ... Standard case, test absolute precision with no rounding
      call pretty('.12321',0,0,4,0+10*0,sout,nout)
      print 333, 0,4,0, '.12321 ->|', sout(1:nout), '|'
C ... Test leading 0, leading '-', roundoff
      call pretty('-.12996',0,0,4,0+10*0,sout,nout)
      print 333, 0,4,0, '-.12996 ->|', sout(1:nout), '|'
C ... ditto, with ndec > 0
      call pretty('-.12996',0,3,4,0+10*0,sout,nout)
      print 333, 3,4,0, '-.12996 ->|', sout(1:nout), '|'
C ... ditto, with ndec > nprec
      call pretty('-.12996',0,5,4,0+10*0,sout,nout)
      print 333, 5,4,0, '-.12996 ->|', sout(1:nout), '|'
C ... ditto, with ndec > number of available digits
      call pretty('-.13',0,5,9,0+10*0,sout,nout)
      print 333, 5,9,0, '-.13 ->|', sout(1:nout), '|'
C ... A funny case
      call pretty('-9.9996',0,0,3,0+10*0,sout,nout)
      print 333, 0,3,0, '-9.9996 ->|', sout(1:nout), '|'
C ... Standard case, test relative precision with rounding
      call pretty('+9.123456',0,0,4,1+10*0,sout,nout)
      print 333, 0,4,1, '+9.123456 ->|', sout(1:nout), '|'
C ... Test relative precision with lots of rounding
      call pretty('-0.12996',0,3,4,1+10*0,sout,nout)
      print 333, 3,4,1, '-.12996 ->|', sout(1:nout), '|'
C ... Test large number with precision to left of decimal
      call pretty('12645678.9098765',0,2,2,1+10*0,sout,nout)
      print 333, 2,2,1, '12345678.9098765 ->|', sout(1:nout), '|'
C ... Test very large number with precision to left of decimal
      call pretty('126456789098765',0,2,12,1+10*0,sout,nout)
      print 333, 2,12,1, '-123456789098765 ->|', sout(1:nout), '|'
C ... Standard case, test leading blank
      call pretty('-9.123456',3,0,4,0+10*0,sout,nout)
      print 333, 0,4,0, '-9.123456 ->|', sout(1:nout), '|'
C ... case when value is 0
      call pretty('0.0',0,0,4,10,sout,nout)
      print 333, 0,4,0, '0.0 ->|', sout(1:nout), '|'
C ... case when value is small but < 0
      call pretty('-0.00001',0,0,4,10,sout,nout)
      print 333, 0,4,0, '-0.00001 ->|', sout(1:nout), '|'
C ... A small number, relative precision
      call pretty('0.0000012345',1,2,2,1,sout,nout)
      print 333, 2,2,1, '0.0000012345 ->|', sout(1:nout), '|'
C ... A small number, relative precision, rounded
      call pretty('0.0000012345',1,4,4,1,sout,nout)
      print 333, 4,4,1, '0.0000012345 ->|', sout(1:nout), '|'

C --- test 'E' representation ---
C ... A standard case
      call pretty('-0.123000E+03',0,0,4,0+10*0,sout,nout)
      print 333, 0,4,0, '-0.123000E+03 ->|', sout(1:nout), '|'
C ... test precision & truncation, absolute precision
      call pretty('-0.45678E-002',0,6,6,0+10*0,sout,nout)
      print 333, 6,6,0, '-0.45678E-002 ->|', sout(1:nout), '|'
C ... test supression of trailing zeros, absolute precision
      call pretty('-0.456780000E-002',0,8,9,0+10*0,sout,nout)
      print 333, 8,9,0, '-0.456780000E-002 ->|', sout(1:nout), '|'
C ... test addition of trailing zeros, absolute precision
      call pretty('-0.45678E-002',0,8,9,0+10*0,sout,nout)
      print 333, 8,9,0, '-0.45678E-002 ->|', sout(1:nout), '|'
C ... test precision & truncation, absolute precision
      call pretty('-0.45678E-002',0,3,3,1+10*0,sout,nout)
      print 333, 3,3,1, '-0.45678E-002 ->|', sout(1:nout), '|'
C ... test supression of trailing zeros, absolute precision
      call pretty('-0.456780000E-002',0,6,7,1+10*0,sout,nout)
      print 333, 6,7,1, '-0.456780000E-002 ->|', sout(1:nout), '|'
C ... test addition of trailing zeros, absolute precision
      call pretty('-0.45678E-002',0,6,7,1+10*0,sout,nout)
      print 333, 6,7,1, '-0.45678E-002 ->|', sout(1:nout), '|'
C ... A funny case
      call pretty('+1.996e3',0,0,3,1+10*0,sout,nout)
      print 333, 0,3,1, '+1.996e3 ->|', sout(1:nout), '|'
C ... Another funny case
      call pretty('+1.E+07',1,1,1,1,sout,nout)
      print 333, 1,1,1, '+1.E+07 ->|', sout(1:nout), '|'
C ... case when value is small but < 0
      call pretty('-0.00001d-1',0,0,4,10,sout,nout)
      print 333, 0,4,0, '-0.00001d-1 ->|', sout(1:nout), '|'
C ... Very big and small numbers
      call pretty('1.2300000+197',0,6,7,1+10*0,sout,nout)
      print 333, 6,7,1, '1.2300000+197 ->|', sout(1:nout), '|'
      call pretty('1.2300000-197',0,6,7,1+10*0,sout,nout)
      print 333, 6,7,1, '1.2300000-197 ->|', sout(1:nout), '|'
      call pretty('1.2499+197',0,3,3,1+10*0,sout,nout)
      print 333, 3,3,1, '1.2499+197 ->|', sout(1:nout), '|'
      call pretty('1.2499-197',0,3,3,1+10*0,sout,nout)
      print 333, 3,3,1, '1.2499-197 ->|', sout(1:nout), '|'

      call cexit(1,1)
      end
