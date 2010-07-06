      subroutine fmain
      implicit none
      character*80 outstr
      integer ip
      double precision xx

      call pshprt(121)
      call bin2a0(10)

C      This points out an SGI compiler bug
C      xx = 999.99999999999977d0
C      outstr = ' '
C      ip = 0
C      call bin2a('g4 ',0,0,xx,4,0,60,outstr,ip)

C --- Test fortran f format, floating point ---
      print *, 'fortran f format'
      ip = 0
      outstr = ' '
      call bin2a('(f12.6)',1,3,1.23d5,4,0,60,outstr,ip)
      call bin2a('(f12.6)',1,3,1.23d1,4,0,60,outstr,ip)
      call bin2a('(f12.6)',1,0,1.23456d-2,4,0,60,outstr,ip)
      call bin2a('(f12.6)',1,3,1.23d-6,4,0,60,outstr,ip)
      call bin2a('(f12.6)',1,3,1.23d-7,4,0,60,outstr,ip)
      call bin2a('(f12.6)',1,0,1.23d-7,4,0,60,outstr,ip)
      call bin2a('(e18.6)',1,0,1.23d-7,4,0,60,outstr,ip)

C --- Test d format, floating point, big numbers ---
      print *, 'd format'
      ip = 0
      outstr = ' '
      call bin2a('d',1,0,1.23d7,4,0,80,outstr,ip)
      call bin2a('d',1,3,1.23d7,4,0,80,outstr,ip)
      call bin2a('d2',1,3,1.23d7,4,0,80,outstr,ip)
C ... note ndec=3 properly -> ndec=0
      call bin2a('d2:1',1,3,1.23d7,4,0,80,outstr,ip)
      call bin2a('d2:1',1,3,1.2499d7,4,0,80,outstr,ip)
      call bin2a('d2:1',1,3,1.2501d7,4,0,80,outstr,ip)

C --- Test d format, floating point, middle size numbers ---
      ip = 0
      outstr = ' '
      call bin2a('d',1,0,1.2349d0,4,0,80,outstr,ip)
      call bin2a('d2',1,3,1.234951d0,4,0,80,outstr,ip)
      call bin2a('d3',1,3,1.234951d0,4,0,80,outstr,ip)
      call bin2a('d4',1,4,1.234951d0,4,0,80,outstr,ip)
C ... note ndec=4 properly -> ndec=precsn
      call bin2a('d2:1',1,5,1.234951d0,4,0,80,outstr,ip)
      call bin2a('d3:1',1,5,1.234951d0,4,0,80,outstr,ip)
      call bin2a('d4:1',1,5,1.2349510d0,4,0,80,outstr,ip)
      call bin2a('d5:1',1,5,1.2349510d0,4,0,80,outstr,ip)
      call bin2a('d6:1',1,5,1.2349510d0,4,0,80,outstr,ip)

C --- Test d format, floating point, small numbers ---
      ip = 0
      outstr = ' '
      call bin2a('d',1,3,1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('d2:10',1,3,1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('d10',1,3,1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('d2:1',1,0,1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('d3:11',1,0,1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('d4:1',1,0,1.2349501d-6,4,0,80,outstr,ip)

C --- Test D format, floating point, big numbers ---
      print *, 'D format'
      ip = 0
      outstr = ' '
      call bin2a('D12',1,1,1.23d7,4,0,80,outstr,ip)
      call bin2a('D12',1,2,1.23d7,4,0,80,outstr,ip)
      call bin2a('D12',1,3,1.23d7,4,0,80,outstr,ip)
      call bin2a('D12',1,4,1.23d7,4,0,80,outstr,ip)

C --- Test D format, floating point, middle size numbers ---
      ip = 0
      outstr = ' '
      call bin2a('D8',1,1,1.234951d0,4,0,80,outstr,ip)
      call bin2a('D8',1,2,1.234951d0,4,0,80,outstr,ip)
      call bin2a('D8',1,3,1.234951d0,4,0,80,outstr,ip)
      call bin2a('D8',1,4,1.234951d0,4,0,80,outstr,ip)
      call bin2a('D8',1,5,1.234951d0,4,0,80,outstr,ip)
      call bin2a('D8',1,6,1.234951d0,4,0,80,outstr,ip)
      call bin2a('D8',1,7,1.2349510d0,4,0,80,outstr,ip)

C --- Test D format, floating point, small numbers ---
      ip = 0
      outstr = ' '
      call bin2a('D12',1,7,1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('D12:1',1,8,1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('D12',1,9,1.2349501d-6,4,0,80,outstr,ip)

C --- Test D format, floating point, small negative numbers ---
      ip = 0
      outstr = ' '
      call bin2a('D12',1,7,-1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('D12:1',1,8,-1.2349501d-6,4,0,80,outstr,ip)
      call bin2a('D12',1,9,-1.2349501d-6,4,0,80,outstr,ip)

C --- Test e format, floating point, big numbers ---
      print *, 'e format'
      ip = 0
      outstr = ' '
      call bin2a('e',1,0,1.23d7,4,0,80,outstr,ip)
      call bin2a('e',1,4,1.23d7,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.23d7,4,0,80,outstr,ip)
      call bin2a('e3',1,3,1.23d7,4,0,80,outstr,ip)
      call bin2a('e4',1,3,1.23d7,4,0,80,outstr,ip)
      call bin2a('e5',1,5,1.2349510d7,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.2499d7,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.2501d7,4,0,80,outstr,ip)
      print *, '---'
      ip = 0
      outstr = ' '
      call bin2a('e-2:0',1,0,1.2349501d5,4,0,80,outstr,ip)
      call bin2a('e1:0',1,0,1.2349501d5,4,0,80,outstr,ip)
      call bin2a('e1:0',1,1,1.2349501d5,4,0,80,outstr,ip)
      call bin2a('e5:0',1,0,1.2349501d5,4,0,80,outstr,ip)
      call bin2a('e5:0',1,4,1.2349501d5,4,0,80,outstr,ip)

C --- Test e format, floating point, middle-size numbers ---
      ip = 0
      outstr = ' '
      call bin2a('e',1,0,1.23d1,4,0,80,outstr,ip)
      call bin2a('e',1,4,1.23d1,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.23d1,4,0,80,outstr,ip)
      call bin2a('e3',1,3,1.23d1,4,0,80,outstr,ip)
      call bin2a('e4',1,3,1.23d1,4,0,80,outstr,ip)
      call bin2a('e5',1,5,1.2349510d1,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.2499d1,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.2501d1,4,0,80,outstr,ip)
      call bin2a('e',1,0,1.23d0,4,0,80,outstr,ip)
      print *, '---'
      ip = 0
      outstr = ' '
      call bin2a('e-2:0',1,0,1.2349501d1,4,0,80,outstr,ip)
      call bin2a('e1:0',1,0,1.2349501d1,4,0,80,outstr,ip)
      call bin2a('e3:0',1,3,1.2349501d1,4,0,80,outstr,ip)
      call bin2a('e4:0',1,4,1.2349501d1,4,0,80,outstr,ip)
      call bin2a('e5:0',1,5,1.2349501d1,4,0,80,outstr,ip)
      call bin2a('e5:0',1,0,1.2349501d1,4,0,80,outstr,ip)

C --- Test e format, floating point, small numbers ---
      ip = 0
      outstr = ' '
      call bin2a('e',1,0,1.23d-5,4,0,80,outstr,ip)
      call bin2a('e',1,4,1.23d-5,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.23d-5,4,0,80,outstr,ip)
      call bin2a('e3',1,3,1.23d-5,4,0,80,outstr,ip)
      call bin2a('e4',1,3,1.23d-5,4,0,80,outstr,ip)
      call bin2a('e5',1,5,1.2349510d-5,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.249d-5,4,0,80,outstr,ip)
      call bin2a('e2',1,3,1.2501d-5,4,0,80,outstr,ip)
      print *, '---'
      ip = 0
      outstr = ' '
      call bin2a('e1:0',1,0,1.2349501d-5,4,0,80,outstr,ip)
      call bin2a('e5:0',1,0,1.2349501d-5,4,0,80,outstr,ip)
      call bin2a('e6:0',1,0,1.2349501d-5,4,0,80,outstr,ip)
      call bin2a('e7:0',1,0,1.2349501d-5,4,0,80,outstr,ip)
      call bin2a('e8:0',1,0,1.2349501d-5,4,0,80,outstr,ip)
      call bin2a('e9:0',1,0,-1.2349501d-5,4,0,80,outstr,ip)

C --- Test g format, floating point, very big numbers ---
      print *, 'g format'
      ip = 0
      outstr = ' '
      call bin2a('g',1,0,1.23d197,4,0,80,outstr,ip)
      call bin2a('g',1,4,1.23d97,4,0,80,outstr,ip)
      call bin2a('g1',1,2,1.2499d97,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2499d97,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2501d97,4,0,80,outstr,ip)
C ... note ndec=3 properly -> ndec=0
      call bin2a('g2:0',1,0,1.23d97,4,0,80,outstr,ip)
      call bin2a('g2:0',1,0,1.249d97,4,0,80,outstr,ip)
      call bin2a('g1:0',1,1,1.250d97,4,0,80,outstr,ip)

C --- Test g format, floating point, big numbers ---
      print *
      ip = 0
      outstr = ' '
      call bin2a('g',1,0,1.23d7,4,0,80,outstr,ip)
      call bin2a('g',1,4,1.23d7,4,0,80,outstr,ip)
      call bin2a('g1',1,2,1.2499d7,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2499d7,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2501d7,4,0,80,outstr,ip)
C ... note ndec=3 properly -> ndec=0
      call bin2a('g2:0',1,0,1.23d7,4,0,80,outstr,ip)
      call bin2a('g2:0',1,0,1.249d7,4,0,80,outstr,ip)
      call bin2a('g1:0',1,1,1.250d7,4,0,80,outstr,ip)

C --- Test g format, floating point, middle numbers ---
      print *
      ip = 0
      outstr = ' '
      call bin2a('g',1,0,1.23d0,4,0,80,outstr,ip)
      call bin2a('g',1,4,1.23d0,4,0,80,outstr,ip)
      call bin2a('g1',1,2,1.2499d0,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2499d0,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2501d0,4,0,80,outstr,ip)
      call bin2a('g2:0',1,0,1.23d0,4,0,80,outstr,ip)
      call bin2a('g2:0',1,0,1.249d0,4,0,80,outstr,ip)
      call bin2a('g1:0',1,1,1.2501d0,4,0,80,outstr,ip)
      call bin2a('g0:0',1,0,1.2501d0,4,0,80,outstr,ip)

C --- Test g format, floating point, small numbers ---
      print *
      ip = 0
      outstr = ' '
      call bin2a('g',1,0,4.567d-99,4,0,80,outstr,ip)
      call bin2a('g',1,0,1.23d-5,4,0,80,outstr,ip)
      call bin2a('g',1,4,1.23d-4,4,0,80,outstr,ip)
      call bin2a('g1',1,2,1.2499d-5,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2499d-5,4,0,80,outstr,ip)
      call bin2a('g2',1,2,1.2501d-5,4,0,80,outstr,ip)
      call bin2a('g2:10',1,0,1.23d-5,4,0,80,outstr,ip)
      call bin2a('g3:10',1,3,1.24996d-5,4,0,80,outstr,ip)
      call bin2a('g4:10',1,4,1.24996d-5,4,0,80,outstr,ip)
      call bin2a('g4:10',1,0,1.24996d-5,4,0,80,outstr,ip)

C --- Test g format, floating point, near or at 0 ---
      print *
      ip = 0
      outstr = ' '
      call bin2a('g',1,0,4.567d-229,4,0,80,outstr,ip)
      call bin2a('g',1,0,0d0,4,0,80,outstr,ip)
      call bin2a('g2',1,0,0d0,4,0,80,outstr,ip)
      call bin2a('g2',1,2,0d0,4,0,80,outstr,ip)
      call bin2a('g2:1',1,2,0d0,4,0,80,outstr,ip)
      call bin2a('g2:0',1,2,0d0,4,0,80,outstr,ip)

C --- Test 'c', 'i' and 'l' representations ---
      print *
      ip = 0
      outstr = ' '
      call bin2a('The',1,0,0,1,0,14,outstr,ip)
      call bin2a(' ',2,0,.true.,0,0,14,outstr,ip)
      call bin2a('rue abc murders',0,0,0,1,0,14,outstr,ip)
      call bin2a(' ',2,0,.true.,0,0,14,outstr,ip)
      call bin2a(' ',1,0,123,2,0,20,outstr,ip)
      call bin2a(' ',1,0,456,2,0,20,outstr,ip)

C --- Test F format, floating point, big numbers ---
      ip = 0
      outstr = ' '
      print *, 'F format, field width = 3'
C     call bin2a0(00)
      call bin2a('F3',1,0,1.23456789d2,4,0,60,outstr,ip)
      call bin2a('F3',1,0,1.23456789d1,4,0,60,outstr,ip)
      call bin2a('F3',1,0,1.23456789d0,4,0,60,outstr,ip)
      call bin2a('F3',1,0,1.23456789d-1,4,0,60,outstr,ip)
      call bin2a('F3',1,0,1.23456789d-2,4,0,60,outstr,ip)
      call bin2a('F3',1,0,1.23456789d-3,4,0,60,outstr,ip)
      print '(1x,a)', outstr(1:ip)

      ip = 0
      outstr = ' '
      print *, 'F format, field width = 7'
      call bin2a('F7',1,0,1.23456789d6,4,0,80,outstr,ip)
      call bin2a('F7',1,0,1.23456789d5,4,0,80,outstr,ip)
      call bin2a('F7',1,0,1.23456789d2,4,0,80,outstr,ip)
      call bin2a('F7',1,0,1.23456789d1,4,0,80,outstr,ip)
      call bin2a('F7',1,0,1.23456789d0,4,0,80,outstr,ip)
      call bin2a('F7',1,0,1.23456789d-1,4,0,80,outstr,ip)
      call bin2a('F7',1,0,1.23456789d-3,4,0,80,outstr,ip)
      call bin2a('F7',1,0,1.23456789d-4,4,0,80,outstr,ip)
      print '(1x,a)', outstr(1:ip)

      ip = 0
      outstr = ' '
      print *, 'F format, field width = 12'
      call bin2a('F12',1,0,1.23456789d11,4,0,80,outstr,ip)
      call bin2a('F12',1,0,1.23456789d11,4,0,80,outstr,ip)
      call bin2a('F12',1,0,1.23456789d10,4,0,80,outstr,ip)
      call bin2a('F12',1,0,1.23456789d2,4,0,80,outstr,ip)
      call bin2a('F12',1,0,1.23456789d-3,4,0,80,outstr,ip)
      call bin2a('F12',1,0,1.23456789d-4,4,0,80,outstr,ip)
      print '(1x,a)', outstr(1:ip)

      print *, 'F format, special cases'
      ip = 0
      outstr = ' '
      call bin2a('F6',1,0,-1.4926251754d-3,4,0,len(outstr),outstr,ip)
      call bin2a('F5',1,0,-1.4926251754d-3,4,0,len(outstr),outstr,ip)
      call bin2a('F4',1,0,-1.4926251754d-3,4,0,len(outstr),outstr,ip)
      call bin2a('F3',1,0,-1.4926251754d-3,4,0,len(outstr),outstr,ip)
      call bin2a('F2',1,0,-1.4926251754d-3,4,0,len(outstr),outstr,ip)
      call bin2a('F6',1,0,-1.9971712021d-4,4,0,len(outstr),outstr,ip)

      end
