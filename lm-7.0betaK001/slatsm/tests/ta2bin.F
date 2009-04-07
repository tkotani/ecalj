      subroutine fmain
      implicit none
      logical a2bin,a2bina,lx,a2d
      character instr*40
      double precision resd(1),rev(10)
      integer i,j,resi(1),i1mach,ival
      logical resl(1)
      real resr(1)

      ival = 0
      call pshpr(100)
      call dpzero(resd,1)
      call iinit(resi,1)
      do  10  i = 1, 10
   10 rev(i) = i
C#ifdef POINTER
      call addsvv('abc',10,ival)
      call lodsvv('abc',ival,0,1,10,rev)
      call addsvv('xyz',8,ival)
      call lodsvv('xyz',ival,0,1,8,rev(3))
      call lodsyv('yy',1,3d0,i)
C#endif

C --- Test a2d ---
      print '(/'' ------- test a2d -------'')'
      j = 1
      call word(' 2' ,1,j,i)
      j = j-1
      lx = a2d(' 2',.false.,j,i-1,' ',resd)
      print '(1x,''2d '',L1,f15.10)', lx, resd(1)

      j = 0
      lx = a2d('2d',.false.,j,-1,'d',resd)
      print '(1x,''2d '',L1,f15.10)', lx, resd(1)
C ... This should generate a parse error
      j = 0
      lx = a2d('    +1.23.4e5 ',.true.,j,-1,'z',resd)
C ... But this should not
      j = 0
      lx = a2d('    +1.23.4e5 ',.true.,j,8,'z',resd)
      print '(1x,''    +1.23.4e5 '',L1,f15.10,i4)', lx, resd(1), j

      j = 1
      lx = a2d('  -3+',.true.,j,-1,'z',resd)
      print '(1x,''  -3+'',L1,f15.10,i4)', lx, resd(1), j
      j = 0
      lx = a2d('    +1.234e5 ',.true.,j,-1,'z',resd)
      print '(1x,''    +1.234e5 '',L1,f15.5,i4)', lx, resd(1), j
      j = 0
      lx = a2d('    +1.234e-5 ',.true.,j,-1,'z',resd)
      print '(1x,''    +1.234e-5 '',L1,1pe20.10,i4)', lx, resd(1), j
      j = 0
C ... This generates a parse error
      lx = a2d('    +1.234e5. ',.true.,j,-1,'z',resd)
C     print '(1x,''    +1.234e5. '',L1,1pe20.10,i4)', lx, resd(1), j
      j = 0
      lx = a2d('.000000000 00000001234 ',.true.,j,-1,'z',resd)
      print '(1x,''.000000000 00000001234 '',L1,1pe20.10,i4)',
     .  lx, resd(1), j
      j = 0
      lx = a2d('.000000000 00000001234 ',.false.,j,-1,'z',resd)
      print '(1x,''.000000000 00000001234 '',L1,1pe20.10,i4)',
     .  lx, resd(1), j

      print '(/'' ------- test a2bin -------'')'
      j = 0
      lx = a2bina(' 3.0 ',resd,4,0,' ',j,-1)
      print '(1x,'' 3.0 '',L1,2(1pe20.10),i4)',lx,resd(1),resd(1)-3,j
C ... This should generate error
      j = 0
      lx = a2bina(' 3.0. ',resd,4,0,' ',j,-1)
C     print *,   ' 3.0. ', lx, resd(1), resd(1)-3, j
C ... But this should not (test jmax)
      j = 0
      lx = a2bina(' 3.0. ',resd,4,0,' ',j,3)
      print '(1x,'' 3.0. '',L1,2(1pe20.10),i4)',lx,resd(1),resd(1)-3,j
C ... This should return 2, since terminator stops at d
      j = 0
      lx = a2bina(' 2d1*2',resd,4,0,'d',j,-1)
      print '(1x,'' 2d1*2'',L1,f15.10)', lx, resd(1)
      j = 0
      lx = a2bina(' 3?1+2*3+4/5:2*(3+cos(pi)*-1) ',resd,4,0,' ',j,-1)
      print '(1x,'' 3?1+2*3+4/5:2*(3+cos(pi)*-1) '',L1,1pe20.10,i4)',
     .  lx, resd(1), j
      j = 0
      lx = a2bina(' 0?1+2*3+4/5:2*(3+cos(pi)*-1) ',resd,4,0,' ',j,-1)
      print '(1x,'' 0?1+2*3+4/5:2*(3+cos(pi)*-1) '',L1,1pe20.10,i4)',
     .  lx, resd(1), j
      j = 0
      lx = a2bina(' 3?2-1:4-1 ',resi,2,0,' ',j,-1)
      print '(1x,'' 3?2-1:4-1 '',L1,i4,i4)', lx, resi, j
      j = 0
      lx = a2bina(' 0?2-1:4-1 ',resi,2,0,' ',j,-1)
      print '(1x,'' 0?2-1:4-1 '',L1,i4,i4)', lx, resi, j
      j = 0
      lx = a2bina(' 3*2 ',resi,2,0,' ',j,-1)
      print '(1x,'' 3*2 '',L1,i4,i4)', lx, resi
      j = 0
      lx = a2bina(' 3==2 ',resi,2,0,' ',j,-1)
      print '(1x,'' 3==2 '',L1,i4,i4)', lx, resi
      j = 0
      lx = a2bina(' 3>2 ',resi,2,0,' ',j,-1)
      print '(1x,'' 3>2 '',L1,i4,i4)', lx, resi
      j = 0
      lx = a2bina(' 3<2 ',resi,2,0,' ',j,-1)
      print '(1x,'' 3<2 '',L1,i4,i4)', lx, resi
      j = 0
      lx = a2bina(' 2==2 ',resi,2,0,' ',j,-1)
      print '(1x,'' 2==2 '',L1,i4,i4)', lx, resi
      j = 0
      lx = a2bina(' 2>2 ',resi,2,0,' ',j,-1)
      print '(1x,'' 2>2 '',L1,i4,i4)', lx, resi
      j = 0
      lx = a2bina(' 2<3 ',resi,2,0,' ',j,-1)
      print '(1x,'' 2<3 '',L1,i4,i4)', lx, resi
      j = 0
      lx = a2bina(' 2<3 ',resl,0,0,' ',j,-1)
      print '(1x,'' 2<3 '',L1,L1)', lx, resl
      j = 0
      lx = a2bina('  0  ',resl,0,0,' ',j,-1)
      print '(1x,''  0  '',L1,L1)', lx, resl
      j = 0
      lx = a2bina('  1  ',resl,0,0,' ',j,-1)
      print '(1x,''  1  '',L1,L1)', lx, resl

C --- Test a2bina-specific ---
      print '(/'' ------- test a2bina-specific -------'')'

      j = 0
      lx = a2bina(' yy>=3 ',resl,0,0,' ',j,-1)
      print '(1x,'' yy>=3 '',L1,L1)', lx, resl

      j = 0
      lx = a2bina(' yy>=3.001 ',resl,0,0,' ',j,-1)
      print '(1x,'' yy>=3.001 '',L1,L1)', lx, resl

      j = 0
      lx = a2bina('yz=3?1+2*3+4/5:2*(3+cos(pi)*-1), x=3, yz/=x ',
     .  resd,4,0,' ',j,-1)
      print '(1x,''yz=3?1+2*3+4/5:2*(3+cos(pi)*-1), x=3, yz+=x '',
     .  L1,1pe20.10,i4)',lx, resd(1), j

      j = 0
      lx = a2bina(' x-=yz-2.0. ',resd,4,0,' ',j,9)
      print '(1x,'' x-=yz-2.0. '',L1,2(1pe20.10),i4)',lx,resd(1),
     .  resd(1)-3,j

      j = 0
      lx = a2bina(' 1,2,x+1 ',resd,4,0,' ',j,-1)
      print '(1x,'' 1,2,x+1 '',L1,2(1pe20.10),i4)',lx,resd(1),
     .  resd(1)-3,j

C --- Test calling of symbolic vectors ---
C#ifdef POINTER
      print '(/'' ------- test vectors -------'')'
      j = 0
      lx = a2bina(' xyz(3) ',resd,4,0,' ',j,-1)
      print '(1x,'' xyz(3) '',L1,1pe20.10,i4)', lx, resd(1), j

      j = 0
      lx = a2bina(' xyz(4)*abc(3) ',resd,4,0,' ',j,-1)
      print '(1x,'' xyz(4)*abc(3) '',L1,1pe20.10,i4)', lx, resd(1), j
C#endif
      end
