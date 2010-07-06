      subroutine fmain
C     implicit none
      integer wksize,ox,oy,oz,ow,lenx,leny,lenz
      parameter(wksize= 440)
      integer w(wksize)
      common /w/ w
      logical cmdopt
      character*40 strn
      call wkinit(wksize)
      call wkprnt(1)

      lenx = 10
      leny = 20
      lenz = lenx-4
      call defrr(oz, lenz)
      call defrr(ox, lenx)
      call defrr(oy, leny)
      call dvset(w(oz),1,lenz,6d0)
      call dvset(w(ox),1,lenx,1d0)
      call dvset(w(oy),1,leny,2d0)
      call prmx('x',w(ox),lenx,lenx,1)
      call prmx('y',w(oy),leny,leny,1)
      call prmx('z',w(oz),lenz,lenz,1)
      print 333, 'stack 0..2: oy,ox,oz=',oy,ox,oz
  333 format(1x,a,4i4)
      call defps3(oy,ox,oz)
      call wkchk('testing links after swap')
      print 333, 'stack 0..2: oz,oy,ox=',oz,oy,ox
      call prmx('x',w(ox),lenx,lenx,1)
      call prmx('y',w(oy),leny,leny,1)
      call prmx('z',w(oz),lenz,lenz,1)

      print *, 'this makes an error unless ow allocated'
      call defrr(ow, leny)
      call dvset(w(ow),1,leny,4d0)
      print 333, 'stack 0..3: ow,oz,oy,ox=',ow,oz,oy,ox
      call defps4(ow,oz,oy,ox)
      print 333, 'stack 0..3: ox,ow,oz,oy=',ox,ow,oz,oy
      call prmx('x',w(ox),lenx,lenx,1)
      call prmx('y',w(oy),leny,leny,1)
      call prmx('z',w(oz),lenz,lenz,1)
      call prmx('w',w(ow),leny,leny,1)


      print *, 'this destroys links and generates error unless fast set'
      call defrr(oz, lenx)
      if (cmdopt('-fast',5,0,strn)) call wkfast(.true.)
      call dpzero(w(ox),3*lenx)
      call defrr(ow, lenx)
      end
      subroutine prmx(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
C     implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2)
      character*(10) fmt, strn*(*), outs*80
      integer i,j,fopna,i1mach
      fmt = '(9f8.4)'
      ifi = 6
*      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
      print *, strn
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(i,j,1), j=1,nc)
      write(ifi,*)
*      call fclose(ifi)
*      outs = 'prm: done writing data '//strn
*      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
*      read(*,'(a80)') outs
*      if (outs .eq. 'q') call rx0('quit in prmx')
      end
