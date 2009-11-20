      subroutine gcontr(z,nrz,nx,ny,cv,ncv,zmax,cs,bitmap,drwc,dparms)
C- draw a contour through equal values of an array.
C ----------------------------------------------------------------
Ci Inputs
Ci   z,nrz,nx,ny: array of the dependent variable; z is dimensioned
Ci         as z(nrz,*).  z is tabulated on a mesh of (nx,ny) elements.
Ci   cv,ncv: values of z for which contours are to be drawn and number
Ci         of values
Ci   dparms: parameters passed to drwc
Ci   zmax: maximum value of z for consideration.  A value of
Ci         z(i,j) > zmax is used as a signal that point and the
Ci         grid line segments radiating from that point to its
Ci         neighbors are to be excluded from the contour.
Ci   cs    contour color
Ci   bitmap: a work area large enough to hold 2*nx*ny*ncv bits.  It is
Ci         accessed by low-level routines, which are described below.
Ci         let j be the number of useful bits in each word of bitmap,
Ci         as determined by the user machine and implementation of
Ci         the bitmap manipulation subprograms described below.  then
Ci         the number of words required for the bitmap is the floor of
Ci         (2*nx*ny*ncv+j-1)/j.
Ci   drwc is a user-provided subroutine used to draw contours; see
Ci         description before the provided routine drwc.
Co Outputs
Co   Contours are plotted
Cv Verbosity
Cv   40: prints out number of elements in mesh, contours
Cr Remarks
Cr   drwc is the user-supplied line drawing subprogram described above.
Cr   drwc may be sensitive to the host computer and to the plot device.
Cr   fill0 is used to fill a bitmap with zeroes.  call fill0 (bitmap,n)
Cr   fills the first n bits of bitmap with zeroes.
Cr   mark1 is used to place a 1 in a specific bit of the bitmap.
Cr   call mark1 (bitmap,n) puts a 1 in the nth bit of the bitmap.
Cr   iget is used to determine the setting of a particular bit in the
Cr   bitmap.  i=iget(bitmap,n) sets i to zero if the nth bit of the
Cr   bitmap is zero, and sets i to one if the nth bit is one.
Cr   fill0, mark1 and iget are machine sensitive.
Cr   ---
Cr   l1 and l2 contain limits used during the spiral search for the
Cr     beginning of a contour.
Cr   i1, i2 and i3 are used for subscript computations during the
Cr     examination of lines from z(i,j) to its neighbors.
Cr   ij stores subcripts used during the spiral search.
Cr   xint is used to mark intersections of the contour under
Cr     consideration with the edges of the cell being examined.
Cr   xy is used to compute coordinates for the drwc subroutine.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer bitmap(1),nrz,nx,ny,ncv
      double precision z(nrz,1),zmax,dparms(*),cv(ncv),cs(0:3)
C Local parameters
      integer i,ibkey,icur,icv,idir,iedge,iflag,ii,imax,imin,ix,
     .        j,jcur,jj,jmax,jmin,k,ks,l,ni,nxidir,iget,iprint,jump
      integer l1(4), l2(4), ij(2), i1(2), i2(2), i3(6)
      double precision cval,dmax,x,y,z1,z2,zz
      double precision xint(4),xy(2)
      external drwc,fill0,iget,iprint,mark1

      equivalence (l2(1),imax), (l2(2),jmax), (l2(3),imin),
     *  (l2(4),jmin)
      equivalence (ij(1),i), (ij(2),j)
      equivalence (xy(1),x), (xy(2),y)
c
      data l1(3) /-1/, l1(4) /-1/
      data i1 /1,0/, i2 /1,-1/, i3 /1,0,0,1,1,0/
c
      l1(1) = nx
      l1(2) = ny
      dmax = zmax
c
c     set the current pen position.  the default position corresponds
c     to z(1,1).
c
      x = 1.0d0
      y = 1.0d0
      call drwc(x, y, 6, dparms, cs)
      icur = max0(1,min0(int(x),nx))
      jcur = max0(1,min0(int(y),ny))
c
c     clear the bitmap
c
      call fill0(bitmap, 2*nx*ny*ncv)
c
c     search along a rectangular spiral path for a line segment having
c     the following properties:
c          1.  the end points are not excluded,
c          2.  no mark has been recorded for the segment,
c          3.  the values of z at the ends of the segment are such that
c              one z is less than the current contour value, and the
c              other is greater than or equal to the current contour
c              value.
c
c     search all boundaries first, then search interior line segments.
c     note that the interior line segments near excluded points may be
c     boundaries.
c
      ibkey = 0
   10 i = icur
      j = jcur
   20 imax = i
      imin = -i
      jmax = j
      jmin = -j
      idir = 0
c     direction zero is +i, 1 is +j, 2 is -i, 3 is -j.
   30 nxidir = idir + 1
      k = nxidir
      if (nxidir.gt.3) nxidir = 0
   40 i = iabs(i)
      j = iabs(j)
      if (z(i,j).gt.dmax) go to 140
      l = 1
c     l=1 means horizontal line, l=2 means vertical line.
   50 if (ij(l).ge.l1(l)) go to 130
      ii = i + i1(l)
      jj = j + i1(3-l)
      if (z(ii,jj).gt.dmax) go to 130
c      assign 100 to jump !takao
      jump=100
c     the next 15 statements (or so) detect boundaries.
   60 ix = 1
      if (ij(3-l).eq.1) go to 80
      ii = i - i1(3-l)
      jj = j - i1(l)
      if (z(ii,jj).gt.dmax) go to 70
      ii = i + i2(l)
      jj = j + i2(3-l)
      if (z(ii,jj).lt.dmax) ix = 0
   70 if (ij(3-l).ge.l1(3-l)) go to 90
   80 ii = i + i1(3-l)
      jj = j + i1(l)
      if (z(ii,jj).gt.dmax) go to 90
c      if (z(i+1,j+1).lt.dmax) go to jump, (100, 280)
      if (z(i+1,j+1).lt.dmax) then
         if(jump==100) goto 100
         if(jump==280) goto 280
      endif
   90 ix = ix + 2
c      go to jump, (100, 280)
       if(jump==100) goto 100
       if(jump==280) goto 280
  100 if (ix.eq.3) go to 130
      if (ix+ibkey.eq.0) go to 130
c     now determine whether the line segment is crossed by the contour.
      ii = i + i1(l)
      jj = j + i1(3-l)
      z1 = z(i,j)
      z2 = z(ii,jj)
      do 120 icv=1,ncv
        if (iget(bitmap,2*(nx*(ny*(icv-1)+j-1)+i-1)+l).ne.0) go to 120
        if (cv(icv).le.dmin1(z1,z2)) go to 110
        if (cv(icv).le.dmax1(z1,z2)) go to 190
  110   call mark1(bitmap, 2*(nx*(ny*(icv-1)+j-1)+i-1)+l)
  120 continue
  130 l = l + 1
      if (l.le.2) go to 50
  140 l = mod(idir,2) + 1
      ij(l) = isign(ij(l),l1(k))
c
c     lines from z(i,j) to z(i+1,j) and z(i,j+1) are not satisfactory.
c     continue the spiral.
c
  150 if (ij(l).ge.l1(k)) go to 170
      ij(l) = ij(l) + 1
      if (ij(l).gt.l2(k)) go to 160
      go to 40
  160 l2(k) = ij(l)
      idir = nxidir
      go to 30
  170 if (idir.eq.nxidir) go to 180
      nxidir = nxidir + 1
      ij(l) = l1(k)
      k = nxidir
      l = 3 - l
      ij(l) = l2(k)
      if (nxidir.gt.3) nxidir = 0
      go to 150
  180 if (ibkey.ne.0) return
      ibkey = 1
      go to 10
c
c     an acceptable line segment has been found.
c     follow the contour until it either hits a boundary or closes.
c
  190 iedge = l
      cval = cv(icv)
      if (ix.ne.1) iedge = iedge + 2
      iflag = 2 + ibkey
      xint(iedge) = (cval-z1)/(z2-z1)
  200 xy(l) = dble(ij(l)) + xint(iedge)
      xy(3-l) = dble(ij(3-l))
      call mark1(bitmap, 2*(nx*(ny*(icv-1)+j-1)+i-1)+l)
      call drwc(x, y, iflag+10*icv, dparms, cs)
      if (iflag.lt.4) go to 210
      icur = i
      jcur = j
      go to 20
c
c     continue a contour.  the edges are numbered clockwise with
c     the bottom edge being edge number one.
c
  210 ni = 1
      if (iedge.lt.3) go to 220
      i = i - i3(iedge)
      j = j - i3(iedge+2)
  220 do 250 k=1,4
        if (k.eq.iedge) go to 250
        ii = i + i3(k)
        jj = j + i3(k+1)
        z1 = z(ii,jj)
        ii = i + i3(k+1)
        jj = j + i3(k+2)
        z2 = z(ii,jj)
        if (cval.le.dmin1(z1,z2)) go to 250
        if (cval.gt.dmax1(z1,z2)) go to 250
        if (k.eq.1) go to 230
        if (k.ne.4) go to 240
  230   zz = z1
        z1 = z2
        z2 = zz
  240   xint(k) = (cval-z1)/(z2-z1)
        ni = ni + 1
        ks = k
  250 continue
      if (ni.eq.2) go to 260
c
c     the contour crosses all four edges of the cell being examined.
c     choose the lines top-to-left and bottom-to-right if the
c     interpolation point on the top edge is less than the
c     interpolation point on the bottom edge.  otherwise, choose the
c     other pair. This method produces the same results if the axes
c     are reversed.  thecontour may close at any edge, but must not
c     cross itself inside any cell.
c
      ks = 5 - iedge
      if (xint(3).lt.xint(1)) go to 260
      ks = 3 - iedge
      if (ks.le.0) ks = ks + 4
c
c     determine whether the contour will close or run into a boundary
c     at edge ks of the current cell.
c
  260 l = ks
      iflag = 1
c      assign 280 to jump
      jump=280
      if (ks.lt.3) go to 270
      i = i + i3(ks)
      j = j + i3(ks+2)
      l = ks - 2
  270 if (iget(bitmap,2*(nx*(ny*(icv-1)+j-1)+i-1)+l).eq.0) go to 60
      iflag = 5
      go to 290
  280 if (ix.ne.0) iflag = 4
  290 iedge = ks + 2
      if (iedge.gt.4) iedge = iedge - 4
      xint(iedge) = xint(ks)
      go to 200
c
      end
      subroutine fill0(bitmap, n)
C- fill the first n bits of bitmap with zeroes.
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   nbpw is the minimum number of significant bits per word used
Cr   by integer arithmetic.  this is usually one less than the
Cr   actual number of bits per word, but an important exception is
Cr   the cdc-6000 series of machines, where nbpw should be 48.
C ----------------------------------------------------------------
C Passed parameters
      integer n,bitmap(1)
C Local parameters
      integer i,loop,nblw,nbpw
c
      data nbpw /31/
c
      loop = n/nbpw
      nblw = mod(n,nbpw)
      if (loop.eq.0) go to 20
      do 10 i=1,loop
        bitmap(i) = 0
   10 continue
   20 if (nblw.ne.0)
     .  bitmap(loop+1) = mod(bitmap(loop+1),2**(nbpw-nblw))
      return
      end
      subroutine mark1(bitmap, n)
C- put a one in the nth bit of bitmap.
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   nbpw is the minimum number of significant bits per word used
Cr   by integer arithmetic.  this is usually one less than the
Cr   actual number of bits per word, but an important exception is
Cr   the cdc-6000 series of machines, where nbpw should be 48.
C ----------------------------------------------------------------
C Passed parameters
      implicit none
      integer bitmap(1), n
C Local parameters
      integer nword,nbpw,i,nbit
c
      data nbpw /31/
c
      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      i = 2**(nbpw-nbit-1)
      bitmap(nword+1) = bitmap(nword+1)+ i*(1-mod(bitmap(nword+1)/i,2))
      return
      end
      function iget(bitmap, n)
C- iget=0 if the nth bit of bitmap is zero, else iget is one.
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   nbpw is the minimum number of significant bits per word used
Cr   by integer arithmetic.  this is usually one less than the
Cr   actual number of bits per word, but an important exception is
Cr   the cdc-6000 series of machines, where nbpw should be 48.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer bitmap(1), n
C Local parameters
      integer nword,iget,nbit,nbpw
      data nbpw /31/
c
      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      iget = mod(bitmap(nword+1)/2**(nbpw-nbit-1),2)
      return
      end
      subroutine drwc(x, y, iflag, lt, cs)
C- do output for gcontr.
C ----------------------------------------------------------------
Ci Inputs
Ci   x,y,iflag (see remarks)
Co Outputs
Co
Cv Verbosity
Cv   >110: prints out all pen movements
Cr Remarks
Cr   let nx = integer part of x, fx = fractional part of x.
Cr   then x should be interpreted such that increases in nx
Cr   correspond to increases in the first subscript of z, and
Cr   fx is the fractional distance from the abscissa corresponding
Cr   to nx to the abscissa corresponding to nx+1,
Cr   and y should be interpreted similarly for the second
Cr   subscript of z.
Cr
Cr   the low-order digit of iflag will have one of the values:
Cr      1 - continue a contour,
Cr      2 - start a contour at a boundary,
Cr      3 - start a contour not at a boundary,
Cr      4 - finish a contour at a boundary,
Cr      5 - finish a closed contour (not at a boundary).
Cr          note that requests 1, 4 and 5 are for pen-down
Cr          moves, and that requests 2 and 3 are for pen-up
Cr          moves.
Cr      6 - set x and y to the approximate 'pen' position, using
Cr          the notation discussed above.  this call may be
Cr          ignored, the result being that the 'pen' position
Cr          is taken to correspond to z(1,1).
Cr   iflag/10 is the contour number.
Cr
Cr   entry drwcin scale the (integral)
Cr   x,y pairs gcontr uses to the actual dimensions of the plot.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x,y,lt(4,1),cs(0:3)
      integer iflag
C Local parameters
      logical dcmpre
      integer nx,ny,ic,jump,nx0,ny0,iprint,iedgst,iedge
      double precision xscale,yscale,xmin0,ymin0,xmax0,ymax0,x1,x2
      double precision xcur,ycur,xmin,xmax,ymin,ymax,xst,yst
      double precision x0,y0,Area
      save xcur, ycur, xmin, xmax, ymin, ymax, nx, ny, xst, yst
      save x0, y0, Area, iedgst
      dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8

      xscale(x) = (xmax-xmin)*(x-1)/(nx-1) + xmin
      yscale(x) = (ymax-ymin)*(y-1)/(ny-1) + ymin

      ic = iflag/10
      jump = mod(iflag,10)

      if (iprint() .ge. 80 .and. jump .ge. 2 .and. jump .le. 3)
     .  print 333, ic, nint(lt(1,ic)),lt(2,ic),lt(3,ic)
  333 format('drwc: drawing contour',i3,'  lt=',i2,2f7.3)

      go to (10, 20, 30, 40, 50, 60), jump
C --- Continue a contour ---
   10 continue
C     if (iprint() .gt. 90) write (*,999) ic, xscale(x), yscale(y)
      call info8(90,0,0,' continue contour%,3i to  %,4;4g   %,4;4g',
     .  ic,xscale(x),yscale(y),0,0,0,0,0)

      call drw(xscale(x),yscale(y))
C     Area = r cross dr/2 = (xdy - ydx)/2
      Area = Area + (x+xcur)/4*(y-ycur) - (y+ycur)/4*(x-xcur)
      go to 70
C --- Start a new contour on the boundary ---
   20 continue
C     if (iprint() .ge. 80) write (*,998) ic, xscale(x), yscale(y)
      call plntyp(nint(lt(1,ic)),nint(lt(2,ic)),lt(3,ic),lt(4,ic),
     .  0d0,0d0)
C     xst,yst are scaled versions of x0,y0
C     Only keep both for historical reasons
      x0 = x
      y0 = y
      xst = xscale(x)
      yst = yscale(y)
      call mve(xscale(x),yscale(y))
      Area = 0
C     Determine what the wall makes the initial boundary
C     1s  digit 0,1,2 for no wall, bottom wall, top wall
C     10s digit 0,1,2 for no wall, left wall, right wall
      iedgst = 0
      if (x .eq. 1) iedgst = iedgst + 10
      if (x .eq. nx) iedgst = iedgst + 20
      if (y .eq. 1) iedgst = iedgst + 1
      if (y .eq. ny) iedgst = iedgst + 2
      call info8(80,0,0,'    start contour%,3i on the '//
     .  '%?#(n==01)#bottom boundary,##%-1j'//
     .  '%?#(n==02)#top boundary,   ##%-1j'//
     .  '%?#(n==10)#left boundary,  ##%-1j'//
     .  '%?#(n==20)#right boundary, ##'//
     .  '  x=%,4;4g   y=%,4;4g',
     .  ic,iedgst,xscale(x),yscale(y),Area,0,0,0)
      go to 70
C --- Start a new contour in the interior ---
   30 continue
      call plntyp(nint(lt(1,ic)),nint(lt(2,ic)),lt(3,ic),lt(4,ic),
     .  0d0,0d0)
C     if (iprint() .ge. 80) write (*,997) ic, xscale(x), yscale(y)
      call info8(80,0,0,'    start contour%,3i in the interior at'//
     .  '  x=%,4;4g   y=%,4;4g',
     .  ic,xscale(x),yscale(y),Area,0,0,0,0)
C     xst,yst are scaled versions of x0,y0
      x0 = x
      y0 = y
      xst = xscale(x)
      yst = yscale(y)
      call mve(xscale(x),yscale(y))
      Area = 0
      iedgst = 0
      go to 70
C --- Close a contour at a boundary ---
   40 continue
C     if (iprint() .ge. 80) write (*,996) ic, xscale(x), yscale(y)

C     Determine what the wall makes the final boundary
C     1s  digit 0,1,2 for no wall, bottom wall, top wall
C     10s digit 0,1,2 for no wall, left wall, right wall
      iedge = 0
      if (x .eq. 1) iedge = iedge + 10
      if (x .eq. nx) iedge = iedge + 20
      if (y .eq. 1) iedge = iedge + 1
      if (y .eq. ny) iedge = iedge + 2

C ... Contour encloses a corner
      if (x .ne. x0 .and. y .ne. y0) then
        call drw(xscale(x),yscale(y))
        Area = Area + (x+xcur)/4*(y-ycur) - (y+ycur)/4*(x-xcur)
C   ... STARTING x at xmin or xmax. ENDING y at ymin or ymax
        if ((dcmpre(xst,xmin) .or. dcmpre(xst,xmax)) .and.
     .      (dcmpre(yscale(y),ymin) .or. dcmpre(yscale(y),ymax))) then
C         Move horizontally to corner
          call drwe(xst,yscale(y),cs)
C         Add to line integral, segment from (x,y) to (x0,y) = -y dx/2
          Area = Area - y/2 * (x0-x)
C         Add to line integral, segment from (x0,y) to (x0,y0) = x0*dy/2
          Area = Area + x0/2 * (y0-y)
C   ... STARTING y at ymin or ymax. ENDING x at xmin or xmax.
        elseif ((dcmpre(yst,ymin) .or. dcmpre(yst,ymax)) .and.
     .      (dcmpre(xscale(x),xmin) .or. dcmpre(xscale(x),xmax))) then
C         Move horizontally to corner
          call drwe(xscale(x),yst,cs)
C         Add to line integral, segment from (x,y) to (x,y0) = x dy/2
          Area = Area + x/2 * (y0-y)
C         Add to line integral, segment from (x,y0) to (x0,y0) = -y0*dx/2
          Area = Area - y0/2 * (x0-x)
        endif
C ... Case starting wall same as ending wall 
      else
        call drwe(xscale(x),yscale(y),cs)
        Area = Area + (x+xcur)/4*(y-ycur) - (y+ycur)/4*(x-xcur)
C       Start, stop on same wall
        if (iedge .eq. iedgst) then
C         Top or bottom wall
          if (iedge .eq. 1 .or. iedge .eq. 2) then
            Area = Area + y/2 * (x-x0)
C         Left or right wall.  Add -x (dy)/2 to r cross dr
          else
            Area = Area - x/2 * (y-y0)
          endif
        else
        endif
      endif
      Area = dabs( Area * (xmax-xmin)/(nx-1) * (ymax-ymin)/(ny-1) )
      call info8(80,0,0,'   finish contour%,3i on the '//
     .  '%?#(n==01)#bottom boundary,##%-1j'//
     .  '%?#(n==02)#top boundary,   ##%-1j'//
     .  '%?#(n==10)#left boundary,  ##%-1j'//
     .  '%?#(n==20)#right boundary, ##'//
     .  '  x=%,4;4g   y=%,4;4g   Area=%,4;4g',
     .  ic,iedge,xscale(x),yscale(y),Area,0,0,0)
      go to 70
C --- Close a new contour in the interior ---
   50 continue
      call drwe(xscale(x),yscale(y),cs)
      Area = Area + (x+xcur)/4*(y-ycur) - (y+ycur)/4*(x-xcur)
      Area = dabs( Area * (xmax-xmin)/(nx-1) * (ymax-ymin)/(ny-1) )
C     if (iprint() .ge. 80) write (*,995) ic, xscale(x), yscale(y), Area
      call info8(80,0,0,'   finish contour%,3i in the interior at'//
     .  '  x=%,4;4g   y=%,4;4g   Area=%,4;4g',
     .  ic,xscale(x),yscale(y),Area,0,0,0,0)

      go to 70
C --- Pen position ---
   60 continue
      if (iprint() .gt. 90) write (*,994)
      x = xcur
      y = ycur
   70 continue
      xcur = x
      ycur = y
      return
C  999 format (' continue contour', i3,' to ', 0p,2g14.6)
C  998 format ('    start contour', i3,' on the boundary at ', 0p,2g14.6)
C  997 format ('    start contour', i3,' in the interior at ', 0p,2g14.6)
C  996 format ('   finish contour', i3,' on the boundary at ', 0p,2g14.6)
C  995 format ('   finish contour', i3,' in the interior at ', 0p,2g14.6,
C     .  '  Area =', g14.6)
  994 format (' request for current pen position')

      entry drwcin(nx0, xmin0, xmax0, ny0, ymin0, ymax0)
      nx = nx0
      ny = ny0
      xmin = xmin0
      xmax = xmax0
      ymin = ymin0
      ymax = ymax0

      end
