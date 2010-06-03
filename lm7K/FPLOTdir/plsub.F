C#define INCH_UNITS
C --- Fortran primitives for plot procedures. ---
Cr This is a collection of routines for drawing lines and symbols.
Cr The hardware dependence is mostly concentrated into an initialization
Cr procedure (pltini) and one (pltam) that draws straight lines, to the
Cr extent it was feasible to do so.  The core routines are:
Cr   pltini(character device, integer opt)
Cr     defines and possibly initializes the current plot device
Cr   pltam(real x, real y, integer opt)
Cr     Moves the pen from the current position to point x-y on the
Cr     current device, in the units appropriate to the plot medium.
Cr     opt=1 plots with pen down; opt=0 plots with pen up.
Cr     All pen moves are done by calling pltam.
Cr     (These units are set in pltini).
Cr   pltdmp(io)
Cr     Dumps the current plot to the output device.  This is hardware-
Cr     dependent.  When it is applicable and io=1, pltdmp prompts for
Cr     a carriage return before executing.
Cr   plsym(xp,yp,zp,np,type,syma,symcs)
Cr     draws a symbol at a collection of points (xp,yp)
Cr
Cr Besides the (hardware-dependent) medium units, there are two sets of
Cr device-independent units.  By default, graphics units range from (0,0)
Cr to (1,1) for the largest square that fits within the plot medium.
Cr User's units define the user's values of x and y at the corners
Cr of the graph, and obviously will depend on the specific plot.
Cr There are several routines for defining the various units
Cr   pltstu(l,r,b,t):
Cr     defines the corners of a graph, in user's units
Cr   pltstp(l,r,b,t):
Cr     restricts (clips) a graph to a portion of the medium's surface
Cr     (in graphics units).  no check to see within (l,r,b,t)
Cr   pltstg(l,r,b,t):
Cr     Change the corners of the medium (in graphics units) from their
Cr     default values (0,0) and (1,1).  Not normally needed.
Cr   plthet(alpha):
Cr     Shears the orthogonal (x,y) axes by the linear transformation:
Cr       x->x
Cr       y->y*sin(alpha) + x*cos(alpha)
Cr      Default alpha is pi/2 (orthogonal x,y)
Cr The following procedures are for drawing straight lines.
C  Arguments (x,y) in user's units, except when stated:
Cr   mve(x,y)  moves the pen to (x,y) from current posn, pen up
Cr   drw(x,y)  draws the pen to (x,y) from current posn, pen down
Cr   mver(x,y) like mve, but motion relative to current x,y
Cr   drwr(x,y) like drw, but motion relative to current x,y
Cr   pposu(x,y) returns the current pen position in (x,y).
Cr   drwe(x,y,cs) draws the last point in a curve.  Area enclosed by
Cr             curve may be be filled in, according to cs.
Cr   drwre(x,y,cs) like drwe, but motion relative to current x,y
Cr   mveg(x,y) like mve, but x,y in GU
Cr   mvegr(x,y) like mveg, but motion relative to current x,y
Cr   pposg(x,y) like pposu, but in GU
Cr
Cr Other routines:
Cr   plcrv  draws a smooth curve for a collection of points.
Cr          For postscript, it writes the file directly instead of
Cr          calling the generic pltam, to preserve continuity of
Cr          dotted, dashed lines, and some other features.
Cr   plcrve closes a curve, possibly shading area enclosed by curve
Cr   frme,ticset,ticnxt draw a frame and tic marks.
Cr   pltsts sets current plot to linear, log scale
Cr   pstr   puts a string at current x,y
Cr   setfnt sets the font for drawing strings.
Cr   plntyp inputs variables to make dashed, broken or dotted lines.
Cr   pltcol inputs variables for color (not implemented)
Cr   plndpi sets the medium resolution (dots per inch)
Cr   key    draws a key
Cr   plstmv sets or reads the plot medium boundaries as symbolic
Cr          variables.  This permits user to redefine, e.g. the edges
Cr          of the page.
Cr
Cr Some Postscript-specific information and routines are:
Cr   psinit writes to the file stream some postscript macros
Cr   pslabl writes a string
Cr For postscript, the "edges" of the default to something smaller than
Cr a 8+1/2 x 11 inch page.  plstmv permits you to override the internal
Cr defaults (now set at (l,r,t,b) = (144,612-72,144,792-72) points).
Cr
Cr The initialization routine does the following:
Cr   1. sets the corners of the medium, in MU
Cr   2. defines scalxc,scalyc converting from MU into centimeters
Cr   3. defines xb,yb, the "line thickness" (for devices where
Cr      this is not adjustable, e.g. Hercules graphics card)
Cr   4. other device-specific initialization (eg calls psinit)
Cr   5. Defines the default graphics units and sets the clip region
Cr      to the largest available square.
Cr
Cv Common block variables:
Cv The variables corresponding to these units are:
Cv   plgl,r,b,t: graphics units demarcating value of x and y at four
Cv               corners of the plotting medium.  They default to (0,0)
Cv               and (1,1) but can be set explicitly.
Cv   plpl,r,b,t: values of x and y, in graphics units, at four corners
Cv               of a plot (useful when a plot is clipped to a portion
Cv               of the entire medium).
Cv   plul,r,b,t: user's units demarcating value of x and y at four
Cv               corners of a plot.
Cv   plml,r,b,t: units of the actual plot medium
Cv               (set automatically by pltini)
Cv   scalxg,scalyg,offxg,offyg (calculated in pltstg), make the
Cv               linear transformation between medium units
Cv               and graphics units:
Cv               xm = (offxg+xg)*scalxg, ym = (offyg+yg)*scalyg
Cv   scalxu,scalyu,offxu,offyu (calculated in pltstu), make the
Cv               linear transformation between medium units
Cv               and user units:
Cv               xm = (offxu+xu)*scalxu, ym = (offyu+yu)*scalyu
Cv   xcm,ycm:    current pen position in medium's units
Cv   xcu,ycu:    current pen position in user's units
Cv   xcc,ycc:    not used
Cv   strtch      the ratio of MU unit lengths in x,y that produces a
Cr               true square as seen by the viewer.  Unit lengths in x,y
Cr               may have differing true dimensions, as they probably
Cr               will when units refer to pixels.
Cv   stream:     file logical unit for plot output
Cr   xb,yb       the "line thickness" in MU (for devices where
Cr               this is not adjustable)
Cr   lintyp,ltbld define current
Cr               line type and boldness.  lintyp:
Cr               0 for none, 1 for solid, 2 for dashed, 3 for dotted
Cr   linl1,linl2,linl3,linl4, sequence of lengths for dashed lines
Cr               (dash, space, dash, space)
Cv   lnprd       scaling of line segments to medium units.  Defined
Cv               so that "typically sized" line segments have size of
Cv               order 1, independently of the medium.
Cv   noclip      true if no clipping has been set
Cv   medium      index defining current medium
Cv   logx,logy   true when user's units refer to a logarithmic scale
Cv   pstat       current pen status: 1 is down, 0 is up
Cv   ndpi,ndwrt  resolution, in dots per inch. number of digits
Cv               needed to reach medium's resolution
C ----------------------------------------------------------------
      subroutine pltini(device,opt,lrotp)
C- Initialization for plot
C ----------------------------------------------------------------
Ci Inputs
Ci   device: device name
Ci   opt:    see Remarks
Ci   lrot:   rotate figure
Co Outputs
Co   Initializes plot medium
Cr Remarks
Cr   pltini can be invoked using a number (opt) or (if 0) a device name.
Cv Variables
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*(*) device
      integer opt
      logical lrotp

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      integer nnames,fopna,idum,i1mach,isw
      double precision transx,transy,scale,tmp

      parameter (nnames=3)
      character*10 names(nnames)
      data names /'ascii','hercules','ps'/

      lrot = lrotp

C --- If no number specified, seek one corresponding to name ---
      medium = opt
      if (medium .eq. 0) then
   10   medium = medium+1
        if (names(medium) .eq. device .or. medium .ge. nnames) goto 20
        goto 10
      endif

   20 goto (110,120,130), medium
      call rx('PLTINI: device not present')

C --- Ascii: write to file unit 'stream' ---
  110 continue
c      stream = i1mach(2)
      if (stream .ne. i1mach(2)) idum = fopna('wplt',stream,0)
C#ifdef INCH_UNITS
C --- Set medium's boundaries  ---
      lnprd = .1d0
      plml = 1.5d0
      plmr = 7.5d0
      plmb = 0d0
      plmt = 9d0
C --- Define size of cm units in terms of medium units ---
      scalxc = 1/2.54d0
      scalyc = 1/2.54d0
      strtch = dabs(scalxc/scalyc)
      offxc = 0d0
      offyc = 0d0
      xb = 1d0/150d0
      yb = xb
      goto 30
C#endif INCH_UNITS
C#ifdefC CM_UNITS
CC --- Set medium's boundaries  ---
C      lnprd = .27
C      plml = 0.
C      plmr = 24.
C      plmb = 0.
C      plmt = 18.
CC --- Define size of cm units in terms of medium units ---
C      scalxc = 1.
C      scalyc = 1.
C      strtch = abs(scalxc/scalyc)
C      offxc = 0.
C      offyc = 0.
C      goto 30
C#endif INCH_UNITS

C --- Hercules graphics card ---
  120 continue
      lnprd = 8
C --- Set medium's boundaries  ---
      plml = 0d0
      plmr = 719d0
      plmb = 317d0
      plmt = 0d0
C --- Define size of cm units in terms of medium units ---
C Calculate scaling by matching physical dimensions of the screen
C and putting (0,0) in cm units at (plml,plmb)
      scalxc = (plmr-plml)/19.5d0
      scalyc = (plmt-plmb)/14.0d0
      strtch = dabs(scalxc/scalyc)
      offxc = plml/scalxc - 0d0
      offyc = plmb/scalyc - 0d0
      goto 30

C --- PostScript ---
  130 continue
      idum = fopna('ps',stream,0)
C --- Set medium's boundaries  ---
      lnprd = 7.2d0
      plml = 144
      plmr = 612-72
      plmb = 144
      plmt = 792-72
      call plstmv(plml,plmr,plmb,plmt,' ',.true.)
      bbox(1) = 9999
      bbox(2) = -9999
      bbox(3) = 9999
      bbox(4) = -9999

C --- Define size of cm units in terms of medium units ---
      scalxc = 1d0/(2.54d0*72d0)
      scalyc = 1d0/(2.54d0*72d0)
      strtch = dabs(scalxc/scalyc)
      offxc = 0d0
      offyc = 0d0
      xb = 1d0/150d0
      yb = xb
      transx = 0
      transy = 0
      scale = 1
      call psinit(stream,isw(lrot),transx,transy,scale)
      goto 30

C --- Setup common to all devices ---
   30 continue
C      if (lrot) then
C        tmp  = plmt
C        plmt = plmr
C        plmr = plmb
C        plmb = plml
C        plml = tmp
C        tmp    = scalxc
C        scalxc = scalyc
C        scalyc = tmp
C        strtch = dabs(scalxc/scalyc)
C        tmp   = offxc
C        offxc = offyc
C        offyc = tmp
C        tmp = xb
C        xb = yb
C        yb = tmp
C      endif

C --- Define GU ---
      if (dabs((plmt-plmb)*strtch) .gt. dabs(plmr-plml)) then
        call pltstg(0d0,1d0,0d0,dabs((plmt-plmb)*strtch/(plmr-plml)))
      else
        call pltstg(0d0,dabs((plmr-plml)/(plmt-plmb)/strtch),0d0,1d0)
      endif

C --- Clip to largest available centered square ---
      plpl = (plgr+plgl)/2 - .5d0
      plpr = (plgr+plgl)/2 + .5d0
      plpb = (plgt+plgb)/2 - .5d0
      plpt = (plgt+plgb)/2 + .5d0
      end
      subroutine plstmv(plml,plmr,plmb,plmt,s,lread)
C- Set or read plot medium boundaries from as symbolic variables
      implicit none
      logical lread
      character *(*) s
      double precision plml,plmr,plmb,plmt
      character*4 sp,st,pnam
      logical a2bin,lxst
      integer i,n,ich,ls,it,ich0,iprint,i1mach
      double precision xx(4)
      data sp /'lbrt'/ st /', '/

      if (lread) then
        lxst = .false.
        do  10  i = 1, 4
          call getsyv('plm'//sp(i:i),xx(i),n)
C         call shosyv(0,0,0,i1mach(2))
          if (n .ne. 0) then
            lxst = .true.
            if (i .eq. 1) plml = xx(i)
            if (i .eq. 2) plmb = xx(i)
            if (i .eq. 3) plmr = xx(i)
            if (i .eq. 4) plmt = xx(i)
          endif
   10   continue
        if (lxst .and. iprint() .ge. 0) call awrit4(
     .  ' plstmv:  plml=%d  plmb=%d  plmr=%d  plmt=%d',
     .  ' ',80,i1mach(2),plml,plmb,plmr,plmt)
      else
        ls = len(s)
C        print *, s
        ich = 0
        call chrpos(s,'=',ls,ich)
        if (ich .ge. ls) goto 22
        ich = ich+1
        ich0 = ich
        do  20  i = 1, 4
          ich = ich0
          call chrps2(s,st,2,ls,ich,it)
C          print *, ich0, ich, s(ich0+1:ich+1)
          if (ich .gt. ich0) then
            if (a2bin(s,xx(i),4,0,st(it:it),ich0,-1)) then
              pnam = 'plm'//sp(i:i)
              call lodsyv(pnam,1,xx(i),n)
C             call lodsyv('plm'//sp(i:i),1,xx(i),n)
              if (st(it:it) .eq. ' ') goto 22
            endif
          else
            if (s(ich0:ich0) .eq. ' ') goto 22
            ich0 = ich0+1
          endif
   20   continue
   22   continue
      endif
      end
      subroutine pltstg(l,r,b,t)
C- set coordinates of graphics units
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision l,r,b,t

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      integer iprint,i1mach

      plgl = l
      plgr = r
      plgb = b
      plgt = t
      lnused = 0

C --- Set the scale and offsets ---
C Calculate scaling by taking the smaller of the x and y scales
C and offset by centering both graphics units in medium units.
      scalxg = dabs((plmr-plml)/(plgr-plgl)*strtch)
      scalyg = dabs((plmt-plmb)/(plgt-plgb))
      if (scalxg .le. scalyg) then
        scalxg = dsign(scalxg,plmr-plml)
        scalyg = dsign(scalxg/strtch,plmt-plmb)
      else
        scalxg = dsign(scalyg*strtch,plmr-plml)
        scalyg = dsign(scalyg,plmt-plmb)
      endif

      offxg = ((plml+plmr)/scalxg - (plgl+plgr))/2
      offyg = ((plmb+plmt)/scalyg - (plgb+plgt))/2

      if (iprint() .ge. 50)
     .call awrit4(' pltstg: medium corners (GU): lb = (%d %d)'//
     .'  rt = (%d %d)',' ',80,i1mach(2),plgl,plgb,plgr,plgt)
      if (iprint() .gt. 100)
     .call awrit4('         scal(x,y)g = %d %d  off(x,y)g = '//
     .'(%d %d)',' ',80,i1mach(2),scalxg,scalyg,offxg,offyg)
      if (iprint() .gt. 100)
     .call awrit4(' pltstg: medium corners (MU): lb = (%d %d)'//
     .'  rt = (%d %d)',' ',80,i1mach(2),plml,plmb,plmr,plmt)

      end
      subroutine pltstu(l,r,b,t,uuytox)
C- set coordinates of user's units
C ----------------------------------------------------------------
Ci Inputs
Ci   l,r,b,t: frame corners (UU)
Ci   uuytox: if zero, does nothing.  Otherwise relative lengths in UU
Ci     are linked to those in MU as:  dy/dx (UU) = uuytox * dy/dx (MU).
Ci     This is accomplished by resizing the current frame to the largest
Ci     size permitted by the above constraint.
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision l,r,b,t,uuytox

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      integer iprint,i1mach,is,scrwid
      double precision xx1,xx2
      logical error
      parameter (scrwid=100)
      character*(scrwid) s

      call lodsyv('ul',1,l,is)
      call lodsyv('ur',1,r,is)
      call lodsyv('ub',1,b,is)
      call lodsyv('ut',1,t,is)
      plul = l
      plur = r
      plub = b
      plut = t
      error = .true.
      if (logx) then
        if (l .le. 0 .or. r .le. 0) goto 99
        plul = dlog(l)
        plur = dlog(r)
      endif
      if (logy) then
        if (b .le. 0 .or. t .le. 0) goto 99
        plub = dlog(b)
        plut = dlog(t)
      endif
      error = .false.

C --- Set initial frame to max size if uuytox spec'd and noclip ---
      if (uuytox .ne. 0 .and. noclip) then
        if (iprint() .ge. 50) then
          call awrit1(' pltstu: resizing frame in x to satisfy '//
     .    'dy/dx(UU) = %1;4g*dy/dx(MU).',' ',80,i1mach(2),uuytox)
          print '(10x,''Frame initially set at ...'')'
        endif
        if (.not. noclip) call pltstp(plpl,plpr,plpb,plpt)
        if (noclip) call pltstp(plgl,plgr,plgb,plgt)
        if (iprint() .ge. 50) print '(10x,''and reset to ...'')'
      endif

C --- Set the scale and offsets in this frame ---
C Calculate by equating points at corners, eg
C scalxu (plul + offxu) = scalxg (plpl + offxg)
      scalxu = scalxg*(plpr-plpl)/(plur-plul)
      scalyu = scalyg*(plpt-plpb)/(plut-plub)
      offxu = (((plpl+plpr)/2 + offxg)*scalxg/scalxu - (plul+plur)/2)
      offyu = (((plpt+plpb)/2 + offyg)*scalyg/scalyu - (plut+plub)/2)
c      if (lrot) then
c        scalxu = scalxg*(plpr-plpl)/(t-b)
c        scalyu = scalyg*(plpt-plpb)/(l-r)
c        offxu = (((plpl+plpr)/2 + offxg)*scalxg/scalxu - (t+b)/2)
c        offyu = (((plpt+plpb)/2 + offyg)*scalyg/scalyu - (l+r)/2)
c      endif

C --- Clip frame to largest permitted by uuytox ---
      if (uuytox .ne. 0) then
        if (scalxu .gt. scalyu*strtch*uuytox) then
          scalxu = scalyu*uuytox
          xx1 = (plur-plul)*scalxu/scalxg/2
          xx2 = (plpr+plpl)/2
          call pltstp(xx2-xx1,xx2+xx1,plpb,plpt)
        else
          scalyu = scalxu/uuytox
          xx1 = (plut-plub)*scalyu/scalyg/2
          xx2 = (plpt+plpb)/2
          call pltstp(plpl,plpr,xx2-xx1,xx2+xx1)
        endif
        scalxu = scalxg*(plpr-plpl)/(plur-plul)
        scalyu = scalyg*(plpt-plpb)/(plut-plub)
        offxu = (((plpl+plpr)/2 + offxg)*scalxg/scalxu - (plul+plur)/2)
        offyu = (((plpt+plpb)/2 + offyg)*scalyg/scalyu - (plut+plub)/2)
      endif

   99 continue
      if (iprint() .ge. 50 .or. error) then
        s = ' '
        is = 0
        call bin2a(' pltstu: frame corners (',1,0,0,1,0,scrwid,s,is)
        if (logx) call bin2a('logx',1,0,0,1,0,scrwid,s,is)
        if (logy) call bin2a('logy',1,0,0,1,0,scrwid,s,is)
        call bin2a('UU ): lb = (%1;4g %1;4g)',1,0,0,1,0,scrwid,s,is)
        call bin2a('  rt = (%1;4g %1;4g)',1,0,0,1,0,scrwid,s,is)
        call awrit4(s(1:is),' ',scrwid,i1mach(2),l,b,r,t)
        if (error) call rx(
     .  'pltstu (abort): log scale incompatible with negative limit')
        if (iprint() .ge. 100)
     .  call awrit4(' pltstu: scal(x,y)u = %1;4g %1;4g  off(x,y)u = '//
     .  '(%1;4g %1;4g)',' ',scrwid,i1mach(2),scalxu,scalyu,offxu,offyu)
      endif

      if (iprint() .ge. 100) call shosyv(0,0,0,i1mach(2))

      end
      subroutine pltstp(l,r,b,t)
C- Bound area of current plot, in GU; redefine UU at these corners
C ----------------------------------------------------------------
Ci Inputs
Ci   l,r,b,t: left, right, bottom, top
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision l,r,b,t
      integer i

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz
      integer iprint,i1mach

      call lodsyv('gl',1,l,i)
      call lodsyv('gr',1,r,i)
      call lodsyv('gb',1,b,i)
      call lodsyv('gt',1,t,i)
      call lodsyv('ml',1,(offxg+l)*scalxg,i)
      call lodsyv('mr',1,(offxg+r)*scalxg,i)
      call lodsyv('mb',1,(offyg+b)*scalyg,i)
      call lodsyv('mt',1,(offyg+t)*scalyg,i)
      plpl = l
      plpr = r
      plpb = b
      plpt = t
      noclip = .false.
      if (iprint() .ge. 30) then
        call awrit4(' pltstp:  frame corners (GU): lb = (%d %d)'//
     .  '  rt = (%d %d)',' ',80,i1mach(2),plpl,plpb,plpr,plpt)
        if (medium .eq. 3)
     .  call awrit4('%% pltstp: frame corners (GU): lb= (%d %d)'//
     .  '  rt = (%d %d)',' ',80,stream,plpl,plpb,plpr,plpt)
      endif

C ... Update plot bounds in user units
C     call pltstu(plul,plur,plub,plut)

      end
      subroutine pltsts(plogx,plogy,io)
C- Set scale (linear or log)
C ----------------------------------------------------------------
Ci Inputs
Ci   plogx,plogy,io
Cr Remarks
Cr   Sets logx,logy from ilogx or ilogy (io=1=> set ilog* from log*)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer io
      logical plogx,plogy

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz
      integer iprint,i1mach

      if (io .eq. 0) then
        logx = plogx
        logy = plogy
        if (iprint() .ge. 30) then
          call awrit2(' pltsts:  set  logx=%l  logy=%l',
     .    ' ',80,i1mach(2),logx,logy)
          if (medium .eq. 3) call awrit2('%% pltsts:  set  logx=%l'//
     .    '  logy=%l',' ',80,stream,logx,logy)
        endif
      else
        plogx = logx
        plogy = logy
      endif

      end
      subroutine plthet(alpha)
C- Set angle for nonorthogonal shear of Cartesian axes
C ----------------------------------------------------------------
Ci Inputs
Ci   alpha : shearing angle for linear transformation of (x,y)
Cr Remarks
Cr   By default, the axes are orthogonal; their angle is pi/2.
Cr   For other values of alpha, (x,y) are sheared by
Cr   the linear transformation:
Cr       x->x
Cr       y->y*sin(alpha) + x*cos(alpha)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision alpha

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz


      cosa = cos(alpha)
      sina = sin(alpha)

      end
      subroutine pltam(x,y,opt,cs)
C- Plot absolute position, medium's units
C ----------------------------------------------------------------
Ci Inputs
Ci   x,y
Ci   opt: -1, return the current pen position
Ci         0, lift pen before move
Ci         1, drop pen before move
Ci         2, drop pen before move, lift after
Ci   cs:   passed to plcrve (for PostScript)
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x,y,cs(0:3)
      integer opt

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local variables
      integer iprint,is
      double precision d,ds,xs,ys,cost,sint,d1mach
      logical tog
      character s*100

      is = 0
      s = ' '
      if (opt .eq. -1) goto 5
      if (opt .ge. 1) goto 2
      goto 3

C --- Draw all line types, Postscript ---
    2 continue
      if (medium .eq. 3) then
        if (pstat .eq. 0) then
          call pbbox(xcm,ycm)
          call awrit4('newpath %1;nd %1;nd moveto',s,100,0,
     .    ndwrt,xcm,ndwrt,ycm)
C          call bin2a('newpath',0,0,0,1,0,80,s,is)
C          call bin2a('(f8.1)',1,0,xcm,4,0,100,s,is)
C          call bin2a('(f8.1)',1,0,ycm,4,0,100,s,is)
C          call bin2a('moveto',1,0,0,1,0,80,s,is)
        endif
C        call bin2a('(f8.1)',1,0,x,4,0,100,s,is)
C        call bin2a('(f8.1)',1,0,y,4,0,100,s,is)
C        call bin2a('lineto',1,0,0,1,0,80,s,is)
        call pbbox(x,y)
        call awrit4('%a %1;nd %1;nd lto',s,100,0,ndwrt,x,ndwrt,y)
        call skpblb(s,100,is)
        is = is+1
        pstat = 1
        goto 5
      endif
      pstat = 1

      goto (10,20) lintyp
      goto 3

C --- Line type 1 Draw, ASCII  ---
   10 continue
      if (medium .eq. 1) then
        write(stream,110) x, y
        d = dsqrt((x-xcm)**2 + (y-ycm)**2) + d1mach(1)
        cost = (x-xcm)/d
        sint = (y-ycm)/d
        if (ltbld .ge. 1) then
          write(stream,310) xcm-xb*sint/1.75d0, ycm+yb*cost/1.75d0
          write(stream,110) x-xb*sint/1.75d0, y+yb*cost/1.75d0
          write(stream,310) xcm+xb*sint/1.75d0, ycm-yb*cost/1.75d0
          write(stream,110) x+xb*sint/1.75d0, y-yb*cost/1.75d0
          if (ltbld .eq. 1) write(stream,310) x, y
        elseif (ltbld .ge. 2) then
          write(stream,310) xcm-xb*sint, ycm+yb*cost
          write(stream,110) x-xb*sint, y+yb*cost
          write(stream,310) xcm+xb*sint, ycm-yb*cost
          write(stream,110) x+xb*sint, y-yb*cost
          write(stream,310) x, y
        endif
  110   format('PA ',2f14.5)
      else
        call rx('PLTAM: unknown medium')
      endif
      goto 5

C --- Line type 2 draw, ASCII ---
   20 continue
      d = dsqrt((x-xcm)**2 + (y-ycm)**2) + d1mach(1)
      cost = (x-xcm)/d
      sint = (y-ycm)/d
      tog = .true.
C Move pen up when tog false, down when tog true
   25 if (tog) then
        ds = min(max(0d0,linl1-linl2-lnused),d)
      else
        ds = min(max(0d0,linl1-lnused),d)
      endif
      xs = ds*cost + xcm
      ys = ds*sint + ycm
C ... Pen down motion
      if (tog .and. ds .gt. 0) then
        if (medium .eq. 1) then
          write(stream,110) xs, ys
          if (ltbld .ge. 1) then
            write(stream,310) xcm-xb*sint/1.75d0, ycm+yb*cost/1.75d0
            write(stream,110) xs-xb*sint/1.75d0, ys+yb*cost/1.75d0
            write(stream,310) xcm+xb*sint/1.75d0, ycm-yb*cost/1.75d0
            write(stream,110) xs+xb*sint/1.75d0, ys-yb*cost/1.75d0
            write(stream,310) xs, ys
          elseif (ltbld .ge. 2) then
            write(stream,310) xcm-xb*sint, ycm+yb*cost
            write(stream,110) xs-xb*sint, ys+yb*cost
            write(stream,310) xcm+xb*sint, ycm-yb*cost
            write(stream,110) xs+xb*sint, ys-yb*cost
            write(stream,310) xs, ys
          endif
        else
          call rx('PLTAM: unknown medium')
        endif
C ... Pen up motion
      else
        if (medium .eq. 1) then
          write(stream,310) xs, ys
        endif
      endif
      xcm = xs
      ycm = ys
      d = d - ds
      lnused = mod(lnused+ds,linl1)
      tog = .not. tog
      if (d .gt. 0) goto 25
      goto 5

C --- Move (pen up) ---
    3 continue
C ... Close previously opened curve
      call plcrve(medium,pstat,-2d0,s,is)
      lnused = 0
      if (medium .eq. 1) then
        write(s,310) x, y
  310   format('MA ',2f14.5)
        is = 31
      elseif (medium .eq. 3) then
      else
        call rx('PLTAM: unknown medium')
      endif
      goto 5

C --- Exit ---
    5 continue
      if (opt .eq. -1) then
        x = xcm
        y = ycm
      else
        xcm = x
        ycm = y
        if (opt .eq. 2) call plcrve(medium,pstat,cs,s,is)
        if (is .gt. 0) write(stream,'(a)') s(1:is)
      endif

    6 continue
      if (iprint() .ge. 120)
     .print 333, x, y, opt, lintyp,linl1,linl2,linl3,linl4
  333 format(' Pltam: x,y=', 2f10.3, '  opt=',i1,
     .'  lintyp=',i2,4f7.3)

      end
      subroutine mve(x,y)
C- Plot absolute position, user's units
C ----------------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   Nonorthogonal shear transformation.  General linear transformation:
Cr
Cr      (xs)   (a  b) (x)      (x)
Cr      (  ) = (    ) ( ) =  R ( )
Cr      (ys)   (c  d) (y)      (y)
Cr
Cr   We seek special shear with the following properties:
Cr
Cr         (1)   (1)         |   (1) |             |   (0) |
Cr       R ( ) = ( )   and   | R ( ) | = 1   and   | R ( ) | = 1
Cr         (0)   (0)         |   (0) |             |   (1) |
Cr
Cr   Also the the angle between xs,ys should be cos(theta). Then
Cr      a=1, c=0, b=cos(theta), d=sin(theta),   or
Cr
Cr      (xs)   (1  cos(theta)) (x)
Cr      (  ) = (             ) ( )
Cr      (ys)   (0  sin(theta)) (y)
Cr
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x,y,cs0(0:3)
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local variables
      integer io,iprint,i1mach
      double precision cs(0:3),xs,ys

      cs(0) = -2d0
      io = 0
      goto 10

      entry drwe(x,y,cs0)
      call dcopy(4,cs0,1,cs,1)
      io = 2
      if (cs(0) .eq. 10) cs(0) = 5
      if (cs(0) .eq. 11) cs(0) = 7
      goto 10

      entry drw(x,y)
      io = 1
   10 continue
      xcu = x
      ycu = y
      goto 30

      entry mver(x,y)
      io = 0
      goto 20

      entry drwr(x,y)
      io = 1
      goto 20

      entry drwre(x,y,cs0)
      call dcopy(4,cs0,1,cs,1)
      io = 2
      if (cs(0) .eq. 10) cs(0) = 5
      if (cs(0) .eq. 11) cs(0) = 7

   20 continue
      if (logx) then
        xcu = xcu*dexp(x)
      else
        xcu = xcu+x
      endif
      if (logy) then
        ycu = ycu*dexp(y)
      else
        ycu = ycu+y
      endif

   30 continue
      xs = xcu
      if (logx) then
        if (xcu .le. 0)
     .  call rx1('x=%g encountered, but log x scale',xcu)
        xs = dlog(xcu)
      endif
      ys = ycu
      if (logy) then
        if (ycu .le. 0)
     .  call rx1('y=%g encountered, but log y scale',ycu)
        ys = dlog(ycu)
      endif
C  ... Nonorthogonal transformation
      xs = (offxu+xs)*scalxu + cosa*(offyu+ys)*scalyu
      ys = sina*(offyu+ys)*scalyu

      if (iprint() .ge. 110) then
        call awrit5(' mve (opt=%i): %1;4g %1;4g(UU)  %d %d(MU)',
     .  ' ',80,i1mach(2),io,xcu,ycu,xs,ys)
      endif

      call pltam(xs,ys,io,cs)
C      if (lrot) then
C        call pltam(ys,xs,io,cs)
C      else
C        call pltam(xs,ys,io,cs)
C      endif
      return

      entry pposu(x,y)
C- Return, update current pen pos'n by scaling from the medium's units
      xcu = xcm/scalxu - offxu
      ycu = ycm/scalyu - offyu
      if (logx) xcu = dexp(xcu)
      if (logy) ycu = dexp(ycu)
      x = xcu
      y = ycu
      end
      subroutine mveg(x,y)
C- Plot absolute position, graphics units
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision x,y
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local variables
      integer io
      double precision xcg,ycg,cs0

      io = 0
      goto 10

      entry drwg(x,y)
      io = 1
   10 continue
      xcg = x
      ycg = y
      goto 30

      entry mvegr(x,y)
      io = 0
      goto 20

      entry drwgr(x,y)
      io = 1
   20 continue
      call pposg(xcg,ycg)
      xcg = x+xcg
      ycg = y+ycg

   30 continue
      cs0 = -2d0
      call pltam((offxg+xcg)*scalxg,(offyg+ycg)*scalyg,io,cs0)
C the following call updates xcu,ycu
      call pposu(xcg,ycg)
      end
      subroutine pposg(x,y)
C- Return absolute position, graphics units
C Passed parameters
      double precision x,y
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz
C- Return current pen position by scaling from the medium's units
      x = xcm/scalxg - offxg
      y = ycm/scalyg - offyg
      end
      subroutine cdrw(x,y,xold,yold,lclip,box,cs)
C- Draw with explicit clips of boundaries
C ----------------------------------------------------------------------
Ci Inputs
Ci   lx     :T if x is on a log scale
Ci   ly     :T if y is on a log scale
Ci   x,y    :draw line to this point
Ci   xold,yold :draw line from this point
Ci   box    :bounding box:
Ci          :xmin=box(1),ymin=box(2),xmax=box(3),ymax=box(4)
Ci   lclip  :T clip line segments outside bounding box
Ci   cs     :if cs(0) = 0, has no effect
Ci          :otherwise, substitute drwe for drw, passing cs
Co Outputs
Cr Remarks
Cb Bugs
Cb   not implemented for log scale
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lclip
      double precision x,y,xold,yold,box(4),cs(0:3)
C ... Local parameters
      logical le1,le2
      integer n
      double precision xe1,ye1,xe2,ye2

      if (lclip) then
        call gtxeff(xold,yold,x,y,box,xe1,ye1,xe2,ye2,le1,le2,n)
C       Both points within clip region
        if (le1 .and. le2) then
          if (cs(0) .eq. 0) then
            call drw(x,y)
          else
            call drwe(x,y,cs)
          endif
C       Only first point within clip region
        elseif (le1) then
          call drw(xe1,ye1)
          call mve(x,y)
C       First point outside, second point on the boundary
        elseif (le2 .and. n .eq. 0) then
          if (cs(0) .eq. 0) then
            call drw(x,y)
          else
            call drwe(x,y,cs)
          endif
C       Only second point within clip region
        elseif (le2) then
          call mve(xe1,ye1)
          if (cs(0) .eq. 0) then
            call drw(x,y)
          else
            call drwe(x,y,cs)
          endif
C       Neither point within clip region
        else
          if (n .eq. 0) goto 99
          call mve(xe1,ye1)
          if (cs(0) .eq. 0) then
            call drw(xe2,ye2)
          else
            call drwe(xe2,ye2,cs)
          endif
          call mve(x,y)
        endif

C     Neither point within clip region
      else
        if (cs(0) .eq. 0) then
          call drw(x,y)
        else
          call drwe(x,y,cs)
        endif
      endif

C     Exit, setting xold,yold to current x,y
   99 continue
      xold = x
      yold = y
      end
      subroutine gtxeff(x1,y1,x2,y2,box,xe1,ye1,xe2,ye2,le1,le2,n)
C- Given a line segment, return a segment lying within a clipped region
C ----------------------------------------------------------------------
Ci Inputs
Ci   lx     :T if x is on a log scale
Ci   ly     :T if y is on a log scale
Ci   x1,y1  :first x-y pair
Ci   x2,y2  :second x-y pair
Ci   box    :bounding box:
Ci          :xmin=box(1),ymin=box(2),xmax=box(3),ymax=box(4)
Ci   xe1,ye1:first point of subsegment
Ci   xe2,ye2:second point of subsegment
Co Outputs
Co   le1    :T if xe1,ye1 and x1,y1 are the same point
Co   le2    :T if xe2,ye2 and x2,y2 are the same point
Co    n     :number of boundary points
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical le1,le2
      double precision box(4),x1,y1,x2,y2,xe1,ye1,xe2,ye2
C ... Local parameters
      logical inside,betw
      integer n,nb
      double precision x,y,fuzz,alpha,beta,xx
      parameter (fuzz=1d-7)
      inside(x,y) = x+fuzz.ge.box(1) .and. x-fuzz.le.box(3) .and.
     .y+fuzz.ge.box(2) .and. y-fuzz.le.box(4)
      betw(x,alpha,beta) =
     .x .ge. min(alpha,beta) .and. x .le. max(alpha,beta)

C     Sanity check
      if (box(3) .le. box(1)) call rx('gtxeff: bad box')
      if (box(4) .le. box(2)) call rx('gtxeff: bad box')

      le1 = inside(x1,y1)
      le2 = inside(x2,y2)
      xe1 = x1
      xe2 = x2
      ye1 = y1
      ye2 = y2
      if (le1 .and. le2) then
        return
      elseif (le1 .or. le2) then
        call gtxff2(1,x1,x2,y1,y2,box,xe1,ye1,n,nb)
        if (n .ne. 1 .and. nb .eq. 0) then
          call rx('bug in gtxeff')
        endif
      else
        call gtxff2(1,x1,x2,y1,y2,box,xe1,ye1,n,nb)
        call gtxff2(2,x1,x2,y1,y2,box,xe2,ye2,n,nb)
        if (n .eq. 0) return
        if (n .ne. 2) call rx('bug in gtxeff')
        if (betw(xe2,x1,xe1)) then
          xx = xe2
          xe2 = xe1
          xe1 = xx
          xx = ye2
          ye2 = ye1
          ye1 = xx
        endif
      endif

      end
      subroutine gtxff2(n,x1,x2,y1,y2,box,xe,ye,k,kb)
C- Return nth occurrence of intersection of box w/ line seg.
C ----------------------------------------------------------------------
Ci Inputs
Ci   n      :nth occurence to seek.  Should be 1 or 2.
Ci   x1,y1  :first x-y pair of line segment
Ci   x2,y2  :second x-y pair of line  segment
Ci   box    :box
Co Outputs
Co   xe     :x coordinate of last occurence
Co   ye     :y coordinate of last occurence
Co   k      :last occurence found.  If successful, k=n. See Remarks
Co   kb     :number of points on the boundary. See Remarks
Cr Remarks
Cr   A point on the box boundary will be returned in (xe,ye),
Cr   but k will not reflect this.  Caller can check by evaluating
Cr   distance between (xe,ye) and box.
Cu Updates
Cu   Added argument kb
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,k,kb
      double precision x1,x2,y1,y2,xe,ye,box(4)
C ... Local parameters
      logical betw
      integer m
      double precision isec,xa,dy,sloi,fuzz,alpha,beta,x,xmin,xmax
      parameter (fuzz=1d-4)
      isec(sloi,xa,dy) = xa + dy*sloi
C     betw => x between alpha and beta and between xmin and xmax
      betw(x,alpha,beta,xmin,xmax) =
     .x .ge. min(alpha,beta) .and. x .le. max(alpha,beta) .and.
     .x .ge. xmin .and. x .le. xmax

      k = 0
      kb = 0

C     Intersection of two horizontals with line segment
      if (y2 .ne. y1) then
        do  m = 2, 4, 2
          ye = box(m)
          xe = isec((x2-x1)/(y2-y1),x1,ye-y1)
          if (betw(xe,x1,x2,box(1),box(3))) then
C         Exclude points coincident with (x1,y1) or (x2,y2)
            if (min(dabs(xe-x1)+dabs(ye-y1),dabs(xe-x2)+dabs(ye-y2))
     .      .gt. fuzz) then
              k = k+1
              if (k .eq. n) return
            else
              kb = kb+1
            endif
          endif
        enddo
      endif

C     Intersection of two verticals with line segment
      if (x2 .ne. x1) then
        do  m = 1, 3, 2
          xe = box(m)
          ye = isec((y2-y1)/(x2-x1),y1,xe-x1)
          if (betw(ye,y1,y2,box(2),box(4))) then
C         Exclude points coincident with (x1,y1) or (x2,y2)
            if (min(dabs(xe-x1)+dabs(ye-y1),dabs(xe-x2)+dabs(ye-y2))
     .      .gt. fuzz) then
              k = k+1
              if (k .eq. n) return
            else
              kb = kb+1
            endif
          endif
        enddo
      endif

C     call rx('bug in gtxeff')

      end
      subroutine plwstr(string,mode)
C- Write string directly to output stream
Ci mode=1: string is name of a file. Insert contents of the file
      implicit none
      character*(*) string
      integer mode
C Local
      integer ifi,fopng,i,iprint,i1mach
      character*1024 s
      logical rdstrn
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

      if (mode .eq. 1) then
        ifi = fopng(string,-1,1)
        s(1:60) = string
        if (iprint() .ge. 40) call awrit0(' plwstr:  '//
     .  'inserting file "'//s(1:60)//'%a"',' ',80,i1mach(2))
        call awrit0('%N%% --- plwstr:  inserting file "'//
     .  s(1:60)//'%a"',' ',80,stream)
   10   if (rdstrn(ifi,s,len(s),.false.)) then
          call skpblb(s,len(s),i)
          write (stream,'(a)') s(1:i+1)
          goto 10
        endif
        call fclose(ifi)
      else
        s(1:60) = string
        if (iprint() .ge. 40) call awrit0(' plwstr:  '//
     .  'inserting string "'//s(1:60)//'%a"',' ',80,i1mach(2))
        call awrit0('%N%% --- plwstr:  inserting string "'//
     .  s(1:60)//'%a"',' ',80,stream)
        write (stream,'(a)') string
      endif

      end
      subroutine pltdmp(io)
C- Dumps plot file
C ----------------------------------------------------------------
Ci Inputs
Ci   io:  0, dump plot immediately; 1, wait for carriage return
Co Outputs
Cr Remarks
Cr  This is only meaningful if plot is dumped to device at one time.
Cr  Hercules card: restore screen to ASCII mode when plot is complete.
Cr  When applicable and io=1, prompt for a carriage return first.
Cr  PostScript file: write out showpage
C ----------------------------------------------------------------
      implicit none
      integer io
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz
      if (medium .eq. 1) then
      elseif (medium .eq. 2) then
      elseif (medium .eq. 3) then
        call awrit4('%%%%BoundingBox: %;1d %;1d %;1d %;1d',' ',
     .  80,stream,bbox(1),bbox(3),bbox(2),bbox(4))
        write (stream,10)
   10   format('showpage')
      endif
C     if (isopen(stream,.false.)) call fclose(stream)
      call fclose(stream)

      end
      subroutine pltcol(color)
C- Set pen color of current medium
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer color
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

      pcol = color
      end
      subroutine plndpi(n)
C- Set number of dots-per-inch
      implicit none
      integer n
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz


      ndpi = n
      ndwrt = int(dlog10(ndpi*2.54d0*scalxc)) + 1

      end

      subroutine plntyp(type,bold,l1,l2,l3,l4)
C- Set line type
C ----------------------------------------------------------------
Ci Inputs
Ci   type: 0, move, pen up
Ci         1, draws simple lines
Ci         2, draws broken lines, pen down for length (l1-l2) followed
Ci            by pen up for length l2.
Cr   l1:   repeat length of the line type.  Scale is
Cr         hardware dependent but defined so that 1 corresponds
Cr         to a visually pleasing period.  Scale can be changed
Cr         in pltini.
Cr   l2:   length for which pen is up during a draw.
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer type,bold
      double precision l1,l2,l3,l4

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      logical err
      integer is
      character s*100
      integer lintys,ltblds
      double precision linl1s,linl2s,linl3s,linl4s
      save lintys,ltblds,linl1s,linl2s,linl3s,linl4s

C --- Scale linetype lengths from "sensible" scales (numbers of order 1) to
C     medium units ---
      lintyp = type
      ltbld  = bold
      linl1 = l1*lnprd
      linl2 = l2*lnprd
      linl3 = l3*lnprd
      linl4 = l4*lnprd
      lnused = 0

   10 if (medium .ne. 3) return
      s = ' '
      is = 0
      if (lintyp .gt. 0) then
        call bin2a('d2',1,0,dble(ltbld)/2,4,0,100,s,is)
        call bin2a('setlinewidth 1 setlinejoin [',1,0,0,1,0,80,s,is)
        if (lintyp .eq. 1) then
          is = is+1
        elseif (lintyp .eq. 2) then
          call bin2a('d2',0,0,linl1,4,0,100,s,is)
          call bin2a('d2',1,0,linl2,4,0,100,s,is)
          if (linl3 .ne. 0) then
            call bin2a('d2',1,0,linl3,4,0,100,s,is)
            call bin2a('d2',1,0,linl4,4,0,100,s,is)
          endif
        elseif (lintyp .eq. 3) then
          is = is-1
          call bin2a('1 setlinecap [0',0,0,0,1,0,80,s,is)
          call bin2a(' ',1,0,ltbld,2,0,100,s,is)
        endif
        call bin2a('] 0 setdash',0,0,0,1,0,80,s,is)
      endif
      if (is .gt. 0) write (stream,*) s(1:is)
      return

      entry plntyr(type,bold,l1,l2,l3,l4,err)
C- Retrieves current line type
      err = lnused .ne. 0
      if (err) return
      type = lintyp
      bold = ltbld
      l1 = linl1/lnprd
      l2 = linl2/lnprd
      l3 = linl3/lnprd
      l4 = linl4/lnprd
      return

      entry plntys
C- Saves current line type
      lintys = lintyp
      ltblds  = ltbld
      linl1s = linl1
      linl2s = linl2
      linl3s = linl3
      linl4s = linl4
      return

      entry plntyg
C- Restores saved line type
      lintyp = lintys
      ltbld  = ltblds
      linl1 = linl1s
      linl2 = linl2s
      linl3 = linl3s
      linl4 = linl4s
      goto 10

      end
      subroutine ticset(llog,mode,mt0,bot,top,ts0,tx0,ch,lmode,
     .tic,mt,mt1,itic,s)
C- Set a tic pattern
C ----------------------------------------------------------------------
Ci Inputs
Ci   llog : T if log scale, F if linear scale
Ci   mode : default tic pattern.  Actual tic pattern is returned in
Ci        : lmode, which see for description.
Ci        : mode>=0 => lmode is assigned to mode
Ci        : mode<0 => ticset chooses a default lmode
Ci   mt0  : default number of tics per major tic.  Actual number is
Ci        : returned in mt, which see for description.
Ci        : mt0<0 => ticset chooses a value for mt.
Ci bot,top: axis bounds
Ci   ts0  : (mode=0 only) default tic spacing.  Uses only ts0(1).
Ci        : ts0(1)<0 => ticset chooses ts
Ci        : (mode>10) vector of positions at which to set tic marks
Ci   tx0  :force tic through tx0
Ci        :tx0 < bot => ticset chooses a value
Ci   ch   :a character of length 1, used for printout
Co Outputs
Co  lmode :specifies tic pattern; see input 'mode' for description
Co        : mode=0: uniformly spaced tics; this is the
Co        :         usual mode for linear scale, but may be
Co        :         appropriate for log scale also.
Co        : ... modes 1..3 are appropriate for log scale:
Co        : mode=1: tics at 1,2,3,4,5,6,7,8,9 * 10^integer
Co        : mode=2: tics at 1,2,5 * 10^integer
Co        : mode=3: tics 1 * 10^integer
Co        :
Co        : mode>10: use supplied vector of tic marks, passed
Co        :        : in array ts0 (mode-10 points)
Co   tic  : position of 1st tic
Co   mt   : number of of tics between major tics, i.e. for tic mark i,
Co        : major tic at any (i+mt1 modulus mt) = 0
Co        : minor tic at any (i+mt1 modulus mt) > 0
Co   mt1  : offset for modulus of mt at 1st tic; see description for mt
Co   ts   : lmode=0, lmode>10; spacing between first and second tics
Co        : other values of lmode: not set (ts has no meaning)
Cr  Remarks
Cr   *tic marks are of two types, major and minor.
Cr    mt and mt1 specify of which type each tic mark belongs.
Cr   *The spacing between tic marks depends on lmode.
Cr    tic is the position of the first tic mark.
Cr    ts specifies spacing between 1st and 2nd mark (lmode=0,10)
Cr    For the log scale modes (lmode=1..3), see description of lmode
Cr    ts0 is a list of all marks (lmode=10)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical llog
      character ch*1, s*(*)
      double precision bot,top,ts0(*),tx0,ts,tic
      integer mode,mt0,mt,mt1,lmode,itic
C ... Local parameters
      double precision botl,topl,l10,width,ticspc,xx,fuzz,fuzfac
      parameter (fuzfac=1d-10)
      integer j,k

C --- Make lmode ---
      lmode = mode
      if (lmode .lt. 0 .and. .not. llog) lmode = 0
      if (lmode .lt. 0 .and. llog) then
        botl = dlog(bot)
        topl = dlog(top)
        l10  = dlog(10d0)
        width = dabs(topl-botl)/l10
        lmode = 3
        if (width .lt. 3d0)  lmode = 2
        if (width .lt. 1d0)  lmode = 1
        if (width .lt. .5d0) lmode = 0
      endif

C --- Make mt ---
      mt = mt0
      if (mt .eq. -1) then
        mt = 1
        if (lmode .eq. 0) mt = 2
        if (lmode .eq. 2) mt = 3
      endif

C     Shift bot,top by fuzz factor
      botl = min(top,bot)
      topl = max(top,bot)
      if (llog) then
        fuzz = (topl/botl)*fuzfac
        botl = botl*(1-fuzz)
        topl = topl*(1+fuzz)
      else
        fuzz = (topl - botl)*fuzfac
        botl = botl - fuzfac
        topl = topl + fuzfac
      endif

C --- ts, tic and mt1 for mode=0 ---
      if (lmode .eq. 0) then
        ts = ts0(1)
        if (ts0(1) .le. 0) ts = ticspc(dabs(top-bot)/10,1)
        if (top .lt. bot) ts = -ts
        tic = tx0
        mt1 = 0
        if (tic .lt. botl) then
          j = int(botl/ts)
          if (j*ts .lt. botl) j = j+1
          tic = j*ts
          mt1 = mod(abs(j),mt)
        else
   11     continue
          if (tic-ts .ge. botl) then
            tic = tic-ts
            mt1 = mt1+1
            goto 11
          endif
          mt1 = mod(mt1,mt)
        endif
        itic = 1
C ...   For printout
        call awrit3('%10z%a  '//ch//'t1=%1;4g ts'//ch//'=%1;4g mt'//
     .  ch//'=%i',s,len(s),0,tic,ts,mt)
        ts0(1) = ts
        return
      endif

C --- ts, tic and mt1 for mode>10 ---
      if (lmode .gt. 10) then
        tic = ts0(1)
        if (mode .eq. 11) then
          ts = 0
        else
          ts = ts0(2) - ts0(1)
        endif
        itic = 1
        k = 1
        if (tx0 .lt. botl .or. tx0 .gt. topl) goto 22
C ...   Determine mt1 by finding point closest to tx0
        xx = dabs(tx0-ts0(1))
        do  20  j = 1, mode-10
          if (dabs(tx0-ts0(j)) .lt. xx) then
            k = j
            xx = dabs(tx0-ts0(j))
          endif
   20   continue
   22   continue
        mt1 = mod(k-1,mt)
        call awrit5('%10z%a  '//ch//'t1=%1;4g mode=%i(%i tics) mt'//
     .  ch//'=%i:%i',s,len(s),0,tic,5,lmode-10,mt,mt1)
        return
      endif

C --- ts, tic, mt1 and itic for mode=1,2,3 ---
      if (lmode .ge. 1 .and. lmode .le. 3) then
        ts = 0
        tic = dmax1(tx0,botl)
        itic = 0
        call ticnxt(lmode,tic,topl,0d0,itic)
        if (tic .lt. botl) call ticnxt(lmode,tic,topl,0d0,itic)
        mt1 = mod(mt-itic+1,mt)
        call awrit3('%10z%a  '//ch//'t1=%1;4g mode=%i mt'//
     .  ch//'=%i',s,len(s),0,tic,lmode,mt)
        return
      endif

C      print *, 'ticset:', j,mt1
C      pause
      call rx('ticset: mode not defined')

      end
      subroutine ticnxt(mode,tic,top,ts0,inew)
C- Finds next tic, depending on mode (see ticset)
Ci mode; tic: last tic; ts: tic spacing (mode=0)
Co inew,tic: next tic, and index to next tic.
      implicit none
      integer mode
      double precision tic,top,ts0(1),d1mach
      integer is,i,inew,nx(3),iprint
      double precision tx(11,3),scale,ss,err
      save scale
      data tx /1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,
     .1d0,2d0,5d0,10d0,20d0,0d0,0d0,0d0,0d0,0d0,0d0,
     .1d0,10d0,100d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/,
     .nx /10,4,2/

      if (mode .eq. 0) then
        tic = tic + ts0(1)
        inew = inew+1
        if (tic .gt. top) inew = 0
        goto 100
      endif

      if (mode .gt. 10) then
        inew = inew+1
        if (inew .gt. mode-10) then
          inew = 0
        else
          tic = ts0(inew)
        endif
        goto 100
      endif

C ... for the logarithmic modes, reset tic to closest table value
      if (inew .eq. 0) then
        scale = nint(dlog10(dabs(tic)))
        if (10**scale .gt. dabs(tic)) scale = scale-1
        ss = dabs(tic)/10**scale

        err = 9d9
        do  10  i = 1, nx(mode)
          if (dabs(ss/tx(i,mode)-1) .gt. err) goto 10
          is = i
          err = dabs(ss/tx(i,mode)-1)
   10   continue
        inew = mod(is-1,nx(mode)+1)
      endif
      inew = inew+1
      if (inew .gt. nx(mode)) then
        inew = 2
        scale = scale+1
      endif
      tic = tx(inew,mode)*10**scale
      if (tic .gt. top*(1+4*d1mach(3))) inew = 0

  100 continue
      if (iprint() .gt. 90) print 333, tic, mode, inew
  333 format(' ticnxt:  tic=',f9.6,' mode=',i2,' inew=',i2)

      end
      function ticspc(width,n)
C- Procedure determines tic spacing given a width, linear axes
C ----------------------------------------------------------------
Ci Inputs
Ci   width
Ci   n:  not used for now
Co Outputs
Co   ticspc
Cr Remarks
Cr
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      double precision width,ticspc
      integer n
C Local parameters
      double precision scale,tmp
C This construction avoids taking integer part of negative number
      scale = dlog10(width)
      if (scale .gt. 0) then
        scale = int(scale)
      else
        scale = -int(-scale) - 1
      endif
      scale = 10d0**(-scale)
C tmp is the width scaled to a number between 1 and 10
      tmp = width*scale
      if (tmp .gt. 7.5d0) then
        tmp = 10
      elseif (tmp .gt. 3.5d0) then
        tmp = 5
      elseif (tmp .gt. 1.5d0) then
        tmp = 2.5d0
      else
        tmp = 1
      endif
      ticspc = int(tmp + .00001d0)/scale
      end
      subroutine frme(tx0,tsx,mtx,rmtx,rmntx,modx,xnum,xnfmt,xlabel,
     .ty0,tsy,mty,rmty,rmnty,mody,ynum,ynfmt,ylabel,
     .sfxdef,np,xp,fxdef,title,lt,cs)
C- draw a frame around the plot and either tic marks or grid lines.
C ----------------------------------------------------------------
Ci Inputs
Ci   ... the following apply to the x-axis
Ci   tx0:  Force a major tic through this point
Ci         tx0 < plul makes default
Ci   tsx:  spacing in user's units between tic marks along the x axis
Ci         tsx <= 0 makes default of approximately 1/10 grid width
Ci   mtx:  number of major tics per tic, x axis
Ci   rmtx: size of major tic, in proportion to plut-plub (.03 is good)
Ci   rmntx:size of minor tic, in proportion to major (.6 is good).
Ci   modx :tic mode (see ticset)
Ci   xnum :1  to number axis on TOP
Ci        :0  to number axis on 
Ci        :-1 to number axis on BOTTOM
Ci   xnfmt:awrite format string for x-axis numbering
Ci   xlabel:label for abscissa
Ci   ... similar variables for y-axis
Ci   --- these variables are a work in progress ---
Ci   sfxdef:a string for mapping absicssa tic marks to a function
Ci   np    :number of points describing function
Ci   xp    :x-points for function
Ci   fxdef :y-points for function
Ci   --- 
Ci   title :
Ci   lt:   line thickness:
Ci         lt(1) is line thickness
Ci         lt(2) = 0 draws both top and bottom axis
Ci                 1 draws bottom only
Ci                 2 draws top only
Ci                 3 draws neither
Ci         lt(3) = 0 draws both left and right axes
Ci                 1 draws left only
Ci                 2 draws right only
Ci                 3 draws neither
Ci         NB: the frame is closed and filled only if lt(2)+lt(3)=0
Ci   cs
Co Outputs
Cr Remarks
Cr   Axes are drawn.
Cu Updates
Cu   19 Feb 03 Added numbering to rhs to top
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer mtx,mty,lt(3),modx,mody,np
      double precision tx0,tsx(1),rmtx,rmntx,ty0,tsy(1),rmty,rmnty,
     .cs(0:3),xp(np),fxdef(np)
      character*(*) xlabel,ylabel,title,sfxdef
      integer xnum,ynum

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      integer i,mt,it,i1mach,ipr,lmode,mt1,side,jx,lgunit,stdo
      double precision xtic,ytic,ticlen,xl,yl,xxtic,dxtic
      character*20 xnfmt,ynfmt
      character*120 s,s2
      double precision l,r,b,t
      logical maktic

C --- Setup ---
      call getpr(ipr)
      stdo = lgunit(1)
      l = plul
      r = plur
      b = plub
      t = plut
      if (logx) then
        l = dexp(plul)
        r = dexp(plur)
      endif
      if (logy) then
        b = dexp(plub)
        t = dexp(plut)
      endif
      if (ipr .ge. 10) then
        call awrit2(' FRME : axes at x=%1;4g %1;4g',s,len(s),0,l,r)
        if (logx) call awrit1('%a (log)',s,len(s),0,l)
        call awrit2('%a  y=%1;4g %1;4g',s,len(s),0,b,t)
        if (logy) call awrit1('%a (log)',s,len(s),0,l)
        call awrit1('%a  bold=%i',s,len(s),0,lt)
        if (lt(2)+lt(3) .gt. 0)
     .  call awrit4('%a: no axes for: '//
     .  '%?#n>0&n<>2# top,##%?#n>1# bottom,##'//
     .  '%?#n>0&n<>2# right,##%?#n>1# left,##%b',
     .  s,len(s),0,lt(2),lt(2),lt(3),lt(3))
        call awrit0('%a',s,len(s),-i1mach(2))
        s2 = '% ---'//s
        if (medium .eq. 3) call awrit0('%a',s2,len(s),-stream)
      endif
      if (medium .eq. 3) write (stream, 333)
  333 format('/max 0 def ')
      call plntys
      call plntyp(1,lt,0d0,0d0,0d0,0d0)

c --- draw frame ---
      call mve (l, b)
      if (lt(2) .eq. 0 .or. lt(2) .eq. 1) then
        call drw (r, b)
      else
        call mve (r, b)
      endif
      if (lt(3) .eq. 0 .or. lt(3) .eq. 2) then
        call drw (r, t)
      else
        call mve (r, t)
      endif
      if (lt(2) .eq. 0 .or. lt(2) .eq. 2) then
        call drw (l, t)
      else
        call mve (l, t)
      endif
      if (lt(2)+lt(3) .eq. 0) then
        call drwe(l, b, cs)
      elseif (lt(3) .eq. 0 .or. lt(3) .eq. 1) then
        call drw (l, b)
C       if (lt(2) .eq. 0 .or. lt(2) .eq. 1) call drw (r, b)
C      else
C        call mve(l, b)
C        call drwre(0d0, 0d0, cs)
      endif
*     write (stream,'(''% we are done'')')

c --- Draw and label tic marks, x axis ---
C ... tic marks for each side
      ytic = b
      ticlen = rmtx*(plut-plub)
      if (logy) ticlen = rmtx*dlog(t/b)
      s = ' '
      do  10  side = 1, 2
C ... Tic mark setup
        s2 = ' '
        call ticset(logx,modx,mtx,l,r,tsx,tx0,'x',lmode,xtic,mt,mt1,it,s2)
        if (side .eq. 1) s(6:len(s)) = s2
        maktic = (lt(2) .eq. 0 .or. lt(2) .eq. side) .and. lmode .ne. 10
        if (.not. maktic) goto 16
C ... for each tic do (mt1 is zero at major tics)
        do  15  i = 1, 999
C   ... for now
          xxtic = xtic
          if (sfxdef(1:1) .ne. ' ') then
            xxtic = xtic
            call polint(fxdef,xp,np,4,xxtic,0d0,0,jx,xtic,dxtic)
            if (ipr .ge. 50) then
              call awrit2(' plsym: abscissa %1:-1,4;4g mapped'//
     .        ' to %1:-1,4;4g',' ',80,stdo,xxtic,xtic)
            endif
          endif
          call mve(xtic, ytic)
          if (mt1 .eq. 0) then
            call drwre(0d0, ticlen, -2d0)
            if (ytic .eq. b .and. xnum .gt. 0) then
              if (medium .eq. 3) then
                call mve(xtic,ytic)
                call pslabl('a',xnfmt,xxtic,' ','l',' ',1,0d0)
              endif
            endif
            if (ytic .eq. t .and. xnum .lt. 0) then
              if (medium .eq. 3) then
                call mve(xtic,ytic)
                call pslabl('A',xnfmt,xxtic,' ','l',' ',1,0d0)
              endif
            endif
            mt1 = mt
          else
            call drwre(0d0, rmntx*ticlen, -2d0)
          endif
   14     continue
          mt1 = mt1-1
          call ticnxt(lmode,xtic,r,tsx,it)
          if (it .eq. 0) goto 16
   15   continue
        call rx('FRME:  too many tics in x')
   16   continue
        ytic = t
        ticlen = -ticlen
   10 continue

c --- Draw and label tic marks, y axis ---
C ... tic marks for each side
      xtic = l
      ticlen = rmty*(plur-plul)
      if (logx) ticlen = rmty*dlog(r/l)
      do  20  side = 1, 2
C ... Tic mark setup
        s2 = ' '
        if (side .eq. 2) s2 = s
        call ticset(logy,mody,mty,b,t,tsy,ty0,'y',lmode,ytic,mt,mt1,it,s2)
        maktic = (lt(3) .eq. 0 .or. lt(3) .eq. side) .and. lmode .ne. 10
        if (.not. maktic) goto 26
C ... for each tic do (mt1 is zero at major tics)
        do  25  i = 1, 999
          call mve(xtic, ytic)
          if (mt1 .eq. 0) then
            call drwre(ticlen, 0d0, -2d0)
            if (xtic .eq. l .and. ynum .gt. 0) then
              if (medium .eq. 3) then
                call mve(xtic,ytic)
                call pslabl('o',ynfmt,ytic,' ','l',' ',1,0d0)
              endif
            endif
            if (xtic .eq. r .and. ynum .lt. 0) then
              if (medium .eq. 3) then
                call mve(xtic,ytic)
                call pslabl('O',ynfmt,ytic,' ','l',' ',1,0d0)
              endif
            endif
            mt1 = mt
          else
            call drwre(rmnty*ticlen, 0d0, -2d0)
          endif
   24     continue
          mt1 = mt1-1
          call ticnxt(lmode,ytic,t,tsy,it)
          if (it .eq. 0) goto 26
   25   continue
        call rx('FRME:  too many tics in y')
   26   continue
        xtic = r
        ticlen = -ticlen
   20 continue

      if (ipr .ge. 30) call awrit0('%a',s2,-len(s),-i1mach(2))
      call plntyg
      if (medium .ne. 3) return
      s = s2
      s(1:5) = '% ---'
      call awrit0('%a',s,len(s),-stream)

C --- Label axes ---
      write (stream,*) '% FRME: Label axes ...'
C ... abcissa
      xl = 0.5d0*(plul + plur)
      yl = plub
      write (stream,130) (offxu+xl)*scalxu,(offyu+yl)*scalyu
      call pslabl('0',' ',0d0,xlabel,'h','f',1,0d0)
      if (xlabel .ne. ' ')
     .call pbbox((offxu+xl)*scalxu,(offyu+yl)*scalyu-3*fsiz)
C ... ordinate
      xl = plul
      yl = 0.5d0*(plut + plub)
      write (stream,230) (offxu+xl)*scalxu,(offyu+yl)*scalyu
      call pslabl('0',' ',0d0,ylabel,'v','r',1,0d0)
      if (ylabel .ne. ' ')
     .call pbbox((offxu+xl)*scalxu-3*fsiz,(offyu+yl)*scalyu)

C --- Add title ---
      write (stream,*) '% FRME: adding title ...'
      xl = plur
      yl = plut
      write (stream,330) (offxu+xl)*scalxu,(offyu+yl)*scalyu
      call pslabl('4',' ',0d0,title,'h','r',1,0d0)
      if (title  .ne. ' ')
     .call pbbox((offxu+xl)*scalxu,(offyu+yl)*scalyu+fsiz)

  130 format(2f8.1,' moveto 0 h 3 mul neg rmoveto')
  230 format(2f8.1,' moveto h 2 mul max add neg 0 rmoveto')
  330 format(2f8.1,' moveto 0 h 0.5 mul rmoveto')
      end
      subroutine plcrv(xp,yp,wp,wp2,wp3,np,xmin,xmax,ymin,ymax,opts,cs,
     .csw,iclos)
C- Draws or continues a smooth curve
C ----------------------------------------------------------------------
Ci Inputs:
Ci   xp,yp,np: set of x-y pairs, and number
Ci   wp,wp2 : color weights, used if csw is set
Ci   xmin..ymax: bounding box
Ci   opts :1s digit
Ci        :if nonzero, start a new curve whenever abscissa decreases
Ci        :10s digit
Ci        :if nonzero, plcrv clips curves to when they go outside
Ci         bounding box
Ci   iclos: 00: continues on existing curve, if one already open;
Ci              otherwise opens a new curve.
Ci          01: closes curve after drawing last point
Ci          10: closes any existing curve before starting new curve
Ci          11: both 10 and 01 above
Co Outputs:
Co   points smoothly connected on the plot stream
Cr Remarks
Cr   When a single color weight wp is present, color determined by
Cr       (1-wp)*cs + wp*csw
Cr   When two color weights are present, color determined by
Cr       (1-wp-wp2)*cs + wp*csw(1) + wp2*csw(2)
Cr   When three color weights are present, color determined by
Cr       (1-wp-wp2-wp3)*cs + wp*csw(1) + wp2*csw(2) + wp3*csw(3)
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer np,iclos,opts
      double precision xp(np),yp(np),wp(np),wp2(np),wp3(np),xmin,xmax,
     .ymin,ymax,xi,yi,cs(0:3),csw(0:3,3)
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz
C Local Variables
      logical lclip,lxbrk
      integer ip,ipr,getdig,i1,is,i1mach
      character s*120, s2*120, fmt*20
      double precision xsmin,ysmin,ysmax,xsmax,csl(0:3),xold,yold
      double precision box(4),cspt(0:3),w1,w2,w3
      save xold,yold

      lxbrk = mod(opts,10) .ne. 0
      lclip = mod(opts/10,10) .ne. 0

      if (lintyp .eq. 0) return
      call awrit2('%x(2f%i.%i,'' lto'')',fmt,20,0,ndwrt+8,ndwrt)
      call getpr(ipr)
C ... Local copy of cs which 10s digit of cs(0) is stripped
      call dcopy(4,cs,1,csl,1)
      if (cs(0) .ge. 10) csl(0) = mod(cs(0),10d0)
C ... Close curve if requested
      is = 0
      s = ' '
      if (getdig(iclos,1,10) .eq. 1) call plcrve(medium,pstat,-2d0,s,is)
      if (is .gt. 0) write(stream,*) s(1:is)
      is = 0
      s = ' '
      if (pstat .eq. 0) then
        call bin2a('plcrv: new',1,0,0,1,0,120,s,is)
        i1 = 2
        xold = xp(1)
        yold = yp(1)
      else
        call bin2a('plcrv: continue',1,0,0,1,0,120,s,is)
        i1 = 1
        call pposu(xold,yold)
      endif
      call awrit3('%a curve, %i pts lt %i bold %i',s,120,0,np,lintyp,
     .ltbld)
      if (cs(0) .gt. -2d0) call awrit1('%a col=%3:1d',s,120,0,cs(1))
      if (csw(0,1) .gt. 0d0)
     .call awrit1('%a colw=%3:1d',s,120,0,csw(1,1))
      if (csw(0,2) .gt. 0d0)
     .call awrit1('%a colw2=%3:1d',s,120,0,csw(1,2))
      if (csw(0,3) .gt. 0d0)
     .call awrit1('%a colw3=%3:1d',s,120,0,csw(1,3))
      if (lintyp .eq. 2 .and. linl3 .eq. 0)
     .call awrit2('%a len=(%d %d)',s,120,0,linl1,linl2)
      if (lintyp .eq. 2 .and. linl3 .ne. 0) call
     .awrit4('%a len=(%d %d %d %d)',s,120,0,linl1,linl2,linl3,linl4)
      if (ipr .ge. 40) call awrit1('%a',s,120,-i1mach(2),1)

C --------------- Postscript -----------------
      if (medium .ne. 3) goto 100
      s2 = '% ---' // s
      call awrit0('%a ---',s2,120,-stream)
C...  Clip beyond plot limits
      if (pstat .eq. 0) then
        xsmin = (offxu+xmin)*scalxu
        ysmin = (offyu+ymin)*scalyu
        ysmax = (offyu+ymax)*scalyu
        xsmax = (offxu+xmax)*scalxu
        if (logx) then
          xsmin = (offxu+dlog(xmin))*scalxu
          xsmax = (offxu+dlog(xmax))*scalxu
        endif
        if (logy) then
          ysmin = (offyu+dlog(ymin))*scalyu
          ysmax = (offyu+dlog(ymax))*scalyu
        endif
        call pbbox(xsmax,ysmax)
        call pbbox(xsmin,ysmin)
        call awrit8('gsave newpath %1;nd %1;nd moveto %1;nd %1;nd '//
     .  'lineto',' ',120,stream,ndwrt,xsmin,ndwrt,ysmin,ndwrt,xsmin,
     .  ndwrt,ysmax)
        call awrit8('              %1;nd %1;nd lineto %1;nd %1;nd '//
     .  'lineto closepath clip',' ',120,stream,ndwrt,xsmax,
     .  ndwrt,ysmax,ndwrt,xsmax,ndwrt,ysmin)
C        call awrit8('gsave newpath %1;1d %1;1d moveto %1;1d %1;1d '//
C     .   'lineto %1;1d %1;1d lineto %1;1d %1;1d lineto closepath clip',
C     .   ' ',120,stream,xsmin,ysmin,xsmin,ysmax,xsmax,ysmax,xsmax,ysmin)

C        if (cs(0) .ge. 10) then
C          call dcopy(4,cs,1,cspt,1)
C          if (csw(0,1) .ge. 10) then
C            call dscal(3,wp,cspt(1),1)
C            call daxpy(3,1-wp,csw(1,1),1,cspt(1),1)
C          endif
C          call awrit1('%3:1;4d setrgbcolor',' ',80,stream,cspt(1))
C        endif
        if (cs(0) .ge. 10) call awrit1('%3:1;4d setrgbcolor',' ',80,
     .  stream,cs(1))
        call mve(xp(1),yp(1))
      endif

      if (lclip) then
        box(1) = xmin
        box(2) = ymin
        box(3) = xmax
        box(4) = ymax
      endif

C --- Draw the curve ---
      if (lxbrk) then
        if (csw(1,1) .ne. 0) call rx('not ready for lxbrk and colw')
        do  1  ip = i1, np
          if (ip .gt. i1) then
            if (xp(ip) .gt. xp(ip-1)) then
              call drw(xp(ip),yp(ip))
            else
              if (ipr .ge. 4) call awrit1(
     .        '    ... break at point %i',' ',80,i1mach(2),ip)
              call mve(xp(ip),yp(ip))
            endif
          else
            call drw(xp(ip),yp(ip))
c         write(stream,fmt)  (offxu+xp(ip))*scalxu,(offyu+yp(ip))*scalyu
          endif
    1   continue
      else
        cspt(0) = 0
        if (csw(0,1) .ne. 0) then
          w1 = (wp(i1)+wp(max(i1-1,1)))/2
          w2 = 0
          w3 = 0
          if (csw(0,2) .ne. 0) then
            w2 = (wp2(i1)+wp2(max(i1-1,1)))/2
          endif
          if (csw(0,3) .ne. 0) then
            w3 = (wp3(i1)+wp3(max(i1-1,1)))/2
          endif
          call dpcopy(cs,cspt,1,4,1-w1-w2)
          call dpadd(cspt,csw(0,1),1,4,w1)
          call dpadd(cspt,csw(0,2),1,4,w2)
          call dpadd(cspt,csw(0,3),1,4,w3)
          cspt(0) = 15
        endif
        if (i1.le.np) call cdrw(xp(i1),yp(i1),xold,yold,lclip,box,cspt)
C        if (csw(0,1) .ne. 0) then
C          call mve(xp(i1),yp(i1))
C        endif
        if (logx .or. logy .or. lrot .or. lclip .or. csw(0,1).ne.0) then
          do  ip = i1+1, np-1
            if (csw(0,1) .ne. 0) then
              w1 = (wp(ip)+wp(max(ip-1,1)))/2
              w2 = 0
              if (csw(0,2) .ne. 0) then
                w2 = (wp2(ip)+wp2(max(ip-1,1)))/2
              endif
              if (csw(0,3) .ne. 0) then
                w3 = (wp3(i1)+wp3(max(i1-1,1)))/2
              endif
              call dpcopy(cs,cspt,1,4,1-w1-w2)
              call dpadd(cspt,csw(0,1),1,4,w1)
              call dpadd(cspt,csw(0,2),1,4,w2)
              call dpadd(cspt,csw(0,3),1,4,w3)
              cspt(0) = 15
            endif
            call cdrw(xp(ip),yp(ip),xold,yold,lclip,box,cspt)
C            if (csw(0,1) .ne. 0) then
C              call mve(xp(ip),yp(ip))
C            endif
          enddo
        else
C     ... Above is more generic, but this way is faster
          do  3  ip = i1+1, np-1
    3     write(stream,fmt)  (offxu+xp(ip))*scalxu,(offyu+yp(ip))*scalyu
        endif
C       Last point
        if (i1+1 .le. np) then
          if (csw(0,1) .ne. 0) then
            w1 = (wp(ip)+wp(max(ip-1,1)))/2
            w2 = 0
            if (csw(0,2) .ne. 0) then
              w2 = (wp2(ip)+wp2(max(ip-1,1)))/2
            endif
            if (csw(0,3) .ne. 0) then
              w3 = (wp3(i1)+wp3(max(i1-1,1)))/2
            endif
            call dpcopy(cs,cspt,1,4,1-w1-w2)
            call dpadd(cspt,csw(0,1),1,4,w1)
            call dpadd(cspt,csw(0,2),1,4,w2)
            call dpadd(cspt,csw(0,3),1,4,w3)
            cspt(0) = 15
          endif
          call cdrw(xp(np),yp(np),xold,yold,lclip,box,cspt)
        endif
      endif
      goto 20

C --- Generic plcrv ---
  100 continue
C In this case, end does not matter
      xi = xp(1)
      yi = yp(1)
      if (ipr .ge. 40) call awrit2('    ... start line %1,6;6d %1,6;6d',
     .' ',80,i1mach(2),xi,yi)
      call mve(xi,yi)
C      if (xi .ge. xmin .and. xi .le. xmax .and.
C     .    yi .ge. ymin .and. yi .le. ymax) then
C        berr = .false.
C      else
C        berr = .true.
C      endif
      do  10  ip = 2, np
C cludge for now...
        xi = xp(ip)
        yi = yp(ip)
        if (xi .ge. xmin .and. xi .le. xmax .and.
     .  yi .ge. ymin .and. yi .le. ymax) then

          if (lxbrk) then
            if (xp(ip) .gt. xp(ip-1)) then
              call drw(xp(ip),yp(ip))
            else
              if (ipr .ge. 40) call awrit1(
     .        '    ... break at point %i',' ',80,i1mach(2),ip)
              call mve(xp(ip),yp(ip))
            endif
          else
            call drw(xp(ip),yp(ip))
          endif
          if (ipr .ge. 50) print 345, 'draw within boundary',xi,yi
  345     format(' plcrv: ',a20,2f12.6)
        else
          call mve(xi,yi)
          if (ipr .ge. 50) print 345, 'mve outside boundary',xi,yi
        endif
   10 continue

C --- Close the curve ---
   20 continue
      if (mod(iclos,10) .eq. 1) then
        is = 0
        s = ' '
        call plcrve(medium,pstat,csl,s,is)
        if (medium .eq. 3) call bin2a('grestore',1,0,0,1,0,120,s,is)
        if (is .gt. 0) write(stream,'(/a)') s(1:is)
      endif

      end
      subroutine plsym(xp,yp,zp,np,stype,syma,symcs)
C- Draw symbol at a collection of (x,y) pairs
C ----------------------------------------------------------------------
Ci Inputs
Ci   xp,yp (zp) x-y pairs where to draw symbols.
Ci   stype: compound integer of 1+10's digit and 100's digit:
Ci          1-10's digit for symbol type: 1-9 for:
Ci                 x, square, diamond, +, polygon, circle, arrow,
Ci                 errbar, timelin, hist
Ci          100's  digit = j>0: use zp in place of syma(j)
Ci   symcs:    color for symbol
Ci          (0) is used by plcrve and specifies fill and color
Ci              0 => no fill
Ci              1 => fill no color
CC         whether it is in color or B&W
CC __cs(1..3) specify the RGB color
CC __cs(4) specifies the symbol bold (used in cymcs only for now)
CC         cs(4)<0 => take value from -lt specification
Ci   syma(1),  size of symbol:
Ci     Other parameters for syma are specific to symbol.
Ci     polygon: (1): size
Ci              (2): number of sides
Ci              (3): angle
Ci     diamond: (2): dimension of height of diamond
Ci     arrow:   tip at (xp,yp), and syma:
Ci              (1,2): tail coordinates relative to head (GU)
Ci              (3): head length (fraction of arrow length)
Ci              (4): head angle rel to tail (degrees)
Ci              (5): head length along axis (fraction of arrow length)
Ci     timelin  zp  bar length (UU)
Ci              (1) endbar height (UU)
Ci              (2) bar thickness
Ci              (3) endbar thickness
Cr Remarks
Cu Updates
Cu   13 May 04  Extra option in errbar
Cu   05 Jun 02  Added histogram symbol
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer stype,np
      double precision syma(8),symcs(0:4),xp(np),yp(np),zp(np)

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local parameters
      double precision x,y,xold,yold
      integer nsym
      parameter (nsym=10)
      integer nb,ib,ipr,ip,is,nsyma(nsym),typex,msyma,ttype
      double precision tpi,phi0,radx,rady,d,dx,dy,hsyma
      character*80 s
      character*7 synam(nsym),nam*(*)
      data synam /'x','square','diamond','+','polygon','circle',
     .'arrow','errbar','timelin','hist'/
      data nsyma /1,1,1,1,2,1,5,2,3,2/
      data tpi /6.283185307179586d0/

      goto 1

      entry plsyn(nam,ttype)
C- return index to name of symbol
      call tokmat(nam,synam,nsym,len(nam),' ',ttype,is,.false.)
      ttype = ttype+1
      return

C --- Setup ---
    1 continue
      typex = mod(stype,100)
      if (typex .eq. 0) return
      msyma = stype/100
      if (msyma .gt. 0) hsyma = syma(msyma)
      call getpr(ipr)
      is = 0
      s = ' '
      call plcrve(medium,pstat,-2d0,s,is)
      if (is .gt. 0) write(stream,*) s(1:is)
      call awrit5('plsym: '//synam(typex)//' %i pts  gs=%1;5d  attr=%'
     .//char(ichar('0')+nsyma(typex))//':1;5d  x,y(1)=%1;4g %1;4g',
     .s,80,0,np,symcs(1),syma,xp,yp)
      call skpblb(s,79,is)
      if (ipr .ge. 40) print *, s(1:is+1)
      if (medium .eq. 3)
     .write (stream,'(/a)') '% --- '//s(1:is+1)//' ---'

      call plntys
      if (symcs(4) .ge. 0) ltbld = nint(symcs(4))
      call plntyp(1,ltbld,0d0,0d0,0d0,0d0)

C --- For each point, do ----
      do  2  ip = 1, np
        if (msyma .gt. 0) syma(msyma) = zp(ip)
        call mve(xp(ip),yp(ip))
        x = .015d0*syma(1)*scalxg/scalxu
        if (syma(2) .ne. 0 .and. iabs(typex) .le. 4) then
          y = .015d0*syma(2)*scalyg/scalyu
        else
          y = .015d0*syma(1)*scalyg/scalyu
        endif
        goto (10,20,30,40,50,60,70,80,90,100), iabs(typex)
        goto 5

C --- Symbol type cross ---
   10   continue
        call drwr(x,y)
        call drwr(-2*x,-2*y)
        call drwr(x,y)
        call drwr(-x,y)
        call drwr(2*x,-2*y)
        call drwr(-x,y)
        goto 5

C --- Symbol type square (rectangle) ---
   20   continue
        call mver(x,y)
        call drwr(-2*x,0d0)
        call drwr(0d0,-2*y)
        call drwr(+2*x,0d0)
        call drwr(0d0,+2*y)
        goto 5

C --- Symbol type diamond ---
   30   continue
        call mver(0d0,y)
        call drwr(-x,-y)
        call drwr(+x,-y)
        call drwr(+x,+y)
        call drwr(-x,+y)
        goto 5

C --- Symbol type +  ---
   40   call drwr(0d0,y)
        call drwr(0d0,-2*y)
        call drwr(0d0,y)
        call drwr(x,0d0)
        call drwr(-2*x,0d0)
        call drwr(x,0d0)
        goto 5

C --- Symbol type poly: syma(1..3): radius, no. sides, angle 1st pt(rad)
   50   continue
        radx = dsqrt(2d0)*x
        rady = dsqrt(2d0)*y
        nb = nint(syma(2))
        phi0 = tpi/(2*nb) + syma(3)
        x = radx*dsin(phi0)
        y =-rady*dcos(phi0)
        call mver(x,y)
        do  55  ib = 1, nb
          xold = x
          yold = y
          x = radx*dsin(ib*tpi/nb+phi0)
          y =-rady*dcos(ib*tpi/nb+phi0)
          call drwr(x-xold,y-yold)
   55   continue
        goto 5

C --- Symbol type cirle, radius a. PostScript only ---
   60   continue
        if (medium .eq. 1) then
          syma(2) = 8
          goto 50
        endif
        if (medium .ne. 3) return
        write (stream,530) xcm,ycm,syma(1)*6
        pstat = 1
        goto 5

C --- Symbol type arrow ---
   70   continue
C ...   Draw tail
        call mver(syma(1)*(plur-plul),syma(2)*(plut-plub))
        call drwe(xp(ip),yp(ip),-2d0)
        if (syma(3) .eq. 0) goto 5
C ...   Position of three arrow points, UU. phi0 is arrow direction.
        dx = syma(1)*(plur-plul)*scalxu
        dy = syma(2)*(plut-plub)*scalyu
        phi0 = datan2(dy,dx)
        d = dsqrt(dx**2+dy**2)
        radx = phi0 + tpi*syma(4)/360
        call drw(d*syma(3)*dcos(radx)/scalxu + xp(ip),
     .  d*syma(3)*dsin(radx)/scalyu + yp(ip))
        call drw(d*syma(5)*dcos(phi0)/scalxu + xp(ip),
     .  d*syma(5)*dsin(phi0)/scalyu + yp(ip))
        radx = phi0 - tpi*syma(4)/360
        call drw(d*syma(3)*dcos(radx)/scalxu + xp(ip),
     .  d*syma(3)*dsin(radx)/scalyu + yp(ip))
        call drw(xp(ip),yp(ip))
        goto 5

C --- Symbol type errbar ---
   80   continue
        d = yp(ip)+zp(ip)*(1+syma(2))
        if (logy) then
          if (d .le. 0) d = yp(ip)
        endif
        call drw(xp(ip),d)
        call drwr(x,0d0)
        call drwr(-2*x,0d0)
        call drwr(x,0d0)
        d = yp(ip)-zp(ip)*(1-syma(2))
        if (logy) then
          d = yp(ip)**2/(yp(ip)+(zp(ip))*(1-syma(2)))
          if (d .le. 0) d = yp(ip)
        endif
        call drw(xp(ip),d)
        call drwr(x,0d0)
        call drwr(-2*x,0d0)
        call drwr(x,0d0)
        goto 5

C --- Symbol type timelin ---
   90   continue
        call plntyp(1,nint(syma(3)),0d0,0d0,0d0,0d0)
        call drwr(0d0,syma(1)/2)
        call drwr(0d0,-syma(1))
        call drwre(0d0,syma(1)/2,symcs)
        call plntyp(1,nint(syma(2)),0d0,0d0,0d0,0d0)
        call drwre(zp(ip),0d0,symcs)
        call plntyp(1,nint(syma(3)),0d0,0d0,0d0,0d0)
        call drwr(0d0,syma(1)/2)
        call drwr(0d0,-syma(1))
        call drwre(0d0,syma(1)/2,symcs)
        goto 5

C --- Symbol type histogram ---
  100   continue
        call mve(xp(ip),syma(2))
C       call drwr(-syma(1)/2*(plur-plul),0)
        call drwr(-syma(1)/2,0)
        call drwr(0d0,yp(ip)-syma(2))
C       call drwr(syma(1)*(plur-plul),0)
        call drwr(syma(1),0)
        call drwr(0d0,syma(2)-yp(ip))
C       call drwr(-syma(1)/2*(plur-plul),0d0)
        call drwr(-syma(1)/2,0d0)
        goto 5

C --- Cleanup for this point ---
    5   continue
        is = 0
        s = ' '
C        if (medium .eq. 3) then
C          is = awrite(' closepath gsave%3:1;4d setrgbcolor '//
C     .      'fill grestore',s,80,0,symcs(1),0,0,0,0,0,0,0)
C        endif
        call plcrve(medium,pstat,symcs,s,is)
        if (is .gt. 0) write(stream,*) s(1:is)
    2 continue
      call plntyg
      if (msyma .gt. 0) syma(msyma) = hsyma

  530 format(' newpath',3f8.2,' 0 360 arc')
      end
      subroutine plcrve(medium,pstat,cs,s,is)
C- Returns string signaling termination of a curve; sets pstat to 0
C ----------------------------------------------------------------------
Ci Inputs
Ci   medium
Ci   pstat  If pstat is not zero, nothing is done
Ci   cs     color scale.  0 flags filling (see Outputs) 1..3 are colors
Co Outputs
Co   pstat is set to 0
Co   Appended to string s (Postscript mode):
Co   cs(0)<1:  s = 'stroke'
Co   cs(0)=1:  s = 'closepath stroke'
Co   cs(0)=2:  s = 'closepath gsave # setgray fill grestore stroke'
Co   cs(0)=3:  s = 'closepath gsave ### setrgbcolor fill grestore stroke'
Co   cs(0)=4:  s = '# setgray stroke'
Co   cs(0)=5:  s = '### setrgbcolor stroke'
Co   cs(0)=6:  s = '# setgray closepath stroke'
Co   cs(0)=7:  s = '### setrgbcolor closepath stroke'
Cu Updates
Cu   13 Feb 01 revised what plcrve does depending on cs(0)
C ----------------------------------------------------------------------
      double precision cs(0:3)
      integer medium,pstat,is,ls
      character s*(*)
      double precision cs0

      if (pstat .ne. 0 .and. medium .eq. 3) then

C        if (cs(0) .ge. 10) then
C          print *, 'plcrve',cs(0)
C        endif
        cs0 = mod(cs(0),10d0)
        if (cs0 .ge. 7) then
          ls = len(s)
          call awrit1(' %3:1;4d setrgbcolor closepath',
     .    s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 .ge. 6) then
          ls = len(s)
          call awrit1(' %1;4d setgray closepath',
     .    s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 .ge. 5) then
          ls = len(s)
          call awrit1(' %3:1;4d setrgbcolor ',
     .    s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 .ge. 4) then
          ls = len(s)
          call awrit1(' %1;4d setgray ',
     .    s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 .ge. 3) then
          ls = len(s)
          call awrit1(' closepath gsave%3:1;4d setrgbcolor '//
     .    'fill grestore',s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 .ge. 2) then
          ls = len(s)
          call awrit1(' closepath gsave %1;4d setgray fill grestore',
     .    s(is+1:ls),ls-is,0,cs(1))
          call skpblb(s,ls,is)
          is = is+1
        elseif (cs0 .ge. 1) then
          call bin2a('closepath',1,0,0,1,0,len(s),s,is)
        endif
        call bin2a('stroke',1,0,0,1,0,len(s),s,is)
      endif
      pstat = 0
      end
      subroutine key(style,xk,yk,xlen,xmin,xmax,ymin,ymax,
     .symtyp,syma,symcs,crvcs,legend,blank)
C- Draw the key's legend at (xk,yk)
C ----------------------------------------------------------------------
Ci Inputs: xk,yk (user units), xlen: length of line to draw, syma,b,typ:
Ci         style: sets style for legend (See Remarks)
Ci                0 let this routine choose legend
Ci                1 or 2: use this style for legend
Ci         symbol, legend: legend
Cr Remarks
Cr   Draws legend containing symbol s. drawn in one of two styles:
Cr     style=1:  ---- s ----
Cr     style=2:  s
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      double precision xk,yk,xlen,syma,symcs(0:3),crvcs(0:3,2)
      integer style,symtyp,blank
      character*60 legend

C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C Local Parameters
      double precision x(2),y(2),xmax,xmin,ymax,ymin
      double precision xm,ym,csw(0:3,3)
      integer i1mach,isw,lstyle
      character*80 s

      lstyle = style
      if (style .eq. 0) lstyle = 1
      if (style .eq. 0 .and. lintyp .eq. 0) lstyle = 2
      call dpzero(csw,12)

      s = ' '
      call awrit5(' key:  x,y=%1;4g,%1;4g(UU)  lt %i  sym %i  len=%1;4g'
     .,s,80,0,xk,yk,lintyp,symtyp,xlen)
      call awrit1('%a  legend='//legend,s,80,0,xk)
      call awrit1('%a',s,80,-i1mach(2),xk)
      if (medium .eq. 3) call awrit1(' %% KEY%a',s,80,-stream,xk)
      x(1) = xk
      x(2) = xk + xlen
      if (logx) x(2) = xk*dexp(xlen)
      y(1) = yk
      y(2) = yk
      if (lstyle .eq. 2) x(2) = x(1)
      call plcrv(x,y,x,x,x,2,xmin,xmax,ymin,ymax,0,crvcs,csw,11)
      if (.not. logx) x(1) = (x(1)+x(2))/2
      if (      logx) x(1) = dsqrt(x(1)*x(2))
      call plsym(x,y,0d0,1,symtyp,syma,symcs)
C      if (lintyp .eq. 0) then
C        x(1) = x(1) + xlen/10
C        x(2) = x(2) - xlen/10
C        call plsym(x,y,0d0,2,symtyp,syma,symcs)
C        x(1) = x(1) - xlen/10
C        x(2) = x(2) + xlen/10
C      else
C        if (.not. logx) x(1) = (x(1)+x(2))/2
C        if (      logx) x(1) = dsqrt(x(1)*x(2))
C        call plsym(x,y,0d0,1,symtyp,syma,symcs)
C      endif
      if (medium .eq. 3) then
        call mve(x(2),y(2))
        call pltam(xm,ym,-1,1d0)
        call awrit3(' %1;1d %1;1d moveto h %?#n#.8#.5# mul h .5 mul neg'
     .  //' rmoveto',' ',120,stream,xm,ym,isw(x(2).eq.x(1)))
C       write (stream,030) xm,ym
        call pslabl('1',' ',0d0,legend,'h','f',blank,0d0)
      endif
  030 format(2f8.1,' moveto h .8 mul h .5 mul neg rmoveto')
      end

      subroutine pslabl(shift,nfmt,val,strn,hvl,fr,blank,rots)
C- Put a PostScript label or string at the current postion
C ----------------------------------------------------------------------
Ci Inputs:
Ci    shift: one of the following
Ci           'a','o' for abscissa, ordinate labelling, 
Ci                   In this case, val is converted into a string
Ci           '0'-'4' for centring or placing in quadrants 1-4 about
Ci                   the current position
Ci           'tx'    for simple text 'show'
Ci           'cc'    or any of the nine pairs of characters 
Ci                   lu lc ld  cu cc cd  ru rc rd
Ci                   1st character justfies left, center, or right
Ci                   2nd character justfies top, middle, or bottom
Ci    strn:  the string to print, according to hvl. A '\' in the string
Ci           forces a line break. Multiple lines are written forwards
Ci           starting at the current point or reversed ending at the
Ci           current point, according to fr
Ci    fr:    'f' or 'r' for whether multiple lines are written forward
Cr           or reversed
Ci    nfmt:  awrite format for val (see hvl='l')
Ci    hvl:   if 'l', then convert val to a string, using awrite format.
Ci           if 'h' or 'v', label is horizontal or vertical
Ci    rots:  rotation of string.
Cr Remarks
Cr  pslabl writes a  PostScript segment assuming currentpoint has been
Cr  previously set. It does NOT take the curent postion from COMMON
Cb Bugs
Cb   This is a hoplessly confusing routine.  Plans are to move some
Cb   of its functionality to pstr.
Cu Updates
Cu   30 Apr 07 bug fix: rots passed to fwrite 
Cu   25 Apr 07 hvl='l' translates fortran aEn format to a.10^n
C ----------------------------------------------------------------------
      implicit none
C Passed Paramters
      character*1 shift*(*),hvl,fr,strn(0:*),nfmt*20
      integer blank
      double precision val,rots
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C local parameters
      integer maxlen,maxlns,isw,i1,i2
      parameter (maxlns=10,maxlen=120)
      character*(maxlen) label,show*6,newl*9,s,ss,line(maxlns),strn2,
     .bs*1
      integer lbrk(maxlns),nlines,ich,count,il,ip,ilen,l1,l2,l3,
     .icop,llen(maxlns),iprint,i1mach
      logical cmdopt

c      call spchar(1,bs)
      bs='\'
      show = 'show'//shift
      if (shift(1:2) .eq. 'tx') show = 'show'
      newl = 'newline'//hvl//fr
      if (hvl .eq. 'l') then
        ip = 0
        count = 0
        label = ' '
C       if (ndec .gt. 4) stop 'PSLABL: more than four decimal places !?'
C       call bin2a('(f10.4)',0,ndec,val,4,count,maxlen,label,ip)
C       call bin2a('d4:0',0,ndec,val,4,count,maxlen,label,ip)
        call bin2a0(0)
        call awrit1(nfmt,label,maxlen,0,val)

C  ...  Convert fortran 'e' format to 10^... and use fwrite
        call chrps2(label,'eEdD',4,len(label),ip,i2)
        if (i2 .ne. 0) then
          call word(label,1,i1,i2)
          ss = label
          ss = label(1:ip) // '^{.}10^{' // label(ip+2:i2) // '}'
          call word(ss,1,i1,i2)
C         Shorten 1.10^... to 10^...
          if (ss(i1:i1+4) .eq. '1^{.}') i1=i1+5
          if (ss(i1:i1+5) .eq. '-1^{.}') then
            i1=i1+5
            ss(i1:i1) = '-'
          endif
          write (stream,200)
  200     format('gsave /font /Symbol def fontsize font choosefont')
          if (shift .eq. 'a'.or. shift .eq. 'A') then
            write (stream,030) xcm,ycm
            call fwrite(stream,blank,'h',0d0,'showcD','newlinevf',
     .      ss(i1:i2),i2-i1+1,label)
          else if (shift .eq. 'o' .or. shift .eq. 'O') then
            write (stream,030) xcm,ycm,' h hshift mul neg 0 rmoveto'
            call fwrite(stream,blank,'h',0d0,'showlc','newlinevf',
     .      ss(i1:i2),i2-i1+1,label)
          endif
          write (stream,300)
  300     format ('grestore')
C  ...  Just use standard awrite
        else
          call skpblb(label,maxlen,ilen)
          ilen = ilen + 1
          if (iprint() .gt. 50) call awrit2(
     .    ' pslabl: label "'//label(1:ilen)//'" at %1;2d,%1;2d (MU)',
     .    s,maxlen,i1mach(2),xcm,ycm)
          call pbbox(xcm,ycm)
          write (ss,031) xcm,ycm,label(1:ilen),show
          call awrit3('%xgsave blank-%a%?#n#on#off#%?#n# %d rotate#%j#',
     .    s,maxlen,0,isw(blank.ne.0),isw(rots.ne.0),rots)
C       s = 'gsave blank-on'
C       if (blank .eq. 0) s = 'gsave blank-off'
          call awrit1('%a '//ss,s,maxlen,0,1)
          call awrit1('%a',s,maxlen,-stream,1)
        endif
      else
        if (hvl .eq. 'v') then
          if (shift .eq. '1') show = 'show2'
          if (shift .eq. '2') show = 'show3'
          if (shift .eq. '3') show = 'show4'
          if (shift .eq. '4') show = 'show1'
        endif
        call skpblb(strn,maxlen,ilen)
        if (ilen .lt. 0) return
        nlines = 1
        ich = 0
        do  1  il = 1, maxlns
C         call chrpos(strn,'\\',ilen,ich)
          call chrpos(strn,bs,ilen,ich)
          lbrk(il) = ich
          if (ich .eq. ilen) goto 2
          ich = ich + 1
          nlines = nlines + 1
    1   continue
    2   continue
        ip = 0
        do  3  il = 1, nlines
          if (nlines .eq. 1) then
            icop = ilen+1
          else
            if (il .eq. 1) then
              icop = lbrk(1)
            elseif (il .eq. nlines) then
              icop = lbrk(nlines) - lbrk(nlines-1)
            else
              icop = lbrk(il)-1 - lbrk(il-1)
            endif
          endif
C         call strcop(line(il),strn(ip),icop,'\\',ich)
          call strcop(line(il),strn(ip),icop,bs,ich)
          llen(il) = ich
          if (iprint() .gt. 50) write (*,100) line(il)(1:llen(il))
          ip = ip + icop+1
    3   continue
        if (fr .eq. 'f') then
          l1 = 1
          l2 = nlines
          l3 = 1
        else
          l1 = nlines
          l2 = 1
          l3 = -1
        endif
        do  4  il = l1, l2, l3
          s = 'gsave blank-'
          call awrit3('%a%?#n#on#off#%?#n# %d rotate#%j#',s,maxlen,0,
     .    blank.ne.0,isw(rots.ne.0),rots)
          if (cmdopt('-plaintext ',11,0,strn2) .or.
     .    shift(1:2) .eq. 'tx') then
            if (hvl .eq. 'h') then
C             write (stream,130) line(il)(1:llen(il)),show,newl
              write (ss,131) line(il)(1:llen(il)),show,newl
            else
C             write (stream,230) line(il)(1:llen(il)),show,newl
              write (ss,231) line(il)(1:llen(il)),show,newl
            endif
            call awrit1('%a '//ss,s,maxlen,0,1)
            call awrit1('%a',s,maxlen,-stream,1)
          else
            call fwrite(stream,blank,hvl,rots,show,newl,
     .      line(il)(1:llen(il)),llen(il),strn2)
          endif
    4   continue
      endif
  030 format(2f8.1,' moveto':/a)
  031 format(2f8.1,' moveto (',a,') ',a6,' grestore')
  130 format('gsave (',a,') ',a6,' grestore ',a9)
  131 format('(',a,') ',a6,' grestore ',a9)
  230 format('gsave 90 rotate (',a,') ',a6,' grestore ',a9)
  231 format('90 rotate (',a,') ',a6,' grestore ',a9)
  100 format(' PSLABL new line: ',a)
      end

      subroutine fwrite(stream,blank,hvl,rots,show,newl,line,ilen,ss)
C- Generalised write using FORTRAN 
C ----------------------------------------------------------------------
Ci Inputs:
Ci   stream,blank
Ci   hvl:   if 'v', rotate string 90 degrees; otherwise:
Ci   rots:  if nonzero, rotate string by angle
Ci   show:  same as shift in pslabl EXCEPT these options are not allowed
Ci          'a','o' 
Ci          'tx'    for simple 'show'
Ci   newl:  controls justification in newlines and should be one of:
Ci          newlinehf, newlinevf, newlinehr, newlinevr
Ci   line:string to be written; ilen its length (including control chars)
Ci   ss: work string
Co Outputs:
Co   
Cr Remarks
Cr   fplot will accept a string with control characters ^ _ ~ to 
Cr   make super- and sub-scripts and greek symbols. The substring to
Cr   be converted is enclosed in { .. }. For example
Cr   ~{m}_{0}=4~{p}.10^{7} is equivalent to
Cr   $\mu$_0=4$\pi$.10$^7$ in TeX. The { .. } may contain any number
Cr   of characters and may also be nested, say, to make sub-sub-scripts.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ilen,maxlen,stream,blank
      parameter (maxlen=120)
      character*(maxlen) show*6,newl*9,line,hvl*1,ss
      double precision rots
C ... Local parameters
      integer i,tlen,iw
      character c*1,cp*1,cn*1,cw*(maxlen)
      character bkslsh*1

c      call spchar(1,bkslsh)
      bkslsh="\"
C --- strip line of its control chars to determine true stringwidth ---
      tlen = 0
      do  i = 1, ilen
        c = line(i:i)
        if (c .ne. '^' .and. c .ne. '_' .and. c .ne. '~' .and.
     .  c .ne. '{' .and. c .ne. '}') then
          tlen = tlen + 1
          ss(tlen:tlen) = c
        endif
      enddo

C --- save the graphics state ---
      write (stream,100)

      if (hvl .eq. 'v') then
        write (stream,101) 90d0
      elseif (rots .ne. 0) then
        write (stream,101) rots
      endif

C --- inform PostScript of the true stringwidth ---
      write (stream,102) ss(1:tlen)

C --- align text according to right/left up/down rules ---
      if (show .eq. 'showru' .or. show .eq. 'show1') then
      elseif (show .eq. 'showcu') then
        write (stream,104)  
      elseif (show .eq. 'showlu' .or. show .eq. 'show4') then
        write (stream,105) 
      elseif (show .eq. 'showrc') then
        write (stream,106) 
      elseif (show .eq. 'showcc' .or. show .eq. 'show0') then
        write (stream,107) 
      elseif (show .eq. 'showlc') then
        write (stream,108) 
      elseif (show .eq. 'showrd' .or. show .eq. 'show2') then
        write (stream,109) 
      elseif (show .eq. 'showld' .or. show .eq. 'show3') then
        write (stream,110) 
      elseif (show .eq. 'showcd') then
        write (stream,111) 
      elseif (show .eq. 'showcD') then
        write (stream,112) 
      endif

C --- blank out area to be written on ---
      if (blank .eq. 1) then
        write (stream,103)
      endif

C --- generalised write ---
      cw = ' '
      iw = 0
      write (stream,121)
      do i = 1, ilen
        c = line(i:i)
        cp = ' '
        if (i .gt. 1) cp = line(i-1:i-1)
        cn = ' '
        if (i .lt. ilen) cn = line(i+1:i+1)
        if ((c .eq. '^' .and. cn .eq. '{') .or.
     .  (c .eq. '_' .and. cn .eq. '{') .or.
     .  (c .eq. '~' .and. cn .eq. '{')) then
C         Flush internal buffer
          if (iw .gt. 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
        elseif (c .eq. '{' .and. cp .eq. '^') then
          write (stream,122) 
        elseif (c .eq. '{' .and. cp .eq. '_') then
          write (stream,123) 
        elseif (c .eq. '{' .and. cp .eq. '~') then
          write (stream,124) 
        elseif (c .eq. '}') then
C         Flush internal buffer
          if (iw .gt. 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
          if (cn .ne. '}') then
            write (stream,125)
          endif
        elseif (c .eq. '(') then
C         Flush internal buffer
          if (iw .gt. 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
          write (stream,129) bkslsh,'050'
        elseif (c .eq. ')') then
C         Flush internal buffer
          if (iw .gt. 0) then
            write (stream,130) cw(1:iw)
            iw = 0
          endif
          write (stream,129) bkslsh,'051'
        else
          iw = iw + 1
          cw(iw:iw) = c
        endif
      enddo
C     Flush internal buffer
      if (iw .gt. 0) then
        write (stream,130) cw(1:iw)
        iw = 0
      endif

      write (stream,200) newl

  100 format ('gsave')
  101 format (f6.1,' rotate')
  102 format ('/w (',a,') stringwidth pop def')
  103 format ('savecurrentpoint gsave w /x exch def'/
     .'newpath xbak ybak h tails mul sub moveto'/
     .'x 0 rlineto 0 h 1.1 mul h tails mul add rlineto x',
     .' neg 0 rlineto'/'closepath 1 setgray fill grestore')
  104 format ('w 2 div neg 0 rmoveto')
  105 format ('w neg 0 rmoveto')
  106 format ('0 h 2 div neg rmoveto')
  107 format ('w 2 div neg h 2 div neg rmoveto')
  108 format ('w neg h 2 div neg rmoveto')
  109 format ('0 h neg rmoveto')
  110 format ('w neg h neg rmoveto')
  111 format ('w 2 div neg h neg rmoveto')
  112 format ('w 2 div neg h 1.7 mul neg rmoveto')
  121 format ('/fntsiznow fontsize def /accshift 0 def /frac 1 def')
  122 format ('/fntsiznow fntsiznow smaller mul def'/
     .'fntsiznow font choosefont'/
     .'/accshift accshift h shift mul frac mul add def'/
     .'0 h shift mul frac mul rmoveto'/
     .'/frac frac 0.7 mul def')
  123 format ('/fntsiznow fntsiznow smaller mul def'/
     .'fntsiznow font choosefont'/
     .'/accshift accshift h shift mul frac mul sub def'/
     .'0 h shift mul frac mul neg rmoveto'/
     .'/frac frac 0.7 mul def')
  124 format ('fontsize /Symbol choosefont')
  125 format ('/fntsiznow fontsize def fontsize font choosefont'/
     .'0 accshift neg rmoveto /accshift 0 def /frac 1 def')
  129 format ('(',a1,a,') show')
  130 format ('(',a,') show')

  200 format ('grestore ',a9)
      end

      subroutine pstr(x,y,rots,string,arg1,arg2,arg3,pos,units,optio)
C- Paint a string in the current font at coordinates x,y
C ----------------------------------------------------------------------
Ci Inputs:
Ci   optio:  a compound of digits specifying options
Ci           1s digit
Ci           1  blank-on
Ci          10s digit
Ci           1  string is a format entering into awrit2
Ci           2  string is a format ...?
Ci         100s digit
Ci           1  string a current coordinates; thus x,y are not used
Ci   x,y:    coordinates at which to paint string
Ci   rots:   rotation angle at which to write
Ci   units:  'u' for x,y in user units, 'm' for medium (cm) units
Ci   string: string to paint, or formatting command; see optio.
Ci   arg1,arg2: optional arguments entering into string format
Ci   pos:    compound of two characters, the first one of
Ci           l,c,r followed by one of u,c,d
Cb   Bugs
Cb     10's and 100's digit of optio never tested.
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      double precision x,y,rots,arg1,arg2,arg3
      character pos*2,units*1,string*(*)
      integer optio
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz
      integer iprint,i1mach,isw,blank,maxlen,ifmt
      parameter (maxlen=120)
      character*(maxlen) s, strn
      double precision xx,yy

      if (medium .ne. 3) return

C --- Postscript label ---
      blank = mod(optio,10)
      ifmt  = mod(optio/10,10)
      if (ifmt .eq. 0) then
        strn = string
      elseif (ifmt .eq. 1) then
        strn = ' '
        call awrit3(string,strn,maxlen,0,arg1,arg2,arg3)
      else
        call rxi('pstr: bad formatting switch',ifmt)
      endif
      call awrit3(' pstr: '//pos//': '//'%?#n#rot=%d #%j#'//
     .'blk=%l  "'//strn,s,80,0,isw(rots.ne.0),rots,blank.ne.0)
      if (mod(optio/100,10) .eq. 0) then
        call awrit0('%a" at',s,80,0)
        xx = x
        yy = y
        if (units .eq. 'm' .and. scalxu*scalyu .ne. 0) then
          xx = x/scalxu - offxu
          yy = y/scalyu - offyu
          if (logx) xx = exp(xx)
          if (logy) yy = exp(yy)
          call awrit2('%a %1;4g,%1;4g(u)',s,80,0,xx,yy)
          xx = x
          yy = y
        elseif (units .ne. 'm') then
          xx = (offxu+x)*scalxu
          yy = (offyu+y)*scalyu
          if (logx) xx = (offxu+dlog(x))*scalxu
          if (logy) yy = (offyu+dlog(y))*scalyu
          call rxx(scalxu*scalyu .eq. 0,'pstr: user''s units not set')
          call awrit2('%a %1;4g,%1;4g(u)',s,80,0,x,y)

        endif
        call awrit2('%a %1,1;1d,%1,1;1d(m)',s,80,0,xx,yy)
      endif
      if (iprint() .ge. 30) call awrit0('%a',s,len(s),-i1mach(2))
      if (medium .eq. 3) call awrit0('%5o%0p%% ---%a',s,len(s),-stream)
      if (mod(optio/100,10) .eq. 0) then
        write (stream,130) xx,yy
        call pbbox(xx,yy)
  130   format(2f8.1,' moveto')
      endif
      call pslabl(pos,' ',0d0,strn,'h','f',blank,rots)
      end
      subroutine pbbox(x,y)
C- Update bbox
      implicit none
C Passed Parameters
      double precision x,y
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

C      if (x .lt. bbox(1) .or. x .gt. bbox(2) .or.
C     .    y .lt. bbox(3) .or. y .gt. bbox(4)) then
C        print *, 'x,y=',x,y
C      endif
      bbox(1) = min(bbox(1),x)
      bbox(2) = max(bbox(2),x)
      bbox(3) = min(bbox(3),y)
      bbox(4) = max(bbox(4),y)


      end

      subroutine setfnt(font,fontsz)
C- Set font in PostScript protocol
C ----------------------------------------------------------------------
Ci Inputs:
Ci   font and font_size
C ----------------------------------------------------------------------
      implicit none

      character*12 font
      integer fontsz,i
C Plot parameters block
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

      fsiz = fontsz
      if (medium .eq. 3) write (stream,130) font,fontsz
      call lodsyv('ch',1,.6d0*fontsz,i)
  130 format(' % FONT: set font and fontsize ..'/
     .'/font ',a12,' def /fontsize ',i3,' def'/
     .'/h fontsize 0.6 mul def % h is about char height'/
     .'fontsize font choosefont')
      end
      subroutine psinit(stream,irot,transx,transy,scale)
C- Initialise PostScript file with defs and defaults
C ----------------------------------------------------------------------
Ci Inputs
Ci   stream:file logical unit
Ci   transx:shift of origin from bottom left corner of paper
Ci   transy:shift of origin from bottom left corner of paper
Ci   scale :shrinks or expands entire plot
Co Outputs
Co   definitions are written to file stream
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer stream,irot
      double precision transx,transy,scale

      write(stream,030)
C --- translation, scaling, rotating ---
      call awrit4('/inch{72 mul}def /cm{72 mul 2.54 div}def' //
     .' %1;1d %1;1d translate  %1;1d %1;1d scale',
     .' ',90,stream,transx, transy, scale, scale)
      if (irot .ne. 0) write(stream,031)
C --- fontsize and definitions ---
      write (stream,230)
      write (stream,231)
      write (stream,240)
      write (stream,241)
      write (stream,250)

  030 format('%!PS-Adobe-2.0 generated by plsub.f ...')
  031 format(' 90 rotate 0 -612 translate')
  130 format(2f8.1,' translate ',2f8.1,' scale')
  230 format(
     .'/tails 0.4 def % set to 0 for blanking chars without tails'/
     .'/lto{lineto}def % shorthand for lineto'/
     .'/blank-on{/b 1 def}def /blank-off{/b 0 def}def blank-on'/
     .'/smaller 0.67 def % sub/superscripts are this much smaller'/
     .'/shift   0.33 def % and are shifted by this much of the char ',
     .'height'/
     .'/hshift  0.10 def % horizontal space, fraction of char ',
     .'height'/'% push size then font onto stack to select font'/
     .'/choosefont {findfont exch scalefont setfont} def'/
     .'% paint object on top of stack at current position')
  231 format(
     .'/debug {dup /str 30 string def str cvs show} def'/
     .'% finding width of longest ordinate label (current width on ',
     .'stack) ...'/'/max 0 def /mxw{/w exch def w max gt {/max w ',
     .'def} if}def'/'% string centring, justification in the four ',
     .'quadrants'/'% and absicca and ordinate label centring ...'/
     .'/showru{show}def'/'/showcu{dup stringwidth pop 2 div neg ',
     .'0 rmoveto show}def'/'/showlu{dup stringwidth pop neg 0 ',
     .'rmoveto show}def'/'/showrc{0 h 2 div neg rmoveto show}def'/
     .'/showcc{dup stringwidth pop 2 div neg h 2 div neg rmoveto ',
     .'show}def'/'/showlc{dup stringwidth pop neg h 2 div neg ',
     .'rmoveto show}def'/'/showrd{0 h neg rmoveto show}def')
  240 format(
     .'/showld{dup stringwidth pop neg h neg rmoveto show}def'/
     .'/showcd{dup stringwidth pop 2 div neg h neg rmoveto show}def'/
     .'/showa{gsave fontsize /Symbol choosefont'/
     .'       dup stringwidth pop 2 div neg h 1.5 mul neg rmoveto ',
     .'show grestore}def'/
     .'/showA{gsave fontsize /Symbol choosefont'/
     .'       dup stringwidth pop 2 div neg h 0.5 mul rmoveto ',
     .'show grestore}def'/
     .'/showo{gsave fontsize /Symbol choosefont'/
     .'       dup stringwidth pop dup mxw h 2 div add neg h 2 div ',
     .'neg rmoveto show grestore}def'/
     .'/showO{gsave fontsize /Symbol choosefont'/
     .'       dup h 2 div h 2 div neg rmoveto show grestore}def'/
     .'% aliases for pslabl ...'/
     .'/show0{showcc}def /show1{show}def /show2{showrd}def ',
     .'/show3{showld}def'/'/show4{showlu}def')
  241 format('% newline macros for ',
     .'horizontal and vertical labelling ...'/
     .'/newlinehf{currentpoint fontsize sub moveto}def'/
     .'/newlinevf{currentpoint exch fontsize add exch moveto}def'/
     .'/newlinehr{currentpoint fontsize add moveto}def'/
     .'/newlinevr{currentpoint exch fontsize sub exch moveto}def'/
     .'/savecurrentpoint {currentpoint /ybak exch def /xbak exch ',
     .'def} def')
  250 format('% end preamble ...'/)

      end


      block data plt
      implicit none
C --- Common variables for plot parameters ---
C plml,r,b,t: mediums's units for left,right,bottom and top coordinates
C             demarcating edge of plot medium.
C plul,r,b,t: user's units for left,right,bottom and top coordinates
C             demarcating edge of plot medium.
C plgl,r,b,t: user's units marking frame of interior of a current plot.
C xcm,ycm:    current pen position in medium's units.
C xcu,ycu:    current pen position in user's units.
C pcol:       current pen color
C stream:     stream where to send plot output (hardware dependent)
      double precision plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox(4),scalxc,scalyc,offxc,offyc,
     .xcm,ycm,strtch,xcu,ycu,xcc,ycc,xb,yb,
     .linl1,linl2,linl3,linl4,lnprd,lnused
      integer stream,pcol,lintyp,ltbld,medium,pstat,ndpi,ndwrt,fsiz
      logical lrot,noclip,logx,logy
      common /pltdat/ plml,plmr,plmb,plmt,plpl,plpr,plpb,plpt,
     .plul,plur,plub,plut,scalxu,scalyu,offxu,offyu,
     .plgl,plgr,plgb,plgt,scalxg,scalyg,offxg,offyg,
     .cosa,sina,bbox,scalxc,scalyc,offxc,offyc,
     .xcm,ycm,xcu,ycu,xcc,ycc,strtch,pcol,stream,xb,yb,
     .lintyp,ltbld,linl1,linl2,linl3,linl4,lnprd,lnused,
     .medium,lrot,noclip,logx,logy,pstat,ndpi,ndwrt,fsiz

      data cosa /0d0/, sina /1d0/
      data plml /0d0/, plmr /1d0/, plmb /0d0/, plmt /1d0/
      data plgl /0d0/, plgr /1d0/, plgb /0d0/, plgt /1d0/
      data plpl /0d0/, plpr /1d0/, plpb /0d0/, plpt /1d0/
      data plul /0d0/, plur /1d0/, plub /0d0/, plut /1d0/
      data xcm  /0d0/, ycm  /0d0/, strtch /1d0/ ndpi /600/ ndwrt /1/
      data lintyp /1/, linl3 /0d0/, linl4 /0d0/ ltbld /2/
      data stream /11/, pstat /0/, noclip /.true./
      data logx /.false./, logy /.false./
      end
