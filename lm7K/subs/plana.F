      subroutine plana(wksize,npadl,npadr,nbaspp,slabl,slat,sspec,ssite,
     .  ves,pnu,qnu)
C- Plane analysis for lmto programs
C ----------------------------------------------------------------------
Ci Inputs
Ci   wksize:size of global work array w
Ci   nbaspp:size of doubly padded basis (layer programs)
Ci          nbaspp = nbas + 2*(nbas(left bulk) + nbas(right bulk))
Ci   slabl :list of species labels.
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat plat2 avw vol
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: z qc lmxa rmt
Ci     Stored:    *
Ci     Passed to: *
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: class spec
Ci     Stored:    clabel pos
Ci     Passed to: spacks spackv
Ci   ves   :electrostatic potential at MT boundaries
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Co Outputs
Cu Updates
Cu   03 Nov 04 update last arg in partok call (need integer)
Cu   31 Oct 03 plana orders planes L-C-R for layer geometry
Cu   01 Mar 02 plana can write to site file
Cu             Altered argument list
Cu   12 Feb 02 Adapted for ASA v6
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer wksize,npadl,npadr,nbaspp
      character*8 slabl(*)
      double precision slat(1),sspec(1),ssite(1),ves(1),pnu(1),qnu(1)
c     double precision pnu(nl,nsp,1),qnu(3,nl,nsp,1),pp(1),bas(3,1)
C Heap:
      integer w(1)
      common /w/ w

C Local parameters
      integer ipc(nbaspp),ips(nbaspp),lmx(nbaspp)
      double precision qc(nbaspp),qt(nbaspp),dq(nbaspp)
      double precision z(nbaspp),bas(3,nbaspp),wsr(nbaspp)
      integer oplan,oiplan,nplane,i1mach,mcont2,recln0,offset,iosw,i,
     .  fext,fopn,iarg(10),iout,j,lgunit,iosits,
     .  oqzpl,oqplan,odqpln,ovharp,ovmadp,opotqz,quit,
     .  k,m,n,ifi,npl,mapc(200),owk,ib,ic,is,iclbas,obas2,onatpl,l,lio
      integer nbas,nclspp,nl,nsp,nglob,idxdn(1),mxint,nbasp,nclasp,nspec
      double precision darg(10),cellen,ddot,avw,carea,vol,hhat(3),hvc
      double precision alat,plat(3,3),plat0(3,3),hdotb,planvc(3,3)
      double precision cphi,sphi,ctheta,stheta,phi,
     .  tmp(3,3),tmp2(3,3),tmp3(3,3),pi,vol0,glat(3,3),r(5),r2(5),
     .  platl(3,3),platr(3,3)
      character*8 nam
      parameter (mcont2=17+2-1, recln0=72)
      integer ncont2(0:mcont2)
      logical sw,a2bin,switch(1)
      character*8 clabl(nbaspp)
      character*80 carg(10)

      data ncont2 /-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/,
     .  iosw/1/

C     Make ipc, ips, clabl
      do  ib = 1, nbaspp
        call upack('site class spec',ssite,ib,ic,is,0,0)
        ipc(ib) = ic
        ips(ib) = is
        call spacks(0,'site clabel',ssite,clabl(ic),ib,ib)
        call upack('spec z qc lmxa rmt',sspec,is,z(ic),qc(ic),lmx(ic),
     .    wsr(ic))
      enddo
      nbas = nglob('nbas')
      nspec = nglob('nspec')
      nbasp  = (nbas+nbaspp)/2
      nclasp = mxint(nbasp,ipc)
      nclspp = mxint(nbaspp,ipc)
C     Unpack basis
      call spackv(10,'site pos',ssite,1,nbaspp,bas)

      nl = nglob('nl')
      nsp = nglob('nsp')
      call upack('lat alat plat plat2 avw vol',slat,alat,plat,planvc,
     .  avw,vol)
      if (ddot(3,planvc(1,3),1,planvc(1,3),1) .eq. 0) then
        call dcopy(9,plat,1,planvc,1)
      endif
      hvc = dsqrt(ddot(3,planvc(1,3),1,planvc(1,3),1))
      call getq(nsp,nl,lmx,nclspp,z,pnu,qnu,0,sspec,qc,qt,dq)

      call dcopy(9,plat,1,plat0,1)
      rewind (fopn('PLAN'))
      call defdr(oplan,nbasp)
      call defi(oiplan,nbasp)
      call defi(onatpl,nbasp)
      call gtplan(npadl,npadr,nbasp,bas,planvc(1,3),ipc,clabl,
     .  nplane,w(oplan),w(oiplan),w(onatpl))

C --- Input from from stdin ---
C Offset in work array so that read ctrl file read at end ...
      offset = wksize*i1mach(6) - 1*recln0 - 100
  100 continue
      print 332
  332 format(/'Menu for plane options.  Type help for help')
      do  110  i = 0, mcont2
        call dpzero(darg,10)
        carg(1) = ' '
        carg(2) = ' '
        call lmplio(recln0,i1mach(1),offset,iosw,ncont2(i),
     .    plat,planvc,iout,iarg,darg,carg)
        if (iout .ne. 0) goto 120
  110 continue
      iosw = 1
      goto 100

  120 continue
      assign 100 to quit
      goto (100,100,3,4,5,6,7,8,9,10,11,12,13,14,15,16), iout

C --- Write file ---
    3 continue
C     print *, 'implementing write'
      do  31  ic = 1, nclspp
        mapc(ic) = ic
   31 continue
      call dscal(nclspp,1/avw,wsr,1)
      call prtbas(nbasp,ipc,clabl,1,nclasp,mapc,-1,1,bas,z,avw,wsr,
     .  ves,switch,nl,pnu,qnu)
      call dscal(nclspp,avw,wsr,1)
      goto quit

C --- Read plane ---
    4 continue
      stop 'need pnu qnu z wsr idxdn for ioplan'
      ifi = fopn('PLAN')
      call ioplan(nbas,nl,nsp,nclspp,plat,
     .  clabl,z,wsr,
     .  pnu,qnu,
     .  idxdn,
     .  ipc,bas,w(oiplan),npl,ifi)

C --- Begin loop that makes various adjustments ---
      do  41  i = 1, 3
        hhat(i) = plat(i,3) / hvc
   41 continue
      do  42  i = 1, npl
        ib = i+nbasp
        ic = i+nclspp
        ipc(ib) = ic

C --- Check for replication of class names ---
c        do  43  j = 1, nclspp
        do  43  j = 1, ic-1
          if (clabl(j) .ne. clabl(ic)) goto 43
          print 432, clabl(ic), j
  432     format('   Class ',a4,' already named in class',i3)
          print 336
#if SVS
  336     format('   Confirm name, or type new name: '\)
#elif APOLLO
  336     format('   Confirm name, or type new name: '$)
#else
  336     format('   Confirm name, or type new name: ')
#endif
          read(*,'(a4)') nam
          if (nam .ne. ' ')  clabl(ic) = nam
   43   continue

C --- Adjust bas to proper projection onto h ---
        hdotb = ddot(3,bas(1,ib),1,hhat,1)
        bas(1,ib) = bas(1,ib) + hhat(1)*(darg(1) - hdotb)
        bas(2,ib) = bas(2,ib) + hhat(2)*(darg(1) - hdotb)
        bas(3,ib) = bas(3,ib) + hhat(3)*(darg(1) - hdotb)
        print 412, clabl(ic), (bas(j,ib), j=1,3)
  412   format(' Atom ',a4,' inserted at:',3f12.5)
   42 continue

      nbas = nbas+npl
      nbasp = nbasp+npl
      nclasp = nclasp+npl
      nclspp = nclspp+npl

C --- Remap class table, finding equivalent names ---
      call defi(owk,nbasp)
      call fixcl(nbasp,nclspp,ipc,clabl,w(owk))
      call rlse(owk)

C --- Make new plane table ---
      call rlse(oplan)
      call defdr(oplan,nbasp)
      call defi(oiplan,nbasp)
      call gtplan(npadl,npadr,nbasp,bas,planvc(1,3),ipc,clabl,
     .  nplane,w(oplan),w(oiplan),w(onatpl))
      goto quit

C --- Write a single plane to plan file  ---
    5 continue
      stop 'need pnu qnu z wsr idxdn for ioplan'
      call ioplan(nbas,nl,nsp,nclspp,plat,
     .   clabl,z,wsr,
     .   pnu,qnu,
     .   idxdn,
     .   ipc,bas,w(oiplan),iarg,-fopn('PLAN'))
      goto quit

C --- Skip over planes  ---
    6 continue
      stop 'need pnu qnu z wsr for ioplan'
      ifi = fopn('PLAN')
      rewind ifi
      do  61  i = 1, iarg(1)-1
        print *, '... skipping past plane ', i
        call ioplan(nbas,nl,nsp,nclspp,plat,
     .    clabl,z,wsr,
     .    pnu,qnu,
     .    idxdn,
     .    ipc,bas,w(oiplan),npl,ifi)
   61 continue
      goto quit

C --- Change extension to name  ---
    7 continue
      call fclose(fopn('PLAN'))
      ifi = fext(carg)
      ifi = fopn('PLAN')
      goto quit

C --- Stretch plat ---
    8 continue
      do  81  i = 1, 3
        hhat(i) = plat(i,3) / hvc
   81 continue
      call daxpy(3,darg(1),hhat,1,plat(1,3),1)
      print 351
      print 350, ((plat0(m,k),m=1,3),(plat(m,k),m=1,3),k=1,3)
  351 format(/13x,'ORG PLAT',26x,'NEW PLAT')
  350 format(3f10.5,5x,3f10.5)
      goto quit

C --- Sum charge in plane ---
    9 continue
      call defdr(oqplan,-nplane)
      call defdr(odqpln,-nplane)
      call defdr(oqzpl,-nplane)
      call defdr(ovharp,-nplane)
      call defdr(ovmadp,-nplane)
      call defdr(opotqz,nplane)
      cellen = alat*hvc
      carea =  avw
      carea =  vol / cellen
      call xxpl(nbasp,alat,ipc,nplane,w(oplan),w(oiplan),ves,iarg,
     .  qt,dq,cellen,carea,w(oqplan),w(odqpln),w(oqzpl),w(ovharp),
     .  w(onatpl),wsr,w(ovmadp),darg,w(opotqz))
      call rlse(oqplan)
      goto quit

C --- Slide planes ---
   10 continue
      do  101  i = 1, 3
        hhat(i) = plat(i,3) / hvc
  101 continue
      do  102  ib = 1, nbasp
        j = w(oiplan+ib-1)
        if (j .ge. iarg(1) .and. j .le. iarg(2)) then
          call daxpy(3,darg(1),hhat,1,bas(1,ib),1)
          print 337, ib, (bas(j,ib), j = 1,3)
  337     format(' Sliding site',i3,' to',3f10.5)
        endif
  102 continue

C --- Make new plane table ---
      call rlse(oplan)
      call defdr(oplan,nbasp)
      call defi(oiplan,nbasp)
      call gtplan(npadl,npadr,nbasp,bas,planvc(1,3),ipc,clabl,
     .  nplane,w(oplan),w(oiplan),w(onatpl))
      goto quit

C --- Rename classes ---
   11 continue
      call defi(owk,nbaspp)
      call fixcl(nbaspp,nclspp,ipc,clabl,w(owk))
      call rlse(owk)
      do  111 ic = 1, nclspp
        print 112, clabl(ic), ic, nclspp
  112   format(' Set name of class ', a4, ', number',i3,' of',i3,' to:')
        read(*,'(a4)') nam
        if (nam .ne. ' ')  clabl(ic) = nam
  111 continue
      goto quit

C --- Shift ves in planes ---
   12 continue
      do  122  ic = 1, nclspp
        ib = iclbas(ic,ipc)
        j = w(oiplan+ib-1)
        if (j .ge. iarg(1) .and. j .le. iarg(2)) then
          print 338, clabl(ic), ves(ic), ves(ic)+darg(1)
  338     format(' Shift ves of class ',a4,' from',f10.5,' to',f10.5)
          ves(ic) = ves(ic)+darg(1)
        endif
  122 continue
      goto quit

C --- Rotate plat ---
   13 continue
      pi = 4*datan(1d0)
      darg(1) = darg(1)*pi
      do  131  i = 1, 3
        hhat(i) = planvc(i,3)/hvc
  131 continue
c      hhat(1) = .6123724
c      hhat(2) = .6123724
c      hhat(3) = .5
      if (dabs(hhat(2))+dabs(hhat(2)) .lt. 1d-8) then
        phi = 0
      else
        phi = datan2(hhat(2),hhat(1))
      endif
      cphi = dcos(phi)
      sphi = dsin(phi)
      ctheta = hhat(3)
      stheta = dsqrt(1-ctheta**2)
      tmp(1,1) = cphi*ctheta
      tmp(2,1) = sphi*ctheta
      tmp(3,1) = -stheta
      tmp(1,2) = -sphi
      tmp(2,2) = cphi
      tmp(3,2) = 0
      tmp(1,3) = cphi*stheta
      tmp(2,3) = sphi*stheta
      tmp(3,3) = ctheta
      print 333, tmp
      ctheta = dcos(darg(1))
      stheta = dsin(darg(1))
      call dpzero(tmp2,9)
      tmp2(1,1) = ctheta
      tmp2(2,1) = stheta
      tmp2(1,2) = -stheta
      tmp2(2,2) = ctheta
      tmp2(3,3) = 1
  333 format(3f10.5)
      call dmpy(tmp,3,1,tmp2,3,1,tmp3,3,1,3,3,3)
      call dmpy(tmp3,3,1,tmp,1,3,tmp2,3,1,3,3,3)
      print 349
  349 format(/9x,'ROTATION MATRIX')
      print 333, tmp2
      call dmpy(tmp2,3,1,plat,3,1,tmp,3,1,3,3,3)
      print 351
      print 350, ((plat(m,k),m=1,3),(tmp(m,k),m=1,3),k=1,3)
      call dcopy(9,tmp,1,plat,1)
      call defdr(obas2,3*nbasp)
      call dmpy(tmp2,3,1,bas,3,1,w(obas2),3,1,3,nbasp,3)
      call dcopy(3*nbasp,w(obas2),1,bas,1)
      call rlse(obas2)
      goto quit

C --- Display projection of basis vectors onto planvc ---
   14 continue
      call dinv33(plat,0,tmp,vol0)
      print *, ' ... Projections in fractions of plat'
      print 341
  341 format(' IB',14x,'BAS',33x,'PROJ')
  342 format(i3,3f10.5,5x,3f10.5)
      do  145  ib = 1, nbasp
        call dmpy(tmp,3,1,bas(1,ib),3,1,tmp2,3,1,3,1,3)
        print 342, ib, (bas(j,ib), j=1,3), (tmp2(j,1), j=1,3)
  145 continue

      print *, ' ... Projection when subtract out planvc'
      darg(1) = 0
C --- Adjust bas to proper projection onto h ---
      do  141  i = 1, 3
        hhat(i) = planvc(i,3)/hvc
  141 continue
      print 341
      do  146  ib = 1, nbasp
        hdotb = ddot(3,bas(1,ib),1,hhat,1)
        tmp(1,1) = bas(1,ib) + hhat(1)*(darg(1) - hdotb)
        tmp(2,1) = bas(2,ib) + hhat(2)*(darg(1) - hdotb)
        tmp(3,1) = bas(3,ib) + hhat(3)*(darg(1) - hdotb)
        print 342, ib, (bas(j,ib), j=1,3), (tmp(j,1), j=1,3)
  146 continue
      goto quit

C --- Atomic positions within range of plane ---
C To circumscribe 3 planes, first 4 classes only, use eg
C pos 6 6 2 const p1<=4*sqrt(2)&p2<=2.01&p3>1&p3<1.5&ic<=4
C and extract with grep pos log.ext | pextract ssssssssss hijb
C Variables are p1,p2,p3,x,y,z,ic,n, ib(first occurrence)
   15 continue
C ... Shift all basis vectors to first cell
      call dinv33(plat,0,glat,vol0)
      print *, ' ... shifting basis vectors to first cell'
      do  151  ib = 1, nbasp
          call dmpy(glat,3,1,bas(1,ib),3,1,r,3,1,3,1,3)
          call dpzero(r2,3)
          do  152  i = 1, 3
          r(i) = -nint(r(i)-.49999d0)
  152     call daxpy(3,r(i),plat(1,i),1,r2,1)
          call daxpy(3,1d0,r2,1,bas(1,ib),1)
          call awrit4('ib=%i %6psh=%3:1,5;5d (%3;0d)%42pnew bas='//
     .      '%3:1,6;6d',carg,80,i1mach(2),ib,r2,r,bas(1,ib))
  151 continue
C ... Generate table of all positions
      call defrr(obas2, 5*(iarg(1)+1)*(iarg(2)+1)*(iarg(3)+1)*nbasp)
      call defrr(owk,     (iarg(1)+1)*(iarg(2)+1)*(iarg(3)+1)*nbasp)
      call dcopy(3,planvc(1,3),1,hhat,1)
      call dscal(3,1/hvc,hhat,1)
      l = 0
      do  154  ib = 1, nbasp
        do  155  i = 0, iarg(1)
        do  155  j = 0, iarg(2)
        do  155  k = 0, iarg(3)
          call dpzero(r(3),3)
          call daxpy(3,dble(i),plat(1,1),1,r(3),1)
          call daxpy(3,dble(j),plat(1,2),1,r(3),1)
          call daxpy(3,dble(k),plat(1,3),1,r(3),1)
          call daxpy(3,1d0,bas(1,ib),1,r(3),1)
          r(1) = ddot(3,hhat,1,r(3),1)
          r(2) = ib
          call dpscop(r,w(obas2),5,1,5*l+1,1d0)
          if (carg(1) .eq. ' ') then
            l = l+1
          else
            call lodsyv('ib',1,dble(ib),n)
            call lodsyv('ic',1,dble(ipc(ib)),n)
            call lodsyv('x',1,r(3),n)
            call lodsyv('y',1,r(4),n)
            call lodsyv('z',1,r(5),n)
            call lodsyv('n',1,r(1),n)
            do  156  n = 1, 3
  156       r2(n) = ddot(3,planvc(1,n),1,r(3),1)/
     .          dsqrt(ddot(3,planvc(1,n),1,planvc(1,n),1))
            call lodsyv('p1',1,r2(1),n)
            call lodsyv('p2',1,r2(2),n)
            call lodsyv('p3',1,r2(3),n)
            n = 0
            if (.not. a2bin(carg,sw,0,0,' ',n,-1)) goto 155
            if (sw) l = l+1
          endif
  155   continue
  154 continue
      call dvshel(5,l,w(obas2),w(owk),0)
C ... Write out table of positions within (lo,hi)
      print 356
  356 format(' ic       ib',13x,'projection',26x,'position')
      write(lgunit(2),357)
  357 format(5x,'ic',7x,'ib',13x,'projection',22x,'position')
      m = 0
      do  158  i = 0, l-1
        call dpscop(w(obas2),r,5,5*i+1,1,1d0)
        j = nint(r(2))
        k = ipc(j)
        if (darg(4).eq.-99 .and. darg(5).eq.-99 .or.
     .     (r(1) .ge. darg(4) .and. r(1) .le. darg(5))) then
          m = m+1
          do  159  n = 1, 3
  159     r2(n) = ddot(3,planvc(1,n),1,r(3),1)/
     .            dsqrt(ddot(3,planvc(1,n),1,planvc(1,n),1))
          print 363, k,clabl(k),j,(r2(n),n=1,3),r(3),r(4),r(5)
  363     format(i3,1x,a4,i4,2(f12.6,2f11.6))
          write(lgunit(2),364) k,clabl(k),j,(r2(n),n=1,3),(r(n),n=3,5)
  364     format('pos',i4,1x,a4,i4,3f11.6,3f10.5)
        endif
  158 continue
      k = nbasp*(iarg(1)+1)*(iarg(2)+1)*(iarg(3)+1)
      call awrit3(' pos: %i points written; %i satisfied constraints'//
     .  ' out of %i',carg,80,i1mach(2),m,l,k)
      call rlse(obas2)
      goto quit

C --- Write to site file ---
   16 continue
C     k will be the number of sites to write to site file (normally nbas)
      k = nbas
C     tmp will be the plat written to site file (normally plat)
      call dcopy(9,plat,1,tmp,1)

C     Some padding sites .. possibly alter k, tmp
      if (nbasp .gt. nbas) then
        call upack2('lat platl platr',slat,platl,platr)
        if (carg(1) .eq. '-ppad') then
          k = nbaspp
          carg(1) = carg(2)
        elseif (carg(1) .eq. '-pad') then
          k = nbasp
          call daxpy(3,-1d0,platl(1,3),1,tmp(1,3),1)
          call daxpy(3,-1d0,platr(1,3),1,tmp(1,3),1)
          carg(1) = carg(2)
        else
          call daxpy(3,-2d0,platl(1,3),1,tmp(1,3),1)
          call daxpy(3,-2d0,platr(1,3),1,tmp(1,3),1)
        endif
      endif
      if (carg(1) .eq. ' ') carg(1) = 'site'

      lio = 1000*(2+4+8+16+32) + 1
      j = iosits(lio,3d0,0,carg,ifi,slabl,alat,tmp,k,nspec,sspec,ssite)
      call fclose(ifi)
      goto quit

      end
      subroutine ioplan(nbas,nl,nsp,nclspp,plat,
     .  clabl,z,wsr,
     .  pnu,qnu,
     .  idxdn,
     .  ipc,bas,iplane,ipl,ifi)
C- I/o of site and class data for a plane of atoms
Ci   ipl write plane i to file (write only)
Co   ipl number of atoms in plane read (read only)

C     implicit none

      integer ifi

C Parameters for plane
      integer iplane(2),ipl

C Parameters and for options, structure, basis
      integer nbas,nl,nsp,nclspp
      double precision plat(9)

C Class parameters
      double precision z(1),wsr(1)
      character*8 clabl(1)
      integer idxdn(nl,1)
      double precision pnu(nl,nsp,1),qnu(3,nl,nsp,1)

C Site parameters
      double precision bas(3,nbas)
      integer ipc(nbas)

C local variables
      integer np, jfi, i, ic, j, k, m, ib
      double precision plat0(9)

      if (ifi .lt. 0) goto 1
C --- Read branch ---
      ic = nclspp
      ib = nbas
      jfi = ifi
      read(jfi,*,end=120,err=120) np, plat0
      do  101  i = 1, 6
        if (dabs(plat(i)-plat0(i)) .lt. 1d-4) goto 101
        print *, 'ioplan:  plat mismatch in plane file'
  101 continue
      do  110  i = 1, np
        ic = ic+1
        ib = ib+1
        read(jfi,335) m, clabl(ic)
        read(jfi,*) z(ic), wsr(ic), (bas(m,ib), m=1,3)
        read(jfi,*) (idxdn(j,ic), j=1,nl)
        read(jfi,*) ((pnu(j,k,ic), j=1,nl), k=1, nsp)
        read(jfi,*) (((qnu(m,j,k,ic), m=1,3), j=1,nl), k=1, nsp)
        print *, 'read site', i, ', class ', clabl(ic)
  110 continue
      ipl = np
      return
  120 ipl=0
      return

C --- Write branch ---
    1 continue
      jfi = -ifi
C Count number of atoms to write
      np = 0
      do  5  i = 1, nbas
        if (iplane(i) .eq. ipl) np = np+1
    5 continue

      write(jfi,*) np, plat
      do  10  i = 1, nbas
        if (iplane(i) .eq. ipl) then
          ic = ipc(i)
          print *, 'writing site', i, ', class ', clabl(ic)
          write(jfi,335) ic, clabl(ic)
  335     format(i2,1x,a4)
          write(jfi,*) z(ic), wsr(ic), (bas(m,i), m=1,3)
          write(jfi,*) (idxdn(j,ic), j=1,nl)
          write(jfi,*) ((pnu(j,k,ic), j=1,nl), k=1, nsp)
          write(jfi,*) (((qnu(m,j,k,ic), m=1,3), j=1,nl), k=1, nsp)
        endif
   10 continue

      end
      subroutine xxpl(nbas,alat,ipc,nplane,plane,iplane,ves,iarg,
     .  qt,dq,cellen,carea,qplan,dqplan,qzpl,vharpl,natpln,rmax,vmadpl,
     .  darg,potqz)
C     implicit none
      integer nbas,ipc(1),nplane,iplane(1),natpln(1),iarg(2)
      double precision plane(1),qt(1),dq(1),qplan(1),qzpl(1),ves(1),
     .  cellen,carea,alat,dqplan(1),vharpl(1),rmax(1),vmadpl(1),
     .  darg(2),potqz(1)

C local variables
      integer i,j,iprint,lgunit,imnarg,imxarg,imnpln
      double precision sumq,field0,pi8,dipole,vintra,
     .  rytoev,qdipol,dipmad,diphar,wfbpln,wfbmad,wfbhar,
     .  wfsmad,wfshar,efpln,efmad,efhar

      parameter (rytoev = 13.605826d0)

      if (iarg(1) .gt. 0 .and. iarg(2) .gt. 0) then
        continue
      elseif (iarg(1) .gt. 0 .or. iarg(2) .gt. 0) then
        call rx('need specify two planes, e.g. "qsum 2 4"')
      endif

      pi8 = 32 * datan(1.d0)

C --- Calculate planar average of Madelung and Hartree potentials ---
      do  10  i = 1, nbas
        qplan(iplane(i)) = qplan(iplane(i)) + qt(ipc(i))
        dqplan(iplane(i)) = dqplan(iplane(i)) + dq(ipc(i))
        vintra = 2*qt(ipc(i)) / rmax(ipc(i))
        vmadpl(iplane(i)) = vmadpl(iplane(i))
     .        + (ves(ipc(i)) - vintra) / natpln(iplane(i))
        vharpl(iplane(i)) = vharpl(iplane(i))
     .        + ves(ipc(i)) / natpln(iplane(i))
   10 continue

C --- Find constant-field term to make potential periodic ---
      sumq = 0.d0
      qdipol = 0.d0
      field0 = 0.d0
      imnarg = min0(iarg(1),iarg(2))
      imxarg = max0(iarg(1),iarg(2))
      imnpln = 1

      do  20  i = 1, nplane
        sumq = sumq + qplan(i)
        if ( (i .ge. imnarg) .and. (i .le. imxarg) )
     .          qdipol = qdipol + qplan(i)
        field0 = field0 - alat*plane(i)*qplan(i)
        if (plane(i) .lt. plane(imnpln)) imnpln = i

        do  15  j = 1, nplane
          if (plane(j) .le. plane(i))
     .      qzpl(i) = qzpl(i) + alat*( plane(i) - plane(j) )*qplan(j)
   15   continue

   20 continue

      field0 = field0 + (cellen + alat*plane(imnpln))*sumq
      field0 = field0 / cellen

C --- Printout ---
      if (iprint() .lt. 20) return

      do  25  i = 1, 2
   25 write(lgunit(i),332) pi8*field0 / carea

C  332 format(//,' Constant field that makes V periodic ',f12.8,//,
C     .          2x,'Plane       qt         z*qt       potqz   ',
C     .          '   vmadpl      vharpl      moment',/)
  332 format(//,' Constant field to make V periodic = ',f12.8,//
     .          2x,'Plane       qt         z*qt       potqz   ',
     .          '   vmadpl      vharpl        z   ',/)


C --- Calculate potential from planar charges ---
      do  30  i = 1, nplane
        potqz(i) = (pi8 / carea)*( -qzpl(i)
     .    + alat*( plane(i) - plane(imnpln) )*field0 )

C        do  26  j = 1, 2
C   26   write(lgunit(j),333) i, qplan(i), qzpl(i), potqz(i),
C     .     vmadpl(i), vharpl(i), dqplan(i)
        do  26  j = 1, 2
   26   write(lgunit(j),333) i, qplan(i), qzpl(i), potqz(i),
     .     vmadpl(i), vharpl(i), plane(i)

  333   format(i5,1x,6f12.6)

   30 continue

C --- Calculate dipole (3 ways) ---
      if (iarg(1) .gt. 0 .and. iarg(2) .gt. 0) then
        dipole = potqz(iarg(2)) - potqz(iarg(1))
        dipmad = vmadpl(iarg(2)) - vmadpl(iarg(1))
        diphar = vharpl(iarg(2)) - vharpl(iarg(1))
      else
        dipole = 0
        dipmad = 0
        diphar = 0
      endif

      do  40  i = 1, 2
   40 write(lgunit(i),337) sumq, iarg(1), iarg(2), qdipol, dipole,
     .  rytoev*dipole, dipmad, rytoev*dipmad, diphar, rytoev*diphar

  337 format(//,' Total Charge = ',f12.7,//,' Dipoles for planes ',
     .  i3,' through ',i3,/,' Charge in these planes = ',f12.7,//,
     .  ' Dipole(planar charges) =     ',f12.8,' = ',f12.8,' eV',
     .  ' (from potqz)'/,
     .  ' Dipole(Madelung potential) = ',f12.8,' = ',f12.8,' eV',
     .  ' (from vmadpl)'/,
     .  ' Dipole(Hartree potential) =  ',f12.8,' = ',f12.8,' eV',
     .  ' (from vharpl)')

      if ( dabs(darg(1)) + dabs(darg(2)) .lt. 1.d-5 ) return

C --- Calculate difference between the vacuum and the
C     bulk Fermi level (or VBM) for the surface (3 ways) ----
      wfbpln = rytoev*( dipole - darg(2) )
      wfbmad = rytoev*( dipmad - darg(2) )
      wfbhar = rytoev*( diphar - darg(2) )

C --- Calculate surface work function - difference
C     between vacuum and the supercell Fermi level (2 ways) ---
      wfsmad = rytoev*( vmadpl(iarg(2)) - darg(1) )
      wfshar = rytoev*( vharpl(iarg(2)) - darg(1) )

C --- Calculate Fermi level relative to
C     the bulk Fermi level (or VBM) (3 ways) ---

      efpln = darg(1) - potqz(iarg(1)) + potqz(iarg(2))
     .  - 0.5d0*vmadpl(iarg(2)) - 0.5d0*vharpl(iarg(2)) - darg(2)
      efmad = darg(1) - vmadpl(iarg(1)) - darg(2)
      efhar = darg(1) - vharpl(iarg(1)) - darg(2)

      do  50  i = 1, 2
   50 write(lgunit(i),338) darg(1), rytoev*darg(1), darg(2),
     .    rytoev*darg(2), wfbpln, wfbmad, wfbhar, wfsmad,
     .    wfshar, efpln, rytoev*efpln, efmad, rytoev*efmad,
     .    efhar, rytoev*efhar

  338 format(' Fermi level of supercell = ',f10.6,' = ',f12.8,' eV',/,
     .  ' Fermi level of bulk =      ',f10.6,' = ',f12.8,' eV',//,
     .  ' Separation between bulk Fermi level (or VBM) and vacuum:',//,
     .  '    PE threshold(planar charges) =     ',f12.8,' eV',/,
     .  '    PE threshold(Madelung potential) = ',f12.8,' eV',/,
     .  '    PE threshold(Hartree potential) =  ',f12.8,' eV',//,
     .  ' Separation between supercell Fermi level and vacuum:',//,
     .  '    Work Function(Madelung potential) = ',f12.8,' eV',/,
     .  '    Work Function(Hartree potential) =  ',f12.8,' eV',//,
     .  ' Supercell Fermi level relative to bulk Fermi level ',
     .  '(or VBM):',//,
     .  '    Fermi level(planar charges) =     ',f12.8,
     .  ' = ',f12.8,' eV',/,
     .  '    Fermi level(Madelung potential) = ',f12.8,
     .  ' = ',f12.8,' eV',/,
     .  '    Fermi level(Hartree potential) =  ',f12.8,
     .  ' = ',f12.8,' eV',/)

      end
      subroutine lmplio(recln0,unit,offset,iosw,ido,
     .  plat,planvc,iout,iarg,darg,carg)
C- Control for plane subroutines

C     implicit none
      double precision plat(3,3),planvc(3,3)

C control parameters
      integer recln
      integer unit,offset,ido,iosw,recln0
      integer mxchr,loop0,nlin,nlist,ctlen
      parameter (mxchr=20,ctlen=120)
      character ctbl(mxchr,2)*(ctlen), cxx*1

C output
      character*80 carg(10)
      integer iout, iarg(10)
      double precision darg(10)

C local variables
      integer recl0
      parameter (recl0=72)
      character*1 recrd(0:1)
      character*(recl0) a
      integer i,i1,iosw2,i0,partok,i1mach
      logical T,F,sw
      double precision cost,ddot,hvc,xx

C --- Common blocks ---
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr

      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio
      common /w/ recrd
      data T /.true./, F /.false./

      recoff = offset
      iout = 0
      hvc = dsqrt(ddot(3,planvc(1,3),1,planvc(1,3),1))

C --- Do according to switch iosw ---
      do  1000  iosw2 = iosw, iosw + iosw/3
      optio = iosw2 - 2*(iosw/3)

      goto (1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     .  170), ido+1
      i = recln(recln0)
      return

    1 if (unit .ne. i1mach(1)) rewind unit
      nrecs = 0
      call rdfiln(unit,'#{}%',0,loop0,nlin,xx,100,
     .  xx,nlist,cxx,ctbl,mxchr,a,recrd(recoff),reclen,nrecs)
c      print *,nrecs, ' lines'
      return

C --- IO ---
   10 call getcat(recrd,'io ',' ',f)
      if (.not. noerr) goto 1000
      i0 = partok(recrd(catbeg),'show=','=',sw,' ',-1,0,0,0)
      if (noerr .and. sw) iosw = 3
      i0 = partok(recrd(catbeg),'verbos=','=',i1,' ',-1,2,0,0)
      if (noerr) call pshprt(i1)
      i0 = partok(recrd(catbeg),'wkp=','=',sw,' ',-1,0,0,0)
      if (noerr .and. sw) call wkprnt(1)
      i0 = partok(recrd(catbeg),'iactiv=','=',sw,' ',-1,0,0,0)
      if (noerr) call initqu(sw)
      goto 1000

   20 call getcat(recrd,'help ',' ',f)
      if (noerr) then
        iosw = 0
        return
      endif
      goto 1000

   30 call getcat(recrd,'write ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) print *,
     .  '  ... writes CLASS, SITE and START to log file'
      if (optio .eq. 1 .and. noerr) iout = 3
      goto 1000

   40 call getcat(recrd,'rplan ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) then
        print *,
     .    '  ... rplan h  reads next plane from plan file and puts at h'
        print *,
     .    '  ... rplan (x,y,z) h  same but puts at h / h.(x,y,z)'
      endif
      call getxyz('rplan ',4,darg,i)
      if (optio .eq. 1 .and. noerr) then
        iout = 4
C --- evaluate h from projection ---
        if (i .gt. 0) then
          cost = plat(i,3)/hvc
          darg(1) = darg(1)/cost
        endif
      endif
      goto 1000

   50 call getcat(recrd,'wplan ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) print *,
     .  '  ... wplan n  writes plane n to plan file'
      i0 = partok(recrd(catbeg),'wplan ',' ',iarg,' ',-1,2,0,1)
      if (optio .eq. 1 .and. noerr) iout = 5
      goto 1000

   60 call getcat(recrd,'splan ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) print *,
     .  '  ... splan n  skips to plane n in plan file'
      i0 = partok(recrd(catbeg),'splan ',' ',iarg,' ',-1,2,0,1)
      if (optio .eq. 1 .and. noerr) iout = 6
      goto 1000

   70 call getcat(recrd,'name ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0)
     .  print *, '  ... name ext changes file extension to ext'
      i = partok(recrd(catbeg),'name ',' ',4,carg,1,1,0,1)
      if (optio .eq. 1 .and. noerr) iout = 7
      goto 1000

   80 call getcat(recrd,'stretch ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) then
        print *,
     .  '  ... stretch h  stretches plat(3) by h, or:'
        print *,
     .  '  ... stretch (x,y,z) h  stretches plat(3) by h / h.(x,y,z)'
      endif
      call getxyz('stretch ',4,darg,i)
      if (optio .eq. 1 .and. noerr) then
        iout = 8
C --- evaluate h from projection ---
        if (i .gt. 0) then
          cost = plat(i,3) / hvc
          darg(1) = darg(1)/cost
        endif
      endif
      goto 1000

   90 call getcat(recrd,'qsum ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) then
        print *, '  ... qsum ip1 ip2 [efcell efblk]'
        print *, '      ',
     .  'calculates dipole betw/ planes ip2 and ip2, work func, and Ef'
      endif
      if (optio .eq. 1 .and. noerr) iout = 9
      i0 = partok(recrd(catbeg),'qsum ',' ',darg,' ',4,4,0,1)
      iarg(1) = idnint(darg(1))
      iarg(2) = idnint(darg(2))
      if (dabs(iarg(1)-darg(1)) .gt. 1d-5 .or.
     .    dabs(iarg(2)-darg(2)) .gt. 1d-5 .and. optio .eq. 1)
     .  call fexit(-1,119,'qsum expects first 2 args to be integer',0)
      darg(1) = darg(3)
      darg(2) = darg(4)
      goto 1000

  100 call getcat(recrd,'slide ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) then
        print *,
     .  '  ... slide n1 n2 (x,y,z) h  slides planes n1..n2 by h'
      endif
      i0 = partok(recrd(catbeg),'slide ',' ',iarg,' ',-2,2,0,1)
      call getxyz('slide ',3,darg,i)
      if (optio .eq. 1 .and. noerr) then
        iout = 10
C --- evaluate h from projection ---
        if (i .gt. 0) then
          cost = plat(i,3) / hvc
          darg(1) = darg(1)/cost
        endif
      endif
      goto 1000

  110 call getcat(recrd,'cnam ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0)
     .  print *, '  ... cnam renames classes'
      if (optio .eq. 1 .and. noerr) iout = 11
      goto 1000

  120 call getcat(recrd,'vshift ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) then
        print *,
     .  '  ... vshift n1 n2 v  shifts potential in planes n1..n2 by v'
      endif
      i0 = partok(recrd(catbeg),'vshift ',' ',darg,' ',-3,4,0,1)
      if (optio .ne. 1) goto 1000
      iarg(1) = nint(darg(1))
      iarg(2) = nint(darg(2))
      if (dabs(iarg(1)-darg(1)) .gt. 1d-5 .or.
     .    dabs(iarg(2)-darg(2)) .gt. 1d-5 .and. optio .eq. 1)
     .  call fexit(-1,119,'vshift expects first 2 args to be integer',0)
      darg(1) = darg(3)
      if (optio .eq. 1 .and. noerr) iout = 12
      goto 1000

  130 call getcat(recrd,'rot ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) print *,
     .  '  ... rot x rotates plat by x (fractions of pi) about vc'
      i0 = partok(recrd(catbeg),'rot ',' ',darg,' ',-1,4,0,1)
      if (optio .eq. 1 .and. noerr) iout = 13
      goto 1000

  140 call getcat(recrd,'proj ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0)
     .  print *,
     .  '  ... displays projection of basis vectors into normal plane'
      if (optio .eq. 1 .and. noerr) iout = 14
      goto 1000

  150 call getcat(recrd,'pos ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) then
        print 345
  345   format(
     .'  ... pos [l m n lo hi]  [const expr] logs atom positions for'/
     .'      lattice vectors (0..l,0..m,0..n) projecting betw. (lo,hi)'/
     .'      onto normal.  Logical expr uses vars x,y,z for position'/
     .'      (eg x<2&y>1) and px,py,pz or n for projection onto normal')
        goto 100
      endif
      if (optio .eq. 1 .and. noerr) iout = 15
      if (optio .eq. 1) darg(4)=-99
      if (optio .eq. 1) darg(5)=-99
      i0 = partok(recrd(catbeg),'pos ',' ',darg,' ',5,4,0,1)
      iarg(1) = idnint(darg(1))
      iarg(2) = idnint(darg(2))
      iarg(3) = idnint(darg(3))
      if (i0 .eq. 4) darg(5) = darg(4)
      i0 = partok(recrd(catbeg),'const ',' ',60,carg,1,1,0,0)
      goto 1000

  160 call getcat(recrd,'wsite ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) then
        print *, '  ... writes site data to site file'
        print *, '      Syntax is : wsite [-pad | -ppad] [filename]'
        print *, '        -pad and -ppad add padded'
        print *, '          (or doubly padded) sites to site file.'
        print *, '      If [filename] is not included, site filename'
        print *, '      defaults to ''site''.'
      endif
C     There seems to be a bug in partok.  Do this way
      i = partok(recrd(catbeg),'wsite ',' ',8,carg,1,1,0,0)
      if (carg(1)(1:1) .eq. '-') then
        i = partok(recrd(catbeg),'wsite ',' ',8,carg(2),2,1,0,0)
      endif

      if (optio .eq. 1 .and. noerr) iout = 16
      goto 1000

  170 call getcat(recrd,'q ',' ',f)
      if (.not. noerr) goto 1000
      if (iosw .eq. 0) print *, '  ... quits program'
      if (optio .eq. 1 .and. noerr) call fexit(-1,119,' ',0)
      goto 1000

 1000 continue
      end
      subroutine getxyz(ccat,nch,h,i)
C- Get distance along a line given projection along component
C     implicit none
      integer i,nch
      character*(*) ccat
      double precision h

C --- iolib functions ---
      integer partok, i0

C --- Common blocks ---
      character*1 recrd(0:1)
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr

      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio
      common /w/ recrd

C local variables
      logical T,F
      data T /.true./, F /.false./


      nchoos = nch
      i = -1
      i0 = partok(recrd(catbeg),'x ',' ',h,' ',-1,4,0,1)
      if (noerr) then
        i=1
        goto 85
      endif
      i0 = partok(recrd(catbeg),'y ',' ',h,' ',-1,4,0,1)
      if (noerr) then
        i=2
        goto 85
      endif
      i0 = partok(recrd(catbeg),'z ',' ',h,' ',-1,4,0,1)
      if (noerr) then
        i=3
        goto 85
      endif
      if (nchoos .ge. 4)
     .  i0 = partok(recrd(catbeg),ccat,' ',h,' ',-1,4,0,1)
   85 continue

      end
      subroutine gtplan(npadl,npadr,nbasp,bas,r,ipc,clabl,nplane,planes,
     .  iplane,natpln)
C- Enumerate all planes in basis, make table of planes and pointers to table
C ----------------------------------------------------------------
Ci Inputs
Ci   nbasp,bas
Ci   r:  vector defining plane normal
Ci   clabl, ipc (for printout only)
Co Outputs
Co   nplane,planes:  number of planes and component along r
Co   iplane:         maps basis index to plane index
Co   natpln:         number of atoms in a plane
Cr Remarks
Cr   planes and iplane must be dimensioned sufficiently large
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer npadl,npadr,nbasp,nplane,iplane(1),natpln(1),ipc(1)
      character*8 clabl(1)
      double precision bas(3,nbasp),r(3),planes(1)
C Local parameters
      integer i,j,ii,iprint,lgunit,l,nbas
      double precision ddot,z,tiny,dsqrt,hold,v(3)
      character*4 cc
      parameter (tiny = 1d-4)

      nbas = nbasp - npadl - npadr

C --- Find all planes, generate pointer table ---
C     and determine number of atoms per plane
      nplane = 0
      do  10  ii = 1, nbasp
        call pplan1(1,ii,nbas,npadl,npadr,i)
        z = ddot(3,bas(1,i),1,r,1) / dsqrt(ddot(3,r,1,r,1))
        do  20  j = 1, nplane
          if (dabs(z-planes(j)) .lt. tiny) then
            iplane(i) = j
            natpln(j) = natpln(j) + 1
            goto 10
          endif
   20   continue
        nplane = nplane+1
        planes(nplane) = z
        iplane(i) = nplane
        natpln(nplane) = 1
   10 continue

C --- Printout ---
      if (iprint() .lt. 30) return
      do  30  i = 1, 2
   30 write(lgunit(i),333) nplane, r
  333 format(/' Gtplan:',i3,' different planes normal to (', 3f9.5,')',
     .  //,'  ib      Class  Plane',6x,
     .  'x-x.h     y-y.h     z-z.h      h       del h',/)

      hold = planes(iplane(nbasp)) - dsqrt(ddot(3,r,1,r,1))
      do  50  ii = 1, nbasp
        call pplan1(1,ii,nbas,npadl,npadr,i)
        z = ddot(3,bas(1,i),1,r,1) / dsqrt(ddot(3,r,1,r,1))
        call dcopy(3,bas(1,i),1,v,1)
        call daxpy(3,-z/dsqrt(ddot(3,r,1,r,1)),r,1,v,1)
        if (i .le. nbas) then
          cc = ' '
        elseif (i .le. nbas+npadl) then
          cc = '(L) '
        else
          cc = '(R) '
        endif
        do  40  l = 1, 2
   40   write(lgunit(l),334) i, cc, ipc(i), clabl(ipc(i)), iplane(i),
     .  (v(j), j=1,3), planes(iplane(i)), planes(iplane(i))-hold
  334   format(i4,a4,i4,':',a4,i4,3x,5f10.5)
        hold = planes(iplane(i))
   50 continue

      end
      subroutine fixcl(nbas,nclspp,ipc,clabl,iwk)
C     implicit none
      integer nbas,nclspp,ipc(16),iwk(16)
      character*8 clabl(8)

C Local variables
      integer ib,ic,i,j,ndel

c      nclspp = 20
c      clabl(17)='a'
c      clabl(18)='b'
c      clabl(19)='c'
c      clabl(20)='d'
c      IPC(2) = 18
c      IPC(3) = 20
c      IPC(14) = 18
c      IPC(15) = 20
c      PRINT *, IPC

C --- Eliminate replication of class names ---
      do  10  ic = 1, nclspp
        do  13  j = 1, ic-1
          if (clabl(j) .ne. clabl(ic)) goto 13
          print 432, clabl(ic), ic, j
  432     format(' fixcl:  reset class ',a4,', number',i3,
     .           ' to earlier class',i3)
          do  20  ib = 1, nbas
            if (ipc(ib) .eq. ic) ipc(ib) = j
   20     continue
          goto 10
   13   continue
   10 continue

C --- Eliminate unused classes ---
      call icopy(nbas,ipc,1,iwk,1)
      call ishell(nbas,iwk)
      ib = 1
   40 ib = ib+1
        if (iwk(ib) .eq. iwk(ib-1) .or. iwk(ib) .eq. iwk(ib-1)+1)
     .    goto 49
        print *, ' fixcl:  discarding unused class(es)',
     .    iwk(ib-1)+1, ' through ', iwk(ib)-1

C --- Collapse tables clabl and ipc --
        ndel =  iwk(ib) - iwk(ib-1) - 1

c        PRINT *, CLABL
        do  45  ic = iwk(ib), nclspp
   45   clabl(ic-ndel) = clabl(ic)
        nclspp = nclspp - ndel
c        PRINT *, CLABL

        do  48  i = 1, nbas
   48   if (ipc(i) .ge. iwk(ib)) ipc(i) = ipc(i)-ndel
        call icopy(nbas,ipc,1,iwk,1)
        call ishell(nbas,iwk)

c        PRINT *, IPC
c        PRINT *, IWK

   49 if (ib .lt. nbas) goto 40
C Highest class is no larger than iwk(nb) ...
      nclspp = min(nclspp,iwk(nbas))

C For debugging, check ipc
      do  50  ib = 1, nbas
   50 if (ipc(ib) .gt. nclspp)
     .    call fexit(-1,119, 'fixcl: bad ipc',0)

c      PRINT *, IPC
      end
      subroutine pplan1(mode,ib,nbas,npadl,npadr,ibp)
C- Return a permuted site index that orders sites L,1:nbas,R
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 for forward permutation, 1 for reverse
Ci   ib    :unpermuted site index
Ci   nbas  :size of basis without padding
Ci   npadl :number of L padding sites 
Ci   npadr :number of R padding sites
Co Outputs
Co   ibp   :forward permutation reorders sites as follows:
Co         :ib  =  1:nbas,1+nbas:npadl+nbas,1+nbas+npadl:npadr+nbas+npadl
Co         :ibp =  npadl+1:npadl+nbas,1:npadl,1:nbas,1:npadr
Co         :reverse permutation not implemented
Co         :ib  =  1:nbas,1+nbas:npadl+nbas,1+nbas+npadl:npadr+nbas+npadl
Co         :ibp =  npadl+1:npadl+nbas,1:npadl,1:nbas,1:npadr
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   31 Oct 03 First created
C ----------------------------------------------------------------------
C     implicit none
      integer mode,ib,nbas,npadl,npadr,ibp

      if (mode .eq. 0) then
        if (ib .le. nbas) then
          ibp = ib+npadl
        elseif (ib .le. nbas+npadl) then
          ibp = ib - nbas
        else
          ibp = ib
        endif
      else
        if (ib .le. npadl) then
          ibp = ib+nbas
        elseif (ib .le. nbas+npadl) then
          ibp = ib - npadl
        else
          ibp = ib
        endif
      endif
      end

      

