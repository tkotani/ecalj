      subroutine ovmin(sovmin,nbas,nbasp,alat,plat,rmax,rmt,dclabl,
     .ips,mode,z,ontab,oiax,pos,iprtbl)

C- Check volume and sphere overlaps, optionally minimizing them
C ----------------------------------------------------------------
Ci Inputs
Ci   sovmin: a set of modifiers, with the syntax
Ci          -mino[:dxmx=#][:xtol=#][:style=#]:site-list
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   rmax  :potential radius, in a.u.
Ci   rmt   :augmentation radius, in a.u.
Ci   dclabl:class name, packed as a real number
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   mode:  vector of length 3 governing how pos shortened (see shorps)
Ci   z     :nuclear charge
Ci   w(ontab):ntab(ib)=# pairs in iax table preceding ib (pairc.f)
Ci   w(oiax):neighbor table containing pair information (pairc.f)
Ci   pos   :basis vectors
Ci   iprtbl: nonzero if to call ovlchk and print table of overlaps
Co Outputs
Co   Sphere overlaps are printed out
Cr Remarks
Cr   rmt(1)  not used now
Cu Updates
Cu   22 Oct 02  weight ES-ES and atom-ES overlaps differently when
Cu              minimizing sphere overlap positions
Cu    9 Dec 98  replace call to frpmin with call to gradzr.
Cu    8 Sep 98  small patches in minimizing algorithm
Cu   24 Nov 97  changed ovmin to call fovlp for fast execution
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nbasp,ontab,oiax,iprtbl
C     integer ntab(nbas+1),iax(niax,1)
      double precision plat(3,3),pos(3,nbasp),rmax(1),rmt(1),z(1),alat
      double precision dclabl(1)
      integer ips(1),mode(3)
      character sovmin*(*)
C static:
      double precision alato,plato(9),xx
      integer nbaso,nbaspo,mxlst,nlst,opos,oz,ormax,oips,modeo(3),
     .ontabo,oiaxo,novl
      parameter (mxlst=256)
      integer ilst(mxlst)
      common /ovstat/ plato,alato,nbaso,nbaspo,nlst,ilst,opos,oz,
     .ormax,oips,modeo,ontabo,oiaxo
C heap:
      integer w(1)
      common /w/ w
C Local parameters
      double precision fovl,xtol,gtol,dxmn,dxmx,ovcall,fovmx
      double precision wk(0:27)
      integer i1mach,isw,op,og,ir,i,j,j1,j2,ls,m,lstyle,
     .iv,parg,olist,nlstc,mxint,nclass,ib,ic,iclbsj,maxit,ipr,n
      character dc*1
      external ovcall,mxint

C --- Print out positions and overlaps ---
      call getpr(ipr)
      if (iprtbl .gt. 0) call ovlchk(nbas,nbasp,pos,alat,rmax,0d0,
     .dclabl,ips,mode,plat,fovmx,xx)
      call fovlp(1,nbas,w(ontab),w(oiax),plat,pos,ips,alat,rmax,z,6d0,
     .1d0,.75d0,.5d0,fovmx,fovl,novl)
      if (novl .eq. 0) novl = 1
      if (ipr .ge. 10 .or. iprtbl .gt. 0)
     .call awrit3('%N OVMIN:     fovl = %;6g   <ovlp> = %;1d%%'//
     .'   max ovlp = %;1d%%',' ',80,
     .i1mach(2),fovl/novl,(fovl/novl)**(1/6d0)*100,fovmx*100)

C --- Minimize overlaps wrt positions in list ---
      if (sovmin .ne. ' ') then
C   ... Default values for gradzr call
        xtol = 2d-4
        gtol = 1d-5
        dxmn = 1d-6
        dxmx = .10d0
        maxit = 20
        isw = 10051

        ls = len(sovmin)
        j1 = 1
        dc = sovmin(j1:j1)
        j1 = j1+1
        lstyle = 0

C   ... Return here to resume parsing for arguments
   40   continue
        call nwordg(sovmin,0,dc//' ',1,j1,j2)

C   ... Parse special arguments
        if (sovmin(j2+1:j2+1) .ne. ' ')  then
          m = j1-1
          i = parg('dxmx=',4,sovmin,m,ls,dc,1,1,iv,dxmx)
          m = j1-1
          i = parg('xtol=',4,sovmin,m,ls,dc,1,1,iv,xtol)
          m = j1-1
          i = parg('style=',2,sovmin,m,ls,dc,1,1,iv,lstyle)
          m = j1-1
          i = parg('maxit=',2,sovmin,m,ls,dc,1,1,iv,maxit)
          j1 = j2+2
          goto 40
        endif

C   ... List of all sites to move
        if (lstyle .gt. 0) then
          nclass = mxint(nbas,ips)
          call defi(olist, nclass)
          call clist(lstyle,sovmin(j1:j2+1),dclabl,z,nclass,nlstc,
     .    w(olist))
          nlst = 0
          do  12  i = 1, nlstc
            ic = w(olist+i-1)
            do  14  j = 1, nbas
              ib = iclbsj(ic,ips,-nbas,j)
              if (ib .lt. 0) goto 12
              nlst = nlst+1
              ilst(nlst) = ib
   14       continue
   12     continue
        elseif (sovmin(j1:j1+1) .eq. 'z ' .or.
     .  sovmin(j1:j1+1) .eq. 'Z ') then
          nlst = 0
          do  10  ib = 1, nbasp
            ic = ips(ib)
            if (z(ic) .eq. 0) then
              nlst = nlst+1
              ilst(nlst) = ib
            endif
   10     continue
        else
          call mkilst(sovmin(j1:),nlst,ilst)
        endif
        call awrit2(' min wrt:  %n:1i',' ',80,i1mach(2),nlst,ilst)
        call awrit3(' setup:     xtol = %,2g   dxmx = %,2g   maxit = %i'
     .  ,' ',80,i1mach(2),xtol,dxmx,maxit)
        if (nlst .le. 0) then
          print *, 'no sites in list ... no minimization'
          return
        endif

C  ...  set up static block for ovcall
        alato = alat
        nbaso = nbas
        nbaspo = nbasp
        ontabo = ontab
        oiaxo = oiax
        call defrr(opos,3*nbasp)
        call dpcopy(pos,w(opos),1,3*nbasp,1d0)
        nclass = mxint(nbas,ips)
        call defrr(oz,nclass)
        call dpcopy(z,w(oz),1,nclass,1d0)
        call defrr(ormax,nbasp)
        call dpcopy(rmax,w(ormax),1,nbasp,1d0)
        call defi(oips,nbasp)
        call icopy(nbasp,ips,1,w(oips),1)
        call icopy(3,mode,1,modeo,1)
        call dpcopy(plat,plato,1,9,1d0)

C  ...  initialization for gradzr call
        n = 3*nlst
        call defrr(op,-10*n)
        call defrr(og,n)
        ir = 0
        do  20  i = 1, nlst
          j = ilst(i)
          call dpscop(w(opos),w(op),3,3*j-2,3*i-2,1d0)
   20   continue
        xx = ovcall(n,0d0,w(op),ir)
        call pshpr(ipr-5)
   22   call gradzr(n,w(op),w,dxmn,dxmx,xtol,gtol,1.0d0,wk,isw,ir)
        xx = ovcall(n,0d0,w(op),ir)
        if (ir .lt. 0) goto 22
        call poppr

C ...   Update positions
        do  30  i = 1, nlst
          j = ilst(i)
          call dpscop(w(op),pos,3,3*i-2,3*j-2,1d0)
   30   continue

C --- Print out updated positions and overlaps ---
        print '(/'' OVMIN:  updated site positions:'')'
        if (iprtbl .gt. 0) call ovlchk(nbas,nbasp,pos,alat,rmax,0d0,
     .  dclabl,ips,mode,plat,fovmx,xx)
        call fovlp(1,nbas,w(ontab),w(oiax),plat,pos,ips,alat,rmax,z,6d0,
     .  1d0,.75d0,.5d0,fovmx,fovl,novl)
        if (novl .eq. 0) novl = 1
        if (ipr .ge. 10)
     .  call awrit3(' minimized: fovl = %;6g   <ovlp> = %;1d%%'//
     .  '   max ovlp = %;1d%%',' ',80,
     .  i1mach(2),fovl/novl,(fovl/novl)**(1/6d0)*100,fovmx*100)
      endif

      end subroutine ovmin

      double precision function ovcall(n,x,p,ir)

C- Generic function call for projection grad fovl in a spec'd direction
Ci x,p,ir see gradzr
C     implicit none
      integer ir,n
      double precision x,p(3*n)
C static:
      double precision alato,plato(9)
      integer nbaso,nbaspo,mxlst,nlst,opos,oz,ormax,oips,modeo(3),oposb,
     .ontabo,oiaxo,novl,novlp,novlm
      parameter (mxlst=256)
      integer ilst(mxlst)
      common /ovstat/ plato,alato,nbaso,nbaspo,nlst,ilst,opos,oz,
     .ormax,oips,modeo,ontabo,oiaxo
C heap:
      integer w(1)
      common /w/ w
C Local
      logical cmdopt
C     character*8 clabl(10)
      integer j,i,ix,ipr,lgunit,novl0
      double precision fovl,ddot,dx,val,fovp,fovm,pos(3),xx,fov0
      character*120 outs
      parameter (dx=1d-4)

C ... Save pos, other initialization
      call getpr(ipr)
      call defrr(oposb,3*nbaspo)
      call dpcopy(w(opos),w(oposb),1,3*nbaspo,1d0)
      call pshpr(0)

      do  12  i = 1, nlst
        j = ilst(i)
        call dpscop(p,w(opos),3,3*i-2,3*j-2,1d0)
   12 continue

      call ovlchk(nbaso,nbaspo,w(opos),alato,w(ormax),0d0,0d0,
     .w(oips),modeo,plato,fovl,xx)
      call fovlp(1,nbaso,w(ontabo),w(oiaxo),plato,w(opos),w(oips),alato,
     .w(ormax),w(oz),6d0,1d0,.75d0,.5d0,xx,fovl,novl)

      if (fovl .eq. 0) then
        print *, 'ovmin: no spheres overlap:'
        call poppr
C        call fovlp(1,nbaso,w(ontabo),w(oiaxo),plato,w(opos),w(oips),
C     .    alato,w(ormax),w(oz),6d0,1d0,.75d0,.5d0,xx,fovl,novl)
        call ovlchk(nbaso,nbaspo,w(opos),alato,w(ormax),0d0,0d0,
     .  w(oips),modeo,plato,fovp,xx)
        if (cmdopt('--wpos=',7,0,outs))
     .  call iopos(.true.,0,outs(8:),nbaso,w(opos))
        call rx('ovmin: no spheres overlap')
      endif

C ... Gradient of fovl wrt pos
      do  20  i = 1, nlst
        j = ilst(i)
        call fovlp(j,j,w(ontabo),w(oiaxo),plato,w(opos),w(oips),
     .  alato,w(ormax),w(oz),6d0,1d0,.75d0,.5d0,xx,fov0,novl0)
        do  20  ix = 1, 3
          val = p(3*i-3+ix)
          call dvset(w(opos),3*j-3+ix,3*j-3+ix,val+dx)
C        call ovlchk(nbaso,nbaspo,w(opos),alato,w(ormax),0d0,0d0,
C     .    w(oips),modeo,plato,fovp,xx)
          call fovlp(j,j,w(ontabo),w(oiaxo),plato,w(opos),w(oips),
     .    alato,w(ormax),w(oz),6d0,1d0,.75d0,.5d0,xx,fovp,novlp)
          call dvset(w(opos),3*j-3+ix,3*j-3+ix,val-dx)
C        call ovlchk(nbaso,nbaspo,w(opos),alato,w(ormax),0d0,0d0,
C     .    w(oips),modeo,plato,fovm,xx)
          call fovlp(j,j,w(ontabo),w(oiaxo),plato,w(opos),w(oips),
     .    alato,w(ormax),w(oz),6d0,1d0,.75d0,.5d0,xx,fovm,novlm)
          call dvset(w(opos),3*j-3+ix,3*j-3+ix,val)
          fovp = fovl + 2*(fovp-fov0)
          fovm = fovl + 2*(fovm-fov0)
          p(n+3*i-3+ix) = dlog(fovp/fovm)/2/dx
*       print *, '... i,j,ix=',i,j,ix,fovp,fovm,p(n+3*i-3+ix)
   20 continue
      ovcall = ddot(n,p(n+1),1,p(2*n+1),1)
      if (ipr .ge. 50) then
        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
     .  'glf.x %;4g',' ',80,lgunit(1),x,fovl/novl,dlog(fovl/novl),
     .  dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
     .  'glf.x %;4g',' ',80,lgunit(2),x,fovl/novl,dlog(fovl/novl),
     .  dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
        do  30  i = 1, nbaspo
          call dpscop(w(opos),pos,3,3*i-2,1,1d0)
          write(lgunit(2),140) pos
  140     format(3f12.6)
   30   continue
        call query('continue',-1,0)
      endif

C      call prmx('grad fovl',p(1+n),n,n,1)
C      call prmx('pos now',w(opos),3,3,nbaspo)

C ... restore pos
      call dpcopy(w(oposb),w(opos),1,3*nbaspo,1d0)
      call rlse(oposb)
      call poppr

      end function ovcall


