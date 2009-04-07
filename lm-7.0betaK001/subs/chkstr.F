      subroutine chkstr(sctrl,slat,sspec,sstr,sarray)
C- Plots screened envelope functions for checking
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbasp nl
Ci     Stored:
Ci     Passed to:
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat avw ocy
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: opos
Ci     Stored:
Ci     Passed to:
Co Outputs
Cl Local variables
Cl   mode  : plotting mode:
Cl         :  1  tabulate screened envelope in a plane
Cl         :  2  tabulate screened envelope on a line
Cl         :  3  check one-center expansion
Cr Remarks
Cr   Routine does one of:
Cr   1.  tabulate screened envelope in a plane for contour plots
Cr       Invoke with --plot:con
Cr   2.  tabulate screened envelope on a line
Cr       Invoke with --plot:line[:v1x,v1y,v1z,v2x,v2y,v2z]
Cr   3.  check one-center expansion
Cr       Invoke with --plot:onec
Cu Updates
Cu   23 Jul 08 (S Lozovoi)) Adapted to l-dependent augmentation radii
Cu             and species-dependent lmax
Cu   06 Aug 06 Revised to work with 2-kappa basis
Cu   19 May 04 Added --plot options
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision sctrl(1),slat(1),sspec(1),sstr(1),sarray(1)
C ... Local parameters
      double precision plat(3,3),alat,avw,ekap(20),ckbas,cksumf,tstmax
      integer npmax,nlmy,nvmax
      parameter (npmax=25,nlmy=16,nvmax=1001)
      double precision slj(25,2),slk(25,2),
     .  vvaln(nlmy,nlmy),vvalj(nlmy,nlmy),
     .  dvaln(nlmy,nlmy),dvalj(nlmy,nlmy)
      double precision reslin(nvmax,12)
      double precision resk(-npmax:npmax,-npmax:npmax),xx,ctr(3)
      equivalence (resk,reslin)
      double precision xp(3,122),xpt(122,3),wp(122),yl(122,nlmy),
     .  rsq(122),v1(6),v2(3),radius,
     .  xmrp(3),fi(10),gi(10),gi0(10),hcr(10)
      equivalence (v1(4),v2)
      integer oalpha,obas,ocy,offR,offRp,oiax,ontab,ortab,ob,os,
     .  oips,olmx,ohcr,ohalp
      integer a2vec,fopn,i,iat,iatp,ifi,ikap,imax,ip,itral,ival,ixi,
     .  j,j1,j2,ll,il,ila,ilb,lmaxw,mode,nbasp,niax,nitab,nkap,nl,nlf,
     .  np,nttab,nvec,lio,nbasf,stdo,nglob,lio23,lio45,getdig,ncplx,
     .  loka,isw,is,nclus,nlma,nlmb,la,lb,nspec,nds,irot
      logical ltmp,iostr,iostr1,cmdopt
      character outs*120,out2*120,dc*1,dc2*1
      parameter (niax=10)
c for error printing
      integer ier1,ier2,ixx
      double precision err

C ... Heap
      integer w(1)
      common /w/ w

      stdo = nglob('stdo')
      call upack('ctrl nbasp nspec nl',sctrl,nbasp,nspec,nl,0,0)
      call upack('lat alat plat avw',slat,alat,plat,avw,0,0)
      call upack('array opos oips',sarray,obas,oips,0,0,0)
      call upack1('lat ocy',slat,ocy)

C ... Read some strux data from from disk
      call info0(10,1,0,' CHKSTR: check or plot strx from file STR ...')
      lio = 0
      ltmp = iostr1(lio,'STR',itral,nttab,nitab,nlf,nbasf,nkap,lmaxw,xx,
     .  ifi)
      lio23 = mod(lio/100,100)
      lio45 = mod(lio/10000,100)
      ncplx = 1+getdig(lio23,0,2)
      loka = isw(lio45 .lt. 2)

      call info8(10,0,0,
     .  ' file contains: '//
     .  '%?#(n==0)#2nd generation ##%-1j'//
     .  '%?#(n==1)#NMTO ##%-1j'//
     .  '%?#(n==2)#Screened ##'//
     .  '%?#n#sdot#strux#'//
     .  '%?#(n<2)#, OKA defs##'//
     .  '%?#(n==0)# (1-center coffs)##'//
     .  '%?#(n>1)#, %-1j%i energies##',
     .  lio45,
     .  getdig(lio23,3,2),
     .  lio45,
     .  getdig(lio23,2,2),
     .  nkap,
     .  nkap,
     .  ekap,
     .  0)
      call isanrg(nbasf,nbasp,nbasp,'file:CHKSTR:','nbas',.true.)
      ckbas = cksumf(w(obas),3*nbasp)
      ltmp = iostr(lio+8,'STR',nl,nbasp,nkap,ekap,itral,ckbas,lmaxw,
     .             nitab,oalpha,oiax,ontab,ob)
      call info5(10,0,0,
     .  ' %i pairs in neighbor table   nlf = %i   E =%n:1;4,4d%N',
     .  w(ontab+nbasp),nlf,nkap,ekap,0)
      nds = nlf**2

C     Setup for rotating to value, slope.  irot=0 for no rotation
c     irot = 0
      irot = 1
      call query('convert to value-slope functions (0=no)',2,irot)
      if (irot .ne. 0) irot = 100
      if (irot .eq. 100) then
        call sp2cls('spec lmxb',sspec,0,1,1,nspec,olmx)
        call sp2cls('spec hcr',sspec,0,nl,1,nspec,ohcr)
        call defdr(ohalp,-nds*nbasp*nkap**2)
        call makalp(nl,nds,nbasp,nkap,ekap,w(ohcr),100*loka+2,
     .    w(oips),w(olmx),w,w(ohalp))
      else
        call defrr(ohalp,1)
      endif

      if (loka .eq. 1) then
        call dscal(nds**2*nitab*nkap**2,-1d0,w(ob),1)
      endif

C     Convert 1-center expansion to strux matrix
      nttab = w(ontab+nbasp)
      call defdr(os,nds**2*nitab*nkap**2)
      call hanb2s(irot+0,nds,nbasp,nkap,w(ontab),w(oiax),w(oalpha),
     .  w(ohalp),w(os),w(ob))
C     1-center also for (val,slo)
      if (irot .eq. 100) then
      call hanb2s(irot+2,nds,nbasp,nkap,w(ontab),w(oiax),w(oalpha),
     .  w(ohalp),w(ob),w)
      endif

C ... Initialize some variables and defaults
      ltmp = .false.
      mode = 0
      tstmax = alat
      v1(1) = -tstmax
      v1(2) = 0
      v1(3) = 0
      v2(1) = tstmax
      v2(2) = 0
      v2(3) = 0
      nvec = 501

C ... get plotting(/testing) mode
      if (cmdopt('--plot',6,0,outs)) then
        out2 = outs(7:)
        dc = out2(1:1)
        if (dc .ne. ' ') then
C     ... Return here to resume parsing for arguments
          j2 = 0
   50     continue
          j2 = j2+1
          if (out2(j2:j2) .eq. dc) goto 50
          j1 = min(len(out2),j2)
          call nwordg(out2,0,dc//' ',1,j1,j2)
          if (j2 .ge. j1) then
            if (.false.) then
            elseif (out2(j1:j1+2) .eq. 'con')  then
              mode = 1
            elseif (out2(j1:j1+3) .eq. 'line')  then
              mode = 0
              dc2 = out2(j1+4:j1+4)
              if (dc2 .ne. ' ') then
                ip = 6
                i = a2vec(out2,len(out2),ip,4,', '//dc,3,2,6,fi,v1)
                if (i .ne. 6) goto 52
                j2 = ip
              endif
            elseif (out2(j1:j1+3) .eq. 'onec')  then
              mode = 2
            else
              goto 52
            endif
            goto 50
   52       continue
            call rxs2('chktr: failed to parse --plot option: "',
     .        out2(1:j2),'%a ..."')
          endif
        endif
      endif

C ... Make rtab, plat in atomic units
      nttab = ival(w(ontab),nbasp+1)
      call defrr(ortab,3*nttab)
      call dpzero(ctr,3)
      call mkrtab(100,alat,plat,w(obas),w(oiax),nttab,ctr,w(ortab))
C     call prmx('rtab',w(ortab),3,3,nttab)
      call dscal(9,alat,plat,1)
      call dscal(3*nbasp,alat,w(obas),1)

C ... Scaling by avw
      if (loka .eq. 0) avw = 1
      call dscal(3*nttab,1/avw, w(ortab),1)
      call dscal(nkap,   avw**2,ekap,   1)

      if (mode .eq. 0) goto 70
      if (mode .eq. 1) goto 60
      if (mode .eq. 2) goto 20
      call rxi('invalid mode',mode)

C --- Contour plot ---
   60 continue
      call query('half length of square (a.u.)',4,tstmax)
      imax = npmax
      call query('number of mesh points=',2,imax)
      iat = 1
   63 call query('atom=',2,iat)
      if (iat .gt. nbasp) then
        print *, 'atom cannot be larger than nbas=',nbasp
        goto 63
      endif
   62 ilb = 3
c  62 ilb = 1
      call query('orbital L to plot=',2,ilb)
      ixi = 1
      offR  = ival(w(ontab),iat)
      nlmb = ival(w(oiax),niax*offR+9)
      if (ilb .gt. nlmb) then
        print *, 'L cannot be larger than nlmb=',nlmb
        goto 62
      endif
   64 ikap = 1
      if (nkap .gt. 1) call query('ikap=',2,ikap)
      if (ikap .gt. nkap) then
        print *, 'atom cannot be larger than nkap=',nkap
        goto 64
      endif
      lb = ll(nlmb)
      is = ival(w(oips),iat)
      call upack('spec hcr',sspec,is,hcr,0,0,0)
      call info5(10,1,0,' Atom %i has hcr=%n:1;4,4d',iat,lb+1,hcr,0,0)
      call info5(0,0,1,
     .  ' plotting plane from (x,y) = %da.u. = %d*alat%N%16f'//
     .  'to   (x,y) = %da.u. = %d*alat',-tstmax,-tstmax/alat,tstmax,
     .  tstmax/alat,0)

C --- Contour plot ---
      do  15  i = -imax, imax
        do  16  j = -imax, imax
          xmrp(1) =  (i*tstmax)/imax
          xmrp(2) =  (j*tstmax)/imax
          xmrp(3) =  1d-4
          call dscal(3,1/avw,xmrp,1)
          call strck(nlf**2,iat,ixi,w(ortab),ekap,ikap,nkap,nbasp,
     .      w(oalpha),w(ohalp),w(oiax),w(ontab),w(os),w(ob),nitab,
     .      w(ocy),xmrp,nlf-1,irot+loka,slk,slj)
C          resk(i,j) = ddot(nkap,slk(1,1),25,w(ohalp),nlf**2*nbasp)
           resk(i,j) = max(min(slk(ilb,1),999d0),-99d0)
   16   continue
   15 continue
      ifi = fopn('PLOT')
      rewind ifi
      call ywrm(0,' ',1,ifi,'(%5,6g)',resk,0,2*imax+1,2*imax+1,2*imax+1)

      return

C --- One-center expansion ---
C     Input radius; direction angles make from random number generator
   20 continue
C     pi = 4d0*datan(1d0)
C     y00 = 1/sqrt(4*pi)
C     call ran1in(12)
C --- Test expansion for R' corresponding to iax(ixi) ---
      print 926
  926 format(/
     .  ' Compare numerically integrated YL expansions of H^a and its'/
     .  ' 1-center expansion on a neighboring sphere surface.  These'/
     .  ' should be equal apart from numerical integration errors.')
      if (irot .eq. 0) print 927
      if (irot .ne. 0) print 929
  927 format(
     .  ' When evaluated at radius=hcr,'/
     .  ' function should be H^0 on (head sphere)'/
     .  ' and vanish on tail spheres.'/)
  929  format(
     .  ' When evaluated at radius=hcr,'/
     .  ' function should be 1 on head sphere (0 for ikap=2) '/
     .  ' and vanish on tail spheres.'/)

c     iat = 2
      iat = 1
   23 call query('atom=',2,iat)
      if (iat .gt. nbasp) then
        print *, 'atom cannot be larger than nbas=',nbasp
        goto 23
      endif
   24 ikap = 1
c  24 ikap = 2
      if (nkap .gt. 1) call query('ikap=',2,ikap)
      if (ikap .gt. nkap) then
        print *, 'ikap cannot be larger than nkap=',nkap
        goto 24
      endif
c  25 ixi = 2
   25 ixi = 1
      call query('ixi= (1 for head, 2 for 1st NN, etc)',2,ixi)
      offR  = ival(w(ontab),iat)
      nclus = ival(w(ontab+1),iat) - offR
      nlma = ival(w(oiax),niax*offR+9)
      la = ll(nlma)
      if (ixi .gt. nclus) then
        print *, 'ixi cannot be larger than nclus=',nclus
        goto 25
      endif
      offRp = ixi-1 + offR
      iatp = ival(w(oiax),niax*offRp+2)
      is = ival(w(oips),iatp)
C ... la, nlma - head of the cluster, lb, nlmb - expansion sphere 
      nlmb = ival(w(oiax),niax*offRp+9)
      lb = ll(nlmb)
      call info5(10,0,0,' Neighbor ixi=%i corresponds to atom %i',
     .  ixi,iatp,0,0,0)
      call upack('spec hcr',sspec,is,hcr,0,0,0)
      call info5(10,0,0,' Atom %i has lmax=%i hcr=%n:1;4,4d',
     .  iatp,lb,lb+1,hcr,0)

c Uncomment next 4 lines and comment line 'radius = hcr(il+1)' further down
c if one wishes to integrate with fixed sphere radius introduced interactively
c     radius = hcr(1)
c  21 continue
c     call query('radius=',4,radius)
c     call info5(10,1,0,' ... using neighbor ixi=%i (jb=%i), radius=%d'
c    .  ,ixi,iatp,radius,0,0)

      call info5(10,1,0,
     .  ' ... using neighbor ixi=%i (jb=%i), radius=%n:1;4,4d'
     .  //' ikap=%i',ixi,iatp,lb+1,hcr,ikap)

C --- One-center expansion by brute force integration ---
C     Integration takes place on sphere R' = iax(2,ixi)
c     call fpiint(-60,0,np,xp,wp)
      call fpiint(-122,0,np,xp,wp)
C     Normalized spherical harmonics for all points on sphere R'
      call dmcpy(xp,1,3,xpt,np,1,np,3)
      call ropyln(np,xpt(1,1),xpt(1,2),xpt(1,3),lb,np,yl,rsq)
C     For each head channel, do
      call dpzero(vvaln,nlmy*nlf**2)
      call dpzero(vvalj,nlmy*nlf**2)
C     Integrate H^a and its one-center expansion on sphere R'
      do il = 0, lb
        radius = hcr(il+1)
        do  ip = 1, np
          xmrp(1) =  radius*xp(1,ip)
          xmrp(2) =  radius*xp(2,ip)
          xmrp(3) =  radius*xp(3,ip)
          call dscal(3,1/avw,xmrp,1)
          call strck(nlf**2,iat,ixi,w(ortab),ekap,ikap,nkap,nbasp,
     .      w(oalpha),w(ohalp),w(oiax),w(ontab),w(os),w(ob),nitab,
c    .      w(ocy),xmrp,nlf-1,irot+loka,slk,slj)
     .      w(ocy),xmrp,la,irot+loka,slk,slj)
C       Integration over the surface of sphere a_R' means that
C       K^a_RL(x on sphere R') = sum_L' vvaln_L'(a_R') Y_L'(a_R')
C       vvaln holds K^a(x-R'); vvalj holds 1-center expansion to K
          do  ilb = il*il+1, (il+1)**2
            do  ila = 1, nlma
              vvaln(ilb,ila) = vvaln(ilb,ila) + 
     .          slk(ila,1)*wp(ip)*yl(ip,ilb)
              vvalj(ila,ilb) = vvalj(ila,ilb) + 
     .          slj(ila,1)*wp(ip)*yl(ip,ilb)
            enddo
          enddo
        enddo
      enddo
C     Calculate corresponding radial derivatives
      call dpzero(dvaln,nlmy*nlf**2)
      call dpzero(dvalj,nlmy*nlf**2)
      do il = 0, lb
        radius = hcr(il+1)
        do  ip = 1, np
          xmrp(1) =  (radius + .0005d0)*xp(1,ip)
          xmrp(2) =  (radius + .0005d0)*xp(2,ip)
          xmrp(3) =  (radius + .0005d0)*xp(3,ip)
          call strck(nlf**2,iat,ixi,w(ortab),ekap,ikap,nkap,nbasp,
     .      w(oalpha),w(ohalp),w(oiax),w(ontab),w(os),w(ob),nitab,
c    .      w(ocy),xmrp,nlf-1,irot+loka,slk(1,1),slj(1,1))
     .      w(ocy),xmrp,la,irot+loka,slk(1,1),slj(1,1))
          xmrp(1) =  (radius - .0005d0)*xp(1,ip)
          xmrp(2) =  (radius - .0005d0)*xp(2,ip)
          xmrp(3) =  (radius - .0005d0)*xp(3,ip)
          call strck(nlf**2,iat,ixi,w(ortab),ekap,ikap,nkap,nbasp,
     .      w(oalpha),w(ohalp),w(oiax),w(ontab),w(os),w(ob),nitab,
c    .      w(ocy),xmrp,nlf-1,irot+loka,slk(1,2),slj(1,2))
     .      w(ocy),xmrp,la,irot+loka,slk(1,2),slj(1,2))
          do  ilb = il*il+1, (il+1)**2
            do  ila = 1, nlma
              dvaln(ilb,ila) = dvaln(ilb,ila) +
     .          (slk(ila,1)-slk(ila,2))/.001d0*wp(ip)*yl(ip,ilb)
              dvalj(ilb,ila) = dvalj(ilb,ila) +
     .          (slj(ila,1)-slj(ila,2))/.001d0*wp(ip)*yl(ip,ilb)
            enddo
          enddo
        enddo
      enddo

c ... make H0
      if (ixi .eq. 1 .and. irot .eq. 0) then
        do il = 1, lb+1
          radius = hcr(il)
          call besslr(ekap(ikap)*(radius/avw)**2,loka,0,il-1,fi,gi0)
          gi(il) = gi0(il)/(radius/avw)**il
        enddo
      endif
      do  ila = 1, nlma
        call info(10,1,0,' Compare numerical YL expns and radial '//
     .    'deriv of H, 1C for ilm=%i, ikap=%i',ila,ikap)
        if (ixi .eq. 1 .and. irot .eq. 0) 
     .  write(stdo,928) 'H0',(0d0,ilb=1,ila-1),
     .    gi(1+ll(ila)),(0d0,ilb=ila+1,min(nlmy,nlmb))
        write(stdo,928) 'Ha',( vvaln(ilb,ila), ilb = 1,min(nlmy,nlmb))
        write(stdo,928) '1c',( vvalj(ilb,ila), ilb = 1,min(nlmy,nlmb))
        write(stdo,928) '- ',(vvaln(ilb,ila)-vvalj(ilb,ila),
     .                    ilb = 1, min(nlmy,nlmb))
        write(stdo,928) 'H''',(dvaln(ilb,ila), ilb = 1,min(nlmy,nlmb))
        write(stdo,928) '1''',(dvalj(ilb,ila), ilb = 1,min(nlmy,nlmb))
        write(stdo,928) '- ',(dvaln(ilb,ila)-dvalj(ilb,ila),
     .                    ilb = 1, min(nlmy,nlmb))
      enddo
  928 format(a,1x,25F9.5)

      call info0(10,1,0,' Summary')
      call prterr(nlmy,nlmb,nlma,vvalj,vvaln,
     .  err,ier1,ier2,xx,ixx,ixx)
      call info5(10,0,0,' Ha vs 1c:    maximum difference %g '//
     .  'at ilm = %i, ilm'' = %i',err,ier2,ier1,0,0)
      call prterr(nlmy,nlmb,nlma,dvalj,dvaln,
     .  err,ier1,ier2,xx,ixx,ixx)
      call info5(10,0,2,' Ha'' vs 1c'':  maximum difference %g '//
     .  'at ilm = %i, ilm'' = %i',err,ier2,ier1,0,0)


      call getqu(ltmp)
      if (.not. ltmp) call fexit(0,0,outs,0)

c     goto 21

C --- Plot screened Hankel along a vector ---
   70 continue

      call info2(0,1,1,' ... plotting line from x =%3:1,d to x =%3:1,d',
     .  v1,v2)

   71 call query('number of points=',2,nvec)
      if (nvec .gt. nvmax) then
        print *, 'number cannot exceed max=',nvmax
        goto 71
      endif
      iat = 1
   73 call query('atom=',2,iat)
      if (iat .gt. nbasp) then
        print *, 'atom cannot be larger than nbas=',nbasp
        goto 73
      endif
      ixi = 1
      offR  = ival(w(ontab),iat)
      nlmb = ival(w(oiax),niax*offR+9)
C   72 ilb = 1
C      call query('orbital L to plot=',2,ilb)
C      if (ilb .gt. nlmb) then
C        print *, 'L cannot be larger than nlmb=',nlmb
C        goto 72
C      endif
   74 ikap = 1
      if (nkap .gt. 1) call query('ikap=',2,ikap)
      if (ikap .gt. nkap) then
        print *, 'atom cannot be larger than nkap=',nkap
        goto 74
      endif
      lb = ll(nlmb)
      is = ival(w(oips),iat)
      call upack('spec hcr',sspec,is,hcr,0,0,0)
      call info5(10,1,0,' Atom %i has hcr=%n:1;4,4d',iat,lb+1,hcr,0,0)

      do  i = 1, nvec
        xmrp(1) = v1(1) + dble(i-1)/max(nvec-1,1)*(v2(1)-v1(1))
        xmrp(2) = v1(2) + dble(i-1)/max(nvec-1,1)*(v2(2)-v1(2))
        xmrp(3) = v1(3) + dble(i-1)/max(nvec-1,1)*(v2(3)-v1(3))
        call dscal(3,1/avw,xmrp,1)
        call strck(nlf**2,iat,ixi,w(ortab),ekap,ikap,nkap,nbasp,
     .    w(oalpha),w(ohalp),w(oiax),w(ontab),w(os),w(ob),nitab,w(ocy),
     .    xmrp,nlf-1,irot+loka,slk,slj)
        reslin(i,1) = avw*xmrp(1)
        reslin(i,2) = avw*xmrp(2)
        reslin(i,3) = avw*xmrp(3)
        do  ilb = 1, min(nlmb,9)
          reslin(i,3+ilb) = slk(ilb,1)
        enddo
      enddo

      call info0(10,0,0,' writing line data to file "plot"')
      ifi = fopn('PLOT')
      rewind ifi
      call ywrm(0,' ',1,ifi,'(%6,6g)',reslin,0,nvmax,nvec,
     .  min(3+nlmb,12))
      call fclose(ifi)
C     call fexit(0,0,outs,0)

c ... write line data to file 76 for gnuplot
c     do i = 1, nvec
c       write(76,'(1x,3f12.6,1x,18g14.5)')
c    .    (reslin(i,ilb),ilb=1, 3+min(nlmb,9))
c     enddo

      end

