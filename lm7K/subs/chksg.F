      subroutine chksg(sctrl,slat,sspec,sstr,sarray)
C- Plots or checks the value-Laplacian envelope functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nl
Ci     Stored:
Ci     Passed to:
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat ocy
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
Cr       Invoke with --pltg:con
Cr   2.  tabulate screened envelope on a line
Cr       Invoke with --pltg:line[:v1x,v1y,v1z,v2x,v2y,v2z]
Cr   3.  check one-center expansion
Cr       Invoke with --pltg:onec
Cb Bugs
Cu Updates
Cu   23 Jul 08 (S. Lozovoi) Adapted from chkstr.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision sctrl(1),slat(1),sspec(1),sstr(1),sarray(1)
C ... Local parameters
      double precision plat(3,3),alat,ckbas,cksumf,tstmax
      integer npmax,nlmy,nvmax,nkap
      parameter (npmax=25, nlmy=16, nvmax=1001, nkap=2)
      double precision slj(2,25),slk(2,25),
     .  vvalk(nlmy,nlmy,nkap),vvalj(nlmy,nlmy,nkap),
     .  v1ck(nlmy,nlmy,nkap),v1cj(nlmy,nlmy,nkap) 
      double precision reslin(nvmax,21)
      double precision resk(-npmax:npmax,-npmax:npmax),xx
      equivalence (resk,reslin)
      double precision xp(3,122),xpt(122,3),wp(122),yl(122,nlmy),
     .  rsq(122),v1(6),v2(3),radius,pix(3),
     .  xmrp(3),fi(10),hcr(10)
      equivalence (v1(4),v2)
      integer obas,ocy,offR,offRp,oiax,ontab,oips,olmx,ohcr
      integer a2vec,fopn,i,iat,iatp,ifi,ikap,ilm,imax,ip,ival,ixi,
     .  j,j1,j2,ll,ila,ilb,mode,nbas,niax,nl,nlf,np,
     .  nttab,nvec,lio,nbasf,stdo,nglob,
     .  is,nclus,nlma,nlmb,lb,nspec,ivl,ivlf,il
      logical ltmp,iosg,iosg1,cmdopt
      character outs*120,out2*120,dc*1,dc2*1
      parameter (niax=10)
      integer nttabg,ii,ier1,ier2,ixx
      integer okmx,orsma,ontabg,oehvl,osg,orsmh
      integer ocg,ojcg,oidxcg         
      double precision err
C     parameter (tiny=1d-8)

C ... Heap
      integer w(1)
      common /w/ w

      stdo = nglob('stdo')
C ... Clebsch-Gordan and related coefficents
      call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)

      call upack('ctrl nbas nspec nl',sctrl,nbas,nspec,nl,0,0)
      call upack('lat alat plat',slat,alat,plat,0,0,0)
      call upack('array opos oips',sarray,obas,oips,0,0,0)
      call upack1('lat ocy',slat,ocy)
      call upack1('str ivl',sstr,ivl)
      call isanrg(ivl,0,2,'ASASTR:','ivl mode (IVFUN)',.true.)

      call sp2cls('spec lmxb',sspec,0,1,1,nspec,olmx)
      call sp2cls('spec hcr',sspec,0,nl,1,nspec,ohcr)
      call sp2cls('spec orbp',sspec,0,nl,1,nspec,orsmh)
      call sp2cls('spec rsma',sspec,0,1,1,nspec,orsma)
      call sp2cls('spec kmxt',sspec,0,1,1,nspec,okmx)
      if (ivl .ne. 0)
     .  call sp2cls('spec ehvl',sspec,0,nl,1,nspec,oehvl)


C ... Read some strux data from from disk
      call info0(10,1,0,
     .  ' CHKSG: check or plot strx from file STRG ...')
      lio = 0
      ltmp = iosg1(lio,'STRG',ivlf,nttab,nlf,nbasf,nttabg,xx,
     .  ifi)

      call isanrg(ivlf,ivl,ivl,'file:CHKSG:','ivl',.true.)
      call isanrg(nbasf,nbas,nbas,'file:CHKSG:','nbas',.true.)
      ckbas = cksumf(w(obas),3*nbas)
      ltmp = iosg(lio+8,'STRG',ivl,nl,nbas,ckbas,
     .             nttabg,oiax,ontab,ontabg,osg)
      call isanrg(nlf,nl,nl,'file:CHKSG:','nl',.true.)
      call isanrg(w(ontab+nbas),nttab,nttab,
     .  'file:CHKSG:','nttab',.true.)
      call info5(10,0,0,
     .  ' %i pairs in neighbor table  out of %i total.'//
     .  ' nl = %i, nbas = %i',
     .  nttabg,nttab,nl,nbas,0)
C     nds = nl**2

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
      nvec = 1001

C ... get plotting(/testing) mode
      if (cmdopt('--pltg',6,0,outs)) then
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
            call rxs2('chksg: failed to parse --pltg option: "',
     .        out2(1:j2),'%a ..."')
          endif
        endif
      endif

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
c     iat = 4
   63 call query('atom=',2,iat)
      if (iat .gt. nbas) then
        print *, 'atom cannot be larger than nbas=',nbas
        goto 63
      endif
   62 continue 
c     ilb = 1
c     ilb = 3
      ilb = 5
      call query('orbital L to plot=',2,ilb)
      offR  = ival(w(ontab),iat)
      nlmb = ival(w(oiax),niax*offR+9)
      if (ilb .gt. nlmb) then
        print *, 'L cannot be larger than nlmb=',nlmb
        goto 62
      endif
   64 continue 
      ikap = 1
c     ikap = 2
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

C ... obtain f(x,y) on a regular mesh 
      call acoord(iat,1,alat,plat,w(obas),w(oiax),w(ontab),pix)
      do  15  i = -imax, imax
        do  16  j = -imax, imax
          xmrp(1) =  (i*tstmax)/imax + pix(1)
          xmrp(2) =  (j*tstmax)/imax + pix(2)
          xmrp(3) =  1d-4 + pix(3)
          call strgck(ivl,nl**2,nl,iat,w(oips),w(orsmh),w(oehvl),
     .      alat,plat,w(obas),w(oiax),w(ontab),w(ontabg),w(ocy),
     .      xmrp,w(osg),slj,slk)
          resk(i,j) = max(min(slj(ikap,ilb),999d0),-99d0)
c         resk(i,j) = max(min(slk(ikap,ilb),999d0),-99d0)
   16   continue
   15 continue
      ifi = fopn('PLOTG')
      rewind ifi
      call ywrm(0,' ',1,ifi,'(%5,6g)',resk,0,2*imax+1,2*imax+1,2*imax+1)
      call fclose(ifi)

      return

C--------------------------------------------------------------------------------
C --- One-center expansion ---
   20 continue
C --- Test expansion for R' corresponding to iax(ixi) ---
      call info0(10,1,0,
     .  ' Compare numerically integrated YL expansions of'//
     .  ' the value-Laplacian functions')
      call info0(10,0,0,  
     .  ' and their 1-center expansion'//
     .  ' on a neighboring sphere surface.')
      call info0(10,1,0,  
     .  ' These should be equal'//
     .  ' apart from numerical integration errors.')
      call info0(10,1,0,  
     .  ' When evaluated at radius=hcr,'//
     .  ' function should be 1 on head sphere (0 for ikap=2)')
      call info0(10,0,1,  
     .  ' and vanish on tail spheres.')

      iat = 1
c     iat = 4
   23 continue                     
      call query('atom=',2,iat)
      if (iat .gt. nbas) then
        print *, 'atom cannot be larger than nbas=',nbas
        goto 23
      endif
   25 continue 
c     ixi = 5
      ixi = 1
      call query('ixi= (1 for head, 2 for 1st NN, etc)',2,ixi)
      offR  = ival(w(ontab),iat)
      nclus = ival(w(ontabg),iat) 
      nlma = ival(w(oiax),niax*offR+9)
      if (ixi .gt. nclus) then
        print *, 'ixi cannot be larger than nclus=',nclus
        goto 25
      endif
      offRp = ixi-1 + offR
      iatp = ival(w(oiax),niax*offRp+2)
      is = ival(w(oips),iatp)
      nlmb = ival(w(oiax),niax*offRp+9)
      lb = ll(nlmb)
c la, nlma - head of the cluster (defines number of U functions)
c lb, nlmb - expansion sphere (defines max L to check)

      call info5(10,0,0,' Neighbor ixi=%i corresponds to atom %i',
     .  ixi,iatp,0,0,0)
      call upack('spec hcr',sspec,is,hcr,0,0,0)
      call info5(10,0,0,' Atom %i has nlmb=%i hcr=%n:1;4,4d',
     .  iatp,nlmb,lb+1,hcr,0)

c Uncomment next 4 lines and comment line 'radius = hcr(il+1)' further down 
c if one wishes to integrate over fixed sphere radius introduced interactively
c  21 continue
c     call query('radius=',4,radius)
c     call info5(10,1,0,' ... using neighbor ixi=%i (jb=%i), radius=%d' 
c    .  ,ixi,iatp,radius,0,0)
      call info5(10,1,0,
     .  ' ... using neighbor ixi=%i (jb=%i), radius=%n:1;4,4d', 
     .  ixi,iatp,lb+1,hcr,0)

C --- One-center expansion by brute force integration ---
C ... Integration takes place on sphere R' = iax(2,ixi)
      call fpiint(-122,0,np,xp,wp)
C ... Normalized spherical harmonics for all points on sphere R'
      call dmcpy(xp,1,3,xpt,np,1,np,3)
      call ropyln(np,xpt(1,1),xpt(1,2),xpt(1,3),lb,np,yl,rsq)
      do ip = 1, np
        do  ilb = 1, nlmb
          yl(ip,ilb) = yl(ip,ilb)*wp(ip)
        enddo
      enddo
      call dpzero(vvalj,nlmy*nlmy*nkap)
      call dpzero(vvalk,nlmy*nlmy*nkap)
C ... Integrate U on sphere R'
      call acoord(iat,ixi,alat,plat,w(obas),w(oiax),w(ontab),pix)
      do il = 0, lb
        radius = hcr(il+1)
        do  ip = 1, np
          do ii = 1, 3
            xmrp(ii) =  radius*xp(ii,ip) + pix(ii)
          enddo 
            call strgck(ivl,nl**2,nl,iat,w(oips),w(orsmh),w(oehvl),
     .        alat,plat,w(obas),w(oiax),w(ontab),w(ontabg),w(ocy),
     .        xmrp,w(osg),slj,slk)
C       Integration over the surface of sphere a_R' means that
C       U^a_RL(x on sphere R') = sum_L' vvalj_L'(a_R') Y_L'(a_R')
C       vvalj holds U^a(x-R'); vvalk holds \lap U^a(x-R')
          do  ilb = il*il+1, (il+1)**2
c           do  ila = 1, nl**2
            do  ila = 1, nlma
              do  ikap = 1, nkap
                vvalj(ilb,ila,ikap) = vvalj(ilb,ila,ikap) + 
     .                         slj(ikap,ila)*yl(ip,ilb)
                vvalk(ilb,ila,ikap) = vvalk(ilb,ila,ikap) + 
     .                         slk(ikap,ila)*yl(ip,ilb)
              enddo
            enddo
          enddo
        enddo
      enddo

C --- Get the one-center expansion through the P_kL route ---
      call strg1c(ivl,nl**2,nlmy,nl,iat,ixi,w(oips),w(orsmh),w(oehvl),
     .  w(ohcr),w(okmx),w(orsma),alat,plat,w(obas),w(oiax),
     .  w(ontab),w(ontabg),w(ocy),w(ocg),w(oidxcg),w(ojcg),w(osg),
     .  v1cj,v1ck)

      do  ikap = 1, nkap
        do  ila = 1, nlma
          call info(10,1,0,' Compare numerical YL (num) and '//
     .      '1C expansion (1c) for ilm=%i, ikap=%i',ila,ikap)
          write(stdo,928) 
     .      'val(num)',(vvalj(ilm,ila,ikap), ilm = 1,min(nlmy,nlmb))
          write(stdo,928) 
     .      'val(1c) ',( v1cj(ilm,ila,ikap), ilm = 1,min(nlmy,nlmb))
          write(stdo,928) 
     .      '   -    ',(vvalj(ilm,ila,ikap)-v1cj(ilm,ila,ikap),
     .                      ilm = 1, min(nlmy,nlmb))

          write(stdo,928) 
     .    'Lap(num)',(vvalk(ilm,ila,ikap), ilm = 1,min(nlmy,nlmb))
          write(stdo,928) 
     .    'Lap(1c) ',( v1ck(ilm,ila,ikap), ilm = 1,min(nlmy,nlmb))
          write(stdo,928) 
     .    '   -    ',(vvalk(ilm,ila,ikap)-v1ck(ilm,ila,ikap),
     .                    ilm = 1, min(nlmy,nlmb))
        enddo
      enddo
  928 format(a,1x,25f9.5)
c 930 format(a,1x,25d9.1)

      do  ikap = 1, nkap
        call info2(10,1,0,' Summary for ikap = %i',ikap,0)
        call prterr(nlmy,nlf**2,nlf**2,vvalj(1,1,ikap),v1cj(1,1,ikap),
     .    err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,0,' value:     maximum difference %g '//
     .    'at ilm = %i, ilm'' = %i',err,ier2,ier1,0,0)
        call prterr(nlmy,nlf**2,nlf**2,vvalk(1,1,ikap),v1ck(1,1,ikap),
     .    err,ier1,ier2,xx,ixx,ixx)
        call info5(10,0,1,' Laplacian: maximum difference %g '//
     .    'at ilm = %i, ilm'' = %i',err,ier2,ier1,0,0)
      enddo

      call getqu(ltmp)
      if (.not. ltmp) call fexit(0,110,outs,0)

      goto 23

C --- Plot value-Laplacian functions along a vector ---
   70 continue

      call info2(0,1,1,' ... plotting line from x =%3:1,d to x =%3:1,d',
     .  v1,v2)

   71 continue
      call query('number of points=',2,nvec)
      if (nvec .gt. nvmax) then
        print *, 'number cannot exceed max=',nvmax
        goto 71
      endif
      iat = 1
   73 continue
      call query('atom=',2,iat)
      if (iat .gt. nbas) then
        print *, 'atom cannot be larger than nbas=',nbas
        goto 73
      endif

      offR  = ival(w(ontab),iat)
      nlmb = ival(w(oiax),niax*offR+9)
      lb = ll(nlmb)
      is = ival(w(oips),iat)
      call upack('spec hcr',sspec,is,hcr,0,0,0)
      call info5(10,1,0,' Atom %i has hcr=%n:1;4,4d',iat,lb+1,hcr,0,0)

      call acoord(iat,1,alat,plat,w(obas),w(oiax),w(ontab),pix)
      do  i = 1, nvec
        do ii = 1, 3
          xmrp(ii) = pix(ii) + v1(ii) + 
     .      dble(i-1)/max(nvec-1,1)*(v2(ii)-v1(ii)) 
        enddo

        call strgck(ivl,nl**2,nl,iat,w(oips),w(orsmh),w(oehvl),
     .    alat,plat,w(obas),w(oiax),w(ontab),w(ontabg),w(ocy),
     .    xmrp,w(osg),slj,slk)

        reslin(i,1) = xmrp(1)
        reslin(i,2) = xmrp(2)
        reslin(i,3) = xmrp(3)
        do  ilb = 1, min(nlmb,9)
          reslin(i,3+ilb) = slj(1,ilb)
          reslin(i,3+min(nlmb,9)+ilb) = slj(2,ilb)
        enddo
c ... fix the printout bug
        do ilb = 4, 2*min(nlmb,9)+3
          if (dabs(reslin(i,ilb)) .lt. 1d-99) reslin(i,ilb)=1d-99
        enddo
      enddo

      call info0(10,0,0,' writing line data to file "plotg"')
      ifi = fopn('PLOTG')
      rewind ifi
      call ywrm(0,' ',1,ifi,'(%6,6g)',reslin,0,nvmax,nvec,
     .  min(3+nlmb,21))
      call fclose(ifi)

c write line data to file fort.77 for gnuplot for both kappa
      do i = 1, nvec
        write(77,'(1x,3f12.6,1x,18g14.5)')
     .    (reslin(i,ilb),ilb=1, 3+2*min(nlmb,9))
      enddo

      end

