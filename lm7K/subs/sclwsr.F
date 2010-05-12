Cdelw1 var ['ontab', '0', 'ontab', 'nbasp+1', 'integer', 'iv_w_', 'sclwsr']
Cdelw1 var ['oiax', '0', 'oiax', '-niax*mxnbr', 'integer', 'iv_w_', 'sclwsr']
Cdelw1 var ['owk', '0', 'owk', '3*mxnbr', 'real(8)', 'rv_w_', 'sclwsr']
      subroutine sclwsr(opt,nbas,nbasp,nspec,alat,plat,basp,ips,modep,
     .  slabl,z,lock,volfac,wsmax,omax1,omax2,wsr)

C- Scales sphere radii to specified volume, satisfying constraints.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   opt  :specifies options in sclwsr.
Ci         :1s digit specifies whether scaling reaches target volfac
Ci         :         is required:
Ci         : 0 meeting volume target is optional
Ci         : 1 meeting volume target is optional but
Ci         :   a warning is printed if not met.
Ci         : 2 meeting volume target is required
Ci         :10s digit concerns treatment of empty spheres
Ci         : 0 ES and other sites are treated symmetrically
Ci         : 1 all sites with z>0 are resized first; then
Ci             all sites are resized.
Ci         : 2 all sites with z>0 are resized first; then
Ci             the ES sites only are resized.
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci   nspec :number of atom species
Ci   alat  :length scale
Ci   plat  :primitive lattice vectors, in units of alat
Ci   basp  :basis vectors (scaled by alat; padded for layer mode)
Ci   ips   :the jth atom belongs to species ips(j)
Ci   modep :specifies which dimensions have periodic b.c.
Ci         :In particular,
Ci         :0 suppresses shifts along plat(j)
Ci         :2 shifts to minimize length of pos
Ci          are used
Ci   slabl :species labels
Ci   z     :nuclear charge, by species
Ci   lock  :constraints specifying which species are locked and
Ci         :which are free to float.
Ci         :On input, lock should be zero or two for each species.
Ci         :Each species for which lock(i)=2 is constrained not to change.
Ci         :The radii for other species are floated.
Ci         :lock(1..nspec) is OVERWRITTEN on output
Ci   volfac:scale until sum of sphere vol = volfac * cell vol or until
Ci         :all sphere radii are constrained (see Remarks, sclws2)
Ci   wsmax :(wsmax>0) a global maximum on the size of MT spheres.
Ci         :No sphere is allowed to exceed wsmax.
Ci         :(wsmax=0) no constraint is imposed.
Ci   omax1 :max. allowed overlap divided by distance (r1+r2-d)/d<omax1
Ci         :omax1(1) is constraint for atom-atom overlaps
Ci         :omax1(2) is constraint for atom-empty-sphere overlaps
Ci         :omax1(3) is constraint for empty-sphere-empty-sphere overlaps
Ci   omax2 :max. allowed overlap divided by radius  (r1+r2-d)/r1<omax2
Ci         :omax2(1) is constraint for atom-atom overlaps
Ci         :omax2(2) is constraint for atom-empty-sphere overlaps
Ci         :omax2(3) is constraint for empty-sphere-empty-sphere overlaps
Co Inputs/Outputs:
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cio        :On input, starting values for wsr
Cio        :On output, final values for wsr
Cl Local:
Cl   dovl1 :maximum overlap divided by distance
Cl   dovl2 :maximum overlap divided by radius
Cl   gamma :(zero) passed to sclws2
Cl   nrspec:number of atoms in each species
Cr Remarks
Cr  Sphere radii are scaled in an iterative procedure.  In any
Cr  iteration, species are divided into those that are 'locked'
Cr  (frozen) and those that are allowed to float.  The largest scaling
Cr  factor is determined for all those species of the latter type, that
Cr  satisfies the constraints (see below).  These species are scaled
Cr  and a new iteration begins.  By construction, each iteration will
Cr  cause at least one new species to be locked; thus, the total number
Cr  of iterations will not exceed the number of species.
Cu Updates
Cu   17 Jan 09  bug fix 10s digit mode=2, no ES
Cu   17 May 02  New 10s digit opt switch gives optional lower priority
Cu              to empty spheres
Cu   22 Apr 02  First created; adapted from Stuttgart LMTO56.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer opt,nbas,nbasp,nspec,ips(nbas),modep(3),lock(nspec)
      double precision alat,basp(3,nbasp),volfac,wsmax,
     .                 omax1(3),omax2(3),plat(3,3),wsr(nspec),z(nspec)
      character*8 slabl(*)
C Local variables:
Cdelw1       integer i,k,owk,nrspec(nspec),ib,is,opt1
      integer:: i , k , nrspec(nspec) , ib , is , opt1
      real(8) ,allocatable :: rv_w_owk(:)

      integer niax,mxcsiz,mxnbr,nttab,ipr,stdo,lgunit,llock(nspec)
Cdelw1       integer ontab,oiax
      integer ,allocatable :: iv_w_ontab(:)
      integer ,allocatable :: iv_w_oiax(:)

      logical les
      double precision dovl1(3),dovl2(3),gamma,range(nbas),wsrs(nspec)
      double precision facr,tiny,avw,vol,avwsr,volnew,volnes,volold,
     .  volsph
      parameter(niax=10,mxcsiz=200,facr=2d0,tiny=1d-5)
C ... Heap
      integer w(1)
      common /w/ w

      integer:: istdo,fopna

C     omax1(1) = -.01
C      wsr(1) = 3.232247d0
C      wsr(2) = 3.232247d0
C      wsr(3) = 2.248243d0
C      wsr(4) = 2.097284d0
C      wsr(5) = 1.647651d0
C      wsr(6) = 1.647651d0

C --- Setup ---
      call getpr(ipr)
      stdo  = lgunit(1)
      gamma = 0
      call maknrs(nbas,ips,ib,nrspec)
      if (ib .gt. nspec) call rx('sclwsr: wrong number of species')
      opt1 = mod(opt/10,10)
      call dcopy(nspec,wsr,1,wsrs,1)
      avw = avwsr(plat,alat,vol,nbas)

C --- Make a neighbor table and adjust llock to freeze ES sites ---
      do  ib = 1, nbas
        is = ips(ib)
        range(ib) = facr * wsr(is)
      enddo
      nbasp = nbas
      mxnbr = mxcsiz*nbas
Cdelw1       call defi(ontab, nbasp+1)
      allocate(iv_w_ontab(nbasp+1))
Ctakao_ZeroClear_NotRequiered        if (nbasp+1<0) iv_w_ontab(:)=0

Cdelw1       call defi(oiax, -niax*mxnbr)
      allocate(iv_w_oiax(abs(-niax*mxnbr)))
      if (-niax*mxnbr<0) iv_w_oiax(:)=0

Cdelw1       call defdr(owk,  3*mxnbr)
      allocate(rv_w_owk(3*mxnbr))
Ctakao_ZeroClear_NotRequiered        if (3*mxnbr<0) rv_w_owk(:)=0.0d0

      nttab = mxnbr
      call pshpr(ipr-20)
Cdelw1       call pairc(1,nbas,nbasp,modep,20,0,alat,plat,basp,basp,range,-1,1,
Cdelw1      .  nttab,w(ontab),w(oiax),w(owk),k)
      call pairc ( 1 , nbas , nbasp , modep , 20 , 0 , alat , plat 
     . , basp , basp , range , - 1 , 1 , nttab , iv_w_ontab , iv_w_oiax 
     . , rv_w_owk , k )

      call poppr

C --- Scale the sphere radii, freezing empty spheres ---
      if (opt1 .ge. 1) then

C   ... Local copy of lock and wsr, adjusting to freeze ES sites
        call icopy(nspec,lock,1,llock,1)
        les = .false.
        do  is  = 1, nspec
          if (z(is) .eq. 0) then
            les = .true.
            llock(is) = 2
            wsr(is) = 0
          endif
        enddo

        if (les) then

C     ... Scale wsr with z=0 sites locked at wsr=0
          volold = volsph(nspec,nrspec,wsr)/vol
          volnes = volfac
Cdelw1           call sclws2(nbas,nspec,alat,plat,basp,slabl,w(oiax),ips,
Cdelw1      .      w(ontab),z,nrspec,omax1,omax2,gamma,wsmax,llock,volnes,
Cdelw1      .      dovl1,dovl2,wsr)
          call sclws2 ( nbas , nspec , alat , plat , basp , slabl , iv_w_oiax 
     .     , ips , iv_w_ontab , z , nrspec , omax1 , omax2 , gamma , wsmax 
     .     , llock , volnes , dovl1 , dovl2 , wsr )

          if (ipr .ge. 30) then
            call awrit2(' SCLWSR:  initial sphere packing = %;1d%%'//
     .        ' scaled to %;1d%% (no empty spheres)',
     .        ' ',120,stdo,100*volold,100*volnes)
          endif

C     ... Restore wsr(Z=0)
          call icopy(nspec,lock,1,llock,1)
          do  is  = 1, nspec
            if (z(is) .eq. 0) wsr(is) = wsrs(is)
          enddo
        endif

C   ... 10s digit opt=2 : freeze wsr(Z>0)
        if (opt1 .ge. 2 .and. les) then
          do  is  = 1, nspec
            if (z(is) .ne. 0) lock(is) = 2
          enddo
        endif

      endif

C --- Scale the sphere radii (2nd pass for les) ---
      volold = volsph(nspec,nrspec,wsrs)/vol
      volnew = volfac
Cdelw1       call sclws2(nbas,nspec,alat,plat,basp,slabl,w(oiax),ips,w(ontab),
Cdelw1      .  z,nrspec,omax1,omax2,gamma,wsmax,lock,volnew,dovl1,dovl2,wsr)
      call sclws2 ( nbas , nspec , alat , plat , basp , slabl , iv_w_oiax 
     . , ips , iv_w_ontab , z , nrspec , omax1 , omax2 , gamma , wsmax 
     . , lock , volnew , dovl1 , dovl2 , wsr )


C --- Printout ---
Cki#error, have return with len(w_varlist)>0 at line 183
Cdelw1 w_varlist remains: ontab oiax owk
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: [ontab oiax owk]
Cdelw1       if (ipr .lt. 10) return
      if ( ipr .lt. 10 ) then
        if (allocated(rv_w_owk)) deallocate(rv_w_owk)
        if (allocated(iv_w_oiax)) deallocate(iv_w_oiax)
        if (allocated(iv_w_ontab)) deallocate(iv_w_ontab)
        return
      endif

C      if (ipr .ge. 10) write(stdo,309)
C     .  vol,100*volold,100*volnew
C  309 format(/' SCLWSR: vol=',f11.3,
C     .  '  sphere fraction=',f5.1,
C     .  '%(initial)  ',f5.1,'%(scaled)')

      if (ipr .ge. 10) then
        call awrit3(' SCLWSR:  vol = %,5g a.u..%3f'//
     .    'Initial sphere packing = %;1d%%  '//
     .    'scaled to %;1d%%',
     .    ' ',120,stdo,vol,100*volold,100*volnew)
      endif

      if (ipr .ge. 30) then
        write(stdo,310)
     .    (omax1(i)*100,i=1,3),(omax2(i)*100,i=1,3),
     .    (dovl1(i)*100,i=1,3),(dovl2(i)*100,i=1,3)
  310   format(1x,'constr omax1=',3f6.1,' %    omax2=',3f6.1,' %',
     .        /1x, 'actual omax1=',3f6.1,' %    omax2=',3f6.1,' %')
c        if (ipr .gt. 30) then
        write(stdo,311)
  311   format(/' spec  name',8x,'old rmax    new rmax     ratio')
        do  i = 1, nspec
          write(stdo,312) i,slabl(i), wsrs(i), wsr(i), wsr(i)/wsrs(i)
  312     format(i4,3x,a,3f12.6)
        enddo
c        endif
      endif

      istdo = fopna('rmt',-1,2)
      do  i = 1, nspec
        write(istdo,"(a,f12.6)") slabl(i), wsr(i)
      enddo
      call fclose(istdo)

      if (dabs(volnew-volfac) .gt. tiny) then
        if (mod(opt,10) .ge. 1) then
          write(stdo,321) int(volfac*100)
  321     format(/' SCLWSR (warning): failed to reach target vol (',
     .      i3,'% of cell vol)')
        endif
        if (mod(opt,10) .eq. 2) then
          call rx('SCLWSR: failed to reach target VOL.  Increase omax.')
        endif
      endif

Cdelw1 rlse name= ontab old_list= owk oiax ontab
Cdelw1 rlse name= ontab new_list= owk oiax
Cdelw1       call rlse(ontab)
      if (allocated(iv_w_ontab)) deallocate(iv_w_ontab)

Cgetarg       end
Cdelw1 w_varlist remains: owk oiax
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: [owk oiax]
      if (allocated(iv_w_oiax)) deallocate(iv_w_oiax)
      if (allocated(rv_w_owk)) deallocate(rv_w_owk)

      end subroutine sclwsr


      subroutine sclws2(nbas,nspec,alat,plat,bas,slabl,iax,ips,ntab,z,
     .  nrspec,omax1,omax2,gamma,wsmax,lock,volfac,dovl1,dovl2,wsr)

C- Enlarges the spheres to reach a specified volume under constraints
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   slabl :species labels (for printout)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   z     :nuclear charge, by species
Ci   nrspec:number of atoms in the ith species
Ci   omax1 :max. allowed overlap divided by distance (s1+s2-d)/d<omax1
Ci   omax2 :max. allowed overlap divided by radius  (s1+s2-d)/s1<omax2
Ci   gamma :a factor that changes scaling wsr from simple multiplicative
Ci         :scaling to a combination of additive + multiplicative scaling
Ci         :That is, in each iteration,
Ci         :scaling is r -> a(r+b) with a*b=gamma*(a-1)*avw
Ci         :gamma>0 tends enlarge small spheres faster than large ones
Ci   wsmax :(wsmax>0) a global maximum on the size of MT spheres.
Ci         :No sphere is allowed to exceed wsmax.
Ci         :(wsmax=0) no constraint is imposed.
Cio Inputs/Outputs:
Cio  lock  :constraints specifying which species are locked and
Cio        :which are free to float.
Cio        :On input, lock should be zero or two for each species.
Cio        :Each species for which lock(i)=2 is constrained not to change.
Cio        :The radii for other species are floated.
Cio        :lock(1..nspec) is OVERWRITTEN on output
Cio  volfac:scale until sum of sphere vol = volfac * cell vol
Cio        :or until all sphere radii are constrained (see Remarks)
Cio        :On input, volfac=target ratio (sum of sphere vol)/(cell vol)
Cio        :On output, volfac=actual ratio
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cio        :On input, starting values.
Cio        :On output, scaled values.
Co Outputs:
Co   dovl1 :maximum overlap divided by distance
Co   dovl2 :maximum overlap divided by radius
Cr Remarks:
Cr  Sphere radii are scaled in an iterative procedure.  In any
Cr  iteration, species are divided into those that are 'locked'
Cr  (frozen) and those that are allowed to float.  The largest scaling
Cr  factor is determined for all those species of the latter type, that
Cr  satisfies the constraints (see below).  These species are scaled
Cr  and a new iteration begins.  By construction, each iteration will
Cr  cause at least one new species to be locked; thus, the total number
Cr  of iterations will not exceed the number of species.
Cr
Cr  Typically sclws2 is called by is a higher-level routine that creates
Cr  the neighbor-table and other necessary arrays.
Cr
Cr  This code was adapted from Stuttgart routine blowup, v LMTO56.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer niax
      parameter (niax=10)
      integer nbas,nspec,iax(niax,*),ntab(*),ips(nbas),nrspec(nspec)
      integer lock(nspec)
      double precision alat,bas(3,nbas),dovl1(3),dovl2(3),wsmax,
     .  volfac,gamma,omax1(3),omax2(3),plat(3,3),wsr(nspec),z(nspec)
      character*8 slabl(nspec)
C Local variables:
      integer ib,jb,is,iclbsj,ip,k1,k2,k3,
     .        kb,ks,kco,kpr,lgunit,nloop,stdo,ipr,npcol,locki
      double precision a,amax,amax1,amax2,amax3,amax4,avw,b,bmax,d,
     .                 dm(0:3),dovlap,dscl(nspec),
     .                 dsclmx,drr2,dr(3),fpi3,gw,opo1,omo2,p,q,r,ratio,
     .                 rik,riko,s,t,tiny,u,v,vol,vola,volb,
     .                 volsph,wsri,wsrk,x,avwsr,Vconst
      logical fin
      character*72 fmt
      parameter(fpi3=4.18879020478639053d0,tiny=1d-5)
C ... External calls
      external dcopy,dpzero,dscal,errmsg,getpr,iinit

      avw  = avwsr(plat,alat,vol,nbas)
      gw   = avw*gamma
      stdo = lgunit(1)
      call getpr(ipr)
C      print *, '!!'
C      ipr = 40
      npcol = 7

C     call dcopy(nspec,1d0,0,dscl(1,0),1)

      fmt = '(''   SPEC:    '',6(1x,a8):/(12x,6(1x,a8)))'
      write(fmt(17:17),'(I1)') npcol
      write(fmt(32:32),'(I1)') npcol
      if (ipr .ge. 50) then
        write(stdo,'('' '')')
        write(stdo,fmt) (slabl(is),is=1,nspec)
      endif
      fmt = '(1x,''init rmt'',1x,6f9.5:/(10x,6f9.5))'
      write(fmt(19:19),'(I1)') npcol
      write(fmt(31:31),'(I1)') npcol
C     print *, fmt
      if (ipr .ge. 50) write(stdo,fmt) (wsr(is),is=1,nspec)

      do  nloop = 1, nspec+1
        amax = 9d9
        do  is = 1, nspec
          if (lock(is) .ne. 2) lock(is) = 0
        enddo

C   --- Lock radii of those spheres with maximum allowed overlap ---
C       and unlock those with radii exceeding maximum allowed.
        do  is = 1, nspec
          wsri = wsr(is)

C     ... If overlap criterion exactly satisfied for any connecting vector,
C         set lock for this species
          do  jb = 1, nrspec(is)
            ib = iclbsj(is,ips,nbas,jb)
            do  kpr = ntab(ib)+2, ntab(ib+1)
              kb = iax(2,kpr)
              k1 = iax(3,kpr)
              k2 = iax(4,kpr)
              k3 = iax(5,kpr)
              ks = ips(kb)
C           ip selects which ommax to use (A-A, A-E, or E-E)
              ip = 2
              if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
              if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
              opo1 = 1+omax1(ip)
              omo2 = 1-omax2(ip)
              wsrk = wsr(ks)
              rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
              rik = rik*alat
C           Set lock if any of these conditions are met:
C             wi + wk - rik = o1*rik
C             wi + wk - rik = o2*wi
C             wi + wk - rik = o2*wk
C             wi = wsmax
              locki = lock(is)
              if (dabs(opo1*rik-wsri-wsrk) .lt. tiny) locki=1
              if (dabs(rik-omo2*wsri-wsrk) .lt. tiny) locki=1
              if (dabs(rik-wsri-omo2*wsrk) .lt. tiny) locki=1
              if (wsmax.gt.0 .and. dabs(wsri-wsmax).lt.tiny) locki=1
              lock(is) = max(locki,lock(is))
            enddo
          enddo
C         print *,'is,lock0=',is,lock(is)

C     ... If overlap criterion exceeded for any connecting vector,
C         unset lock for this species
          if (lock(is) .ne. 2) then
            do  jb = 1, nrspec(is)
              ib = iclbsj(is,ips,nbas,jb)
              do  kpr = ntab(ib)+2, ntab(ib+1)
                kb = iax(2,kpr)
                k1 = iax(3,kpr)
                k2 = iax(4,kpr)
                k3 = iax(5,kpr)
                ks = ips(kb)
                wsrk = wsr(ks)
                rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                rik = rik*alat
C           ip selects which ommax to use (A-A, A-E, or E-E)
                ip = 2
                if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
                if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
                opo1 = 1+omax1(ip)
                omo2 = 1-omax2(ip)
C           Unset lock if any of these conditions are exceeded
C             wi + wk - rik < o1*rik
C             wi + wk - rik < o2*wi
C             wi + wk - rik < o2*wk
C             wi > wsmax
                if (opo1*rik-wsri-wsrk .lt. -tiny) lock(is)=0
                if (rik-omo2*wsri-wsrk .lt. -tiny) lock(is)=0
                if (rik-wsri-omo2*wsrk .lt. -tiny) lock(is)=0
                if (wsmax.gt.0 .and. wsri-wsmax.gt.tiny) lock(is)=0
              enddo
            enddo
          endif
        enddo
C       print *, 'lock', (lock(is), is=1,nspec)

C  --- Find amax=largest allowed scaling for unlocked species ---
        do  is = 1, nspec
          if (lock(is) .eq. 0) then
            riko = -1
            kco = -1
            wsri = wsr(is)
            do  jb = 1, nrspec(is)
              ib = iclbsj(is,ips,nbas,jb)
              do  kpr = ntab(ib)+2, ntab(ib+1)
                kb = iax(2,kpr)
                k1 = iax(3,kpr)
                k2 = iax(4,kpr)
                k3 = iax(5,kpr)
                ks = ips(kb)
                rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                rik = rik*alat
                ip = 2
                if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
                if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
                opo1 = 1+omax1(ip)
                omo2 = 1-omax2(ip)
                if (dabs(rik-riko).gt.tiny .or. kco.ne.ks) then
                  wsrk = wsr(ks)
                  riko = rik
                  kco = ks
                  amax1 = 9d9
                  amax2 = 9d9
                  amax3 = 9d9
                  amax4 = 9d9
C               If second site ks is locked, can only scale site i
                  if (lock(ks) .ne. 0) then
                    amax1 = (opo1*rik-wsrk+gw)/(wsri+gw)
                    if (omo2 .gt. 0d0)
     .            amax2 = (rik-wsrk+omo2*gw)/(omo2*wsri+omo2*gw)
                    amax3 = (rik-omo2*wsrk+gw)/(wsri+gw)
C               If neither site is locked, both will scale
                  else
                    amax1 = (opo1*rik+gw+gw)/(wsri+wsrk+gw+gw)
                    if (wsrk+omo2*wsri+(1+omo2)*gw .gt. 0d0)
     .            amax2 = (rik+(1+omo2)*gw)/
     .                  (wsrk+omo2*wsri+(1+omo2)*gw)
                    if (wsri+omo2*wsrk+(1+omo2)*gw .gt. 0d0)
     .            amax3 = (rik+(1+omo2)*gw)/
     .                  (wsri+omo2*wsrk+(1+omo2)*gw)
                  endif
                  if (wsmax .gt. 0) then
                    amax4 = wsmax/wsri
                  endif
                  amax = dmin1(amax,amax1,amax2,amax3,amax4)
C               print *, 'kc,kbas,kpr,amax',ks,kb,kpr,amax
                endif
              enddo
            enddo
          endif
        enddo
        bmax = (1d0-1/amax)*gw

C   ... Determine what new volume will be after scaling with a,b
        vola = 0d0
        volb = 0d0
        do  ib = 1, nbas
          is = ips(ib)
          if (lock(is) .eq. 0) then
            volb = volb + (amax*wsr(is)+bmax)**3
          else
            vola = vola + wsr(is)**3
          endif
        enddo
        vola = vola * fpi3
        volb = volb * fpi3

C   --- Case scaling will lead to new volume > final volume ---
        if (vol*volfac .lt. vola+volb) then
          call dpzero(dm,4)
          do  ib = 1, nbas
            is = ips(ib)
            if (lock(is) .eq. 0) then
              a = wsr(is)+gw
C             For numerical reasons distinguish cases
              if (dabs(gamma) .gt. 1d0) then
                b = wsr(is)
              else
                b = -gw
              endif
C             fpi3*dm(0) = sum of sphere volumes not constrained
              dm(0) = dm(0) +       b*b*b
              dm(1) = dm(1) + 3d0 * a*b*b
              dm(2) = dm(2) + 3d0 * a*a*b
              dm(3) = dm(3) +       a*a*a
            endif
          enddo
C         Vconst + a**3*Vuncst = Vtarget; Vconst + Vuncst = vola
C         If Vconst > Vtarget, constraint cannot be satisified
          Vconst = vola - fpi3*dm(3)
          if (Vconst .gt. vol*volfac) then
            call fexit2(-1,111,' Exit -1 : SCLWSR: constrained '//
     .        'sphere vol (%;0d) exceeds target vol=%;0d',
     .        Vconst,vol*volfac)
          endif
          if (dabs(dm(3)) .gt. tiny) then
            r =  dm(2) / dm(3)
            s =  dm(1) / dm(3)
            t = (dm(0) - (vol*volfac-vola)/fpi3) / dm(3)
            p = s - r*r/3d0
            q = 2d0*r*r*r/27d0 - r*s/3d0 + t
            d = p*p*p/27d0 + q*q/4d0
            u = (dsqrt(d)-q/2d0)**(1d0/3d0)
            v = -p/u/3d0
            x = u+v-r/3d0
            if (dabs(gamma) .gt. 1d0) then
              amax = x+1d0
              bmax = x*gw/amax
            else
              amax = x
              bmax = (1d0-1d0/amax)*gw
            endif
            if (ipr .ge. 100) then
              write(stdo,300)'R S T',r,s,t
              write(stdo,300)'P Q  ',p,q
              write(stdo,300)'  D  ',d
              write(stdo,300)' U V ',u,v
              write(stdo,300)' AMAX',amax
              write(stdo,300)' BMAX',bmax
              write(stdo,300)' -------------------------'
            endif
          endif
        endif

C   --- Scale unlocked spheres by w <- a*w+b ---
        fin = .true.
        do  is = 1, nspec
          dscl(is) = 1d0
          if (lock(is) .eq. 0) then
            dsclmx  = amax + bmax/wsr(is)
            wsr(is) = dsclmx*wsr(is)
            dscl(is)= dsclmx
            fin = .false.
          endif
        enddo

C       ratio = new volume / old volume
        ratio = volsph(nspec,nrspec,wsr)/vol
C       fin=T if volume change is small or all spheres are locked
        fin = fin .or. dabs(ratio-volfac).lt.tiny .or. nloop.eq.nspec+1

        fmt = '(1x,''iter:'',i3,1x,6f9.5:/(10x,6f9.5))'
C       print *, fmt
        write(fmt(19:19),'(I1)') npcol
        write(fmt(31:31),'(I1)') npcol
C       print *, fmt
        if (ipr .ge. 50) write(stdo,fmt) nloop,(dscl(is),is=1,nspec)

C   ... Last iteration
        if (fin) then
          call dpzero(dovl1,3)
          call dpzero(dovl2,3)
          do  is = 1, nspec
            wsri = wsr(is)
            do  jb = 1, nrspec(is)
              ib = iclbsj(is,ips,nbas,jb)
              do  kpr = ntab(ib)+2, ntab(ib+1)
                kb = iax(2,kpr)
                k1 = iax(3,kpr)
                k2 = iax(4,kpr)
                k3 = iax(5,kpr)
                ks = ips(kb)
                ip = 2
                if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
                if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
                wsrk = wsr(ks)
                rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                rik = rik*alat
                dovlap = wsri+wsrk-rik
                dovl1(ip) = dmax1(dovl1(ip),dovlap/rik)
                dovl2(ip) = dmax1(dovl2(ip),dovlap/wsri,dovlap/wsrk)
C              if (dovl1(ip).gt.omax1(ip)+tiny .or.
C     .            dovl2(ip).gt.omax2(ip)+tiny)
C     .          write(stdo,'('' SCLWS2: warning dovl>omax'')')
              enddo
            enddo
          enddo
          if (ipr .ge. 50) then
            fmt = '(1x,'' new rmt'',1x,6f9.5:/(10x,6f9.5))'
            write(fmt(19:19),'(I1)') npcol
            write(fmt(31:31),'(I1)') npcol
            write(stdo,fmt) (wsr(is),is=1,nspec)
          endif
          volfac = ratio
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
          return
        endif
      enddo
      call rx('sclws2: this cannot happen')

  300 format(6x,a,3f13.7)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine sclws2


