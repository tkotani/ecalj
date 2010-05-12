Cdelw1 var ['oiwtmp', '0', 'oiwtmp', '-nirrmx', 'integer', 'iv_w_', 'findes']
Cdelw1 var ['orirr', '0', 'orirr', '3*nirrmx', 'real(8)', 'rv_w_', 'findes']
Cdelw1 var ['oirr', '0', 'oirr', '-nmesh', 'integer', 'iv_w_', 'findes']
Cdelw1 var ['ores', '0', 'ores', 'nirr', 'real(8)', 'rv_w_', 'findes']
Cdelw1 var ['olini', '0', 'olini', '-nirr', 'integer', 'iv_w_', 'findes']
Cdelw1 var ['obas1', '1', 'obas1', '3*(nbas+nsymop)', 'real(8)', 'rv_w_', 'findes']
Cdelw1 var ['oicls1', '1', 'oicls1', 'nbas+nsymop', 'integer', 'iv_w_', 'findes']
Cdelw1 var ['owk', '0', 'owk', '5*nbas', 'real(8)', 'rv_w_', 'findes']
Cdelw1 var ['oiclas', '1', 'oiclas', 'nbas', 'integer', 'iv_w_', 'findes']
Cdelw1 var ['obas', '-1', 'undef', '0', 'undef', 'x', 'findes']
Cdelw1 undel [['call', 'undef', 'wref', 'redef'], 'oiclas', 'findes']
Cdelw1 undel [['call', 'undef', 'wref'], 'obas', 'findes']
Cdelw1 undel [['wref'], 'obas1', 'findes']
Cgetarg...info...           structure ['findes', 'ctrl', 'sctrl']
Cgetarg...info...           use_to_add ['findes', 'm_struc_def']
C     This series of subroutines was adapted from the Stuttgart package
C     to facilitate finding empty spheres on a lattice.
      subroutine findes(sctrl,alat,dclabl,nbas,nclass,nl,
     .  nrclas,mxclas,nsymop,obas,oiclas,plat,symopm,symopv,lock,wsr,z)

      use m_struc_def  !Cgetarg

C- Finds position and radius of possible empty spheres
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   nrxyz :no of divisions made along each lattice vector
Ci   nsymop:number of symmetry operations
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Cio Inputs/Outputs:
Cio  dclabl:class labels, stored as real numbers
Cio  nbas  :number of atoms in the basis
Cio        :number of empty spheres is added
Cio  nclass:number of classes, atoms in same class are symmetry-related
Cio        :number of empty spheres classes is added
Cio  nl    :number of l's
Cio  nrclas:number of atoms in the ith class
Cio  obas  :offset to the heap for basis vectors (scaled by alat)
Cio        :position of empty spheres is added in bas
Cio  oiclas:offset to the heap for iclass
Cio        :iclass of empty spheres is added
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cio        :wsr should be dimensioned wsr(mxclas)
Cio        :On input the first nclass values must be set
Cio  z     :nuclear charge.
Cio        :z should be dimensioned z(mxclas)
Cio        :On input the first nclass values must be set
Cio  lock  :switches that lock sphere radii.
Cio        :lock should be dimensioned lock(mxclas)
Cio        :On input the first nclass values must be set
Cl Local variables
Cl   facvol:sum of sphere volumes = facvol * cell volume
Cl   omax1:max. allowed overlap divided by distance (s1+s2-d)/d<omax1
Cl   omax2:max. allowed overlap divided by radius  (s1+s2-d)/s1<omax2
Cl   rmaxes:maximum radius for empty sphere (in atomic units)
Cl   rmines:minimum radius for empty sphere (in atomic units)
Cr Remarks:
Cr   This subroutines searches possible positions for empty
Cr   (better: interstitial spheres). The algorithm is purely geometric,
Cr   the potential is not used at all.
Cr   Starting from a muffin-tin geometry, the program searches for
Cr   that position where one can put the largest empty sphere.
Cr
Cr   The symmetry of the crystal is preserved, so one needs on input
Cr   all symmetry-operations of the crystal.
Cr
Cr   In order to consider all atoms, the basis vectors must lie in
Cr   the WS-primitive cell :  -0.5 < b1,b2,b3 < 0.5  !!
Cr   with: bas(m) = b1*plat(m,k)+b2*plat(m,2)+b3*plat(m,3)
Cr   A call to ORDBAS garantees that this condition is fulfilled.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nbas,nclass,nl,nrclas(*),nsymop,mxclas,lock(nclass)
Cdelw1 do not change  [['call', 'undef', 'wref'], 'obas', 'findes']
Cdelw1 do not change  [['call', 'undef', 'wref', 'redef'], 'oiclas', 'findes']
      integer obas,oiclas
Cgetarg       double precision sctrl(1),alat,plat(3,3),
Cgetarg      .  symopm(*),symopv(*),wsr(nclass),z(nclass)
      real(8):: alat , plat(3,3) , symopm(*) , symopv(*) , wsr(nclass) 
     . , z(nclass)
      type(s_ctrl)::sctrl

      double precision dclabl(nclass)
C Local variables:
      integer ic,idamax,ierr,iprint,jmax(3),jmin(3),
     .  lmx,nb1,nirr,nirrmx,nmesh,modep(3)
      integer stdo,nglob,nrxyz(3),isw
Cdelw1 do not change  [['wref'], 'obas1', 'findes']
Cdelw1       integer obas1,oicls1,oirr,olini,ores,orirr,oiwtmp,owk
      integer:: obas1
      integer ,allocatable :: iv_w_oicls1(:)
      integer ,allocatable :: iv_w_oirr(:)
      integer ,allocatable :: iv_w_olini(:)
      real(8) ,allocatable :: rv_w_ores(:)
      real(8) ,allocatable :: rv_w_orirr(:)
      integer ,allocatable :: iv_w_oiwtmp(:)
      real(8) ,allocatable :: rv_w_owk(:)

      double precision avstep,bases(3),
     .  ratio,resmx,platcp(9),qlat(9),tiny,vol,volsph,wmax,wsres,
     .  wsrc(mxclas),facvol,rmaxes,rmines,omax1(3),omax2(3),wsmax
      character*8 clabl(mxclas)
      parameter(tiny=1.d-5)
C Heap
      integer w(1)
      common/w/w
C External calls:
C      external  addes,bigges,chkes,cpusrt,cpustp,dcopy,defdr,defi,
C     .          defl,deflmx,dinv33,dscall,errmsg,fillat,renam,
C     .          icopy,idamaxl,iprint,mdeses,rlse,
C     .          rsmesh,sclwsr,symes,volsph
C Intrinsic functions:
C      intrinsic  dabs,dmax1,dmin1,dsqrt,max0

C      call cpusrt('FES ',.false.)

C --- Setup ---
      stdo = nglob('stdo')
Cgetarg...info...           integer :: nrxyz(3) , integer(8) :: sctrl%nesabc 3
Cgetarg...info...           real(8) :: rmines , real(8) :: sctrl%rmines 1
Cgetarg...info...           real(8) :: rmaxes , real(8) :: sctrl%rmaxes 1
Cgetarg       call upack('ctrl nesabc rmines rmaxes',sctrl,nrxyz,rmines,rmaxes,
Cgetarg      .  0,0)

      i_copy_size=size(sctrl%nesabc)
      call i8icopy(i_copy_size,sctrl%nesabc,1,nrxyz,1)
      rmines=sctrl%rmines
      rmaxes=sctrl%rmaxes

Cgetarg...info...           integer :: modep(3) , integer(8) :: sctrl%modep 3
Cgetarg       call upack('ctrl modep',sctrl,modep,0,0,0,0)

      i_copy_size=size(sctrl%modep)
      call i8icopy(i_copy_size,sctrl%modep,1,modep,1)

Cgetarg...info...           real(8) :: facvol , real(8) :: sctrl%sclwsr 1
Cgetarg...info...           real(8) :: omax1(3) , real(8) :: sctrl%omax1 3
Cgetarg...info...           real(8) :: omax2(3) , real(8) :: sctrl%omax2 3
Cgetarg...info...           real(8) :: wsmax , real(8) :: sctrl%wsrmax 1
Cgetarg       call upack('ctrl sclwsr omax1 omax2 wsrmax',sctrl,facvol,omax1,
Cgetarg      .  omax2,wsmax,0)

      facvol=sctrl%sclwsr
      i_copy_size=size(sctrl%omax1)
      call dcopy(i_copy_size,sctrl%omax1,1,omax1,1)
      i_copy_size=size(sctrl%omax2)
      call dcopy(i_copy_size,sctrl%omax2,1,omax2,1)
      wsmax=sctrl%wsrmax

      if (facvol .eq. 0) facvol = 1
      call info2(10,1,0,' --- FINDES : find new empty spheres'//
     .    '%?#n#  Constrain rmax<=%d## --- ',isw(wsmax.ne.0),wsmax)
      call info2(30,0,1,'%5fUse:  rmines = %d  rmaxes = %d  ',rmines,
     .  rmaxes)
      do  ic = 1, nclass
        call r8tos8(dclabl(ic),clabl(ic))
      enddo

C --- Scale by alat ---
      call dscal(nclass,1.d0/alat,wsr   ,1)
      call dscal(1     ,1.d0/alat,rmaxes,1)
      call dscal(1     ,1.d0/alat,rmines,1)

C --- Use the most compact unit cell ---
      call dcopy(9,plat(1,1),1,platcp,1)
      call dinv33(platcp,1,qlat,vol)
      vol = dabs(vol)*alat**3
Cdelw1 warning(1) , probably  oiclas  is not defined yet at linenumber= 146
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 146
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
      call ordbas(nbas,plat,w(obas),w(oiclas))

C --- Get irreducible mesh points ---
      nmesh = nrxyz(1)*nrxyz(2)*nrxyz(3)
      nirrmx = 1000 + 4*nmesh/nsymop/3
C ... nirrmx should be nmesh but this requires to much space
      ierr=0
Cdelw1    10 call defi (oiwtmp,-nirrmx)
 10   allocate(iv_w_oiwtmp(abs(-nirrmx)))
      if (-nirrmx<0) iv_w_oiwtmp(:)=0

Cdelw1       call defdr(orirr ,3*nirrmx)
      allocate(rv_w_orirr(3*nirrmx))
Ctakao_ZeroClear_NotRequiered        if (3*nirrmx<0) rv_w_orirr(:)=0.0d0

Cdelw1       call defi (oirr  ,-nmesh)
      allocate(iv_w_oirr(abs(-nmesh)))
      if (-nmesh<0) iv_w_oirr(:)=0

Cdelw1       call rsmesh(avstep,ierr,w(oirr),w(oiwtmp),jmax,jmin,nrxyz(1),
Cdelw1      .            nrxyz(2),nrxyz(3),nirr,nirrmx,nsymop,platcp,qlat,
Cdelw1      .            w(orirr),symopm,symopv)
      call rsmesh ( avstep , ierr , iv_w_oirr , iv_w_oiwtmp , jmax 
     . , jmin , nrxyz ( 1 ) , nrxyz ( 2 ) , nrxyz ( 3 ) , nirr , nirrmx 
     . , nsymop , platcp , qlat , rv_w_orirr , symopm , symopv )

      if (ierr.ne.0) then
Cdelw1 rlse name= oiwtmp old_list= oiwtmp orirr oirr
Cdelw1 rlse name= oiwtmp new_list= (None)
Cdelw1         call rlse(oiwtmp)
        if (allocated(iv_w_oirr)) deallocate(iv_w_oirr)
        if (allocated(rv_w_orirr)) deallocate(rv_w_orirr)
        if (allocated(iv_w_oiwtmp)) deallocate(iv_w_oiwtmp)

        nirrmx = nirrmx + nirrmx
        goto 10
      endif
c     write(stdo,301)avstep*alat

Cdelw1 rlse name= orirr old_list= (None)

Cki#error ERROR, try to release name= orirr ,but list=null at linenumber= 166 list= (None)
Cdelw1       call rlse(orirr)

Cdelw1       call defrr(orirr,3*nirr)
      allocate(rv_w_orirr(3*nirr))
Ctakao_ZeroClear_NotRequiered        if (3*nirr<0) rv_w_orirr(:)=0.0d0

Cdelw1       call defrr(ores,   nirr)
      allocate(rv_w_ores(nirr))
Ctakao_ZeroClear_NotRequiered        if (nirr<0) rv_w_ores(:)=0.0d0

Cdelw1       call defi (olini, -nirr)
      allocate(iv_w_olini(abs(-nirr)))
      if (-nirr<0) iv_w_olini(:)=0


c      call cpustp('MESH',.true.,.true.)

C --- Initialize res and lini ---
Cdelw1       call mdeses(alat,w(ores),w(olini),nirr,nsymop,platcp,qlat,
Cdelw1      .            w(orirr),rmines-avstep,symopm,symopv)
      call mdeses ( alat , rv_w_ores , iv_w_olini , nirr , nsymop , 
     . platcp , qlat , rv_w_orirr , rmines - avstep , symopm , symopv 
     . )

      ratio = volsph(nclass,nrclas,wsr)*(alat**3)/vol

C --- Repeat until no more ES found ---
      nb1 = 1
      do while (.true.)

        if (iprint().ge.20) write(stdo,306)
        call info5(20,1,1,
     .    ' FINDES:  Find next largest available ES to reach '//
     .    'target vol frac = %d.'//
     .    '%N%10fCurrent vol frac = %;3,3d, nbas = %i, nclass = %i',
     .    facvol,ratio,nbas,nclass,0)

C   --- Check whether the volume can be reached for given OMAX ---
C       wsrc is only used as a work array; except contents printed last iter
        call dpcopy(wsr,wsrc,1,nclass,alat)
Cdelw1 warning(1) , probably  oiclas  is not defined yet at linenumber= 192
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 192
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
        call sclwsr(0,nbas,nbas,nclass,alat,plat,w(obas),w(oiclas),
     .    modep,clabl,z,lock,facvol,wsmax,omax1,omax2,wsrc)

        ratio = volsph(nclass,nrclas,wsrc)/vol
C       Exit, if volume target met
        if (dabs(ratio-facvol) .lt. tiny) then
          call info2(20,1,0,' FINDES: reached target, sphere '//
     .      'vol = %d * cell vol',facvol,0)
          goto 20
        endif

C   ... Fill unit cell with potential new sites
Cdelw1 warning(1) , probably  oiclas  is not defined yet at linenumber= 204
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 204
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
Cdelw1 warning(1) , probably  oiclas  is not defined yet at linenumber= 204
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 204
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
Cdelw1         call fillat(alat,w(obas),w(oiclas),jmax,jmin,nb1,nbas,nirr,
Cdelw1      .    nrxyz,w(olini),platcp,qlat,w(ores),w(orirr),
Cdelw1      .    rmaxes,rmines-avstep,wsr)
        call fillat ( alat , w ( obas ) , w ( oiclas ) , jmax , jmin 
     .   , nb1 , nbas , nirr , nrxyz , iv_w_olini , platcp , qlat , rv_w_ores 
     .   , rv_w_orirr , rmaxes , rmines - avstep , wsr )

        nb1 = nbas+1
C   ... Find biggest empty sphere in unit cell with center on mesh-point
Cdelw1 warning(1) , probably  oiwtmp  is not defined yet at linenumber= 209
Cdelw1         call bigges(alat,bases,w(oiwtmp),nirr,w(olini),w(ores),resmx,
Cdelw1      .              w(orirr))
        call bigges ( alat , bases , iv_w_oiwtmp , nirr , iv_w_olini 
     .   , rv_w_ores , resmx , rv_w_orirr )

C   ... Move empty sphere to high-symmetry point, leave the mesh
Cdelw1 warning(1) , probably  oiclas  is not defined yet at linenumber= 212
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 212
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
        if (resmx.gt.tiny)
     .   call symes(alat,w(obas),bases,clabl,w(oiclas),nbas,nrxyz,
     .               nsymop,platcp,qlat,symopm,symopv,wsr,wsres)
        if (resmx.le.tiny .or. wsres.lt.rmines) then
          call info2(20,0,0,' FINDES: Largest empty sphere found '//
     .      'smaller than rmines = %d',rmines*alat,0)
          goto 20
        endif
        if (nclass+1 .gt. mxclas) then
          call info2(20,0,0,' FINDES: number of classes exceeds '//
     .      'mxclas = %i',mxclas,0)
          goto 20
        endif
        wsres = dmin1(wsres,rmaxes)
        wmax = wsr(idamax(nclass,wsr,1))
C   ... Check possible overlap w/ atoms outside parallelipiped
        call chkes(alat,bases,platcp,qlat,wmax,wsres)
C   ... Increase number of classes by 1
        nclass = nclass+1
        ic     = nclass
        wsr(ic) = wsres
        z(ic)   = 0.d0
        call deflmx(lmx,wsr(ic),z(ic))
        nl = max0(nl,lmx+1)
        if (wsr(ic)*alat .lt. 0.5d0) then
          call info2(20,0,0,' FINDES (warning) new ES (class %i)'//
     .      ' radius is very small: r = %d',nclass,wsr(ic)*alat)
        endif
C   ... Find label for new empty sphere
        call renam(clabl,ic,nclass,z)
        if (iprint().ge.20)
     .    write(stdo,303) clabl(ic),wsr(ic)*alat,nclass

C   ... Add all symmetry-related spheres
Cdelw1 do not change obas1 because of ['wref']
        call defdr(obas1 ,3*(nbas+nsymop))
Cdelw1         call defi (oicls1,   nbas+nsymop)
        allocate(iv_w_oicls1(nbas+nsymop))
Ctakao_ZeroClear_NotRequiered          if (nbas+nsymop<0) iv_w_oicls1(:)=0

Cdelw1 warning(1) , probably  oiclas  is not defined yet at linenumber= 248
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  oiclas  is not defined yet at linenumber= 248
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1         call icopy(nbas,w(oiclas),1,w(oicls1),1)
        call icopy ( nbas , w ( oiclas ) , 1 , iv_w_oicls1 , 1 )

Cdelw1 do not change obas1 because of ['wref']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 249
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
        call dcopy(3*nbas,w(obas),1,w(obas1),1)
Cdelw1 do not change obas1 because of ['wref']
Cdelw1 do not change obas1 because of ['wref']
Cdelw1         call addes2(w(obas1),bases,clabl,w(oicls1),nbas,nclass,nrclas,
Cdelw1      .             nsymop,platcp,qlat,symopm,symopv)
        call addes2 ( w ( obas1 ) , bases , clabl , iv_w_oicls1 , nbas 
     .   , nclass , nrclas , nsymop , platcp , qlat , symopm , symopv 
     .   )

Cdelw1 do not change obas1 because of ['wref']
        call redfrr(obas1,3*nbas)
Cdelw1 do not change obas1 because of ['wref']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 253
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
        obas = obas1
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
        call defi(oiclas,nbas)
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1         call icopy(nbas,w(oicls1),1,w(oiclas),1)
        call icopy ( nbas , iv_w_oicls1 , 1 , w ( oiclas ) , 1 )


      enddo
   20 continue
C     call rlse(oiwtmp)

C --- Cleanup ---
C ... Rescale wsr,rmines,rmaxes
      call dscal(nclass,alat,wsr   ,1)
      call dscal(1     ,alat,rmines,1)
      call dscal(1     ,alat,rmaxes,1)
Cdelw1       call defdr(owk,5*nbas)
      allocate(rv_w_owk(5*nbas))
Ctakao_ZeroClear_NotRequiered        if (5*nbas<0) rv_w_owk(:)=0.0d0

Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 267
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
      call ordbas(nbas,plat,w(obas),w(oiclas))
Cdelw1 rlse name= owk old_list= orirr ores olini obas1 oicls1 oiclas owk
Cdelw1 rlse name= owk new_list= orirr ores olini obas1 oicls1 oiclas
Cdelw1       call rlse(owk)
      if (allocated(rv_w_owk)) deallocate(rv_w_owk)

C ... Printout
      call info0(20,0,0,' Writing new basis to file poses ...')
Cdelw1 do not change oiclas because of ['call', 'undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  obas  is not defined yet at linenumber= 271
Cdelw1 do not change obas because of ['call', 'undef', 'wref']
      call prpos(nbas,nclass,clabl,w(oiclas),w(obas),wsrc,z)

      call rx0('findes')

  306 format(/72('='))
  303 format(/' FINDES:  new class: ',a8,', WSR=',f7.5,', class no',i4)

Cgetarg       end
Cdelw1 w_varlist remains: orirr ores olini obas1 oicls1 oiclas
Cdelw1 w_varlistundel: [oiclas obas obas1]
Cdelw1 w_varlist (undel), remains: [obas1 oiclas]
Cdelw1 w_varlist (del), remains: [orirr ores olini oicls1]
Cdelw1 not deallocate oiclas because of [call undef wref redef]
Cdelw1 not deallocate obas1 because of [wref]
      if (allocated(iv_w_oicls1)) deallocate(iv_w_oicls1)
      if (allocated(iv_w_olini)) deallocate(iv_w_olini)
      if (allocated(rv_w_ores)) deallocate(rv_w_ores)
      if (allocated(rv_w_orirr)) deallocate(rv_w_orirr)

      end subroutine findes


      subroutine symes(alat,bas,bases,clabl,iclass,nbas,nrxyz,nsymop,
     .                 plat,qlat,symopm,symopv,wsr,wsres)

C- Tries to move bases to a high-symmetry point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Cio  bases :basis vectors of new empty sphere (scaled by alat)
Ci   clabl :name of the different inequivalent atoms
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   nrxyz :no of divisions made along each primitive lattice vector
Ci   nsymop:number of symmetry operations as supplied by the generators
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (scaled by alat)
Co Outputs:
Cio  bases :basis vectors of new empty sphere (scaled by alat)
Co   wsres :Wigner-Seitz radius of new empty sphere
Cr Remarks
Cr   facmrg:decides when two empty spheres should be merged
Cr   facovl:decides how much empty spheres should overlap
Cr   First all symmetry-related positions of bases are generated.
Cr   Then the empty sphere is replaced by the average of all positions
Cr   within a sphere with radius wsres*facmrg centered around bases.
Cr   The parameter which decides when to merge the spheres: facmrg, is
Cr   set to 1.5 but this number may be changed.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer iclass(*),nbas,nrxyz(3),nsymop
      double precision alat,bas(3,*),bases(3),plat(3,3),qlat(3,3),
     .                 symopm(3,3,*),symopv(3,*),wsr(*),wsres
      character*8 clabl(*)
C Local variables:
      integer ibas,ibasm,ii,isop,iprint,j(3),k(3),ll1,ltmax,lunit,
     .        m,mm,n,nrxyzz(3),nid,nit,nnz,nsbas
      double precision d,db,dbas(3),desesm,desatm,dm1,
     .                 dm1mx,dm2,d3nrm2,desat,deses,facmrg,facovl,fmwsr,
     .                 r(3),rb(3,3),sdbas(3),tiny,tolwsr
      logical ldum,lready
      parameter(ltmax=1,ll1=ltmax*2+1)
      parameter(facmrg=2.d0,facovl=1.d0,tiny=.1d-5,tolwsr=.1d-2)
C External calls:
      external  daxpy,dcopy,dpzero,distes,d3nrm2,getirr,icopy,iprint,
     .          lunit,mdesat,mdeses,popprt,pshprt
C Intrinsic functions:
      intrinsic  dmax1,dmin1,mod
C Statement functions:
      mm(ii,m)=ltmax-(mod(ii,ll1**m)-mod(ii,ll1**(m-1)))/ll1**(m-1)

      if (iprint().ge.60) write(*,300)facmrg,facovl

      call icopy(3,nrxyz,1,nrxyzz,1)

C --- Optimize position of empty sphere (especially low symmetry points)
      nit=0
      lready=.false.
      do while (.not.lready)
        nit=nit+1
        nnz=0
        db=1.d0
        do m = 1, 3
          nrxyz(m)=nrxyz(m)*2
          do n = 1, 3
            rb(n,m) = plat(n,m)/nrxyz(m)
C           qb(n,m) = qlat(n,m)*nrxyz(m)
          enddo
          db=db*d3nrm2(rb(1,m))
        enddo
        db=alat*db**(1.d0/6.d0)
        call getirr(k,nrxyz,qlat,bases)
        dm1mx=0.d0
        do ii=0,ll1**3-1
          do m=1,3
            j(m) = k(m)+mm(ii,m)
            j(m) = mod(j(m),nrxyz(m))
            if (j(m).lt.0) j(m)=j(m)+nrxyz(m)
            if (j(m).ge.nrxyz(m)/2) j(m)=j(m)-nrxyz(m)
          enddo
          do m=1,3
            r(m) = j(1)*rb(m,1) + j(2)*rb(m,2) + j(3)*rb(m,3)
          enddo
          call mdesat(alat,bas,desat,ibas,iclass,nbas,plat,r,wsr)
          call mdeses(alat,deses,ldum,1,nsymop,plat,qlat,r,0.d0,
     .                symopm,symopv)
          dm1=dmin1(desat,deses)
          if (dm1.gt.desat*0.5d0) nnz=nnz+1
          if (dm1.gt.dm1mx) then
            dm1mx=dm1
            desesm=deses
            desatm=desat
            ibasm=ibas
            call dcopy(3,r,1,bases,1)
          endif
          if(iprint().ge.100)
     .      write(*,311)ii,r,desat*alat,deses*alat
        enddo
        lready=db.lt.tolwsr.or.nnz.lt.3
        if (nit.eq.1.or.lready) call pshpr(iprint()+10)
        if (iprint().ge.30) write(*,308)nit,bases,dm1mx*alat
        if (iprint().ge.50)
     .   write(*,305)desesm*alat,clabl(iclass(ibasm)),desatm*alat
        if (iprint().ge.60) write(*,307)db,tolwsr,nnz
        if (nit.eq.1.or.lready) call poppr
      enddo

C --- Consider symmetry-related positions
      fmwsr=facmrg*facovl*desatm+tiny+tiny
      nsbas=1
      nid=1
      call dpzero(sdbas,3)
      if (iprint().ge.50) write(*,301)
      do isop=2,nsymop
        call distes(bases,d,dbas,isop,plat,qlat,symopm,symopv)
        if (d.lt.fmwsr) then
          call daxpy(3,1.d0,dbas,1,sdbas,1)
          nsbas=nsbas + 1
          if (d.lt.tiny) nid = nid + 1
        endif
        if (iprint().ge.50)
     .    write(*,302)isop,dbas,d*alat,fmwsr*alat
      enddo
C --- Merge the ES
      if (nsbas/nid.ne.1) then
        if (iprint().ge.30) write(*,303)bases,nid,nsymop/nid
        call dcopy(3,bases,1,bases,1)
        call daxpy(3,-1.d0/nsbas,sdbas,1,bases,1)
C ----- Now examine merged sphere, get deses,desat
        call mdesat(alat,bas,desat,ibas,iclass,nbas,plat,bases,wsr)
        call mdeses(alat,deses,ldum,1,nsymop,plat,qlat,bases,0.d0,
     .              symopm,symopv)
        dm2=dmin1(desat,deses)
        if (iprint().ge.30) write(*,304)bases,nsbas/nid
        if (iprint().ge.40)
     .    write(*,305)deses*alat,clabl(iclass(ibas)),desat*alat
        if (dm2.lt.dm1mx) then
          call daxpy(3,+1.d0/nsbas,sdbas,1,bases,1)
          if (iprint().ge.40) write(*,306)
        endif
      else
        dm2=0.d0
      endif

      wsres=dmax1(dm1mx,dm2)*facovl

      call icopy(3,nrxyzz,1,nrxyz,1)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

300   format( ' SYMES : FACMRG=',f6.4,',  FACOVL=',f6.4)
301   format(/' SYMES :ISOP',16x,'DBAS',16x,'DIST',5x,'FMWSR',
     .       /8x,57('-'))
302   format(8x,i2,3x,3f10.6,1x,2f10.6)
303   format( ' SYMES : BEFORE: ',3f10.6,',  degeneracy: ',i2,
     .       /9x,'Number of equivalent points in the unit cell:',i2)
304   format( ' SYMES : AFTER:  ',3f10.6,', averaged over',i2,' points')
305   format( ' SYMES : shortest distance to ES of same class :',f10.7,
     .      /'       : shortest distance to another atom ',a8,':',f10.7)
306   format( ' SYMES : spheres not merged')
307   format(' SYMES : DB=',f7.5,', TOLWSR=',f7.5,', NNZ=',i2)
308   format( ' SYMES : STEP ',i2,':',3f10.6,'  WSR=',f8.6)
311   format(2x,i3,3f11.6,2f10.7)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine symes


      subroutine getirr(j,nrxyz,qlat,r)

C- gets coordinates j of irreducible mesh point
C-----------------------------------------------------------------------
Ci Inputs:
Ci   nrxyz :no of divisions made along each primitive lattice vector
Ci   r     :actual coordinate
Ci   qlat  :primitive translation vectors in reciprocal space
Co Outputs:
Co   j     :if r is a mesh point then rirr(m) = sum_k  rb(m,k)*j(k)
Co         :where rirr and r are connected by a lattice vector.
C-----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nrxyz(*),j(3)
      double precision qlat(3,*),r(3)
C Local variables:
      integer iprint,k,lunit
      double precision tiny,x(3)
      logical namp
      parameter(tiny=.1d-5)
C External calls:
      external iprint,lunit
C Intrinsic functions:
      intrinsic dabs,idnint,mod

      namp=.false.
      do k=1,3
        x(k) = (r(1)*qlat(1,k)+r(2)*qlat(2,k)+r(3)*qlat(3,k))*nrxyz(k)
        j(k) = idnint(x(k))
        if (dabs(x(k)-j(k)).gt.tiny) namp=.true.
        j(k) = mod(j(k),nrxyz(k))
        if (j(k).lt.0) j(k)=j(k)+nrxyz(k)
      enddo
      if (iprint().ge.70) write(*,300)r,j
      if (namp.and.iprint().ge.50)
     .  write(*,400)r,x,idnint(x(1)),idnint(x(2)),idnint(x(3))
ct     .  write(lunit(2),400)r,x,idnint(x(1)),idnint(x(2)),idnint(x(3))
300   format( ' GETIRR: VECTOR  ',3f10.6,'  -->',3i6)
400   format(/' GETIRR: the point:',3f8.4,' not a mesh point',
     .       /'         X=',3f8.2,
     .       /'         J=',3(i5,3x))
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine getirr

      subroutine mdesat(alat,bas,desat,ibas,iclass,nbas,plat,r,wsr)

C- Finds the closest atom to a given mesh-point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   bas   :basis vectors (scaled by alat)
Ci   ibas  :index of closest atom
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (scaled by alat)
Co Outputs:
Cio  desat  :distance from mesh-point to nearest sphere
Cr Remarks:
Cr Contrary to fillat the distance is even evaluated if it is smaller
Cr than rmines.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer ibas,iclass(*),nbas
      double precision alat,bas(3,*),desat,plat(*),r(3),wsr(*)
C Local variables:
      integer i1,i2,i3,iprint,jbas
      double precision c01,c02,c03,c11,c12,c13,c21,c22,c23,
     .                 d1mach,d2,dx,dy,dz,r1,r2,r3,res,rpw2,wsri
C External calls:
      external d1mach,iprint
C Intrinsic functions:
      intrinsic dsqrt

      r1=r(1)
      r2=r(2)
      r3=r(3)
      desat=d1mach(2)
      do jbas=1,nbas
        rpw2=d1mach(2)
        wsri=wsr(iclass(jbas))
        c01=bas(1,jbas)-r1
        c02=bas(2,jbas)-r2
        c03=bas(3,jbas)-r3
        do i1=-1,1
          c11=c01+i1*plat(1)
          c12=c02+i1*plat(2)
          c13=c03+i1*plat(3)
          do i2=-1,1
            c21=c11+i2*plat(4)
            c22=c12+i2*plat(5)
            c23=c13+i2*plat(6)
            do i3=-1,1
              dx=c21+i3*plat(7)
              dy=c22+i3*plat(8)
              dz=c23+i3*plat(9)
              d2=dx*dx+dy*dy+dz*dz
              if (rpw2.gt.d2) then
                res=dsqrt(d2)-wsri
                rpw2=res+wsri
                rpw2=rpw2*rpw2
                if (res.lt.desat) ibas=jbas
              endif
            enddo
          enddo
        enddo
        desat=dmin1(desat,res)
        if (iprint().ge.70) write(*,300) ibas,res*alat,desat*alat
      enddo
300   format(' MDESAT: IBAS=',i3,', RES=',f11.6,', DESAT=',f11.6)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine mdesat


      subroutine distes(bases,d,dbas,isop,plat,qlat,symopm,symopv)

C- Calculates the distance between symmetry-related sites
C ----------------------------------------------------------------------
Ci Inputs:
Ci   bases :basis vectors of new empty sphere (scaled by alat)
Ci   isop  :number of symmetry opration
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Co Outputs:
Ci   d     :distance between symmetry-related sites
Cr Remarks:
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer isop
      double precision bases(*),d,dbas(*),plat(*),qlat(*),symopm(9,*),
     .                 symopv(3,*)
C Local variables:
      integer m
      double precision bast(3),d3nrm2
C External calls:
      external  daxpyl,dmpy,d3nrm2,reducv
C Intrinsic functions:
      intrinsic  dsqrt

      call dmpy(symopm(1,isop),3,1,bases,3,1,bast,3,1,3,1,3)
      call daxpy(3,1.d0,symopv(1,isop),1,bast,1)
      do m=1,3
        dbas(m)=bases(m)-bast(m)
      enddo
      call reducv(plat,qlat,dbas,2)
      d=dsqrt(d3nrm2(dbas))

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine distes


      subroutine chkes(alat,bases,plat,qlat,wmax,wsres)

C- Checks possible overlap of empty sphere and atoms outside
C- the parallelipiped.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bases :basis vectors of new empty sphere (scaled by alat)
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive reciprocal lattice vectors
Ci   wmax  :largest Wigner-Seitz radius (scaled by alat)
Ci   wsres :Wigner-Seitz radius of empty sphere (scaled by alat)
Cr Remarks:
Cr   the minimal distance of the empty sphere to the faces of the
Cr   parallelipiped is calculated.
Cr   If this distance is smaller the the sum of the empty sphere radius
Cr   and the largest radius it could be possible that there is overlap
Cr   of the empty sphere and another sphere outside the parallelipiped.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision alat,bases(3),plat(3,3),qlat(3,3),wmax,wsres
C Local variables:
      integer i,iprint,k,lunit,m
      double precision ddot,dmin,qn(3,3),vec(3)
      character*288 messg
C External calls:
      external  ddot,errmsg,iprint,lunit,nrmliz
C Intrinsic functions:
      intrinsic  dabs,dmin1

      dmin=1.d20
      call nrmliz(3,qlat,qn)
      do k=1,3
        do i=-1,1,2
          do m=1,3
            vec(m)= 1.5d0*i*plat(m,k)-bases(m)
          enddo
          dmin=dmin1(dmin,dabs(ddot(3,vec,1,qn(1,k),1)))
        enddo
      enddo
      if (wsres+wmax-dmin.gt.0.d0) then
        write(messg,400)dmin*alat,wsres*alat,wmax*alat,
     .                     (dmin-wsres-wmax)*alat
        call errmsg(messg,1)
      elseif (iprint().ge.50) then
        write(*,300)dmin*alat,wsres*alat,wmax*alat,
     .                     (dmin-wsres-wmax)*alat
      endif

300   format(/' CHKES : DMIN=',f8.5,', WSR=',f8.5,', WMAX=',f8.5,
     .       /'         DMIN-WSR-WMAX=',f8.5)
400   format( ' CHKES : DMIN=',f8.5,', WSR=',f8.5,', WMAX=',f8.5,
     .        '|        DMIN-WSR-WMAX=',f8.5,
     .        '|        possible overlap of empty sphere',
     .                ' with atom out of range',
     .        '|        use most compact unit cell.$')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine chkes


      subroutine deflmx(lmx,wsr,z)

C- Returns default value of lmx for given nuclear charge and wsr
C ----------------------------------------------------------------------
Ci Inputs:
Ci   wsr   :Wigner-Seitz sphere radius (in atomic units)
Ci   z     :nuclear charge
Co Outputs:
Co   lmx   :lmx(j) = maximum l for atom j
Cr Remarks:
Cr The minimal lmx depends on the radius and is:
Cr
Cr       lmx=1  for 0.0 to 1.5
Cr       lmx=2  for 1.5 to 3.5
Cr       lmx=3  for 3.5 to 5.5
Cr       lmx=4  for 5.5 to 7.5
Cr       ....
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer lmx
      double precision wsr,z
C Local variables:
      integer iz,h,rb
      parameter (h=1,rb=37)
      character*72 messg
C External calls:
      external errmsg
C Intrinsic functions:
      intrinsic idnint,int,max0

      lmx=int(0.5d0*wsr+1.25d0)

      iz=idnint(z)
      if (wsr.lt.0.d0.or.wsr.gt.10.d0) then
        write(messg,400)wsr
        call errmsg(messg,1)
      endif

      if (iz.ge.rb) then
        lmx=max0(3,lmx)
      elseif (iz.ge.h) then
        lmx=max0(2,lmx)
      endif

400   format(' DEFLMX: bad Wigner Seitz radius WSR=',f10.6,'$')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine deflmx


      subroutine renam(clabl,ic,nclass,z)

C- Makes a new label for a class
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ic    :class label
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   z     :nuclear charge
Cio Inputs/Outputs:
Cio  clabl :name of the different inequivalent atom
Cio         nclass clabl's are needed on input
Cio         clabl(ic) is changed on output
Cr Remarks:
Cr A new clabl is found for class ic.
Cr If class ic is the only one with that nuclear charge, then
Cr the label will be H, He, ....
Cr If there are more than one classes of the nclass classes with
Cr the same nuclear charge, the label of class ic will be XXn with
Cr XX= H, He, ...  and n is the smallest nonzero integer which is not
Cr yet used in a label of a class ic < jc.
Cw Warning: if RENAME is called for class ic, it should also be called
Cw for classes ic+1, ..., nclass (else it could happen that different
Cw classes have the same label).
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer ic,nclass
      double precision z(*)
      character*8  clabl(*)
C Local variables:
      integer iat,iz,jc,len
      logical exist,lonly
      character*8 fmt
C External calls:
      external zclabl
C Intrinsic functions:
      intrinsic idnint

      iz=idnint(z(ic))
      lonly=.true.
      do jc=1,nclass
        if (ic.ne.jc.and.idnint(z(jc)).eq.iz) lonly=.false.
      enddo
      call zclabl(-1,clabl(ic),iz)
      len=2
      if (clabl(ic)(2:2).eq.' ') len=1
      iat=1
      exist=.true.
      do while (exist)
        if (iat.ne.1.or..not.lonly) then
          call getfmt(' ',iat,1,2,1,3,fmt)
          write(clabl(ic)(len+1:),fmt=fmt)iat
        endif
        exist=.false.
        do jc = 1,ic-1
          if (clabl(ic).eq.clabl(jc)) exist=.true.
        enddo
        iat=iat+1
      enddo
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine renam


      subroutine zclabl(iopt,clabl,iz)

C- Returns name of compound
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iopt  :-1  z    =input  ; clabl=output
Ci           1  clabl=input  ; z    =output
Ci Inputs/Outputs:
Cio  clabl :name of chemical formula
Cio  iz    :nuclear charge
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer iopt,iz
      character*(*) clabl
C Local variables:
      integer icha,lunit
      character*2 aclabl(0:100)
      character*8 ccopy
C External calls:
      external  chcase,errmsg,lunit
C Intrinsic functions:
      intrinsic  ichar
C Data statements:
      data aclabl/'E ',
     .            'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     .            'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     .            'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     .            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     .            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     .            'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     .            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     .            'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     .            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     .            'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/

      if (iopt.eq.-1) then
        if (iz.lt.0.or.iz.gt.100) then
ct          write(lunit(2),400)iz
          write(*,400)iz
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
          return
        endif
        clabl=aclabl(iz)
      elseif (iopt.eq.1) then
        ccopy=clabl
        call  chcase( 1,1,ccopy(1:1))
        call  chcase(-1,1,ccopy(2:2))
        icha=ichar(ccopy(2:2))
        if (icha.lt.97.or.icha.gt.122)ccopy(2:2)=' '
        do iz=0,100
          if (aclabl(iz).eq.ccopy(1:2)) then
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
            return
          endif
        enddo
ct        write(lunit(2),401)clabl
        write(*,401)clabl
        iz=-1
      else
        call errmsg(' ZCLABL: bad iopt.$',5)
      endif

400   format(/' ZCLABL: bad Z=',i5,' set to -1')
401   format(/' ZCLABL: could not find atom ',a8,' in my list,',
     .       ' Z set to -1')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine zclabl


      subroutine getfmt(chr,res,nelts,icast,nwmax,ndigmx,fmt)

C- Supplies the format for writing a token
C ----------------------------------------------------------------------
Ci Inputs:
Ci   chr   :array of character
Ci   res   :array of logical,integer,real or double precision
Ci   nelts :number of elements in the array.
Ci   icast :integer controlling the type of the token(s):
Ci         :0=,logical, 1=char, 2=int, 3=real, 4=double
Ci   nwmax :maximum repetition number
Ci   ndigmx:maximum possible number of digits.
Co Outputs:
Co   fmt   :format descriptor
Co          EX: (11f12.5) (3 a 6)
Cr Remarks:
C     For a character token, nrdig returns the number of nonblank
C     characters in an array of strings.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer icast,ndigmx,nelts,nwmax
      integer res(*)
      character*1 chr(*)
      character*(*) fmt
C Local variables:
      integer chrpos,i,icount,idig,k,lenchr,n,nrdig,nw,resi
      double precision d1mach,dpres,zero0,zero,resd,tiny
      real resr
      logical resl
      parameter (tiny=1.d-6)
C External calls:
      external  chrpos,d1mach,tocast
C Intrinsic functions:
      intrinsic dabs,dlog10,dnint,iabs,max0,min0

      n=0
      nrdig=0
      zero0=0.1d0**ndigmx
      do icount=1,iabs(nelts)
        call tocast(res,res,res,res,resl,resi,resr,resd,1,icount,icast)
        i=1
        if (icast.eq.1) then
          lenchr=res(1)
          k=1+lenchr*(icount-1)
          i=chrpos(chr(k),' ',1,lenchr)-1
        elseif (icast.eq.2) then
          if (iabs(resi).ne.0)   i=alog10(.5+iabs(resi))+1
          if (resi.lt.0) i=i+1
        elseif (icast.ge.3) then
          if (icast.eq.3) dpres=dble(resr)
          if (icast.eq.4) dpres=resd
          zero=zero0
          do idig=0,ndigmx-1
            if (dabs(dpres-dnint(dpres)).lt.zero) goto 10
            dpres=dpres*10.0d0
            zero=zero*10.0d0
          enddo
  10      continue
          if (icast.eq.3) dpres=dble(resr)
          if (icast.eq.4) dpres=resd
          if (idig.ne.0) i = 0
          if (dabs(dpres).ge.1.d0-tiny) i=dlog10(tiny+dabs(dpres))+1
          if (dpres.lt.-zero0) i=i+1
        endif
        nrdig=max0(nrdig,idig)
        n=max0(n,i)
      enddo
      nw=min0(nelts,nwmax)
      if(icast.eq.0) then
        write(fmt,300) nw,'l',1
      elseif(icast.eq.1)then
        write(fmt,300) nw,'a',n
      elseif(icast.eq.2)then
        write(fmt,300) nw,'i',n
      elseif(icast.eq.3) then
        write(fmt,301) nw,'e',n+nrdig+1,nrdig
      elseif(icast.eq.4) then
        write(fmt,301) nw,'f',n+nrdig+1,nrdig
      endif
300   format('(',i2,a1,i2,')')
301   format('(',i2,a1,i2,'.',i2,')')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine getfmt


      subroutine tocast(l,i,r,d,l1,i1,r1,d1,m,n,icast)

C ----------------------------------------------------------------------
Ci Inputs:
Ci   l     :array of logical
Ci   i     :array of integer
Ci   r     :array of real
Ci   d     :array of double precision
Co Outputs:
Ci   l1    :array of logical
Ci   i1    :array of integer
Ci   r1    :array of real
Ci   d1    :array of double precision
Cr Remarks:
C     For a character token, nrdig returns the number of nonblank
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer n,m,icast
      real    r(*),r1(*)
      integer i(*),i1(*)
      logical l(*),l1(*)
      double precision d(*),d1(*)

      if (icast.eq.0) then
        l1(m)=l(n)
      elseif (icast.eq.2) then
        i1(m)=i(n)
      elseif (icast.eq.3) then
        r1(m)=r(n)
      elseif (icast.eq.4) then
        d1(m)=d(n)
      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine tocast


      subroutine addes2(bas,bases,clabl,iclass,nbas,nclass,nrclas,
     .             nsymop,plat,qlat,symopm,symopv)

C-Adds the new empty sphere basis atom and its symmetry related atoms
C ----------------------------------------------------------------------
Ci Inputs:
Cio  bas   :basis vectors (scaled by alat)
Cio  bases :basis vectors of new empty sphere (scaled by alat)
Ci   clabl :name of the different inequivalent atom
Ci   iclass:the jth atom belongs to class iclass(j)
Cio  nbas  :number of atoms in the basis
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   nsymop:number of symmetry operations as supplied by the generators
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Co   symopm:symmetry operation matrix
Co   symopv:symmetry operation vector (scaled by alat)
Co Outputs:
Cio  bas   :basis vectors (scaled by alat)
Cio         on output list has been completed by the new positions
Cio  nbas  :number of atoms in the basis
Co   nrclas:number of atoms in the ith class
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer iclass(*),nbas,nclass,nrclas(*),nsymop
      double precision plat(*),qlat(*),bas(3,*),bases(3),symopm(3,3,*),
     .                 symopv(3,*)
      character*8 clabl(*)
C Local variables:
      integer i,isop,iprint,lunit,nbasnw,ibas,jbas
      double precision bas1(3)
      logical latrel,ladd
      character*72 messg
C External calls:
      external  daxpy,dcopy,dmpy,errmsg,iprint,latrel,lunit,reducv

C     ifi = fopn('espos')

      nbasnw=nbas
      do isop=1,nsymop
        call dmpy(symopm(1,1,isop),3,1,bases,3,1,bas1,3,1,3,1,3)
        call daxpy(3,1.d0,symopv(1,isop),1,bas1,1)
        jbas=0
        ladd=.true.
        do while (jbas.lt.nbasnw.and.ladd)
          jbas=jbas+1
          if (latrel(qlat,bas1,bas(1,jbas))) then
            ladd=.false.
            if (iclass(jbas).ne.nclass) then
              write(messg,400)jbas,iclass(jbas)
              call errmsg(messg,4)
            endif
          endif
        enddo
        if (ladd) then
          nbasnw = nbasnw+1
          iclass(nbasnw) = nclass
          call reducv(plat,qlat,bas1,2)
          call dcopy(3,bas1,1,bas(1,nbasnw),1)
        endif
      enddo

ct      if (iprint().ge.10) then
ct        write(lunit(1),300)nbas,nbasnw,nbasnw-nbas
ct        write(lunit(1),301)(clabl(iclass(ibas)),iclass(ibas),
ct     .                         (bas(i,ibas),i=1,3),ibas=nbas+1,nbasnw)
ct      endif

      write(*,300)nbas,nbasnw,nbasnw-nbas
      write(*,301)(clabl(iclass(ibas)),iclass(ibas),
     .                      (bas(i,ibas),i=1,3),ibas=nbas+1,nbasnw)
ct        write(ifi,301) (clabl(iclass(ibas)),iclass(ibas),
ct     .                      (bas(i,ibas),i=1,3),ibas=nbas+1,nbasnw)
      nrclas(nclass)=nbasnw-nbas
      nbas=nbasnw

300   format(/' ADDES : The basis has been extended from ',i3,
     .         ' to ',i3,' atoms.',
     .       /9x,'The ',i2,' new positions are: ')
301   format(9x,'ATOM=',a8,' ICLASS=',i3,1x,'POS=',3f9.5)

400   format(' ADDES : empty sphere sits on atom ',i3,', class',i3,'$')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine addes2


      logical function latrel(qlat,pos1,pos2)

C- Checks if two positions differ by a lattice vector
C ----------------------------------------------------------------------
Ci Inputs:
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   pos1  :first position
Ci   pos2  :second position
Co Outputs:
Co   latrel:.true. if positions differ by a lattice vector
Cr Remarks:
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision qlat(3,*),pos1(*),pos2(*)
C Local variables:
      integer i,m
      double precision diff(3),tolpos,vdiff
      parameter(tolpos=4.d-3)
C Intrinsic functions:
      intrinsic  dabs,dnint

      latrel=.false.
      do i=1,3
        diff(i)=pos1(i)-pos2(i)
      enddo
      do m=1,3
        vdiff=diff(1)*qlat(1,m)+diff(2)*qlat(2,m)+diff(3)*qlat(3,m)
        vdiff=dabs(vdiff-dnint(vdiff))
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        if (vdiff.gt.tolpos) return
      enddo
      latrel=.true.
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end function latrel


      subroutine prpos(nbas,nclass,clabl,iclass,pos,wsr,z)

C     implicit none
      integer nbas,iclass(*),i,ifi,ibas,fopn,nclass,icls
      double precision  pos(3,nbas)
      double precision wsr(*),z(*)
      character*(8) clabl(*)

      ifi = fopn('poses')
      do ibas = 1,nbas
        write(ifi,301) clabl(iclass(ibas)),(pos(i,ibas),i=1,3)
      end do

      write(ifi,*)

      do icls = 1,nclass
        write(ifi,302) clabl(icls),z(icls),wsr(icls)
      end do

  301 format(9x,'ATOM=',a8, 1x,'POS=',3f12.7)
  302 format(9x,'ATOM=',a8, 1x,'Z = ',f5.1,' R = ',f6.3)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine prpos


      subroutine bigges(alat,bases,iwtmp,nirr,lini,res,resmx,rirr)

C- Finds biggest emtpy sphere in the irreducible part of the unit cell
C ----------------------------------------------------------------------
Ci Inputs :
Ci   alat  :length scale
Ci   iwtmp :weight of mesh-point
Ci   lini  :logical bitmap, T if point closer than rmines to a sphere
Ci   n1-3  :no of divisions made along each lattice vector
Ci   nirr  :number of irreducible mesh points
Ci   res  :square of radius of empty sphere (in units of mesh length)
Ci   rirr  :coordinates of irreducible mesh points
Co Outputs:
Co   bases :basis vectors of new empty sphere (scaled by alat)
Co   resmx:largest res, zero if all lini are .true.
Cr Remarks:
Cr   only those points with lini=.false. are considered
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer iwtmp(*),nirr
      double precision alat,bases(3),rirr(3,*),res(*),resmx
      logical lini(*)
C Local variables:
      integer iirr,iprint,iw,iwtmp0,k,lunit
      double precision r2
C External calls:
      external dcopyl,iprint,lunit

ct      if (iprint().ge.50) write(lunit(1),301)
      if (iprint().ge.50) write(*,301)
      resmx=0.d0
      iwtmp0=48
      do iirr=1,nirr
        if(.not.lini(iirr)) then
          iw=iwtmp(iirr)
          r2=res(iirr)
          if (r2.gt.resmx.or.(r2.eq.resmx.and.iw.lt.iwtmp0)) then
            call dcopy(3,rirr(1,iirr),1,bases,1)
            resmx=r2
            iwtmp0=iw
          endif
          if (iprint().ge.60) then
            if (r2.ge.resmx.or.iprint().ge.65)
     .       write(*,302)iirr,(rirr(k,iirr),k=1,3),r2*alat,iw
ct     .       write(lunit(1),302)iirr,(rirr(k,iirr),k=1,3),r2*alat,iw
          endif
        endif
      enddo

      if (resmx.ne.0.d0) then
ct        if (iprint().ge.40) write(lunit(1),300)bases,resmx*alat,iwtmp0
        if (iprint().ge.40) write(*,300)bases,resmx*alat,iwtmp0
      else
ct        if (iprint().ge.40) write(lunit(1),303)resmx*alat
        if (iprint().ge.40) write(*,303)resmx*alat
      endif
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

 300  format(/' BIGGES: BASES=  ',3f10.6,'  WSR=',f8.6,' IWT=',i2)
 301  format(/' BIGGES: IIRR',15x,'RIRR',16x,'RES   IWTMP',
     .       /8x,51('-'))
 302  format(7x,i6,3f10.6,2x,f7.5,1x,i4)
 303  format(/' BIGGES:  RESMX=',f8.6)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine bigges


      subroutine fillat(alat,bas,iclass,jmax,jmin,nb1,nbas,nirr,nrxyz,
     .                  lini,plat,qlat,res,rirr,rmaxes,rmines,wsr)

C- Fills the unit cell with atoms; finds the interstitial region
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   irr   :points to irreducible mesh-point
Ci   jmax  :upper limit for coordinates of irreducible points
Ci   jmin  :lower limit for coordinates of irreducible points
Ci   nrxyz :no of divisions made along each lattice vector
Ci   nb1   :filling with atoms nb1,...,nbas
Ci   nbas  :number of atoms in the basis
Ci   nirr  :number of irreducible mesh points
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   rirr  :coordinates of irreducible mesh points
Ci   rmaxes:maximum radius for empty sphere (scaled by alat)
Ci   rmines:minimum radius for empty sphere (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (scaled by alat)
Co Outputs:
Cio  lini  :logical bitmap, T if point closer than rmines to a sphere
Cio  res   :distance from mesh-point to nearest sphere
Cr Remarks:
Cr   Each irreducible mesh point it checked whether it lies
Cr   closer or not to a sphere than rmines.
Cr
Cr   In order to consider all atoms, the basis vectors must lie in
Cr   the WS-primitive cell i.e. with:
Cr
Cr          bas(m) = b1*plat(m,k)+b2*plat(m,2)+b3*plat(m,3)
Cr
Cr   b1,b2,b3 must lie between -0.5 and +0.5 !!!
Cr   Further the unit cell must be compact (be sure that plat was passed
Cr   trough CPPLAT).
Cr   If these two conditions are fullfilled, it is sufficient to sweep
Cr   the 3*3*3=27 unit cells around the central unit cell.
Cr
Cr   The irreducible points come from RSMESH have m-coorrdinates:
Cr
Cr          rirr(m) = g1*plat(m,k)+g2*plat(m,2)+g3*plat(m,3)
Cr
Cr   and gm fulfills:
Cr
Cr       -0.5 <= jmin(m)/nrxyz(m) <= gm <= jmax(m)/nrxyz(m) <= 0.5
Cr
Cr   Atoms with a distance to this paralleliped larger than rmaxes
Cr   are not considered (ins=false).
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer iclass(*),jmax(3),jmin(3),nrxyz(*),nb1,nbas,nirr
      double precision alat,bas(3,*),plat(3,*),qlat(3,*),res(*),
     .                 rirr(3,*),rmaxes,rmines,wsr(*)
      logical lini(*)
C Local variables:
      integer i(3),ibas,ic,ii,iirr,inz,iprint,isi,iz,ll1,ltmax,lunit,
     .        m,mm,stdo,nglob
      double precision a(3),c(3),c1,c2,c3,d3nrm2,d2,d12,dx,dy,dz,f(3),
     .                 pn(3,3),qn(3,3),rpw2,vec1(3),
     .                 vec3(3),wsri,wprmn2,wprmx2
      logical ins
      parameter(ltmax=1,ll1=ltmax*2+1)
C External calls:
      external  cross,d3nrm2,iprint,lunit,nrmliz
C Intrinsic functions:
      intrinsic  dsqrt,iabs,mod
C Statement functions:
      mm(ii,m)=ltmax-(mod(ii,ll1**m)-mod(ii,ll1**(m-1)))/ll1**(m-1)

      stdo=nglob('stdo')

      if (iprint().ge.70.and.iprint().lt.80) write(stdo,300)
      if (iprint().ge.80) write(stdo,302)
      if (iprint().ge.70.and.iprint().lt.80) write(stdo,300)
      if (iprint().ge.80) write(stdo,302)

C --- normalize plat and qlat, and store into pn and qn
      call nrmliz(3,plat,pn)
      call nrmliz(3,qlat,qn)

      do ibas=nb1,nbas
c       call reducv(plat,qlat,bas(1,ibas),2)
        ic = iclass(ibas)
        wsri=wsr(ic)
        wprmx2=(wsri+rmaxes)*(wsri+rmaxes)
        wprmn2=(wsri+rmines)*(wsri+rmines)
        if (iprint().ge.70.and.iprint().lt.80)
     .    write(stdo,301) ibas,ic,(bas(m,ibas),m=1,3),wsr(ic)*alat
ct        write(stdo,301) ibas,ic,(bas(m,ibas),m=1,3),wsr(ic)*alat
        do ii=0,ll1**3-1
          do m=1,3
            i(m)=mm(ii,m)
            if (i(m).lt.0) then
              f(m)=jmin(m)/dfloat(nrxyz(m))
            elseif (i(m).eq.0) then
              f(m)=0.5d0*(jmin(m)+jmax(m))/dfloat(nrxyz(m))
            else
              f(m)=jmax(m)/dfloat(nrxyz(m))
            endif
          enddo
          do m = 1, 3
            c(m)=bas(m,ibas)+i(1)*plat(m,1)+i(2)*plat(m,2)+i(3)*plat(m,3)
            a(m)=f(1)*plat(m,1)+f(2)*plat(m,2)+f(3)*plat(m,3)
            vec1(m)=c(m)-a(m)
          enddo

          isi=iabs(i(1))+iabs(i(2))+iabs(i(3))
          if (isi.eq.2) then
            do m=1,3
              if (i(m).eq.0) iz=m
            enddo
            call cross(vec1,pn(1,iz),vec3)
            d12=d3nrm2(vec3)
          elseif (isi.eq.1) then
            do m=1,3
              if (iabs(i(m)).eq.1)  inz=m
            enddo
            d12=vec1(1)*qn(1,inz)+vec1(2)*qn(2,inz)+vec1(3)*qn(3,inz)
            d12=d12*d12
          elseif (isi.eq.3) then
            d12=d3nrm2(vec1)
          elseif (isi.eq.0) then
            d12=0.d0
          endif
          ins=(d12.lt.wprmx2)
C ------- ins is false when atom ibas in unit-cell corresponding to
C ------- ii can never overlap with an empty sphere in the irreducible
C ------- volume
          if (ins) then
            c1=c(1)
            c2=c(2)
            c3=c(3)
            do iirr=1,nirr
              if (.not.lini(iirr)) then
                dx=c1-rirr(1,iirr)
                dy=c2-rirr(2,iirr)
                dz=c3-rirr(3,iirr)
                d2=dx*dx+dy*dy+dz*dz
                if (d2.lt.wprmn2) then
                  lini(iirr)=.true.
                  res(iirr)=0.d0
                else
                  rpw2=res(iirr)+wsri
                  rpw2=rpw2*rpw2
                  if (d2.lt.rpw2) res(iirr)=dsqrt(d2)-wsri
                endif
              endif
            enddo
          endif
          if (iprint().ge.80)
     .      write(stdo,303)ibas,ic,c,i,d12*alat,wprmx2*alat,ins
ct        write(stdo,303)ibas,ic,c,i,d12*alat,wprmx2*alat,ins
        enddo
      enddo

300   format(/,' FILLAT: fill unit cell:',
     .    /,9x,'IBAS',2x,'ICLASS',12x,'BAS',13x,'WSR',/,9x,50('-'))
301   format(8x,i3,4x,i3,3x,3f8.5,3x,f8.5)
302   format(/,' IBAS  ICLASS',12x,'BAS',9x,
     .           'I1  I2  I3    D12 WPRMX2 INS',/,9x,63('-'))
303   format(1x,i3,4x,i3,3x,3f7.4,3(2x,i2),2x,f6.3,1x,f5.3,2x,l1)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine fillat


      subroutine mdeses(alat,deses,lini,nirr,nsymop,plat,qlat,rirr,
     .                  rmines,symopm,symopv)

C- Initializes lini and deses, taking into account the symmetry.
C-----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci         :needed only for output
Ci   nirr  :number of irreducible points
Ci   nsymop:number of symmetry operations
Ci   plat  :primitive lattice vectors
Co   rirr  :irreducible point
Ci   rmines:minimum radius for empty sphere (scaled by alat)
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Co Outputs:
Co   deses :half of distance from mesh-point to nearest equivalent point
Co         :deses = distance empty sphere - empty sphere
Co         :zero if lini=.true.
Co   lini  :logical bitmap, T if point closer than rmines to a sphere
Co         (see remarks)
Cr Remarks:
Cr
Cr Due to the given symmetry, if an empty sphere is put into the
Cr crystal, empty spheres must be put on all symmetry-related positions.
Cr If the minimum radius of the empty spheres is given,
Cr and if the spheres should not overlap,
Cr then parts of the unit cell can be excluded as possible positions
Cr for empty spheres because an empty sphere would overlap with
Cr an symmetry-related sphere.
C-----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nirr,nsymop
      double precision alat,plat(9),qlat(9),deses(*),rirr(3,*),rmines,
     .                 symopv(3,*),symopm(9,*)
      logical lini(*)
C Local variables:
      integer iirr,isop,iprint,lunit,m,nexcl
      double precision d1,d2,d3,df2min,diff2,d3nrm2,frm2,
     .                 p2min,pmin(3),r1,r2,r3,mat(9,48),tiny,tiny2,
     .                 v1,v2,v3
      parameter(tiny=1.d-6)
C External calls:
      external d3nrm2,gtpmin,iprint,lunit
C Intrinsic functions:
      intrinsic dmin1,dnint,dsqrt

      nexcl=0
      tiny2=tiny*tiny
      frm2=4.d0*rmines*rmines

C --- get shortest lattice vector
      call gtpmin(plat,pmin)
      p2min=d3nrm2(pmin)

C --- mat_ij=symopm_ij - delta_ij
      call dcopy(9*nsymop,symopm,1,mat,1)
      do isop = 1, nsymop
        call daxpy(3,-1.d0,1.d0,0,mat(1,isop),4)
      enddo

      do iirr = 1, nirr
        df2min=p2min
        r1=rirr(1,iirr)
        r2=rirr(2,iirr)
        r3=rirr(3,iirr)
        do isop = 2, nsymop
C ------- d = r - rrot
C ------- where r is transformed to rrot by symmetry operation number is
C ------- v = d in units of plat and always between -.5 and +.5
C ------- dsqrt(diff2) is then the distance (scaled by alat) from r
C ------- to rrot (minus eventually some lattice vectors)

          d1=mat(1,isop)*r1+mat(4,isop)*r2+mat(7,isop)*r3+symopv(1,isop)
          d2=mat(2,isop)*r1+mat(5,isop)*r2+mat(8,isop)*r3+symopv(2,isop)
          d3=mat(3,isop)*r1+mat(6,isop)*r2+mat(9,isop)*r3+symopv(3,isop)
          v1=d1*qlat(1)+d2*qlat(2)+d3*qlat(3)
          v2=d1*qlat(4)+d2*qlat(5)+d3*qlat(6)
          v3=d1*qlat(7)+d2*qlat(8)+d3*qlat(9)
          v1=v1-dnint(v1)
          v2=v2-dnint(v2)
          v3=v3-dnint(v3)
          d1=v1*plat(1)+v2*plat(4)+v3*plat(7)
          d2=v1*plat(2)+v2*plat(5)+v3*plat(8)
          d3=v1*plat(3)+v2*plat(6)+v3*plat(9)
          diff2=d1*d1+d2*d2+d3*d3
          if (diff2.gt.tiny2) then
            if (diff2.lt.frm2) then
              lini(iirr)=.true.
              deses(iirr)=0.d0
              nexcl=nexcl+1
              goto 10
            endif
            df2min=dmin1(df2min,diff2)
          endif
        enddo
        deses(iirr)=dsqrt(df2min)*0.5d0
10      continue
      enddo

C --- Printout
      if (nirr.gt.1) then
c        if (iprint() .ge. 30) write(lunit(1),300) nirr-nexcl,nirr
c        if (iprint() .ge .50) write(lunit(1),301) alat*dsqrt(p2min)
        write(*,300) nirr-nexcl,nirr
        write(*,301) alat*dsqrt(p2min)
        if (iprint() .ge. 60) then
c          write(lunit(1),302)
          write(*,302)
          do iirr = 1, nirr
            if (lini(iirr)) then
c              write(lunit(1),303)iirr,(rirr(m,iirr),m=1,3),rmines*alat
              write(*,303)iirr,(rirr(m,iirr),m=1,3),rmines*alat
            else
c              write(lunit(1),304)iirr,(rirr(m,iirr),m=1,3),
              write(*,304)iirr,(rirr(m,iirr),m=1,3),
     .                           deses(iirr)*alat
            endif
          enddo
        endif
      endif

300   format(/' MDESES:',i6,' mesh points left from',i6)
301   format(8x,' shortest lattice vector',f10.5)
302   format(/' MDESES:  NIRR',7x,'Rx',8x,'Ry',8x,'Rz',8x,'RES0',
     .       /10x,47('-'))
303   format(7x,i5,2x,3f10.4,2x,'<',f9.6)
304   format(7x,i5,2x,3f10.4,1x,f11.6)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine mdeses


      subroutine gtpmin(plat,pmin)

C- gets shortest nozero lattice vector
C ----------------------------------------------------------------------
Ci Inputs:
Ci   plat  :primitive lattice vectors
Co Outputs:
Co   pmin  :shortest nozero lattice vector
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision pmin(3),plat(9)
C Local variables:
      integer i,i1,i2,i3,j,k
      double precision d(3),d1mach,d2,d2min,d3nrm2,rmax
C External calls:
      external d1mach,d3nrm2,latlim
C Intrinsic functions:
      intrinsic dsqrt,dmin1

      d2min=d1mach(2)
      rmax=dsqrt(dmin1(d3nrm2(plat(1)),d3nrm2(plat(4)),d3nrm2(plat(7))))
      call latlim(plat,rmax,i1,i2,i3)
      do i = -i1, i1
        do j = -i2, i2
          do k = -i3, i3
            d(1)=i*plat(1)+j*plat(4)+k*plat(7)
            d(2)=i*plat(2)+j*plat(5)+k*plat(8)
            d(3)=i*plat(3)+j*plat(6)+k*plat(9)
            d2 = d3nrm2(d)
            if (d2.lt.d2min.and.(i.ne.0.or.j.ne.0.or.k.ne.0)) then
              d2min=d2
              call dcopy(3,d,1,pmin,1)
            endif
          enddo
        enddo
      enddo
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine gtpmin


      double precision function d3nrm2(r)

C-  'norm'-square
C ----------------------------------------------------------------------
Ci Inputs:
Ci   r     :vector
Co Outputs:
Co   d3nrm2:norm-square
Cr Remarks:
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision r(3)

      d3nrm2=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end function d3nrm2


      subroutine ordbas(nbas,plat,bas,iclass)

C- Orders basis in ascending order of iclass
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas  :number of atoms in the basis
Ci   plat  :primitive lattice vectors
Co Inputs/Outputs:
Cio  bas   :basis vectors
Cio  iclass:the jth atom belongs to class iclass(j)
Cw Work arrays:
Cw   wk    :work array of length 5*nbas
Cr Remarks:
Cr   Orders bas so that the first positions correspond to iclass=1,
Cr   followed by those corresponding to iclass=2, ....
Cr   Within each class, sites are also ordered
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nbas,iclass(nbas)
      double precision bas(3,nbas),plat(3,3)
C Local variables:
      integer ib,ibp
      double precision danrm2,qlat(3,3),vol,wk(5,nbas)
      integer iprm(nbas)

      call info0(40,1,0,' ORDBAS:  shorten and order basis')

      call dinv33(plat,1,qlat,vol)
      do  ib = 1, nbas
C       call shorbz(bas(1,ib),bas(1,ib),plat,qlat)
        call reducv(plat,qlat,bas(1,ib),2)
        wk(1,ib) = dble(iclass(ib))
        wk(2,ib) = danrm2(bas(1,ib))
        call dcopy(3,bas(1,ib),1,wk(3,ib),1)
      enddo

      call dvheap(5,nbas,wk,iprm,1d-8,1)
      do  ib = 1, nbas
        ibp = iprm(ib)
        iclass(ib) = idnint(wk(1,ibp))
        call dcopy(3,wk(3,ibp),1,bas(1,ib),1)
      enddo

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine ordbas


      subroutine rsmesh(avstep,ierr,irr,iwtmp,jmax,jmin,n1,n2,n3,nirr,
     .                  nirrmx,nsymop,plat,qlat,rirr,symopm,symopv)

C-Divides the unit cell into microcells and finds irreducible part
C-----------------------------------------------------------------------
Ci Inputs:
Ci   n1-3  :no of divisions made along each primitive lattice vector
Ci   nirrmx:maximum number of irreducible points
Ci   nsymop:number of symmetry operations
Ci   plat  :primitive lattice vectors
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Cio Inputs/Outputs
Co   irr   :irr(i1,i2,i3) points to corresponding irreducible point
Co Outputs:
Co   avstep:average step length
Co   ierr  :1 if nirr.gt.nirrmx
Co   iwtmp :weight of mesh-point
Co   jmax  :upper limit for coordinates of irreducible points
Co   jmin  :lower limit for coordinates of irreducible points
Co   nirr  :number of irreducible points
Co   rirr  :irreducible point
Cr Remarks:
Cr  The lattice is divided into nrxyz(1)*nrxyz(2)*nrxyz(3)
Cr  microcells which are parallelipipeds with 8 corners.
Cr  The corners are nodes of the real space mesh in the whole
Cr  lattice unit cell. Some of these will be symmetry-related leaving
Cr  nirr irreducible k-points.
Cr  These are returned in rirr(3,j) j = 1,nirr; for each corner defined
Cr  by the triple (i1,i2,i3), irr(i1,i2,i3) points to the corresponding
Cr  vector in rirr.
Cr  It is possible to operate with integer algebra, if the mesh
Cr  is invariant under the symmetry operations considered.
C-----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer ierr,n1,n2,n3,irr(0:n1-1,0:n2-1,0:n3-1),iwtmp(*),jmin(3),
     .        jmax(3),nirr,nirrmx,nsymop,stdo,nglob
      double precision avstep,plat(3,3),qlat(3,3),symopv(3,*),
     .                 symopm(9,*),rirr(3,*)
C Local variables:
      integer i1,i2,i3,iirr,iprint,isop,ism(9,48),isv(3,48),j1,
     .        j2,j3,jj,k,k1,k2,k3,kk,m,n0,nnc,nrxyz(3),n12,n22,n32
      double precision dasum,diff(12),dmxdev,dnrm2,qb(3,3),rb(3,3),
     .        step(3),tiny,toldev,tmp1(9),tmp2(9)
      logical lcomp(48)
      character*144 messg
      parameter(toldev=1d0,tiny=1.d-3)
C External calls:
      external  dmpy
C               ,errmsg,iprint
C Intrinsic functions:
      intrinsic  dmax1,dsqrt,idnint,max0,min0,mod

      nrxyz(1)=n1
      nrxyz(2)=n2
      nrxyz(3)=n3
      n12=n1/2
      n22=n2/2
      n32=n3/2
      n0=n1*n2*n3

      stdo=nglob('stdo')

C --- Mesh-size ...
      if (iprint().ge.20) write(stdo,298)nrxyz
      do m=1,3
        step(m)=dnrm2(3,plat(1,m),1)/nrxyz(m)
      enddo

      avstep=(step(1)*step(2)*step(3))**(1.d0/3.d0)
      dmxdev=dmax1(step(1)/avstep,step(2)/avstep,step(3)/avstep)-1.d0
      if (dmxdev.gt.toldev .and.iprint().ge.30) then
        write(messg,402)step
        call errmsg(messg,0)
      endif

      do m = 1, 3
        jmin(m)= nrxyz(m)
        jmax(m)=-nrxyz(m)
        do k = 1, 3
          rb(m,k) = plat(m,k)/nrxyz(k)
          qb(m,k) = qlat(m,k)*nrxyz(k)
        enddo
      enddo
      if (iprint().ge.60) then
        write(stdo,300)((rb(m,k),m=1,3),(qb(m,k),m=1,3),k=1,3)
      endif

C --- check the compatiblity of nrxyz with symgrp

      nnc=0
ct      write(*,*) 'nsymop',nsymop
      do isop=1,nsymop
C         write(*,*) 'isop=',isop
        call dmpy(qb,1,3,symopm(1,isop),3,1,tmp1,3,1,3,3,3)
        call dmpy(tmp1,3,1,rb,3,1,tmp2,3,1,3,3,3)
        do m=1,9
          ism(m,isop)=idnint(tmp2(m))
          diff(m)=tmp2(m)-ism(m,isop)
        enddo

        call dmpy(qb,1,3,symopv(1,isop),3,1,tmp1,3,1,3,1,3)

        do m=1,3
          isv(m,isop)=idnint(tmp1(m))
          diff(9+m)=tmp1(m)-isv(m,isop)
        enddo
        lcomp(isop)=dasum(12,diff,1).lt.tiny
        if (.not.lcomp(isop)) then
c        if (iprint().ge.50) write(lunit(2),400)nrxyz,qb,(symopm(k,isop),
c     .                                    k=1,9),rb,tmp2,(tmp1(k),k=1,3)
          write(*,296)
          nnc=nnc+1
C          write(*,*) 'isop=',isop
        endif
      enddo
ct      write(*,*) nnc

      if (nnc.ne.0.and.ierr.eq.0) then
        write(messg,401)nnc,nsymop
        call errmsg(messg,0)
      endif
      write(*,297)

C ---

      nirr = 0
      do i3 = 0, n3-1
        do i2 = 0, n2-1
          do i1 = 0, n1-1
            if (irr(i1,i2,i3) .eq. 0) then
              nirr = nirr+1
              if (nirr.gt.nirrmx) then
c                if (iprint().ge.70)
c     .            call errmsg(' RSMESH: increase nirrmx.$',-1)
                ierr=1
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
                return
              endif
              jj=n0
              do isop = 1, nsymop
                if (lcomp(isop)) then
                  k1=isv(1,isop)+ism(1,isop)*i1+ism(4,isop)*i2
     .                          +ism(7,isop)*i3
                  k2=isv(2,isop)+ism(2,isop)*i1+ism(5,isop)*i2
     .                          +ism(8,isop)*i3
                  k3=isv(3,isop)+ism(3,isop)*i1+ism(6,isop)*i2
     .                          +ism(9,isop)*i3
                  k1 = mod(k1,n1)
                  k2 = mod(k2,n2)
                  k3 = mod(k3,n3)
                  if (k1.lt.0) k1=k1+n1
                  if (k2.lt.0) k2=k2+n2
                  if (k3.lt.0) k3=k3+n3
                  if (irr(k1,k2,k3).eq.0) then
                    iwtmp(nirr)=iwtmp(nirr)+1
                    irr(k1,k2,k3) = nirr
C ----------------- now determine j1,j2,j3
C ----------------- Attention: i1 in [0,n1-1] but k1 in [-n1/2,n1/2-1]
                    if (k1.ge.n12) k1=k1-n1
                    if (k2.ge.n22) k2=k2-n2
                    if (k3.ge.n32) k3=k3-n3
                    kk=k1+k2+k3
                    if (kk.lt.jj) then
                      jj=kk
                      j1=k1
                      j2=k2
                      j3=k3
                    endif
                  endif
                endif
              enddo
              jmin(1)=min0(j1,jmin(1))
              jmax(1)=max0(j1,jmax(1))
              jmin(2)=min0(j2,jmin(2))
              jmax(2)=max0(j2,jmax(2))
              jmin(3)=min0(j3,jmin(3))
              jmax(3)=max0(j3,jmax(3))
              rirr(1,nirr) = j1*rb(1,1) + j2*rb(1,2) + j3*rb(1,3)
              rirr(2,nirr) = j1*rb(2,1) + j2*rb(2,2) + j3*rb(2,3)
              rirr(3,nirr) = j1*rb(3,1) + j2*rb(3,2) + j3*rb(3,3)
c              if (iprint().ge.70)
c     .          write(stdo,301)i1,i2,i3,nirr,(rirr(k,nirr),k=1,3)
            endif
          enddo
        enddo
      enddo

      do m=1,3
        if (iprint().ge.70) write(stdo,302)m,m,jmin(m),jmax(m)
      enddo

      if (iprint() .ge. 30) write(stdo,303) nirr,n0,n1,n2,n3
      if (iprint() .ge. 60) then
        write(stdo,304)
        do iirr = 1, nirr
          write(stdo,305) iirr,(rirr(m,iirr),m=1,3),iwtmp(iirr)
        enddo
      endif

ct     do m=1,3
ct        write(*,302)m,m,jmin(m),jmax(m)
ct     enddo

ct      write(*,303)nirr,n0,n1,n2,n3
ct      write(*,304)
C      do iirr = 1, nirr
C         write(*,305)iirr,(rirr(m,iirr),m=1,3),iwtmp(iirr)
C      enddo

      ierr = 0

  296 format(T10,'nrxyz not compatible with symgrp')
  297 format(T10,'Compatiblity of nrxyz with symgrp -- ok')
  298 format(/' RSMESH: NRXYZ= ',3i4)
  300 format(/' RSMESH:',11x,' RB ',31x,' QB '/3(3x,3f10.5,5x,3f11.5/))
C 301 format(3i5,i8,2x,3f6.3)
  302 format( ' RSMESH: JMIN',i1,', JMAX',i1,'=',2i4)
  303 format(8x,i6,' irreducible mesh-points from',i8,
     .  ' (',3i4,' )')
  304 format(/' RSMESH:  NIRR',7x,'Rx',8x,'Ry',8x,'Rz',6x,'IWTMP',
     .  /9x,46('-'))
  305 format(9x,i5,2x,3f10.4,5x,i2)

C  400 format(' RSMESH: WARNING: NRXYZ=',3i4,
C     .       ' not compatible with symmetry:',
C     .      //24x,3f9.4,/8x,'QB=QLAT*NRXYZ= ',1x,3f9.4,/24x,3f9.4/,
C     .       /24x,3f9.4,/8x,'SYMOPM=        ',1x,3f9.4,/24x,3f9.4/,
C     .       /24x,3f9.4,/8x,'RB=PLAT/NRXYZ  ',1x,3f9.4,/24x,3f9.4/,
C     .       /24x,3f9.4,/8x,'QB^t*SYMOPM*RB=',1x,3f9.4,/24x,3f9.4/,
C     .                  /8x,'QB^t*SYMOPV    ',1x,3f9.4/,
C     .       /' and the last 12 values should be integer values.'/)
  401 format(' RSMESH: mesh not invariant under ',i2,
     .       ' symmetry operations|        (from ',i2,')$')
  402 format(' RSMESH: bad NRXYZ:',
     .       ' different step length along lattice vectors:',
     .       '|',8x,'STEP1=',f7.4,',  STEP1=',f7.4,',  STEP1=',f7.4,'$')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine rsmesh


      subroutine reducv(plat,qlat,vec,flag)

C- Removes real space lattice vectors from a given real space vector
C ----------------------------------------------------------------------
Ci Inputs:
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   plat  :primitive translation vectors in real space
Ci          they must fulfil: sum_k plat_ki * qlat_kj = delta_ij
Cio  vec   :real space test-vector (cartesian components)
Ci   flag  :controls form of vec at output
Ci          1: vec is given as coordinates of plat
Ci          2: vec is given in cartesian coordinates (vec in unit-cell)
Ci          3: vec is given in cartesian coordinates (vec in WS-cell)
Ci          4: vec is given in cartesian coordinates and is a multiple
Ci             of plat/12 -> all symopv fulfill this condition for
Ci             a standard choice of the origin.
Co Outputs:
Cio  vec   :real space test-vector
Co          vector with shortest length after removal of lattice vectors
Cr Remarks:
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer flag
      double precision plat(9),qlat(9),vec(3)
C Local variables:
      integer i,i1,i2,i3,j,k,m,n1,n2,n3
      double precision danrm2,d1mach,tolpos,v2,vec2(3),vdiff(3),vmin,
     .       v(3),v01,v02,v03,v11,v12,v13,v21,v22,v23,v31,v32,v33,x,y,z
      parameter(tolpos=4.d-3)
C External calls:
      external danrm2,d1mach,dcopy,dpzero
C Intrinsic functions:
      intrinsic dnint

      vdiff(1)=vec(1)*qlat(1)+vec(2)*qlat(2)+vec(3)*qlat(3)
      vdiff(2)=vec(1)*qlat(4)+vec(2)*qlat(5)+vec(3)*qlat(6)
      vdiff(3)=vec(1)*qlat(7)+vec(2)*qlat(8)+vec(3)*qlat(9)
      vdiff(1)=vdiff(1)-dnint(vdiff(1)-d1mach(3))
      vdiff(2)=vdiff(2)-dnint(vdiff(2)-d1mach(3))
      vdiff(3)=vdiff(3)-dnint(vdiff(3)-d1mach(3))
      if (flag.eq.1) then
        call dcopy(3,vdiff,1,vec,1)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
      endif

      if (flag.eq.4) then
        do i=1,3
          v(i)=dnint(12.d0*vdiff(i))/12.d0
          if (dabs(vdiff(i)-v(i)).lt.tolpos) vdiff(i)=v(i)
        enddo
      endif

C --- Now the coordinates of the reduced vectors in the basis plat
C --- are all between -.5 and + .5
      vec(1)=vdiff(1)*plat(1)+vdiff(2)*plat(4)+vdiff(3)*plat(7)
      vec(2)=vdiff(1)*plat(2)+vdiff(2)*plat(5)+vdiff(3)*plat(8)
      vec(3)=vdiff(1)*plat(3)+vdiff(2)*plat(6)+vdiff(3)*plat(9)

      if (flag.eq.3) then
C ----- Try shortening by adding basis vectors
        vmin=d1mach(2)
        v01=vec(1)
        v02=vec(2)
        v03=vec(3)
        do i1=-1,1
          v11=v01+i1*plat(1)
          v12=v02+i1*plat(2)
          v13=v03+i1*plat(3)
          do i2=-1,1
            v21=v11+i2*plat(4)
            v22=v12+i2*plat(5)
            v23=v13+i2*plat(6)
            do i3=-1,1
              v(1)=v21+i3*plat(7)
              v(2)=v22+i3*plat(8)
              v(3)=v23+i3*plat(9)
              v2=danrm2(v)
              if (v2.lt.vmin) then
                vmin=v2
                call dcopy(3,v,1,vec,1)
              endif
            enddo
          enddo
        enddo
      endif
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine reducv


