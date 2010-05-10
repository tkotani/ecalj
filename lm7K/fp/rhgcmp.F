Cdelw1 var ['oyl', '0', 'oyl', 'ng*nlmtop', 'real(8)', 'rv_w_', 'rhgcmp']
Cdelw1 var ['og2', '0', 'og2', 'ng', 'real(8)', 'rv_w_', 'rhgcmp']
Cdelw1 var ['og', '0', 'og', 'ng*3', 'real(8)', 'rv_w_', 'rhgcmp']
Cdelw1 var ['ogv', '-1', 'undef', '0', 'undef', 'x', 'rhgcmp']
Cdelw1 var ['oiv', '0', 'oiv', 'ng*3', 'integer', 'iv_w_', 'rhgcmp']
Cdelw1 var ['ocs', '0', 'ocs', 'ng', 'real(8)', 'rv_w_', 'rhgcmp']
Cdelw1 var ['osn', '0', 'osn', 'ng', 'real(8)', 'rv_w_', 'rhgcmp']
Cdelw1 var ['oqkl', '1', 'oqkl', 'nlm*(kmax+1)', 'real(8)', 'rv_w_', 'rhgcmp']
Cdelw1 undel [['undef', 'wref'], 'ogv', 'rhgcmp']
Cgetarg...info...           structure ['rhgcmp', 'lat', 'slat']
Cgetarg...info...           structure ['rhgcmp', 'spec', 'sspec']
Cgetarg...info...           structure ['rhgcmp', 'site', 'ssite']
Cgetarg...info...           use_to_add ['rhgcmp', 'm_struc_def']
      subroutine rhgcmp(mode,ib1,ib2,ssite,sspec,slat,orhat,kmax,ng,cg)
       
       use m_struc_def  !Cgetarg

C- Adds density of compensating gaussians to FT list
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : a compound of digits specifying what is to be included
Ci         : in the expansion coefficients
Ci         : 1s   digit = 1 add local density rho1-rho2
Ci         :              2 add local density rho1
Ci         :              3 add local density rho2
Ci         : 10s  digit = 1 add core density rhoc
Ci         :              2 add -1 * core density from sm-hankel
Ci         :                in the local density, restoring it
Ci         :                by adding the sm-hankel to the FT mesh
Ci         :              3 combination 1+2
Ci         : 100s digit = 1 add -1 * nuclear density Z delta(r)
Ci         :                In this mode, Z is smoothed into the G_kL
Ci         :              2 add -1 * nuclear density Z delta(r)
Ci         :                In this mode, Z is incporporated directly
Ci         :                in a PW expansion (Z is not smoothed).
Ci         :
Ci         :Examples:
Ci         :mode=130 include the core, the core tail and nuclear charges
Ci         :         This should make the system charge-neutral.
Ci         :mode=131 Like mode=130, but exclude nuclear charge.
Ci         :         The system should have net charge sum_z
Ci         :mode=2   Exclude all core charges, i.e. gaussian (qcorg-z)
Ci         :  and qcorh from the foca Hankel density.
Ci         :  The system should have the valence charge.
Ci         :3 Like 0, but include nuclear charge -Z delta(r) 
Ci         :  directly in a PW expansion (Z is not smoothed).
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: rhogkl
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl rg
Ci     Stored:    *
Ci     Passed to: corprm rhogkl
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat qlat nabc ng ogv okv vol
Ci     Stored:    *
Ci     Passed to: *
Ci   w(orhat):vector of offsets to local site density arrays
Ci   ng    :number of G-vectors
Co Outputs
Co   cg    :FT of local densities is added to cg, depending on mode.
Cr Remarks
Cr   The local charges inside each augmentation sphere
Cr   (including -1 * the core tail) are smoothed by expanding
Cr   in a  G_kL expansion for k=0..kmax.  The latter is
Cr   subsequently converted into a PW expansion.
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   23 Oct 01 rhgcmp now expands local densities in
Cu             GkL for k=0..kmax, l=1..nlml for each site
Cu             Recovers old rhgcmp for kmax=0.  New argument list.
Cu   09 Feb 01 Added mode
Cu   30 May 00 Adapted from nfp rho_gcomp.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ib1,ib2,ng,kmax
      integer orhat
Cgetarg       double precision ssite(1),sspec(1),slat(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex cg(ng)
C ... Local parameters
      integer ib,is,iv0,igetss,lmxl,ltop,n1,n2,n3,ng1,nglob,nlm,
     .  nlmtop,nspec,ngabc(3),lfoc,modgkl
Cdelw1 do not change  [['undef', 'wref'], 'ogv', 'rhgcmp']
Cdelw1       integer oqkl,ocs,og,og2,ogv,oiv,okv,osn,oyl,nsp
       integer:: ogv , okv , nsp 
       real(8) ,allocatable :: rv_w_oqkl(:)
       real(8) ,allocatable :: rv_w_ocs(:)
       real(8) ,allocatable :: rv_w_og(:)
       real(8) ,allocatable :: rv_w_og2(:)
       integer ,allocatable :: iv_w_oiv(:)
       real(8) ,allocatable :: rv_w_osn(:)
       real(8) ,allocatable :: rv_w_oyl(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg,
     .  vol,z,q0(3),df(0:20),plat(3,3),qlat(3,3),tau(3)
C ... External calls
      external corprm,defi,defrr,poppr,pshpr,rhgcm2,rhogkl,rlse,stdfac,
     .         suphas,suphs0,suylg,tcn,tcx,upack
C ... Heap
      integer w(1)
      common /w/ w
      data q0 /0d0,0d0,0d0/

      call tcn('rhgcmp')
      call stdfac(20,df)
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: qlat(3,3) , real(8) :: slat%qlat 9
Cgetarg       call upack('lat alat plat qlat',slat,alat,plat,qlat,0,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       i_copy_size=size(slat%qlat) 
       call dcopy(i_copy_size,slat%qlat,1,qlat,1) 

Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           integer :: ng1 , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           integer :: okv , integer(8) :: slat%okv 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat nabc ng ogv okv vol',slat,ngabc,ng1,ogv,okv,vol)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
       ng1=slat%ng
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 115
Cdelw1 do not change ogv because of ['undef', 'wref']
       ogv=slat%ogv
       okv=slat%okv
       vol=slat%vol

      nspec = nglob('nspec')
      nsp   = nglob('nsp')
      modgkl = mode
      if (mode .ge. 200) modgkl = mod(mode,100)
C      if (mode .eq. 0) then
C        modgkl = 131
C      elseif (mode .eq. 1 .or. mode .eq. 3) then
C        modgkl = 31
C      elseif (mode .eq. 2) then
C        modgkl = 1
C      endif
C      call isanrg(mode,0,3,'rhgcmp:','mode', .true.)

C --- Set up help arrays ---
      ltop = 0
      do  is = 1, nspec
Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
         lmxl = int(sspec(is)%lmxl) 

        ltop = max0(ltop,lmxl)
      enddo
      nlmtop = (ltop+1)**2
Cdelw1       call defrr (oyl,  ng*nlmtop)
       allocate(rv_w_oyl(ng*nlmtop))
       if (ng*nlmtop<0) rv_w_oyl(:)=0.0d0

Cdelw1       call defrr (og2,  ng)
       allocate(rv_w_og2(ng))
       if (ng<0) rv_w_og2(:)=0.0d0

Cdelw1       call defrr (og,   ng*3)
       allocate(rv_w_og(ng*3))
       if (ng*3<0) rv_w_og(:)=0.0d0

Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 144
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 144
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call suylg(ltop,alat,ng,w(ogv),w(og),w(og2),w(oyl))
       call suylg ( ltop , alat , ng , w ( ogv ) , rv_w_og , rv_w_og2 
     . , rv_w_oyl ) 

Cdelw1 rlse name= og old_list= oyl og2 og 
Cdelw1 rlse name= og new_list= oyl og2 
Cdelw1       call rlse (og)
       if (allocated(rv_w_og)) deallocate(rv_w_og)


Cdelw1       call defi (oiv, ng*3)
       allocate(iv_w_oiv(ng*3))
       if (ng*3<0) iv_w_oiv(:)=0

Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 148
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 148
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call suphs0(plat,ng,w(ogv),w(oiv))
       call suphs0 ( plat , ng , w ( ogv ) , iv_w_oiv ) 

Cdelw1       call defrr (ocs, ng)
       allocate(rv_w_ocs(ng))
       if (ng<0) rv_w_ocs(:)=0.0d0

Cdelw1       call defrr (osn, ng)
       allocate(rv_w_osn(ng))
       if (ng<0) rv_w_osn(:)=0.0d0


C --- Loop over sites ---
      iv0 = 0
      do  ib = ib1, ib2
Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg         call upack('site spec pos',ssite,ib,is,tau,0,0)
         
         is=ssite(ib)%spec
         i_copy_size=size(ssite(ib)%pos) 
         call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg         call upack('spec lmxl rg',sspec,is,lmxl,rg,0,0)
         
         lmxl=sspec(is)%lmxl
         rg=sspec(is)%rg

        if (lmxl .eq. -1) goto 10
        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        if (mode .eq. 2) cofh = 0
        nlm = (lmxl+1)**2
Cdelw1         call suphas(q0,tau,ng,w(oiv),n1,n2,n3,qlat,w(ocs),w(osn))
         call suphas ( q0 , tau , ng , iv_w_oiv , n1 , n2 , n3 , qlat 
     .   , rv_w_ocs , rv_w_osn ) 


Cdelw1         call defrr(oqkl,nlm*(kmax+1))
         allocate(rv_w_oqkl(nlm*(kmax+1)))
         if (nlm*(kmax+1)<0) rv_w_oqkl(:)=0.0d0

        call pshpr(0)
Cdelw1         call rhogkl(ib,ib,nsp,modgkl,ssite,sspec,orhat,kmax,w(oqkl))
         call rhogkl ( ib , ib , nsp , modgkl , ssite , sspec , orhat 
     .   , kmax , rv_w_oqkl ) 

        call poppr
Cdelw1         call rhgcm2(vol,rg,rfoc,ceh,cofh,kmax,
Cdelw1      .    mod(mode/10,10).ge.2,w(oqkl),nlm,ng,
Cdelw1      .    w(og2),w(oyl),w(ocs),w(osn),cg)
         call rhgcm2 ( vol , rg , rfoc , ceh , cofh , kmax , mod ( mode 
     .   / 10 , 10 ) .ge.2 , rv_w_oqkl , nlm , ng , rv_w_og2 , rv_w_oyl 
     .   , rv_w_ocs , rv_w_osn , cg ) 

Cdelw1         if (mode .ge. 200) call rhgcm3(-z,vol,ng,w(ocs),w(osn),cg)
         if ( mode .ge. 200 ) call rhgcm3 ( - z , vol , ng , rv_w_ocs 
     .   , rv_w_osn , cg ) 

Cdelw1 rlse name= oqkl old_list= oyl og2 oiv ocs osn oqkl 
Cdelw1 rlse name= oqkl new_list= oyl og2 oiv ocs osn 
Cdelw1         call rlse(oqkl)
         if (allocated(rv_w_oqkl)) deallocate(rv_w_oqkl)


        iv0 = iv0+nlm
   10   continue
      enddo

Cdelw1 rlse name= oyl old_list= oyl og2 oiv ocs osn 
Cdelw1 rlse name= oyl new_list= (None)
Cdelw1       call rlse (oyl)
       if (allocated(rv_w_osn)) deallocate(rv_w_osn)
       if (allocated(rv_w_ocs)) deallocate(rv_w_ocs)
       if (allocated(iv_w_oiv)) deallocate(iv_w_oiv)
       if (allocated(rv_w_og2)) deallocate(rv_w_og2)
       if (allocated(rv_w_oyl)) deallocate(rv_w_oyl)

      call tcx('rhgcmp')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ogv]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine rhgcmp 


      subroutine rhgcm2(vol,rg,rfoc,ceh,cofh,kmax,lcor,qkl,nlm,ng,g2,yl,
     .  cs,sn,cg)

C- Convert G_kL expansion of function centered at a site to PW's
C     implicit none
C ... Passed parameters
      integer ng,nlm,kmax
      logical lcor
      double precision ceh,cofh,rfoc,rg,vol,qkl(0:kmax,nlm)
      double precision g2(ng),yl(ng,1),cs(ng),sn(ng)
      double complex cg(ng)
C ... Local parameters
      integer i,ilm,l,ll,lmxl,m,k
      double precision aa,cfoc,cvol,gam,gamf,pi,y0,fac,sqkl
      double complex phase,cc

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (nlm .eq. 0) return
      lmxl = ll(nlm)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      gam = 0.25d0*rg*rg
      gamf = 0.25d0*rfoc*rfoc
      cvol = 1d0/vol
      cfoc = -4d0*pi*y0/vol
      do  i = 1, ng
        phase = dcmplx(cs(i),sn(i))
        aa = dexp(-gam*g2(i))*cvol
        cc = aa*phase*(0d0,1d0)
        ilm = 0
        do  l = 0, lmxl
          cc = cc*(0d0,-1d0)
          do m = -l,l
            ilm = ilm+1
            fac = 1
            sqkl = 0
            do  k = 0, kmax
              sqkl = sqkl + qkl(k,ilm)*fac
              fac = -g2(i)*fac
            enddo
            cg(i) = cg(i) + sqkl*cc*yl(i,ilm)
          enddo
        enddo

        if (lcor) then
          aa = cfoc*dexp(gamf*(ceh-g2(i)))/(ceh-g2(i))
          cg(i) = cg(i) + cofh*aa*phase
        endif

      enddo

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine rhgcm2 

      subroutine rhgcm3(z,vol,ng,cs,sn,cg)

C- PW expansion of Z * delta(r)
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :size of delta-function
Ci   vol   :cell volume
Ci   ng    :number of G-vectors
Ci   cs    :cos(-p*G)
Ci   sn    :cos(-p*G)
Co Outputs
Co   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   26 Oct 01 
C ----------------------------------------------------------------------

C     implicit none
C ... Passed parameters
      integer ng
      double precision z,vol,cs(ng),sn(ng)
      double complex cg(ng)
C ... Local parameters
      integer i
      double complex phase

      do  i = 1, ng
        phase = dcmplx(cs(i),sn(i))
        cg(i) = cg(i) + z*phase/vol
      enddo

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine rhgcm3 


