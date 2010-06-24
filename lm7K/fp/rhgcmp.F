Cchp1 cmdvar= real(8),rv_p_,ogv
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
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex cg(ng)
C ... Local parameters
      integer ib,is,iv0,igetss,lmxl,ltop,n1,n2,n3,ng1,nglob,nlm,
     .nlmtop,nspec,ngabc(3),lfoc,modgkl
Cchp1       integer:: ogv , okv , nsp
       integer:: okv , nsp 
       real(8),pointer :: rv_p_ogv(:)

      real(8) ,allocatable :: qkl_rv(:)
      real(8) ,allocatable :: cs_rv(:)
      real(8) ,allocatable :: g_rv(:)
      real(8) ,allocatable :: g2_rv(:)
      integer ,allocatable :: iv_iv(:)
      real(8) ,allocatable :: sn_rv(:)
      real(8) ,allocatable :: yl_rv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg,
     .vol,z,q0(3),df(0:20),plat(3,3),qlat(3,3),tau(3)
C ... External calls
      external corprm,defi,defrr,poppr,pshpr,rhgcm2,rhogkl,rlse,stdfac,
     .suphas,suphs0,suylg,tcn,tcx,upack
C ... Heap
      integer w(1)
      common /w/ w
      data q0 /0d0,0d0,0d0/

      call tcn('rhgcmp')
      call stdfac(20,df)

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)


      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      ng1=slat%ng
Cchp1       ogv=slat%ogv
       rv_p_ogv => slat%rv_p_ogv 

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
        lmxl = int(sspec(is)%lmxl)

        ltop = max0(ltop,lmxl)
      enddo
      nlmtop = (ltop+1)**2
      allocate(yl_rv(ng*nlmtop))
Ctakao_ZeroClear_NotRequiered        if (ng*nlmtop<0) yl_rv(:)=0.0d0

      allocate(g2_rv(ng))
Ctakao_ZeroClear_NotRequiered        if (ng<0) g2_rv(:)=0.0d0

      allocate(g_rv(ng*3))
Ctakao_ZeroClear_NotRequiered        if (ng*3<0) g_rv(:)=0.0d0

Cchp1       call suylg ( ltop , alat , ng , w ( ogv ) , g_rv , g2_rv 
Cchp1      ., yl_rv )
       call suylg ( ltop , alat , ng , rv_p_ogv , g_rv , g2_rv , yl_rv 
     . ) 


      if (allocated(g_rv)) deallocate(g_rv)


      allocate(iv_iv(ng*3))
Ctakao_ZeroClear_NotRequiered        if (ng*3<0) iv_iv(:)=0

Cchp1       call suphs0 ( plat , ng , w ( ogv ) , iv_iv )
       call suphs0 ( plat , ng , rv_p_ogv , iv_iv ) 


      allocate(cs_rv(ng))
Ctakao_ZeroClear_NotRequiered        if (ng<0) cs_rv(:)=0.0d0

      allocate(sn_rv(ng))
Ctakao_ZeroClear_NotRequiered        if (ng<0) sn_rv(:)=0.0d0


C --- Loop over sites ---
      iv0 = 0
      do  ib = ib1, ib2

        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pos)
        call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1)


        lmxl=sspec(is)%lmxl
        rg=sspec(is)%rg

        if (lmxl .eq. -1) goto 10
        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        if (mode .eq. 2) cofh = 0
        nlm = (lmxl+1)**2
        call suphas ( q0 , tau , ng , iv_iv , n1 , n2 , n3 , qlat 
     .  , cs_rv , sn_rv )


        allocate(qkl_rv(nlm*(kmax+1)))
Ctakao_ZeroClear_NotRequiered          if (nlm*(kmax+1)<0) qkl_rv(:)=0.0d0

        call pshpr(0)
        call rhogkl ( ib , ib , nsp , modgkl , ssite , sspec , orhat 
     .  , kmax , qkl_rv )

        call poppr
        call rhgcm2 ( vol , rg , rfoc , ceh , cofh , kmax , mod ( mode 
     .  / 10 , 10 ) .ge.2 , qkl_rv , nlm , ng , g2_rv , yl_rv 
     .  , cs_rv , sn_rv , cg )

        if ( mode .ge. 200 ) call rhgcm3 ( - z , vol , ng , cs_rv 
     .  , sn_rv , cg )

        if (allocated(qkl_rv)) deallocate(qkl_rv)


        iv0 = iv0+nlm
   10   continue
      enddo

      if (allocated(sn_rv)) deallocate(sn_rv)
      if (allocated(cs_rv)) deallocate(cs_rv)
      if (allocated(iv_iv)) deallocate(iv_iv)
      if (allocated(g2_rv)) deallocate(g2_rv)
      if (allocated(yl_rv)) deallocate(yl_rv)

      call tcx('rhgcmp')
      end subroutine rhgcmp


      subroutine rhgcm2(vol,rg,rfoc,ceh,cofh,kmax,lcor,qkl,nlm,ng,g2,yl,
     .cs,sn,cg)

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

      end subroutine rhgcm3


