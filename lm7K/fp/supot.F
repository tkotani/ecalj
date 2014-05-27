      subroutine supot(mode,sctrl,slat,spot)
      use m_struc_def  
C- Initialization for potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make Madelung matrix for monopoles (ASA)
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: nbas nbasp nspin lcd,4
Ci   slat  :struct containing information about the lattice
Ci     Elts read: alat vol opos awald nkd nkq odlv oqlv plat nabc gmax
Ci                nsgrp osymgr oag
Ci     Stored:    ng ogv okv gmax nabc oips0 obgv
Ci   spot  :struct containing information about the potential
Ci     Stored:    omad osmrho osmpot oorhat
Co Outputs
Cr Remarks
Cr   Permanent arrays are created to handle the potential, depending
Cr   on the mode
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_pot)::spot
C ... Local parameters
      integer:: nbas , nbasp , nbaspp , nsp , nkd , nkq , igets , ngabc(3) 
     .             , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ng , ngrp , iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(3,3),gmax,xx
      integer ::iwdummy
      real(8):: wdummy(3)=0d0
      call info0(50,1,0,' supot : potential setup ... ')
      nbas=sctrl%nbas
      nbasp=sctrl%nbasp
      nsp=sctrl%nspin
      alat=slat%alat
      vol=slat%vol
      awald=slat%awald
      nkd=slat%nkd
      nkq=slat%nkq
      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
        allocate(spot%rv_a_omad(nbaspp**2))
        spot%rv_a_omad(:)=0.0d0
        call madmat ( nbaspp , slat%rv_a_opos , awald , alat , vol , slat%rv_a_odlv
     .      , nkd , slat%rv_a_oqlv , nkq , spot%rv_a_omad )
      endif
C --- Setup for FT charge density, potential representation ---
      if ( iand(4,int(sctrl%lcd)) .ne. 0 ) then
        alat = slat%alat
        plat = slat%plat
        ngabc= slat%nabc
        gmax = slat%gmax
C   ... Generate energy cutoff gmax or n1..n3
c        call pshpr(iprint()-0)
       call mshsiz(alat,plat,0,gmax,ngabc,ngmx)
c       print *,' xxxx mshsiz=',ngmx
        call fftz30(n1,n2,n3,k1,k2,k3)
C   ... Make list of lattice vectors within cutoff
        allocate(slat%rv_a_ogv(abs(ngmx*3)))
        if (allocated(slat%iv_a_okv)) deallocate(slat%iv_a_okv)
        allocate(slat%iv_a_okv(abs(ngmx*3)))
c        call gvlist ( alat , plat , iwdummy , n1 , n2 , n3 , gmax , 8, ngmx , ng , slat%iv_a_okv , slat%rv_a_ogv , xx , xx )
        call gvlst2(alat, plat, wdummy, n1,n2,n3, 0d0,gmax,0,8, ngmx, ng, slat%iv_a_okv, slat%rv_a_ogv, xx, xx)
        if (ng .ne. ngmx) then
          print *,' gmax,ng ngmx=',gmax,ng,ngmx
          call rx('supot: bug in gvlst2')
        endif
        slat%ng = ng
        slat%gmax = gmax
        slat%nabc = ngabc !it was call ii8copy(i_copy_size=3,ngabc,1,slat%nabc,1)
        if (allocated(spot%zv_a_osmrho)) deallocate(spot%zv_a_osmrho)
        allocate(spot%zv_a_osmrho(k1*k2*k3*nsp))
        spot%zv_a_osmrho(:)=0d0
        if (allocated(spot%zv_a_osmpot)) deallocate(spot%zv_a_osmpot)
        allocate(spot%zv_a_osmpot(k1*k2*k3*nsp))
        spot%zv_a_osmpot(:)=0.0d0
        if (allocated(spot%sv_a_oorhat)) deallocate(spot%sv_a_oorhat)
        allocate(spot%sv_a_oorhat(3*nbaspp))
c        do i=1,3*nbaspp
c        enddo
C   ... Setup for symmetrization
        allocate(slat%iv_a_oips0(ng))
        slat%iv_a_oips0(:)=0.0d0
        allocate(slat%zv_a_obgv(ng))
        slat%zv_a_obgv(:)=0.0d0
        ngrp=slat%nsgrp
        call sgvsym ( ngrp , slat%rv_a_osymgr , slat%rv_a_oag , ng , slat%rv_a_ogv , slat%iv_a_oips0
     .       , slat%zv_a_obgv )
      endif
      end subroutine supot
