      module m_mkqp
      integer,allocatable,protected :: iv_a_oidtet(:), iv_a_oipq(:), iv_a_ostar(:)
      real(8),allocatable,protected :: rv_p_oqp (:), rv_a_owtkp(:)
      contains
      subroutine m_mkqp_init() !lnoirr,lreduc,lgstar)
      use m_struc_def 
      use m_mksym,only:    rv_a_osymgr,lat_nsgrp,lat_npgrp
      use m_lmfinit,only: bz_nabc,bz_lshft,bz_nkp,bz_ntet,ctrl_lmet,bz_lmet,ctrl_ldos
      use m_lattic,only: lat_plat
!! Set up k-points and related quantities for BZ integration
!! ----------------------------------------------------------------------
!! gettet: = bz_lmet/=0 .or. ctrl_ldos/=0 given at setup
!!            T read or generate tetrahedra corners, if
!!             tetrahedron integration set
!!   lnoirr: =F (T suppress generation of inequivalent tetrahedra)
!!   lreduc: xxx 0 do not save array ipq
!!         : =1 save array ipq
!!         : xxx -1 ignore symmetry operations, make qp for full BZ.
!!   lgstar: nozero, generate igstar according to bzmesh, which see
!!         : xxx 0 igstar is not made
!!         : xxx 2 igstar contains inverse mapping of ipq
!!         : =-2 igstar contains group ops rotating irreducible
!!         :    to to full BZ.
!! Outputs protected
!! idtet, iqp, qp, star, wtkp : for the Brillouin Zone
      implicit none
      logical gettet
      integer::  lgstar=-2,lreduc=1, i_copy_size,i_data_size
      logical:: lgors,ltet,lnoirr=.false.,llshft(3),lipq !lsx,
      integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
     .       , lshft(3) , lpbc , ntet , i , iprint , igets 
      integer,allocatable :: iv_a_owk(:)
      integer,allocatable :: iv_a_tmp(:)
      double precision plat(3,3),qlat(3,3),vol
      character outs*80
      integer ::iwdummy
c      data prgnam/'mkqp'/
C ... Setup
      gettet = bz_lmet/=0 .or. ctrl_ldos/=0
      ntet = 0
c      stdo = globalvariables%stdo
      nkxyz=bz_nabc
      lshft=bz_lshft
      plat =lat_plat
      nsgrp=lat_nsgrp
      npgrp=lat_npgrp
      lpbc = 0
      ltet = gettet .and. iand(2,int(ctrl_lmet)) .ne.0
      lipq = lreduc .eq. 1
      if (lreduc .eq. -1) npgrp = 1
C not now ... q-points from BZMAP, C ... Read qp from disk
c      if ( iand(2,int(sbz%lopt)) .ne.0 ) call rx('recompile mkqp with BZMAP option')
c      if ( iand(1,int(sbz%lio)) .ne.0 ) call rx('not support reading QPTS')
C ... Make the qp list from bzmesh
      mxkp = nkxyz(1)*nkxyz(2)*nkxyz(3)
      if (lgstar .eq. 0) then
        allocate(iv_a_ostar(1))
        iv_a_ostar(:)=0
      else
        allocate(iv_a_ostar(mxkp+1))
        iv_a_ostar(:)=0
      endif
      allocate(rv_a_owtkp(mxkp))
      rv_a_owtkp(:)=0.0d0
c#if TRSYMM
c      call dvset ( rv_a_owtkp , 1 , 1 , dble ( nsgrp ) )
c#endif
      if (allocated(rv_p_oqp)) deallocate(rv_p_oqp)
      allocate(rv_p_oqp(3*mxkp))
      if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
      allocate(iv_a_oipq(6*mxkp))
      call dinv33(plat,1,qlat,vol)
c$$$C   ... Restrict BZ to two dimensions
c$$$      if (lpbc .eq. 1) then
c$$$        outs = ' ' // prgnam
c$$$        if (nkxyz(3) .gt. 1 .and. iprint() .ge. 10) then
c$$$          write(stdo,*) ' '
c$$$          call awrit2('%a (warning): nk3=%i, shft3=%i; reset to 1,0',
c$$$     .     outs,80,-stdo,nkxyz(3),lshft)
c$$$        endif
c$$$        lshft(3)=0
c$$$        nkxyz(3) = 1
c$$$        call projql(qlat)
c$$$      endif
      do  10  i = 1, 3
        llshft(i) = lshft(i) .ne. 0
 10   continue
      call ivset ( iv_a_ostar , 1 , 1 , lgstar )
      call info0(20,0,0,' ')
      call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz
     . ( 3 ) , llshft , rv_a_osymgr, npgrp , iv_a_oipq , rv_p_oqp , !lat_rv_a_osymgr
     . rv_a_owtkp , nkp , mxkp , iv_a_ostar , lpbc )
C   ... Sanity check
c        if (lsx .and. lshft(1)+lshft(2)+lshft(3) .gt. 0) call
c     .    rx('MKQP:  shifted BZ mesh not allowed with SX')
C   ... Remake, with exactly allocated arrays, ipq on top
      if (allocated(rv_a_owtkp)) deallocate(rv_a_owtkp)
      allocate(rv_a_owtkp(nkp))
      rv_a_owtkp(:)=0.0d0
c#if TRSYMM
c      call dvset ( rv_a_owtkp , 1 , 1 , dble ( nsgrp ) )
c#endif
      if (allocated(rv_p_oqp)) deallocate(rv_p_oqp)
      allocate(rv_p_oqp(abs(3*nkp)))
      if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
      allocate(iv_a_oipq(abs(6*mxkp)))
      call pshpr(0)
      call ivset ( iv_a_ostar , 1 , 1 , lgstar )
      call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz
     . ( 3 ) , llshft , rv_a_osymgr, npgrp , iv_a_oipq , rv_p_oqp , !lat_
     . rv_a_owtkp , nkp , mxkp , iv_a_ostar , lpbc )
      call poppr
C   ... Hang on to ipq if (a) call tetirr below or (b) lipq
      if (.not. ((ltet .and. .not. lnoirr) .or. lipq)) then
        if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
        if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
        if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
      endif
C   ... Generate inequivalent tetrahedra
      if (ltet .and. .not. lnoirr) then
        allocate(iv_a_oidtet(abs(mxkp*30)))
        if (mxkp*30<0) iv_a_oidtet(:)=0
C     ... Save a copy of ipq since tetirr destroys this array
        if (lipq) then
          allocate(iv_a_owk(abs(6*mxkp))); iv_a_owk(:)=iv_a_oipq(:)
          iv_a_oipq(:)=0
          call icopy ( 6 * mxkp , iv_a_owk , 1 , iv_a_oipq , 1 )
        endif
        call tetirr ( qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz ( 3 )
     .   , iv_a_oipq , ntet , iv_a_oidtet )

        i_data_size=size(iv_a_oidtet); allocate(iv_a_tmp(i_data_size))
        iv_a_tmp=iv_a_oidtet; deallocate(iv_a_oidtet)
        i_data_size=min(i_data_size,ntet*5); allocate(iv_a_oidtet(ntet*5))
        iv_a_oidtet(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)

C     ... Unless ipq still needed, throw it away
        if (lipq) then
          iv_a_oipq(:)=iv_a_owk(:); deallocate(iv_a_owk)
        else
          if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
          if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
          if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
        endif
      endif
C --- Pack new info into structures ---
      bz_nkp   = nkp
      bz_nabc = nkxyz
      bz_ntet  = ntet
      if (.not.allocated(iv_a_oidtet))then
        write(*,*)'warnig in mkqp: iv_a_oidtet is not allocated. but continue.'
      endif
c$$$C --- Write q-points to disc ---
c$$$      if ( iand(2,int(sbz%lio)) .ne.0 ) then
c$$$        nfilqp = fopna('QPTS',-1,0)
c$$$        call getqp ( 2 , - nfilqp , nkp , nkxyz , lshft , ntet , sbz%rv_p_oqp
c$$$     .     , sbz%rv_a_owtkp , sbz%iv_a_oidtet )
c$$$        call fclose(nfilqp)
c$$$      endif
      end subroutine m_mkqp_init
      end module m_mkqp

