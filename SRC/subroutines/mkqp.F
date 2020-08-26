      subroutine mkqp(sctrl,slat,gettet,lnoirr,lreduc,lgstar)
      use m_struc_def  !Cgetarg
      use m_globalvariables
      use m_work,only:iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet, !sbz
     &     rv_a_osymgr
      use m_lmfinit,only: bz_nabc,bz_lshft,bz_nkp,bz_ntet

!! Set up k-points and related quantities for BZ integration
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: lpgf,lmet,lsx
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat npgrp osymgr
Ci     Stored:    *
Ci     Passed to: *
Ci   gettet: T read or generate tetrahedra corners, if
Ci             tetrahedron integration set
Ci   lnoirr: T suppress generation of inequivalent tetrahedra
Ci   lreduc: 0 do not save array ipq
Ci         : 1 save array ipq
Ci         :-1 ignore symmetry operations, make qp for full BZ.
Ci   lgstar: nozero, generate igstar according to bzmesh, which see
Ci         : 0 igstar is not made
Ci         : 2 igstar contains inverse mapping of ipq
Ci         :-2 igstar contains group ops rotating irreducible
Ci         :   to to full BZ.
Cio Inputs/Outputs
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Cio    Elts read: nkabc lshft lopt lio
Cio    Stored:    nkp nkabc oqp owtkp ostar ntet oidtet oipq
!!taken from lmsuit
      implicit none
      logical gettet
      integer lgstar,lreduc,i_copy_size,i_data_size
      type(s_ctrl)::sctrl
      
c      integer :: sbz
c      type(s_bz)::sbz

      type(s_lat)::slat
C ... Local parameters
      logical lgors,ltet,lnoirr,llshft(3),lipq !lsx,
      integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
     .       , lshft(3) , lpbc , ntet , i , iprint , igets , stdo , nglob
      integer,allocatable :: iv_a_owk(:)
      integer,allocatable :: iv_a_tmp(:)
      double precision plat(3,3),qlat(3,3),vol
      character prgnam*4, outs*80
      integer ::iwdummy
      data prgnam/'mkqp'/
C ... Setup
      ntet = 0
      stdo = globalvariables%stdo
      nkxyz=bz_nabc
      lshft=bz_lshft
      plat =slat%plat
      nsgrp=slat%nsgrp
      npgrp=slat%npgrp
      lpbc = 0
      ltet = gettet .and. iand(2,int(sctrl%lmet)) .ne.0
      if(ltet) then
        sctrl%lmet = ior(int(sctrl%lmet),2)               !ltet=T --> lmet2=on (i think)
      else
        sctrl%lmet = sctrl%lmet - iand(int(sctrl%lmet),2) !ltet=F --> lmet2=off (i think)
      endif  

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
      allocate(rv_a_owtkp(abs(-mxkp)))
      rv_a_owtkp(:)=0.0d0
#if TRSYMM
      call dvset ( rv_a_owtkp , 1 , 1 , dble ( nsgrp ) )
#endif
      if (allocated(rv_p_oqp)) deallocate(rv_p_oqp)
      allocate(rv_p_oqp(3*mxkp))
      if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
      allocate(iv_a_oipq(6*mxkp))
      call dinv33(plat,1,qlat,vol)
C   ... Restrict BZ to two dimensions
      if (lpbc .eq. 1) then
        outs = ' ' // prgnam
        if (nkxyz(3) .gt. 1 .and. iprint() .ge. 10) then
          write(stdo,*) ' '
          call awrit2('%a (warning): nk3=%i, shft3=%i; reset to 1,0',
     .     outs,80,-stdo,nkxyz(3),lshft)
        endif
        lshft(3)=0
        nkxyz(3) = 1
        call projql(qlat)
      endif
      do  10  i = 1, 3
        llshft(i) = lshft(i) .ne. 0
 10   continue
      call ivset ( iv_a_ostar , 1 , 1 , lgstar )
      call info0(20,0,0,' ')
      call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz
     . ( 3 ) , llshft , rv_a_osymgr, npgrp , iv_a_oipq , rv_p_oqp , !slat%rv_a_osymgr
     . rv_a_owtkp , nkp , mxkp , iv_a_ostar , lpbc )
C   ... Sanity check
c        if (lsx .and. lshft(1)+lshft(2)+lshft(3) .gt. 0) call
c     .    rx('MKQP:  shifted BZ mesh not allowed with SX')
C   ... Remake, with exactly allocated arrays, ipq on top
      if (allocated(rv_a_owtkp)) deallocate(rv_a_owtkp)
      allocate(rv_a_owtkp(abs(-nkp)))
      if (-nkp<0) rv_a_owtkp(:)=0.0d0
#if TRSYMM
      call dvset ( rv_a_owtkp , 1 , 1 , dble ( nsgrp ) )
#endif
      if (allocated(rv_p_oqp)) deallocate(rv_p_oqp)
      allocate(rv_p_oqp(abs(3*nkp)))
      if (3*nkp<0) rv_p_oqp(:)=0.0d0
      if (allocated(iv_a_oipq)) deallocate(iv_a_oipq)
      allocate(iv_a_oipq(abs(6*mxkp)))
      if (6*mxkp<0) iv_a_oipq(:)=0
      call pshpr(0)
      call ivset ( iv_a_ostar , 1 , 1 , lgstar )
      call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz
     . ( 3 ) , llshft , rv_a_osymgr, npgrp , iv_a_oipq , rv_p_oqp , !slat%
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
      end subroutine mkqp


