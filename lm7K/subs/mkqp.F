      subroutine mkqp(sctrl,sbz,slat,gettet,lnoirr,lreduc,lgstar)

      use m_struc_def  !Cgetarg
      use m_globalvariables
C- Set up k-points and related quantities for BZ integration
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
Cl Local variables
Cl   lipq  :T save array ipq
Cr Remarks
Cu  27 Jun 08 Adapt to new getqp.f
Cu  15 Sep 02 Can use sign of wgt to flag which irr points contain
Cu            equivalent points from time-reversal symmetry
Cu  21 Jul 02 Bug fix in second call to bzmesh
Cu   2 Feb 01 revised code to be consistent with comments (lreduc=0,1)
Cr   9 Oct 00 New lreduc, replacing lipq
Cr   6 Jan 98 (MvS) Split lnoirr into lnoirr+lipq options.
Cr  19 Nov 97 (WRL) added lpgf option, projecting qp to 2D
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical gettet
      integer lgstar,lreduc
      type(s_ctrl)::sctrl
      type(s_bz)::sbz
      type(s_lat)::slat

C ... Local parameters
      logical lgors,ltet,lnoirr,llshft(3),lipq !lsx,
             integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
     .       , lshft(3) , lpbc , ntet , i , iprint , igets , stdo , nglob 
ckino 12月.12.2011:               real(8),pointer :: rv_p_osymgr(:) =>NULL()

ckino Dec.8.2011:              integer,pointer :: iv_p_ogstar(:) =>NULL()

ckino Dec.8.2011:             integer,pointer :: iv_p_oidtet(:) =>NULL()

ckino Dec.14 manual correction          integer,pointer :: iv_p_owk(:) =>NULL()
          integer,allocatable :: iv_a_owk(:) 

         integer,pointer :: iv_p_oipq(:) =>NULL()

ckino Dec.8.2011:          real(8),pointer :: rv_p_owtkp(:) =>NULL()

ckino Dec.8.2011:         real(8),pointer :: rv_p_oqp(:) =>NULL()

       integer,allocatable :: iv_a_tmp(:)

      double precision plat(3,3),qlat(3,3),vol
      character prgnam*4, outs*80


C ... Heap
       integer ::iwdummy 

      data prgnam/'mkqp'/

C ... Setup
      ntet = 0
cki      oidtet = 1
ckino Dec.8.2011:        nullify(iv_p_oidtet)
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo

      i_copy_size=size(sbz%nkabc)
      call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
      i_copy_size=size(sbz%lshft)
      call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)


      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      nsgrp=slat%nsgrp
      npgrp=slat%npgrp
ckino 12月.12.2011:         rv_p_osymgr => slat%rv_p_osymgr 


c       lpgf = int(sctrl%lpgf(1))

C     lpbc = 0 for kmesh in 3 dimensions, 1 kmesh in 2 dimensions
      lpbc = 0
c      if (lpgf .gt. 0) lpbc = 1
      ltet = gettet .and. iand(2,int(sctrl%lmet)) .ne.0


      call lsets_bitop_i8(sctrl%lmet,1,ltet,2)

c       lsx = iand(1,int(sctrl%lsx)) .ne.0 .or. iand(1,int(sctrl%lscr)) .ne.0

      lipq = lreduc .eq. 1
      if (lreduc .eq. -1) npgrp = 1

C ... q-points from BZMAP
      if ( iand(2,int(sbz%lopt)) .ne.0 ) then

        call rx('recompile mkqp with BZMAP option')
C ... Read qp from disk
      elseif ( iand(1,int(sbz%lio)) .ne.0 ) then

        call info0(30,0,0,' MKQP:   reading data from file QPTS ...')
        nfilqp = fopna('QPTS',-1,1)

         call getqp ( 0 , nfilqp , nkp , nkxyz , lshft , ntet , iwdummy 
     .   , iwdummy , iwdummy ) 

        if (ltet) then
          if (ntet .eq. 0)
     .    call rx('tetrahedron method specd but no tet weights given')
        else
          ntet = 0
        endif
ckino Dec.8.2011:           allocate(rv_p_oqp(abs(3*nkp)))
         allocate(sbz%rv_p_oqp(abs(3*nkp)))
ckino Dec.8.2011:           if (3*nkp<0) rv_p_oqp(:)=0.0d0
         if (3*nkp<0) sbz%rv_p_oqp(:)=0.0d0

ckino Dec.8.2011:           allocate(rv_p_owtkp(abs(-nkp)))
         allocate(sbz%rv_p_owtkp(abs(-nkp)))
ckino Dec.8.2011:           if (-nkp<0) rv_p_owtkp(:)=0.0d0
         if (-nkp<0) sbz%rv_p_owtkp(:)=0.0d0

        if (ntet .gt. 0) then
ckino Dec.8.2011:             allocate(iv_p_oidtet(abs(ntet*5)))
           allocate(sbz%iv_p_oidtet(abs(ntet*5)))
ckino Dec.8.2011:             if (ntet*5<0) iv_p_oidtet(:)=0
           if (ntet*5<0) sbz%iv_p_oidtet(:)=0

        else
ckino Dec.8.2011:             allocate(iv_p_oidtet(abs(1)))
           allocate(sbz%iv_p_oidtet(abs(1)))
ckino Dec.8.2011:             if (1<0) iv_p_oidtet(:)=0
           if (1<0) sbz%iv_p_oidtet(:)=0

        endif
ckino Dec.8.2011:             call getqp ( 2 , nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
           call getqp ( 2 , nfilqp , nkp , nkxyz , lshft , ntet , sbz%rv_p_oqp 
ckino Dec.8.2011:       .     , rv_p_owtkp , iv_p_oidtet ) 
ckino Dec.8.2011:       .     , sbz%rv_p_owtkp , iv_p_oidtet ) 
     .     , sbz%rv_p_owtkp , sbz%iv_p_oidtet ) 




        call fclose(nfilqp)
        if (iprint() .ge. 20) call
     .  awrit1(' MKQP:   read %i qp from disc',' ',80,stdo,nkp)
        call rxx(ltet,'tet. integration with non-standard k-mesh')
        if (lgstar .ne. 0) then
          call rx('mkqp: lgstar not allowed with user supplied k-mesh')
        endif
        call rx('done')
C ... Make the qp list from bzmesh
      else
        mxkp = nkxyz(1)*nkxyz(2)*nkxyz(3)
        if (lgstar .eq. 0) then
ckino Dec.8.2011:             allocate(iv_p_ogstar(abs(-1)))
           allocate(sbz%iv_p_ostar(abs(-1)))
ckino Dec.8.2011:             if (-1<0) iv_p_ogstar(:)=0
           if (-1<0) sbz%iv_p_ostar(:)=0

        else
ckino Dec.8.2011:             allocate(iv_p_ogstar(abs(-mxkp-1)))
           allocate(sbz%iv_p_ostar(abs(-mxkp-1)))
ckino Dec.8.2011:             if (-mxkp-1<0) iv_p_ogstar(:)=0
           if (-mxkp-1<0) sbz%iv_p_ostar(:)=0

        endif
ckino Dec.8.2011:           allocate(rv_p_owtkp(abs(-mxkp)))
         allocate(sbz%rv_p_owtkp(abs(-mxkp)))
ckino Dec.8.2011:           if (-mxkp<0) rv_p_owtkp(:)=0.0d0
         if (-mxkp<0) sbz%rv_p_owtkp(:)=0.0d0

#if TRSYMM
ckino Dec.8.2011:           call dvset ( rv_p_owtkp , 1 , 1 , dble ( nsgrp ) ) 
         call dvset ( sbz%rv_p_owtkp , 1 , 1 , dble ( nsgrp ) ) 

#endif
ckino Dec.8.2011:           allocate(rv_p_oqp(abs(3*mxkp)))
         allocate(sbz%rv_p_oqp(abs(3*mxkp)))
ckino Dec.8.2011:           if (3*mxkp<0) rv_p_oqp(:)=0.0d0
         if (3*mxkp<0) sbz%rv_p_oqp(:)=0.0d0

         allocate(iv_p_oipq(abs(6*mxkp)))
         if (6*mxkp<0) iv_p_oipq(:)=0

        call dinv33(plat,1,qlat,vol)

C   ... Restrict BZ to two dimensions
        if (lpbc .eq. 1) then
          outs = ' ' // prgnam
          if (nkxyz(3) .gt. 1 .and. iprint() .ge. 10) then
            write(stdo,*) ' '
            call awrit2('%a (warning): nk3=%i, shft3=%i; reset to 1,0',
     .      outs,80,-stdo,nkxyz(3),lshft)
          endif
          lshft(3)=0
          nkxyz(3) = 1
          call projql(qlat)
        endif

        do  10  i = 1, 3
          llshft(i) = lshft(i) .ne. 0
   10   continue
ckino Dec.8.2011:           call ivset ( iv_p_ogstar , 1 , 1 , lgstar ) 
         call ivset ( sbz%iv_p_ostar , 1 , 1 , lgstar ) 

        call info0(20,0,0,' ')
             call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
ckino Dec.8.2011:       .       ( 3 ) , llshft , rv_p_osymgr , npgrp , iv_p_oipq , rv_p_oqp , 
ckino 12月.12.2011:       .       ( 3 ) , llshft , rv_p_osymgr , npgrp , iv_p_oipq , sbz%rv_p_oqp , 
     .       ( 3 ) , llshft , slat%rv_p_osymgr, npgrp , iv_p_oipq , sbz%rv_p_oqp , 
ckino Dec.8.2011:       .       rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 
ckino Dec.8.2011:       .       sbz%rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 
     .       sbz%rv_p_owtkp , nkp , mxkp , sbz%iv_p_ostar , lpbc ) 






C   ... Sanity check
c        if (lsx .and. lshft(1)+lshft(2)+lshft(3) .gt. 0) call
c     .    rx('MKQP:  shifted BZ mesh not allowed with SX')

C   ... Remake, with exactly allocated arrays, ipq on top
ckino Dec.8.2011:           if (associated(rv_p_owtkp)) deallocate(rv_p_owtkp)
         if (associated(sbz%rv_p_owtkp)) deallocate(sbz%rv_p_owtkp)


ckino Dec.8.2011:           allocate(rv_p_owtkp(abs(-nkp)))
         allocate(sbz%rv_p_owtkp(abs(-nkp)))
ckino Dec.8.2011:           if (-nkp<0) rv_p_owtkp(:)=0.0d0
         if (-nkp<0) sbz%rv_p_owtkp(:)=0.0d0

#if TRSYMM
ckino Dec.8.2011:           call dvset ( rv_p_owtkp , 1 , 1 , dble ( nsgrp ) ) 
         call dvset ( sbz%rv_p_owtkp , 1 , 1 , dble ( nsgrp ) ) 

#endif
ckino Dec.8.2011:           allocate(rv_p_oqp(abs(3*nkp)))
         allocate(sbz%rv_p_oqp(abs(3*nkp)))
ckino Dec.8.2011:           if (3*nkp<0) rv_p_oqp(:)=0.0d0
         if (3*nkp<0) sbz%rv_p_oqp(:)=0.0d0

         allocate(iv_p_oipq(abs(6*mxkp)))
         if (6*mxkp<0) iv_p_oipq(:)=0

        call pshpr(0)
ckino Dec.8.2011:           call ivset ( iv_p_ogstar , 1 , 1 , lgstar ) 
         call ivset ( sbz%iv_p_ostar , 1 , 1 , lgstar ) 

             call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
ckino Dec.8.2011:       .       ( 3 ) , llshft , rv_p_osymgr , npgrp , iv_p_oipq , rv_p_oqp , 
ckino 12月.12.2011:       .       ( 3 ) , llshft , rv_p_osymgr , npgrp , iv_p_oipq , sbz%rv_p_oqp , 
     .       ( 3 ) , llshft , slat%rv_p_osymgr, npgrp , iv_p_oipq , sbz%rv_p_oqp , 
ckino Dec.8.2011:       .       rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 
ckino Dec.8.2011:       .       sbz%rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 
     .       sbz%rv_p_owtkp , nkp , mxkp , sbz%iv_p_ostar , lpbc ) 





        call poppr

C   ... Hang on to ipq if (a) call tetirr below or (b) lipq
        if (.not. ((ltet .and. .not. lnoirr) .or. lipq)) then
           if (associated(iv_p_oipq)) deallocate(iv_p_oipq)


cki          oipq = -999
             nullify(iv_p_oipq)
        endif

C   ... Generate inequivalent tetrahedra
        if (ltet .and. .not. lnoirr) then
ckino Dec.8.2011:             allocate(iv_p_oidtet(abs(mxkp*30)))
           allocate(sbz%iv_p_oidtet(abs(mxkp*30)))
ckino Dec.8.2011:             if (mxkp*30<0) iv_p_oidtet(:)=0
           if (mxkp*30<0) sbz%iv_p_oidtet(:)=0

C     ... Save a copy of ipq since tetirr destroys this array
          if (lipq) then
cki            owk = oipq
ckino Dec.14.2011 manual correction     iv_p_owk => iv_p_oipq
             allocate(iv_a_owk(abs(6*mxkp))); iv_a_owk(:)=iv_p_oipq(:)
ckino Dec.14.2011 manual correction             allocate(iv_p_oipq(abs(6*mxkp)))
             
ckino Dec.14.2011 manual correction             if (6*mxkp<0) iv_p_oipq(:)=0
              iv_p_oipq(:)=0

ckino Dec.14.2011 manual correction              call icopy ( 6 * mxkp , iv_p_owk , 1 , iv_p_oipq , 1 ) 
              call icopy ( 6 * mxkp , iv_a_owk , 1 , iv_p_oipq , 1 ) 


          endif
            call tetirr ( qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz ( 3 ) 
ckino Dec.8.2011:       .      , iv_p_oipq , ntet , iv_p_oidtet ) 
     .      , iv_p_oipq , ntet , sbz%iv_p_oidtet ) 


ckino Dec.8.2011:             i_data_size=size(iv_p_oidtet); allocate(iv_a_tmp(i_data_size))
           i_data_size=size(sbz%iv_p_oidtet); allocate(iv_a_tmp(i_data_size))
ckino Dec.8.2011:             iv_a_tmp=iv_p_oidtet; deallocate(iv_p_oidtet)
           iv_a_tmp=sbz%iv_p_oidtet; deallocate(sbz%iv_p_oidtet)
ckino Dec.8.2011:             i_data_size=min(i_data_size,ntet*5); allocate(iv_p_oidtet(ntet*5))
           i_data_size=min(i_data_size,ntet*5); allocate(sbz%iv_p_oidtet(ntet*5))
ckino Dec.8.2011:             iv_p_oidtet(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)
           sbz%iv_p_oidtet(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)

C     ... Unless ipq still needed, throw it away
          if (lipq) then
cki            oipq = owk
ckino Dec.14.2011 manual correction            iv_p_oipq => iv_p_owk
               iv_p_oipq(:)=iv_a_owk(:); deallocate(iv_a_owk) 
          else
cki            call defps2(oipq,oidtet)
             if (associated(iv_p_oipq)) deallocate(iv_p_oipq)


cki            oipq = -999
              nullify(iv_p_oipq)
          endif
        endif
      endif

C --- Pack new info into structures ---
C     call prmx('qp',w(oqp),3,3,nkp)

      sbz%nkp=nkp
      i_copy_size=size(sbz%nkabc)
      call ii8copy(i_copy_size,nkxyz,1,sbz%nkabc,1)
ckino Dec.8.2011:         sbz%rv_p_oqp => rv_p_oqp 

ckino Dec.8.2011:         sbz%rv_p_owtkp => rv_p_owtkp 

ckino Dec.8.2011:         sbz%iv_p_ostar => iv_p_ogstar 



      sbz%ntet=ntet
ckino Dec.8.2011:         if (.not.associated(iv_p_oidtet))then
       if (.not.associated(sbz%iv_p_oidtet))then
ckino Dec.8.2011:            write(*,*)'warnig in mkqp: iv_p_oidtet is not allocated. but continue.'
          write(*,*)'warnig in mkqp: sbz%iv_p_oidtet is not allocated. but continue.'
       endif
ckino Dec.8.2011:         sbz%iv_p_oidtet => iv_p_oidtet 


cki      if ( oipq .gt. 0 ) then

         sbz%iv_p_oipq => iv_p_oipq 

cki      endif


C --- Write q-points to disc ---
      if ( iand(2,int(sbz%lio)) .ne.0 ) then

        nfilqp = fopna('QPTS',-1,0)
ckino Dec.8.2011:             call getqp ( 2 , - nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
           call getqp ( 2 , - nfilqp , nkp , nkxyz , lshft , ntet , sbz%rv_p_oqp 
ckino Dec.8.2011:       .     , rv_p_owtkp , iv_p_oidtet ) 
ckino Dec.8.2011:       .     , sbz%rv_p_owtkp , iv_p_oidtet ) 
     .     , sbz%rv_p_owtkp , sbz%iv_p_oidtet ) 



        call fclose(nfilqp)
      endif

      end subroutine mkqp


