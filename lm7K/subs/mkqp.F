Cdwdef1 have_common_w_in mkqp , #_of_w_access= 1
Cchp1 cmdvar= real(8),rv_p_,osymgr
Cchp1 cmdvar= integer,iv_p_,ogstar
Cchp1 cmdvar= integer,iv_p_,ostar
Cchp1 cmdvar= integer,iv_p_,oidtet
Cchp1 cmdvar= integer,iv_p_,owk
Cchp1 cmdvar= integer,iv_p_,oipq
Cchp1 cmdvar= real(8),rv_p_,owtkp
Cchp1 cmdvar= real(8),rv_p_,oqp
      subroutine mkqp(sctrl,sbz,slat,gettet,lnoirr,lreduc,lgstar)

      use m_struc_def  !Cgetarg

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
Cchp1       integer fopna,mxkp,nfilqp,nkp,nkxyz(3),nsgrp,npgrp,lshft(3),
Cchp1      .lpbc,ntet,i,oidtet,oqp,owtkp,ogstar,oipq,owk,osymgr,
Cchp1      .iprint,igets,stdo,nglob
Cchp1        integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
Cchp1      . , lshft(3) , lpbc , ntet , i , oidtet , owtkp , ogstar , oipq 
Cchp1      . , owk , osymgr , iprint , igets , stdo , nglob 
Cchp1         integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
Cchp1      .  , lshft(3) , lpbc , ntet , i , oidtet , ogstar , oipq , owk , 
Cchp1      .  osymgr , iprint , igets , stdo , nglob 
Cchp1          integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
Cchp1      .   , lshft(3) , lpbc , ntet , i , oidtet , ogstar , owk , osymgr 
Cchp1      .   , iprint , igets , stdo , nglob 
Cchp1           integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
Cchp1      .    , lshft(3) , lpbc , ntet , i , oidtet , ogstar , osymgr , iprint 
Cchp1      .    , igets , stdo , nglob 
Cchp1            integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
Cchp1      .     , lshft(3) , lpbc , ntet , i , ogstar , osymgr , iprint , igets 
Cchp1      .     , stdo , nglob 
Cchp1             integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
Cchp1      .      , lshft(3) , lpbc , ntet , i , osymgr , iprint , igets , stdo 
Cchp1      .      , nglob 
             integer:: fopna , mxkp , nfilqp , nkp , nkxyz(3) , nsgrp , npgrp 
     .       , lshft(3) , lpbc , ntet , i , iprint , igets , stdo , nglob 
             real(8),pointer :: rv_p_osymgr(:)

            integer,pointer :: iv_p_ogstar(:)

           integer,pointer :: iv_p_oidtet(:)

          integer,pointer :: iv_p_owk(:)

         integer,pointer :: iv_p_oipq(:)

        real(8),pointer :: rv_p_owtkp(:)

       real(8),pointer :: rv_p_oqp(:)

       integer,allocatable :: iv_a_tmp(:)

      double precision plat(3,3),qlat(3,3),vol
      character prgnam*4, outs*80


C ... Heap
Cdwdef1       integer w(1)
Cdwdef1       common /w/ w
       integer ::iwdummy 

      data prgnam/'mkqp'/

C ... Setup
      ntet = 0
cki      oidtet = 1
      nullify(iv_p_oidtet)
      stdo = nglob('stdo')

      i_copy_size=size(sbz%nkabc)
      call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
      i_copy_size=size(sbz%lshft)
      call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)


      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      nsgrp=slat%nsgrp
      npgrp=slat%npgrp
Cchp1       osymgr=slat%osymgr
       rv_p_osymgr => slat%rv_p_osymgr 


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

Cdwdef1         call getqp(0,nfilqp,nkp,nkxyz,lshft,ntet,w,w,w)
         call getqp ( 0 , nfilqp , nkp , nkxyz , lshft , ntet , iwdummy 
     .   , iwdummy , iwdummy ) 

        if (ltet) then
          if (ntet .eq. 0)
     .    call rx('tetrahedron method specd but no tet weights given')
        else
          ntet = 0
        endif
Cchp1         call defdr(oqp,3*nkp)
         allocate(rv_p_oqp(abs(3*nkp)))
         if (3*nkp<0) rv_p_oqp(:)=0.0d0

Cchp1         call defdr(owtkp,-nkp)
         allocate(rv_p_owtkp(abs(-nkp)))
         if (-nkp<0) rv_p_owtkp(:)=0.0d0

        if (ntet .gt. 0) then
Cchp1           call defi(oidtet,ntet*5)
           allocate(iv_p_oidtet(abs(ntet*5)))
           if (ntet*5<0) iv_p_oidtet(:)=0

        else
Cchp1           call defi(oidtet,1)
           allocate(iv_p_oidtet(abs(1)))
           if (1<0) iv_p_oidtet(:)=0

        endif
Cchp1         call getqp(2,nfilqp,nkp,nkxyz,lshft,ntet,w(oqp),w(owtkp),
Cchp1      .  w(oidtet))
Cchp1          call getqp ( 2 , nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
Cchp1      .   , w ( owtkp ) , w ( oidtet ) ) 
Cchp1           call getqp ( 2 , nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
Cchp1      .    , rv_p_owtkp , w ( oidtet ) ) 
           call getqp ( 2 , nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
     .     , rv_p_owtkp , iv_p_oidtet ) 




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
Cchp1           call defi(ogstar,-1)
           allocate(iv_p_ogstar(abs(-1)))
           if (-1<0) iv_p_ogstar(:)=0

        else
Cchp1           call defi(ogstar,-mxkp-1)
           allocate(iv_p_ogstar(abs(-mxkp-1)))
           if (-mxkp-1<0) iv_p_ogstar(:)=0

        endif
Cchp1         call defdr(owtkp,-mxkp)
         allocate(rv_p_owtkp(abs(-mxkp)))
         if (-mxkp<0) rv_p_owtkp(:)=0.0d0

#if TRSYMM
Cchp1         call dvset(w(owtkp),1,1,dble(nsgrp))
         call dvset ( rv_p_owtkp , 1 , 1 , dble ( nsgrp ) ) 

#endif
Cchp1         call defdr(oqp,3*mxkp)
         allocate(rv_p_oqp(abs(3*mxkp)))
         if (3*mxkp<0) rv_p_oqp(:)=0.0d0

Cchp1         call defi(oipq,6*mxkp)
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
Cchp1         call ivset(w(ogstar),1,1,lgstar)
         call ivset ( iv_p_ogstar , 1 , 1 , lgstar ) 

        call info0(20,0,0,' ')
Cchp1         call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,
Cchp1      .  w(osymgr),npgrp,w(oipq),w(oqp),w(owtkp),nkp,mxkp,w(ogstar),
Cchp1      .  lpbc)
Cchp1          call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .   ( 3 ) , llshft , w ( osymgr ) , npgrp , w ( oipq ) , rv_p_oqp 
Cchp1      .   , w ( owtkp ) , nkp , mxkp , w ( ogstar ) , lpbc ) 
Cchp1           call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .    ( 3 ) , llshft , w ( osymgr ) , npgrp , w ( oipq ) , rv_p_oqp 
Cchp1      .    , rv_p_owtkp , nkp , mxkp , w ( ogstar ) , lpbc ) 
Cchp1            call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .     ( 3 ) , llshft , w ( osymgr ) , npgrp , iv_p_oipq , rv_p_oqp 
Cchp1      .     , rv_p_owtkp , nkp , mxkp , w ( ogstar ) , lpbc ) 
Cchp1             call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .      ( 3 ) , llshft , w ( osymgr ) , npgrp , iv_p_oipq , rv_p_oqp 
Cchp1      .      , rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 
             call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
     .       ( 3 ) , llshft , rv_p_osymgr , npgrp , iv_p_oipq , rv_p_oqp , 
     .       rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 






C   ... Sanity check
c        if (lsx .and. lshft(1)+lshft(2)+lshft(3) .gt. 0) call
c     .    rx('MKQP:  shifted BZ mesh not allowed with SX')

C   ... Remake, with exactly allocated arrays, ipq on top
Cchp1         call rlse(owtkp)
         if (associated(rv_p_owtkp)) deallocate(rv_p_owtkp)


Cchp1         call defdr(owtkp,-nkp)
         allocate(rv_p_owtkp(abs(-nkp)))
         if (-nkp<0) rv_p_owtkp(:)=0.0d0

#if TRSYMM
Cchp1         call dvset(w(owtkp),1,1,dble(nsgrp))
         call dvset ( rv_p_owtkp , 1 , 1 , dble ( nsgrp ) ) 

#endif
Cchp1         call defdr(oqp,3*nkp)
         allocate(rv_p_oqp(abs(3*nkp)))
         if (3*nkp<0) rv_p_oqp(:)=0.0d0

Cchp1         call defi(oipq,6*mxkp)
         allocate(iv_p_oipq(abs(6*mxkp)))
         if (6*mxkp<0) iv_p_oipq(:)=0

        call pshpr(0)
Cchp1         call ivset(w(ogstar),1,1,lgstar)
         call ivset ( iv_p_ogstar , 1 , 1 , lgstar ) 

Cchp1         call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,
Cchp1      .  w(osymgr),npgrp,w(oipq),w(oqp),w(owtkp),nkp,mxkp,w(ogstar),
Cchp1      .  lpbc)
Cchp1          call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .   ( 3 ) , llshft , w ( osymgr ) , npgrp , w ( oipq ) , rv_p_oqp 
Cchp1      .   , w ( owtkp ) , nkp , mxkp , w ( ogstar ) , lpbc ) 
Cchp1           call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .    ( 3 ) , llshft , w ( osymgr ) , npgrp , w ( oipq ) , rv_p_oqp 
Cchp1      .    , rv_p_owtkp , nkp , mxkp , w ( ogstar ) , lpbc ) 
Cchp1            call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .     ( 3 ) , llshft , w ( osymgr ) , npgrp , iv_p_oipq , rv_p_oqp 
Cchp1      .     , rv_p_owtkp , nkp , mxkp , w ( ogstar ) , lpbc ) 
Cchp1             call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
Cchp1      .      ( 3 ) , llshft , w ( osymgr ) , npgrp , iv_p_oipq , rv_p_oqp 
Cchp1      .      , rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 
             call bzmesh ( plat , qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz 
     .       ( 3 ) , llshft , rv_p_osymgr , npgrp , iv_p_oipq , rv_p_oqp , 
     .       rv_p_owtkp , nkp , mxkp , iv_p_ogstar , lpbc ) 





        call poppr

C   ... Hang on to ipq if (a) call tetirr below or (b) lipq
        if (.not. ((ltet .and. .not. lnoirr) .or. lipq)) then
Cchp1           call rlse(oipq)
           if (associated(iv_p_oipq)) deallocate(iv_p_oipq)


cki          oipq = -999
             nullify(iv_p_oipq)
        endif

C   ... Generate inequivalent tetrahedra
        if (ltet .and. .not. lnoirr) then
Cchp1           call defi(oidtet,mxkp*30)
           allocate(iv_p_oidtet(abs(mxkp*30)))
           if (mxkp*30<0) iv_p_oidtet(:)=0

C     ... Save a copy of ipq since tetirr destroys this array
          if (lipq) then
cki            owk = oipq
             iv_p_owk => iv_p_oipq
Cchp1             call defi(oipq,6*mxkp)
             allocate(iv_p_oipq(abs(6*mxkp)))
             if (6*mxkp<0) iv_p_oipq(:)=0

Cchp1             call icopy(6*mxkp,w(owk),1,w(oipq),1)
Cchp1              call icopy ( 6 * mxkp , w ( owk ) , 1 , iv_p_oipq , 1 ) 
              call icopy ( 6 * mxkp , iv_p_owk , 1 , iv_p_oipq , 1 ) 


          endif
Cchp1           call tetirr(qlat,nkxyz(1),nkxyz(2),nkxyz(3),w(oipq),ntet,
Cchp1      .    w(oidtet))
Cchp1            call tetirr ( qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz ( 3 ) 
Cchp1      .     , iv_p_oipq , ntet , w ( oidtet ) ) 
            call tetirr ( qlat , nkxyz ( 1 ) , nkxyz ( 2 ) , nkxyz ( 3 ) 
     .      , iv_p_oipq , ntet , iv_p_oidtet ) 


Cchp1           call redfi(oidtet,ntet*5)
           i_data_size=size(iv_p_oidtet); allocate(iv_a_tmp(i_data_size))
           iv_a_tmp=iv_p_oidtet; deallocate(iv_p_oidtet)
           i_data_size=min(i_data_size,ntet*5); allocate(iv_p_oidtet(ntet*5))
           iv_p_oidtet(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)

C     ... Unless ipq still needed, throw it away
          if (lipq) then
cki            oipq = owk
            iv_p_oipq => iv_p_owk
          else
cki            call defps2(oipq,oidtet)
Cchp1             call rlse(oipq)
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
Cchp1       sbz%oqp=oqp
       sbz%rv_p_oqp => rv_p_oqp 

Cchp1       sbz%owtkp=owtkp
       sbz%rv_p_owtkp => rv_p_owtkp 

Cchp1       sbz%ostar=ogstar
       sbz%iv_p_ostar => iv_p_ogstar 



      sbz%ntet=ntet
Cchp1       sbz%oidtet=oidtet
       sbz%iv_p_oidtet => iv_p_oidtet 


cki      if ( oipq .gt. 0 ) then

Cchp1         sbz%oipq=oipq
         sbz%iv_p_oipq => iv_p_oipq 

cki      endif


C --- Write q-points to disc ---
      if ( iand(2,int(sbz%lio)) .ne.0 ) then

        nfilqp = fopna('QPTS',-1,0)
Cchp1         call getqp(2,-nfilqp,nkp,nkxyz,lshft,ntet,w(oqp),w(owtkp),
Cchp1      .  w(oidtet))
Cchp1          call getqp ( 2 , - nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
Cchp1      .   , w ( owtkp ) , w ( oidtet ) ) 
Cchp1           call getqp ( 2 , - nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
Cchp1      .    , rv_p_owtkp , w ( oidtet ) ) 
           call getqp ( 2 , - nfilqp , nkp , nkxyz , lshft , ntet , rv_p_oqp 
     .     , rv_p_owtkp , iv_p_oidtet ) 



        call fclose(nfilqp)
      endif

      end subroutine mkqp


