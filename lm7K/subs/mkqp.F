Cdelw1 var ['oqp', '1', 'oqp', '3*nkp', 'real(8)', 'rv_w_', 'mkqp']
Cdelw1 var ['owtkp', '1', 'owtkp', '-nkp', 'real(8)', 'rv_w_', 'mkqp']
Cdelw1 var ['oidtet', '2', 'oidtet', 'ntet*5', 'integer', 'iv_w_', 'mkqp']
Cdelw1 var ['ogstar', '2', 'ogstar', '-1', 'integer', 'iv_w_', 'mkqp']
Cdelw1 var ['oipq', '1', 'oipq', '6*mxkp', 'integer', 'iv_w_', 'mkqp']
Cdelw1 var ['osymgr', '-1', 'undef', '0', 'undef', 'x', 'mkqp']
Cdelw1 var ['owk', '-1', 'undef', '0', 'undef', 'x', 'mkqp']
Cdelw1 undel [['undef', 'wref'], 'osymgr', 'mkqp']
Cdelw1 undel [['undef', 'wref'], 'owk', 'mkqp']
Cdelw1 undel [['wref'], 'oipq', 'mkqp']
Cdelw1 undel [['wref'], 'oidtet', 'mkqp']
Cdelw1 undel [['wref'], 'oqp', 'mkqp']
Cdelw1 undel [['wref'], 'owtkp', 'mkqp']
Cdelw1 undel [['wref'], 'ogstar', 'mkqp']
Cgetarg...info...           structure ['mkqp', 'bz', 'sbz']
Cgetarg...info...           structure ['mkqp', 'lat', 'slat']
Cgetarg...info...           structure ['mkqp', 'ctrl', 'sctrl']
Cgetarg...info...           use_to_add ['mkqp', 'm_struc_def']
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
Cgetarg       double precision sctrl(1),sbz(1),slat(1)
       type(s_ctrl)::sctrl
       type(s_bz)::sbz
       type(s_lat)::slat

C ... Local parameters
      logical lgors,ltet,lnoirr,llshft(3),lipq !lsx,
Cdelw1 do not change  [['wref'], 'oidtet', 'mkqp']
Cdelw1 do not change  [['wref'], 'oqp', 'mkqp']
Cdelw1 do not change  [['wref'], 'owtkp', 'mkqp']
Cdelw1 do not change  [['wref'], 'ogstar', 'mkqp']
Cdelw1 do not change  [['wref'], 'oipq', 'mkqp']
Cdelw1 do not change  [['undef', 'wref'], 'owk', 'mkqp']
Cdelw1 do not change  [['undef', 'wref'], 'osymgr', 'mkqp']
      integer fopna,mxkp,nfilqp,nkp,nkxyz(3),nsgrp,npgrp,lshft(3),
     .  lpbc,ntet,i,oidtet,oqp,owtkp,ogstar,oipq,owk,osymgr,
     .  iprint,igets,stdo,nglob
      double precision plat(3,3),qlat(3,3),vol
      character prgnam*4, outs*80


C ... Heap
      integer w(1)
      common /w/ w
      data prgnam/'mkqp'/

C ... Setup
      ntet = 0
Cdelw1 warning(1) , probably  oidtet  is not defined yet at linenumber= 70
Cdelw1 do not change oidtet because of ['wref']
      oidtet = 1
      stdo = nglob('stdo')
Cgetarg...info...           integer :: nkxyz(3) , integer(8) :: sbz%nkabc 3
Cgetarg...info...           integer :: lshft(3) , integer(8) :: sbz%lshft 3
Cgetarg       call upack2('bz nkabc lshft',sbz,nkxyz,lshft)
       
       i_copy_size=size(sbz%nkabc) 
       call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1) 
       i_copy_size=size(sbz%lshft) 
       call i8icopy(i_copy_size,sbz%lshft,1,lshft,1) 

Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           integer :: nsgrp , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: npgrp , integer(8) :: slat%npgrp 1
Cgetarg...info...           integer :: osymgr , integer(8) :: slat%osymgr 1
Cgetarg       call upack('lat plat nsgrp npgrp osymgr',slat,plat,nsgrp,npgrp,
Cgetarg      .  osymgr,0)
       
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       nsgrp=slat%nsgrp
       npgrp=slat%npgrp
Cdelw1 warning(1) , probably  osymgr  is not defined yet at linenumber= 92
Cdelw1 do not change osymgr because of ['undef', 'wref']
       osymgr=slat%osymgr

Cgetarg       lpgf = igets('ctrl lpgf',sctrl)
c       lpgf = int(sctrl%lpgf(1)) 

C     lpbc = 0 for kmesh in 3 dimensions, 1 kmesh in 2 dimensions
      lpbc = 0
c      if (lpgf .gt. 0) lpbc = 1
Cgetarg       ltet  = gettet .and. lgors('ctrl lmet,2',sctrl)
       ltet = gettet .and. iand(2,int(sctrl%lmet)) .ne.0 

Cgetarg       call lsets('ctrl lmet',sctrl,ltet,2)
       
       call lsets_bitop_i8(sctrl%lmet,1,ltet,2)

Cgetarg       lsx  = lgors('ctrl lsx,1',sctrl) .or. lgors('ctrl lscr,1',sctrl)
c       lsx = iand(1,int(sctrl%lsx)) .ne.0 .or. iand(1,int(sctrl%lscr)) .ne.0 

      lipq = lreduc .eq. 1
      if (lreduc .eq. -1) npgrp = 1

C ... q-points from BZMAP
Cgetarg       if (lgors('bz lopt,2',sbz)) then
       if ( iand(2,int(sbz%lopt)) .ne.0 ) then 

        call rx('recompile mkqp with BZMAP option')
C ... Read qp from disk
Cgetarg       elseif (lgors('bz lio,1',sbz)) then
       elseif ( iand(1,int(sbz%lio)) .ne.0 ) then 

        call info0(30,0,0,' MKQP:   reading data from file QPTS ...')
        nfilqp = fopna('QPTS',-1,1)

        call getqp(0,nfilqp,nkp,nkxyz,lshft,ntet,w,w,w)
        if (ltet) then
          if (ntet .eq. 0)
     .      call rx('tetrahedron method specd but no tet weights given')
        else
          ntet = 0
        endif
Cdelw1 do not change oqp because of ['wref']
        call defdr(oqp,3*nkp)
Cdelw1 do not change owtkp because of ['wref']
        call defdr(owtkp,-nkp)
        if (ntet .gt. 0) then
Cdelw1 do not change oidtet because of ['wref']
          call defi(oidtet,ntet*5)
        else
Cdelw1 do not change oidtet because of ['wref']
          call defi(oidtet,1)
        endif
Cdelw1 do not change oqp because of ['wref']
Cdelw1 do not change owtkp because of ['wref']
Cdelw1 do not change oidtet because of ['wref']
        call getqp(2,nfilqp,nkp,nkxyz,lshft,ntet,w(oqp),w(owtkp),
     .    w(oidtet))

        call fclose(nfilqp)
        if (iprint() .ge. 20) call
     .    awrit1(' MKQP:   read %i qp from disc',' ',80,stdo,nkp)
        call rxx(ltet,'tet. integration with non-standard k-mesh')
        if (lgstar .ne. 0) then
          call rx('mkqp: lgstar not allowed with user supplied k-mesh')
        endif
        call rx('done')
C ... Make the qp list from bzmesh
      else
        mxkp = nkxyz(1)*nkxyz(2)*nkxyz(3)
        if (lgstar .eq. 0) then
Cdelw1 do not change ogstar because of ['wref']
          call defi(ogstar,-1)
        else
Cdelw1 do not change ogstar because of ['wref']
          call defi(ogstar,-mxkp-1)
        endif
Cdelw1 do not change owtkp because of ['wref']
        call defdr(owtkp,-mxkp)
#if TRSYMM
Cdelw1 do not change owtkp because of ['wref']
        call dvset(w(owtkp),1,1,dble(nsgrp))
#endif
Cdelw1 do not change oqp because of ['wref']
        call defdr(oqp,3*mxkp)
Cdelw1 do not change oipq because of ['wref']
        call defi(oipq,6*mxkp)
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
   10   llshft(i) = lshft(i) .ne. 0
Cdelw1 do not change ogstar because of ['wref']
        call ivset(w(ogstar),1,1,lgstar)
        call info0(20,0,0,' ')
Cdelw1 do not change oqp because of ['wref']
Cdelw1 do not change owtkp because of ['wref']
Cdelw1 do not change ogstar because of ['wref']
Cdelw1 do not change oipq because of ['wref']
Cdelw1 warning(1) , probably  osymgr  is not defined yet at linenumber= 183
Cdelw1 do not change osymgr because of ['undef', 'wref']
        call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,
     .    w(osymgr),npgrp,w(oipq),w(oqp),w(owtkp),nkp,mxkp,w(ogstar),
     .    lpbc)

C   ... Sanity check
c        if (lsx .and. lshft(1)+lshft(2)+lshft(3) .gt. 0) call
c     .    rx('MKQP:  shifted BZ mesh not allowed with SX')

C   ... Remake, with exactly allocated arrays, ipq on top
Cdelw1 rlse name= owtkp old_list= oqp owtkp oidtet ogstar oipq 
Cdelw1 rlse name= owtkp new_list= oqp 
Cdelw1 not deallocate oipq because of [wref]
Cdelw1 not deallocate ogstar because of [wref]
Cdelw1 not deallocate oidtet because of [wref]
Cdelw1 not deallocate owtkp because of [wref]
Cdelw1         call rlse(owtkp)
         call rlse(owtkp)

Cdelw1 do not change owtkp because of ['wref']
        call defdr(owtkp,-nkp)
#if TRSYMM
Cdelw1 do not change owtkp because of ['wref']
        call dvset(w(owtkp),1,1,dble(nsgrp))
#endif
Cdelw1 do not change oqp because of ['wref']
        call defdr(oqp,3*nkp)
Cdelw1 do not change oipq because of ['wref']
        call defi(oipq,6*mxkp)
        call pshpr(0)
Cdelw1 warning(1) , probably  ogstar  is not defined yet at linenumber= 200
Cdelw1 do not change ogstar because of ['wref']
        call ivset(w(ogstar),1,1,lgstar)
Cdelw1 do not change oqp because of ['wref']
Cdelw1 do not change owtkp because of ['wref']
Cdelw1 warning(1) , probably  ogstar  is not defined yet at linenumber= 201
Cdelw1 do not change ogstar because of ['wref']
Cdelw1 do not change oipq because of ['wref']
Cdelw1 warning(1) , probably  osymgr  is not defined yet at linenumber= 201
Cdelw1 do not change osymgr because of ['undef', 'wref']
        call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,
     .    w(osymgr),npgrp,w(oipq),w(oqp),w(owtkp),nkp,mxkp,w(ogstar),
     .    lpbc)
        call poppr

C   ... Hang on to ipq if (a) call tetirr below or (b) lipq
        if (.not. ((ltet .and. .not. lnoirr) .or. lipq)) then
Cdelw1 rlse name= oipq old_list= oqp owtkp oipq 
Cdelw1 rlse name= oipq new_list= oqp owtkp 
Cdelw1 not deallocate oipq because of [wref]
Cdelw1           call rlse(oipq)
           call rlse(oipq)

Cdelw1 warning(1) , probably  oipq  is not defined yet at linenumber= 209
Cdelw1 do not change oipq because of ['wref']
          oipq = -999
        endif

C   ... Generate inequivalent tetrahedra
        if (ltet .and. .not. lnoirr) then
Cdelw1 do not change oidtet because of ['wref']
          call defi(oidtet,mxkp*30)
C     ... Save a copy of ipq since tetirr destroys this array
          if (lipq) then
Cdelw1 warning(1) , probably  oipq  is not defined yet at linenumber= 217
Cdelw1 do not change oipq because of ['wref']
Cdelw1 warning(1) , probably  owk  is not defined yet at linenumber= 217
Cdelw1 do not change owk because of ['undef', 'wref']
            owk = oipq
Cdelw1 do not change oipq because of ['wref']
            call defi(oipq,6*mxkp)
Cdelw1 do not change oipq because of ['wref']
Cdelw1 warning(1) , probably  owk  is not defined yet at linenumber= 219
Cdelw1 do not change owk because of ['undef', 'wref']
            call icopy(6*mxkp,w(owk),1,w(oipq),1)
          endif
Cdelw1 do not change oidtet because of ['wref']
Cdelw1 do not change oipq because of ['wref']
          call tetirr(qlat,nkxyz(1),nkxyz(2),nkxyz(3),w(oipq),ntet,
     .      w(oidtet))
Cdelw1 do not change oidtet because of ['wref']
          call redfi(oidtet,ntet*5)
C     ... Unless ipq still needed, throw it away
          if (lipq) then
Cdelw1 do not change oipq because of ['wref']
Cdelw1 warning(1) , probably  owk  is not defined yet at linenumber= 226
Cdelw1 do not change owk because of ['undef', 'wref']
            oipq = owk
          else
Cdelw1 do not change oidtet because of ['wref']
Cdelw1 do not change oipq because of ['wref']
            call defps2(oipq,oidtet)
Cdelw1 rlse name= oipq old_list= oqp owtkp oidtet oipq 
Cdelw1 rlse name= oipq new_list= oqp owtkp oidtet 
Cdelw1 not deallocate oipq because of [wref]
Cdelw1             call rlse(oipq)
             call rlse(oipq)

Cdelw1 warning(1) , probably  oipq  is not defined yet at linenumber= 230
Cdelw1 do not change oipq because of ['wref']
            oipq = -999
          endif
        endif
      endif

C --- Pack new info into structures ---
C     call prmx('qp',w(oqp),3,3,nkp)
Cgetarg...info...           integer :: nkp , integer(8) :: sbz%nkp 1
Cgetarg...info...           integer :: nkxyz(3) , integer(8) :: sbz%nkabc 3
Cgetarg...info...           integer :: oqp , integer(8) :: sbz%oqp 1
Cgetarg...info...           integer :: owtkp , integer(8) :: sbz%owtkp 1
Cgetarg...info...           integer :: ogstar , integer(8) :: sbz%ostar 1
Cgetarg       call pack5('bz nkp nkabc oqp owtkp ostar',sbz,
Cgetarg      .  nkp,nkxyz,oqp,owtkp,ogstar)
       
       sbz%nkp=nkp 
       i_copy_size=size(sbz%nkabc) 
       call ii8copy(i_copy_size,nkxyz,1,sbz%nkabc,1) 
Cdelw1 do not change oqp because of ['wref']
       sbz%oqp=oqp 
Cdelw1 do not change owtkp because of ['wref']
       sbz%owtkp=owtkp 
Cdelw1 warning(1) , probably  ogstar  is not defined yet at linenumber= 250
Cdelw1 do not change ogstar because of ['wref']
       sbz%ostar=ogstar 

Cgetarg...info...           integer :: ntet , integer(8) :: sbz%ntet 1
Cgetarg...info...           integer :: oidtet , integer(8) :: sbz%oidtet 1
Cgetarg       call pack2('bz ntet oidtet',sbz,ntet,oidtet)
       
       sbz%ntet=ntet 
Cdelw1 do not change oidtet because of ['wref']
       sbz%oidtet=oidtet 

Cgetarg...info...           integer :: oipq , integer(8) :: sbz%oipq 1
Cgetarg       if (oipq .gt. 0) call pack1('bz oipq',sbz,oipq)
Cdelw1 warning(1) , probably  oipq  is not defined yet at linenumber= 261
Cdelw1 do not change oipq because of ['wref']
       if ( oipq .gt. 0 ) then 
       
Cdelw1 warning(1) , probably  oipq  is not defined yet at linenumber= 263
Cdelw1 do not change oipq because of ['wref']
       sbz%oipq=oipq 
       endif 


C --- Write q-points to disc ---
Cgetarg       if (lgors('bz lio,2',sbz)) then
       if ( iand(2,int(sbz%lio)) .ne.0 ) then 

        nfilqp = fopna('QPTS',-1,0)
Cdelw1 do not change oqp because of ['wref']
Cdelw1 do not change owtkp because of ['wref']
Cdelw1 do not change oidtet because of ['wref']
        call getqp(2,-nfilqp,nkp,nkxyz,lshft,ntet,w(oqp),w(owtkp),
     .    w(oidtet))
        call fclose(nfilqp)
      endif

Cgetarg       end
Cdelw1 w_varlist remains: oqp owtkp oidtet 
Cdelw1 w_varlistundel: [osymgr owk oipq oidtet oqp owtkp ogstar]
Cdelw1 w_varlist (undel), remains: [oqp owtkp oidtet]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate oidtet because of [wref]
Cdelw1 not deallocate owtkp because of [wref]
Cdelw1 not deallocate oqp because of [wref]
       end subroutine mkqp 


