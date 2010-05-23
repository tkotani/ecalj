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
      oidtet = 1
      stdo = nglob('stdo')

      i_copy_size=size(sbz%nkabc)
      call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
      i_copy_size=size(sbz%lshft)
      call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)


      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      nsgrp=slat%nsgrp
      npgrp=slat%npgrp
      osymgr=slat%osymgr

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

        call getqp(0,nfilqp,nkp,nkxyz,lshft,ntet,w,w,w)
        if (ltet) then
          if (ntet .eq. 0)
     .      call rx('tetrahedron method specd but no tet weights given')
        else
          ntet = 0
        endif
        call defdr(oqp,3*nkp)
        call defdr(owtkp,-nkp)
        if (ntet .gt. 0) then
          call defi(oidtet,ntet*5)
        else
          call defi(oidtet,1)
        endif
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
          call defi(ogstar,-1)
        else
          call defi(ogstar,-mxkp-1)
        endif
        call defdr(owtkp,-mxkp)
#if TRSYMM
        call dvset(w(owtkp),1,1,dble(nsgrp))
#endif
        call defdr(oqp,3*mxkp)
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
          llshft(i) = lshft(i) .ne. 0
   10   continue
        call ivset(w(ogstar),1,1,lgstar)
        call info0(20,0,0,' ')
        call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,
     .    w(osymgr),npgrp,w(oipq),w(oqp),w(owtkp),nkp,mxkp,w(ogstar),
     .    lpbc)

C   ... Sanity check
c        if (lsx .and. lshft(1)+lshft(2)+lshft(3) .gt. 0) call
c     .    rx('MKQP:  shifted BZ mesh not allowed with SX')

C   ... Remake, with exactly allocated arrays, ipq on top
        call rlse(owtkp)

        call defdr(owtkp,-nkp)
#if TRSYMM
        call dvset(w(owtkp),1,1,dble(nsgrp))
#endif
        call defdr(oqp,3*nkp)
        call defi(oipq,6*mxkp)
        call pshpr(0)
        call ivset(w(ogstar),1,1,lgstar)
        call bzmesh(plat,qlat,nkxyz(1),nkxyz(2),nkxyz(3),llshft,
     .    w(osymgr),npgrp,w(oipq),w(oqp),w(owtkp),nkp,mxkp,w(ogstar),
     .    lpbc)
        call poppr

C   ... Hang on to ipq if (a) call tetirr below or (b) lipq
        if (.not. ((ltet .and. .not. lnoirr) .or. lipq)) then
          call rlse(oipq)

          oipq = -999
        endif

C   ... Generate inequivalent tetrahedra
        if (ltet .and. .not. lnoirr) then
          call defi(oidtet,mxkp*30)
C     ... Save a copy of ipq since tetirr destroys this array
          if (lipq) then
            owk = oipq
            call defi(oipq,6*mxkp)
            call icopy(6*mxkp,w(owk),1,w(oipq),1)
          endif
          call tetirr(qlat,nkxyz(1),nkxyz(2),nkxyz(3),w(oipq),ntet,
     .      w(oidtet))
          call redfi(oidtet,ntet*5)
C     ... Unless ipq still needed, throw it away
          if (lipq) then
            oipq = owk
          else
            call defps2(oipq,oidtet)
            call rlse(oipq)

            oipq = -999
          endif
        endif
      endif

C --- Pack new info into structures ---
C     call prmx('qp',w(oqp),3,3,nkp)

      sbz%nkp=nkp
      i_copy_size=size(sbz%nkabc)
      call ii8copy(i_copy_size,nkxyz,1,sbz%nkabc,1)
      sbz%oqp=oqp
      sbz%owtkp=owtkp
      sbz%ostar=ogstar


      sbz%ntet=ntet
      sbz%oidtet=oidtet

      if ( oipq .gt. 0 ) then

        sbz%oipq=oipq
      endif


C --- Write q-points to disc ---
      if ( iand(2,int(sbz%lio)) .ne.0 ) then

        nfilqp = fopna('QPTS',-1,0)
        call getqp(2,-nfilqp,nkp,nkxyz,lshft,ntet,w(oqp),w(owtkp),
     .    w(oidtet))
        call fclose(nfilqp)
      endif

      end subroutine mkqp


