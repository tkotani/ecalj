#define USEOPTIONAL 1
      module m_seneinput
      implicit none
      integer(4),allocatable:: npair(:,:), nlat(:,:,:,:)
      integer(4)::npairmx
      real(8),allocatable:: qsmesh(:,:,:,:)
      real(8),allocatable:: qwgt(:,:,:)
      end module

!! == Give \Sigma-Vxc for given qp from sfz (this does not work well)==
      subroutine seneinterp2(qp,nbas,nk1,nk2,nk3,plat, !npair,nlat,qwgt,npairmx,qsmesh,
     i  sfz,ldim,ndimh, 
     o  sene)
      use m_hamindex,only: napwmx, igv2,napwk,igv2rev,ibasindex,getikt
      use m_seneinput,only: npair,nlat,qwgt,npairmx,qsmesh !they are read only
      implicit none
      integer(4):: nbas,nk1,nk2,nk3,ib1,ib2,ii,ik1,ik2,ik3,ldim,ndimh,
     &  ikc(3,0:1,0:1,0:1), igshift(3,0:1,0:1,0:1),
     &  j1x,j2x,ix1,ix2,ix3,ikt,ig,nkk(3),igvecc0(3),iktm,igvecc(3),j1,j2
c     &  npair(nbas,nbas),nlat(3,npairmx,nbas,nbas),
      complex(8):: sfz(nk1,nk2,nk3,ndimh,ndimh), 
     &  pi2i,sene(ndimh,ndimh),sss
      complex(8),allocatable:: wgtq(:,:,:,:,:)
      real(8)::rrr(3),plat(3,3) !,qsmesh(3,nk1,nk2,nk3)
      logical:: oncewrite,debug=.false.
      integer(4),allocatable::igx(:,:,:,:)
      real(8):: wgc(0:1,0:1,0:1),qpara(3),qmm(3),qqq(3),www,qp(3) !,qwgt(npairmx,nbas,nbas)

      if(debug) print *,' seneinterp2: qp=',qp,nk1,nk2,nk3
      sene=0d0
! --- Generate interpolation weight wgtq.
      allocate( wgtq(nbas,nbas,nk1,nk2,nk3) )
      pi2i = 2d0*4d0*datan(1d0)*(0d0,1d0)
      wgtq = 0d0
      do ib1=1,nbas
        do ib2=1,nbas
          do ii = 1,npair(ib1,ib2)
            rrr = matmul(plat, nlat(:,ii,ib1,ib2))
! wgtq is interpolation weight for MTO.
            do ik1=1,nk1
              do ik2=1,nk2
                do ik3=1,nk3
                  wgtq(ib1,ib2,ik1,ik2,ik3) = wgtq(ib1,ib2,ik1,ik2,ik3) 
     &      + qwgt(ii,ib1,ib2)*exp( -pi2i*sum( (qp-qsmesh(:,ik1,ik2,ik3))*rrr ) )
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
! ---  MTO part of interpolated sigm(i,j)
      do j1=1,ldim
        do j2=1,ldim
          ib1= ibasindex(j1)
          ib2= ibasindex(j2)
          do ik1=1,nk1
            do ik2=1,nk2
              do ik3=1,nk3
                sene(j1,j2)= sene(j1,j2) + wgtq(ib1,ib2,ik1,ik2,ik3)*sfz(ik1,ik2,ik3,j1,j2) !,isp)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wgtq)
!!--- return when no APW
      if(ldim==ndimh) then
        return
      endif

c      ikt=getikt(qp) ! index for qp
c      if(napwk(ikt)==0) then
c      if(oncewrite(4)) print *,' napw=0 skip napw part '
c        return
c      endif


!!--- APW part of interpolarion sigm(i,j)
!!... range check qqq is in [ikc1/nk1, (ikc1+1)/nk1],... qpara is [0,1],[0,1],[0,1]
      print *,' --- goto qparam--- '
      call qparam(qp,qsmesh,nk1,nk2,nk3,plat, ikc,qpara,igshift)
      ! eight corners of a cube -->
      !  do ix1=0,1
      !  do ix2=0,1
      !  do ix3=0,1
      !    qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
      !  enddo
      !  enddo
      !  enddo
!!... interpolation weight at eight corners of a cube
      call interparam(qpara, wgc)
      print *,' --- end of interparam ---'
!!... obtain G index for eight corners of a cube. qp is in the cube.
      allocate( igx(ndimh,0:1,0:1,0:1) )

c      ikt=getikt(qp) ! index for qp
c      print *,'ikt qp=',ikt,qp
      ix1=0
      ix2=0
      ix3=0
      qmm = qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
      ikt=getikt(qmm) ! index for qp
      print *, 'test xxxxxxxxx333333333333 qp  ',qp
      print *, 'test xxxxxxxxx333333333333 qmm ',qmm
      print *, 'test xxxxxxxxx333333333333 qp-qmm ',qp-qmm
      print *, sum(plat(:,1)*(qp-qmm))
      print *, sum(plat(:,2)*(qp-qmm))
      print *, sum(plat(:,3)*(qp-qmm))
!!    
      nkk= (/nk1,nk2,nk3/)
      do ig=1,napwk(ikt) !qp+G
        do ix1=0,1
          do ix2=0,1
            do ix3=0,1
              igvecc0 = igv2(:,ig,ikt) + igshift(:,ix1,ix2,ix3)
              qmm = qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
              iktm = getikt(qmm) !index for qmm
              qqq = matmul(qp-qmm,plat)
              if(ig==1) then
                write(6,"('qmm qqq=', 13f13.4)")qmm,qqq
              endif  
              igvecc = igvecc0 + qqq + 10000 + 1d-8 ! +10000 is to make igvecc+qqq positve. and then subtract -10000
              igvecc = igvecc - 10000              !
              igx(ig,ix1,ix2,ix3) = igv2rev(igvecc(1),igvecc(2),igvecc(3),iktm)
              !igx at eight corners. If missing, 9999999 !see suham.F
              ! if(ix1==0.and.ix2==0.and.ix3==0) then
              !   write(6,"(a,3f8.3)")'=', igvecc0 + qqq- igvecc
              !   print *,'iktm qmm=',iktm,qmm
              !   write(6,"(a,20i8)")'ix1 ix2 ix3 igvecc=', igvecc,igshift(:,ix1,ix2,ix3),igv2(:,ig,ikt)
              !   write(6,"(a,10i8)")'ix1 ix2 ix3 ig igx=', ix1,ix2,ix3, ig, igx(ig,ix1,ix2,ix3)
              ! endif
            enddo
          enddo
        enddo
      enddo

!!... fill APW part
      print *,'goto sss section'
      do ix1=0,1
      do ix2=0,1
      do ix3=0,1
        print *, 'ix1,ix2,ix3 wgc =',ix1,ix2,ix3, wgc(ix1,ix2,ix3)
      enddo
      enddo
      enddo

      do j1=1,ndimh
        do j2=1,ndimh
          if(j1<=ldim.and.j2<=ldim) cycle
          sss=0d0
          www=0d0
          do ix1=0,1
          do ix2=0,1
          do ix3=0,1
            !if(wgc(ix1,ix2,ix3)==0d0) cycle
               if(j1<=ldim) then
                  j1x= j1
                else
                  j1x= igx(j1-ldim,ix1,ix2,ix3)+ldim  !G1
                endif
                if(j2<=ldim) then
                  j2x= j2
                else
                  j2x= igx(j2-ldim,ix1,ix2,ix3)+ldim  !G2
                endif
                if(j1x == 999999+ldim) goto 888
                if(j2x == 999999+ldim) goto 888
                ik1= ikc(1,ix1,ix2,ix3)
                ik2= ikc(2,ix1,ix2,ix3)
                ik3= ikc(3,ix1,ix2,ix3)
                sss = sss+ wgc(ix1,ix2,ix3)* sfz(ik1,ik2,ik3,j1x,j2x) !,isp)
                www = www+ wgc(ix1,ix2,ix3)
 888        continue
          enddo                                         ! sigma at eight corners
          enddo
          enddo
          sene(j1,j2)= sene(j1,j2) + sss/www
c         print *,'888 loop out', ix1,ix2,ix3,j1-ldim,j2-ldim
        enddo
      enddo
      deallocate(igx)
 1019 continue
      end subroutine seneinterp2



      subroutine seneinterp(nbas,nsp,ndimh,slat,sham,sbz,ifis,rsrnge,
     &nk1,nk2,nk3,sfz) !,qin,isp,sene)
      use m_struc_def  !Cgetarg
      use m_hft2rs
      use m_globalvariables
      use m_sigflg
C- takao. obtain sene at given q (qin).
C
Ci sbz, and others
Co sfz
C ------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,ifis,ndimh,lwsig
C     integer osig(nbas),otau(nbas),oppi(nbas)
      type(s_lat)::slat
      type(s_ham)::sham
      type(s_bz)::sbz

C     double precision ssite(1),sspec(1)
      double precision rsrnge
C ... Local parameters
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nsp,nglob,ival,ldham(16),hreal,lrsig,
     .moditp,nttabs,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .ndhrs,j1,k1,k2,k3,iq1,nl,nspc,mxorb,nqsig,modsgp,nmin,nmax,
     .fopna,kcplx,lrot,iprint,stdo,lssym,ledit,niax,nscnts,
     .ndims,ndimz,iq,n123(4),lcore,lhigh
      parameter (niax=10)
c      equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
      integer:: ohrss
      real(8),pointer :: rv_p_oqp(:) =>NULL()
ckino Jan.4.2012
      real(8),allocatable,target :: rv_a_oqp(:)
      integer ,allocatable :: gstar_iv(:)
      integer ,allocatable :: ipq_iv(:)
      real(8) ,allocatable :: wgt_rv(:)
      complex(8) ,allocatable :: wk_zv(:)
      complex(8) ,allocatable :: sigm_zv(:)
      complex(8) ,allocatable :: sigm2_zv(:)
      real(8) ,allocatable :: delt_rv(:)
      integer ,allocatable :: istb2_iv(:)

      double precision tmp(3)
      real(8),allocatable:: evls(:),evlz(:),sigii(:)
      complex(8),allocatable:: siglda(:,:),z(:,:),sigo(:,:)
C     For offset q mesh
      integer is(3),lshft(3),ifac(3),lqoffo
      double precision rb(3,3),qb(3,3),qoffi(3),qoffo(3)
C MPI
c      logical mlog
c      integer procid,mpipid,master
C ... for mixing sigma files: Use alf1*sigm + alf2*sigm1
      double precision alf(2)
C     integer os
      double precision plat(3,3),qp(3),tolq,rsstol,sigp(10)  ,qlat(3,3)
      double precision emin,emax,asig,bsig,efit,qoff(3),rotm(3,3),ddot
      equivalence (emin,sigp(3)),(emax,sigp(5))
      equivalence (asig,sigp(6)),(bsig,sigp(7)),(efit,sigp(8))
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      double precision qk
      integer jj1,jj2,jj3,k
      integer ::iwdummy
      real(8):: eseavr
      integer:: napw_in
      complex(8)::sfz(nk1,nk2,nk3,ndimh,ndimh,nsp)
      real(8):: qin(3)
ckino delete integer(4) def.      integer(4):: ndhamx,debugmode
      integer:: ndhamx,debugmode
      complex(8):: sene(ndimh,ndimh),zv_dummy(1)
      real(8):: rv_dummy(1)
      integer:: i_copy_size,nhrss,sigswitch
      real(8),allocatable::hrss(:)
      logical:: l_dummy_isanrg, isanrg, debug=.true.
!!
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .(jj2*ifac(2)-1)*qb(k,2) +
     .(jj3*ifac(3)-1)*qb(k,3)
c
      call tcn('seneinterp')
      print *
      print *,'seneinterp: startxxx'
      mode=12 !lrsig=12
      lwsig=0

C --- Setup and printout ---
c      procid = mpipid(1)
c      master = 0
c      mlog = cmdopt('--mlog',6,0,outs)

C     tolq is allowed fuzz in qp
      tolq = 1d-6
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      ledit = 0
      lrsig = mod(mode,10)
      lwsig = 0
      moditp = mod(mode/10,10)
      alf(1) = 1
      alf(2) = 0
      nscnts = 0
      lqoffo = 0
      lfbzout = .false.

C     Switch flagging whether input file has symops or not
      hreal = 0
      if (lrsig .ge. 4) hreal = 1
      if (lrsig .eq. 0) return

      ltrans = mod(mode/1000,10) .eq. 0
      kcplx = 1
      if (mod(mode/1000,10) .eq. 1) kcplx = 0
      lphase = mod(mode/1000,10) .eq. 0
      lssym = mod(mode/10000,10)
Ckino isanrg is logical function,       call isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)
      l_dummy_isanrg=isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)

C ... Command-line switches for input sigma and defaults
#ifdef USEOPTIONAL
      call sigflg (dc= ' ' , sio='rsig' , strn=' ' 
     . , lfbz=lfbzin  , lonesp=lonesp , qpoff=qoffi  )
#else
      call sigflg ( ' ' , 'rsig' , ' ' , iwdummy , iwdummy , iwdummy 
     . , lfbzin , iwdummy , lonesp , qoffi , iwdummy )
#endif

      if (cmdopt('--rsig',6,0,outs)) then
        out2 = outs(7:)
        dc = out2(1:1)
#ifdef USEOPTIONAL
        call sigflg ( dc=dc , sio='rsig' , strn=out2 
     .   , lfbz=lfbzin  ,lonesp= lonesp , qpoff=qoffi  )
#else
        call sigflg ( dc , 'rsig' , out2 , iwdummy , iwdummy , iwdummy 
     .   , lfbzin , iwdummy , lonesp , qoffi , iwdummy )
#endif

      endif
#if 1
      if ( lfbzin .or. mod(lssym,2).eq.1 ) then
        lfbzin = .true.
      else
        lfbzin = .false.
      endif
#else
      lfbzin = lfbzin .or. mod(lssym,2)
#endif

      if (cmdopt('--mixsig=',9,0,outs)) then
        j = 0
        i = parg('--mixsig=',4,outs,j,len(outs),
     .  ', ',2,2,ifac,alf)
        if (i .ne. 1 .and. i .ne. 2) call rx
     .  ('rdsigm: failed to parse arguments to --mixsig=')
      endif

      call info(10,1,0,' RDSIGM: read file sigm and create '//
     .'%?#n#REAL#COMPLEX# sigma(R) by FT ...',hreal,0)

c      i_copy_size=size(sham%ldham)
c      call i8icopy(i_copy_size,sham%ldham,1,ldham,1)
ckino Dec.22.2011:         iv_p_ooffh => sham%iv_p_ooffh
ckino Dec.22.2011:         iv_p_oiprmb => sham%iv_p_oindxo
c      i_copy_size=size(sham%sigp)
c      call dcopy(i_copy_size,sham%sigp,1,sigp,1)
      ldham =sham%ldham
      sigp  =sham%sigp
      rsstol=sham%rsstol

      modsgp = nint(sigp(1))
      nmin   = nint(sigp(2))
      nmax   = nint(sigp(4))

!! info blocks
      call info2(30,0,0,
     .'%9fSigm will be approximated by:  '//
     .'%?#(n==0)#Simple bloch sum##%-1j'//
     .'%?#(n==1)#diagonal Sigma for high and low states##%-1j'//
     .'%?#(n==2)#Perturbation for higher kappa blocks##%-1j'//
     .'%?#(n==3)#Interpolation from known points##%-1j'//
     .' ',moditp,0)
      if (moditp .eq. 1 .or. moditp .eq. 2) then
        if( sigswitch()==0 ) then
        call info5(30,0,0,
     .  '%9fApproximate sigma '//
     .  '%?#(n<0)#for energies E(lda)<%d; and %-2j#%-1j#'//
     .  '%?#(n>0)#for states up to %-1jn=%i; and ##%j'//
     .  '%?#(n<=0)#for energies E(lda)>%d%-2j#%-1j#'//
     .  '%?#(n>0)#for states above %-1jn=%i##%j',
     .  nmin,emin,nmax,emax,0)
        call info5(30,0,0,'%9fFor high states '//
     .  '%?#(n==0)#Sigii > %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
     .  '%?#(n==1)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
     .  '%?#(n==3)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
     .  '%?#(n==2)#%,3;4d < Sigii < %,3;4d%-2j##%-1j'//
     .  ' ',modsgp,asig,bsig,0,0)
        if (nmax .gt. 0) call info5(30,0,0,
     .  '%9fLinear fit to diagonal sigma for E(lda)>%d'//
     .  '%?#(n==3)# and n<%i'//
     .  ' ',efit,modsgp,nmax,0,0)
        if (nmax .eq. 0) call info5(30,0,0,
     .  '%9fLinear fit to diagonal sigma for E(lda)>%d'//
     .  '%?#(n==3)# and E(lda)<%d'//
     .  ' ',efit,modsgp,emax,0,0)
        endif
Ckino isanrg is logical function,         call isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
        l_dummy_isanrg=isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
      endif
      if (lfbzin .or. ddot(3,qoffi,1,qoffi,1).ne.0) then
        call info5(30,0,0,
     .  '%?#(n==1)#%9fFile sigm in saved in FBZ: '//
     .  'symmetrization suppressed##'//
     .  '%?#(n==1)#%N%9fFile k-mesh is offset by:%3:2g#%j#'//
     .  ' ',isw(lfbzin),isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi,0,0)
      endif
      if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
        call info5(30,0,0,
     .  '%9fUse for sigma: %d*(file sigm)%?#n# + '//
     .  '%d*(file sigm1)',alf(1),isw(alf(2).ne.0),alf(2),0,0)
      endif

!! --- Read header of sigma(orbital basis) from file ---
C ... Read k-mesh parameters
c      if (procid .eq. master) then
C       Require file contains sigm(orbital basis)
        call iosigh(0,nscnts,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
     .  lshft(3),ifis)
        print *,'nscnts111=',nscnts
        if (nscnts .ne. 0 .and. nscnts .ne. 5) then
          call iosigh(2,0,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
     .    lshft(3),ifis)
        endif
c      endif
!!

c      call mpibc1(i,1,2,.false.,'rdsigm','nsp')
c      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
c      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
c      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')

      if (i .eq. nsp) then
        lsplts = .false.
      elseif (i .gt. nsp) then
        call rx('rdsigm: sigm file spin polarized but calculation is not')
      else
        lsplts = .true.
      endif

C     Check for consistency in remaining ndimh, nqp and set file sigma file pointer past header
ctakao think this is not necessary.

c comment out this block, since ndimh check can stop here when ndimh=ldim in PMT case. aug2012takao
c      if (procid .eq. master) then 
c        rewind ifis
c        if (lfbzin) nqp = nk1*nk2*nk3
c        call iosigh(3,nscnts,i,ndimh,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
c      endif

c      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
c      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
c      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
      call info5(30,0,0,'%9fsigm file has %i irreducible QP: nk = ( %i %i %i )  shift=%3:1l',nqp,nk1,nk2,nk3,lshft)
      if (lsplts) call info0(30,0,0,'%9f(warning) sigm file not spin pol .. splitting spins') !sanity check

ctakao feb2012
c      i_copy_size=size(slat%plat)
c      call dcopy(i_copy_size,slat%plat,1,plat,1)
ckino Dec.22.2011:         rv_p_opos => slat%rv_p_opos
ckino Dec.22.2011:         iv_p_oistab => slat%iv_p_oistab
ckino Dec.22.2011:         rv_p_og => slat%rv_p_osymgr
ckino Dec.22.2011:         rv_p_oag => slat%rv_p_oag
C     Get leading dimension of hrs
ckino Dec.22.2011:         call offshp ( iv_p_oiprmb , 1 , nbas , 2 , 0 , ndimh , 0 , ndhrs
ckino Dec.27.2011:         call offshp ( sham%iv_p_oindxo , 1 , nbas , 2 , 0 , ndimh , 0 , ndhrs

      if(debug) print *,'seneinterp: debug EEEEEEEEEE'
      plat = slat%plat
      nsgrp= slat%nsgrp
      call offshp ( sham%iv_a_oindxo , 1 , nbas , 2 , 0 , ndimh , 0 , ndhrs )
C     Setup for FFT
c      call fftz30(nk1,nk2,nk3,k1,k2,k3) !takao comment out because k1 k2 k3 are supplied at the begining of routine.
c      if (nk1.ne.k1 .or. nk2.ne.k2 .or. nk3.ne.k3) call rx('rdsigm: not ready for FFT w/ dimensions ne no div.')
      k1=nk1
      k2=nk2
      k3=nk3

C ... Number of group operations for input file sigma
      if (mod(lssym,2) .eq. 0) nsgrps = nsgrp
      if (mod(lssym,2) .eq. 1) nsgrps = 1
      if (lfbzin) nsgrps = 1
C ... Make is,ifac,qb,qlat,qoff
      do i = 1, 3
        llshft(i) = lshft(i) .ne. 0
      enddo
      call pshpr(0)
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
      call poppr
      do i = 1, 3
        qoff(i) = qk(i,1,1,1) + qoffi(i)
      enddo
      if(debug) print *,'seneinterp: debug FFFFFFF 2222'
C ... Setup for r.s. hamiltonian: allocate and create iax,ntab; allocate hrs
      i = 1000 + 10*hreal + 0
C      if(debug) print *, '!!'; i = 1000*0 + 10*hreal + 0
      if (mod(lssym,4) .ge. 2) i = 10*hreal + 0
      if (ltrans) i = i + 10000
      print *,' rrrrrrrrrr go to hft2rs i=',i
      call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , 1 ,
     .    nsp , nbas , slat%rv_a_osymgr , slat%rv_a_oag , nsgrp , rsrnge , sham%iv_a_ooffH
     .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , slat%rv_a_opos , sham%iv_a_ontabs ,
     .    sham%iv_a_oiaxs , ndhrs , rv_dummy )
      if(debug) print *,'seneinterp: debug FFFFFFF size iv_p_ontabls ndhrs=',size(sham%iv_a_ontabs),ndhrs

ckino Dec.14.2011:         nttabs = ival ( iv_p_ontabs , nbas + 1 )
ckino Dec.27.2011:         nttabs = ival ( sham%iv_a_ontabs , nbas + 1 )
C     Allocate memory for hrs
ckino Dec.22.2011:           allocate(rv_p_ohrs(ndhrs**2*nttabs*nsp))
ckino Dec.28.2011:           allocate(sham%rv_p_ohrs(ndhrs**2*nttabs*nsp))
ckino Dec.22.2011:           rv_p_ohrs(:)=0d0
ckino Dec.28.2011:           sham%rv_p_ohrs(:)=0d0
ckino Dec.22.2011:           allocate(rv_p_ohrs(2*ndhrs**2*nttabs*nsp))
ckino Dec.28.2011:           allocate(sham%rv_p_ohrs(2*ndhrs**2*nttabs*nsp))
ckino Dec.22.2011:           rv_p_ohrs(:)=0d0
ckino Dec.28.2011:           sham%rv_p_ohrs(:)=0d0
ckino Dec.22.2011:         allocate(rv_p_oqsig(abs(3*mxkp)))
ckino Jan.4.2012
ckino Dec.22.2011:         if (3*mxkp<0) rv_p_oqsig(:)=0.0d0
ckino Jan.4.2012       allocate(rv_p_oqp(abs(3*mxkp)))

      nttabs = ival ( sham%iv_a_ontabs , nbas + 1 )
      if(debug) print *,'nttabs=',nttabs
      if(hreal .eq. 1) then
        if (allocated(sham%rv_a_ohrs)) deallocate(sham%rv_a_ohrs)
        allocate(sham%rv_a_ohrs(ndhrs**2*nttabs*nsp))
        sham%rv_a_ohrs(:)=0d0
      else
        if (allocated(sham%rv_a_ohrs)) deallocate(sham%rv_a_ohrs)
        allocate(sham%rv_a_ohrs(2*ndhrs**2*nttabs*nsp))
        sham%rv_a_ohrs(:)=0d0
      endif

C --- Get list of irreducible k-points, and ipq and gstar arrays ---
      mxkp = nk1*nk2*nk3
      if (associated(sham%rv_p_oqsig)) deallocate(sham%rv_p_oqsig)
      allocate(sham%rv_p_oqsig(abs(3*mxkp)))
      if (3*mxkp<0) sham%rv_p_oqsig(:)=0.0d0
      allocate(rv_a_oqp(abs(3*mxkp)))
      rv_p_oqp=>rv_a_oqp
ckino end
      if (3*mxkp<0) rv_p_oqp(:)=0.0d0
      allocate(gstar_iv(abs(-mxkp-1)))
      if (-mxkp-1<0) gstar_iv(:)=0
      gstar_iv = - 2
      allocate(ipq_iv(mxkp))
      allocate(wgt_rv(abs(-mxkp)))
      if (-mxkp<0) wgt_rv(:)=0.0d0
      print *,'seneinterp: debug HHHHH'
      call info(20,1,0,' q-points in full BZ where sigma calculable ...',0,0)
      call bzmesh(plat, qb, nk1, nk2, nk3, llshft, iwdummy, 0 , ipq_iv, sham%rv_p_oqsig, wgt_rv, nqsig, mxkp, 0, 0 )
ckino Dec.22.2011:       .  , 0 , ipq_iv , rv_p_oqsig , wgt_rv , nqsig , mxkp , 0 , 0 ) 
ckino Dec.22.2011:         sham%rv_p_oqsig => rv_p_oqsig
ckino Dec.22.2011:          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og
ckino Dec.28.2011:          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_p_osymgr
      print *,'seneinterp: debug IIIIIIII'
      sham%nqsig=nqsig
      call dpzero ( wgt_rv , mxkp )
      call info2(20,0,0,' Irr. qp for which sigma is calculated ...',0,0)
      call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_a_osymgr 
     .  , nsgrps , ipq_iv , rv_p_oqp , wgt_rv , nqps , mxkp , gstar_iv, 0 )

Ckino isanrg is logical function,       call isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
      l_dummy_isanrg=isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
C     call prmx('q for which sigma is read',w(oqp),3,3,nqp)
      call info2(30,0,0,'%?#(n==1)#%1fFile sigm has k-offset=%3:2g '//
     .'relative to the above#%j# ',isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi)

C ... titus
CC     List of sites
C      call getoffH(w(ooffH),offH,n0H,nkap0,nbas)
C      ib1 = 1
C      ib2 = nbas
C      lidim = offH(4,1,nbas+1)
C
C      print*,'lidim',lidim
C
C      do 100 ib = ib1, ib2
C      rewind ifis
C      call iosigh(0,nscnts,nsp1,ndimh,nk1,nk2,nk3,nqp,
C     .  lshft(1),lshft(2),lshft(3),ifis)
C      nlmaa = offH(4,1,ib+1) - offH(4,1,ib)
C      print*,'nlmaa',nlmaa
C      print*,'offH(4,1,ib)',offH(4,1,ib)
C      print*,'offH(4,1,ib+1)',offH(4,1,ib+1)
C ... titus



C --- Generate hrs = sigma(T) from file sigma(k) ---
C ... Loop over spins, keeping sig in full BZ for only one spin
c      if (procid .eq. master) then
        do  isp = 1, nsp
C     If sigma file not spin polarized, use sigma from spin 1
          if (isp .eq. 2 .and. lsplts) then
            call iosigh(3,0,1,ndimh,nk1,nk2,nk3,nqp,
     .      lshft(1),lshft(2),lshft(3),ifis)
          endif

C     If to mix in some other sigma file, set up ifis2
          if (alf(2) .ne. 0) then
            if (isp .eq. 1 .or. isp .eq. 2 .and. lsplts) then
              ifis2 = fopna('sigm1',-1,4)
              call iosigh(3,0,nsp,ndimh,nk1,nk2,nk3,nqp,
     .        lshft(1),lshft(2),lshft(3),ifis2)
            endif
          endif

cccccccccccccccccccccccccccccccccccccccccccc
          print *,' goto 8888 isp nsp=',isp,nsp
          goto 8888
cccccccccccccccccccccccccccccccccccccccccccc







C ... File sigma in irr BZ -> sigma in full BZ by rotations
c      call defcc(osfz,k1*k2*k3*ndimh**2)

C ... titus
C      call defcc(osfz,k1*k2*k3*ndimh*nlmaa)
C      call defcc(osfz1,k1*k2*k3*ndimh*nlmaa)
C ... titus
          allocate(wk_zv(ndimh**2))

Changenglob          nl = nglob('nl')
          nl = globalvariables%nl
          nspc = 1
          allocate(sigm_zv(ndimh**2))

ccccccccccccccccccccccccc
ctakaox
c      call defcc(oss,ndimh**2)
ccccccccccccccccccccccccc
cccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do  iq1 = 1, nqps
            read(ifis) qp,  eseavr !june2009  add eseavr
            if(debugmode()>0) print *,' uuu2 iq1 eseavr=',iq1,eseavr
            call dpscop(qp,tmp,3,1,1,1d0)
            call dpsadd ( tmp , rv_p_oqp , 3 , 1 , 3 * iq1 - 2 , - 1d0 )

            call dpsadd(tmp,qoffi,3,1,1,-1d0)

            !tmp= qp - w((oqp(3*(iq1-1)+1:3*(iq1-1)+3))+qoffi)
            if (.not. latvec(1,tolq,plat,tmp)) then
              call dpscop ( rv_p_oqp , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 )

              if (lssym .ge. 4) then
                call info5(20,0,0,' rdsigm (warning) qp %i : '//
     .          'expected%3;8,4D read%3;8,4D',iq1,tmp,qp,0,0)
                call dpscop(tmp,qp,3,1,1,1d0)
              else
                print 456, iq1,tmp,qp
  456           format(' rdsigm: error on reading qp no',i4/
     .          ' Expected qp=',3f12.6/
     .          '     File qp=',3f12.6)
                call rx(' incompatible q-mesh')
              endif
            endif
C       File read sigma
            print *,' nscnts =',nscnts
            if (nscnts == 5) call dpdump(tmp,1,ifis)
            call dpdump ( sigm_zv , ndimh * * 2 * 2 , ifis )


C        if (isp .eq. 2) then
C          print *, iq1,isp
C         call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
C        endif

cccccccccccccccccccccccccccccccccccccc
ctakaox !see also bndfp.F and hambls.F
c        if(iq1==1) then
c          print *,'reading ovlmat'
c          ifio = fopna('ovlmat',-1,4)
c        endif
c        call dpdump(w(oss),ndimh**2*2,ifio)
c        print *,'read ovlmat for iq1=',iq1
c        call divctakao(w(osigm),w(oss),ndimh)
ccccccccccccccccccccccccccccccccccccc
c

c ... scalced sigm. alf(1:2)
            if (alf(1) .ne. 1d0) then
              call dscal ( ndimh * * 2 * 2 , alf ( 1 ) , sigm_zv , 1 )

            endif
            if (alf(2) .ne. 0d0) then
              read(ifis2) qp
              call dpscop(qp,tmp,3,1,1,1d0)
              call dpsadd ( tmp , rv_p_oqp , 3 , 1 , 3 * iq1 - 2 , - 1d0 )

              call dpsadd(tmp,qoffi,3,1,1,-1d0)
              if (abs(tmp(1))+abs(tmp(2))+abs(tmp(3)) .gt. tolq) then
                call dpscop ( rv_p_oqp , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 )

                if (lssym .ge. 4) then
                  call info5(20,0,0,' rdsigm (warning) qp %i : '//
     .            'expected%3;8,4D read%3;8,4D (file sigm1)',iq1,tmp,qp,
     .            0,0)
                else
                  call rx(' file sigm1 has incompatible irr mesh')
                endif
              endif
C         call dpsadd(qp,w(oqp),3,1,3*iq1-2,1d0)
              allocate(sigm2_zv(ndimh**2))

              call dpdump ( sigm2_zv , ndimh * * 2 * 2 , ifis2 )

              call daxpy ( ndimh * * 2 * 2 , alf ( 2 ) , sigm2_zv , 1 , 
     .        sigm_zv , 1 )

              if (allocated(sigm2_zv)) deallocate(sigm2_zv)

            endif

c ... rotate sigm. w(osigm) at qp ---> w(ofbz) at stars of qp.
C       call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
            write(6,"(a,13f13.5)")'goto hamfbk3 xx input qp=', qp
c        napw_in=0 ! for test

ctakao feb2012
c            i_copy_size=size(slat%qlat)
c            call dcopy(i_copy_size,slat%qlat,1,qlat,1)


            call hamfb3k ( qp , iq1 , nk1 , nk2 , nk3 , k1 , k2 , k3 , ipq_iv
c     .      , napw_in , ndimh , ndimh , ndimh , qb , plat , qlat , ifac , 
     .      , napw_in , ndimh , ndimh , ndimh , qb , ifac , 
     .      gstar_iv , sigm_zv , sfz ( 1 , 1 , 1 , 1 , 1 , isp ) )

            if(debugmode()>0) write(6,"(a,3f13.5)")'end of hamfbk3'
c
c        i = 100*kcplx + 00
c        if (lphase) i = i+40
c        print *,' hamfb3: i =',i
c        call hamfb3(nbas,nl,w(ooffH),w(oiprmb),i,w(opos),iq1,nk1,
c     .    nk2,nk3,k1,k2,k3,w(oipq),w(oistab),w(og),w(oag),
c     .    w(ogstar),ifac,ndimh,ndimh,ndimh,nspc,qb,w(osigm),w(owk),
c     .    w(owk),w(osfz))

          enddo








ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 8888     continue
          print *,' 8888 continue'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccc
ctakaox !see also bndfp.F and hambls.F
c      print *,' end of iq1 loop'
c      call fclose(ifio)
c      call rlse(oss)
c      print *,' end of iq1 loop xxxx'
cccccccccccccccccccccccccccccccccccccccccccc


C      print *, 'one element of sigma(k)'
C      call pvtrod(w(osfz),k1,k2,k3,ndimh)

C ... FT sfz and copy sfz(T) to hrs
          i = 10*hreal + 1
          if (ltrans) i = i + 10000
C ... titus
C      call hft2rs(i,nk1,nk2,nk3,k1,k2,k3,qoff,isp,nsp,nbas,w(og),w(oag)
C     .  ,nsgrps,rsrnge,w(ooffH),ib,ib,nlmaa,ndimh,w(osfz),plat,
C     .  w(opos),ontabs,
C     .  oiaxs,ndhrs,w(ohrs))
C ... titus
          call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , isp
     .         , nsp , nbas , slat%rv_a_osymgr , slat%rv_a_oag , nsgrps , rsrnge , sham%iv_a_ooffH 
     .         , 1 , nbas , ndimh , ndimh , sfz ( 1 , 1 , 1 , 1 , 1 , isp )
     .         , plat , slat%rv_a_opos ,sham%iv_a_ontabs , sham%iv_a_oiaxs , ndhrs , sham%rv_a_ohrs )
          print *,' 8888 out hft2rs i ndimh ndhrs=',i,ndimh,ndhrs
C      print *, 'print out sigma(T)'
C      call pvtrof(w(ohrs),ndhrs,nsp,nk1,nk2,nk3,w(oiaxs),1,
C     .  nttabs,w(osfz),ndimh)
C      stop

C ... End loop over spins
c      call rlse(osfz)
        enddo
C ... titus
C 100  continue
C     exit of MPI master loop
c      endif

      if(debug) print *,' a8 111111111111111'
      if (allocated(sigm_zv)) deallocate(sigm_zv)
      if (allocated(wk_zv)) deallocate(wk_zv)
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq_iv)) deallocate(ipq_iv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)
      sham%ndhrs=ndhrs

      if(debug) print *,' a8 22222222222222222'

ckino Dec.14.2011:         sham%iv_p_ontabs => iv_p_ontabs
ckino Dec.14.2011:         sham%iv_a_ontabs => iv_a_ontabs
ckino Dec.27.2011:        call move_alloc(from= iv_a_ontabs ,to= sham%iv_a_ontabs )
ckino Dec.14.2011:         sham%iv_p_oiaxs => iv_p_oiaxs
ckino Dec.14.2011:         sham%iv_a_oiaxs => iv_a_oiaxs
ckino Dec.27.2011:        call move_alloc(from= iv_a_oiaxs ,to= sham%iv_a_oiaxs )
ckino Dec.22.2011:         sham%rv_p_ohrs => rv_p_ohrs

C ... Check that FT replicates sigma(k) at file qp to tolerance rsstol
C     rsstol = 1d-6
C     print *, '!!'
C     if (alf(1) .ne. 1 .or. alf(2) .ne. 0 .or. .true.) then
      if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
        call info0(20,0,0,' skipping check of Bloch summed sigma'//
     .  ' (sigma scaled) ')
      else
        j1 = 0
        if (lssym .ge. 4) j1 = 2
c        if (procid .eq. master) then
          call dpzero(tmp,3)
          print *,'a8 22222222222 goto chksgr'
          call chksgr ( j1 , ltrans , kcplx , plat , nsp , ndimh , ifis
     .      , hreal , sham%iv_a_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_a_ohrs , ndhrs 
     .      , rsstol , i , rv_p_oqp , nbas , 0 , rotm , iwdummy )
c        endif
      endif
      if(debug) print *,' a8 end of chksgr 3333333333'
      if ( .not. ( cmdopt ( '--wsig' , 6 , 0 , outs ) .or.cmdopt ( '-wsig' , 5 , 0 , outs ) ) ) then
        if (associated(rv_p_oqp)) nullify(rv_p_oqp)
        if (allocated(rv_a_oqp)) deallocate(rv_a_oqp)
      endif

      if(debug) print *,' a8 3333333333'

C ... Broadcast R.S. sigma before symmetrization
c      if (hreal .eq. 1) then
c        call mpibc1( sham%rv_a_ohrs , ndhrs * * 2 * nttabs * nsp, 4 , mlog, 'rdsigm' , 'hrs' )
c      else
c        call mpibc1( sham%rv_a_ohrs , 2 * ndhrs * * 2 * nttabs * nsp ,4, mlog , 'rdsigm' , 'hrs' )
c      endif

C --- Symmetrize hrs ---
C     Best to use nsgrp here instead of nsgrps.
C     However, neighbor table was generated for nsgrps, and may be
C     increase if nsgrp>nsgrps.  This branch will fail in that case.
      if (mod(lrsig,4) .ge. 2 .and. nsgrp .gt. 1) then
cC       Allocate memory for hrs
c        if (hreal .eq. 1) then
c          call defdr(ohrss,ndhrs**2*nttabs*nsp)
cCchp1           call dcopy(ndhrs**2*nttabs*nsp*1,w(ohrs),1,w(ohrss),1)
c           call dcopy ( ndhrs * * 2 * nttabs * nsp * 1 , rv_p_ohrs , 1 ,
c     .     w ( ohrss ) , 1 )
c        else
c          call defdc(ohrss,ndhrs**2*nttabs*nsp)
cCchp1           call dcopy(ndhrs**2*nttabs*nsp*2,w(ohrs),1,w(ohrss),1)
c           call dcopy ( ndhrs * * 2 * nttabs * nsp * 2 , rv_p_ohrs , 1 ,
c     .     w ( ohrss ) , 1 )
c        endif
ckino Dec.22.2011:          nhrss= size(rv_p_ohrs)
ckino Dec.28.2011:          nhrss= size(sham%rv_p_ohrs)
        nhrss= size(sham%rv_a_ohrs)
        allocate(hrss(nhrss))
ckino Dec.22.2011:          hrss = rv_p_ohrs
ckino Dec.28.2011:          hrss = sham%rv_p_ohrs
        hrss = sham%rv_a_ohrs
Changenglob        mxorb = nglob('mxorb')
        mxorb = globalvariables%mxorb
C       10s digit distributes ri-rj pairs to avg (ri-rj),(rj-ri) pairs
C       Use 10*2 for hermitian matrices, 10*1 for symmetric ones
        i = 1-hreal + 10*2 + 100*(1-hreal) + 100000
        allocate(istb2_iv(nsgrp*nbas))

ckino Dec.22.2011:           call istbpm ( iv_p_oistab , nbas , nsgrp , istb2_iv )
ckino Dec.28.2011:           call istbpm ( slat%iv_p_oistab , nbas , nsgrp , istb2_iv )
Changenglob        nl = nglob('nl')
ckino Dec.22.2011:                call rsmsym ( i , plat , mxorb , iv_p_oiprmb , ndimh , nbas ,
ckino Dec.27.2011:                call rsmsym ( i , plat , mxorb , sham%iv_p_oindxo , ndimh , nbas ,
ckino Dec.14.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_p_ontabs , iv_p_oiaxs
ckino Dec.14.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_p_ontabs , iv_a_oiaxs
ckino Dec.22.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_a_ontabs , iv_a_oiaxs
ckino Dec.27.2011:       .        slat%rv_p_opos , nl , nsp , 1 , nttabs , iv_a_ontabs , sham%iv_a_oiaxs
ckino Dec.27.2011:       .        slat%rv_p_opos , nl , nsp , 1 , nttabs , sham%iv_a_ontabs , sham%iv_a_oiaxs
ckino Jan.04.2012:       .        slat%rv_p_opos , nl , nsp , 1 , nttabs , sham%iv_a_ontabs , sham%iv_a_oiaxs 
ckino Dec.22.2011:       .        , rv_p_og , istb2_iv , nsgrp , ndhrs , hrss , rv_p_ohrs )
ckino Dec.22.2011:       .        , slat%rv_p_osymgr , istb2_iv , nsgrp , ndhrs , hrss , rv_p_ohrs )
ckino Dec.28.2011:       .        , slat%rv_p_osymgr , istb2_iv , nsgrp , ndhrs , hrss , sham%rv_p_ohrs )
ckino Dec.28.2011:       .        , slat%rv_a_osymgr , istb2_iv , nsgrp , ndhrs , hrss , sham%rv_p_ohrs ) 
        call istbpm ( slat%iv_a_oistab , nbas , nsgrp , istb2_iv )
        nl = globalvariables%nl

       if(debug) print *,' a8 44444444444444'
        call rsmsym ( i , plat , mxorb , sham%iv_a_oindxo , ndimh , nbas ,
     .        slat%rv_a_opos , nl , nsp , 1 , nttabs , sham%iv_a_ontabs , sham%iv_a_oiaxs
     .        , slat%rv_a_osymgr , istb2_iv , nsgrp , ndhrs , hrss , sham%rv_a_ohrs )
        if (allocated(istb2_iv)) deallocate(istb2_iv)
        deallocate(hrss)
C        Debugging
C        ifis2 = fopna('out',-1,0)
C        print *, 'hello',ndhrs,nttabs
C        call ywrm(0,'sig',2,ifis2,'(9f20.10)',w(ohrs),1,ndhrs**2,
C     .    ndhrs**2,nttabs*nsp)
C        call rx0('done')

      if(debug) print *,' a8 55555555555555'

C   ... Check how well symmetrized FT replicates sigma(k) at file qp
c        if (procid .eq. master) then
          rsstol = 0
          if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
            call info0(20,0,0,' skipping check of symmetrized sigma'//
     .      ' (sigma scaled) ')
          else
            call chksgr ( 0 , ltrans , kcplx , plat , nsp , ndimh , ifis
ckino Dec.14.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs
ckino Dec.22.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs
ckino Dec.22.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs
ckino Dec.27.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs
ckino Dec.27.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs
ckino Dec.28.2011:       .          , hreal , sham%iv_a_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .          , hreal , sham%iv_a_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_a_ohrs , ndhrs 
     .          , rsstol , i , 0 , nbas , 0 , rotm , iwdummy )
c          endif
        endif
C   ... Broadcast R.S. sigma after symmetrization
c        if (hreal .eq. 1) then
c          call mpibc1 ( sham%rv_a_ohrs , ndhrs * * 2 * nttabs * nsp , 4 , mlog, 'rdsigm' , 'hrs' )
c        else
c          call mpibc1 ( sham%rv_a_ohrs , 2 * ndhrs * * 2 * nttabs * nsp , 4, mlog , 'rdsigm' , 'hrs' )
c        endif
      endif
      if(debug) print *,' a8 6666666666666'

C ... lwsig=0 modes that internally transform sigma, save and exit
      if (lwsig .eq. 0) then
        if (cmdopt('--wsig',6,0,outs) .or. cmdopt('-wsig',5,0,outs)) then
c          if (procid .eq. master) then
            ifis2 = fopna('sigm2',-1,4)
            rsstol = 0
            out2 = outs(7:)
            if (outs(1:5) .eq. '-wsig') out2 = outs(6:)
            dc = out2(1:1)
            call sigflg(dc,'wsig',out2,ledit,lwsig,lrot,lfbzout,lnwmsh,
     .      lonesp,qoffo,rots)

            if (lonesp .eq. 1) then
              if (nsp .eq. 1) call rx(
     .        'onesp option nonsensical unless sigma file is spin pol')
ckino Dec.22.2011:                     call siged ( 1 , nbas , nsp , ndhrs , plat , rv_p_opos , ndimh
ckino Jan.04.2012:                     call siged ( 1 , nbas , nsp , ndhrs , plat , slat%rv_p_opos , ndimh
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_p_oiaxs , rv_p_ohrs
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino Dec.22.2011:       .             , iv_p_oiprmb , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino Dec.22.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs
ckino Dec.28.2011:       .             , sham%iv_a_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
              call siged ( 1 , nbas , nsp , ndhrs , plat , slat%rv_a_opos , ndimh
     .           , sham%iv_a_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_a_ohrs )
            endif
            if (ledit .eq. 1) then
ckino Dec.22.2011:                     call siged ( 0 , nbas , nsp , ndhrs , plat , rv_p_opos , ndimh
ckino Jan.04.2012:                     call siged ( 0 , nbas , nsp , ndhrs , plat , slat%rv_p_opos , ndimh
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_p_oiaxs , rv_p_ohrs
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino Dec.22.2011:       .             , iv_p_oiprmb , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino Dec.22.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs
ckino Dec.28.2011:       .             , sham%iv_a_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
              call siged ( 0 , nbas , nsp , ndhrs , plat , slat%rv_a_opos , ndimh
     .             , sham%iv_a_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_a_ohrs )
            endif

C       Handle case q-mesh for file output differs from file input
C       ltmp = true in cases where file symops artificially reduced
C        ltmp = mod(lssym,2) .eq. 1 .and. lwsig .ne. 3
C        if (lnwmsh .or. ltmp) then
c              i_copy_size=size(sbz%nkabc)
c              call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
c              i_copy_size=size(sbz%lshft)
c              call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)
ckino Jan.4.2012 needless operation               rv_p_oqp => sbz%rv_p_oqp
ckino Jan.4.2012               allocate(rv_p_oqp(abs(3*mxkp)))
            if (lnwmsh) then
              lshft=sbz%lshft
              nkxyz=sbz%nkabc
              nqp=sbz%nkp
              mxkp = nk1*nk2*nk3
              allocate(rv_a_oqp(abs(3*mxkp)))
              rv_p_oqp=>rv_a_oqp
              if (3*mxkp<0) rv_p_oqp(:)=0.0d0
              allocate(gstar_iv(abs(-mxkp-1)))
              if (-mxkp-1<0) gstar_iv(:)=0
              gstar_iv = - 2
              allocate(ipq_iv(mxkp))
              allocate(wgt_rv(abs(-mxkp)))
              if (-mxkp<0) wgt_rv(:)=0.0d0
              do  i = 1, 3
                llshft(i) = lshft(i) .ne. 0
              enddo
              call dpzero ( wgt_rv , mxkp )
              if (lfbzout) then
                call info0(20,1,0,
     .          ' Creating sigma for new mesh, full BZ (file sigm2) ...')
                call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
     .            , 0 , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , 0 , 0 )
              else
                call info0(20,1,0,
     .          ' Creating sigma for new mesh (file sigm2) ...')
ckino Dec.22.2011:                    call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og
ckino Dec.28.2011:                    call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_p_osymgr
                call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_a_osymgr 
     .            , nsgrp , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , gstar_iv 
     .            , 0 )
              endif
            elseif (lfbzout) then
cki              oqp = oqsig
ckino Dec.22.2011:                rv_p_oqp => rv_p_oqsig
              rv_p_oqp => sham%rv_p_oqsig
              nqp = nqsig
              call info2(20,1,0,
     .        ' Creating sigma for full BZ, %i qp (file sigm2) ...',nqp,0)
            else
              call info2(20,0,0,
     .        ' Writing sigma for irr BZ, %i qp (file sigm2) ...',nqp,0)
            endif
C   ... Setup rotation matrix
            if (rots .ne. ' ') then
              call a2rotm(rots,.false.,0,rotm)
              if (iprint() .ge. 20) then
                call info0(20,0,0,' Rotate sigma by rotation matrix:')
                write (stdo,350) ((rotm(i,j),j=1,3),i=1,3)
  350           format(3f11.6)
              endif
              lrot = lrot + 1
            endif

C   ... Setup phase matrix
            if (lrot .ge. 2) then
              allocate(delt_rv(abs(-3*nbas)))
              if (-3*nbas<0) delt_rv(:)=0.0d0
              call info0(0,0,0,' Phase shift of sigma.  '//
     .        'Translation vectors from shorps file:')
              call iopos ( .false. , - 1 , 'shorps' , nbas , delt_rv )
            endif
C       call prmx('q for which sigm is written',w(oqp),3,3,nqp)
            call iosigh(0,nscnts,nsp,ndimh,nk1,nk2,nk3,nqp,
     .      lshft(1),lshft(2),lshft(3),-ifis2)
            if (lqoffo .ne. 0) then
C         call prmx('qp',w(oqp),3,3,nqp)
              do  i = 1, nqp
                call dmsadd ( rv_p_oqp , 1 , qoffo , 1 , 1 , 3 , 1 , 1 , 3 * 
     .           i - 2 , 1 , 1d0 )

              enddo
C         call prmx('qp',w(oqp),3,3,nqp)
            endif
            call chksgr ( 1 , ltrans , kcplx , plat , nsp , ndimh , ifis2
ckino Dec.14.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs
ckino Dec.22.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs
ckino Dec.22.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs
ckino Dec.27.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs
ckino Dec.27.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs
ckino Dec.28.2011:       .          , hreal , sham%iv_a_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .          , hreal , sham%iv_a_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_a_ohrs , ndhrs 
     .          , rsstol , nqp , rv_p_oqp , nbas , lrot , rotm , delt_rv )
            call fclose(ifis2)
            call rx0('done writing sigma, file sigm2')
C      elseif (lwsig .eq. 1 .or. lwsig .eq. 2) then
C        if (lfbzout) call info0(20,0,0,' ')
C        if (lfbzout) goto 99
c          endif
        endif
      endif

C     Broadcast entire sham structure (shouldn't be needed)
C     call mpibc1(sham,i,4,mlog,'rdsigm','sham')
C     Broadcast neighbor table (shouldn't be needed)
C     call mpibc1(w(ontabs),nbas+1,2,mlog,'rdsigm','ntabs')
C     call mpibc1(w(oiaxs),niax*nttabs,2,mlog,'rdsigm','iaxs')


      if(debug) print *,' a8 999999999111111111111111'

C ... Transformations done by caller: write qp list to disk
   99 continue
      if (lwsig .ne. 0) then

        call info0(30,0,0,' ')

        if (lnwmsh) then

          i_copy_size=size(sbz%nkabc)
          call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
          i_copy_size=size(sbz%lshft)
          call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)
          nqp=sbz%nkp
          rv_p_oqp => sbz%rv_p_oqp


          call info2(20,0,0,
     .    ' Transform sigma on new k-mesh:  '//
     .    'nk=%3:1i  shft=%3:1l',nkxyz,lshft)
        endif

        if (lfbzout) then
          mxkp = nk1*nk2*nk3
          call info(20,0,0,' rdsigm: use full BZ (%i q-points)...',mxkp,0)
ckino Jan.4.2012           allocate(rv_p_oqp(abs(3*mxkp)))
          allocate(rv_a_oqp(abs(3*mxkp)))
          rv_p_oqp=>rv_a_oqp
          if (3*mxkp<0) rv_p_oqp(:)=0.0d0
          allocate(ipq_iv(mxkp))
          allocate(wgt_rv(abs(-mxkp)))
          if (-mxkp<0) wgt_rv(:)=0.0d0
          i_copy_size=size(slat%plat)
          call dcopy(i_copy_size,slat%plat,1,plat,1)

C   ... Make is,ifac,qb,qlat
          do   i = 1, 3
            llshft(i) = lshft(i) .ne. 0
          enddo
          call pshpr(0)
          call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
          call poppr
          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
     .      , 0 , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , 0 , 0 )
C       call pack2('bz nkp oqp',sbz,nqp,oqp)
          if (allocated(delt_rv)) deallocate(delt_rv)
          if (allocated(wgt_rv)) deallocate(wgt_rv)
        endif
C ... Write points to disk, shifting by qoffo
        if (lqoffo .ne. 0) then
C       call prmx('qp',w(oqp),3,3,nqp)
          do  i = 1, nqp
            call dmsadd ( rv_p_oqp , 1 , qoffo , 1 , 1 , 3 , 1 , 1 , 3 * 
     .       i - 2 , 1 , 1d0 )
          enddo
C       call prmx('qp',w(oqp),3,3,nqp)
        endif
        call info0(30,0,0,' RDSIGM: writing data to file QPTS ...')
        ifiz = fopna('qpts',-1,0)
        call getqp ( 1 , - ifiz , nqp , nkxyz , lshft , 0 , rv_p_oqp 
     .    , iwdummy , iwdummy )
        if (lfbzout) then
          if (allocated(ipq_iv)) deallocate(ipq_iv)
        endif
      endif
      call tcx('seneinterp')
      if (allocated(gstar_iv)) deallocate(gstar_iv)
      end subroutine seneinterp

!!--------------------------------------------------------------------------------------------
      subroutine rdsigm2(nbas_dummy,nsp,ndimh,slat,sham,sbz,ifis, 
     & nk1,nk2,nk3,ldim,qsmesh,sfz     ,qsfz,nqps) !,qsfz,qqps are for test
      use m_struc_def,only:      s_lat,s_ham,s_bz  
      use m_globalvariables,only: globalvariables
      use m_hamindex,only : getikt,napwk,plat
!! == Expand self-energy (read by ifis) to all the q point on mesh. In developing. ==
!! nbas is in this structure
!! input
!!    ifis:  file hundle for self-energy file sigm. only at irreducible q point.
!! output
!!   complex(8)::sfz(nk1,nk2,nk3,ndimh,ndimh,nsp):  self-energy (\Sigma-Vxc) for  all the q points on mesh.
!!   real(8):: qsmesh(3,nk1,nk2,nk3)
c  We have to clean up this routine. The purpose of this routine is "read sigm file and expand it in full BZ".
c  Not do more than that. (in future, we do scaling of simga in bndfp.F.
c  Many un-used local variables are contained.
c  Especially qsmesh (regular q mesh for self-energy.) is very problematic. It should be given at a place, and then
c  it should be used somewhere else.
c
c  original text below are not so meanigful. they are just for hints.
c
C- Read sigm(k) from file and generate sigm(R)
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do not read self-energy; just exit
Ci         :1 read and FT sigma(k) to make sigma(T)
Ci         :2 symmetrize sigma(T)
Ci         :4 take the real part of sigma(T) only
Ci         :10s digit (used here only for printout)
Ci         :0 Simple bloch sum of sigma
Ci         :1 approx high- and low-energy sigma with diagonal
Ci         :2 perturbation approach, diagonalizing exactly only
Ci         :  the first kappa block.
Ci         :3 Linear interpolation of sigma from neighboring p
Ci         :  sig(q) = sum_i wt_i A+_i sig(qpi) A_i
Ci         :  where A_i = z^LDA_qi (z_q^LDA)^-1
Ci         :100s digit specifies number of interpolation points
Ci         :    for interpolation mode 3.
Ci         :  0 -> use default (4 points)
Ci         :1000s digit specifies ASA conventions:
Ci         :      poke sigm(k)_RL,R'L' -> sigm(T)_RL,R'L'
Ci         :      (10000s digit of hft2rs)
Ci         :10000s digit
Ci         :  1 specifies file sigm stored with no group operations
Ci         :    sigm is stored at k-points in the full BZ.
Ci         :    Equivalent to --rsig:fbz
Ci         :  2 do not force R.S. neighbor table to conform to
Ci         :    symmetry
Ci         :  4 do not force file qp to match those generated by
Ci              bzmesh.  Use qp from bzmesh.
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ndimh :hamiltonian dimension
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat opos nsgrp oistab osymgr oag
Ci     Stored:
Ci     Passed to:
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: nkabc lshft nkp oqp
Ci         : (used only when rdsigm writes sigma on a new mesh)
Ci   ifis  :file logical unit for sigma
Ci   rsrnge:maximum length for connecting vectors in sigma(R)
Cio Inputs/Outputs:
Cio  sham  :struct for parameters defining hamiltonian; see routine uham
Cio    Elts read: ldham ooffH oindxo sigp rsstol
Cio    Stored:    nqsig oqsig ndhrs ontabs oiaxs ohrs
Cio Outputs:
Co   lwsig :0  No special transformations of sigma
Co         :1  Mode transforms sigma from orbital into LDA basis.
Co         :   rdsigm reads sigma from file sigm in orbital basis;
Co         :   Calling program generates and stores transformed sigma
Co         :2  Similar to lwsig=1, except
Co             low- and high- energy blocks replaced by diagonal parts
Co         :-1 Mode transforms sigm from LDA to orbital basis
Co         :   (Inverse operation of lwsig=1 or lwsig=2).
Co         :   This mode requires both sigma and LDA eigenvectors go be
Co         :   stored on disk in files 'sigm' and 'evec'
Co         :   rdsigm generates and stores transformation in file 'sigm2'
Co         :3  Returns lwsig=3 to flag calling program.  It should generate and
Co         :   store LDA eigenvalues and eigenvectors. No sigm file is read.
Co         :4  Returns lwsig=4 to flag calling program.  It should generate and
Co         :   store eigenvalues and eigenvectors.
Co         :5  Writes sigm as in lwsig=2 with low- and high- energy blocks
Co         :   replaced by diagonal parts, but writes sigm in orbital bas.
Cl Local variables
Cl   sigp  :parameters for approximating self-energy sigma.  sigma
Cl         :is approximated by its diagonal part sigii for energies
Cl         :below a low-energy cutoff (specified nmin or emin) and
Cl         :above a low-energy cutoff (specified nmax or emax).
Cl         : arg 1: specifies how to set diagonal part sigii
Cl         :        for states above the high-energy cutoff nmax or emax
Cl         :        0 constrain sigii to be > asig+bsig*e
Cl         :        1 constrain sigii to be = asig+bsig*e
Cl         :        2 constrain sigii to be > asig and < bsig
Cl         :        3 constraint same as case 1.
Cl         :          arg1=3 differs in that the least-squares fit to
Cl         :          sigii (for informational purposes only, to help
Cl         :          estimate asig and bsig) is done for states between
Cl         :          efit and nmax or emax
Cl         : arg 2: nmin : sigma for states 1..nmin are approximated by sigii
Cl         : arg 3: emin : (used only if nmin<0)
Cl         :             : sigma for levels e<emin are approximated by sigii
Cl         : arg 4: nmax : sigma for levels i>nmax are approximated by
Cl         :               sigii AND constrained according to mode
Cl         : arg 5: emax : (used only if nmax<=0)
Cl         :             : sigma for levels e<emax are approximated by
Cl         :               sigii AND constrained according to mode
Cl         : arg 6: asig : constraint used to approximate
Cl         :               sigii = asig + E * bsig  or
Cl         :               asig < sigii < bsig
Cl         : arg 7: bsig : constraint used to approximate
Cl         :               sigii = asig + E * bsig  or
Cl         :               asig < sigii < bsig
Cl         : arg 8: efit : (mode 3) energy minimium
Cl                         for fitting asig and bsig (not used here)
Cl   Note  :sigp takes a similar, but slightly different meaning when
Cl         :used in conjuction with sigm(LDA)->sigm(orbital), i.e. lwsig=-1
Cl         :In this case, sigp applies to the new basis, whose
Cl         :eigenvalues and eigenvectors are stored in 'evec.'
Cl         :Arguments are used in rotevs, with the following relations:
Cl         : arg1 (mode) has same meaning as before
Cl         : arg2 (nmin) takes meaning of lcore.
Cl         :       nmin>0 => sigm does not contain this block,  but the
Cl         :                 new basis does.  Use in emin (arg3) in
Cl         :                 place of sigm for this block .
Cl         :       nmin<0 => new basis does not contain this block
Cl         : arg3 (emin) If nmin>0, use emin for diag. sigma, this block
Cl         : arg4, arg5 (nmax,emax) could play the role of lhigh in
Cl         :       rotevs.  However, they are not used.  lhigh is
Cl         :       determined from constraint ndims+lccore+lhigh = ndimz
Cl
Cl   hreal :0 allow sigm(R) to be complex
Cl         :1 assume sigm(R) is real
Cl  ltrans :specifies indexing of sigm(T)_RL,R'L'L
Cl         : ltrans = F taken from sigm(k)_RL,R'L'
Cl         : ltrans = T taken from sigm(k)_RL,R'L'
Cl  kcplx  :0 sigm(k) has real, imaginary separated
Cl         :1 sigm(k) is in complex*16 format:
Cl         :2 sigm(k) has real, imaginary separated by columns
Cl  lphase :T if phase convention phi = q * [(g R_j + a) - R_i]
Cl         :  for rotations should be scaled by -1
Cl  lssym  :10000s digit mode
Cl  lfbzin :flags whether input self-energy file has suppressed symops
Cl  nscnts :compound of switches containing file contents of sigma
Cl  qoffi  :k-mesh offset for input sigma file
Cl  qoffo  :k-mesh offset for generated sigma file
Cr Remarks
Cb Bugs
Cb   Routine should be able to accomodate offset qp case (ifac<>0)
Cb   Routine should be able to accomodate case when a different
Cb   list of irreducible qp is used than the one generated by
Cb   bzmesh.
Cu Updates
Cu   27 Jun 08 Redesigned transformation of sigma to new basis
Cu   15 May 07 Parallelized symmetrizer
Cu   24 Jan 07 New option --wsig:onesp
Cu   20 Jan 07 Bug fix, --mixsig, spin polarized case
Cu   16 Jan 07 New option lwsig=-1; updated iosigh
Cu   24 Jul 06 MPI broadcast of sigma
Cu   24 Feb 05 Rotation of sigma matrix now in hrs
Cu             Switch lwsig to flag transformed sigm to be generated
Cu   27 Jan 05 New 40000s digit mode (allow qp mismatch)
Cu   20 Dec 04 New --wsig:rot and --wsig:phase switches
Cu             New 10000s digit mode
Cu   22 Nov 04 New --wsig:newkp switch
Cu   25 Sep 04 rdsigm can read sigm(q) file for offset q-mesh
Cu             rdsigm can write sigm(q) for mesh sbz->nkabc
Cu    5 Apr 04 bug fix for --wsig switch.  New --mixsig switch
Cu   15 Feb 04 Enable rdsigm to read a linear combination of
Cu             data from file ifis and data from file 'sigm1'
Cu    4 Jan 04 rdsigm checks bloch sum against ham->rsstol
Cu   10 Oct 03 rdsigm can spin-split a non-spin-polarized sigma
Cu   23 Sep 03 rdsigm can now read sigm fit to ASA hamiltonian
Cu   11 Jul 03 enable write of sigma to full BZ
Cu   24 May 03 Redesigned switches in accord with new interpolation
Cu   09 May 03 Added symmetrization of sigma
Cu   10 Jan 03 Some adaptations for Sergey's interpolation
Cu   14 Aug 02 Added option to orthogonalize sigm.
Cu   27 Jul 02 first created
C ------------------------------------------------------------
      implicit none
      integer mode,ifis,ndimh,lwsig, nbas_dummy!nbas,
      type(s_lat)::slat
      type(s_ham)::sham
      type(s_bz)::sbz
      double precision rsrnge
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nsp,nglob,ival,ldham(16),hreal,lrsig,
     .moditp,nttabs,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .ndhrs,j1,k1,k2,k3,iq1,nl,nspc,mxorb,nqsig,modsgp,nmin,nmax,
     .fopna,kcplx,lrot,iprint,stdo,lssym,ledit,niax,nscnts,
     .ndims,ndimz,iq,n123(4),lcore,lhigh
      parameter (niax=10)
c      equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
      integer:: ohrss , osigm2 , odelt , oistb2
      integer ,allocatable :: gstar_iv(:)
      integer ,allocatable :: ipq_iv(:)
      real(8) ,allocatable :: qp_rv(:)
      real(8) ,allocatable :: wgt_rv(:)
      complex(8) ,allocatable :: wk_zv(:)
      complex(8) ,allocatable :: sigm_zv(:,:)
      double precision tmp(3)
      real(8),allocatable:: evls(:),evlz(:),sigii(:)
      complex(8),allocatable:: siglda(:,:),z(:,:),sigo(:,:)
C     For offset q mesh
      integer is(3),lshft(3),ifac(3) !,lqoffo
      double precision rb(3,3),qb(3,3) !,qoffi(3),qoffo(3)
C MPI
c      logical mlog
c      integer procid,mpipid,master
C ... for mixing sigma files: Use alf1*sigm + alf2*sigm1
      double precision alf(2)
C     integer os
      double precision qp(3),tolq,rsstol,sigp(10)  !plat(3,3),qlat(3,3)
      double precision emin,emax,asig,bsig,efit,rotm(3,3),ddot !,qoff(3)
      equivalence (emin,sigp(3)),(emax,sigp(5))
      equivalence (asig,sigp(6)),(bsig,sigp(7)),(efit,sigp(8))
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
c      double precision qk
      integer jj1,jj2,jj3,k
C.....titus
C      integer ib,ib1,ib2,nlmaa,lidim,ldh,osfz1,nsp1
C      integer nkap0,n0H
C      parameter  (nkap0=3,n0H=5)
C      integer offH(n0H,nkap0,nbas)
C.....titus
      integer ::iwdummy
      real(8):: eseavr,qsmesh(3,nk1,nk2,nk3)
      integer:: i1,i2,i3,ikt,ldim,napw_in,debugmode
      complex(8)::sfz(nk1,nk2,nk3,ndimh,ndimh,nsp)
      integer:: i_copy_size
      logical:: isanrg, l_dummy_isanrg,debug=.false.


ccccccccccccccccccccccccccccccccc
      real(8):: qsfz(3,*)
      integer:: ii1,ii2,ii4
      sfz=1d99

C     --- Read sigma(orbital basis) from file ---
      call tcn('rdsigm2')
      print *
      print *,'rdsigm2:'
      mode=12
      lwsig=0
c      procid = mpipid(1)
c      master = 0
c      mlog = cmdopt('--mlog',6,0,outs)
      tolq = 1d-6
      stdo = globalvariables%stdo
      ledit = 0
      lrsig = 2  !mod(mode,10)
      lwsig = 0
      moditp = 1 !mod(mode/10,10)
      alf(1) = 1
      alf(2) = 0
      nscnts = 0
      lfbzout = .false.
      hreal = 0 
      ltrans=.true.
      kcplx = 0
      lphase = .true.
      lssym  = 0
      lonesp = 0
      lfbzin = .false.

c      if (procid .eq. master) then !   Require file contains sigm(orbital basis)
        call iosigh(0,nscnts,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2), lshft(3),ifis)
C         if (nscnts .ne. 0 .and. nscnts .ne. 5) then
C         call iosigh(2,0,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
C      .    lshft(3),ifis)
C         endif
c      endif
      lshft=0
c      call mpibc1(i,1,2,.false.,'rdsigm','nsp')
c      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
c      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
c      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
      if (i == nsp) then
        lsplts = .false.
      elseif (i > nsp) then
        call rx('rdsigm: sigm file spin polarized but calculation is not')
      else
        lsplts = .true.
      endif
      if (lsplts) call info0(30,0,0,'%9f(warning) sigm file not spin pol .. splitting spins')

c      if (procid .eq. master) then
        rewind ifis
c        if (lfbzin) nqp = nk1*nk2*nk3
        call iosigh(3,nscnts,i,ndimh,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
c      endif
c      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
c      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
c      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
      call info5(30,0,0,'%9fsigm file has %i irreducible QP: nk ='//
     .' ( %i %i %i )  shift=%3:1l',nqp,nk1,nk2,nk3,lshft)
      nsgrp=slat%nsgrp
      k1=nk1
      k2=nk2
      k3=nk3
      nsgrps = nsgrp

C ... Make is,ifac,qb,qlat,qoff
      do i = 1, 3
        llshft(i) = lshft(i) .ne. 0
      enddo
c      call pshpr(0)
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
c      call poppr

c.takao qmesh
      do i1=1,nk1
        do i2=1,nk2
          do i3=1,nk3
            qsmesh(:,i1,i2,i3) = (i1*ifac(1)-1)*qb(:,1) +
     &      (i2*ifac(2)-1)*qb(:,2) +
     &      (i3*ifac(3)-1)*qb(:,3)
          enddo
        enddo
      enddo
      if(debug) print *,'rdsigm2: debug111111'

C --- Get list of irreducible k-points, and ipq and gstar arrays ---
      mxkp = nk1*nk2*nk3
      if (associated(sham%rv_p_oqsig)) deallocate(sham%rv_p_oqsig)
      allocate(sham%rv_p_oqsig(abs(3*mxkp)))
      if (3*mxkp<0) sham%rv_p_oqsig(:)=0.0d0
      allocate(qp_rv(3*mxkp))
      allocate(gstar_iv(abs(-mxkp-1)))
      if (-mxkp-1<0) gstar_iv(:)=0
      gstar_iv = - 2
      allocate(ipq_iv(mxkp))
      allocate(wgt_rv(abs(-mxkp)))
      if (-mxkp<0) wgt_rv(:)=0.0d0
      if(debug) print *,'rdsigm2: debug2222'
      call info(20,1,0,' q-points in full BZ where sigma calculable ...',0,0)
      call bzmesh ( plat, qb, nk1,nk2,nk3, llshft,iwdummy,0, ipq_iv,sham%rv_p_oqsig, wgt_rv, nqsig, mxkp, 0, 0)
      sham%nqsig=nqsig
      call dpzero ( wgt_rv , mxkp )
      call info2(20,0,0,' Irr. qp for which sigma is calculated ...',0,0)
      call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_a_osymgr 
     . , nsgrps , ipq_iv , qp_rv , wgt_rv , nqps , mxkp , gstar_iv , 0 )
      l_dummy_isanrg=isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)

C --- Generate hrs = sigma(T) from file sigma(k) ---
C ... Loop over spins, keeping sig in full BZ for only one spin
c      if (procid .eq. master) then
        do  isp = 1, nsp
C     If sigma file not spin polarized, use sigma from spin 1
          if (isp .eq. 2 .and. lsplts) then
            call iosigh(3,0,1,ndimh,nk1,nk2,nk3,nqp, lshft(1),lshft(2),lshft(3),ifis)
          endif
C     If to mix in some other sigma file, set up ifis2
          if (alf(2) .ne. 0) then
            if (isp .eq. 1 .or. isp .eq. 2 .and. lsplts) then
              ifis2 = fopna('sigm1',-1,4)
              call iosigh(3,0,nsp,ndimh,nk1,nk2,nk3,nqp, lshft(1),lshft(2),lshft(3),ifis2)
            endif
          endif

          allocate(wk_zv(ndimh**2))
Changenglob          nl = nglob('nl')
          nl = globalvariables%nl
          nspc = 1
          allocate(sigm_zv(ndimh,ndimh))
          do  iq1 = 1, nqps
            read(ifis) qp,  eseavr !june2009  add eseavr
            print *,'readin sigm qp isp=',qp,isp

            if(debugmode()>0) print *,' uuu2 iq1 eseavr=',iq1,eseavr
            call dpscop(qp,tmp,3,1,1,1d0)
            call dpsadd ( tmp , qp_rv , 3 , 1 , 3 * iq1 - 2 , - 1d0 )
c        call dpsadd(tmp,qoffi,3,1,1,-1d0)
            !tmp= qp - w((oqp(3*(iq1-1)+1:3*(iq1-1)+3))+qoffi)
            if (.not. latvec(1,tolq,plat,tmp)) then
              call dpscop ( qp_rv , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 )
              if (lssym .ge. 4) then
                call info5(20,0,0,' rdsigm (warning) qp %i : '//
     .          'expected%3;8,4D read%3;8,4D',iq1,tmp,qp,0,0)
                call dpscop(tmp,qp,3,1,1,1d0)
              else
                print 456, iq1,tmp,qp
  456           format(' rdsigm: error on reading qp no',i4/
     .          ' Expected qp=',3f12.6/
     .          '     File qp=',3f12.6)
                call rx(' incompatible q-mesh')
              endif
            endif
            call dpdump ( sigm_zv , ndimh * * 2 * 2 , ifis )

c ... scalced sigm. alf(1:2)
            if (alf(1) .ne. 1d0) then
              call dscal ( ndimh * * 2 * 2 , alf ( 1 ) , sigm_zv , 1 )
            endif
c ... rotate sigm. w(osigm) at qp ---> w(ofbz) at stars of qp.
C       call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
            if(ldim==ndimh) then
              napw_in=0
              ikt=-9999
            else
              ikt = getikt(qp)
              napw_in= napwk(ikt)
            endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            ii1=mod(iq1,2)+1
c            ii2=mod(iq1/2,2)+1
c            ii4=mod(iq1/4,2)+1
c            print *,' iq1 ii4,ii2,ii1 =',iq1,ii4,ii2,ii1
c            sfz(ii4,ii2,ii1,1:ndimh,1:ndimh,isp)= sigm_zv(:,:)
c            qsfz(:,iq1)=qp
c            cycle
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            write(6,"(a,13f13.5)")' Goto hamfb3k  xxx input qp=', qp
            call hamfb3k ( qp , iq1 , nk1 , nk2 , nk3 , k1 , k2 , k3 , ipq_iv
     .      , napw_in , ndimh , ndimh , ndimh , qb , ldim , 
     .      ifac , gstar_iv , sigm_zv , sfz ( 1 , 1 , 1 , 1 , 1 , isp  ) )
            if(debugmode()>0) write(6,"(a,3f13.5)")'end of hamfbk3'
          enddo
          deallocate(sigm_zv)
          deallocate(wk_zv)
C ... End loop over spins
        enddo
c      endif ! exit of MPI master loop
      call tcx('rdsigm2')
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq_iv)) deallocate(ipq_iv)
      if (allocated(qp_rv)) deallocate(qp_rv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)
      end subroutine rdsigm2




