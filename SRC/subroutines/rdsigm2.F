      module m_rdsigm2
      public:: m_rdsigm2_init, getsenex, senex, dsene !getsenex returns the self-energy term.

      private
      complex(8),allocatable,protected:: sene(:,:),senex(:,:)
      real(8), allocatable,protected ::  rv_a_ohrs (:)
      integer, allocatable,protected :: iv_a_oiaxs(:)
      integer, allocatable,protected :: iv_a_ontabs(:)
      integer,protected:: ndimsig,nspsigm,ham_ndhrs,ham_nqsig
      integer,private :: nsizeohrs,nsizeoiaxs,nsizeontabs,noqsig
      real(8) , allocatable,private ::  rv_p_oqsig (:)  !not yet protected (rdsigm2 rewrite this)
      contains
      
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine getsenex(qp,isp,ndimh,ovlm)
!! this return self-energy senex, which is      
      implicit none
      integer:: isp,ndimh
      real(8):: qp(3)
      complex(8),allocatable:: ovlmtoi(:,:),ovliovl(:,:) 
      complex(8):: ovlm(ndimh,ndimh)
      call getsene(qp,isp) 
      allocate( ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
      ovlmtoi = ovlm(1:ndimsig,1:ndimsig)!,ispc)
      call matcinv(ndimsig,ovlmtoi)
      ovliovl = matmul(ovlmtoi,ovlm(1:ndimsig,1:ndimh))
      deallocate(ovlmtoi)
      allocate(senex(ndimh,ndimh))
      senex = matmul(transpose(dconjg(ovliovl)), matmul(sene,ovliovl))
      deallocate(ovliovl,sene)
      end subroutine
!sssssssssssssssssssssssssssssssssssssssss      
      subroutine dsene()
      deallocate(senex)
      end subroutine
!sssssssssssssssssssssssssssssssssssssssss      
      subroutine getsene(qp,isp) !(lrsig,qp,isp) !call after sigmainit
      use m_lmfinit,only: plat,mxorb,nbas,nl,iprmb
c      use m_suham,only:  iprmb
      integer:: i, mode, isp, ispsigm,iwdummy,ival!lrsig
      real(8):: hreal,qp(3)
      logical:: debug=.false.
      if(debug) print *,'ndimsig=',ndimsig
      allocate(sene(ndimsig,ndimsig))
c      i = lrsig*10
      hreal = 0
c      if (mod(i/10,10) >= 4) hreal = 1 !this is lrsig in hambles.
      mode = 100000 + 4000 + 40*(1-hreal) + 000 + 10 !hreal=1 -> sig(T) is real
      if(isp>nspsigm) then
         ispsigm = nspsigm
      else
         ispsigm=isp
      endif
!! bloch: we need hrs(means real space) and iaxs(pair table) for FFT.
      call bloch (mode,qp,nl,plat,mxorb,iprmb,1
     &           , ival ( iv_a_ontabs , nbas + 1 )
     &           , iv_a_oiaxs , rv_a_ohrs, ham_ndhrs, ispsigm, nspsigm , ndimsig , ndimsig 
     &           , 0 , ndimsig, 0 , ndimsig , 0 , sene, iwdummy, iwdummy )
               if(debug) print *,'end of bloch sene'
      end subroutine getsene
!sssssssssssssssssssssssssssssssssssssssss      
      
! takao expected to replace bloch to be a simple version seneinterp2. But not yet...
c     call seneinterp2(qp,nbas,nk1,nk2,nk3,plat, !npair,nlat,qwgt,npairmx,qsmesh,
c     i        sfz(:,:,:,:,:,isp),ldim,ndimsig, 
c     o        sene) 
c$$$C--- This works OK for lmf mode when noshrbz=.true. because no interpolation is needed.
c$$$C    Check comment out "sene=0" and "if(j1<=ldim.and.j2<=ldim) cycle" or not ---it depends on what you want to test.
c$$$c          allocate(senex(ndimh,ndimh))
c$$$          allocate(sene(ndimh,ndimh))
c$$$c          senex=0d0
c$$$          write(6,"(a,3f13.5)")'    qp=',qp
c$$$          do ik1=1,nk1
c$$$          do ik2=1,nk2
c$$$          do ik3=1,nk3
c$$$             write(6,"(a,3f13.5)")' qsmesh=',qsmesh(1:3,ik1,ik2,ik3)
c$$$             if(sum(abs(qp-qsmesh(1:3,ik1,ik2,ik3)))<1d-6) goto 1013
c$$$          enddo
c$$$          enddo
c$$$          enddo
c$$$          stop 'xxxxx not qp=qsmesh xxxx'
c$$$ 1013    continue
c$$$          do j1=1,ndimh !ldim+napw
c$$$          do j2=1,ndimh !ldim+napw
c$$$c            if(j1<=ldim.and.j2<=ldim) cycle
c$$$c            senex(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$            sene(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$          enddo
c$$$          enddo
c$$$c          print *,'ldim:abs and diff sene-senex=',sum(abs(sene(1:ldim,1:ldim)))
c$$$c     &      ,sum(abs(sene(1:ldim,1:ldim)-senex(1:ldim,1g:ldim)))
c$$$c          print *,'ndimh:abs and diff sene-senex=',sum(abs(sene))
c$$$c     &      ,sum(abs(sene-senex))
c$$$c          print *, 'end of test mode of sene newsigmasw=',newsigmasw()
c$$$c          deallocate(senex)
c$$$c          stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
c
c$$$C--- This works OK for lmf mode when noshrbz=.true. for PMT
c$$$C    Check comment out "sene=0" and "if(j1<=ldim.and.j2<=ldim) cycle" or not ---it depends on what you want to test.
c$$$c          allocate(senex(ndimh,ndimh))
c$$$          allocate(sene(ndimh,ndimh))
c$$$c          senex=0d0
c$$$          write(6,"(a,3f13.5)")'    qp=',qp
c$$$          do ik1=1,nqps
c$$$             write(6,"(a,i3,3f13.5)")' ik1 qsfz=',ik1,qsfz(1:3,ik1)
c$$$             iq1=ik1
c$$$            ii1=mod(iq1,2)+1
c$$$            ii2=mod(iq1/2,2)+1
c$$$            ii4=mod(iq1/4,2)+1
c$$$            print *,' iq1 ii4,ii2,ii1 =',iq1,ii4,ii2,ii1
c$$$             print *,' xxxxxxxxxxx this shows shorbz problem... xxxxxxxxx'
c$$$             if(sum(abs(qp-qsfz(1:3,ik1)))<1d-6) goto 1013
c$$$          enddo
c$$$          stop 'xxxxx not qp=qsmesh xxxx'
c$$$ 1013    continue
c$$$          do j1=1,ndimh !ldim+napw
c$$$          do j2=1,ndimh !ldim+napw
c$$$c            if(j1<=ldim.and.j2<=ldim) cycle
c$$$c            senex(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$            sene(j1,j2)= sfz(ii4,ii2,ii1,j1,j2,isp)
c$$$          enddo
c$$$          enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
!sssssssssssssssssssssssssssssssssssssssss      
      subroutine seneinterp(nbas,nsp,ndimh,ifis,rsrnge, nk1,nk2,nk3,sfz)
      use m_hamindex,only: offH
      use m_mksym,only: rv_a_osymgr ,iv_a_oistab,rv_a_oag,lat_nsgrp
      use m_hft2rs
      use m_sigflg
      use m_lmfinit,only: bz_lshft,ham_rsstol, mxorb,stdo,nl,ham_ldham,iprmb
      use m_mkqp,only: bz_nabc,bz_nkp
      use m_lattic,only: lat_plat,rv_a_opos
      use m_MPItk,only: master_mpi
      use m_ext,only:sname
      implicit none
!! seneinterp is used in hambls.F case (store real-space sigma to ham_rv_a_ohrs)
      integer mode,nbas,ifis,ndimh,lwsig
      double precision rsrnge
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nsp,nglob,ival,ldham(16),hreal,lrsig,
     .moditp,nttabs,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .ndhrs,j1,k1,k2,k3,iq1,nspc,nqsig,!modsgp,nmin,nmax,
     .kcplx,lrot,iprint,lssym,ledit,niax,nscnts,
     .ndims,ndimz,iq,n123(4),lcore,lhigh
      parameter (niax=10)
      integer:: ohrss
      real(8),allocatable :: rv_p_oqp(:)
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
      double precision plat(3,3),qp(3),tolq,rsstol,qlat(3,3),rotm(3,3),qoff(3),ddot
! Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      double precision qk
      integer jj1,jj2,jj3,k
      integer ::iwdummy
      integer:: napw_in
      complex(8)::sfz(nk1,nk2,nk3,ndimh,ndimh,nsp)
      real(8):: qin(3)
      integer:: ndhamx,debugmode
      complex(8):: sene(ndimh,ndimh),zv_dummy(1)
      real(8):: rv_dummy(1)
      integer:: i_copy_size,nhrss,sigswitch,ibas,ifile_handle
      real(8),allocatable::hrss(:)
      logical:: l_dummy_isanrg, isanrg, debug=.false.
!!
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) + (jj2*ifac(2)-1)*qb(k,2) +
     .(jj3*ifac(3)-1)*qb(k,3)
c
      call tcn('seneinterp')
      if(master_mpi) write(stdo,*)'seneinterp: startxxx'
      mode=12 !lrsig=12
      lwsig=0
      tolq = 1d-6
      ledit = 0
      lrsig = mod(mode,10)
      lwsig = 0
      moditp = mod(mode/10,10)
      nscnts = 0
      lqoffo = 0
      lfbzout = .false.
      hreal = 0 !real or complex
      if (lrsig .eq. 0) return
      ltrans = mod(mode/1000,10) .eq. 0
      kcplx = 1
      if (mod(mode/1000,10) .eq. 1) kcplx = 0
      lphase = mod(mode/1000,10) .eq. 0
      lssym = mod(mode/10000,10)
      l_dummy_isanrg=isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)
      
!! Command-line switches for input sigma and defaults
      call sigflg (dc= ' ' , sio='rsig' , strn=' ', lfbz=lfbzin  , lonesp=lonesp , qpoff=qoffi  )
      if (cmdopt('--rsig',6,0,outs)) then
        out2 = outs(7:)
        dc = out2(1:1)
        call sigflg ( dc=dc , sio='rsig' , strn=out2, lfbz=lfbzin  ,lonesp= lonesp , qpoff=qoffi  )
      endif
      if ( lfbzin .or. mod(lssym,2).eq.1 ) then
        lfbzin = .true.
      else
        lfbzin = .false.
      endif
      if (cmdopt('--mixsig=',9,0,outs)) then
         call rx('rdsigm: not support --mixsig anymore')
      endif
      call info(10,1,0,' RDSIGM: read file sigm and create '//
     .'%?#n#REAL#COMPLEX# sigma(R) by FT ...',hreal,0)
      ldham =ham_ldham
      rsstol=ham_rsstol

      call info2(30,0,0,
     .'%9fSigm will be approximated by:  '//
     .'%?#(n==0)#Simple bloch sum##%-1j'//
     .'%?#(n==1)#diagonal Sigma for high and low states##%-1j'//
     .'%?#(n==2)#Perturbation for higher kappa blocks##%-1j'//
     .'%?#(n==3)#Interpolation from known points##%-1j'//
     .' ',moditp,0)

      if (lfbzin .or. ddot(3,qoffi,1,qoffi,1).ne.0) then
        call info5(30,0,0,'%?#(n==1)#%9fFile sigm in saved in FBZ: '//
     .  'symmetrization suppressed##'//
     .  '%?#(n==1)#%N%9fFile k-mesh is offset by:%3:2g#%j#'//
     .  ' ',isw(lfbzin),isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi,0,0)
      endif

!! --- Read header of sigma(orbital basis) from file ---
!! ... Read k-mesh parameters
!!     Require file contains sigm(orbital basis)
      call iosigh(0,nscnts,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
      if (nscnts .ne. 0 .and. nscnts .ne. 5) then
         call iosigh(2,0,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
      endif
      if (i .eq. nsp) then
        lsplts = .false.
      elseif (i .gt. nsp) then
        call rx('rdsigm: sigm file spin polarized but calculation is not')
      else
        lsplts = .true.
      endif

      call info5(30,0,0,
     &   '%9fsigm file has %i irreducible QP: nk = ( %i %i %i )  shift=%3:1l',nqp,nk1,nk2,nk3,lshft)
      if (lsplts) call info0(30,0,0,
     &   '%9f(warning) sigm file not spin pol .. splitting spins') !sanity check
      if(debug) print *,'seneinterp: debug EEEEEEEEEE'
      plat = lat_plat
      nsgrp= lat_nsgrp
      call offshp ( iprmb , 1 , nbas , 2 , 0 , ndimh , 0 , ndhrs )
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
      if (mod(lssym,4) .ge. 2) i = 10*hreal + 0
      if (ltrans) i = i + 10000
      call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , 1 ,
     .    nsp , nbas , rv_a_osymgr , rv_a_oag , nsgrp , rsrnge , offH
     .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , rv_a_opos ,  iv_a_ontabs ,
     .    iv_a_oiaxs , ndhrs , rv_dummy )
      nttabs = ival ( iv_a_ontabs , nbas + 1 )
      if(debug) print *,'seneinterp: FFFFFFF size iv_p_ontabls ndhrs=',size(iv_a_ontabs),ndhrs,nttabs
      if(hreal .eq. 1) then
        if (allocated(rv_a_ohrs)) deallocate(rv_a_ohrs)
        allocate(rv_a_ohrs(ndhrs**2*nttabs*nsp))
        rv_a_ohrs(:)=0d0
      else
        if (allocated(rv_a_ohrs)) deallocate(rv_a_ohrs)
        allocate(rv_a_ohrs(2*ndhrs**2*nttabs*nsp))
        rv_a_ohrs(:)=0d0
      endif

C --- Get list of irreducible k-points, and ipq and gstar arrays ---
      mxkp = nk1*nk2*nk3
      if (allocated(rv_p_oqsig)) deallocate(rv_p_oqsig)
      allocate(rv_p_oqsig(abs(3*mxkp)))
      allocate(rv_p_oqp(abs(3*mxkp)))
      rv_p_oqp(:)=0d0
      allocate(gstar_iv(mxkp+1))
      gstar_iv(:) = 0
      gstar_iv = - 2
      allocate(ipq_iv(mxkp))
      allocate(wgt_rv(mxkp))
      wgt_rv=0d0
      call info(20,1,0,' q-points in full BZ where sigma calculable ...',0,0)
      call bzmesh(plat, qb, nk1, nk2, nk3, llshft, iwdummy, 0 , ipq_iv, rv_p_oqsig, wgt_rv, nqsig, mxkp, 0, 0 )
      ham_nqsig=nqsig
      wgt_rv=0d0
      call info2(20,0,0,' Irr. qp for which sigma is calculated ...',0,0)
      call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_a_osymgr 
     .  , nsgrps , ipq_iv , rv_p_oqp , wgt_rv , nqps , mxkp , gstar_iv, 0 )
      l_dummy_isanrg=isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
      call info2(30,0,0,'%?#(n==1)#%1fFile sigm has k-offset=%3:2g '//
     .'relative to the above#%j# ',isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi)

!! --- Generate hrs = sigma(T) from file sigma(k) ---
!! ... Loop over spins, keeping sig in full BZ for only one spin
!!     if (procid .eq. master) then
      do 1201 isp = 1, nsp
          if(debug) print *,' seneinterp:do 1201 start isp lsplts=',isp,lsplts
          if(debug) print *,'    ndimh nk1nk2nk3,nqp=', ndimh,nk1,nk2,nk3,nqp
          if(debug) print *,'    lshft,ifis=',lshft(1),lshft(2),lshft(3),ifis
!!     If sigma file not spin polarized, use sigma from spin 1
          if (isp .eq. 2 .and. lsplts) then
             rewind ifis
             read(ifis)
          endif
C ... FT sfz and copy sfz(T) to hrs
          i = 10*hreal + 1
          if (ltrans) i = i + 10000
          call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , isp
     .         , nsp , nbas , rv_a_osymgr , rv_a_oag , nsgrps , rsrnge , offH      
     .         , 1 , nbas , ndimh , ndimh , sfz ( 1 , 1 , 1 , 1 , 1 , isp )
     .         , plat , rv_a_opos ,iv_a_ontabs , iv_a_oiaxs , ndhrs , rv_a_ohrs )
          print *,' seneinterp:do 1201 end i ndimh ndhrs=',i,ndimh,ndhrs
 1201  continue
      if(debug) print *,' exit 1201 loop 111111111111111'
      if (allocated(sigm_zv)) deallocate(sigm_zv)
      if (allocated(wk_zv)) deallocate(wk_zv)
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq_iv)) deallocate(ipq_iv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)
      ham_ndhrs=ndhrs
      j1 = 0
      if (lssym .ge. 4) j1 = 2
      call dpzero(tmp,3)
      if(debug) print *,'a8 22222222222 goto chksgr'
!     ! --- chksgr read again sigm file---> but it is not for generating new quantities, just for check.
      call chksgr ( j1 , ltrans , kcplx , plat , nsp , ndimh , ifis
     .     , hreal , iprmb , nttabs , iv_a_oiaxs , rv_a_ohrs , ndhrs 
     .     , rsstol , i , rv_p_oqp , nbas , 0 , rotm , iwdummy )
      if(debug) print *,' a8 end of chksgr 3333333333'
      if ( .not. ( cmdopt ( '--wsig' , 6 , 0 , outs ) .or.cmdopt ( '-wsig' , 5 , 0 , outs ) ) ) then
         if (allocated(rv_p_oqp)) deallocate(rv_p_oqp)
      endif
      if(debug) print *,' a8 3333333333'
      
C --- Symmetrize hrs ---
C     Best to use nsgrp here instead of nsgrps.
C     However, neighbor table was generated for nsgrps, and may be
C     increase if nsgrp>nsgrps.  This branch will fail in that case.
      if (mod(lrsig,4) .ge. 2 .and. nsgrp .gt. 1) then
        nhrss= size(rv_a_ohrs)
        allocate(hrss(nhrss))
        hrss = rv_a_ohrs
C       10s digit distributes ri-rj pairs to avg (ri-rj),(rj-ri) pairs
C       Use 10*2 for hermitian matrices, 10*1 for symmetric ones
        i = 1-hreal + 10*2 + 100*(1-hreal) + 100000
        allocate(istb2_iv(nsgrp*nbas))
        call istbpm ( iv_a_oistab , nbas , nsgrp , istb2_iv )
        call rsmsym ( i , plat , mxorb , iprmb , ndimh , nbas ,
     .          rv_a_opos , nl , nsp , 1 , nttabs , iv_a_ontabs , iv_a_oiaxs
     .        , rv_a_osymgr , istb2_iv , nsgrp , ndhrs , hrss , rv_a_ohrs )
        if (allocated(istb2_iv)) deallocate(istb2_iv)
        deallocate(hrss)
C   ... Check how well symmetrized FT replicates sigma(k) at file qp
c        if (procid .eq. master) then
          rsstol = 0
            call chksgr ( 0 , ltrans , kcplx , plat , nsp , ndimh , ifis
     .          , hreal , iprmb , nttabs , iv_a_oiaxs , rv_a_ohrs , ndhrs 
     .          , rsstol , i , 0 , nbas , 0 , rotm , iwdummy )
C   ... Broadcast R.S. sigma after symmetrization
      endif
      if(debug) print *,' a8 6666666666666'
C ... Transformations done by caller: write qp list to disk
   99 continue
      if (lwsig .ne. 0) then
        call info0(30,0,0,' ')
        if (lnwmsh) then
           nkxyz=bz_nabc !sbz
           lshft=bz_lshft
           nqp=bz_nkp
          call info2(20,0,0,
     .    ' Transform sigma on new k-mesh:  '//
     .    'nk=%3:1i  shft=%3:1l',nkxyz,lshft)
        endif
        if (lfbzout) then
          mxkp = nk1*nk2*nk3
          call info(20,0,0,' rdsigm: use full BZ (%i q-points)...',mxkp,0)
          allocate(rv_p_oqp(abs(3*mxkp)))
          allocate(ipq_iv(mxkp))
          allocate(wgt_rv(mxkp))
          wgt_rv=0d0
          i_copy_size=size(lat_plat)
          call dcopy(i_copy_size,lat_plat,1,plat,1)
C   ... Make is,ifac,qb,qlat
          do   i = 1, 3
            llshft(i) = lshft(i) .ne. 0
          enddo
          call pshpr(0)
          call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
          call poppr
          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
     .      , 0 , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , 0 , 0 )
          if (allocated(delt_rv)) deallocate(delt_rv)
          if (allocated(wgt_rv)) deallocate(wgt_rv)
        endif
C ... Write points to disk, shifting by qoffo
        if (lqoffo .ne. 0) then
          do  i = 1, nqp
            call dmsadd ( rv_p_oqp , 1 , qoffo , 1 , 1 , 3 , 1 , 1 , 3 * 
     .       i - 2 , 1 , 1d0 )
          enddo
        endif
        if(master_mpi.and.iprint()>30) write(stdo,"(a)")' RDSIGM: writing data to file QPTS ...'
        ifiz=ifile_handle()
        open(ifiz,file='qpts.'//trim(sname))
        call getqp(1, -ifiz, nqp , nkxyz , lshft , 0 , rv_p_oqp, iwdummy , iwdummy )
        close(ifiz)
        if (lfbzout.and.allocated(ipq_iv)) deallocate(ipq_iv)
      endif
      nsizeohrs  = size( rv_a_ohrs )
      nsizeoiaxs = size( iv_a_oiaxs )
      nsizeontabs= size( iv_a_ontabs)
      noqsig = size(     rv_p_oqsig)
      call tcx('seneinterp')
      if (allocated(gstar_iv)) deallocate(gstar_iv)
      end subroutine seneinterp
      
!ssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine senebroadcast(procid,master,mlog)
      integer:: procid,master
      logical:: mlog
      call mpibc1(nsizeohrs,1,2,mlog,'bndfp','nsizeohrs')
      call mpibc1(nsizeoiaxs,1,2,mlog,'bndfp','nsizeoiaxs')
      call mpibc1(nsizeontabs,1,2,mlog,'bndfp','nsizeontabs')
      call mpibc1(noqsig,1,2,mlog,'bndfp','noqsig')
      if(procid/=master) then
         print *,'nnnnnn=',nsizeohrs,nsizeoiaxs,nsizeontabs,noqsig
         allocate( rv_a_ohrs(nsizeohrs) )
         allocate( iv_a_oiaxs(nsizeoiaxs))       
         allocate( iv_a_ontabs(nsizeontabs))       
         allocate( rv_p_oqsig(noqsig))
      endif   
      call mpibc1(iv_a_ontabs,nsizeontabs,2,mlog,'bndfp','ontabs')
      call mpibc1(iv_a_oiaxs ,nsizeoiaxs,2,mlog,'bndfp','oiaxs')
      call mpibc1(rv_a_ohrs,  nsizeohrs, 4,mlog,'bndfp','ohrs')
      call mpibc1(rv_p_oqsig, noqsig,4,mlog,'bndfp','noqsig')
      end subroutine senebroadcast

      
      
      
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine m_rdsigm2_init()
      use m_lmfinit,only : rsrnge=>ham_rsrnge,nbas,pwmode=>ham_pwmode
      use m_hamindex,only: symops_af,napwmx
      use m_MPItk,only: procid,master
      use m_ext,only:sname
      use m_lmfinit,only: ham_ldham
#if MPI | MPIK
      include "mpif.h"
#endif
      complex(8),allocatable:: sfz(:,:,:,:,:,:),wgtq(:,:,:,:,:),sfzx(:,:,:,:,:,:),sfzmto(:,:,:,:,:,:)
      integer:: ldim,ierr,ifi,ndimh_dummy,nk1,nk2,nk3,nqps,ifile_handle,ifis2,ik1,ik2,ik3,is
      logical:: laf,mlog,cmdopt, mtosigmaonly,cmdopt0
      character strn*120
      real(8),allocatable:: qsmesh2(:,:,:,:)
      ldim  = ham_ldham(1)
      ndimsig= ldim             !if(mtosigmaonly()) mode. Dimension of sigm is the size of MTOs.
      if(procid==master) then
         write(6,*)' --- goto sigmainit: rdsigm2 initialization section ---',rsrnge
         ifi=ifile_handle()
         open(ifi,file='sigm.'//trim(sname),form='unformatted') !
         read(ifi,err=9995,end=9995) nspsigm,ndimh_dummy,nk1,nk2,nk3
         laf=allocated(symops_af) !jun2015takao !reserved for future
         if(laf) nspsigm=2      ! we need recover laf mode
         allocate( qsmesh2(3,nk1,nk2,nk3) )
         if(mtosigmaonly()) then
            allocate(sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm))
         else
            call rx('mtosigmaonly()=T is needed in the current version sep2012')
         endif
         rewind ifi
         if (mod(pwmode,10) ==0) then
            write(6,"(a,2i5)")" --- ldim(=dim of lmto)= ",ldim
         elseif(napwmx>0) then 
            write(6,"(a,2i5)")" --- ldim(=dim of lmto), napwmx= ",ldim,napwmx
         endif   
!     ! sfz is the self-energy roteted in the full BZ.
         call rdsigm2(nbas,nspsigm,ifi,
     &        nk1,nk2,nk3,ldim,qsmesh2,sfz,  nqps, mtosigmaonly(),ndimsig) 
         close(ifi) 
!! --- write sig_fbz for SYMGRP e
c            if(.not.cmdopt('--wsig_fbz',10,0,strn).and.oncewrite(10)) then
         if(cmdopt0('--wsig_fbz')) then
           open(newunit=ifis2,file='sigm_fbz.'//trim(sname),form='unformatted')
           write(6,"(a)")' Writing sigm_fbz.* for SYMGRP e --wsig_fbz'
           write(ifis2) nspsigm,ndimsig,nk1,nk2,nk3,nk1*nk2*nk3,0,0,0
           do is=1,nspsigm
           do ik1=1,nk1
           do ik2=1,nk2
           do ik3=1,nk3
              write(ifis2) qsmesh2(1:3,ik1,ik2,ik3),is
              write(ifis2) sfz(ik1,ik2,ik3,1:ndimsig,1:ndimsig,is)
           enddo
           enddo
           enddo
           enddo
           close(ifis2) !call fclose(ifis2)
         endif
c
          
!! mesh generator.  Unused yet...->in future, we will use this.
c$$$            allocate(npair(nbas,nbas),pos(3,nbas))
c$$$            call dcopy ( nbas * 3 , rv_a_opos , 1 , pos , 1 )
c$$$            npairmx=nk1*nk2*nk3*2
c$$$            if (allocated(nlat)) deallocate(nlat)
c$$$            if (allocated(nqwgt)) deallocate(nqwgt)
c$$$            do !this loop is just to determined npairmx
c$$$              npairmx= npairmx + (nk1*nk2*nk3+1)*.5 !+1 added Feb2014 for 1x1x1
c$$$              allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
c$$$              call gennlat(pos,nbas,plat,nk1,nk2,nk3,npairmx,ok,npair,nlat,nqwgt)
c$$$              if(ok) exit
c$$$              deallocate( nlat, nqwgt )
c$$$            enddo
c$$$            deallocate(pos,npair)

!! -- We use ndimsig=ldim. This means sigm is projected on MTOs.
!!    The main input of the subroutine seneinterp is sfz. 
!!    The seneinterp stores real-space sigma to ham_rv_a_ohrs (real space representation), 
!!    which is used in the following 'call bloch' to obtain Sigm-Vxc for any q point by interpolation.
!! NOTE: sigm file is read again within seneinterp-chksgr, and compared with FFT-sfz.
            ifi = ifile_handle()
            open(ifi,file='sigm.'//trim(sname),form='unformatted')
            call seneinterp(nbas,nspsigm,ndimsig,ifi,rsrnge, 
     &        nk1,nk2,nk3,sfz)  !NOTE: output is stored in rv_a_ohrs in m_seneinterp
            close(ifi) 
            deallocate(sfz,qsmesh2) !,nlat,nqwgt)
      endif                  !procid==master
      if(cmdopt0('--wsig_fbz')) call rx0('end of --wsig_fbz mode')
#ifdef MPIK
      mlog     = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
      call senebroadcast(procid,master,mlog) !broadcast self-energy to all procid
      call mpibc1(ham_ndhrs,1,2,mlog,'bndfp','ndhrs')
      call mpibc1(nspsigm,1,2,mlog,'bndfp','nspsigm')
      call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
      if(procid==master) write(6,*)' --- end of sigmainit: rdsigm2 initialization section ---'
      return

 9995 continue
      call rx('sigmainit: readin error of sigm file')
      end subroutine m_rdsigm2_init
      
!!--------------------------------------------------------------------------------------------
      subroutine rdsigm2(nbas_dummy,nsp,ifis, 
     &     nk1,nk2,nk3,ldim,qsmesh,sfz     ,nqps,mtosigmaonly,ndimsig) !,qsfz,qqps are for test
      use m_mksym,only: rv_a_osymgr,lat_nsgrp
      use m_hamindex,only : getikt,napwk,plat,   symops_af,ngrp_original,ngrpaf,symops,ngrp
      use m_lmfinit,only: nl,stdo
!! == Expand self-energy (read by ifis) to all the q point on mesh. In developing. ==
!! nbas is in this structure
!! input
!!    ifis:  file hundle for self-energy file sigm. only at irreducible q point.
!! output
!!   complex(8)::sfz(nk1,nk2,nk3,ndimsig,ndimsig,nsp):  self-energy (\Sigma-Vxc) for  all the q points on mesh.
!!   real(8):: qsmesh(3,nk1,nk2,nk3)
!!  Self-energy (\Sigma-Vxc) is read from ifis file.
!!  It is stored into sigm_zv(ndimsig_r,ndimsig_r), which is rotated to be sfz in the full BZ by hamfb3k.
!!  ndimsig<=ndimsig_r 
!!   * ndimsig<=ndimsig_r is because I expect compatibility with current hqpe_sc where ndimsig_r= nlmto+max(napw)
!!     This should be corrected near future (written in 20sep2012).
!!
!!  We have to clean up this routine. The purpose of this routine is "read sigm file and expand it in full BZ".
!!  Not do more than that. (in future, we do scaling of simga in bndfp.F.
!!  Many un-used local variables are contained.
!!  Especially qsmesh (regular q mesh for self-energy.) is very problematic. It should be given at a place, and then
!!  it should be used somewhere else.


c---  original text below are not so meanigful. they are just for hints.
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
Cio    Elts read: ldham offH oindxo sigp rsstol
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
      integer mode,ifis,ndimsig_r,lwsig, nbas_dummy!nbas,
c      type(s_lat)::slat
c      type(s_ham)::sham
c      type(s_bz)::sbz
      double precision rsrnge
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nsp,nglob,ival,ldham(16),hreal,lrsig,
     .moditp,nttabs,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .ndhrs,j1,k1,k2,k3,iq1,nspc,mxorb,nqsig, !modsgp,nmin,nmax,
     .kcplx,lrot,iprint,lssym,ledit,niax,nscnts,
     .ndims,ndimz,iq,n123(4),lcore,lhigh
      parameter (niax=10)
c      equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
      integer:: ohrss , osigm2 , odelt , oistb2
      integer ,allocatable :: gstar_iv(:)
      integer ,allocatable,target :: ipq(:,:,:)
      real(8) ,allocatable :: qp_rv(:,:)
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
c      double precision alf(2)
C     integer os
      double precision qp(3),tolq,rsstol,rotm(3,3),ddot  !plat(3,3),qlat(3,3)
c      double precision emin,emax,asig,bsig,efit !,qoff(3)
c      equivalence (emin,sigp(3)),(emax,sigp(5)),sigp(10)
c      equivalence (asig,sigp(6)),(bsig,sigp(7)),(efit,sigp(8))
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
      real(8):: qsmesh(3,nk1,nk2,nk3) !eseavr,
      integer:: i1,i2,i3,ikt,ldim,napw_in,debugmode
      integer::  ndimsig !sep2012
      complex(8)::sfz(nk1,nk2,nk3,ndimsig,ndimsig,nsp)
      integer:: i_copy_size,ix
      logical:: isanrg, l_dummy_isanrg,debug=.false.,mtosigmaonly,laf

      real(8):: qir(3),diffq(3),platt(3,3)
      integer:: ii1,ii2,ii4,ispr,iaf,ig
      character(300)::aaa
      integer,allocatable,target:: ipqaf(:,:,:)
       integer,pointer:: ipq_pointer(:,:,:)
C     --- Read sigma(orbital basis) from file ---
      call tcn('rdsigm2')
      print *
      print *,'rdsigm2:'
      sfz=1d99
      laf=allocated(symops_af) !jun2015takao
      
      mode=12
      lwsig=0
c      procid = mpipid(1)
c      master = 0
c      mlog = cmdopt('--mlog',6,0,outs)
      tolq = 1d-6
c      stdo = globalvariables%stdo
      ledit = 0
      lrsig = 2  !mod(mode,10)
      lwsig = 0
      moditp = 1 !mod(mode/10,10)
c      alf(1) = 1
c      alf(2) = 0
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
c      print *,'iiiiiiiiiiii',i
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
         print *,'i,nsp=',i,nsp
        call rx('rdsigm: sigm file spin polarized but ctrl is not')
      else
        lsplts = .true.
      endif
      if (lsplts) call info0(30,0,0,'%9f(warning) sigm file not spin pol .. splitting spins')

c      if (procid .eq. master) then
        rewind ifis
c        if (lfbzin) nqp = nk1*nk2*nk3

!! takao add iosig optio=0 to get ndimsig_r.
        call iosigh(0,nscnts,i,ndimsig_r,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
!!

        call iosigh(3,nscnts,i,ndimsig_r,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
c      endif
c      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
c      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
c      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
      call info5(30,0,0,'%9fsigm file has %i irreducible QP: nk ='//
     .' ( %i %i %i )  shift=%3:1l',nqp,nk1,nk2,nk3,lshft)
      nsgrp=lat_nsgrp
      k1=nk1
      k2=nk2
      k3=nk3
      nsgrps = nsgrp

      print *,' lat_nsgrp=',lat_nsgrp
C ... Make is,ifac,qb,qlat,qoff
      do i = 1, 3
        llshft(i) = lshft(i) .ne. 0
      enddo
c      call pshpr(0)
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
c      call poppr


C --- Get list of irreducible k-points, and ipq and gstar arrays ---
      mxkp = nk1*nk2*nk3
      if (allocated(rv_p_oqsig)) deallocate(rv_p_oqsig)
      allocate(rv_p_oqsig(abs(3*mxkp)))
      allocate(qp_rv(3,mxkp))
      allocate(gstar_iv(mxkp+1))
      gstar_iv(:)=0
      gstar_iv = - 2
      allocate(ipq(nk1,nk2,nk3))
      allocate(wgt_rv(abs(-mxkp)))
      if (-mxkp<0) wgt_rv(:)=0.0d0
      if(debug) print *,'rdsigm2: debug2222'
      call info(20,1,0,' q-points in full BZ where sigma calculable ...',0,0)
      call bzmesh ( plat, qb, nk1,nk2,nk3, llshft,iwdummy,0, ipq,rv_p_oqsig, wgt_rv, nqsig, mxkp, 0, 0)
      ham_nqsig=nqsig
      call dpzero ( wgt_rv , mxkp )
      call info2(20,0,0,' Irr. qp for which sigma is calculated ...',0,0)
      call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_a_osymgr 
     . , nsgrps , ipq, qp_rv , wgt_rv , nqps , mxkp , gstar_iv , 0 )
      l_dummy_isanrg=isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)

      platt=transpose(plat)
!! qmesh
      do i1=1,nk1
        do i2=1,nk2
          do i3=1,nk3
            qsmesh(:,i1,i2,i3) = 
     &      (i1*ifac(1)-1)*qb(:,1) +
     &      (i2*ifac(2)-1)*qb(:,2) +
     &      (i3*ifac(3)-1)*qb(:,3)
          enddo
        enddo
      enddo
      if(debug) print *,'rdsigm2: debug111111'
      if(laf) then
        if(iprint()>10) write(6,*)'rdsimg2: AF mode, mapping from irr points to regular mesh point'
        allocate(ipqaf(nk1,nk2,nk3))
        ipqaf=0
        do i1=1,nk1
        do i2=1,nk2
        do i3=1,nk3
c          print *,'qmesh=',qsmesh(:,i1,i2,i3)
          do 1111 iq1=1,nqps
            qir = qp_rv(:,iq1)
c            print *,'qir=',qir
            do ig= ngrp_original+1,ngrp !only AF symmetry (equivalent with symops_af)
              call rangedq( matmul(platt,(qsmesh(:,i1,i2,i3) - matmul(symops(:,:,ig),qir))), diffq)
c              print *,'symops(:,:,ig) diffq',symops(:,:,ig),diffq
              if(sum(abs(diffq))<1d-6) then
               ipqaf(i1,i2,i3) = iq1    !iq1 is pointer to the irreducible q point = qp_rv(:,iq1)
               goto 1122
              endif
            enddo
 1111     continue
          write(aaa,"(3i5,3f13.5)") i1,i2,i3,qsmesh(:,i1,i2,i3)
          call rx('rdsigm2: 1111 loop can not find ipqaf'//trim(aaa))
 1122     continue
        enddo
        enddo
        enddo
      endif


C --- Generate hrs = sigma(T) from file sigma(k) ---
C ... Loop over spins, keeping sig in full BZ for only one spin
c      if (procid .eq. master) then
      do  isp = 1, nsp
C     If sigma file not spin polarized, use sigma from spin 1
          if (isp==2 .and. lsplts) then
            call iosigh(3,0,1,ndimsig_r,nk1,nk2,nk3,nqp, lshft(1),lshft(2),lshft(3),ifis)
          endif

c need to fix followings if necessary. Maybe it is better to clean iosigh(too complicated).
cC     If to mix in some other sigma file, set up ifis2
c          if (alf(2) .ne. 0) then
c            if (isp .eq. 1 .or. isp .eq. 2 .and. lsplts) then
c              ifis2 = fopna('sigm1',-1,4)
c              call iosigh(3,0,nsp,ndimsig_r,nk1,nk2,nk3,nqp, lshft(1),lshft(2),lshft(3),ifis2)
c            endif
c          endif

          allocate(wk_zv(ndimsig_r**2))
c          nl = globalvariables%nl
          nspc = 1
          allocate(sigm_zv(ndimsig_r,ndimsig_r))

!! we now not expect ordered data on qp feb2013
          do iq1 = 1, nqps
c            print *,' iq1 isp=',iq1,isp
!! look for a tag qp in sigm, where qp=qp_rv(:,iq1)  for given iq1
            do ix=0,1
              do 
                read(ifis,end=468) qp, ispr ! ispr is added dec2013
                read(ifis) sigm_zv
c                print *,'      qp ispr ix=',qp,ispr,ix,ndimsig_r
                if(laf) then
                  if(ispr==2) cycle 
                else
                  if(ispr/=isp) cycle
                endif  
                if(sum(abs(matmul(transpose(plat),qp-qp_rv(:,iq1))))<tolq) goto 460
c                   call dpdump ( sigm_zv , ndimsig_r * * 2 * 2 , ifis )
c                   call dpscop(qp,tmp,3,1,1,1d0)
c                   call dpsadd ( tmp , qp_rv , 3 , 1 , 3 * iq1 - 2 , - 1d0 )
c                   qdiff= qp - qp_rv(:,iq1)
cbug               if (latvec(1,tolq,plat,tmp)) goto 460               ! bug before dec2013
c                   if (latvec(1,tolq,plat, qp - qp_rv(:,iq1))) goto 460 ! bugfix isp==ispr dec2013
              enddo
 468          continue
              rewind ifis
              read(ifis)
            enddo
            write(aaa,"(i5,3f13.5)") iq1, qp_rv(:,iq1)
            call rx(' rdsigm2: read error. In sigm, we did not find iq= '//trim(aaa))
 460        continue

c ... scalced sigm. alf(1:2)
c            if (alf(1) .ne. 1d0) then
c               call rx('not use alf. use ScaledSigma in ctrl file')
cc              call dscal ( ndimsig_r * * 2 * 2 , alf ( 1 ) , sigm_zv , 1 )
c            endif
c ... rotate sigm. w(osigm) at qp ---> w(ofbz) at stars of qp.
C       call zprm('sigm',2,w(osigm),ndimsig_r,ndimsig_r,ndimsig_r)
            if(mtosigmaonly.or.ldim==ndimsig_r) then
              napw_in=0
              ikt=-9999
            else
              ikt = getikt(qp)
              napw_in= napwk(ikt)
            endif
c            print *,'ikt qp mtosigmalnly napw_in=',ikt,qp,mtosigmaonly, napw_in
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$            ii1=mod(iq1,8)+1
c$$$            ii2=mod(iq1/2,8)+1
c$$$            ii4=mod(iq1/4,8)+1
c$$$            print *,' iq1 ii4,ii2,ii1 =',iq1,ii4,ii2,ii1
c$$$c            sfz(ii4,ii2,ii1,1:ndimsig_r,1:ndimsig_r,isp)= sigm_zv(:,:)
c$$$            sfz(ii4,ii2,ii1,1:ndimsig,1:ndimsig,isp)= sigm_zv(:,:)
c$$$            qsfz(:,iq1)=qp
c$$$            cycle
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            iaf=0
            ipq_pointer => ipq
            if(laf.and.isp==1) then
              iaf=1
            elseif(laf.and.isp==2) then
              iaf=2
              ipq_pointer => ipqaf
            endif
            write(6,"(a,2i5,' ',13f13.5)")' rdsigm2:Goto hamfb3k  xxx input isp,iaf,qp=', isp,iaf,qp
            if(iprint()>60) write(6,"(a,13f13.5)")' rdsigm2:Goto hamfb3k  xxx input qp=', qp
            call hamfb3k ( qp , iq1 , nk1 , nk2 , nk3 , k1 , k2 , k3 , ipq_pointer,
     .       napw_in , ndimsig , ndimsig , ndimsig , qb , ldim , 
     .       ifac , gstar_iv , sigm_zv(1:ndimsig,1:ndimsig) , iaf, sfz(1,1,1,1,1,isp))
            if(debugmode()>0) write(6,"(a,3f13.5)")'end of hamfbk3'
          enddo
          deallocate(sigm_zv)
          deallocate(wk_zv)
C ... End loop over spins
        enddo
c      endif ! exit of MPI master loop
      call tcx('rdsigm2')
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq)) deallocate(ipq)
      if (allocated(qp_rv)) deallocate(qp_rv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)
      print *,'xxxxxxx end rdsigm2 xxxx'
      end subroutine rdsigm2
      end module m_rdsigm2



      


!!--- xxxx -----future version... hamfb3k, bloch will be replaced -----------------
      module m_seneinput_futureversion_unusedyet
      implicit none
      integer(4),allocatable:: npair(:,:), nlat(:,:,:,:)
      integer(4)::npairmx
      real(8),allocatable:: qsmesh(:,:,:,:)
      real(8),allocatable:: qwgt(:,:,:)
      contains
!! == Give \Sigma-Vxc for given qp from sfz (this does not work well)==
      subroutine seneinterp2(qp,nbas,nk1,nk2,nk3,plat, !npair,nlat,qwgt,npairmx,qsmesh,
     i  sfz,ldim,ndimx, 
     o  sene)
      use m_hamindex,only: napwmx, igv2,napwk,igv2rev,ibasindex,getikt
c      use m_seneinput,only: npair,nlat,qwgt,npairmx,qsmesh !they are read only
      implicit none
      integer(4):: nbas,nk1,nk2,nk3,ib1,ib2,ii,ik1,ik2,ik3,ldim,ndimx,
     &  ikc(3,0:1,0:1,0:1), igshift(3,0:1,0:1,0:1),
     &  j1x,j2x,ix1,ix2,ix3,ikt,ig,nkk(3),igvecc0(3),iktm,igvecc(3),j1,j2
c     &  npair(nbas,nbas),nlat(3,npairmx,nbas,nbas),
      complex(8):: sfz(nk1,nk2,nk3,ndimx,ndimx), 
     &  pi2i,sene(ndimx,ndimx),sss
      complex(8),allocatable:: wgtq(:,:,:,:,:)
      real(8)::rrr(3),plat(3,3) !,qsmesh(3,nk1,nk2,nk3)
      logical:: oncewrite,debug=.false.
      integer(4),allocatable::igx(:,:,:,:)
      real(8):: wgc(0:1,0:1,0:1),qpara(3),qmm(3),qqq(3),www,qp(3) !,qwgt(npairmx,nbas,nbas)
      call tcn('seneinterp2')
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

!!--- return when no APW part of self-energy.
      if(ldim==ndimx) then
        return
      endif

c      ikt=getikt(qp) ! index for qp
c      if(napwk(ikt)==0) then
c      if(oncewrite(4)) print *,' napw=0 skip napw part '
c        return
c      endif


cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c follwings are not so meaningful, although it returns correct Sigma on regular mesh points.
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
      allocate( igx(napwmx,0:1,0:1,0:1) )

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

      do j1=1,ndimx
        do j2=1,ndimx
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
      call tcx('seneinterp2')
      end subroutine seneinterp2
      end module 
