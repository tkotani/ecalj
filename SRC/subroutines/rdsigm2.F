      module m_rdsigm2
      public:: M_rdsigm2_init, Getsenex, Dsene, senex !getsenex returns the self-energy term.
      public:: ndimsig,sene
      private
      complex(8),allocatable,protected:: sene(:,:),senex(:,:)
      complex(8), allocatable,protected :: hhrs (:,:,:,:)
      integer, allocatable,protected :: iaxs(:,:)
      integer, allocatable,protected :: ntabs(:)
      integer,protected:: ndimsig,nspsigm,ndhrs,ham_nqsig
      integer,private :: nsizeohrs,nsizeoiaxs,nsizeontabs,noqsig,nttabs
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
      deallocate(ovliovl)
      end subroutine
!sssssssssssssssssssssssssssssssssssssssss      
      subroutine dsene()
      deallocate(senex,sene)
      end subroutine
!sssssssssssssssssssssssssssssssssssssssss      
      subroutine getsene(qp,isp) !(lrsig,qp,isp) !call after sigmainit
      use m_lmfinit,only: mxorb,nbas,nl,iprmb
      use m_lattic,only: plat=>lat_plat
      integer:: i,  isp, ispsigm,iwdummy!,ival!lrsig
      real(8):: qp(3)
      logical:: debug=.false.
      if(debug) print *,'ndimsig=',ndimsig
      allocate(sene(ndimsig,ndimsig))
      if(isp>nspsigm) then
         ispsigm = nspsigm
      else
         ispsigm=isp
      endif
!! bloch: we need hrs(means real space) and iaxs(pair table) for FFT.
      call bloch (qp,nl,plat,mxorb,iprmb,ntabs(nbas+1)
     &     , iaxs , hhrs, ndhrs, ispsigm, nspsigm ,
     & ndimsig , sene)
      if(debug) print *,'end of bloch sene'
      end subroutine getsene
      
!sssssssssssssssssssssssssssssssssssssssss      
      subroutine seneinterp(nbas,nspsigm,ndimsig,rsrnge, nk1,nk2,nk3,sfz)
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
!! seneinterp is used in hambls.F case (store real-space sigma to ham_hhrs)
      integer mode,nbas,ifis,ndimsig!,lwsig
      double precision rsrnge
      logical llshft(3),cmdopt,ltrans,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nspsigm,nglob,ldham(16),lrsig,
     .moditp,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .j1,k1,k2,k3,iq1,nspc,nqsig,!modsgp,nmin,nmax,
     .lrot,iprint,ledit,niax,nscnts,
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
      real(8),allocatable:: evls(:),evlz(:),sigii(:)
      complex(8),allocatable:: siglda(:,:),z(:,:),sigo(:,:)
C     For offset q mesh
      integer is(3),lshft(3),ifac(3)!,lqoffo
      double precision rb(3,3),qb(3,3),qoffi(3),qoffo(3)
      double precision plat(3,3),qp(3),tolq,rsstol,qlat(3,3),rotm(3,3),qoff(3),ddot
! Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      double precision qk
      integer jj1,jj2,jj3,k
      integer ::iwdummy
      integer:: napw_in
      complex(8)::sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm)
      real(8):: qin(3)
      integer:: ndhamx,debugmode
      complex(8):: sene(ndimsig,ndimsig),zv_dummy(1)
      real(8):: rv_dummy(1)
      integer:: i_copy_size,nhrss,sigswitch,ibas,ifile_handle
      real(8),allocatable::hrss(:)
      logical:: l_dummy_isanrg, isanrg, debug=.false.,cmdopt0
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) + (jj2*ifac(2)-1)*qb(k,2) +
     .(jj3*ifac(3)-1)*qb(k,3)
      call tcn('seneinterp')
c      open(ifis,file='sigm.'//trim(sname),form='unformatted')
      if(master_mpi) write(stdo,*)'seneinterp: startxxx'
      mode=12 !lrsig=12
      tolq = 1d-6
      ledit = 0
      lrsig = 2  ! mod(mode,10)
      moditp = 1 ! mod(mode/10,10)
      nscnts = 0
      ltrans = .true. !mod(mode/1000,10) .eq. 0
      ldham =ham_ldham
      rsstol=ham_rsstol
      plat = lat_plat
      nsgrp= lat_nsgrp
c     call offshp ( iprmb , 1 , nbas , 2 , 0 , ndimsig , 0 , ndhrs )
      ndhrs=ndimsig
      k1=nk1
      k2=nk2
      k3=nk3
C ... Number of group operations for input file sigma
      nsgrps = nsgrp
C ... Make is,ifac,qb,qlat,qoff
      do i = 1, 3
        llshft(i) = .false. !lshft(i) .ne. 0
      enddo
      call pshpr(0)
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
      call poppr
      if(debug) print *,'seneinterp: debug FFFFFFF 2222'
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
      write(stdo,*)' q-points in full BZ where sigma calculable ...'
      call bzmesh(plat,qb,nk1,nk2,nk3,llshft,iwdummy,0,ipq_iv,rv_p_oqsig,wgt_rv,nqsig,mxkp, 0, 0 )
      ham_nqsig=nqsig
      wgt_rv=0d0
      write(stdo,*)' Irr. qp for which sigma is calculated ...'
      call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_a_osymgr 
     .  , nsgrps , ipq_iv , rv_p_oqp , wgt_rv , nqps , mxkp , gstar_iv, 0 )
! Create pair table iax,ntab
      call hft2rsinit ( nk1 , nk2 , nk3 , k1 , k2 , k3 ,  1 ,
     .    nspsigm , nbas , rv_a_osymgr , rv_a_oag , nsgrp , rsrnge , offH
     .    , plat , rv_a_opos ,  ntabs ,
     .    iaxs,ndimsig)! , ndhrs , rv_dummy )
      nttabs = ntabs( nbas + 1 )
      if(debug) print *,'seneinterp: FFFFFFF size iv_p_ontabls ndhrs=',size(ntabs),ndhrs,nttabs
      if (allocated(hhrs)) deallocate(hhrs)
      allocate(hhrs(ndhrs,ndhrs,nspsigm,nttabs))
      hhrs=0d0
!! --- Generate hrs = sigma(T) from file sigma(k) ---
!! ... Loop over spins, keeping sig in full BZ for only one spin
!!     if (procid .eq. master) then
      do 1201 isp = 1, nspsigm
          if(debug) print *,' seneinterp:do 1201 start isp =',isp
          if(debug) print *,'    ndimsig nk1nk2nk3,nqp=', ndimsig,nk1,nk2,nk3,nqps
C ... FT sfz and copy sfz(T) to hrs (real space)
          call hft2rs ( nk1 , nk2 , nk3 , k1 , k2 , k3 , isp
     .         , nspsigm , nbas , rv_a_osymgr , rv_a_oag , nsgrps , rsrnge , offH      
     .         , 1 , nbas , ndimsig, sfz ( 1 , 1 , 1 , 1 , 1 , isp )
     .         , plat , rv_a_opos ,ntabs , iaxs , ndhrs , hhrs )
          print *,' seneinterp:do 1201 end i ndimsig ndhrs=',i,ndimsig,ndhrs
 1201  continue
      if(debug) print *,' exit 1201 loop 111111111111111'
      if (allocated(sigm_zv)) deallocate(sigm_zv)
      if (allocated(wk_zv)) deallocate(wk_zv)
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq_iv)) deallocate(ipq_iv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)
c      ham_ndhrs=ndhrs
      call tcx('seneinterp')
      end subroutine seneinterp
      
c      if (allocated(gstar_iv)) deallocate(gstar_iv)
c$$$      j1 = 0
c$$$      if(debug) print *,'a8 22222222222 goto chksgr'
c$$$!!--- chksgr read again sigm file---> but it is not for generating new quantities, just for check.
c$$$      call chksgr ( j1 , ltrans , plat , nspsigm , ndimh , ifis
c$$$     .     , iprmb , nttabs , iv_a_oiaxs , hhrs , ndhrs 
c$$$     .     , rsstol , i , rv_p_oqp , nbas , 0 , rotm , iwdummy )
c$$$      if(debug) print *,' a8 end of chksgr 3333333333'
c$$$C --- Symmetrize hrs ---
c$$$C     Best to use nsgrp here instead of nsgrps.
c$$$C     However, neighbor table was generated for nsgrps, and may be
c$$$C     increase if nsgrp>nsgrps.  This branch will fail in that case.
c$$$      if (mod(lrsig,4) .ge. 2 .and. nsgrp .gt. 1) then
c$$$        nhrss= size(hhrs)
c$$$        allocate(hrss(nhrss))
c$$$        hrss = hhrs
c$$$C       10s digit distributes ri-rj pairs to avg (ri-rj),(rj-ri) pairs
c$$$C       Use 10*2 for hermitian matrices, 10*1 for symmetric ones
c$$$        i = 1 + 10*2 + 100*1 + 100000
c$$$        allocate(istb2_iv(nsgrp*nbas))
c$$$        call istbpm ( iv_a_oistab , nbas , nsgrp , istb2_iv )
c$$$        call rsmsym ( i , plat , mxorb , iprmb , ndimh , nbas ,
c$$$     .          rv_a_opos , nl , nspsigm , 1 , nttabs , iv_a_ontabs , iv_a_oiaxs
c$$$     .        , rv_a_osymgr , istb2_iv , nsgrp , ndhrs , hrss , hhrs )
c$$$        if (allocated(istb2_iv)) deallocate(istb2_iv)
c$$$        deallocate(hrss)
c$$$C   ... Check how well symmetrized FT replicates sigma(k) at file qp
c$$$c        if (procid .eq. master) then
c$$$          rsstol = 0
c$$$            call chksgr ( 0 , ltrans ,  plat , nspsigm , ndimh , ifis
c$$$     .          , iprmb , nttabs , iv_a_oiaxs , hhrs , ndhrs 
c$$$     .          , rsstol , i , 0 , nbas , 0 , rotm , iwdummy )
c$$$C   ... Broadcast R.S. sigma after symmetrization
c$$$      endif
c$$$      if(debug) print *,' a8 6666666666666'
c$$$C ... Transformations done by caller: write qp list to disk
c$$$   99 continue
c      call tcx('seneinterp')
c      if (allocated(gstar_iv)) deallocate(gstar_iv)
c      close(ifis)
c      end subroutine seneinterp
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine m_rdsigm2_init()
      use m_lmfinit,only : rsrnge=>ham_rsrnge,nbas,pwmode=>ham_pwmode
      use m_hamindex,only: symops_af,napwmx
      use m_MPItk,only: procid,master
      use m_ext,only:sname
      use m_lmfinit,only: ham_ldham
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
!!    The seneinterp stores real-space sigma to ham_hhrs (real space representation), 
!!    which is used in the following 'call bloch' to obtain Sigm-Vxc for any q point by interpolation.
!! NOTE: sigm file is read again within seneinterp-chksgr, and compared with FFT-sfz.
            ifi = ifile_handle()
            call seneinterp(nbas,nspsigm,ndimsig,rsrnge, 
     &        nk1,nk2,nk3,sfz)  !NOTE: output is stored in hhrs in m_seneinterp
c            close(ifi) 
            deallocate(sfz,qsmesh2) !,nlat,nqwgt)
      endif                  !procid==master
      if(cmdopt0('--wsig_fbz')) call rx0('end of --wsig_fbz mode')
c#ifdef MPIK
      mlog     = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
      call senebroadcast(procid,master,mlog) !broadcast self-energy to all procid
      call mpibc1(ndhrs,1,2,mlog,'bndfp','ndhrs')
      call mpibc1(nspsigm,1,2,mlog,'bndfp','nspsigm')
c#endif
      if(procid==master) write(6,*)' --- end of sigmainit: rdsigm2 initialization section ---'
      return
 9995 continue
      call rx('sigmainit: readin error of sigm file')
      end subroutine m_rdsigm2_init
     
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine rdsigm2(nbas_dummy,nspsigm,ifis, 
     &     nk1,nk2,nk3,ldim,qsmesh,sfz     ,nqps,mtosigmaonly,ndimsig) !,qsfz,qqps are for test
      use m_mksym,only: rv_a_osymgr,lat_nsgrp
      use m_hamindex,only : getikt,napwk,   symops_af,ngrp_original,ngrpaf,symops,ngrp
      use m_lmfinit,only: nl,stdo
      use m_lattic,only: plat=>lat_plat
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
Cl         :xxx1 assume sigm(R) is real
Cl  ltrans :specifies indexing of sigm(T)_RL,R'L'L
Cl         : ltrans = F taken from sigm(k)_RL,R'L'
Cl         : ltrans = T taken from sigm(k)_RL,R'L'
Cl  sigm(k) is in complex*16 format:
Cl  lphase :T if phase convention phi = q * [(g R_j + a) - R_i]
Cl         :  for rotations should be scaled by -1
Cl  lssym  :10000s digit mode
Cl  lfbzin :flags whether input self-energy file has suppressed symops
Cl  nscnts :compound of switches containing file contents of sigma
Ccccl  qoffi  :k-mesh offset for input sigma file
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
      double precision rsrnge
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nspsigm,nglob,ldham(16),lrsig,
     .moditp,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .j1,k1,k2,k3,iq1,nspc,mxorb,nqsig, !modsgp,nmin,nmax,
     .lrot,iprint,lssym,ledit,nscnts,
     .ndims,ndimz,iq,n123(4),lcore,lhigh
      integer,parameter::niax=10
      integer:: ohrss , osigm2 , odelt , oistb2
      integer ,allocatable :: gstar_iv(:)
      integer ,allocatable,target :: ipq(:,:,:)
      real(8) ,allocatable :: qp_rv(:,:)
      real(8) ,allocatable :: wgt_rv(:)
      complex(8) ,allocatable :: wk_zv(:)
      complex(8) ,allocatable :: sigm_zv(:,:)
      real(8),allocatable:: evls(:),evlz(:),sigii(:)
      complex(8),allocatable:: siglda(:,:),z(:,:),sigo(:,:)
      integer is(3),lshft(3),ifac(3) !,lqoffo
      double precision rb(3,3),qb(3,3) !,qoffi(3),qoffo(3)
      double precision qp(3),tolq,rsstol,rotm(3,3),ddot  !plat(3,3),qlat(3,3)
      integer jj1,jj2,jj3,k
      integer ::iwdummy
      real(8):: qsmesh(3,nk1,nk2,nk3) !eseavr,
      integer:: i1,i2,i3,ikt,ldim,napw_in,debugmode
      integer::  ndimsig !sep2012
      complex(8)::sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm)
      integer:: i_copy_size,ix
      logical:: isanrg, l_dummy_isanrg,debug=.false.,mtosigmaonly,laf
      real(8):: qir(3),diffq(3),platt(3,3)
      integer:: ii1,ii2,ii4,ispr,iaf,ig,nsp_,ndimh_,nk1_,nk2_,nk3_,nqp_
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
      tolq = 1d-6
      ledit = 0
      lrsig = 2  !mod(mode,10)
      lwsig = 0
      moditp = 1 !mod(mode/10,10)
      nscnts = 0
      lfbzout = .false.
      ltrans=.true.
      lphase = .true.
      lssym  = 0
      lonesp = 0
      lfbzin = .false.
      lshft=0
c$$$      if (i == nsp) then
c$$$        lsplts = .false.
c$$$      elseif (i > nsp) then
c$$$         print *,'i,nsp=',i,nsp
c$$$        call rx('rdsigm: sigm file spin polarized but ctrl is not')
c$$$      else
c$$$        lsplts = .true.
c$$$      endif
c$$$      if (lsplts) call info0(30,0,0,'%9f(warning) sigm file not spin pol .. splitting spins')
!! takao add iosig optio=0 to get ndimsig_r.
      rewind ifis
      read(ifis) nsp_,ndimsig_r,nk1_,nk2_,nk3_,nqp_
c      call iosigh(0,nscnts,i,ndimsig_r,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
c      call iosigh(3,nscnts,i,ndimsig_r,nk1,nk2,nk3,nqp,lshft(1),lshft(2),lshft(3),ifis)
      write(stdo,"(' sigm file has ',i5,' irreducible QP: nk =',3i5)") nqp,nk1,nk2,nk3
      nsgrp=lat_nsgrp
      k1=nk1
      k2=nk2
      k3=nk3
      nsgrps = nsgrp
      print *,' lat_nsgrp=',lat_nsgrp
C ... Make is,ifac,qb,qlat,qoff
      do i = 1, 3
        llshft(i) = .false. !lshft(i) .ne. 0
      enddo
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
C --- Get list of irreducible k-points, and ipq and gstar arrays ---
      mxkp = nk1*nk2*nk3
      if (allocated(rv_p_oqsig)) deallocate(rv_p_oqsig)
      allocate(rv_p_oqsig(abs(3*mxkp)))
      allocate(qp_rv(3,mxkp))
      allocate(gstar_iv(mxkp+1))
      gstar_iv = -2
      allocate(ipq(nk1,nk2,nk3))
      allocate(wgt_rv(mxkp))
      wgt_rv(:)=0.0d0
      if(debug) print *,'rdsigm2: debug2222'
      call info(20,1,0,' q-points in full BZ where sigma calculable ...',0,0)
      call bzmesh(plat, qb, nk1,nk2,nk3, llshft,iwdummy,0, ipq,rv_p_oqsig, wgt_rv, nqsig, mxkp, 0, 0)
      ham_nqsig=nqsig
      call dpzero ( wgt_rv , mxkp )
      call info2(20,0,0,' Irr. qp for which sigma is calculated ...',0,0)
      call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_a_osymgr 
     . , nsgrps , ipq, qp_rv , wgt_rv , nqps , mxkp , gstar_iv , 0 )
      if(nqps/=nqp_) call rx('nqps/=nqp_ from sigm')
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
      do  isp = 1, nspsigm
         if (isp==2.and.nsp_==1 ) then !If sigm file not spin polarized, use sigm from spin 1
            rewind ifis
            read(ifis)
c            call iosigh(3,0,1,ndimsig_r,nk1,nk2,nk3,nqp, lshft(1),lshft(2),lshft(3),ifis)
          endif
          allocate(wk_zv(ndimsig_r**2))
          nspc = 1
          allocate(sigm_zv(ndimsig_r,ndimsig_r))
!! we now not expect ordered data on qp feb2013
          do iq1 = 1, nqps
!! look for a tag qp in sigm, where qp=qp_rv(:,iq1)  for given iq1
            do ix=0,1
              do 
                read(ifis,end=468) qp, ispr ! ispr is added dec2013
                read(ifis) sigm_zv
                if(laf) then
                  if(ispr==2) cycle 
                else
                  if(ispr/=isp) cycle
                endif  
                if(sum(abs(matmul(transpose(plat),qp-qp_rv(:,iq1))))<tolq) goto 460
              enddo
 468          continue
              rewind ifis
              read(ifis)
            enddo
            write(aaa,"(i5,3f13.5)") iq1, qp_rv(:,iq1)
            call rx(' rdsigm2: read error. In sigm, we did not find iq= '//trim(aaa))
 460        continue
            if(mtosigmaonly.or.ldim==ndimsig_r) then
              napw_in=0
              ikt=-9999
            else
              ikt = getikt(qp)
              napw_in= napwk(ikt)
            endif
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
      call tcx('rdsigm2')
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq)) deallocate(ipq)
      if (allocated(qp_rv)) deallocate(qp_rv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)
      print *,'xxxxxxx end rdsigm2 xxxx'
      end subroutine rdsigm2
!!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine hamfb3k(qin,iq,nk1,nk2,nk3,
c    i k1,k2,k3,ipq,napw_in,ndimh,ldima,ldimb,qb,plat,qlat,ldim,ifac,igstar,
     i k1,k2,k3,ipq,napw_in,ndimh,ldima,ldimb,qb,ldim,ifac,igstar,
     i hq,iaf, !iaf at jun2015
     o gfbz)
!! Generate gfbz (full BZ) from hq (1st BZ).
Ci q iq: q and index for q
Ci   nk1,nk2,nk3:  no. divisions for the qp in 3 recip. latt. vecs
Ci   k1,k2,k3: leading dimensions of gfbz
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci-temp  igstar:contains group operation needed to rotated a k-point
Ci          to its irreducible one (bzmesh.f) !this is just for check.
Ci   ndimh: dim of self energy hq at qin
Ci   ldima :dimensions gfbz; also the number of rows in gfbz to fill.
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   ldimb :dimensions gfbz; also the number of columns in gfbz to fill
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   hq    : Self-energy (or non-local potential for this iq
Ci read date in m_hamindex .
Co Outputs
Co   gfbz  : For those qp in star iq, hq stored
C ----------------------------------------------------------------------
      implicit none
      integer:: nk1,nk2,nk3,k1,k2,k3,ipq(*),igstar(0:*),ndimh,ldima,ldimb,napw_in,debugmode
      real(8)::    qin(3),qb(3,3) !,plat(3,3),qlat(3,3)
      complex(8):: hq(ndimh,ndimh),gfbz(k1,k2,k3,ldima,ldimb)
      integer:: i,i1,i2,i3,ig,iq,iq1,is,j,jj1,jj2,jj3,js,k,nl,ierr,ifac(3),j1,j2,ik1,ik2,ik3,isp,ldim,iaf
      real(8):: q1(3),qk
      character(200)::aaa
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .(jj2*ifac(2)-1)*qb(k,2) +
     .(jj3*ifac(3)-1)*qb(k,3)
      call tcn('hamfb3k')
      if(debugmode()>0) print *, 'hamfb3k: start...'
      iq1 = 0
      do  i3 = 1, nk3
        do  i2 = 1, nk2
          do  i1 = 1, nk1
            iq1 = iq1+1
            if(debugmode()>0) print *,'iq iq1 ipq(iq1)',iq,iq1,ipq(iq1)
!! ipq(iq1) ist gives a table to point irreducible point.
!!          q1,iq1 is target on regular mesh <--- qin,iq is irreducible points; this mapping is by rotsig.
!!          iq=ipq(iq1) shows iq for given iq1.
            if (ipq(iq1) .ne. iq) cycle !this must make things efficient
            q1(1) = qk(1,i1,i2,i3)
            q1(2) = qk(2,i1,i2,i3)
            q1(3) = qk(3,i1,i2,i3)
            if(debugmode()>0) print *,q1
c        ig = igstar(iq1)
c        write(6,"(a,i4)") 'symops',ig
c        write(6,"(3f10.5)")symops(1,:,ig)
c        write(6,"(3f10.5)")symops(2,:,ig)
c        write(6,"(3f10.5)")symops(3,:,ig)
            if(debugmode()>0) write(6,"(a,3f13.5)")' input          qin = ', qin
            if(debugmode()>0) write(6,"(a,3f13.5)")' target a        q1 = ', q1 ! qin = g(:,:,ig)^{-1}*q1
c        write(6,"(a,2i3, 3f13.5)")' target b   qp = ',iq1,ig, matmul(symops(:,:,ig),qin)-q1
c       write(6,"(a,3f13.5)")
            call rotsig(qin,q1,ndimh,napw_in,ldim,hq,gfbz(i1,i2,i3,:,:),ierr,iaf)
            if(ierr/=0) write(aaa,"(' qin=',3f13.6,' q1=',3f13.6)") qin,q1
            if(ierr/=0) call rx('hamfb3: rotsig do not map qin to q1;'//trim(aaa))
          enddo
        enddo
      enddo
      if(debugmode()>0) print *, 'hamfb3k: end...'
      call tcx('hamfb3k')
      end
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!! See fpgw/gwsrc/readeigen.F. 
      subroutine rotsig(qin,qout,ndimh,napw_in,ldim,sig,sigout,ierr,iaf) !iaf is added @jun2015takao
      use m_hamindex,only: symops,invgx,miat,tiat,shtvg,qlat,plat,dlmm,ngrp,norbmto,
     &  ibastab,ltab,ktab,offl,offlrev,getikt,igv2,igv2rev,napwk,
     &  ngrp_original
      implicit none
!!=== sigma rotator. sig at qin to sig at qout. ===
!! obtain sigout for qout.
!! a little confusing since qin=symops(qout), and basis rotation. Need to clean up.
!! Both of q and qtarget should be in qq table(in m_hamindex) which is given by gen_hamindex 
!! (read from QGpsi).
!! Used idea is   <base(qin)|sigma(qin)|base(qin)> = <g(base)|sigma |g(base)>.   
!! where qin=g(qout).  qtarget=qin= g(q=qout)
      integer   :: ig,ndimh,napw_in,ibaso,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr,igg,ikt2,ikt,l,ibas,ig2,k,ix
      real(8):: qin(3),qout(3)
      real(8)   :: q(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3) !plat(3,3),qlat(3,3)
      complex(8):: phase,img=(0d0,1d0),img2pi, sig(ndimh,ndimh),sigout(ndimh,ndimh)
      complex(8),allocatable:: sigx(:,:)
      integer :: ldim,debugmode,iaf,ngini,ngend
      character(300)::aaa
      img2pi=2*4d0*datan(1d0)*img
      ierr=1
      platt=transpose(plat) !this is inverse of qlat
!! find symops(3,3,ig),
      if(debugmode()>0) write(6,"('rotsig: qin qout=',3f9.4,x,3f9.4)") qin,qout
      qtarget= qin
      q      = qout  ! we find q
! we try to find qtrget = symops(igx) * q          (this means qin = symops(igc) qout).
!
      if(iaf==1) then !iaf mechanism is adde @jun2015takao for antiferro materials such as Nd2CuO4.
        ngini = 1
        ngend = ngrp_original
      elseif(iaf==2) then
        ngini = ngrp_original + 1
        ngend = ngrp
      else
        ngini = 1
        ngend = ngrp
      endif
!
      do igx=ngini,ngend
        if(debugmode()>0) print *, 'ddd=',matmul(platt,(qtarget-matmul(symops(:,:,igx),q)))
        call rangedq(   matmul(platt,(qtarget-matmul(symops(:,:,igx),q))), qx)
        if(sum(abs(qx))<1d-6) then
          igg=igx
          if(debugmode()>0) then
            print *,'ttt: q      =',q
            print *,'ttt: qtarget=',qtarget
            print *,'ttt: matmul q =',matmul(symops(:,:,igx),q)
            print *,'ttt: rotsig: OK! igg=',igg
            print *
          endif
          goto 1012
        endif
      enddo
      write(aaa,"(a,3f7.3,2x,3f7.3)")' rotsig: qtarget is not a star of q',q,qtarget
      call rx(trim(aaa))
      print *
      return
 1012 continue
      allocate(sigx(ndimh,ndimh))
      sigx=0d0
c      nlmto = ndimh-napw_in
      nlmto=ldim
      if(debugmode()>0) then
        print *,' tttt: invgx =',invgx(igg),shtvg(:,igg)
        print *,' tttt: ntorb napwin',norbmto,ndimh,napw_in,nlmto
      endif
!! mto part ---
c  write(6,"(3i3,2x,2i5,3x,a)") ib,l,k, offl(iorb,ib)+1,offl(iorb,ib)+2*l+1,trim(spid)
      if(nlmto/=0 )then
        ibaso=-999
        do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          !Index of Hamiltonian is devided into these blocks.
          ibas = ibastab(iorb)
          if(ibas/=ibaso) phase = exp( -img2pi*sum(qtarget*tiat(:,ibas,igg)) )
          ibaso=ibas
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+ 1
          iend1 = offl(iorb)+ 2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+ 1
          iend2 = offlrev(miat(ibas,igg),l,k)+ 2*l+1
          do ix=1,ndimh 
            sigx(ix,init1:iend1)= matmul(sig(ix,init2:iend2),dlmm(-l:l,-l:l,l,igg))*phase
          enddo
        enddo
      endif
!! apw part ------------
      if(napw_in/=0) then
        write(*,*) ' Probably OK-->Remove this stop to use this branch.'
     &   //' But need to confirm two apw sections in this routines.(phase factors) '
     &   //' Idea of this routine: <i|\sigma|j>_qout= <g(i)|\sigma|g(i)>_qin, where qin=g(qout)'
        stop 'abort'
        ikt  = getikt(q)    !index for q
        ikt2 = getikt(qtarget) !index for qtarget
        print *,' rotsig ikt ikt2=',ikt,ikt2
        if(napw_in /= napwk(ikt) ) then
          call rx('rotsig: napw_in /= napw(ikt)')
        endif
        do ig = 1,napw_in
          qpg = q + matmul( qlat(:,:),igv2(:,ig,ikt))      !q+G
          qpgr = matmul(symops(:,:,igg),qpg)               !rotated q+G
c          call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
          nnn=nint(matmul(platt,qpgr-qtarget))
          print *,ig,'nnn  ikt2=',nnn,ikt2
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2)
          phase= exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
          do ix=1,ndimh
            sigx(ix,nlmto+ig) = sig(ix,nlmto+ig2) * phase
          enddo
        enddo
      endif
      if(debugmode()>0) print *,' apw part end 111: ikt ikt2=',ikt,ikt2
!! mto part ------
      if(nlmto/=0) then
        ibaso=-999
        do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          !Index of Hamiltonian is devided into these blocks.
          ibas = ibastab(iorb)
          if(ibas/=ibaso) phase = exp( img2pi*sum(qtarget*tiat(:,ibas,igg)) )
          ibaso=ibas
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+ 1
          iend1 = offl(iorb)+ 2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+ 1
          iend2 = offlrev(miat(ibas,igg),l,k)+ 2*l+1
          do ix=1,ndimh
          sigout(init1:iend1,ix)= phase * matmul(transpose(dlmm(-l:l,-l:l,l,igg)),sigx(init2:iend2,ix))
          enddo
        enddo
      endif
      if(debugmode()>0) print *,' end of 2nd mto part q=',q
!! apw part ------------
      if(napw_in/=0) then
        ikt  = getikt(q)    !index for q
        ikt2 = getikt(qtarget) !index for qtarget
        if(debugmode()>0) print *,' rotsig 111 ikt ikt2=',ikt,ikt2
        if(napw_in /= napwk(ikt) ) then
          call rx('rotsig: napw_in /= napw(ikt)')
        endif
        do ig = 1,napw_in
          qpg = q + matmul( qlat(:,:),igv2(:,ig,ikt))      !q+G
          qpgr = matmul(symops(:,:,igg),qpg)               !rotated q+G
c          call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
          nnn=nint( matmul(platt,qpgr-qtarget))
          if(debugmode()>0) print *,ig,'nnn  ikt2=',nnn,ikt2
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2)
          phase=exp(img2pi*sum(qpgr*shtvg(:,igg))) 
          do ix=1,ndimh
            sigout(nlmto+ig,ix) =   sigx(nlmto+ig2,ix) * phase
          enddo
        enddo
        if(debugmode()>0) print *,' apw part end 222: ikt ikt2=',ikt,ikt2
      endif
      ierr=0
      if(debugmode()>0) print *,' goto deallcate sigx'
      deallocate(sigx)
      if(debugmode()>0) print *,' end of rotsig'
      end subroutine rotsig

!ssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine senebroadcast(procid,master,mlog)
      integer,parameter::niax=10
      integer:: procid,master
      logical:: mlog
      nsizeohrs  = size( hhrs )*2
      nsizeoiaxs = size( iaxs )
      nsizeontabs= size( ntabs)
      noqsig = size(     rv_p_oqsig)
      call mpibc1_int(ndhrs,1,'senebroadcast:ndhrs')
      call mpibc1_int(nspsigm,1,'senebroadcast:nspsigm')
      call mpibc1_int(nttabs,1,'senebroadcast:nttabs')
      
      call mpibc1(nsizeoiaxs,1,2,mlog,'bndfp','nsizeoiaxs')
      call mpibc1(nsizeontabs,1,2,mlog,'bndfp','nsizeontabs')
      call mpibc1(noqsig,1,2,mlog,'bndfp','noqsig')
      if(procid/=master) then
         print *,'nnnnnn=',nsizeoiaxs,nsizeontabs,noqsig,ndhrs,ndhrs,nspsigm,nttabs
         allocate(hhrs(ndhrs,ndhrs,nspsigm,nttabs))
         allocate( iaxs(niax,nsizeoiaxs/niax))       
         allocate( ntabs(nsizeontabs))       
         allocate( rv_p_oqsig(noqsig))
      endif   
      call mpibc1(ntabs,nsizeontabs,2,mlog,'bndfp','ntabs')
      call mpibc1_int(iaxs,     size(iaxs),'senebroadcast:iaxs')
      call mpibc1_complex(hhrs, size(hhrs),'senebroadcase:hhrs')
      call mpibc1(rv_p_oqsig, noqsig,4,mlog,'bndfp','noqsig')
      end subroutine senebroadcast
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      end module m_rdsigm2




      







      
c$$$!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c$$$      subroutine chksgr(opt,ltrans,plat,nsp,ndimh,ifis, !hreal, !kcplx
c$$$     .iprmb,nttabs,iaxs,hrs,ndhrs,tols,nqp,qplst,nbas,lrot,rotm,
c$$$     .delT)
c$$$      use m_lmfinit,only: nl,mxorb
c$$$      use m_hamindex,only : symops_af
c$$$C- Compare Bloch summed sigma against file sigma(k)
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   opt   :1s digit
c$$$Ci         :0 calculate sigma at list of qp in file,
c$$$Ci         :  compare against file sigma
c$$$Ci         :1 calculate sigma at list of qp, write to file sigma
c$$$Ci         :  NB: caller must write header to sigma file before
c$$$Ci         :  calling chksgr!
c$$$Ci         :2 like opt=0, but use qplst instead of qp
c$$$Ci   plat  :primitive lattice vectors, in units of alat
c$$$Ci   nsp   :2 for spin-polarized case, otherwise 1
c$$$Ci   ndimh :dimension of hamiltonian
c$$$Ci   ifis  :file logical unit for sigma
c$$$Ci   hreal :0 sigm(R) is complex
c$$$Cixxx         :1 sigm(R) is real
c$$$Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
c$$$Ci   nttabs:number of site pairs in sigm(R)
c$$$Ci   iaxs  :neighbor table for sigm(R)
c$$$Ci   hrs   :sigm(R)  <=================================== important input to be checked
c$$$Ci   ndhrs :leading dimensions of hrs
c$$$Ci   tols  :deviations in sigma > tols causes chksgr to abort
c$$$Ci         :tols=0 -> chksgr never aborts
c$$$Ci   qplst :(opt=0) not used
c$$$Ci         :(opt=1) list of qp at which to make sigma
c$$$Ci   lrot  :(opt=0) not used
c$$$Ci         :(opt=1) 0 -> do not rotate sigma matrix
c$$$Ci         :        1 -> rotate sigma matrix by rotm
c$$$Ci         :        2 -> additionally scale sigma matrix by phase
c$$$Ci         :             owing to change in shortening of basis
c$$$Ci         :             vectors
c$$$Ci   rotm  :(opt=0) not used
c$$$Ci         :(opt=1) list of qp at which to make sigma
c$$$Ci   delT  :(only used if lrot=2)
c$$$Ci         :change in lattice translation vectors on rotation
c$$$Ci         :to add phase shift to sigm.
c$$$Cio Inputs/Outputs
c$$$Cio  nqp   :Input, opt=0:  number of qp for which sigma was checked
c$$$Cio        :Output, opt=1: number of qp for which to make sigma
c$$$Cio  (opt=1) sigma written to logical unit ifis
c$$$Cl Local variables
c$$$Cl   sq    :(opt=0) file's contents of sigma, for comparison
c$$$Cl         :(opt=1) temporary array in rotating sigma
c$$$Cl   sq1   :bloch-summed sigma
c$$$Cr Remarks
c$$$Cr
c$$$Cu Updates
c$$$Cu   20 Dec 04 enable rotation of sigma before writing
c$$$Cu   11 Jul 03 enable write of sigma
c$$$Cu   09 May 03 First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer opt,ifis,ndimh,nsp,nttabs,ndhrs,nqp,lrot
c$$$      integer niax,n0,nkap0
c$$$      parameter (niax=10, n0=10, nkap0=3)
c$$$      integer iprmb(*),iaxs(niax,*),nbas
c$$$      double precision plat(3,3),tols,qplst(3,nqp),rotm(3,3),
c$$$     .delT(3,nbas)
c$$$      double complex hrs(ndhrs,ndhrs,nsp,nttabs)
c$$$      logical ltrans
c$$$C ... Local parameters
c$$$      logical lsplts
c$$$      integer i,j,isp,stdo,iq1,ipr,lgunit,isw,opt0,opt01,
c$$$     .nscnts
c$$$      integer ib,jb,is,norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
c$$$C     these are read by sigma header
c$$$      integer nk1,nk2,nk3,lshft(3),offh(nbas),ndim(nbas),offi,offj,
c$$$     .ndi,ndj
c$$$C     integer os
c$$$      double precision qp(3),errmxx,errmx,err,qpr(3),xx
c$$$      double complex hwk(ndhrs,ndhrs)
c$$$      complex(8),allocatable:: sq(:,:),sq1(:,:),uz(:,:)
c$$$      integer ::iwdummy
c$$$      integer ii,jj
c$$$      integer:: ndimhr !Read from sigm. This can be different from ndimh in the case of PMT, 
c$$$                       ! that is, ndimh< ndimhr=nlmto+napw !apr2012takao
c$$$      complex(8),allocatable:: sqr(:,:)
c$$$      integer:: lb1,lb2 !mar2014takao
c$$$      logical:: laf !jun2015takao
c$$$
c$$$      call getpr(ipr)
c$$$      laf=allocated(symops_af) !jun2015takao
c$$$      stdo = lgunit(1)
c$$$      opt01 = mod(mod(opt,10),2)
c$$$      opt0  = mod(opt,10)
c$$$      nscnts = 0
c$$$      allocate(sq(ndimh,ndimh),sq1(ndimh,ndimh))
c$$$C --- Rotate hrs if switch set ---
c$$$      if (opt01 .ne. 0 .and. mod(lrot,2) .ne. 0) then
c$$$        allocate(uz(ndimh,ndimh))
c$$$        call rothrm(2,ndimh,iprmb,rotm,1,nbas,ndimh,uz,uz)
c$$$        do  ib = 1, nbas
c$$$          call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,j,offl,ndim(ib))
c$$$          offh(ib) = j+1
c$$$        enddo
c$$$C       For each pair and spin, rotate hrs to U(ib) hrs U(jb)+
c$$$        do is = 1, nttabs
c$$$          ib = iaxs(1,is)
c$$$          jb = iaxs(2,is)
c$$$          offi = offh(ib)
c$$$          offj = offh(jb)
c$$$          ndi = ndim(ib)
c$$$          ndj = ndim(jb)
c$$$          do  isp = 1, nsp
c$$$C           U1 sigm(unrotated) U2+
c$$$            if (ltrans) then
c$$$C             ltrans: use (U1 h+ U2+)+ = (U1 (U2 h)+)+ = U2 h U1+
c$$$              call zgemm('N','N',ndj,ndi,ndj,(1D0,0D0),Uz(offj,offj),
c$$$     .        ndimh,hrs(1,1,isp,is),ndhrs,(0D0,0D0),hwk,ndhrs)
c$$$C             U sigm(unrotated) U+
c$$$              call zgemm('N','C',ndj,ndi,ndi,(1D0,0D0),hwk,ndhrs,
c$$$     .        Uz(offi,offi),ndimh,(0D0,0D0),hrs(1,1,isp,is),ndhrs)
c$$$            else
c$$$C             Make U1 h U2+
c$$$              call zgemm('N','N',ndi,ndj,ndi,(1D0,0D0),Uz(offi,offi),
c$$$     .        ndimh,hrs(1,1,isp,is),ndhrs,(0D0,0D0),hwk,ndhrs)
c$$$C             U sigm(unrotated) U+
c$$$              call zgemm('N','C',ndi,ndj,ndj,(1D0,0D0),hwk,ndhrs,
c$$$     .        Uz(offj,offj),ndimh,(0D0,0D0),hrs(1,1,isp,is),ndhrs)
c$$$            endif
c$$$
c$$$          enddo
c$$$        enddo
c$$$      endif
c$$$C ... Read header data
c$$$      lsplts = .false.
c$$$      if (opt01 .eq. 0) then
c$$$        call iosigh(0,nscnts,i,ndimhr,nk1,nk2,nk3,nqp,
c$$$     .  lshft(1),lshft(2),lshft(3),ifis)
c$$$        if (nscnts .ne. 0 .and. nscnts .ne. 5) then
c$$$          call iosigh(2,0,i,ndimhr,nk1,nk2,nk3,nqp,
c$$$     .    lshft(1),lshft(2),lshft(3),ifis)
c$$$        endif
c$$$        if (i .eq. nsp) then
c$$$          lsplts = .false.
c$$$        elseif (i .gt. nsp) then
c$$$          call rx(
c$$$     .    'rdsigm: sigm file spin polarized but calculation is not')
c$$$        else
c$$$          lsplts = .true.
c$$$        endif
c$$$      endif
c$$$      do  isp = 1, nsp
c$$$        if(laf.and.isp==2) cycle !jun2015takao
c$$$        errmxx = 0
c$$$C     If sigma file not spin polarized, use sigma from spin 1
c$$$        if (isp .eq. 2 .and. lsplts) then
c$$$          call iosigh(1,nscnts,1,ndimhr,nk1,nk2,nk3,nqp,
c$$$     .    lshft(1),lshft(2),lshft(3),ifis)
c$$$        endif
c$$$C ... For each q, compare Bloch sum of hrs against file's s(q)
c$$$C     or write to file
c$$$        do  iq1 = 1, nqp
c$$$          if (opt01 .eq. 0) then
c$$$            read(ifis) qp
c$$$            if (opt0 .eq. 2) call dcopy(3,qplst(1,iq1),1,qp,1)
c$$$            if (nscnts .eq. 5) call dpdump(xx,1,ifis)
c$$$            allocate(sqr(ndimhr,ndimhr))
c$$$            if(ndimhr<ndimh) call rx('chksgr :ndimhr<ndimh')
c$$$            call dpdump(sqr,ndimhr**2*2,ifis)
c$$$            sq = sqr(1:ndimh,1:ndimh) !!! readin sq
c$$$            deallocate(sqr)
c$$$            call ztoyy(sq,ndimh,ndimh,ndimh,ndimh,1,1) !kcplx,1)
c$$$          else
c$$$            call dcopy(3,qplst(1,iq1),1,qp,1)
c$$$          endif
c$$$C     Call bloch : 104010=perm orb, transpose, c*16
c$$$c          if(hreal/=0) call rx('going to bloch. We allow only hreal=0 but not')
c$$$c          i = 100000 + 4000 + 40*(1-hreal) + 10
c$$$          if (.not. ltrans) i = i - 4000
c$$$          call bloch (qp , nl , plat , mxorb , iprmb , 1 , nttabs 
c$$$     .     , iaxs , hrs , ndhrs , isp , nsp , ndimh , ndimh , 0 , ndimh 
c$$$     .     , 0 , ndimh , 0 , sq1 , iwdummy , iwdummy )
c$$$          errmx = 0
c$$$C       Compare Bloch sum of hrs against file's s(q)
c$$$          if (opt01 .eq. 0) then
c$$$            ii = 1
c$$$            jj = 1
c$$$            do  j = 1, ndimh
c$$$              do  i = 1, ndimh
c$$$                err = abs(sq(i,j)-sq1(i,j))  !Comparison here.
c$$$                if (err .gt. errmx) then
c$$$                  ii = i
c$$$                  jj = j
c$$$                  if (err .gt. tols .and. tols .gt. 0) then
c$$$                    print 357, i,j, sngl(err),sq1(i,j),sq(i,j)
c$$$  357               format(/
c$$$     .              '   i   j',6x,'diff',14x,'bloch sum',16x,'file value'/
c$$$     .              2i4,f12.6,2x,2f12.6,2x,2f12.6)
c$$$                    call rx1('rdsigm: Bloch sum deviates more than'//
c$$$     .              ' allowed tolerance (tol=%g)',tols)
c$$$                  endif
c$$$                endif
c$$$                errmx = max(errmx,err)
c$$$              enddo
c$$$            enddo
c$$$            if (ipr .ge. 45) write(stdo,332) iq1,errmx,ii,jj
c$$$  332       format(' comparing s(iq) to interpolated sq for iq=',i4,
c$$$     .      '  errmx=',1pe8.1:' at i,j=',2i4)
c$$$            errmxx = max(errmxx,errmx)
c$$$          else
c$$$C ...  Here we had 'Rotation of sigma(k)" but no longer needed, because hrs(T) was rotated now.
c$$$C   ... Write to s(q)
c$$$C         Additional scaling by phase shifts
c$$$            if (lrot .ge. 2) then
c$$$              call dgemm('T','N',3,1,3,1d0,plat,3,qp,3,0d0,qpr,3)
c$$$              call rothph(02,qpr,delT,ndimh,iprmb,1,nbas,ndimh,sq1)
c$$$            endif
c$$$C         Rotate qp to new coordinate system
c$$$            if (mod(lrot,2) .ne. 0) then
c$$$C           In-line multiply avoids bug in DEC fort compiler
c$$$              do  i = 1, 3
c$$$                qpr(i) = rotm(i,1)*qp(1) +
c$$$     .          rotm(i,2)*qp(2) +
c$$$     .          rotm(i,3)*qp(3)
c$$$              enddo
c$$$              call dcopy(3,qpr,1,qp,1)
c$$$            endif
c$$$C         Write qp, sigm for this qp
c$$$            write(ifis) qp
c$$$            call dpdump(sq1,ndimh**2*2,-ifis)
c$$$C         call zprm('sigm',2,sq1,ndimh,ndimh,ndimh)
c$$$          endif
c$$$        enddo
c$$$        if (opt01 .eq. 0)
c$$$     .  call info5(20,0,0,' check FT s(R) against s(q) ... '//
c$$$     .  'maximum error = %;2g%?#n>0# < tol (%;2g)#%j#%?#n==2# spin 2##',
c$$$     .  errmxx,isw(tols.ne.0),tols,isp,0)
c$$$      enddo
c$$$      deallocate(sq,sq1)
c$$$      if (opt01 .ne. 0 .and. mod(lrot,2) .ne. 0) deallocate(uz)
c$$$cccccccccccccccccccccccccc
c$$$      contains
c$$$      subroutine iosigh(optio,mode,nsp,ndimh,nk1,nk2,nk3,nqp, lshft1,lshft2,lshft3,ifi)
c$$$C- Read/write header information in sigma file
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   optio :(read only) 1s digit
c$$$Ci         :0 read header data into
c$$$Ci            nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
c$$$Ci         :1 File data must passed matched value for
c$$$Ci            nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
c$$$Ci         :2 file mode must passed matched value
c$$$Ci         :These digits may be taken in any combination
c$$$Ci   ifi   :File logical unit: ifi>0 for read; ifi<0 for write
c$$$Cio Inputs/Outputs
c$$$Cio  mode  : switches with information about the contents of sigma
c$$$Cio        :1s digit:
c$$$Cio        :0 file contains sigma in the orbital basis
c$$$Cio        :1 file contains sigma in the lda basis, and LDA
c$$$Cio           eigenvalues
c$$$Cio        :2 file contains sigma in the lda basis, and LDA
c$$$Cio           eigenvalues
c$$$Cio        :3 file consists of LDA eigenvectors and eigenvalues
c$$$Cio        :4 file consists of LDA eigenvectors and eigenvalues
c$$$Cio        :5 file contains sigma in the orbital basis, and LDA
c$$$Cio           eigenvalues
c$$$Cio  nsp   :2 for spin-polarized case, otherwise 1
c$$$Cio  ndimh :dimension of hamiltonian
c$$$Cio  nk1..nk3 : number of divisions in 1st,2nd,3rd Qlat
c$$$Cio  nqp   :number of irr. k-points
c$$$Cio  lshft1:1 if offset BZ in Qlat(1), otherwise 0
c$$$Cio  lshft2:1 if offset BZ in Qlat(2), otherwise 0
c$$$Cio  lshft3:1 if offset BZ in Qlat(3), otherwise 0
c$$$Cl Local variables
c$$$Cr Remarks
c$$$Cr   Reads or writes header file for sigma, with sanity checks
c$$$Cr   sigm file contents: (1s digit mode = 0..3,5,6)
c$$$Cr    1.  header
c$$$Cr    2.  qp offset (version 2 and later)
c$$$Cr    ... 3 and 4 are repeated for nqp points (nqp is contained in header)
c$$$Cr    3.  qp
c$$$Cr    4.  sigm(qp) or ham(qp)
c$$$Cr   sigm file contents: (1s digit mode = 1 or 3)
c$$$Cr    1.  header
c$$$Cr    2.  qp offset (version 2 and later)
c$$$Cr    ... 3..5 are repeated for nqp points (nqp is contained in header)
c$$$Cr    3.  qp
c$$$Cr    4.  eval(LDA)
c$$$Cr    5.  sigm(LDA) (mode 1) or evec(LDA) (mode 3)
c$$$Cu Updates
c$$$Cu   27 Jun 08 Revised mode, and input checking.  Old optio=1 -> optio=3
c$$$C ----------------------------------------------------------------------
c$$$      implicit none
c$$$C ... Passed parameters
c$$$      integer nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,ifi,optio, mode
c$$$C ... Local parameters
c$$$      integer jfi,i1,i2,i3,i4,i5,i6,i7,i8,i9,ivsn,modef
c$$$      logical:: isanrg, l_dummy_isanrg
c$$$      if (ifi .lt. 0) call rx('iosigh,ifi<0')
c$$$c$$$        ivsn = 1
c$$$c$$$C       if (mode .eq. 0) ivsn = 0
c$$$c$$$        jfi = -ifi
c$$$c$$$        rewind jfi
c$$$c$$$C       Write 0 => file contains version number (versions 1 and later)
c$$$c$$$        if (ivsn .ne. 0) write(jfi) 0, ivsn, mode
c$$$c$$$        write(jfi) nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
c$$$      if (optio .eq. 0) then
c$$$        rewind ifi
c$$$C       Read version, checking whether header is pre- version 1
c$$$        read(ifi,err=99,end=99) i1,ivsn,mode
c$$$C       File header is pre-version 1; set version to zero
c$$$        if (i1 .ne. 0) then
c$$$          ivsn = 0
c$$$          mode = 0
c$$$          rewind ifi
c$$$        endif
c$$$        read(ifi,err=99,end=99)
c$$$     .  nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
c$$$      else
c$$$        rewind ifi
c$$$C       Read version, checking whether header is pre- version 1
c$$$        read(ifi,err=99,end=99) i1,ivsn,modef
c$$$C       File header is pre-version 1; set version, modef to zero
c$$$        if (i1 .ne. 0) then
c$$$          ivsn = 0
c$$$          modef = 0
c$$$          rewind ifi
c$$$        endif
c$$$C       Assign mode or check with file correspondence
c$$$        if (mod(optio/2,2) .ne. 0) then
c$$$Ckino isanrg is logical function,           call isanrg(modef,mode,mode,'iosigh:','file''s mode',.true.)
c$$$          l_dummy_isanrg=isanrg(modef,mode,mode,'iosigh:','file''s mode',.true.)
c$$$        else
c$$$          mode = modef
c$$$        endif
c$$$C       Assign header parameters or check with file correspondence
c$$$        if (mod(optio,2) .ne. 0) then
c$$$          read(ifi,err=99,end=99) i1,i2,i3,i4,i5,i6,i7,i8,i9
c$$$          l_dummy_isanrg=isanrg(i1,nsp,nsp,'iosigh:','file''s nsp',.true.)
c$$$          l_dummy_isanrg=isanrg(i2,ndimh,ndimh,'iosigh:','file''s ndimh',.true.)
c$$$          l_dummy_isanrg=isanrg(i3,nk1,nk1,'iosigh:','file''s nk1',.true.)
c$$$          l_dummy_isanrg=isanrg(i4,nk2,nk2,'iosigh:','file''s nk2',.true.)
c$$$          l_dummy_isanrg=isanrg(i5,nk3,nk3,'iosigh:','file''s nk3',.true.)
c$$$          l_dummy_isanrg=isanrg(i6,nqp,nqp,'iosigh:','file''s nqp',.true.)
c$$$          l_dummy_isanrg=isanrg(i7,lshft1,lshft1,'iosigh:','file''s lshft1',.true.)
c$$$          l_dummy_isanrg=isanrg(i8,lshft2,lshft2,'iosigh:','file''s lshft2',.true.)
c$$$          l_dummy_isanrg=isanrg(i9,lshft3,lshft3,'iosigh:','file''s lshft3',.true.)
c$$$        else
c$$$          read(ifi,err=99,end=99)
c$$$     .    nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
c$$$        endif
c$$$      endif
c$$$      return
c$$$   99 continue
c$$$      call rx('iosigh: failed to read header from file')
c$$$      end subroutine iosigh
c$$$!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c$$$      subroutine rothph(opt,q,tbas,ldima,iprmb,ib1,ib2,ldH,H)
c$$$C- Extra scaling of hamiltonian by a phase
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   opt   :1s digit:
c$$$Ci         :0 scale on left by phase
c$$$Ci         :1 scale on right by conjugate of phase
c$$$Ci         :2 scale both left and right
c$$$Ci         :10s digit
c$$$Ci         :0 phase = exp(2*pi*tbas*q)
c$$$Ci         :1 phase = exp(-2*pi*tbas*q)
c$$$Ci   q     :k-point for current h
c$$$Ci         :Note: q and tbas must be in reciprocal units
c$$$Ci   tbas  :set of basis vector shifts; see Remarks
c$$$Ci         :Note: q and tbas must be in reciprocal units
c$$$Ci   ldima :points to last orbital in  downfolding block of h to rotate
c$$$Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
c$$$Ci   ib1   :H  generated for (ib1..ib2) site subblock of h
c$$$Ci   ib2   : -//-
c$$$Ci   ldh  :leading dimension of H
c$$$Co Outputs
c$$$Co   H    :H is scaled by a q-dependent phase; see Remarks
c$$$Cl Local variables
c$$$Cr Remarks
c$$$Cr   when the hamiltonian is constructed by Ewald sums, and a basis
c$$$Cr   vector is shortened by dpos, the hamiltonin is scaled by
c$$$Cr   a q-dependent phase.
c$$$Cr   hamiltonian rows belonging to that site are scaled by
c$$$Cr      phase  = exp(2*pi*i * q * dpos)
c$$$Cr   hamiltonian columns belonging to that site are scaled by
c$$$Cr      phase  = exp(-2*pi*i * q * dpos)
c$$$Cb Bugs
c$$$Cb   This routine has only been checked for ib1=1 and ib2=nbas
c$$$Cb   This routine assumes hamiltonian is contiguous block betw/ (ib1..ib2)
c$$$Cb   This routine uses off from orbl and assigns ldimp=0.
c$$$Cb     These assume only lower block is being rotated
c$$$Cu Updates
c$$$Cu   17 Dec 04 First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer opt,ldH,iprmb(1),ldima,ib1,ib2
c$$$      double precision q(3),tbas(3,ib2)
c$$$      double complex H(ldH,ldH)
c$$$C ... Local parameters
c$$$      logical nphase
c$$$      integer n0,nkap0,i,ib,j,j1,j2,ndimi,norb,off,ipr,PRDBG
c$$$      double precision twopi,sp
c$$$      double complex phase
c$$$      parameter (n0=10,nkap0=3,PRDBG=80)
c$$$      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
c$$$
c$$$      if (q(1) .eq. 0 .and. q(2) .eq. 0 .and. q(3) .eq. 0) return
c$$$
c$$$      twopi = 8d0*datan(1d0)
c$$$      nphase = .false.
c$$$      call getpr(ipr)
c$$$      if (mod(opt/10,10) .ne. 0) nphase = .true.
c$$$
c$$$C ... Offsets to starting and ending regions of h(ib1..ib2)
c$$$      call orbl(ib1,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)
c$$$      j1 = off+1
c$$$      call orbl(ib2,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)
c$$$      j2 = off + ndimi
c$$$
c$$$C --- For each site ibas, construct poke rotation matrix rotating Ylm(ibas)
c$$$      do  ib = ib1, ib2
c$$$
c$$$C   ... Row, column offsets to h
c$$$        call orbl(ib,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)
c$$$
c$$$C   ... Phase
c$$$        sp = twopi*(tbas(1,ib)*q(1) + tbas(2,ib)*q(2) + tbas(3,ib)*q(3))
c$$$        if (nphase) sp = -sp
c$$$        phase = dcmplx(dcos(sp),dsin(sp))
c$$$        if (ipr .ge. PRDBG) then
c$$$          print 333, ib,sp*360/twopi,q
c$$$  333     format(' rothph: ib=',i4,'  phase=',f7.1,' deg.  q=',3f7.1)
c$$$        endif
c$$$
c$$$C   --- For all orbitals belonging to site ib, h->phase*h and/or h->h*phase+
c$$$        if (mod(opt,10) .eq. 0 .or. mod(opt,10) .eq. 2) then
c$$$          do  j = j1, j2
c$$$            do  i = off+1, off+ndimi
c$$$              h(i,j) = phase*h(i,j)
c$$$            enddo
c$$$          enddo
c$$$        endif
c$$$        if (mod(opt,10) .eq. 1 .or. mod(opt,10) .eq. 2) then
c$$$          phase = dconjg(phase)
c$$$          do  i = off+1, off+ndimi
c$$$            do  j = j1, j2
c$$$              h(j,i) = h(j,i)*phase
c$$$            enddo
c$$$          enddo
c$$$        endif
c$$$      enddo
c$$$      end subroutine rothph
c$$$cccccccccccccccccccccccccccccccccc      
c$$$      end subroutine chksgr
c$$$!!sssssssssssssssssssssssssssssssssssss
c$$$      subroutine rothrm(opt,ldima,iprmb,rotm,ib1,ib2,ldU,U,Uz)
c$$$C- Generate rotation matrix for hamiltonian
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   opt   :1s digit
c$$$Ci         :  0 use real harmonics; output U is real
c$$$Ci         :  1 use true spherical harmonics; output U is complex
c$$$Ci         :  2 use real harmonics, but return U as complex
c$$$Cixxx   ldmpa :offset to first orbital in downfolding block of h to rotate
c$$$Ci   ldima :points to last orbital in  downfolding block of h to rotate
c$$$Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
c$$$Ci   rotm  :3x3 cartesian rotation matrix
c$$$Ci   ib1   :U  generated for (ib1..ib2) site subblock of h
c$$$Ci   ib2   : -//-
c$$$Ci   ldU  :leading dimension of U
c$$$Co Outputs
c$$$Co   U,Uz :Rotation matrix for (ib1..ib2) subblock of h
c$$$Co        :U is the real form (opt=0); Uz is the complex form.
c$$$Co        :Only of of (U,Uz) is returned.
c$$$Cl Local variables
c$$$Cl   rYL  :rotation matrix that rotates Y_lm
c$$$Cr Remarks
c$$$Cr   This routine generates rotation matrix U for the hamiltonian matrix
c$$$Cr   which consists of blocks which transform in the row and colum indices
c$$$Cr   as Y_lm.
c$$$Cr
c$$$Cr   The unrotated and rotated hamiltonians h and hr are related by :
c$$$Cr      hr = U h U+
c$$$Cr
c$$$Cr   See also rothph, below
c$$$Cb Bugs
c$$$Cb   This routine has only been checked for ib1=1 and ib2=nbas
c$$$Cu Updates
c$$$Cu   17 Dec 04 First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer opt,ldU,iprmb(1),ldmpa,ldima,ib1,ib2
c$$$      double precision rotm(3,3),U(ldU,ldU)
c$$$      double complex Uz(ldU,ldU)
c$$$C ... Local parameters
c$$$      integer n0,nkap0,nlx,nl2,i,ibas,iorb,j,li,
c$$$     .nlm,nlmi,norb,offri,offi
c$$$      parameter (n0=10,nkap0=3,nlx=9,nl2=nlx*nlx)
c$$$      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
c$$$      double precision xx,rYL(nl2,nl2)
c$$$      double complex rYLz(nl2,nl2)
c$$$
c$$$      if (opt .ne. 0 .and. opt .ne. 2)
c$$$     .call rx('rothrm not ready for this opt')
c$$$
c$$$C ... Rotation matrix for spherical harmonics
c$$$      call ylmrtg(nl2,rotm,rYL)
c$$$C     call prmx('rothrm: rYL',rYL,nl2,nl2,nl2)
c$$$
c$$$      if (opt .eq. 1) then
c$$$        call ylmrtg(nl2,rotm,rYL)
c$$$        call s2sph(1,nlx,nlx,rYL,nl2,nl2,nl2,nl2,rYLz)
c$$$C       call yprm('rYL',3,rYL,nl2*nl2,nl2,nl2,nl2)
c$$$      endif
c$$$
c$$$C ... Initialize U
c$$$      if (opt .eq. 0) then
c$$$        call dpzero(U,ldU*ldU)
c$$$      else
c$$$        call dpzero(Uz,ldU*ldU*2)
c$$$      endif
c$$$
c$$$C --- For each (ibas,l) block, poke rotation matrix rotating Ylm(block)
c$$$      do  ibas = ib1, ib2
c$$$
c$$$C   ... Row, column offsets to h
c$$$        call orbl(ibas,0,ldima,iprmb,norb,ltab,ktab,xx,offl,nlm)
c$$$
c$$$C   --- For each (ibas,l) block, poke rYL(i,j) into U(i,j) ---
c$$$C   ... loop over pairs of orbital blocks
c$$$        do  iorb = 1, norb
c$$$          li = ltab(iorb)
c$$$          nlmi = 2*li + 1
c$$$          offi = offl(iorb) !- ldmpa
c$$$          offri = li**2
c$$$          if (li+1 .gt. nlx)
c$$$     .    call rxi('rothrm: increase nlx, need nlx=',(li+1)**2)
c$$$
c$$$C     ... loop over pairs of m channels within this block
c$$$          if (opt .eq. 0) then
c$$$            do  j = 1, nlmi
c$$$              do  i = 1, nlmi
c$$$                U(i+offi,j+offi) = rYL(i+offri,j+offri)
c$$$              enddo
c$$$            enddo
c$$$          elseif (opt .eq. 1) then
c$$$            do  j = 1, nlmi
c$$$              do  i = 1, nlmi
c$$$                Uz(i+offi,j+offi) = rYLz(i+offri,j+offri)
c$$$              enddo
c$$$            enddo
c$$$          else
c$$$            do  j = 1, nlmi
c$$$              do  i = 1, nlmi
c$$$                Uz(i+offi,j+offi) = rYL(i+offri,j+offri)
c$$$              enddo
c$$$            enddo
c$$$          endif
c$$$        enddo
c$$$      enddo
c$$$C      if (opt .eq. 0) then
c$$$C        call prmx('U',U,ldU,ldU,ldU)
c$$$C      else
c$$$C        call zprm('Uz',2,Uz,ldU,ldU,ldU)
c$$$C      endif
c$$$      end subroutine rothrm
c$$$
c$$$      
c$$$      
c$$$! takao expected to replace bloch to be a simple version seneinterp2. But not yet...
c$$$c     call seneinterp2(qp,nbas,nk1,nk2,nk3,plat, !npair,nlat,qwgt,npairmx,qsmesh,
c$$$c     i        sfz(:,:,:,:,:,isp),ldim,ndimsig, 
c$$$c     o        sene) 
c$$$c$$$C--- This works OK for lmf mode when noshrbz=.true. because no interpolation is needed.
c$$$c$$$C    Check comment out "sene=0" and "if(j1<=ldim.and.j2<=ldim) cycle" or not ---it depends on what you want to test.
c$$$c$$$c          allocate(senex(ndimh,ndimh))
c$$$c$$$          allocate(sene(ndimh,ndimh))
c$$$c$$$c          senex=0d0
c$$$c$$$          write(6,"(a,3f13.5)")'    qp=',qp
c$$$c$$$          do ik1=1,nk1
c$$$c$$$          do ik2=1,nk2
c$$$c$$$          do ik3=1,nk3
c$$$c$$$             write(6,"(a,3f13.5)")' qsmesh=',qsmesh(1:3,ik1,ik2,ik3)
c$$$c$$$             if(sum(abs(qp-qsmesh(1:3,ik1,ik2,ik3)))<1d-6) goto 1013
c$$$c$$$          enddo
c$$$c$$$          enddo
c$$$c$$$          enddo
c$$$c$$$          stop 'xxxxx not qp=qsmesh xxxx'
c$$$c$$$ 1013    continue
c$$$c$$$          do j1=1,ndimh !ldim+napw
c$$$c$$$          do j2=1,ndimh !ldim+napw
c$$$c$$$c            if(j1<=ldim.and.j2<=ldim) cycle
c$$$c$$$c            senex(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$c$$$            sene(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$c$$$          enddo
c$$$c$$$          enddo
c$$$c$$$c          print *,'ldim:abs and diff sene-senex=',sum(abs(sene(1:ldim,1:ldim)))
c$$$c$$$c     &      ,sum(abs(sene(1:ldim,1:ldim)-senex(1:ldim,1g:ldim)))
c$$$c$$$c          print *,'ndimh:abs and diff sene-senex=',sum(abs(sene))
c$$$c$$$c     &      ,sum(abs(sene-senex))
c$$$c$$$c          print *, 'end of test mode of sene newsigmasw=',newsigmasw()
c$$$c$$$c          deallocate(senex)
c$$$c$$$c          stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
c$$$c
c$$$c$$$C--- This works OK for lmf mode when noshrbz=.true. for PMT
c$$$c$$$C    Check comment out "sene=0" and "if(j1<=ldim.and.j2<=ldim) cycle" or not ---it depends on what you want to test.
c$$$c$$$c          allocate(senex(ndimh,ndimh))
c$$$c$$$          allocate(sene(ndimh,ndimh))
c$$$c$$$c          senex=0d0
c$$$c$$$          write(6,"(a,3f13.5)")'    qp=',qp
c$$$c$$$          do ik1=1,nqps
c$$$c$$$             write(6,"(a,i3,3f13.5)")' ik1 qsfz=',ik1,qsfz(1:3,ik1)
c$$$c$$$             iq1=ik1
c$$$c$$$            ii1=mod(iq1,2)+1
c$$$c$$$            ii2=mod(iq1/2,2)+1
c$$$c$$$            ii4=mod(iq1/4,2)+1
c$$$c$$$            print *,' iq1 ii4,ii2,ii1 =',iq1,ii4,ii2,ii1
c$$$c$$$             print *,' xxxxxxxxxxx this shows shorbz problem... xxxxxxxxx'
c$$$c$$$             if(sum(abs(qp-qsfz(1:3,ik1)))<1d-6) goto 1013
c$$$c$$$          enddo
c$$$c$$$          stop 'xxxxx not qp=qsmesh xxxx'
c$$$c$$$ 1013    continue
c$$$c$$$          do j1=1,ndimh !ldim+napw
c$$$c$$$          do j2=1,ndimh !ldim+napw
c$$$c$$$c            if(j1<=ldim.and.j2<=ldim) cycle
c$$$c$$$c            senex(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$c$$$            sene(j1,j2)= sfz(ii4,ii2,ii1,j1,j2,isp)
c$$$c$$$          enddo
c$$$c$$$          enddo
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
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
!     !
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
!!== Obtain q parameter for given qin(3) ==
!! This find eight corners of a cube (each corner is one of the qsmesh). qin is within the cube.
      subroutine qparam(qin,qsmesh,nk1,nk2,nk3,plat,  ikc,qpara,igshift)
      implicit none
      intent(in)::      qin,qsmesh,nk1,nk2,nk3,plat  
      intent(out)::                                   ikc,qpara,igshift
!!i qsmesh(3,ik1,ik2,ik3); we assume it orders well;this routine works well even for regular mesh or off-gamma mesh.
!!r Not sophsticated well... At first, we need to clarify the definition and generality of qsmesh itself.
      integer::nk1,nk2,nk3, nkk(3),ik1,ik2,ik3,i,j,k,ikcc(3),ix,
     &nx(0:1),ny(0:1),nz(0:1), mx(0:1), my(0:1), mz(0:1),ix1,ix2,ix3,
     &igshift(3,0:1,0:1,0:1), ikc(3,0:1,0:1,0:1)
      real(8):: qsmesh(3,nk1,nk2,nk3), qpara(3),qqq(3),qin(3),plat(3,3)
      real(8):: eps=1d-5, qmesh(3), qmeshm(3,0:1,0:1,0:1),qmin(3),qmax(3),qx(3),qqqm(3)
!! Get projected q.-> ik1, ik2 ik3 for given q.
      qqq = matmul(qin,plat)
      nkk= (/nk1,nk2,nk3/)
      ikcc = qqq*nkk + 100*nkk !(1) ---> Takao suppose adding 100*nkk is enough to make ikcc positive.
      ikcc = ikcc    - 100*nkk !(2) Note that this conversion to integer, "(1) and (2)",
      !    do not equilvalent with ikc=qqq*nkk when qqq*nkk is negative.
      !
      qpara = nkk*qqq - ikcc
      ikcc= mod(ikcc+100*nkk,nkk)
      qqqm= (qpara + ikcc)/dble(nkk) !positive
      print *,'nkkk=',nkk
      write(6,"(a,3d13.5,3d13.5,3i5)")'  qin   =',qin
      write(6,"(a,3d13.5,3d13.5,3i5)")'  qqqm   =',qqqm
      write(6,"(a,3d13.5,3d13.5,3i5)")'  qpara =',qpara
      write(6,"(a,3i5)")              '  ikcc  =',ikcc
c qmesh is based on on qlat coordinate. Is this too time-consuming (maybe negligible)?
c         But maybe too general or
c      allocate( qmesh(3,nk1,nk2,nk3) )
      do ik1 = nk1,1,-1
        qmesh = matmul(qsmesh(:,ik1,1,1),plat)
        !print *,' aaa=',qmesh(1),qqqm(1)
        if( qmesh(1)-eps < qqqm(1) ) then
          nx(0)= ik1
          nx(1)= mod(ik1+1-1,nk1)+1
          mx(0) = 0
          mx(1) = (ik1+1- nx(1))/nk1
          go to 1011
        endif
      enddo
      stop ' no exit 1xx'
 1011 continue
      do ik2 = nk2,1,-1
        qmesh = matmul(qsmesh(:,1,ik2,1),plat)
        if( qmesh(2)- eps < qqqm(2) ) then
          ny(0)= ik2
          ny(1)= mod(ik2+1-1,nk2)+1
          my(0)= 0
          my(1)= (ik2+1- ny(1))/nk2
          goto 1012
        endif
      enddo
      stop 'no exit 2xx'
 1012 continue
      do ik3 = nk3,1,-1
        qmesh = matmul(qsmesh(:,1,1,ik3),plat)
        if( qmesh(3)- eps < qqqm(3)  ) then
          nz(0)= ik3
          nz(1)= mod(ik3+1-1,nk3)+1
          mz(0)= 0
          mz(1)= (ik3+1- nz(1))/nk3
          goto 1013
        endif
      enddo
      stop 'no exit 3xx'
 1013 continue

      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            ikc(1:3,ix1,ix2,ix3)     = (/nx(ix1),ny(ix2),nz(ix3)/)
            igshift(1:3,ix1,ix2,ix3) = (/mx(ix1),my(ix2),mz(ix3)/)
          enddo
        enddo
      enddo
ccccccccccccccccccccccccccccccccc
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3i3,3x,4i3)") 'ix1 ix2 ix3 ikc=',ix1,ix2,ix3,ikc(1:3,ix1,ix2,ix3)
          enddo
        enddo
      enddo
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3i3,3x,4i3)") '        igshift=',ix1,ix2,ix3,igshift(1:3,ix1,ix2,ix3)
          enddo
        enddo
      enddo
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3i3,3x,4i3)") '    ikc+nkk*igshift=',ix1,ix2,ix3,ikc(1:3,ix1,ix2,ix3) + nkk*igshift(1:3,ix1,ix2,ix3)
          enddo
        enddo
      enddo
c check
      qmin  = matmul( qsmesh(:,nx(0),ny(0),nz(0)), plat )
      qx  = matmul( qsmesh(:,nx(1),ny(0),nz(0)), plat ) + igshift(1:3,1,0,0) !/dble(nk1)
      qmax(1)=qx(1)
      qx = matmul( qsmesh(:,nx(0),ny(1),nz(0)), plat )  + igshift(1:3,0,1,0) !/dble(nk2)
      qmax(2)=qx(2)
      qx = matmul( qsmesh(:,nx(0),ny(0),nz(1)), plat )  + igshift(1:3,0,0,1) !/dble(nk3)
      qmax(3)=qx(3)
!! region check      
      do ix=1,3
        if( qmin(ix)-eps <= qqqm(ix) .and. qqqm(ix) <= qmax(ix)+eps ) then
        else
           write(6,*)' error regioncheck: ix  =',ix
           write(6,*)' error regioncheck: qqqm=',qqqm
           write(6,*)' error regioncheck: qmin=',qmin
           write(6,*)' error regioncheck: qmax=',qmax
           call rx('error in region check!')
        endif
      enddo
      return
      end
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine interparam(q,wgt)
      implicit none
      intent(in)::          q
      intent(out)::           wgt
C- Calculate interpolarion weights when values at corners of a cube is given.
C  The cube is devided into 6 piramids, whose center are the center of a cube.
Ci q(3) is within a cube, that is, q(i) \in (0,1) for i=1,3.
C     At first, we judge which piramids include the q vector.
C  Then we further devide the piramids into four tetrahedron.
C  We judge which tetrehedron includes the q vector.
C  The interpolation weight is calculated from the linear-interpolarion formular from its four corners.
C- cube
      integer:: ip,ix1,ix2,ix3
      real(8):: q(3),wgt(0:1,0:1,0:1),xx(3),x2d(2),wc
      real(8):: wgtx(0:1,0:1),maxx
      integer:: iz(1),ix
      xx = 2d0*q-(/1d0,1d0,1d0/)
c weight for a center.
      ix = maxloc(abs(xx),dim=1)
      !ix=iz(1)
      maxx = abs(xx(ix))
      if(xx(ix)>=0) ip=1
      if(xx(ix)<0)  ip=0
      print *,' xx=',xx
      print *,' ix=',ix
      print *,' maxx=',maxx
c
      wgt=1d0/8d0*(1-maxx) !center weight
      if(ix==1) then
        x2d=(/xx(2),xx(3)/)
        call prodx2d(x2d,wgtx,wc)
        wgt(ip,:,:) = wgt(ip,:,:) + wgtx + 1/4d0 * (maxx + wc)
        ! face-center weight --> 1/4d0*(maxx-1d0) + 1/4d0*wc
      elseif(ix==2)  then
        x2d=(/xx(1),xx(3)/)
        call prodx2d(x2d,wgtx,wc)
        wgt(:,ip,:) = wgt(:,ip,:) + wgtx + 1/4d0 * (maxx + wc)
      else
        x2d=(/xx(1),xx(2)/)
        call prodx2d(x2d,wgtx,wc)
        wgt(:,:,ip) = wgt(:,:,ip) + wgtx + 1/4d0 * (maxx + wc)
      endif
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3f8.3,a,3i3,f8.3)") '  q=', q,'  wgt=',ix1,ix2,ix3,wgt(ix1,ix2,ix3)
          enddo
        enddo
      enddo
      end
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine prodx2d(x,wgt,wc)
      implicit none
      intent(in)::       x
      intent(out)::        wgt,wc
      real(8):: x(2),wgt(0:1,0:1),wc,wpp,wpm
      real(8):: vpp(2),vpm(2)
      vpp=(/.5d0,.5d0/)
      vpm=(/.5d0,-.5d0/)
      wpp = sum(vpp*x)
      wpm = sum(vpm*x)
      print *,' wpp wpm=',wpp,wpm
      wgt=0d0
      if(wpp>=0) then
        wgt(1,1)=  wpp
        wc= -wpp
      else
        wgt(0,0)= -wpp
        wc=  wpp
      endif
      if(wpm>=0) then
        wgt(1,0)=  wpm
        wc= wc-wpm
      else
        wgt(0,1)= -wpm
        wc= wc+wpm
      endif
      end
      end module m_seneinput_futureversion_unusedyet
