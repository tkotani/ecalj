      module m_rdsigm2
      public:: M_rdsigm2_init, Getsenex, Dsene, senex !getsenex returns the self-energy term.
      public:: ndimsig, sene
      integer,protected,public :: nk1,nk2,nk3
      complex(8),pointer,public:: hrr(:,:,:,:,:,:)
      private
      complex(8),allocatable,target :: sfz(:,:,:,:,:,:)
      complex(8),allocatable,protected:: sene(:,:),senex(:,:)
      complex(8), allocatable,protected :: hhrs (:,:,:,:)
      integer, allocatable,protected :: iaxs(:,:)
      integer, allocatable,protected :: ntabs(:)
      integer,protected:: ndimsig,nspsigm,ndhrs,ham_nqsig
      real(8) , allocatable,private ::  rv_p_oqsig (:)  !not yet protected (rdsigm2 rewrite this)
      contains
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine getsenex(qp,isp,ndimh,ovlm)
!! Return self-energy senex at qp,isp    
      implicit none
      integer:: isp,ndimh,ispsigm
      real(8):: qp(3)
      complex(8),allocatable:: ovlmtoi(:,:),ovliovl(:,:) 
      complex(8):: ovlm(ndimh,ndimh)
      allocate(sene(ndimsig,ndimsig))
      if(isp>nspsigm) then
         ispsigm = nspsigm
      else
         ispsigm=isp
      endif
      call bloch2(qp,ispsigm,sene)
      allocate( ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
      ovlmtoi = ovlm(1:ndimsig,1:ndimsig)
      call matcinv(ndimsig,ovlmtoi)
      ovliovl = matmul(ovlmtoi,ovlm(1:ndimsig,1:ndimh))
      deallocate(ovlmtoi)
      allocate(senex(ndimh,ndimh))
      senex = matmul(transpose(dconjg(ovliovl)), matmul(sene,ovliovl))
      deallocate(ovliovl)
      end subroutine
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine dsene()
      deallocate(senex,sene)
      end subroutine
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine m_rdsigm2_init()
      use m_lmfinit,only : rsrnge=>ham_rsrnge,nbas,pwmode=>ham_pwmode,ham_ldham
      use m_hamindex,only: symops_af,napwmx
      use m_MPItk,only: procid,master
      use m_ext,only:sname
      integer:: ldim,ierr,ifi,ndimh_dummy,nqps,ifile_handle,ifis2,ik1,ik2,ik3,is,iset
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
!! Get self-energy sfz in the full BZ.
         call rdsigm2(nbas,nspsigm,ifi,
     &        nk1,nk2,nk3,ldim,qsmesh2,  nqps, mtosigmaonly(),ndimsig) 
         close(ifi) 
!! Write sig_fbz 
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
!FFT hrr is on regular mesh points. For bloch2
         call fftz3(sfz,nk1,nk2,nk3,nk1,nk2,nk3,ndimsig**2*nspsigm,iset,+1) !+1 backward
         hrr => sfz ! rename sfz as hrr
      endif                  !procid==master
      if(cmdopt0('--wsig_fbz')) call rx0('end of --wsig_fbz mode')
      mlog     = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
      call mpibc1_int(nspsigm,1,'senebr:nspsigm')
      call mpibc1_int(nk1,1,'senebr:nk1')
      call mpibc1_int(nk2,1,'senebr:nk2')
      call mpibc1_int(nk3,1,'senebr:nk3')
      if(procid/=master) allocate(sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm))
      call mpibc1_complex(sfz, size(sfz),'senebroadcase:sfz')
      if(procid/=master) hrr => sfz    ! rename sfz as hrr
      if(procid==master) write(6,*)' --- end of sigmainit: rdsigm2 initialization section ---'
      return
 9995 continue
      call rx('sigmainit: readin error of sigm file')
      end subroutine m_rdsigm2_init
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine rdsigm2(nbas_dummy,nspsigm,ifis, 
     &     nk1,nk2,nk3,ldim,qsmesh    ,nqps,mtosigmaonly,ndimsig) !,qsfz,qqps are for test sfz
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
      implicit none
      integer mode,ifis,ndimsig_r,lwsig, nbas_dummy!nbas,
      double precision rsrnge
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nspsigm,nglob,ldham(16),lrsig,
     .nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
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
c      complex(8)::sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm)
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
      tolq = 1d-6
      ledit = 0
      nscnts = 0
      ltrans=.true.
      lshft=0
      rewind ifis
      read(ifis) nsp_,ndimsig_r,nk1_,nk2_,nk3_,nqp_
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
      end subroutine rdsigm2
!!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine hamfb3k(qin,iq,nk1,nk2,nk3,
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
     .     (jj2*ifac(2)-1)*qb(k,2) +
     .     (jj3*ifac(3)-1)*qb(k,3)
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
            if(debugmode()>0) write(6,"(a,3f13.5)")' input          qin = ', qin
            if(debugmode()>0) write(6,"(a,3f13.5)")' target a        q1 = ', q1 ! qin = g(:,:,ig)^{-1}*q1
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
      end module m_rdsigm2
