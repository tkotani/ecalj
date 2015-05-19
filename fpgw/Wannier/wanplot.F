       program wanplot
!! == Wannier function plot. Wannier function is expanded in the PW (spacial mesh). ==
!! NOTE: Because os lazyness, not yet MPI. In cases, it may be useful...
!! === Usage: ===
!!  We first need to run maxloc generation as
!!   echo 0|hbasfp0
!!   echo 1|$nfpgw/hmaxloc   >lmaxloc1
!!   $nfpgw/hpsig            >lpsig
!!   echo 2|$nfpgw/huumat    >luumat2
!!   echo 2|$nfpgw/hmaxloc   >lmaxloc2
!!
!! === Remarks ===
!!   iclass means equivalent sites.
!!   Use Xcrysden to see xsf file. Need to figure out automatic controll.
!!   Currently a little inconvenient.
!!   plot Wannier functions. Do this after Wannier matrix dnk is generated.
!!
!! ==== History ====
!! takao modified jul2014 from calc_wannier6 by H.Kino.
!! 080603, For plotting arbitary region, from calc_wannier
!! 071114, Takashi Miyake, from calc_wfn.F
!! -------------------------------------------------------------------------
c      program calc_wannier6
c      use m_LMTO
c      use m_MLWF
c      use m_wfrho_abc 
      use keyvalue,only: getkeyvalue
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval,readcphi,readgeig
      use m_read_bzdata,only: read_bzdata, nqbz,nqibz,nqbzw,nteti,ntetf
     &     ,n1,n2,n3,qbas,ginv,qbasmc,qbz,wbz,qibz,wibz,qbzw,idtetf,ib1bz,idteti
     &     ,nstar,irk,nstbz,ngrp2=>ngrp,qibz_r,nqibz_r
      use m_DATA4GW,only: read_data4gw,set_mnla,iclass,nclass,zz,alat,nbas,nsp,plat,ldim2,bas
      use m_QG,only: read_qg,ngp
      use m_cubeformat 
      use m_xsfformat 
      use m_expand_mesh 
      implicit none
      integer :: nwf,iko_ix,iko_fx,nband_wfn
      real(8) :: a,qdum(3)
      integer :: nq_wfn,tlat(3),nfac
      integer,allocatable :: bindx_wfn(:),bindx_wfn0(:)
      real(8),allocatable:: qreg(:,:),q_wfn(:,:)
      complex(8),allocatable :: 
     & phipw0(:,:,:,:,:,:), phiaug0(:,:,:,:,:,:),phitot0(:,:,:,:,:,:),
     & phipw(:,:,:,:,:,:),  phiaug(:,:,:,:,:,:), phitot(:,:,:,:,:,:),
     & phipw1(:,:,:,:,:,:), phiaug1(:,:,:,:,:,:),phitot1(:,:,:,:,:,:),
     & wanpw(:,:,:,:,:),  wanaug(:,:,:,:,:),wantot(:,:,:,:,:),
     & wanpw1(:,:,:,:,:), wanaug1(:,:,:,:,:),wantot1(:,:,:,:,:),
     &  dnk(:,:,:,:)
      integer:: i_rini(3),i_rfin(3)
      integer :: npw(3),mesh(3),mesh0(3),mesh1(3),meshrange(2,3)
      integer :: iq,ib,  ngpmx
      integer:: i,id,j,ifi,iqbz2
      integer :: isp,iqbz,ikp,iwf,nprecb,mrecb,mrece,nlmtot,nqbzt, mrecg,ifhbed,iopen,nband
      real(8):: r_rini0(3),r_rfin0(3),ang,alat_ang
      real(8):: r_rini(3), r_rfin(3),r
      real(8):: r_rini1(3), r_rfin1(3)
      logical :: lrho,flag
      character(10):: vis_unit='none'
      character(20):: outputformat ='none'
      character(20)::inputfile='GWinput'
      real(8)::quu(3),det,qlat(3,3),q(3)
      complex(8),allocatable :: cphi(:,:,:,:),geig(:,:,:,:),geig2(:,:),cphi2(:,:)
      character*4::fname
      logical:: debug=.false.,vis_skip
!-----------------------------------------------------
      call getkeyvalue(inputfile,'vis_skip',  vis_skip,  default=.false. ) 
      if(vis_skip) call rx0('wanplot: we found vis_skip on. Do nothing and Quit!')
!! Readin all data
      call read_qg()
      call read_BZDATA()
      call read_DATA4GW()
      call set_mnla()
      call dinv33x (plat,qlat)  !it was dinv33(plat,1,qlat) by Ferdi
      do isp=1,nsp
         ifi = 1000
         if (isp.eq.1) fname='MLWU'
         if (isp.eq.2) fname='MLWD'
         open(ifi,file=fname,form='unformatted',status='old', action='read')
         read(ifi)nqbz,nwf,iko_ix,iko_fx
         if (isp.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp),qreg(3,nqbz))
         do iqbz = 1,nqbz
            read(ifi)iqbz2,q(1:3)
            if(debug)write(6,"(i5,3f13.5)") iqbz,q(:)
            qreg(:,iqbz) = q
            read(ifi)dnk(iko_ix:iko_fx,1:nwf,iqbz,isp)
         enddo
         close(ifi)
      enddo                     ! isp
      write(6,*)'read end of MLWU/D ...'
!!
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      call init_readeigen(ginv,nsp,nband,mrece) !initialization of readEigen
      call init_readeigen2(mrecb,ldim2,mrecg) !initialize m_readeigen

!! replace geig and cphi with those for WF. Converted by the dnk matrix.
!! geig2 = cphi2 = 0 for ik > nqbz
      call readngmx('QGpsi',ngpmx)
      write(6,*)'ngpmx nband ldim2=',ngpmx,nband,ldim2
      write(6,*)'nwf   nqbz  nsp  =',nwf,nqbz,nsp
      allocate(geig2(ngpmx,nband))
      allocate(cphi2(ldim2,nband))
      allocate(geig(ngpmx,nwf,nqbz,nsp))
      allocate(cphi(ldim2,nwf,nqbz,nsp))
      geig = 0d0
      cphi = 0d0
      do ikp = 1,nqbz
      do isp = 1,nsp
c         write(6,"('ikp isp=',2i5,3f12.5)") ikp,isp,qreg(:,ikp)
         call readcphi(qreg(:,ikp),ldim2,isp, quu, cphi2)
         if(sum(abs(qreg(:,ikp)-quu))>1d-6) call rx('wanplot: mmlf111eeeee')
         call readgeig(qreg(:,ikp),ngpmx,isp, quu, geig2)
         if(sum(abs(qreg(:,ikp)-quu))>1d-6) call rx('wanplot: mmlf222eeeee')
!! may2015 use zaxpy. This can avoid bug? when wanplot by qsub.
         do iwf = 1,nwf
         do ib  = iko_ix,iko_fx
            call zaxpy(ngp(ikp),dnk(ib,iwf,ikp,isp),geig2(1,ib),1,geig(1,iwf,ikp,isp),1)
            call zaxpy(ldim2,   dnk(ib,iwf,ikp,isp),cphi2(1,ib),1,cphi(1,iwf,ikp,isp),1)
c           geig(:,iwf,ikp,isp) = geig(:,iwf,ikp,isp) + geig2(:,ib)*dnk(ib,iwf,ikp,isp)
c           cphi(:,iwf,ikp,isp) = cphi(:,iwf,ikp,isp) + cphi2(:,ib)*dnk(ib,iwf,ikp,isp)
         enddo ! ib
         enddo ! iwf
      enddo ! isp
      enddo ! ikp
      deallocate(geig2,cphi2,dnk)
      write(6,*) '### ib,bas(1:3,ib) ############'
      do ib=1,nbas
         write(*,"(i5,3f12.6)")ib,bas(1:3,ib)
      enddo

!! NOTE: nq_wfn = nqbz
      nq_wfn = nqbz  
      allocate(q_wfn(3,nq_wfn))
      q_wfn(1:3,1:nqbz) = qbz(1:3,1:nqbz)

!! == Readin vis_* settings in GWinput
      call getkeyvalue(inputfile,'vis_wan_band_n',nband_wfn,default=nwf)
      if (nband_wfn.gt.nband) call rx('wanplot: nband_wfn.gt.nband !')
      write(6,"(a,2i5)") '### nq_wfn, nband_wfn =',nq_wfn,nband_wfn
      allocate(bindx_wfn0(nband_wfn),bindx_wfn(nband_wfn)) 
      do ib=1,nband_wfn
        bindx_wfn0(ib) = ib
      enddo
      call getkeyvalue(inputfile,'vis_wan_band_id',bindx_wfn,default=bindx_wfn0,size=nband_wfn) 
      do ib=1,nband_wfn
        write(*,"(a,2i5)") 'ib bndinx=',ib,bindx_wfn(ib) 
      enddo
!! 
      call getkeyvalue(inputfile, 'vis_wan_tvec',   tlat,size=3,default=(/0,0,0/))
      write(*,"(a,3i5)")'### tlat',tlat
      call getkeyvalue(inputfile, 'vis_wan_interpolation', nfac,default=1)  !FFT
      vis_unit='abc'
      write(6,*)' CAUTION: range of ubond and lbound are in abc(cell) unit'
c      call getkeyvalue(inputfile, 'vis_wan_unit', vis_unit)
c      if ( trim(vis_unit).ne.'abc') then
c         write(*,*) 'support only vis.wan.unit=abc'
c         stop 'support only vis.wan.unit=abc'
c      endif
      call calc_npw(nfac,npw)
      call getkeyvalue(inputfile,'vis_wan_mesh',  mesh0,  size=3, default=(/10,10,10/) ) 
                                !mesh size 0:mesh0(1),...
      call getkeyvalue(inputfile,'vis_wan_lbound',r_rini0,size=3,default=(/-1d0,-1d0,-1d0/))!lower bound
      call getkeyvalue(inputfile,'vis_wan_ubound',r_rfin0,size=3,default=(/1d0,1d0,1d0/))   !upper bound
      write(*,*) ' mesh=',mesh0
      write(*,*) ' lbound=',r_rini0
      write(*,*) ' ubound=',r_rfin0
c      call getkeyvalue(inputfile,'vis_wan_outputformat',outputformat,default='xsf') 

      do i=1,3
        i_rini(i)= floor(r_rini0(i))
        i_rfin(i)= ceiling(r_rfin0(i))
      enddo
!! for plot mesh
      r_rini= i_rini
      r_rfin = i_rfin
      do i=1,3
         mesh(i)= (i_rfin(i)-i_rini(i))*mesh0(i)
      enddo 
      write(6,*)'mmm: i_rfin=',i_rini
      write(6,*)'mmm: i_rini=',i_rfin
      write(6,*)'mmm:   mesh=',mesh

      allocate(phipw0 (mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,   nband_wfn,nq_wfn,nsp))
      allocate(phiaug0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,  nband_wfn,nq_wfn,nsp))
      allocate(phitot0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,  nband_wfn,nq_wfn,nsp))

!! == Generate phi, which is the real-space rep. of the Bloch functions (on mesh0, real mesh points). 
!! Time consuming part.
      call calc_phiall_abc2(nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     i     npw,mesh0,nsp,nband,ldim2,ngpmx,  
     i     geig,cphi,nwf,
     o     phipw0,phiaug0,phitot0)
      write(6,*)'sumchk 000 =',sum(abs(phipw0)),sum(abs(phiaug0)),sum(abs(phitot0))

      allocate(phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,   nband_wfn,nq_wfn,nsp))
      allocate(phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,  nband_wfn,nq_wfn,nsp))
      allocate(phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,  nband_wfn,nq_wfn,nsp))

!! phipw,phiaug,phitot are the real-space rep. on extended mesh points (mesh).
      call expand_mesh(plat,
     i    nq_wfn,nband_wfn,q_wfn,nsp,
     i    i_rini,i_rfin, 
     o    mesh0, phipw0,phiaug0,phitot0, 
     o    mesh, phipw,phiaug,phitot )
c      write(6,*)'sumchk 22222=',sum(abs(phipw)),sum(abs(phiaug)),sum(abs(phitot))

!! from Bloch to Wannier
      allocate(wanpw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp) )
      allocate(wanaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))
      allocate(wantot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))

      call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phipw,       wanpw )
      call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phiaug,      wanaug)
      call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phitot,      wantot)
      write(*,*) 'mesh in b2w'
      write(*,*) 'mesh=',mesh
      write(*,*) 'r_rini=',r_rini
      write(*,*) 'r_rfin=',r_rfin
      write(*,*)' '

!!..... calculate range  from rini0 and rfin0 !abc (fractional basis only)
!!    r = rini + (rfin-rini)*(i-1)/mesh
      meshrange=0
      r_rini1=0
      r_rfin1=0
      do id=1,3
         do i=1,mesh(id)+1
            r =  r_rini(id) + (r_rfin(id)-r_rini(id))*(i-1)/mesh(id)
            if ( r> r_rini0(id) ) then 
               j=i-1
               meshrange(1,id)=j
               r_rini1(id) = r_rini(id) + (r_rfin(id)-r_rini(id))*(j-1)/mesh(id)
               exit
            endif
         enddo
         do i=1,mesh(id)+1
            r =  r_rini(id) + (r_rfin(id)-r_rini(id))*(i-1)/mesh(id)
            if ( r> r_rfin0(id) .or.i==mesh(id)+1) then 
               j=i
               meshrange(2,id)=j 
               r_rfin1(id) = r_rini(id) + (r_rfin(id)-r_rini(id))*(j-1)/mesh(id)
               exit
            endif
         enddo
      enddo
      print *,'meshrange2=',meshrange(2,:) !upper limits
      print *,'meshrange1=',meshrange(1,:) !lower limits

      mesh1(:)=meshrange(2,:)-meshrange(1,:)
      allocate(wanpw1 (mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
      allocate(wanaug1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
      allocate(wantot1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
      allocate(phipw1 (mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))
      allocate(phiaug1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))
      allocate(phitot1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))

      write(*,*)'range in inputfile'
      write(*,*) 'rini0=',r_rini0
      write(*,*) 'rfin0=',r_rfin0
      write(*,*)' '
      write(*,*)'cutted mesh'
      write(*,*) 'mesh=',mesh1
      write(*,*) 'rini=',r_rini1
      write(*,*) 'rfin=',r_rfin1
      write(*,*)' '

      wanpw1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) =
     .     wanpw(meshrange(1,1):meshrange(1,1)+mesh1(1), 
     .     meshrange(1,2):meshrange(1,2)+mesh1(2),
     .     meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )
      wanaug1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) =
     .     wanaug(meshrange(1,1):meshrange(1,1)+mesh1(1),
     .     meshrange(1,2):meshrange(1,2)+mesh1(2),
     .     meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )
      wantot1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) =
     .     wantot(meshrange(1,1):meshrange(1,1)+mesh1(1),
     .     meshrange(1,2):meshrange(1,2)+mesh1(2),
     .     meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )

      phipw1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) =
     .     phipw(meshrange(1,1):meshrange(1,1)+mesh1(1),
     .     meshrange(1,2):meshrange(1,2)+mesh1(2),
     .     meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
      phiaug1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) =
     .     phiaug(meshrange(1,1):meshrange(1,1)+mesh1(1),
     .     meshrange(1,2):meshrange(1,2)+mesh1(2),
     .     meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
      phitot1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) =
     .     phitot(meshrange(1,1):meshrange(1,1)+mesh1(1),
     .     meshrange(1,2):meshrange(1,2)+mesh1(2),
     .     meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
c rini -> r_rini1
c rfin -> r_rfin1
c mesh -> mesh1

!! Dump phi and wannier functions -------------------
      write(*,*) '-- dump phi(bloch funciton) and wan(Wannier funciton) --'
      qdum = 0.0d0
      ang = 0.529177d0
      alat_ang=alat*ang
!! opendx and cube need to be fixed...
c$$$      if (outputformat.eq.'opendx') then
c$$$         call wfn2dx_2(alat_ang,plat,nsp,1,nband_wfn,qdum,bindx_wfn,
c$$$     &        mesh1,r_rini1,r_rfin1,wanpw1,wanaug1,wantot1)
c$$$         call crystal2dx_2(alat_ang,plat,r_rini1,r_rfin1,
c$$$     &        nbas,bas,nclass,iclass,zz)
c$$$      else if (outputformat.eq.'cube') then
c$$$         call wrt_cube(
c$$$     i        'wan',
c$$$     i        alat,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
c$$$c     i     mesh,rini,rfin,phipw,phiaug,phitot  ! for bloch orbital
c$$$     i        mesh1,r_rini1,r_rfin1,wanpw1,wanaug1,wantot1, ! for wannier function 
c$$$     i        nbas,bas,nclass,iclass,zz )
c$$$         call wrt_cube(
c$$$     i        'phi',
c$$$     i        alat,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
c$$$     i        mesh1,r_rini1,r_rfin1,phipw1,phiaug1,phitot1, ! for bloch orbital
c$$$c     i     mesh,rini,rfin,wanpw,wanaug,wantot,  ! for wannier function 
c$$$     i        nbas,bas,nclass,iclass,zz )
c$$$      else                      !--- if(outputformat.eq.'xsf') then, default
         write(6,*) 'Writing xsf (Xcrysden) file...'
         call wrt_xsf(
     i        'wan',vis_unit,
     i        alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
     i        mesh1,r_rini1,r_rfin1,wanpw1,wanaug1,wantot1, ! for wannier function 
     i        nbas,bas,nclass,iclass,zz )
         call wrt_xsf(
     i        'phi',vis_unit,
     i        alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
     i        mesh1,r_rini1,r_rfin1,phipw1,phiaug1,phitot1, ! for bloch orbital
     i        nbas,bas,nclass,iclass,zz )
c$$$      endif 

! dump rho or not ----------
      lrho=.false. 
      write(6,"(a,l)") 'dump rho? [T/F] (if needed set lrho=T in wanplot.F)=',lrho
      if (lrho) call calc_rho_2(alat_ang,nq_wfn,nband_wfn,mesh,r_rini,r_rfin,wanpw,wanaug,wantot)
!      if (lrho) call calc_rho(nq_wfn,nband_wfn,npw,phipw,phiaug,phitot)
!      if (lrho) call calc_rho(1,nband_wfn,npw,wanpw,wanaug,wantot)
c      call cputid(0)
      call rx0('wanplot: ok')
      end program wanplot

!-----------------------------------------------------
      subroutine calc_npw(nfac, npw)
      use m_QG,only: ngvecp,qqqa,nqnum,ngp
      use m_DATA4GW,only: alat,plat
      implicit none
c input
      integer :: nfac
c output
      integer :: npw(3)
c local
      integer :: iq,ig,id,itmp(3),ntmp(3)
      double precision :: pi,gtmp(3),gcutmax,gcuttmp,at(3,3),g(3,3)
      logical :: debug=.false.
      write(*,"(a)") '--- calc_npw ---'
      call mytranspose(plat,At,3,3)
      call myinv3(At,G)
      pi=4.0d0*atan(1.0d0)
      ntmp(1:3)=0
      do iq=1,nqnum
        gcutmax=-1.0d0
        do ig=1,ngp(iq)
          call mymatvec(G,dble(ngvecp(1:3,ig,iq)),gtmp,3,3)
          gtmp(1:3)=gtmp(1:3)+qqqa(1:3,iq)
          gtmp(1:3)=gtmp(1:3)*2.0d0*pi/alat
          gcuttmp=sqrt(sum(gtmp(1:3)**2))
          if (gcutmax.lt.gcuttmp) gcutmax=gcuttmp
          do id=1,3
            itmp(id)=abs(ngvecp(id,ig,iq))
            if (ntmp(id).lt.itmp(id)) ntmp(id)=itmp(id)
          enddo
        enddo
        if(debug) write(*,"(a,2i5,f10.5)") '# iq ngp gcutmax= ',iq,ngp(iq),gcutmax
      enddo
!      npw(1:3)=2*ntmp(1:3)+2
      npw(1:3)=nfac*ntmp(1:3)+2
      write(*,"(a,3i6)") '# npw(1:3)=',npw(1:3)
      end subroutine calc_npw

!-----------------------------------------------------
c Linear interpolation of gx/r
      double precision function calc_gxr(r,l,n,ic,isp)
c      use m_LMTO      
      use m_DATA4GW,only: bb,nr,aa,alat,gx
      implicit none
c input
      double precision :: r
      integer :: l,n,ic,isp
c local
      double precision :: r1,r2
      integer :: ir
      ir=1+int(log(r/bb(ic)+1.0d0)/aa(ic))
      if (ir.lt.1) stop 'ir.lt.1'
      if (ir.gt.nr(ic)-1) stop 'ir.gt.nr(ic)-1'
      r1=bb(ic)*(exp((ir-1)*aa(ic))-1d0)
      r2=bb(ic)*(exp((ir  )*aa(ic))-1d0)
      if (r1.gt.r) stop 'r1.gt.r'
      if (r2.le.r) stop 'r2.le.r'
      calc_gxr=(r-r2)/(r1-r2)*gx(ir,l,n,ic,isp) 
     &     + (r-r1)/(r2-r1)*gx(ir+1,l,n,ic,isp)
      calc_gxr=calc_gxr/(r+1.0d-20)
      end function calc_gxr

!-----------------------------------------------------
      subroutine b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,npw,
     &                 phi,wan)
!! Make Wannier functions from Bloch functions in real space representation.
      use m_DATA4GW,only: bb,nr,aa,alat,nsp,plat
      implicit none
      integer :: nq_wfn,nband_wfn,npw(3),bindx_wfn(nband_wfn),tlat(3)
      double precision :: q_wfn(3,nq_wfn),tvec(3),phase,pi,rtmp(3)
      double complex :: 
     &     phi(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp),
     &     wan(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nsp)
     &    ,ephase
      integer :: iq,isp
c debug:
c      wan(:,:,:,:,1) = phi(:,:,:,:,2,1)
c      return
      pi = 4.0d0*atan(1.d0)
      rtmp(:) = dble(tlat(:))
      call mymatvec(plat,rtmp,tvec,3,3)          
      tvec(1:3)=alat*tvec(1:3)
      wan = (0.0d0,0.0d0)
      do isp = 1,nsp
      do iq = 1,nq_wfn
        phase=2.0d0*pi/alat*sum(q_wfn(1:3,iq)*tvec(1:3))
        ephase=dcmplx(cos(phase),-sin(phase))
        wan(:,:,:,:,isp) = wan(:,:,:,:,isp) + phi(:,:,:,:,iq,isp)*ephase
      enddo ! iq 
      enddo ! isp
      wan = wan / dble(nq_wfn)
      end subroutine b2w

!--------------------------------------------------------
      subroutine chkinv33(a,b)
      implicit none
      integer :: i,j
      double precision :: a(3,3),b(3,3),c(3,3),r,eps
      eps = 1.d-6
!      c = matmul(a,b)
      do i=1,3
      do j=1,3
         c(i,j) = sum(a(:,i)*b(:,j))
      enddo
      enddo
      do i=1,3
         c(i,i) = c(i,i)-1.0d0
      enddo
      do i=1,3
      do j=1,3
         r = abs(c(i,j))
         if (r.gt.eps) stop 'chkinv33 error'
      enddo
      enddo
      end subroutine chkinv33
!--------------------------------------------------------
      subroutine calc_rho_2(alat_ang,nq_wfn,nband_wfn,mesh,rini,rfin,
     &   phipw,phiaug,phitot)
c      use m_LMTO
      use m_DATA4GW,only: bb,nr,aa,alat,nsp,plat
      implicit none
c input
      integer :: nq_wfn,nband_wfn,mesh(3)
      double precision :: alat_ang,rini(3),rfin(3)
      double complex :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

      double complex,allocatable :: rhopw(:,:,:),
     &     rhoaug(:,:,:),rhotot(:,:,:)

c
      integer :: isp,iq,ib,i1,i2,i3
      double precision :: nel

      double precision :: qdum(3),vol
      integer :: bindxdum
      write(*,*) '--- calc_rho ---'
      call mydet3(plat,VOL)
      VOl=abs(VOL)*alat**3

c Allocate rho
      allocate(rhopw(mesh(1)+1,mesh(2)+1,mesh(3)+1))
      allocate(rhoaug(mesh(1)+1,mesh(2)+1,mesh(3)+1))
      allocate(rhotot(mesh(1)+1,mesh(2)+1,mesh(3)+1))

      rhopw(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0
      rhoaug(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0
      rhotot(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0

      do isp=1,nsp
      do iq=1,nq_wfn
      do ib=1,nband_wfn
        do i3=1,mesh(3)+1
        do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
          rhopw(i1,i2,i3)=rhopw(i1,i2,i3)+
     &         abs(phipw(i1,i2,i3,ib,iq,isp))**2
          rhoaug(i1,i2,i3)=rhoaug(i1,i2,i3)+
     &         abs(phiaug(i1,i2,i3,ib,iq,isp))**2
          rhotot(i1,i2,i3)=rhotot(i1,i2,i3)+
     &         abs(phitot(i1,i2,i3,ib,iq,isp))**2

        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      rhopw(:,:,:)=
     &     rhopw(:,:,:)/dble(nq_wfn)
      rhoaug(:,:,:)=
     &     rhoaug(:,:,:)/dble(nq_wfn)
      rhotot(:,:,:)=
     &     rhotot(:,:,:)/dble(nq_wfn)


      nel=0.0d0
      do i3=1,mesh(3)
      do i2=1,mesh(2)
      do i1=1,mesh(1)
        nel=nel+rhotot(i1,i2,i3)
      enddo
      enddo
      enddo
      nel=nel*dble(3-nsp)*VOL/dble(mesh(1)*mesh(2)*mesh(3))
      write(*,*) 'nel = ',nel

      qdum(1:3)=0.0d0
      bindxdum=0
      call wfn2dx_2(alat_ang,plat,1,1,1,qdum,bindxdum,
     &     mesh,rini,rfin,rhopw,rhoaug,rhotot)

      deallocate(rhopw)
      deallocate(rhoaug)
      deallocate(rhotot)
      end subroutine calc_rho_2

