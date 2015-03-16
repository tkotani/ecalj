      module m_zmel
!> Get zmel <MPB psi|psi> 
!! We do "call get_zmelt". Then we have matrix elements zmel (exchange=F for correlation)
!! or zmeltt (exchange=T). 
!! 
!! input

!! basis data for crystal structure.
      use m_genallcf_v3,only: 
     &  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     &  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw,
     &  alat,ef, diw,dw,delta,deltaw,esmr,symgrp,clabl,iclass,nindxv,nindxc,ncwf,
     &  invg, il, in, im, ilnm, nlnm, ilv,inv,imv,  ilnmv, nlnmv,
     &  ilc,inc,imc,  ilnmc, nlnmc,nindx,konf,icore,ncore,
     &  occv,unoccv ,occc,unoccc, nocc, nunocc, plat, pos,z,ecore, freq, symgg,
     &  done_genallcf_v3  
!! Get eigenfuncitons. cphi is coefficients of MTO+lo part, geig is IPW parts.
      use m_readeigen,only: readcphi,readgeig
!! Basic data set to get zmel*
      use m_rdpp,only:       
     &  nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr,
     &  done_rdpp
      use m_read_bzdata,only: 
     &  nqbz,nqibz,  qbas,ginv,qbz,qibz,wbz,
     &  done_read_bzdata
      use keyvalue,only: getkeyvalue
!!
      integer,parameter:: NULL=-99999
      integer,allocatable :: miat(:,:)
      real(8),allocatable :: tiat(:,:,:),shtvg(:,:)
      integer:: nband=NULL,ngcmx=NULL,ngpmx=NULL,ntq=NULL
      integer,allocatable :: itq(:)
      integer,allocatable:: imdim(:)
      integer:: nmmax=NULL
      real(8),allocatable:: ppbir(:,:,:)
!! input 
      complex(8),allocatable,target :: ppovlz(:,:)
!! output
      complex(8),allocatable :: zmel(:,:,:),zmeltt(:,:,:)
!! local
      real(8),private:: qbasinv(3,3),q_bk(3)=1d10,qk_bk(3)=1d0
      logical,private:: init=.true.
c      complex(8),allocatable :: zcousq(:,:)
      complex(8),allocatable,private :: cphiq(:,:), cphim(:,:)
      real(8),allocatable,private :: rmelt(:,:,:),cmelt(:,:,:)
      integer,private::kxold=-9999
!!
      contains
!!
!! <phi(q,nqmax,ispq) |phi(q-rkvec,nctot+nmmax,ispm) MPB(rkvec,ngb)>
!!
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine get_zmelt(exchange,q,kx,  kvec,irot,rkvec,kr,isp,
     &   ngc,ngb,nmmax,nqmax)
!! kvec is in the IBZ, rk = Rot_irot(kvec)
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
      implicit none
      logical:: exchange
      integer:: kx,kr,isp,ngc,ngb,nmmax,nqmax,irot,ispq,ispm,nmini,nqini
      real(8) ::  quu(3),q(3), kvec(3),rkvec(3)
      ispq = isp
      ispm = isp
      nmini=1
      nqini=1
      call get_zmelt2(exchange,  
     &   kvec,irot,rkvec,ngc,ngb,  !MPB for       MPB_rkvec
     &   nmini,nmmax,ispm,         !middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq)      !end state for phi_q
      end
c ntqxx--->nqmax
c nbmax -->nmmax
!! --------------------
!!note: for usual correlation mode, I think nctot=0
!!note: in self-energy mode, we use
!!       nstate=nctot+nmmax
!!       allocate(zmeltt(ngb,    nstate,  nqmax))
!!           < MPB     phi       | phi > (but true matrix elements are for <phi|phi MPB> (complex conjugate).
!!             rkvec   q-rkvec   |  q
 !                      cphim    | cphiq
!                       ispm     | ispq
!                    nmini:nmmax | nqini:nqmax
!
!!--- For dielectric funciton, we use irot=1 kvec=rkvec=q
!!              q      rkvec     | q + rkvec  
!                      nkmin:nt0 | nkqmin:ntp0
!                         occ    | unocc     
!                      (nkmin=1)
!                      (cphi_k  | cphi_kq !in x0kf)
!
!!     rkvec= rk(:,k)-qq  ! <phi(q+rk,nqmax)|phi(rk,nctot+nmmax)  MPB(q,ngb )>
!!     kvec = rk(:,k)-qq  ! k
!!
      subroutine get_zmelt2(exchange,
     &   kvec,irot,rkvec,ngc,ngb,  !MPB for       MPB_rkvec
     &   nmini,nmmax,ispm,         !middle for    phi_{q-rkvec}
     &    q,nqini,nqmax,ispq)      !end state for phi_q
!! \parameter all inputs
!! \parameter output=rmelt,clemt  matrix <MPB psi|psi>
      implicit none
      logical:: exchange
      integer:: invr,nxx,itp,irot,isp,kr,no,nmmax,ngc,ngb,nqmax,nbcut
      integer:: iatomp(natom),nmini,nqini
      real(8) :: symope(3,3),shtv(3),tr(3,natom),qk(3),det
     &  , quu(3),q(3), kvec(3),rkvec(3),wtt
c      complex(8) :: cphiq(nlmto,nband), cphim(nlmto,nband)
      complex(8),allocatable :: zzzmel(:,:,:),zw (:,:)
      integer:: nt
      real(8),allocatable :: drealzzzmel(:,:,:), dimagzzzmel(:,:,:) ,ppb(:)
      logical:: debug=.false.
      complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)

      complex(8):: expikt(natom)
      integer:: it,ia,kx,verbose,nstate
      logical:: oncew
      real(8),parameter::epsd=1d-6

      integer::ispq,ispm,iii,itps
c      complex(8), allocatable :: zwz4(:,:),zwz44(:,:),pomat(:,:), zwzs(:)
c      complex(8), intent(in) :: pomatr(nnmx,nomx,nkpo)
c      allocate(expikt(natom))
c      integer(4) ::nbcut

!TIME0 
      if(verbose()>80) debug=.true.
      if(debug) write(*,*) 'get_zmel in m_zmel: start'
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )

      if(.not.done_genallcf_v3) call rx('m_zmel: not yet call genallcf_v3')
      if(.not.done_rdpp)        call rx('m_zmel: not yet call rdpp')
      if(.not.done_read_bzdata) call rx('m_zmel: not yet call read_bzdata')
      if(init) then
         call dinv33(qbas,0,qbasinv,det)
         allocate( cphiq(nlmto,nband), cphim(nlmto,nband))
         init=.false.
      endif

      if(sum(abs(q-q_bk))>epsd) then
         call readcphi(q, nlmto,ispq, quu, cphim )
         cphiq(1:nlmto,1:ntq) = cphim(1:nlmto,itq(1:ntq))
         q_bk=q
      endif
      if(debug) write(*,*) 'get_zmel in m_zmel: 11111111'
      allocate( rmelt(ngb, nctot+nmmax, nqmax), ! nstate= nctot+nband
     &  cmelt(ngb, nctot+nmmax, nqmax))
      if(debug) write(*,*) 'get_zmel in m_zmel: 22222222'

c      do      it = 1,ntq
c         itp        = itq(it)
c         cphiq(1:nlmto,it) = cphim(1:nlmto,itp)
cc            write(*,*)'svvvv ',it, itp, sum(cphiq(:,it))
c      enddo
c               write(*,*)'srrrrr 1c',sum(cphiq(:,1:ntq)),ntq
c
c      cphiq(1:nlmto,1:ntq) = cphim(1:nlmto,itq(1:ntq))

!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      qk =  q - rkvec        
      if(sum(abs(qk-qk_bk))>epsd) then
         call readcphi(qk, nlmto,ispm, quu, cphim)
         qk_bk=qk
      endif
      if(debug) write(*,*) 'get_zmel in m_zmel: 3333333333'

      allocate( ppb(nlnmx*nlnmx*mdimx*nclass))
c      call getsrdpp2( nclass,nl,nxx)
!! Rotate atomic positions invrot*R = R' + T
      invr  =  invg(irot)       !invrot (irot,invg,ngrp)
      tr    = tiat(:,:,invr)
      iatomp= miat(:,invr)
      symope= symgg(:,:,irot)
      shtv  = matmul(symope,shtvg(:,invr))

c$$$!! ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
c$$$               call ppbafp_v2 (irot,ngrp,isp,nspin,
c$$$     i              il,in,im,nlnm, !w(i_mnl),
c$$$     d              nl,nn_,nclass,nlnmx,
c$$$     i              mdimx,lx,nx,nxx, !Bloch wave    
c$$$     i              cgr, nl-1,  !rotated CG
c$$$     i              ppbrd,      !radial integrals
c$$$     o              ppb)
      ppb=ppbir(:,irot,ispq)
      if(debug) write(6,*) ' sxcf: 1100c',ncore,nqmax
!TIME1 "end of ppbafp_v2"
!TIME0
!! ====================================================
!!     matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)>
!!     including the phase factor exp(ik.T)
!!     B(rot*k,r) = B(k,invrot*r)
!! ====================================================
!! phase factors expikt(ia) is for exp(ik.T(R))

      do ia = 1,natom
         expikt(ia) = exp(img *tpi* sum(kvec*tr(:,ia)) ) 
      end do
!!   core 
      nt   = nctot + nmmax      ! = nstate for the case of correlation
      allocate( zzzmel(nbloch,nt,nqmax))
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      write(*,*)'srrrrrrr 1',icore,ncore
c$$$      write(*,*)'srrrrrrr 2',nqmax,iclass
c$$$      write(*,*)'srrrrrrr 3',sum(ppb)
c$$$      write(*,*)'srrrrrrr 4',nlnmv,nlnmc,sum(nblocha),sum(imdim),sum(iatomp)
c$$$      write(*,*)'srrrrrrr 5',mdimx,nlmto,nbloch,nlnmx,nt,natom,nclass,nl,nnc
c$$$      write(*,*)'srrrrr 1',sum(cphim(1:nlmto,1:ntq))
c$$$      write(*,*)'srrrrr 1',sum(cphiq(1:nlmto,1:ntq))
c$$$      write(*,*)'srrrrr 1',sum(ppb)
c$$$      write(*,*)'srrrrr 1',sum(expikt)
c$$$      write(*,*)'srrrrr 1',sum(zzzmel)
c$$$ccccccccccccccccccccccccccccccccccc
      call psicb_v2  (icore,ncore,nqmax,iclass,
     i  dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i  cphiq(1,nqini), 
     i  ppb,  
c     i              nlnmv,nlnmc,mdim,
     i  nlnmv,nlnmc,nblocha,
     i  imdim,iatomp,
     d  mdimx,nlmto,nbloch,nlnmx,nt,nqmax,natom,nclass,
     d  nl,nnc,     
     o  zzzmel)
      if(debug)  write(6,*) ' sxcf_fal2sc: goto psi2bc1'
      if(debug) write(6,'("4 zzzqqbbb111 ",3d13.5)') sum(abs(zzzmel)),sum(zzzmel)


!!   valence
      call psi2b_v2  (nmmax, nqmax,iclass,
     i  dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i  cphim(1,nmini),                !occ    q-rk nband
     i  cphiq(1,nqini),                 !unocc  q    ntq
     i  ppb,  
     i  nlnmv,nlnmc,nblocha,nctot,
     i  imdim,iatomp,
     d  mdimx,nlmto,nbloch,nlnmx, nband, nt,nqmax,
     d  natom,nclass,
     o  zzzmel)
      if(verbose()>50) call timeshowx("4 after psi2bc1")
      if(debug) then
         write(6,"('sumof zmel abszmel=',4d23.16)")sum(zzzmel),sum(abs(zzzmel) )
      end if
      if(debug) write(6,'("4 zzzqqbbb111 ",3d13.5)')sum(abs(zzzmel)),sum(zzzmel)

!TIME1 "bfore psi2b_v2"
!TIME0
!! -- IPW part.
      if(debug) write(6,*) ' sxcf_fal1: goto drvmelp2 xxx111'
      allocate(drealzzzmel(nbloch,nt,nqmax),dimagzzzmel(nbloch,nt,nqmax))
      drealzzzmel=dreal(zzzmel)
      dimagzzzmel=dimag(zzzmel)
      deallocate(zzzmel)
!    qk =  q - rkvec   
      itps = nqini
      call drvmelp2( q,             nqmax, ! q in FBZ
     i  qk,  nmmax,    ! q-rk
     i  kvec,        ! k in IBZ for mixed product basis. rk = symope(kvec) 
     i  ispq,ispm,ginv,
     i  ngc,ngcmx, ngpmx,nband,itq, 
     i  symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i  drealzzzmel, dimagzzzmel, nbloch, nt,nctot,
     o  rmelt,cmelt, 
     i  itps)
      if(debug) write(6,*) ' sxcf_fal1: end of drvmelp2 sum rmelt cmelt',sum(rmelt),sum(cmelt)
      deallocate(drealzzzmel,dimagzzzmel)
      if(verbose()>50) call timeshowx("5 after drvmelp")
      if(nbcut/=0.and.(.not.exchange)) then
         do it= nctot+1,nctot+min(nbcut,nmmax)
            rmelt(:, it,:) =0d0
            cmelt(:, it,:) =0d0
         enddo
      endif
!TIME1 "after drvmelp2"

!! NOTE:=========================================
!! zmelt = rmelt(igb(rkvec), iocc(q), iunocc(q-rkvec)) + i* cmelt
!! iunocc: band index at target  q.
!! iocc:   band index at intermediate vector qk = q - rkvec 
!! igb: index of mixed product basis       at rkvec (or written as rk)
!!   igb=1,ngb 
!!   ngb=nbloch+ngc  ngb: # of mixed product basis
!!                   nbloch: # of product basis (within MTs)
!!                   ngc: # of IPW for the Screened Coulomb interaction.
!!                   igc is for given 
!! See readgeig in drvmelp2.
!! ==================================================
c-------------------------------------------------------------------------
c$$$!! smbasis 
c$$$!! smbasis ---need to fix this
!!  Read pomatr
c$$$      if(smbasis()) then  !this smbasis if block is from hsfp0.sc.m.F
c$$$         write(6,*)' smooth mixed basis : augmented zmel'
c$$$         ifpomat = iopen('POmat',0,-1,0) !oct2005
c$$$         nkpo = nqibz+nq0i
c$$$         nnmx=0
c$$$         nomx=0
c$$$         do ikpo=1,nkpo
c$$$            read(ifpomat) q_r,nn_,no,iqx !readin reduction matrix pomat
c$$$            if(nn_>nnmx) nnmx=nn_
c$$$            if(no>nomx) nomx=no
c$$$            allocate( pomat(nn_,no) )
c$$$            read(ifpomat) pomat
c$$$            deallocate(pomat)
c$$$         enddo
c$$$         isx = iclose("POmat")
c$$$         ifpomat = iopen('POmat',0,-1,0) !oct2005
c$$$         allocate( pomatr(nnmx,nomx,nkpo),qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
c$$$         do ikpo=1,nkpo
c$$$            read(ifpomat) qrr(:,ikpo),nn_,no,iqx !readin reduction matrix pomat
c$$$            nnr(ikpo)=nn_
c$$$            nor(ikpo)=no
c$$$            read(ifpomat) pomatr(1:nn_,1:no,ikpo)
c$$$         enddo
c$$$         isx = iclose("POmat")
c$$$         write(6,*)"Read end of POmat ---"
c$$$      endif
c-------------------------------------
c$$$               if(smbasis()) then !
c$$$                  ntp0= nqmax
c$$$                  nn= nnr(kx)
c$$$                  no= nor(kx)
c$$$                  allocate( pomat(nn,no) )
c$$$                  pomat= pomatr(1:nn,1:no,kx)
c$$$                  if( sum(abs(kvec-qrr(:,kx)))>1d-10 .and.kx <= nqibz ) then
c$$$                     call rx( 'qibz/= qrr')
c$$$                  endif
c$$$                  if(no /= ngb.and.kx <= nqibz) then
c$$$!!     A bit sloppy check only for kx<nqibz because qibze is not supplied...
c$$$                     write(6,"(' q  ngb  ',3d13.5,3i5)")  kvec,ngb
c$$$                     write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
c$$$                     call rx( 'x0kf_v2h: POmat err no/=ngb')
c$$$                  endif
c$$$                  if(timemix) call timeshow("xxx2222 k-cycle")
c$$$                  ngb = nn      ! Renew ngb !!!
c$$$                  allocate ( zmel  (nn, nctot+nmmax, ntp0) )
c$$$                  call matm( pomat, dcmplx(rmelt,cmelt), zmel, 
c$$$     &                 nn, no, (nctot+nmmax)*ntp0 )
c$$$                  deallocate(rmelt, cmelt)
c$$$                  allocate( rmelt(ngb, nctot+nmmax, ntp0), !ngb is reduced.
c$$$     &                 cmelt(ngb, nctot+nmmax, ntp0) )
c$$$                  rmelt = dreal(zmel)
c$$$                  cmelt = dimag(zmel)
c$$$                  deallocate(zmel,pomat)
c$$$               else
c$$$                  nn=ngb
c$$$                  no=ngb
c$$$               endif

c      if( oncew() ) then
c         write(6,"('ngb nn no=',3i6)") ngb,nn,no
c      endif
c               if(timemix) call timeshow("22222 k-cycle")
      if(allocated(zzzmel))deallocate(zzzmel) !rmel,cmel)
      if(debug) write(6,*) ' sxcf: goto wtt'
      if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)
!! === End of zmelt ; we now have matrix element zmelt= rmelt + img* cmelt === 

      if(exchange) then
         if(debug) write(*,*) 'exchange mode 0000'
         allocate( zmel (ngb, nctot+nmmax, nqmax))
         zmel = dcmplx (rmelt,cmelt)
         if(debug) write(*,*) 'exchange mode 1111'
         deallocate(rmelt,cmelt)
         if(debug) then
            do it=1,nctot+nmmax
               write(6,"('wwwwwsc ',i5,2f10.4)") it,sum(abs(zmel(:,it,1)))
            enddo
            write(*,*) 'eeeeeeeeeeeee end of wwwwsc',nctot,nmmax
            write(6,*)'sumcheck ppovlz=',sum(abs(ppovlz(:,:)))
c            write(6,*)'sumcheck vcoud=',sum(abs(vcoud(:)))
         endif
c                  if(.not.newaniso() ) then
c                     allocate(vcoult(1:ngb,1:ngb),z1p(ngb,nctot+nmmax,nqmax))
c                     vcoult= transpose(vcoul)
c                     call matm( vcoult, zmel, z1p, ngb,ngb,(nctot+nmmax)*nqmax )
c                     deallocate(vcoult) !, zmel)
c                  endif  
!TIME1 "end of matm vcoult"
!TIME0

c$$$                  if(newaniso()) then
c                    write(*,'(a,5I10)')'kino: nqmax,nctot+nmmax,ngb=',nqmax,nctot+nmmax,ngb
         allocate(zmeltt(nctot+nmmax,nqmax,ngb))
c#ifdef USE_GEMM_FOR_SUM
         if(verbose()>39) then
            write(*,*)'info: USE GEMM FOR SUM (zmeltt=zmel*ppovlz) in sxcf_fal2.sc.F'
            write(*,*)'zgemmsize',nqmax*(nctot+nmmax),ngb,ngb
            write(*,*)'size ,zmel',size(zmel,dim=1),size(zmel,dim=2),size(zmel,dim=3)
            write(*,*)'size ,ppovlz',size(ppovlz,dim=1),size(ppovlz,dim=2)
            write(*,*)'size ,zmeltt',size(zmeltt,dim=1),size(zmeltt,dim=2),size(zmeltt,dim=3)
         endif
         call flush(6)	
         call zgemm('T','N',nqmax*(nctot+nmmax),ngb,ngb,(1d0,0d0),
     .     zmel,ngb,ppovlz,ngb,(0d0,0d0),zmeltt,nqmax*(nctot+nmmax) )
         deallocate(zmel)

      else
!! correlation case. Get zmel
c$$$         if(newaniso()) then
         if(debug) write(*,*) 'correlation mode 0000'
         nstate = nctot + nmmax ! = nstate for the case of correlation
         allocate(zmeltt(ngb, nstate, nqmax))
         zmeltt= dcmplx (rmelt,-cmelt) !zmeltt= <itp|it,ib>
         deallocate(rmelt,cmelt)
!! zmel(igb,it*itp) = C(ppovlz)*N(zmeltt(:,it*itp))
!!  C means Hermitian conjugate, N means normal
!! http://www.netlib.org/lapack/lapack-3.1.1/html/zgemm.f.html
         allocate( zmel (ngb, nstate, nqmax) )
         if(debug) write(6,'("4 zzzppp222aaa ",3d13.5)') sum(abs(zmeltt)),sum(zmeltt)
         call zgemm('C','N',ngb, nstate*nqmax,ngb,(1d0,0d0),
     .     ppovlz, ngb, zmeltt,ngb, (0d0,0d0),zmel,ngb)
         deallocate(zmeltt)
         if(debug) write(*,*)'zz000 nqmax,ngb,nstate ',nqmax,ngb,nstate
         if(debug) write(*,*)'zz000 sumchk zmel ',sum(abs(zmel(1:ngb,1:nstate,1:nqmax)))
         if(debug) write(*,*) 'correlation mode end'
!TIME1 "matmul zmel ppovlz"
      endif         
      end subroutine get_zmelt2
      end module m_zmel

      subroutine timeshowx(info)
      character*(*) :: info
      write(*,'(a,$)')info
      call cputid(0)
      end

