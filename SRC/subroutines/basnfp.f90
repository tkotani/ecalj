subroutine basnfp_v2(nocc,nunocc,nindx, lmxax,nn,nrx,nrofi,r,aa,bb,ic, & !Generate product basis siwhin MT
       phitoto,phitotr,nsp,nclass, cutbase,lcutmx,ixx,alat,nc_max)
  use m_keyvalue,only:getkeyvalue
  use m_ll,only: ll
  use m_read_bzdata,only: Read_bzdata, q0i,nq0i,wqt=>wt
  use m_genallcf_v3,only: lmxa
  !i phitotr: atomic radial functions. raw
  !r       This is avereged as for spin (stored in phiav) and used in the construction of product basis.
  !i phitoto: atomic radial functions. ortogonalized
  !r       Coefficents of eigenfunctions based on this orthogonalized phi.
  !r       See cphix, which is written to CphiGeig in rdata4gw_v2.m.f
  !r       Roughly speaking, the eigen function psi = \sum cphix * phitoto
  !r
  !r    core part is the same in phitoto and phitotr
  !i   indexes, and cutoff conditions to set up product basis.
  !     nocc(l,n)   = 0 or 1
  !     nunocc(l,n) = 0 or 1
  !o output files
  !o   BASFP//atom     :  product basis. used in hvccfp0.
  !o   PPBRD_V2_//atom : <phi phi |B> radial integral. Read in hx0fp0 hsfp0 through rdpp_v2.
  !o
  implicit none
  character(8) :: xtxx
  integer,parameter:: nxxmx=300, npradmx=300
  integer:: nn,nrx, lmxax, nocc(0:lmxax,nn), nunocc(0:lmxax,nn), nindx(0:lmxax), &
       npr(0:lmxax,nn,0:lmxax,nn),nprpd(0:lmxax,nn,0:lmxax,nn), &
       iprad,l1,l2,n1,n2,lx,nx, ifix,icx
  real(8):: phi(nrx,0:lmxax,nn) 
  integer:: nrofi, i,kmax,nx1,nx2,ib1,ib2,iadd,ibx, &
       nxx( 0:2*(lmxax) ),nprad, nprod, &
       nxxold( 0:2*(lmxax) ), ngmxx( 0:2*(lmxax) ),nodnum,nbasen
  real(8)::    aa,bb,tot, hnrofi,rmax, &
       r(nrx),sig,sig0,ovv,rxx(nrofi) ,epsx,sxx,rrt, &
       alpha,polinta,cutbase(0:2*(lmxax)),adist,rax, plgndr
  real(8):: screen(nrofi),screent(nrofi),kappa,axx,epp   ,f0(nrofi)
  real(8),allocatable:: ovvv(:,:),ovvi(:,:),sc(:,:),oso(:,:), &
       rkp (:,:), rkm (:,:), &
       rkp0(:,:), rkm0(:,:), &
       rkpt(:,:), rkmt(:,:), &
       scrnmt(:,:,:),ovmt(:,:,:)  , eb(:)
  real(8),allocatable:: zz(:,:,:),rprodtc(:,:),ovvc(:,:,:),oc(:,:,:),wk(:),tc(:)
  real(8),parameter:: fpi = 4d0*3.14159265358979323846d0
  real(8),parameter::  pi =     3.14159265358979323846d0
  integer:: nb,ngmx,ig,ificrb,kmx,isx,k,ic,nblocha
  integer,allocatable:: iwk(:)
  integer,allocatable:: iprlc(:),lprc(:),ibo(:)
  real(8),allocatable::  rprodx(:,:,:),rprodx2(:,:,:),rprod(:,:)
  logical :: newbase2
  real(8) :: bbase(nrofi),absqg2,aaa,aaa12,rphiphi(nrofi)
  real(8),allocatable ::phij(:),psij(:)
  integer::lxx,nxxx,ir,n,l,ierr,   lcutmx,lxxm
  integer::nsp,nclass,ifppb,isp,ip1,ip2
  real(8) :: phitoto(nrx,0:lmxax,nn,nclass, nsp)
  real(8) :: phitotr(nrx,0:lmxax,nn,nclass, nsp)
  real(8),allocatable :: ppbrd(:,:,:,:,:,:),absqg2x(:) !,wqt(:), q0i(:,:)
  real(8),allocatable :: phiav(:,:,:)
  real(8) :: rnormphi(0:lmxax,nn) ,alat,sss
  integer :: irx,noo,nuu      ,ix !,nq0i, neps, nq0ix
  integer:: zvztest,nlmlsp
  real(8):: ovvs,rrr(nrx)
  !---for mode 8
  integer:: nr_r,nlml_r,nsp_r,ifv,ibas,ilmx,isp1,isp2,kxx,ixx,nxx_i
  real(8),allocatable:: rofi_r(:),rho1(:,:,:),rspin(:),den(:,:),r11(:)
  real(8):: spinvec,spinvec0,sumc(2),sqrtfpi,const
  character(3):: charnum3
  logical:: valmt=.false.
  logical :: newaniso,addbasnew, read_bzdata_done=.false. !smbasis,
  character(15) :: rprodf
  integer :: ifprodmt, l1l2p,l1l2m,inn !smbasiscut,
  integer :: lxlnln(0:2*(lmxax), 0:lmxax,nn,0:lmxax,nn) ,nzz
!  real(8),allocatable:: prodmt (:,:,:)
  real(8):: prr(3),derie,derie2,ddd,derie3
  integer:: nc_max(0:lmxax),nnn !, smbasis_case
  character(len=100):: recxxx
  character(len=160):: recxxx2
  integer:: npbasmax(0:2*(lmxax)),ifinin,iax,izz,naxx,verbose
  integer,allocatable:: ipx(:,:)
  real(8):: bb1,bb1s,aa_in,bb_in
  integer:: nl2m1,nrofi_in,lxx_
  print *,' basnfp_v2: ********** start ******** nrofi=',nrofi
  allocate(phiav(nrx,0:lmxax,nn),rprod(1:nrx, npradmx),ipx(0:2*(lmxax),nxxmx))
  ! To make the sign of phitotr safer; but I may already added something for keeping sss=1d0, right?
  do l1 = 0, lmxax
     do n1 = 1, nindx(l1)
        if(nsp==1) then
           phiav(1:nrx, l1, n1) = phitotr(1:nrx, l1, n1, ic, 1)
        else
           sss = dsign(1d0, phitotr(2,l1,n1,ic, 1)*phitotr(2,l1,n1,ic, 2))
           phiav(1:nrx, l1, n1) = &
                ( phitotr(1:nrx, l1, n1, ic, 1) + sss* phitotr(1:nrx, l1, n1, ic, 2) )/2d0
        endif
     enddo
   enddo
   
  lxx = 2*(lmxax)
  ptestixx4: if(ixx==4) then  ! ccccccccc start of ptest mode cccccccccccccccccccccccccccccccccccccccccccc
    if( .NOT. read_bzdata_done) then
      call read_bzdata()
      read_bzdata_done=.true.
    endif
    print *, ' *** rprodx is given by Bessel.***'
    allocate(phij(0:lxx),psij(0:lxx))
    allocate( absqg2x(nq0i) ) 
    nzz = max(2,nq0i)
    ix=0
    do i=1,nq0i
      if(wqt(i)==0d0 ) then
        ix=ix+1
        absqg2x(ix) =sum( (2*pi/alat *q0i(1:3,i))**2) !nq0i ---> q0i
        if(ix>1) call addd(absqg2x(ix),absqg2x,ix-1)
      endif
    enddo
    fac2l: block
      integer:: fac2l(0:lxx)
      fac2l(0) = 1d0
      do  lx = 1, lxx
        fac2l(lx) = fac2l(lx-1) * (lx+lx-1)
      enddo
      nxxx = ix
      iprad = 0
      nxx = 0
      do lx = 0, lxx
        if(lx >lcutmx) cycle  ! Lmax cutoff for product basis
        do n1 = 1, nxxx
          !           if(smbasis() .AND. lx>  smbasiscut() ) cycle
          iprad  = iprad + 1
          nxx(lx)= nxx(lx)+1
          ipx(lx, nxx(lx)) = iprad
          absqg2 = absqg2x(n1)
          do ir =1,nrofi
            lxxm=max(2,lx)
            call bessl(absqg2*r(ir)**2,lxxm,phij,psij)!lmin must be larger than 1,right?
            phij(0:lxxm) = phij(0:lxxm)*fac2l(0:lxxm)*0.5d0 !Andersen factor
            !              psij(0:lxx) = psij(0:lxx)/fac2l(0:lxx)       !Andersen factor
            rprod(ir,iprad) = phij(lx)* r(ir) **(lx +1)
          enddo
          print *,' sumchk rprod=',lx,n1,sum(abs(rprod(1:nrofi,iprad)))
          if(.false.) then !verbose()>60) then
            write(3100+ic,"(' -- -- -- ',3i3,' --- ' )") lx,n1
            do ir =1,nrofi
              write(3100+ic,"(d13.5,2x,2d18.8)") r(ir), rprod(ir,iprad)
            enddo
          endif
        enddo
      enddo
    end block fac2l
    nprad=iprad
    print *, ' *** TEST nprad=',nprad
    print *, ' nxx =',nxx(0:lxxm)
    goto 1212
  endif ptestixx4
  
  normcheck: do l1 = 0, lmxax
     do n1 = 1, nindx(l1)
        call gintxx(phiav(1,l1,n1),phiav(1,l1,n1),aa,bb,nrofi,ovv)
        write(6,"(' norm check for phi='2i3,d13.6)") l1,n1,ovv
        if(abs(ovv) <1d-10 ) ovv  = 0d0
        rnormphi(l1,n1) = ovv
     enddo
  enddo normcheck
   
  ! product basis construction
  lxlnln=0
  npr    = 0
  nprpd  = 0
  nxx    = 0
  nxxold = 0
  iprad  = 0
  rprod  = 0d0
  do 203 l1 = 0, lmxax
    do 202 n1 = 1, nindx(l1)  ; noo = nocc  (l1,n1)
      do 201 l2 = 0, lmxax
        do 20 n2 = 1, nindx(l2)  ; nuu = nunocc(l2,n2)
          !!  phiav * phiav ----
          if( noo>=1 .AND. nuu>=1 ) then
            if( npr(l1,n1,l2,n2) == 0 ) then
              iprad = iprad + 1
              npr (l1,n1,l2,n2) = iprad
              npr (l2,n2,l1,n1) = iprad
              do lx = abs(l1-l2), l1+l2
                if(mod(lx+l1+l2,2)==1) cycle
                if(lx >lcutmx) cycle  ! Lmax cutoff for product basis
                write(6,"(' ---  lx l1 l2 n1 n2 iprad= ',6i5)") lx,l1,l2, n1, n2,iprad
                nxx(lx) = nxx(lx) + 1
                ipx(lx, nxx(lx)) = iprad
                lxlnln(lx,l1,n1,l2,n2)=iprad
                lxlnln(lx,l2,n2,l1,n1)=iprad
              enddo
              rprod(1,iprad) = 0d0
              rprod(2:nrofi,iprad) = phiav(2:nrofi,l1,n1)*phiav(2:nrofi,l2,n2)/r(2:nrofi) ! phi = u = r \phi
              !          call gintxx(phiav(1,l1,n1), phiav(1,l2,n2),aa,bb,nrofi, sss )
              !          write(6,"(' normchk phiav*phiav=',4i3,d13.6)") l1,n1,l2,n2,sss
            endif
          endif
20      enddo
201   enddo
202 enddo
203 enddo

  !! Additional product bais for smbais()=T.
!   if(smbasis()) then ! Dec2005
!      !$$$!!--- CASE1. Add r^l r^(l+2) anyway. ---
!      if(smbasis_case()==1) then
!         do 110 lx  = 0, 2*(lmxax)
!            if(lx > smbasiscut() ) cycle
!            do inn = 1, 2
!               iprad = iprad + 1
!               rprod(1,iprad) = 0d0
!               if(inn==2) nnn = lx + 1
!               if(inn==1) nnn = lx + 3
!               rprod(2:nrofi,iprad)= r(2:)**nnn
!               nxx(lx) = nxx(lx) + 1
!               ipx(lx, nxx(lx)) = iprad
!               write(6,"('sm -- lx iprad nnn nxx= ',6i5)")lx, iprad,nnn,nxx(lx)
!            enddo
! 110     enddo
!      elseif(smbasis_case()==2) then
!         !!--- CASE2 ---
!         do 111 lx  = 0, 2*(lmxax)
!            if(lx > smbasiscut() ) cycle
!            nxx_i =nxx(lx)
!            do inn = nxx_i+1, 2
!               iprad = iprad + 1
!               rprod(1,iprad) = 0d0
!               if(inn==2) nnn = lx + 1
!               if(inn==1) nnn = lx + 3
!               rprod(2:nrofi,iprad)= r(2:)**nnn
!               nxx(lx) = nxx(lx) + 1
!               ipx(lx, nxx(lx)) = iprad
!               write(6,"('sm -- lx iprad nxx= ',6i5)")lx, iprad,nxx(lx)
!            enddo
! 111     enddo
!      elseif(smbasis_case()==3) then
!         !!--- CASE3 ---
!         !! GaAs local-orbital tests suggest the above choice looks better for Ga 3d core.
!         !!  High priority to low priority
!         !! Priority ordering
!         !!     1. smaller l_1 +l_2,
!         !!     2. smaller |l1-l2|
!         !!     3.  phi*phi,  phi*phidot, phi*local
!         do 112 lx  = 0, 2*(lmxax)
!            if(lx>  smbasiscut() ) cycle
!            inn=0
!            l1= lx/2
!            l2=  lx-l1   !l1<=l2
!            !       Priority ordering is phi*phi,  phi*phidot  !, phi*local
!            do 120 n1 = nc_max(l1) +1,nc_max(l1) +1  !nindx(l1)
!               do 130 n2 = nc_max(l2) +1, nindx(l2)
!                  inn = inn+1
!                  if( lxlnln(lx,l1,n1,l2,n2)==0) then
!                     if( npr(l1,n1,l2,n2) == 0 ) then
!                        iprad = iprad + 1
!                        npr (l1,n1,l2,n2) = iprad
!                        npr (l2,n2,l1,n1) = iprad
!                        rprod(1,iprad) = 0d0
!                        rprod(2:nrofi,iprad)=phiav(2:,l1,n1)*phiav(2:,l2,n2)/r(2:) ! phi = u = r \phi
!                     else
!                        iprad = npr(l1,n1,l2,n2)
!                     endif
!                     nxx(lx) = nxx(lx) + 1
!                     ipx(lx, nxx(lx)) = iprad
!                     lxlnln(lx,l1,n1,l2,n2)= iprad
!                     lxlnln(lx,l2,n2,l1,n1)= iprad
!                     write(6,"('sm -- lx l1 l2 n1 n2 iprad= ',6i5)")lx, l1, l2,n1, n2,iprad
!                  endif
!                  if(inn==2) cycle
! 130           enddo
! 120        enddo

! 112     enddo
!      endif
     !      print *, '======= basis set up =============='
     !      do lx  = 0, 2*(lmxax)
     !      do nx = 1, nxx(lx) ;  ib1 = ipx(lx,nx)
     !        nprod = nprod + 2*lx+1
     !        write(6,"('lx nx number ipx=',4i4)") lx, nx, 2*lx+1, ipx(lx,nx)
     !      enddo
     !      enddo
     !      print *, '--- test end ---------'
     !      return
!  endif
  !! === Add a PW like product basis  ===
  ! Here you can add any functions in the same way.
      !  if(smbasis() .AND. smbasis_case()==1) goto 113
      
  if( .NOT. addbasnew() .AND. zvztest()/=2) then !Add one s type product basis.
     iprad  = iprad  + 1 ! radial function index
     nxx(0) = nxx(0) + 1 ! nxx(l) is the number index of the radial functions for l.
     ipx(0, nxx(0)) = iprad ! the radial function index for each l and nxx(l).
     rprod(1:nrofi,iprad) = r(1:nrofi)
     ! this is better (but larger ProductBasis). Used for the check on Coulomb matrix at q->0.
  elseif(zvztest()/=2) then !followings may be a compromise to give reasoably good enough <M|v|M>.
     do lx = 0,min(lcutmx,1) ! may2015)  !(=0,lcutmx before feb2012) !(=0,1  !feb2012)
        iprad  = iprad  + 1 ! radial function index
        nxx(lx) = nxx(lx) + 1
        ipx(lx, nxx(lx)) = iprad
        rprod(1:nrofi,iprad) = r(1:nrofi)*r(1:nrofi)**lx
        if(lx>=2) cycle !if(lx>=2) cycle !to reduce number of basis basis.
        iprad  = iprad  + 1 ! radial function index
        nxx(lx) = nxx(lx) + 1
        ipx(lx, nxx(lx)) = iprad
        rprod(1:nrofi,iprad) = r(1:nrofi)*r(1:nrofi)**(lx+1)
     enddo
  endif
113 continue
  nprad = iprad
  print *,' number of radial basis: nprad=', nprad
1212 continue

  ! sanity check
  if(maxval(nxx(0:lxx))>nxxmx) call rx(' basnfp: nxx >nxxmx --- Enlarge nxxmx in basnfp.f')
  if(nprad > npradmx)call rx( ' basnfp: nprad > npradmx --- Enlarge npradmx in basnfp.f')
  !! nprod
  nprod = 0
  do lx = 0, lxx !2*(lmxax)
     do nx = 1, nxx(lx)
        nprod = nprod + 2*lx+1
        write(6,"('lx nx number ipx=',4i4)") lx, nx, 2*lx+1, ipx(lx,nx)
     enddo
  enddo
  print *,' *** total number of product basis nprod=', nprod
  !! ovmt
  kmax = lxx !2*(lmxax)
  allocate( ovmt  (nprad,nprad,0:lxx)) !2*(lmxax))  )
  do lx  = 0, lxx !2*(lmxax)
     do nx1 = 1, nxx(lx) ;  ib1 = ipx(lx,nx1)
        do nx2 = 1, nxx(lx) ;  ib2 = ipx(lx,nx2)
           call gintxx(rprod(1,ib1),rprod(1,ib2),aa,bb,nrofi, ovv)
           ovmt(ib1,ib2, lx) = ovv !radial overlap matrix of rprod
        enddo
     enddo
  enddo
  !! orthonormal basis functions.--------------------
  allocate(rprodx(nrofi,maxval(nxx(0:lxx)),0:lxx))
  do lx = 0, lxx
     nb = nxx(lx)
     if(nb==0) cycle
     print *
     allocate( ovvc(nb,nb,2), zz(nb,nb,2),eb(nb),ibo(nb) )
     ovvc = 0d0;  zz = 0d0; eb = 0d0
     do ib1 = 1,nb
        do ib2 = 1,nb
           ovvc(ib1,ib2,1) = ovmt(ipx(lx,ib1),ipx(lx,ib2),lx)
        enddo
     enddo
     call rss(nb, ovvc(:,:,1), eb, zz(:,:,1),ierr) !rs==>rss 2022-6-13
     if(ierr/=0) call rx( ' basnfp: rs error')
     ibx=0
     do ib1 = nb,1,-1
        write(6,"(a,i5,d13.6,a,d13.6)")'    ib eb=',ib1,eb(ib1),'  ecut=',cutbase(lx)
        if(eb(ib1)<cutbase(lx)) then !this is mainly used now 2022jan
           cycle
        endif
        !          endif
        ibx = ibx+1
        ibo(ibx) = ib1
        rprodx(1:nrofi,ibx,lx) = matmul( rprod(1:nrofi,ipx(lx,1:nb)), zz(1:nb,ib1,1) )
        rprodx(1:nrofi,ibx,lx) = rprodx(1:nrofi,ibx,lx)/sqrt(eb(ib1))
     enddo
     nb     =ibx
     nxx(lx)=ibx
     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! debug
     !      if(lx==2) then
     !        ibx =1
     !        print *,' ibx=', ibx, eb(ibx)
     !        do ib1 =1,nb
     !          call gintxx(rprod(1,ipx(lx,ib1)),rprodx(1,ibx,lx)
     !     &              ,aa,bb,nrofi,ovv)
     !          ovv = ovv*sqrt(eb(ibx))
     !          write(6,"(' zz xxx',4d16.7)")
     !     &    zz(ib1,ibx,1), ovv,   !zz is normalised to 1.
     !     &    ovv/ ( zz(ib1,ibx,1)*eb(ibx) )
     !        enddo
     !        print *,' norm=',sum(zz(1:nb,ibx,1)**2)
     !      stop
     !      endif
     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     do ib1=1,nb
        do ib2=1,nb
           call gintxx(rprodx(1,ib1,lx),rprodx(1,ib2,lx),aa,bb,nrofi,ovv)
           if(ib1==ib2 ) then
              write(6,"('    Diag ibx ovv=',i3,d23.16,' eb= ',d18.10,' nod= ',i3)") &
                   ib1,ovv, eb(ibo(ib1)), nodnum(rprodx(1,ib1,lx),nrofi)
           elseif( ovv>1d-10 ) then
              write(6,"('      offdiag ib1 ib2 ovv=',2i3,d23.16,' eb=',d18.10)") &
                   ib1,ib2,ovv, eb(ibo(ib1))
           endif
        enddo
     enddo
     write(6,"('    *** lx  =',i4,'*** Used nb =',i4)") lx,nb
     deallocate( ovvc, zz, eb,ibo)
  enddo

  !- Reserve rprodx
  !      real(8):: crbase(nrx,kmxx,nclass)
  !      integer:: kmxx,kmx(nclass),iprlc(kmxx,nclass),lprc(kmxx,nclass)
  !----------------------------------------------
  kmx = sum( nxx(0:lxx) )
  allocate( iprlc(kmx), lprc(kmx) )
  k = 0
  iprlc(1:kmx) = 0
  do lx = 0, lxx !2*(lmxax)
     do nx = 1, nxx(lx)
        k = k + 1
        if( k==1 ) then
           iprlc(k) = 0
        else
           iprlc(k) = iprlc(k-1) + 2*lprc(k-1)+1
        endif
        lprc (k) = lx   ! irdc(k,ic) = k
     enddo
  enddo
  if(k/=kmx) call rx( ' basnfp: k/=kmx')
  if(kmx/=0 ) then
     nblocha = iprlc (kmx) + 2*lprc(kmx)+1
  elseif(kmx==0 ) then
     nblocha = 0
  endif
  
!!!!!!!!!!!!!!!!!!!!!!----------------------
!  newlxx: do l=0,lxx
!    if(nxx(l)>0) lxx_=l
!  enddo newlxx
!  lxx = lxx_ !lxx is ic dependent 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  print *,' Write BASFP.* reserve rprodx...'
  write(6,"(' basnfp: BASFP... kmx nblocha=',2i5)") kmx,nblocha
  open(newunit=ificrb,file='__BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10) ))
  write(ificrb,"(4i6,2d24.16)") lxx, kmx, nblocha,nrofi,aa,bb
  write(ificrb,"(i5)") nxx(0:lxx)
  k = 0
  do lx = 0, lxx
     do nx = 1, nxx(lx)
        k = k + 1
        write(ificrb,"(3i5)"   ) k,iprlc(k),lprc(k)
        write(ificrb,"(d23.15)") (rprodx(i,nx,lx),i=1,nrofi)
     enddo
  enddo
  close(ificrb)
  
!  deallocate(rprodx,iprlc,lprc)
!  
!  Write and Read even for iread==0, in order to have the exact match on PPB* when iread=1 for the same BASFP* file.
!  print *,' read rprodx...'
!  open(newunit=ificrb,file='__BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10) ))
!  read(ificrb,"(4i6,2d24.16)") lxx, kmx,nblocha,nrofi_in,aa_in,bb_in
!  read(ificrb,"(i5)") nxx(0:lxx)
!  allocate(rprodx(nrofi,maxval(nxx(0:lxx)),0:lxx))
!  allocate( iprlc(kmx), lprc(kmx) )
!  k = 0
!  do lx = 0, lxx !2*(lmxax)
!     do nx = 1, nxx(lx)
!        k = k + 1
!        read(ificrb,"(3i5)"   ) k,iprlc(k),lprc(k)
!        read(ificrb,"(d23.15)") (rprodx(i,nx,lx),i=1,nrofi)
!     enddo
!  enddo
!  close(ificrb)

  ! Calculate radial matrix elements.
  print *,' Calculate radial matrix elements...'
  newlxx: do l=0,lxx
    if(nxx(l)>0) lxx_=l
  enddo newlxx
  lxx=lxx_ !lxx is ic dependent 
  allocate( ppbrd(0:lmxax,nn,0:lmxax,nn, 0:lxx, maxval(nxx(0:lxx)) ) )
  ppbrd =.9999d99 !for safe
  open(newunit=ifppb,file='__PPBRD_V2_'//char( 48+ic/10 )//char(48+mod(ic,10)),form='unformatted')
  write(ifppb) nblocha, lxx
  write(ifppb) nxx(0:lxx)
  ppbrddo: do isp= 1,nsp
    isp1=isp
    isp2=isp
    if(ixx==8 .AND. isp==1) then
      isp1=1
      isp2=2
    elseif(ixx==8 .AND. isp==2) then
      isp1=2
      isp2=1
    endif
    do  lx = 0, lxx
      do  nx = 1, nxx(lx)
        do  l1 = 0, lmxa(ic) !lmxax
          do  n1 = 1, nindx(l1)
            do  l2 = 0, lmxa(ic) !lmxax
              do  n2 = 1, nindx(l2)
                if(lx <abs(l1-l2) .OR. l1+l2<lx) cycle
                rphiphi(1)       = 0d0
                rphiphi(2:nrofi) = phitoto(2:nrofi,l1,n1,ic,isp2) *phitoto(2:nrofi,l2,n2,ic,isp1)/r(2:nrofi) ! phi = u = r \phi
                call gintxx(rprodx(1,nx,lx), rphiphi,aa,bb,nrofi, ppbrd(l1,n1,l2,n2,lx,nx) )
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    write(ifppb) ppbrd
  enddo ppbrddo
  deallocate(ppbrd)
  close(ifppb)

  ! prodmt proddmt (value and slope of the product at MT). !oct2005
  !    Stored into the tail of PPBRD* files.
  ! if(smbasis()) then
  !    allocate( prodmt (2,maxval(nxx(0:2*(lmxax))),0:2*(lmxax)))
  !    prodmt = 1d10 !for safe
  !    do lx = 0, 2*(lmxax)
  !       do nx = 1, nxx(lx)
  !          prodmt(1,nx,lx)= rprodx(nrofi,nx,lx)/r(nrofi)
  !          prodmt(2,nx,lx)= derie2(r, rprodx(1:nrofi,nx,lx)/r(1:nrofi), nrofi)
  !       enddo
  !    enddo
  !    open(newunit=ifprodmt,file='PRODMT_'//charnum3(ic),form='unformatted')
  !    write(ifprodmt) nl
  !    write(ifprodmt) maxval(nxx(0:2*(lmxax)))
  !    write(ifprodmt) nxx(0:2*(lmxax))
  !    write(ifprodmt) prodmt
  !    close(ifprodmt)
  !    deallocate(prodmt)
  ! endif

  ! --- MixSpin= <rho_up - rho_down | B> matrix calculation. May2005
  ! 1  Suppose  "ibas==iclass"--- it is already checked in hbasfp0.m.f
  ! 2  ValMT.* is written with subroutine savemtval(ib,rho1,rofi,nr,nlml,nsp)
  !r    in fp/locpot.f just befor locpt2 in lmto (lmf).
  ixx8if: if(ixx==8) then
    sqrtfpi = sqrt(fpi)
    ibas=ic
    if(valmt) then
      open(newunit=ifv,file='ValMT.'//charnum3(ibas)//'.chk',form='unformatted')
      read(ifv) nr_r,nlml_r,nsp_r
      write(6,"('readin nr nlml nsp=',3i5)") nr_r,nlml_r,nsp_r
      allocate(rofi_r(nr_r),rho1(nr_r,nlml_r,nsp_r),rspin(nrofi),den(nrofi,nsp),r11(nrofi))
      r11(1:nrofi)= 1d0
      read(ifv) rofi_r, rho1
      close(ifv)
    else
      open(newunit=ifv,file='rhoMT.'//xtxx(ibas),form='unformatted',status='old',err=1031)
      goto 1032
1031  continue
      call rx( 'rhoMT by locpot-wrhomt. open error')
1032  continue
      read (ifv) nr_r
      allocate(rofi_r(nr_r))
      read (ifv) rofi_r
      read (ifv) nr_r,nlmlsp,kxx,kxx,nsp_r
      write(6,*)' rho1 xxx=', nr_r,nlmlsp,kxx,kxx,nsp_r
      nlml_r = nlmlsp/nsp_r
      allocate( rho1(nr_r,nlml_r,nsp_r),rspin(nrofi),den(nrofi,nsp),r11(nrofi))
      r11(1:nrofi)= 1d0
      read (ifv) rho1
    endif
    if(nsp_r/=nsp) call rx( " ReadinError: ValMT: nspr/= nsp")
    if(nsp/=2    ) call rx( " This mode is only for nsp==2")
    rho1= sqrtfpi*rho1  !rho1 is not including sqrt(fpi) Right?
    open(newunit=ifv,file='MixSpin.'//charnum3(ibas))
    write(ifv,"(2i10,' ! ibas, max l of product basis' )") ibas,lxx !2*(lmxax)
    write(ifv,"(i10,'           ! nxx(lx)'  )") nxx(0:lxx)
    do ilmx = 1,(lxx+1)**2
      lx = ll(ilmx)
      if(ilmx <=nlml_r) then
        rspin(1) = 0d0  !rspin = rho^{true spin density} * r
        do ir =2,nrofi
          den(ir,1)=  polinta(r(ir), rofi_r,rho1(:,ilmx,1),nr_r)
          den(ir,2)=  polinta(r(ir), rofi_r,rho1(:,ilmx,2),nr_r)
          rspin(ir)  = (den(ir,1) -den(ir,2) )/r(ir)
        enddo
        den(1,1:2)=0d0
      else
        rspin=0d0
        den=0d0
      endif
      ! den = 4 pi r^2 * rho_true(r)  !  rspin = 4 pi r * rho_true(r)
      if( nxx(lx)/=0) then
        do isp=1,nsp
          call gintxx( den(1,isp), r11, aa,bb,nrofi, sumc(isp) )
        enddo
        write(6,"(' charge: ilm charge=',i5,2f13.6)") ilmx,sumc(1:nsp)
      endif
      bb1s=0d0
      do nx = 1, nxx(lx)
        call gintxx( rprodx(1,nx,lx), rspin,aa,bb,nrofi,  spinvec0 )
        spinvec = spinvec0/sqrtfpi
        if(lx==0) then ! const = <1|B> where 1 is normalized within the sphere 2007
          call gintxx( rprodx(1,nx,lx), r,aa,bb,nrofi, const )
        else
          const=0d0
        endif
        const= const *  sqrtfpi !/((fpi/3d0)*r(nrofi)**3)
        ! Now spinvec = <B_I(\bfr) | m_true(\bfr) >
        if(abs(spinvec)<1d-10 ) spinvec=0d0
        write(ifv,"(     2i6,d24.16,2x,f13.10,2x,f13.10,d24.16,' ! I=(ilm, nx), <spin|B(I)>, intrho(1:nsp) <1|B(I)>')") &
             ilmx, nx, spinvec, sumc(1:nsp),const
        write(6,"('ttt:',2i6, d24.16,2x, 2f14.10,d24.16,' ! I=(ilm, nx), <spin|B(I)>, intrho(1:nsp) <1|B(I)>')") &
             ilmx, nx, spinvec, sumc(1:nsp),const
      enddo
    enddo
    deallocate(rofi_r,rho1,rspin)
    close(ifv)
  endif ixx8if !ixc==8 end
  if (allocated(rprod))   deallocate(rprod)
  if (allocated(phiav))   deallocate(phiav)
  if (allocated(phij))    deallocate(phij)
  if (allocated(absqg2x)) deallocate(absqg2x)
  if (allocated(ovmt))    deallocate(ovmt)
  if (allocated(rprodx))  deallocate(rprodx)
  if (allocated(rprodx2)) deallocate(rprodx2)
  if (allocated(iprlc))   deallocate(iprlc)
  if (allocated(lprc))    deallocate(lprc)
  return
end subroutine basnfp_v2

real(8) function derie2 (x,y,n)
  intent(in)::x,y,n
  !     return derivative at n
  integer:: n,i
  real(8):: x(n), y(n), dydi(n),ii(n),dydx(n)
  real(8),external::polinta
  do i=1,n-1
     dydx(i) = (y(i+1)-y(i))/(x(i+1)-x(i))
     ii(i) = i + 0.5d0
  enddo
  derie2 = polinta(dble(n), ii,dydx,n-1)
END function derie2

subroutine addd(a1,a,n)
  ! a1 is choosed so that it is not in agreement with other a(n).
  implicit none
  integer:: n,i
  real(8):: eps=0.1,a1,a(n)
880 continue
  do i=1,n
     if(abs(a1-a(i))< 0.2) then
        a1=a1+eps
        goto 880
     endif
  enddo
end subroutine addd

! taken from basn.f
pure integer function nodnum(f,n)
  intent(in)::f,n
  integer::i,n
  real(8):: f(n)
  nodnum=0
  do i=2,n-1
     if(f(i)*f(i+1)<0) nodnum=nodnum+1
  enddo
END function nodnum
