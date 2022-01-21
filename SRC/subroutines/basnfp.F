      subroutine basnfp_v2 (nocc,nunocc,nindx,      !phi,
     &    nl,nn,nrx,nrofi,r,aa,bb,ic, 
     i    phitoto,phitotr,nsp,nclass,
     i    cutbase,lcutmx,ixx,alat,nc_max) !,iread
      use m_keyvalue,only:getkeyvalue
      use m_lldata,only: ll
c takao kotani Apr 2002.
c gives an index for the allowed product basis
c A new routine by t.kotani mod. from indxbas of fa.
c
Ci phitotr: atomic radial functions. raw
Cr       This is avereged as for spin (stored in phiav) and used in the construction of product basis.
Ci phitoto: atomic radial functions. ortogonalized
Cr       Coefficents of eigenfunctions based on this orthogonalized phi.
Cr       See cphix, which is written to CphiGeig in rdata4gw_v2.m.f
Cr       Roughly speaking, the eigen function psi = \sum cphix * phitoto
Cr
cr    core part is the same in phitoto and phitotr
Ci   indexes, and cutoff conditions to set up product basis.
c     nocc(l,n)   = 0 or 1
c     nunocc(l,n) = 0 or 1
Co output files
Co   BASFP//atom     :  product basis. used in hvccfp0.
Co   PPBRD_V2_//atom : <phi phi |B> radial integral. Read in hx0fp0 hsfp0 through rdpp_v2.
Co
      implicit none
      character*8 xtxx
      integer(4),parameter:: nxxmx=300, npradmx=300
      integer(4):: nl,nn,nrx, 
     &  nocc(0:nl-1,nn), nunocc(0:nl-1,nn), nindx(0:nl-1),
     &  npr(0:nl-1,nn,0:nl-1,nn),nprpd(0:nl-1,nn,0:nl-1,nn),
     &  iprad,l1,l2,n1,n2,lx,nx, ifix,icx
      real(8):: phi(nrx,0:nl-1,nn) !nrx = maxval(nr(1:nclass))
c
      integer(4):: nrofi, i,kmax,nx1,nx2,ib1,ib2,iadd,mintc,mxintc,ibx,
     &    nxx( 0:2*(nl-1) ),nprad, nprod,
     &    nxxold( 0:2*(nl-1) ), ngmxx( 0:2*(nl-1) ),nodnum,nbasen
      real(8)::    aa,bb,tot, hnrofi,rmax,
     &     r(nrx),sig,sig0,ovv,rxx(nrofi) ,epsx,sxx,rrt,
     &     alpha,polinta,cutbase(0:2*(nl-1)),adist,rax, plgndr
      real(8):: screen(nrofi),screent(nrofi),kappa,axx,epp   ,f0(nrofi)
      real(8),allocatable:: ovvv(:,:),ovvi(:,:),sc(:,:),oso(:,:),
     &          rkp (:,:), rkm (:,:),
     &          rkp0(:,:), rkm0(:,:),
     &          rkpt(:,:), rkmt(:,:),
     &          scrnmt(:,:,:),ovmt(:,:,:)  , eb(:)
      real(8),allocatable:: zz(:,:,:),rprodtc(:,:) !,oxx(:,:)
     &        ,ovvc(:,:,:),oc(:,:,:),wk(:),tc(:) !,ovvco(:,:,:)
      real(8),parameter:: fpi = 4d0*3.14159265358979323846d0
      real(8),parameter::  pi =     3.14159265358979323846d0
      integer(4):: nb,ngmx,ig,ificrb,kmx,isx,k,ic,nblocha
      integer(4),allocatable:: iwk(:)
      character*7 filename
      character*11 filenamep
      integer(4),allocatable:: iprlc(:),lprc(:),ibo(:)
      real(8),allocatable::  rprodx(:,:,:),rprodx2(:,:,:),rprod(:,:)
      logical newbase2
      real(8) :: bbase(nrofi),absqg2,aaa,aaa12,rphiphi(nrofi)
      real(8),allocatable ::phij(:),psij(:)
      integer(4)::lxx,nxxx,ir,n,l,ierr,   lcutmx
      integer(4)::nsp,nclass,ifppb,isp,ip1,ip2
      real(8) :: phitoto(nrx,0:nl-1,nn,nclass, nsp)
      real(8) :: phitotr(nrx,0:nl-1,nn,nclass, nsp)
      real(8),allocatable :: ppbrd(:,:,:,:,:,:)
     &  ,wqt(:), q0i(:,:),absqg2x(:)
      real(8),allocatable :: phiav(:,:,:)
      real(8) :: rnormphi(0:nl-1,nn) ,alat,sss
      integer(4) :: irx,noo,nuu      ,nq0i, neps, nq0ix,ix
      integer(4):: zvztest,nlmlsp
      real(8):: ovvs,rrr(nrx)
C---for mode 8
      integer(4):: nr_r,nlml_r,nsp_r,ifv,ibas,ilmx,isp1,isp2,kxx,ixx,nxx_i
      real(8),allocatable:: rofi_r(:),rho1(:,:,:),rspin(:),den(:,:),r11(:)
      real(8):: spinvec,spinvec0,sumc(2),sqrtfpi,const
      character(3):: charnum3
      logical:: valmt=.false.
      logical :: smbasis,newaniso,addbasnew
      character*15 rprodf
      integer(4) :: ifprodmt, smbasiscut,l1l2p,l1l2m,inn
      integer(4) :: lxlnln(0:2*(nl-1), 0:nl-1,nn,0:nl-1,nn) ,nzz
      real(8),allocatable::
     &   prodmt (:,:,:)
      real(8):: prr(3),derie,derie2,ddd,derie3
      integer(4):: nc_max(0:nl-1),nnn, smbasis_case
      character(len=100):: recxxx
      character(len=160):: recxxx2
      integer(4):: npbasmax(0:2*(nl-1)),ifinin,iax,izz,naxx,verbose
      integer(4),allocatable:: ipx(:,:)
      real(8):: bb1,bb1s,aa_in,bb_in
      integer(4):: nl2m1,nrofi_in,ifq0p,ifile_handle
c-------------------------------------------
      print *,' basnfp_v2: ********** start ******** nrofi=',nrofi
c      if(iread==1) goto 2001
      allocate(phiav(nrx,0:nl-1,nn),rprod(1:nrx, npradmx),ipx(0:2*(nl-1),nxxmx))
! To make the sign of phitotr safer; but I may already added something for keeping sss=1d0, right?
      do l1 = 0, nl-1
        do n1 = 1, nindx(l1)
          if(nsp==1) then
            phiav(1:nrx, l1, n1) = phitotr(1:nrx, l1, n1, ic, 1)
          else
            sss = dsign(1d0, phitotr(2,l1,n1,ic, 1)*phitotr(2,l1,n1,ic, 2))
            phiav(1:nrx, l1, n1) = 
     &      ( phitotr(1:nrx, l1, n1, ic, 1) + sss* phitotr(1:nrx, l1, n1, ic, 2) )/2d0
          endif
        enddo
      enddo
c$$$! ptest. Better mesh points for ptest because of its behevior near r=rmax.
c$$$      if(.false.) then
c$$$        aa = 0.2d0*aa
c$$$        bb = r(nrofi)/( exp(aa*(nrofi-1))-1d0)
c$$$        r(1)=0d0
c$$$        do i =2,nrofi
c$$$          r(i)=bb*( exp(aa*(i-1))-1d0)
c$$$        enddo
c$$$      endif
      if(ixx==4) goto 1212
!-norm check
      do l1 = 0, nl-1
      do n1 = 1, nindx(l1)
          call gintxx(phiav(1,l1,n1),phiav(1,l1,n1),aa,bb,nrofi,ovv)
          write(6,"(' norm check for phi='2i3,d13.6)") l1,n1,ovv
          if(abs(ovv) <1d-10 ) ovv  = 0d0
          rnormphi(l1,n1) = ovv
      enddo
      enddo
!! product basis construction
      lxlnln=0
      npr    = 0
      nprpd  = 0
      nxx    = 0
      nxxold = 0
      iprad  = 0
      rprod  = 0d0
      do 203 l1 = 0, nl-1
      do 202 n1 = 1, nindx(l1)  ; noo = nocc  (l1,n1)
      do 201 l2 = 0, nl-1
      do 20 n2 = 1, nindx(l2)  ; nuu = nunocc(l2,n2)
!!  phiav * phiav ----
        if( noo>=1 .and. nuu>=1 ) then
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
            rprod(2:nrofi,iprad)=phiav(2:,l1,n1)*phiav(2:,l2,n2)/r(2:) ! phi = u = r \phi
c          call gintxx(phiav(1,l1,n1), phiav(1,l2,n2),aa,bb,nrofi, sss )
c          write(6,"(' normchk phiav*phiav=',4i3,d13.6)") l1,n1,l2,n2,sss
          endif
        endif
   20 continue
 201  continue
 202  continue
 203  continue
      
!! Additional product bais for smbais()=T.
      if(smbasis()) then ! Dec2005
c$$$!!--- CASE1. Add r^l r^(l+2) anyway. ---
        if(smbasis_case()==1) then
          do 110 lx  = 0, 2*(nl-1)
            if(lx > smbasiscut() ) cycle
            do inn = 1, 2
              iprad = iprad + 1
              rprod(1,iprad) = 0d0
              if(inn==2) nnn = lx + 1
              if(inn==1) nnn = lx + 3
              rprod(2:nrofi,iprad)= r(2:)**nnn
              nxx(lx) = nxx(lx) + 1
              ipx(lx, nxx(lx)) = iprad
              write(6,"('sm -- lx iprad nnn nxx= ',6i5)")lx, iprad,nnn,nxx(lx)
            enddo
 110      continue
        elseif(smbasis_case()==2) then
!!--- CASE2 ---
          do 111 lx  = 0, 2*(nl-1)
            if(lx > smbasiscut() ) cycle
            nxx_i =nxx(lx)
            do inn = nxx_i+1, 2
              iprad = iprad + 1
              rprod(1,iprad) = 0d0
              if(inn==2) nnn = lx + 1
              if(inn==1) nnn = lx + 3
              rprod(2:nrofi,iprad)= r(2:)**nnn
              nxx(lx) = nxx(lx) + 1
              ipx(lx, nxx(lx)) = iprad
              write(6,"('sm -- lx iprad nxx= ',6i5)")lx, iprad,nxx(lx)
            enddo
 111      continue
        elseif(smbasis_case()==3) then
!!--- CASE3 ---
!! GaAs local-orbital tests suggest the above choice looks better for Ga 3d core.
!!  High priority to low priority
!! Priority ordering
!!     1. smaller l_1 +l_2,
!!     2. smaller |l1-l2|
!!     3.  phi*phi,  phi*phidot, phi*local
          do 112 lx  = 0, 2*(nl-1)
            if(lx>  smbasiscut() ) cycle
            inn=0
            l1= lx/2
            l2=  lx-l1   !l1<=l2
!       Priority ordering is phi*phi,  phi*phidot  !, phi*local
            do 120 n1 = nc_max(l1) +1,nc_max(l1) +1  !nindx(l1)
              do 130 n2 = nc_max(l2) +1, nindx(l2)
                inn = inn+1
                if( lxlnln(lx,l1,n1,l2,n2)==0) then
                  if( npr(l1,n1,l2,n2) == 0 ) then
                    iprad = iprad + 1
                    npr (l1,n1,l2,n2) = iprad
                    npr (l2,n2,l1,n1) = iprad
                    rprod(1,iprad) = 0d0
                    rprod(2:nrofi,iprad)=phiav(2:,l1,n1)*phiav(2:,l2,n2)/r(2:) ! phi = u = r \phi
                  else
                    iprad = npr(l1,n1,l2,n2)
                  endif
                  nxx(lx) = nxx(lx) + 1
                  ipx(lx, nxx(lx)) = iprad
                  lxlnln(lx,l1,n1,l2,n2)= iprad
                  lxlnln(lx,l2,n2,l1,n1)= iprad
                  write(6,"('sm -- lx l1 l2 n1 n2 iprad= ',6i5)")lx, l1, l2,n1, n2,iprad
                endif
                if(inn==2) goto 121
 130          continue
 120        continue
 121        continue
 112      continue
        endif
c      print *, '======= basis set up =============='
c      do lx  = 0, 2*(nl-1)
c      do nx = 1, nxx(lx) ;  ib1 = ipx(lx,nx)
c        nprod = nprod + 2*lx+1
c        write(6,"('lx nx number ipx=',4i4)") lx, nx, 2*lx+1, ipx(lx,nx)
c      enddo
c      enddo
c      print *, '--- test end ---------'
c      return
      endif
!! === Add a PW like product basis  ===
c Here you can add any functions in the same way.
      if(smbasis().and.smbasis_case()==1) goto 113
!!      
      if(.not.addbasnew().and.zvztest()/=2) then !Add one s product basis.
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
c ---
 113  continue
      nprad = iprad
      print *,' number of radial basis: nprad=', nprad


ccccccccccc start of ptest cccccccccccccccccccccccccccccccccccccccccccc
 1212 continue
      if(ixx==4) then
        print *, ' *** TEST case ***  rprodx is given by Bessel.'
        lxx   = 2*(nl-1)
c      lxx  = 4
        allocate(phij(0:lxx),psij(0:lxx))
c q near zero test
c This section is just to set nxxx and absqg2x(1:nxxx).
        write(6,*) '--- readin Q0P -------'
        ifq0p=ifile_handle()
        open (ifq0p,file='Q0P')
        read (ifq0p,"(i5)") nq0i
        allocate( wqt(1:nq0i),q0i(1:3,1:nq0i),absqg2x(nq0i) )
c     nq0ix = nq0i
        nzz = max(2,nq0i)
        if( smbasis() ) then
          write(6,"(' smbasis case' )")
          ix = 0
          do i=1,nzz !nq0i
            ix=ix+1
            if(i<= nq0i) then
              read (ifq0p, * ) wqt(i),q0i(1:3,i)
              absqg2x(ix) =sum( (2*pi/alat *q0i(1:3,i))**2) !bug 30jan2005 ---q0i(q:3,nq0i)
            else
              absqg2x(ix) = absqg2x(1)*i*2
            endif
            if(ix>1) call addd(absqg2x(ix),absqg2x,ix-1)
          enddo
ccccccccccccccccccccccc
c       absqg2x(1) =(2*pi/alat*4.18871)**2
c       absqg2x(2) =(2*pi/alat*4.33184)**2
cccccccccccccccccccccc
          do i=1,nzz !nq0i
            write(6,"('smbasis case i absqg =',i5,f13.5 )")
     &            i, sqrt(absqg2x(i)/(2*pi/alat)**2)
          enddo
        else
          ix=0
          do i=1,nq0i
            read (ifq0p, * ) wqt(i),q0i(1:3,i)
            if(wqt(i)==0d0 ) then
              ix=ix+1
              absqg2x(ix) =sum( (2*pi/alat *q0i(1:3,i))**2) !nq0i ---> q0i
              if(ix>1) call addd(absqg2x(ix),absqg2x,ix-1)
            endif
          enddo

        endif
c      neps = nq0i - nq0ix
        nxxx = ix
        deallocate(wqt,q0i)
        close(ifq0p)
c      write(6,*)"----- read end of Q0P -----"
ccccccccccccccccccccccc q near zero test end
c     nxxx= 4
c      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i),absqg2x(nq0i) )
c      absqg2x(1) = 3.66707D+00
c      absqg2x(2) = 2.94331D+00
c      absqg2x(3) = 0.628174D+01
c      absqg2x(4) = 0.703180D+01
ccccccccccccccccccccccccccccccccccccccccccccccccccc
        iprad = 0
        nxx = 0

        do lx = 0, lxx
          if(lx >lcutmx) cycle  ! Lmax cutoff for product basis
          do n1 = 1, nxxx
            if(smbasis().and.lx>  smbasiscut() ) cycle
            iprad  = iprad + 1
            nxx(lx)= nxx(lx)+1
            ipx(lx, nxx(lx)) = iprad
            absqg2 = absqg2x(n1)
            do ir =1,nrofi
              call bessl2x(absqg2*r(ir)**2,max(2,lx),phij,psij) !second argument must be larger than 1.
              rprod(ir,iprad) = phij(lx)* r(ir) **(lx +1)
            enddo
            print *,' sumchk rprod=',lx,n1,sum(abs(rprod(1:nrofi,iprad)))
cccccccccccccccccccccc
c        if(lx==8) then
            if(verbose()>60) then
              write(3100+ic,"(' -- -- -- ',3i3,' --- ' )") lx,n1
              do ir =1,nrofi
                write(3100+ic,"(d13.5,2x,2d18.8)")   
     &     r(ir), rprod(ir,iprad)
              enddo
            endif
c        endif
cccccccccccccccccccc
          enddo
        enddo
        nprad=iprad
        print *, ' *** TEST nprad=',nprad
        print *, ' nxx =',nxx(0:lxx)
      endif
ccccccccccc end of ptest ccccccccccccccccccccccccccccccccccccccccccccccc

      
! sanity check
      if(maxval(nxx(0:2*(nl-1)))>nxxmx) call rx(' basnfp: nxx >nxxmx --- Enlarge nxxmx in basnfp.f')
      if(nprad > npradmx)call rx( ' basnfp: nprad > npradmx --- Enlarge npradmx in basnfp.f')

!! nprod
      nprod = 0
      do lx = 0, 2*(nl-1)
        do nx = 1, nxx(lx)
          nprod = nprod + 2*lx+1
          write(6,"('lx nx number ipx=',4i4)") lx, nx, 2*lx+1, ipx(lx,nx)
        enddo
      enddo
      print *,' *** total number of product basis nprod=', nprod
!! ovmt
      kmax = 2*(nl-1)
      allocate( ovmt  (nprad,nprad,0:2*(nl-1))  )
      do lx  = 0, 2*(nl-1)
        do nx1 = 1, nxx(lx) ;  ib1 = ipx(lx,nx1)
          do nx2 = 1, nxx(lx) ;  ib2 = ipx(lx,nx2)
            call gintxx(rprod(1,ib1),rprod(1,ib2),aa,bb,nrofi, ovv)
            ovmt(ib1,ib2, lx) = ovv !radial overlap matrix of rprod 
          enddo
        enddo
      enddo
!! orthonormal basis functions.--------------------
      allocate(rprodx(nrofi,maxval(nxx(0:2*(nl-1))),0:2*(nl-1)))
      do lx = 0, 2*(nl-1)
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
        call rs(nb, ovvc(:,:,1), eb, zz(:,:,1),ierr)
        if(ierr/=0) call rx( ' basnfp: rs error')
c$$$! PMASMAX option. rarely used now.
c$$$        npbasmax=100
c$$$        if(ixx==0.or.ixx==8) then
c$$$          call getkeyvalue("GWinput","<PBASMAX>", unit=ifinin,status=naxx,errstop='off')
c$$$          if(naxx<0) goto 1250
c$$$          do izz = 1, naxx
c$$$            read(ifinin,"(a100)",end=1100) recxxx
c$$$            recxxx2=recxxx//
c$$$     &      " 100 100 100 100 100 100 100 100 100 100 100 100 100 100"
c$$$            read(recxxx2,*) iax, npbasmax(0:2*(nl-1))
c$$$            if(iax==ic) then
c$$$              write(6,"('<PBASMAX> gives l npbas=',2i3)") lx, npbasmax(lx)
c$$$              goto 1200
c$$$            endif
c$$$          enddo
c$$$ 1100     continue
c$$$          npbasmax=100
c$$$ 1200     continue
c$$$          close(ifinin)
c$$$ 1250     continue  
c$$$        endif
        ibx=0
        do ib1 = nb,1,-1
          write(6,"(a,i5,d13.6,a,d13.6)")'    ib eb=',ib1,eb(ib1),'  ecut=',cutbase(lx)
c          if(npbasmax(lx)/=100) then 
c            if(ibx+1>npbasmax(lx)) then
c              if(smbasis()) then
c                if(ibx+1>2) cycle
c              else
c                cycle
c              endif
c            endif
c          else
          if(smbasis().and.ib1>=nb-1) then
              continue
          elseif(eb(ib1)<cutbase(lx)) then !this is mainly used now 2022jan
              cycle
          endif
c          endif
          ibx = ibx+1
          ibo(ibx) = ib1
          rprodx(1:nrofi,ibx,lx) = matmul( rprod(1:nrofi,ipx(lx,1:nb)), zz(1:nb,ib1,1) )
          rprodx(1:nrofi,ibx,lx) = rprodx(1:nrofi,ibx,lx)/sqrt(eb(ib1))
        enddo
        nb     =ibx
        nxx(lx)=ibx
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c debug        
c      if(lx==2) then
c        ibx =1
c        print *,' ibx=', ibx, eb(ibx)
c        do ib1 =1,nb
c          call gintxx(rprod(1,ipx(lx,ib1)),rprodx(1,ibx,lx)
c     &              ,aa,bb,nrofi,ovv)
c          ovv = ovv*sqrt(eb(ibx))
c          write(6,"(' zz xxx',4d16.7)")
c     &    zz(ib1,ibx,1), ovv,   !zz is normalised to 1.
c     &    ovv/ ( zz(ib1,ibx,1)*eb(ibx) )
c        enddo
c        print *,' norm=',sum(zz(1:nb,ibx,1)**2)
c      stop
c      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do ib1=1,nb
        do ib2=1,nb
            call gintxx(rprodx(1,ib1,lx),rprodx(1,ib2,lx),aa,bb,nrofi,ovv)
            if(ib1==ib2 ) then
              write(6,"('    Diag ibx ovv=',i3,d23.16,' eb= ',d18.10,' nod= ',i3)")
     &              ib1,ovv, eb(ibo(ib1)), nodnum(rprodx(1,ib1,lx),nrofi)
            elseif( ovv>1d-10 ) then
              write(6,"('      offdiag ib1 ib2 ovv=',2i3,d23.16,' eb=',d18.10)")
     &              ib1,ib2,ovv, eb(ibo(ib1))
            endif
        enddo
        enddo
        write(6,"(' *** lx  *** Used nb =',2i5)") lx,nb
        deallocate( ovvc, zz, eb,ibo)
      enddo

cccccccccccccc test2 ccccccccccccccccccccccccccccccccccccccccccc
      if(.false.) then
        allocate(rprodx2(nrx,nxxx,0:lxx))
        do n  =1,nxxx
          if(n==1) absqg2 = 2.524974**2
c       if(n==2) absqg2 = 2.598177**2
          if(n==2) absqg2 =  .612396**2
          do ir =1,nrofi
            call bessl2x(absqg2*r(ir)**2, lxx, phij, psij)
            do l = 0, lxx
              rprodx2(ir,n,l) = phij(l)* r(ir) **(l +1 )
            enddo
          enddo
        enddo
        do l   = 0, lxx
c       rprodx2(1:nrofi,1,l)= rprodx2(1:nrofi,1,l)  + rprodx2(1:nrofi,2,l)
          n = 1
          call gintxx(rprodx2(1,n,l),rprodx2(1,n,l)
     &    ,aa,bb,nrofi, aaa )
          aaa = 1d0/sqrt(aaa)
          rprodx2(1:nrofi,n,l)= aaa*rprodx2(1:nrofi,n,l)
          if(nxxx==1) cycle
          n1=1
          n2=2
          call gintxx(rprodx2(1,n1,l),rprodx2(1,n2,l)
     &    ,aa,bb,nrofi, aaa12 )
          rprodx2(1:nrofi,n2,l) = rprodx2(1:nrofi,n2,l)
     &      - aaa12*rprodx2(1:nrofi,n1,l)
          n = 2
          call gintxx(rprodx2(1,n,l),rprodx2(1,n,l)
     &    ,aa,bb,nrofi, aaa )
          aaa = 1d0/sqrt(aaa)
          rprodx2(1:nrofi,n,l)= aaa*rprodx2(1:nrofi,n,l)
        enddo
      endif
ccccccccccccccend of test2 ccccccccccccccccccccccccccccccccccccccccccc
      
C- Reserve rprodx
C      real(8):: crbase(nrx,kmxx,nclass)
C      integer(4):: kmxx,kmx(nclass),iprlc(kmxx,nclass),lprc(kmxx,nclass)
c----------------------------------------------
      kmx = sum( nxx(0:2*(nl-1)) )
      allocate( iprlc(kmx), lprc(kmx) )
      k = 0
      iprlc(1:kmx) = 0
      do lx = 0, 2*(nl-1)
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
!------------------------
      print *,' reserve rprodx...'
      filename = 'BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10) )
      open(newunit=ificrb,file=trim(filename))
      write(ificrb,"(4i6,2d24.16)") 2*(nl-1), kmx, nblocha,nrofi,aa,bb
      write(ificrb,"(i5)") nxx(0:2*(nl-1))
      k = 0
      do lx = 0, 2*(nl-1)
        do nx = 1, nxx(lx)
          k = k + 1
          write(ificrb,"(3i5)"   ) k,iprlc(k),lprc(k)
          write(ificrb,"(d23.15)") (rprodx(i,nx,lx),i=1,nrofi)
        enddo
      enddo
      close(ificrb)
      deallocate(rprodx,iprlc,lprc)
c 2001 continue
!  Write and Read even for iread==0,
!  in order to have the exact match on PPB* when iread=1 for the same BASFP* file.
      print *,' read rprodx...'
      filename = 'BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10) )
      open(newunit=ificrb,file=trim(filename))
      read(ificrb,"(4i6,2d24.16)") nl2m1, kmx,nblocha,nrofi_in,aa_in,bb_in
      if( nl /= nl2m1/2 +1)   call rx( 'wrong 2*(nl-1) in readin BASNFP')
      if(nrofi_in /= nrofi)   call rx( 'nrofi_in/=nrofi in readin BASNFP')
      if(abs(aa-aa_in)>1d-12) call rx( 'aa_in/=aa in readin BASNFP')
      if(abs(bb-bb_in)>1d-12) call rx( 'bb_in/=bb in readin BASNFP')
      write(6,"(' basnfp: BASFP... kmx nblocha=',2i5)") kmx,nblocha
      read(ificrb,"(i5)") nxx(0:2*(nl-1))
      allocate(rprodx(nrofi,maxval(nxx(0:2*(nl-1))),0:2*(nl-1)))
      allocate( iprlc(kmx), lprc(kmx) )
      k = 0
      do lx = 0, 2*(nl-1)
        do nx = 1, nxx(lx)
          k = k + 1
          read(ificrb,"(3i5)"   ) k,iprlc(k),lprc(k)
          read(ificrb,"(d23.15)") (rprodx(i,nx,lx),i=1,nrofi)
        enddo
      enddo
      close(ificrb)
c$$$      if(iread==1) then
c$$$        print *,' read rprodx...'
c$$$        filename = 'BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10) )
c$$$        ificrb   = iopen ( filename,0,-1,0)
c$$$        read(ificrb) nl2m1, kmx,nblocha,nrofi_in,aa_in,bb_in
c$$$        if( nl /= nl2m1/2 +1) stop 'wrong 2*(nl-1) in readin BASNFP'
c$$$        if(nrofi_in /= nrofi) stop 'nrofi_in/=nrofi in readin BASNFP'
c$$$        if(abs(aa-aa_in)>1d-12) stop 'aa_in/=aa in readin BASNFP'
c$$$        if(abs(bb-bb_in)>1d-12) stop 'bb_in/=bb in readin BASNFP'
c$$$        write(6,"(' basnfp: BASFP... kmx nblocha=',2i5)") kmx,nblocha
c$$$        read(ificrb) nxx(0:2*(nl-1))
c$$$        allocate(rprodx(nrofi,maxval(nxx(0:2*(nl-1))),0:2*(nl-1)))
c$$$        allocate( iprlc(kmx), lprc(kmx) )
c$$$        k = 0
c$$$        do lx = 0, 2*(nl-1)
c$$$        do nx = 1, nxx(lx)
c$$$        k = k + 1
c$$$        read(ificrb) k,iprlc(k),lprc(k)
c$$$        read(ificrb) (rprodx(i,nx,lx),i=1,nrofi)
c$$$        enddo
c$$$        enddo
c$$$        isx = iclose(filename)
c$$$      else
c$$$        print *,' reserve rprodx...'
c$$$        write(6,"(' basnfp: BASFP... kmx nblocha=',2i5)") kmx,nblocha
c$$$        filename = 'BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10) )
c$$$        ificrb   = iopen ( filename,0,-1,0)
c$$$        write(ificrb) 2*(nl-1), kmx, nblocha,nrofi,aa,bb
c$$$        write(ificrb) nxx(0:2*(nl-1))
c$$$        k = 0
c$$$        do lx = 0, 2*(nl-1)
c$$$        do nx = 1, nxx(lx)
c$$$          k = k + 1
c$$$          write(ificrb) k,iprlc(k),lprc(k)
c$$$          write(ificrb) (rprodx(i,nx,lx),i=1,nrofi)
c$$$        enddo
c$$$        enddo
c$$$        isx = iclose(filename)
c$$$      endif

! Calculate radial matrix elements.
      print *,' Calculate radial matrix elements...'
      allocate( ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1)
     &           ,maxval(nxx(0:2*(nl-1))) ) )
      ppbrd =.9999999999999d99 !for safe
      filenamep = 'PPBRD_V2_'//char( 48+ic/10 )//char(48+mod(ic,10))
      open(newunit=ifppb,file=trim(filenamep),form='unformatted')
      write(ifppb) nblocha, 2*(nl-1), nxx(0:2*(nl-1))
      do isp= 1,nsp
        isp1=isp
        isp2=isp
        if(ixx==8.and.isp==1) then
          isp1=1
          isp2=2
        elseif(ixx==8.and.isp==2) then
          isp1=2
          isp2=1
        endif
        do  lx = 0, 2*(nl-1)
        do  nx = 1, nxx(lx)
        do  l1 = 0, nl-1
        do  n1 = 1, nindx(l1)
        do  l2 = 0, nl-1
        do  n2 = 1, nindx(l2)
          if(lx <abs(l1-l2) .or. l1+l2<lx) cycle
          rphiphi(1)       = 0d0
          rphiphi(2:nrofi) = phitoto(2:nrofi,l1,n1,ic,isp2)
     &                      *phitoto(2:nrofi,l2,n2,ic,isp1)/r(2:) ! phi = u = r \phi
          call gintxx(rprodx(1,nx,lx), rphiphi,aa,bb,nrofi, ppbrd(l1,n1,l2,n2,lx,nx) )
        enddo  
        enddo  
        enddo  
        enddo  
        enddo  
        enddo  
        write(ifppb) ppbrd
      enddo
      deallocate(ppbrd)
      close(ifppb)
      
! prodmt proddmt (value and slope of the product at MT). !oct2005
!    Stored into the tail of PPBRD* files.
      if(smbasis()) then
        allocate( prodmt (2,maxval(nxx(0:2*(nl-1))),0:2*(nl-1)))
        prodmt = 1d10 !for safe
        do lx = 0, 2*(nl-1)
          do nx = 1, nxx(lx)
            prodmt(1,nx,lx)= rprodx(nrofi,nx,lx)/r(nrofi)
            prodmt(2,nx,lx)= derie2(r, rprodx(1:nrofi,nx,lx)/r(1:nrofi), nrofi)
          enddo
        enddo
        filenamep = 'PRODMT_'//charnum3(ic)
        open(newunit=ifprodmt,file=trim(filenamep),form='unformatted')
        write(ifprodmt) nl
        write(ifprodmt) maxval(nxx(0:2*(nl-1)))
        write(ifprodmt) nxx(0:2*(nl-1))
        write(ifprodmt) prodmt
        close(ifprodmt)
        deallocate(prodmt)
      endif
      
! --- MixSpin= <rho_up - rho_down | B> matrix calculation. May2005
!r1  Suppose  "ibas==iclass"--- it is already checked in hbasfp0.m.f
!r2  ValMT.* is written with subroutine savemtval(ib,rho1,rofi,nr,nlml,nsp)
!r    in fp/locpot.f just befor locpt2 in lmto (lmf).
      if(ixx==8) then
        sqrtfpi = sqrt(fpi)
        ifv = ifile_handle()
        ibas=ic
        if(valmt) then
          open(ifv,file='ValMT.'//charnum3(ibas)//'.chk',form='unformatted')
          read(ifv) nr_r,nlml_r,nsp_r
          write(6,"('readin nr nlml nsp=',3i5)") nr_r,nlml_r,nsp_r
          allocate(rofi_r(nr_r),rho1(nr_r,nlml_r,nsp_r),rspin(nrofi),den(nrofi,nsp),r11(nrofi))
          r11(1:nrofi)= 1d0
          read(ifv) rofi_r, rho1
          close(ifv)
        else
          open(ifv,file='rhoMT.'//xtxx(ibas),form='unformatted',status='old',err=1031)
          goto 1032
 1031     continue 
          call rx( 'rhoMT by locpot-wrhomt. open error')
 1032     continue
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
        open(ifv,file='MixSpin.'//charnum3(ibas))
        write(ifv,"(2i10,' ! ibas, max l of product basis' )") ibas,2*(nl-1)
        write(ifv,"(i10,'           ! nxx(lx)'  )") nxx(0:2*(nl-1))
        do ilmx = 1, (2*(nl-1)+1)**2
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
! const = <1|B> where 1 is normalized within the sphere 2007
            if(lx==0) then
              call gintxx( rprodx(1,nx,lx), r,aa,bb,nrofi, const )
c             call gintxx( r, r,aa,bb,nrofi, const )
            else
              const=0d0
            endif
            const= const *  sqrtfpi !/((fpi/3d0)*r(nrofi)**3)
! Now spinvec = <B_I(\bfr) | m_true(\bfr) >
            if(abs(spinvec)<1d-10 ) spinvec=0d0
            write(ifv,"(     2i6,d24.16,2x,f13.10,2x,f13.10,d24.16
     &       ' ! I=(ilm, nx), <spin|B(I)>, intrho(1:nsp) <1|B(I)>')")
     &       ilmx, nx, spinvec, sumc(1:nsp),const
            write(6,"('ttt:',2i6, d24.16,2x, 2f14.10,d24.16
     &       ' ! I=(ilm, nx), <spin|B(I)>, intrho(1:nsp) <1|B(I)>')")
     &       ilmx, nx, spinvec, sumc(1:nsp),const
          enddo
        enddo
        deallocate(rofi_r,rho1,rspin)
        close(ifv)
      endif !ixc==8 end
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
      end

C ----------------------------------------------------------------
      subroutine bessl2x(y,lmax,phi,psi)
cnote --- what the difference from nfpsrc/bessl? : I think essentially the same.
C- Radial part of Bessel functions
C ----------------------------------------------------------------
Ci Inputs
Ci   Y = E * R**2;  lmax
Co Outputs
Co   phi:  first (lmax+1) spherical bessel functions / r^l
Co         for e -> 0 returns phi(l) = 1/(2l+1)!!
Co   psi:  first (lmax+1) spherical hankel functions * r^(l+1)
Co         for e -> 0 returns psi(l) = (2l-1)!!
Cr Remarks
Crxxx   Andersen's definition in the limit E->0:
Crxxx   bessel phi(OKA) is phi * (2l-1)!!/2  and
Crxxx   hankel psi(OKA) is psi / (2l-1)!!, making
Crxxx   H(OKA) = r^-l-1 and  J(OKA) = r^l/(2(2l+1))
C ----------------------------------------------------------------
C     implicit none
      integer lmax
      double precision y
      double precision phi(lmax+1),psi(lmax+1)
      integer i,isn,j1,j2,k,l,lmux,lmuxp1,lmuxp2,lp1,nf,tlp1,tmp1,tmp2
      double precision dt,dt1,dt2,exppr,my,srmy,t,t1,tol
      double precision dum(420)
cC#ifdef OKA
cC A table of (2l-1)!!
c      integer fac2l(10)
c      data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
cC#endif
      real(8):: fac2l
      lmux = max0(lmax,2)
      if (lmux .gt. 9) call rx( 'bessl2: lmax gt 9')
      tol = 1.d-8
      my = -y
      l = lmux
    1 tlp1 = l+l+1
      i = 1
      do  2  k = 3, tlp1, 2
        i = i*k
    2 continue
      t1 = 1.d0/dble(i)
      dt = 1.d0
      t = 1.d0
      i = 0
      do  3  k = 1, 10000
        i = i+2
        dt1 = i
        dt2 = i+tlp1
        dt = dt*my/(dt1*dt2)
        t = t+dt
        if (dabs(dt) .lt. tol) goto 4
    3 continue
      goto 10
    4 if (l .lt. lmux) goto 5
      dum(1) = t1*t
      l = lmux-1
      goto 1
    5 dum(2) = t1*t
      tmp1 = lmux + lmux + 1
      tmp2 = tmp1 + 1
      nf = tmp1
      do  6  k = 3, tmp2
        nf = nf-2
        dum(k) = nf*dum(k-1) - y*dum(k-2)
    6 continue
      lmuxp1 = lmux+1
      lmuxp2 = lmux+2
      isn = -1
      do  7  k = 1, lmuxp1
        isn = -isn
        j1 = lmuxp2-k
        j2 = lmuxp1+k
        phi(k) = dum(j1)
        psi(k) = dum(j2)*isn
    7 continue
      if (y .ge. 0d0) goto 40
C ------- NEGATIVE ENERGY CASE ----------
      srmy = dsqrt(-y)
      psi(2) = 1.d0+srmy
      psi(1) = 1.d0
      if (lmux .lt. 2) goto 23
      tlp1 = 1
      do  21  lp1 = 3, lmuxp1
        tlp1 = tlp1+2
        psi(lp1) = tlp1*psi(lp1-1) - y*psi(lp1-2)
   21 continue
   23 exppr = 1.d0/dexp(srmy)
      do  22  lp1 = 1, lmuxp1
        psi(lp1) = psi(lp1)*exppr
   22 continue
C -------- EXIT --------
   40 continue
C#ifdef OKA
      do  42  lp1 = 1, lmuxp1
        phi(lp1) = (phi(lp1)*fac2l(lp1))/2
        psi(lp1) =  psi(lp1)/fac2l(lp1)
   42 continue
C#endif
      return
   10 write(*,11) y
   11 format(' BESSL2: power series not convergent, E*r**2=',e12.4)
      call rx( '')
      end

      real(8) function fac2l(i)
cC A table of (2l-1)!!
c     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
      integer:: i,l
      logical,save::  init=.true.
      real(8),save:: fac2lx(101)
      if(init) then
        fac2lx(1)=1d0
        do l=1,100
          fac2lx(l+1)=fac2lx(l)*(2*l-1)
        enddo
      endif
      fac2l=fac2lx(i)
      end

      real(8) function derie (x,y)
      implicit none
      real(8) :: x(3),y(3),deri1,deri2,xm1,xm2,xx,deriei,dxdi,dydi
c      xm1  = (x(2)+x(1))/2d0
c      deri1= (y(2)-y(1))/(x(2)-x(1))
c      xm2  = (x(3)+x(2))/2d0
c      deri2= (y(3)-y(2))/(x(3)-x(2))
c      xx = x(3)
c      deriei = deri1 + (deri2-deri1)/(xm2-xm1) *(xx - xm1)
cc dxdi at end
cc      dxdi1 = x(2) - x(1)
cc      dxdi2 = x(3) - x(2)
      dxdi = x(3) - x(2) + .5d0*(x(3)- 2*x(2) +x(1))
      dydi = y(3) - y(2) + .5d0*(y(3)- 2*y(2) +y(1))
      derie = dydi/dxdi
      end

      real(8) function derie3 (x,y)
      implicit none
      real(8) :: x(3),y(3),deri1,deri2,xm1,xm2,xx,deriei,dxdi,dydi
      xm1  = (x(2)+x(1))/2d0
      deri1= (y(2)-y(1))/(x(2)-x(1))
      xm2  = (x(3)+x(2))/2d0
      deri2= (y(3)-y(2))/(x(3)-x(2))
      xx = x(3)
      derie3 = deri1 + (deri2-deri1)/(xm2-xm1) *(xx - xm1)
      end

      real(8) function derie2 (x,y,n)
c     return derivative at n
      integer:: n,i
      real(8):: x(n), y(n), dxdi(n),dydi(n),polinta,ii(n),dydx(n)
      do i=1,n-1
        dydx(i) = (y(i+1)-y(i))/(x(i+1)-x(i))
        ii(i) = i + 0.5d0
      enddo
      derie2 = polinta(dble(n), ii,dydx,n-1)
      end

      subroutine addd(a1,a,n)
! a1 is choosed so that it is not in agreement with other a(n).
      implicit none
      integer(4):: n,i
      real(8):: eps=0.1,a1,a(n)
 880  continue
      do i=1,n
        if(abs(a1-a(i))< 0.2) then
          a1=a1+eps
          goto 880
        endif
      enddo
      end

c$$$      character(3) function charnum3n(num)
c$$$      integer(4) ::num
c$$$      charnum3n=''
c$$$      charnum3n = char(48+mod(num,10))
c$$$      if(num>9)  charnum3n=char(48+mod(num/10,10))//charnum3n
c$$$      if(num>99) charnum3n=char(48+mod(num/100,10))//charnum3n
c$$$      end
c$$$
