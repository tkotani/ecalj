program hvccfp0
  !- Coulomb matrix. <f_i | v| f_j>_q. ------------------------
  ! output
  !    VCCFP : the coulomb matrix vcoul(nblochpmx,nblochpmx) for all qibz.
  !-------------------------------------------------------------
  ! int
  !    strx: structure constant for e=0 (means 1/|r-r'| )

  use m_genallcf_v3,only: Genallcf_v3,  alat,nbas=>natom, bas=>pos
  use m_hamindex,only:    Readhamindex,plat,qlat
  use m_read_bzdata,only: Read_bzdata, ginv,nqbz,qbz,nqibz,qibz,nq0i,wqt=>wt,q0i,nq0iadd
  use m_readqg,only:   readqg,readngmx
  use m_mpi,only: MPI__hx0fp0_rankdivider2,mpi__task,MPI__Initialize,MPI__Finalize,mpi__root, &
       MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,mpi__rank,mpi__size, &
       mpi__ranktab,MPI__consoleout,mpi__iend,mpi__iini,mpi__getrange
  use m_readgwinput,only: ReadGwinputKeys, &
       keeppositivecou
  use m_lgunit,only: m_lgunit_init
  
  implicit none
  integer :: ifvcfpout,ifhvccfp,is,  lmxcg,if1011,if3011, ifplane,ngpmx, ngcmx, nblochpmx, nbloch,&
       ibas,ic,lxx,nxx,nrx,l,n,k,isx,kdummy, nkdmx,nkqmx,lmax,nkdest,nkrest,ngp,ngc,nlxx,i,lnjcg,lnxcg, &
       nkd,nkq ,ibas1,ibas2,nlx1,nlx2, iqibz 
  real(8) ::  q(3),p(3),voltot, pi,fpi,tripl,alat0,epsx, tol,as,tpiba,qb0(3,3),vol0,rdist0,qdist0,radd,qadd, &
       a0,awald,alat1,tol1,r0,q0,awald0,qg(3),   absqg2,aaa,aaa12
  integer,allocatable :: jcg(:),indxcg(:), lx(:),kmx(:),nblocha(:),nr(:),ificrb(:), &
       nx(:,:),ngvecp(:,:),ngvecc(:,:),ngvecci(:,:,:),iqibzx(:)
  real(8),allocatable :: rmax(:), cg(:),rprodx(:,:,:,:),dlv(:,:),qlv(:,:),work(:),ngcn(:), &
       rojb(:,:,:), sgbb(:,:,:,:),aa(:),bb(:),rofit(:),phi(:),psi(:)
  complex(8) ,allocatable :: vcoul(:,:),geig(:,:),strx(:,:,:,:), &
       sgpb(:,:,:,:),sgpp(:,:,:,:), fouvb(:,:,:,:),fouvp(:,:,:,:),vcoul0(:,:), &
       s(:,:),sd(:,:),rojp(:,:,:) , vcoulnn(:,:)
  character*7,allocatable :: filename(:)
  character(20) :: xxt
  complex(8):: phasep,img=(0d0,1d0)
  integer(4)::ir,ig1,n1,n2
  complex(8),allocatable :: hh(:,:),oox(:,:),ooxi(:,:),oo(:,:),zz(:,:),zzr(:)
  real(8),allocatable    :: eb(:)
  complex(8),allocatable :: matp(:),matp2(:)
  complex(8) :: xxx,trwv
  integer(4) :: ngb,nev,nmx,iqx,ipl1,ipl2,igx1,igx2
  logical :: checkeig
  logical:: besseltest=.false. !test
  real(8) :: sss1,sss2,dnorm
  complex(8),allocatable:: gbvec(:), ppovl(:,:), b0mat(:)
  integer(4) ::igc,igc0,ifgb0vec,ifgb0vec1,ix, iy
  integer(4) :: iqxini, iqxend,imode
  logical :: allochk=.true.
  complex(8),allocatable:: hh1(:,:),oo1(:,:)
  integer(4)::  nqnumc,ifiqgc 
  real(8):: qqq(3),QpGcut_Cou       !,qq(3)
  integer(4),allocatable:: ngvecc0(:,:)
  integer(4):: ngc0
  real(8):: quu(3) !ginv(3,3),
  real(8),allocatable :: rkpr(:,:,:),rkmr(:,:,:),rofi(:,:)
  real(8):: eee,eees, q_org(3),screenfac
  integer(4):: ifvcfporg,nqbz_in,nblochpmx_in
  complex(8),allocatable:: vcoul_org(:,:)
  logical :: smbasis,debug=.false.,smbb
  integer(4)   :: ifprodmt,nl_r,lx_,nxx_r,nxdim,ibl1,nn,no,ngbnew, &
       nmatch,ifpmatch,nmatch_q,ifpmatch_q,m,ifpomat,nbln,ibln,ngb_in,nnr,igc2
  character(3) :: charnum3
  character(10) :: i2char
  character(11):: filenamep
  integer(4),allocatable:: nx_r(:), ibl(:,:,:,:),imatcho(:),imatchn(:),imatcho_q(:),imatchn_q(:)
  real(8),allocatable:: prodmt(:,:,:,:),rdmatch(:,:,:,:)
  complex(8),allocatable:: ppmt(:,:,:,:),pmat(:,:),pomat(:,:),oon(:,:)
  complex(8):: pval,pslo,phasex
  real(8)::absqq,qqx(3), epsmx,aaaa
  integer(4):: nnmx ,ngcnn,ngbo
  integer(4):: ifgb0vec_a,ifgb0vec_b , ifvcoud,idummy
  logical:: is_mix0vec,wvcc !,newaniso
  character(128):: vcoudfile
  real(8),allocatable:: wqfac(:),qbzwww(:,:)
  integer:: ifiwqfac,iqbz,iqbzx,nnn,ixyz,ifq0p,ifile_handle,incwfin
  character(128) :: ixcc
  logical:: cmdopt2
  character(20):: outs=''
  !!------------------------------------------------------------
  call MPI__Initialize()
  call M_lgunit_init()
  pi  = 4d0*datan(1d0)
  fpi = 4d0*pi
  if( mpi__root) write(6,"(' mode=0,3,202 (0 and 3 give the same results for given bas)' )")
  if(cmdopt2('--job=',outs)) then
     read(outs,*) imode
  elseif( mpi__root ) then
     read(5,*) imode
  endif
  call MPI__Broadcast(imode)
  write(ixcc,"('.mode=',i4.4)")imode
  call MPI__consoleout('hvccfp0'//trim(ixcc))
  call cputid (0)
  if(imode==202 ) then
     write(6,*)' hvccfp0: imode=',imode
  elseif(imode==0) then
  elseif(imode==3) then
  else
     call rx( 'hvccfp0: now hvccfp0 support just normal mode=0 3 202 101')
  endif
  !! jun2020, We read data not thru HVCCIN
  incwfin= 0                !use ForX0 for core in GWIN
  call Genallcf_v3(incwfin) !in module m_genallcf_v3
  call Readhamindex()
  call Read_bzdata()
  call readngmx('QGcou',ngcmx)
  allocate(ngvecc(3,ngcmx))
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  write(6,*)' voltot=',voltot
  call ReadGWinputKeys() !Readin dataset in GWinput

  !! Readin BASFP//atom. The product basis functions.
  if(allochk) write(*,*)'allocte(lx,kmx,nblocha,nr,aa,bb,filename,ificrb'
  allocate(lx(nbas),kmx(nbas),nblocha(nbas), &
       nr(nbas),aa(nbas),bb(nbas),filename(nbas), &
       ificrb(nbas),rmax(nbas) )
  do ibas = 1,nbas
     ic = ibas !
     filename(ibas)= 'BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10))
     open(newunit=ificrb(ibas),file=trim(filename(ibas)))
     read(ificrb(ibas),"(4i6,2d24.16)") &
          lx(ibas), kmx(ibas), nblocha(ibas), nr(ibas),aa(ibas),bb(ibas)
     !! rmax here jun2020
     rmax(ibas) = bb(ibas)*(exp((nr(ibas)-1)*aa(ibas))-1d0)
  enddo
  lxx = maxval(lx)
  if(allochk) write(*,*) 'allocate( nx(0:lxx,nbas) )'
  allocate( nx(0:lxx,nbas) )
  do ibas = 1,nbas
     read(ificrb(ibas),"(i5)") nx(0:lx(ibas),ibas)
  enddo
  nxx = maxval(nx)
  nrx = maxval(nr)
  if(allochk) write(*,*) 'allocate( rprodx(nrx,nxx,0:lxx,nbas) )'
  allocate( rprodx(nrx,nxx,0:lxx,nbas) )
  do ibas = 1,nbas
     do l = 0, lx(ibas)
        do n = 1, nx(l,ibas)
           read(ificrb(ibas),"(3i5)"   ) k, kdummy,kdummy
           read(ificrb(ibas),"(d23.15)") (rprodx(i,n,l,ibas),i=1,nr(ibas))
        enddo
     enddo
  enddo
  do ibas = 1,nbas
     close(ificrb(ibas))
  enddo


  !! Bessel TEST block -------------------------------------
  if(besseltest) then
     write(6,*)
     write(6,*)
     write(6,*) ' *** TEST case ***  rprodx is given by Bessel.'
     ! c test G, corresponding <q+G|v|q+G> should be exact. e.g. ig1=1 and  ig1=35 for iqx=2
     ! c You can change these values for tests. ccccccccccccc
     iqx  = 2
     igx1 = 1
     igx2 = 35
     write(6,"(' iqx=',i3,' ig1 ig2=',2i3)") iqx,igx1,igx2
     write(6,"(a)") &
          ' <q+G|v|q+G> for the corresponding iqx ig1 ig2 should be exact!'
     write(6,"(a)") ' See fort.196'
     write(6,"(a)") &
          ' Errors will be from the radial function integrals !!!'
     write(6,"(a)") ' You can slso so similar test from hbasfp0.'
     write(6,"(a)") ' See test1 in basnfp0.'
     if(allochk) write(*,*) 'deallocate(rprodx,nx)'
     deallocate(rprodx,nx)
     tpiba=8.d0*datan(1.d0)/alat
     lx = 4
     nr = nr(1)
     aa = aa(1)
     bb = bb(1)
     lxx = maxval(lx)
     if(allochk) write(*,*)'allocate( nx(0:lxx,nbas) )'
     allocate( nx(0:lxx,nbas) )
     kmx= 1
     nx = 2
     nxx = maxval(nx)
     nblocha= nxx *(lxx+1)**2
     nrx = maxval(nr)
     if(allochk) write(*,*)'allocate(rprodx,rofi ,phi,psi) '
     allocate(rprodx(nrx,nxx,0:lxx,nbas),rofit(nrx) &
          ,phi(0:lxx),psi(0:lxx))
     rofit(1) = 0d0
     do ir   = 1, nrx
        rofit(ir) = bb(1)*( exp(aa(1)*(ir-1)) - 1d0)
     enddo
     do n = 1, nxx
        if(n==1) ig1 = igx1
        if(n==2) ig1 = igx2
        qg(1:3) = &
             tpiba * (qibz(1:3,iqx)+ matmul(qlat, ngvecci(1:3,ig1,iqx)))
        absqg2  = sum(qg(1:3)**2)
        do ir =1,nrx
           call besslggg(absqg2*rofit(ir)**2,lxx,phi,psi)
           do ibas=1,nbas
              do l = 0, lx(ibas)
                 rprodx(ir,n,l,ibas) = phi(l)* rofit(ir) **(l +1 )
              enddo
           enddo
        enddo
     enddo
     ! --- orthogonalized rprodx.
     do ibas=1,nbas
        do l = 0, lx(ibas)
           rprodx(1:nr(ibas),1,l,ibas)= &
                rprodx(1:nr(ibas),1,l,ibas) &
                + rprodx(1:nr(ibas),2,l,ibas)
           n = 1
           call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas) &
                ,aa(ibas),bb(ibas),nr(ibas), aaa )
           aaa = 1d0/sqrt(aaa)
           rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
           if(nxx==1) cycle
           n1=1
           n2=2
           call gintxx(rprodx(1,n1,l,ibas),rprodx(1,n2,l,ibas) &
                ,aa(ibas),bb(ibas),nr(ibas), aaa12 )
           rprodx(1:nr(ibas),n2,l,ibas) = rprodx(1:nr(ibas),n2,l,ibas) &
                - aaa12*rprodx(1:nr(ibas),n1,l,ibas)
           n = 2
           call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas) &
                ,aa(ibas),bb(ibas),nr(ibas), aaa )
           aaa = 1d0/sqrt(aaa)
           rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
        enddo
     enddo
  endif
  !! Bessel test block end ccccccccccccccccccccccccccccccc


  nbloch    = sum(nblocha)
  nblochpmx = nbloch + ngcmx

  ! --- CG coefficienets. <LM3|lm1 lm2> (confusing, use subroutine clebsh)
  ! inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
  ! Injcg = indxcg(inxcg) to indxcg(inxcg)-1
  ! cg(inxcg)  : = <lm3|lm1 lm2>
  ! jcg(lnjcg) : = lm3
  lmxcg = lxx
  call scg_sizechk(lmxcg,lnjcg,lnxcg) !(lmax,c,cindx,js)
  write(6,*)'scg_sizechk= ',lnjcg,lnxcg
  !      if (lmxcg .le. 6) then
  !        lnjcg = 6500
  !        lnxcg = 1300
  !      else if (lmxcg .le. 8) then
  !        lnjcg = 22700
  !        lnxcg = 3400
  !      else if (lmxcg .le. 10) then
  !        lnjcg = 62200
  !        lnxcg = 7400
  !      else
  !        call rxi('setcg: cannot handle lmxcg=',lmxcg)
  !      endif
  !      if(allochk)
  !     &  write(*,*) 'allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))'
  allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
  call scg(lmxcg,cg,indxcg,jcg)
  if(allochk) write(6,*)' end of scg: cg coefficients generated.'

  ! --- Get real-space vectors and reciprocal-space vectors for Ewald sum.
  ! defaults values for ewald sum
  !      call lattc(awald0,tol,alat,alat,plat0,gx,gy,gz,gam,plat,qlat,
  !     .   lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx,w(owork))
  !- taken from lattc.f

  !! Apr2021 gomi. we still have a problem of posititve definiteness of the
  !!      Coulomb matrix. Probably because of this routine
  ! default values ok?
  awald0 = 2d0   ! See p_lat_0
  tol    = 1d-14 ! It was 1d-9 before aug2019. Sakakibara had a problem when he treat slab model of KF
  ! for 18 atoms per cell. --> The lowest eigenvalue of coulomb matrix becomes negative.
  nkdmx  = 8000  !    800->8000 to treat GdN. !31oct2019
  nkqmx  = 8000  !    nkdmx and nqkmx should be determined automatically in future.
  lmax   = 2*lxx  !lxx or lmax=6 ???

  vol0= abs(tripl(plat,plat(1,2),plat(1,3)))
  as   = awald0
  alat1= alat
  tpiba=8.d0*datan(1.d0)/alat
  call cross(plat(1,2),plat(1,3),qb0)
  call cross(plat(1,3),plat(1,1),qb0(1,2))
  call cross(plat(1,1),plat(1,2),qb0(1,3))
  qb0(1:3,1:3) = qb0(1:3,1:3)/vol0

  rdist0=vol0**(1.d0/3.d0)
  qdist0=1.d0/rdist0
  radd=.7*rdist0
  qadd=.7*qdist0
  a0=as/rdist0
  awald=a0/alat
  tol1= tol*alat**(lmax+1)
  allocate(dlv(3,nkdmx), qlv(3,nkqmx), work(max0(nkdmx,nkqmx)) )
  call lctoff(a0,vol0,lmax,tol1,r0,q0)
  nkdest =4.18879*(r0+radd)**3/vol0+.5
  nkrest =4.18879*(q0+qadd)**3*vol0+.5
  write(6,340) as,tol,lmax,awald,vol0,alat1,nkdest,nkrest
340 format(/' lattc:  as=',f6.3,'   tol=',1p,e9.2,'   lmax=',i2, &
       '   awald=',0p,f7.4,'   v0=',f10.3/' alat1=',f9.5, &
       '   estimates:   nkd',i6,'   nkr',i6)
  call lgen(plat,r0+radd,nkd,nkdmx,dlv,work)
  write(6,342) r0,r0*alat,radd,nkd
342 format('  r0=',f9.4,'   rc=',f9.4,'   radd=',f9.4,'   nkd=', i7)
  call lgen(qb0,q0+qadd,nkq,nkqmx,qlv,work)
  write(6,341) q0,q0*tpiba,qadd,nkq
341 format('  q0=',f9.4,'   qc=',f9.4,'   qadd=',f9.4,'   nkr=', i7)
  deallocate(work)
  !!
  eee=screenfac()           !takao feb2012

  !! for eps_lmf and epsPP_lmf mode,
  !! even the small eee=1d-4 can affect to dielectric function near q=0 when its values is large as one-hundred or more.
  !! Thus we set eee=0d0 to avoid this.
  if(imode==202) then !
     eee=0d0
  endif
  write(6,"(' Coulomb is exp(sqrt(-eee)*r)/r. eee=',d13.6,d13.6)")eee

  !--- bessel and hankel for the expansion of exp(-r/r_0)/r.
  ! bessel and hankel is renomarized so that its behaves as r^l and r^{-l-1} near r=0.
  !  rkpr means r^l*r for e=0 (r0c =infinity) case
  allocate(rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas),rofi(nrx,nbas))
  do ibas=1,nbas
     call genjh(eee,nr(ibas),aa(ibas),bb(ibas),lx(ibas), nrx,lxx, &
          rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas))
  enddo

  !--- onsite integrals <j(e=0)|B> and <B|v(onsite)|B>
  allocate( rojb(nxx, 0:lxx, nbas), sgbb(nxx,  nxx,  0:lxx, nbas))
  do ibas = 1,nbas
     call mkjb_4( lxx, lx(ibas),nxx, nx(0:lxx,ibas), &
          aa(ibas),bb(ibas), nr(ibas), nrx, &
          rprodx(1,1,0,ibas), &
          rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas), &
          rojb(1,0,ibas), sgbb(1,1,0,ibas))
  enddo

  !----------------
  !--- coulomb matrix for each q = qibz
  !----------------
  nlxx= (lxx+1)**2
  allocate(ngvecc0(3,ngcmx))
  call readqg('QGcou',(/0d0,0d0,0d0/),  quu,ngc0, ngvecc0)
  deallocate(ngvecc0)
  ngb = nbloch + ngc0
  if(allochk) write(*,*) 'allocate( vcoul)'
  allocate( vcoul(nblochpmx,nblochpmx) )
  vcoul  = 0d0

  !... q near zero
  !$$$      write(6,*) '--- readin Q0P -------'
  !$$$      open (newunit=ifq0p,file='Q0P')
  !$$$      read (ifq0p,*) nq0i,idummy,nq0iadd
  !$$$      if(allochk) write(*,*)'allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )'
  !$$$      allocate( wqt(1:nq0i+nq0iadd),q0i(1:3,1:nq0i+nq0iadd) )
  !$$$      do i=1,nq0i+nq0iadd
  !$$$         read (ifq0p,*) wqt(i),q0i(1:3,i)
  !$$$      enddo
  !$$$      write (6,"(d13.5,3x, 3d13.5)" )( wqt(i),q0i(1:3,i),i=1,nq0i+nq0iadd)
  !$$$      close(ifq0p)

  !      wvcc=.false.
  write(6,'(a)') " Mix0vec.XXX is not empty only when" &
       //" the corresponding q is in Q0P with zero weight."
  iqxend = nqibz + nq0i+nq0iadd
  !      if(wvcc) open(newunit=ifvcfpout,file='VCCFP')
  open(newunit=ifgb0vec, file="Mix0vec")
  open(newunit=ifgb0vec1,file="Mix0vec1")

  if(imode==202) then
     iqxini= nqibz + 1
  else
     iqxini = 1
  endif
  if(imode==0) then
     iqxini=1
  endif
  write(6,*)'iqxini iqxend=',iqxini,iqxend
  if(abs(sum(qibz(:,1)**2))/=0d0) call rx( 'hvccfp0: sum(q**2)==0d0')
  !     if(wvcc) write(ifvcfpout) nqbz, nblochpmx


  !... Readin PRODMT into prodmt. oct2005
  smbb = smbasis() !smabsis =F usually, may need fixing in future
  write(6,*) ' smooth mixed basis=',smbb
  if(smbasis()) then
     allocate( prodmt(2,nxx,0:lxx,nbas))
     allocate( nx_r(0:lxx))
     do ibas =1,nbas
        !          filenamep = 'PRODMT_'//charnum3(ibas)
        !          ifprodmt  = iopen(filenamep,0,-1,0)
        open(newunit=ifprodmt,file=trim('PRODMT_'//charnum3(ibas)))
        read(ifprodmt) nl_r
        if( 2*(nl_r-1) /= lxx ) then
           write(6,*) 2*(nl_r-1),lxx
           call rx( '2*nl_r-1 /= lxx ')
        endif
        read(ifprodmt) nxx_r
        write(6,"(' nxx =',100i3)")nxx_r
        if(nxx_r>nxx) call rx( 'nxx_r>nxx')
        read(ifprodmt) nx_r(0:lxx)
        write(6,"(' nx_r=',100i3)") nx_r(0:lxx)
        lx_ = lx(ibas)
        if(sum(abs(nx(0:lx_,ibas)-nx_r(0:lx_))) /=0) then
           write(6,*)' debug: nx  =',nx(0:lx_,ibas)
           write(6,*)' debug: nx_r=',nx_r(0:lx_)
           call rx( 'nx /=nx_r')
        endif
        read(ifprodmt) prodmt(1:2, 1:nxx_r, 0:lxx, ibas)
        write(6,*)' sumcheck prodmt=',sum(abs(prodmt(:,:,:,ibas)))
        close(ifprodmt) !isx = iclose(filenamep)
     enddo
     !!... Check write for radial part of the product basis
     if( .FALSE. ) then
        do ibas= 1,1 !1,nbas
           do l   =  0,lx(ibas)
              if1011=ifile_handle()
              open(if1011,file='ProdOld_ibas'//charnum3(ibas)//'_l'//charnum3(l))
              nxdim = nx(l,ibas)
              do ix=1,nxdim
                 write(if1011,"(' -- -- -- ',3i3,' --- ' )") ix,l,ibas
                 do ir =1,nr(ibas)
                    write(if1011,"(d13.5,2x,2d18.8)") &
                         rofi(ir,ibas), rprodx(ir,ix,l,ibas) &
                         , rprodx(ir,ix,l,ibas) /rofi(ir,ibas)
                 enddo
              enddo
              close(if1011)
           enddo
        enddo
        !       stop 'text end'
     endif
     allocate( rdmatch(nxx,nxx,0:lxx,nbas) )
     do ibas= 1, nbas
        do l   = 0, lx(ibas)
           nxdim = nx(l,ibas)
           if(nxdim<=1)write(6,*)'hvccfp0:smbasis case error nxdim <=1'
           !       pval  = prodmt(1, 1:nxdim, l,ibas)
           !       pslo  = prodmt(2, 1:nxdim, l,ibas)
           !       prod(r, inew) = \sum_iold rrmat(inew,iold) * prod(r,iold)
           write(6,"('goto mkradmatch ibas lnxdim =',3i4)")ibas,l,nxdim
           call mkradmatch(prodmt(1:2, 1:nxdim, l,ibas), nxdim, &
                rdmatch(1:nxdim,1:nxdim,l,ibas) )
        enddo
     enddo
     ! index (mx,nx,lx,ibas) ordering: taken from voul_4
     allocate(ibl(-lxx:lxx,nxx,0:lxx,nbas))
     ibl1 = 0
     ibl=999999
     do ibas= 1, nbas
        do l   = 0, lx(ibas)
           do n   = 1, nx(l,ibas)
              do m   = -l, l
                 ibl1  = ibl1 + 1
                 ibl(m,n,l,ibas) = ibl1
                 !       write(6,*)ibl1,n,l,m,lmbl(ibl1)
              enddo
           enddo
        enddo
     enddo
     if(ibl1/= nbloch) then
        write(6,*)' ibl1 nbloch',ibl1, nbloch
        call rx( ' hvccfp0:smbasis mode  error ibl1/= nbloch')
     endif
     ! index (mx,nx,lx,ibas) ordering
     nnr = 2 ! =2 new
     ! =0 equivalence with original mixed basis
     write(6,*)' sss:nbas lx=',nbas,lx(1:nbas)
     nbln=0
     do ibas= 1, nbas
        do l   = 0, lx(ibas)
           write(6,"('sss: nx=',3i4)") ibas,l,nx(l,ibas)
           if(nx(l,ibas)<=0) cycle
           if(nx(l,ibas)==1) call rx( 'nx(l,ibas) =1')
           nbln = nbln + (2*l+1)*(nx(l,ibas)-nnr)
        enddo
     enddo
     allocate( pmat(nbloch+ngcmx, nbln+ngcmx) )
     pmat=0d0
     ibln = 0
     do ibas= 1, nbas
        do l   = 0, lx(ibas)
           do nn  = nnr+1, nx(l,ibas) !nn=1 and nn=2 corresponds to non-zero val sol
              do m   = -l, l
                 ibln = ibln +1
                 nxdim = nx(l,ibas)
                 pmat( ibl(m,1:nxdim,l,ibas), ibln) &
                      =  rdmatch(1:nxdim, nn, l,ibas)
              enddo
           enddo
        enddo
     enddo
     !... Store matting matrix (imatchn,imatcho,pmatch)
     !        ifpomat   = iopen('POmat',0,-1,0)
     open(newunit=ifpomat,file='POmat',form='unformatted')
     !       write(6,*)'ttt= sumchk pmat(b)=',sum(abs(pmat(1:nbloch, 1:nbln)))
  endif

  !! Vcoud file, which contains E(\nu,I), given in PRB81,125102
  !! == main loop for iqx ==
  call MPI__getRange( mpi__iini, mpi__iend, iqxini, iqxend )
  do 1001 iqx = mpi__iini, mpi__iend ! q=(0,0,0) is omitted!
     write(6,"('#### do 1001 start iqx=',5i5)")iqx,nqibz
     vcoudfile='Vcoud.'//i2char(iqx)  !this is closed at the end of do 1001
     open(newunit=ifvcoud,file=trim(vcoudfile),form='unformatted')
     if(iqx > nqibz) then  !       iq = 1
        q  = q0i(:,iqx-nqibz)
     else                  !       iq = iqx
        q  = qibz(:,iqx)
     endif
     !        if(abs(sum(q))<1d-8) cycle
     if(imode==202 .AND. abs(sum(q))<1d-8) cycle
     !! q+G vector
     call readqg('QGcou',q,  quu,ngc, ngvecc ) !qq-->q
     ngb = nbloch + ngc  !it was ngcnn(iq)
     write(6,'(" iqx q ngc =",i5,3f10.4,i5)') iqx,q,ngc
     !! strxq: structure factor.
     allocate( strx(nlxx,nbas,nlxx,nbas))
     do ibas1 =1,nbas
        do ibas2 =1,nbas
           p = bas(:,ibas2)-bas(:,ibas1)
           phasep =exp(img*2*pi*sum(q*p))
           nlx1 = (lx(ibas1)+1)**2
           nlx2 = (lx(ibas2)+1)**2
           if(allochk) write(*,*) 'allocate( s(nlx1,nlx2))'
           allocate( s(nlx1,nlx2),sd(nlx1,nlx2)) !kino add sd----but sd is dummy
           call strxq(1,eee,q,p,nlx1,nlx2,nlx1,alat,voltot, &
                awald,nkd,nkq,dlv,qlv, &
                cg,indxcg,jcg, &
                s,sd)
           strx(1:nlx1,ibas1,1:nlx2,ibas2) = fpi*s      !!! *phasep
           if(allochk) write(*,*)'deallocate( s )'
           deallocate( s,sd )
        enddo
     enddo
     !!  onsite integrals <j(e=0)|P^(q+G)_L> and <B|v(onsite)|B>
     if(allochk) write(*,*)'allocate(rojp,sgpb,fouvb)'
     allocate(  rojp(ngc,      nlxx, nbas), &
          sgpb(ngc, nxx, nlxx, nbas), &
          fouvb(ngc, nxx, nlxx, nbas))
     do ibas = 1,nbas
        if(allochk) write(6,*)' --- goto mkjp_4',ibas
        call mkjp_4(q,ngc, ngvecc, alat, qlat, &
             lxx, lx(ibas),nxx, nx(0:lxx,ibas), &
             bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas), &
             nr(ibas), nrx, rprodx(1,1,0,ibas), &
             eee, rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas), &
             rojp(1,1,ibas),  sgpb(1,1,1,ibas), &
             fouvb(1,1,1,ibas))
     enddo

     !--- the Coulomb matrix
     if(allochk) write(6,*)' goto vcoulq_4'
     call vcoulq_4(q, nbloch, ngc, &
          nbas, lx,lxx, nx,nxx, &
          alat, qlat, voltot, ngvecc, &
          strx, rojp,rojb, sgbb,sgpb, fouvb, nblochpmx, bas,rmax, &
          eee, aa,bb,nr,nrx,rkpr,rkmr,rofi, &
          vcoul)
     if(allochk) write(6,*)' end of vcoulq_4'
     deallocate( strx, rojp,sgpb,fouvb)

     !----check write
     trwv = 0d0
     do i = 1,nbloch
        trwv = trwv + vcoul(i,i)
     enddo
     write(6,'(" vcoul trwi=",i6,2d22.14)') iqx,trwv
     write(6,'("### sum vcoul(1:ngb,      1:ngb) ",2d22.14,2x,d22.14)') &
          sum(vcoul(1:ngb,1:ngb)), sum(abs(vcoul(1:ngb,1:ngb)))
     write(6,'("### sum vcoul(1:nbloch,1:nbloch) ",2d22.14,2x,d22.14)') &
          sum(vcoul(1:nbloch,1:nbloch)),sum(abs(vcoul(1:nbloch,1:nbloch)))
     write(6,*)
1101 continue
     ngbo=ngb

     !!.. Generate ppmt mattix oct2005 .......................
     !!   smbasis option is not used now-- we may recover in future
     if(smbasis()) then
        allocate( ppmt(2,(lxx+1)**2,nbas,ngc) )
        ppmt = 0d0
        call mkppmt(alat,plat,qlat, q, &
             ngc, ngvecc, &
             rmax, nbas,  bas, lx, lxx, &
             ppmt) ! ppmt contains value and slove of e(i q+G r) at MT boundaries.
        ! ppmt(2,lmxaa,nbas)
        write(6,*) 'nbln ngc',nbln,ngc

        !... Matching matrix pmtch. ppmt and prodmt
        pmat(:, nbln+1:nbln+ngc)=0d0
        do igc=1,ngc
           pmat(nbloch+igc, nbln+igc) = 1d0
           do ibas= 1, nbas
              do l  =  0, lx(ibas)
                 do m  = -l, l
                    pval= ppmt(1, l**2 + l+1 +m, ibas,igc)
                    pslo= ppmt(2, l**2 + l+1 +m, ibas,igc)
                    do n = 1,nx(l,ibas)
                       if(n==1 .AND. debug) write(6,"('ttt2: ')")
                       pmat(ibl(m,n,l,ibas), nbln+igc) &
                            =  rdmatch(n,1,l,ibas) * pval &
                            +  rdmatch(n,2,l,ibas) * pslo
                       if(debug .AND. abs(pmat(ibl(m,n,l,ibas), nbln+igc))/=0d0) &
                            write(6,"('ttt2: i1 i2 pmat=',2i5,2d13.5)") &
                            ibl(m,n,l,ibas), nbln+igc, pmat(ibl(m,n,l,ibas), nbln+igc)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        deallocate(ppmt)
        nn = nbln  +ngc ! number for new smooth mixed basis.
        no = nbloch+ngc ! number for original size of mixed basis.
        if(debug) write(6,*) 'end of pmat'

        !... oo(no,no). The original overlap matrix.
        allocate( pomat(nn,no) )
        allocate( ppovl(ngc,ngc),oo(no,no))
        call mkppovl2(alat,plat,qlat, &
             ngc,  ngvecc, &
             ngc,  ngvecc, &
             nbas, rmax, bas, &
             ppovl)
        oo = 0d0
        do ipl1 = 1,nbloch
           oo(ipl1,ipl1) = 1d0
        enddo
        do ix= 1,ngc
           do iy= 1,ngc
              oo(nbloch+ix, nbloch+iy) = ppovl(ix,iy)
           enddo
        enddo
        if(debug) write(6,*) 'end of oo'
        !... oon(nn,nn) is the overlap matrix with new basis
        allocate(oon(nn,nn))
        oon = matmul( dconjg(transpose(pmat(1:no,1:nn))) &
             ,matmul(oo,pmat(1:no,1:nn)) )
        if( .FALSE. ) then
           open(newunit=if3011,file='oontest'//i2char(iqx))
           do ix=nbln+1,nn
              igc=ix-nbln
              qqx(1:3) = (q(1:3)+ matmul(qlat, ngvecc(1:3,igc)))
              absqq  =  sqrt(sum(qqx(1:3)**2))
              do iy=nbln+ 1,nn
                 igc2=iy-nbln
                 if(ix==iy) then
                    write(if3011,"('on : ',2i8,3i3,2x,3i3,f13.5,3x,2f20.10)") &
                         ix,iy, ngvecc(1:3,igc),ngvecc(1:3,igc2),absqq, oon(ix,iy)
                 else
                    write(if3011,"('off:', 2i8,3i3, 2f20.10)")ix,iy, &
                         ngvecc(1:3,igc)-ngvecc(1:3,igc2),  oon(ix,iy)
                 endif
              enddo
           enddo
           close(if3011)
        endif
        !... Generat pomat
        !    zmelt_new(K, ij) = \sum_I pomat(K,I)* zmelt(I, ij)
        !    means <psi_i psi_j | K> where |K> denote new mixed basis.
        !  See sxcf_fal2 and x0kf.
        !   Be carefull its transpose procedure---it is a little confusing...
        call pmatorth(oo,oon, pmat(1:no,1:nn), no, nn, &
             pomat)
        if( iqx <= nqibz ) deallocate(oon)
        deallocate(ppovl,oo)
        !... Store matching matrix
        write(ifpomat) q,nn,no,iqx
        write(ifpomat) pomat
        deallocate(pomat)
     endif

     !! == Write out VCCFP ==
     if(debug) write(6,*) 'write out vcoul'
     if(smbasis()) then
        ngb= nn
        allocate(vcoulnn(ngb,ngb))
        vcoulnn= matmul(transpose(dconjg(pmat(1:no,1:nn))) &
             ,matmul(vcoul(1:no,1:no),pmat(1:no,1:nn)))
        vcoul(1:ngb,1:ngb)= vcoulnn
        deallocate(vcoulnn)
     endif
     !        if(wvcc) then
     !          write(ifvcfpout) ngb
     !          write(ifvcfpout) vcoul(1:ngb,1:ngb),q
     !        endif
     write(6,"(' ngc ngb/ngbo=',6i6)") ngc,ngb,ngbo

     !! Mix0vec ---------------------------------
     !! diagonalize the Coulomb matrix
     !        if(.true.) then
     if(allochk) write(*,*) 'allocate( ppovl(ngc,ngc))'
     allocate( oo(ngb,ngb) )
     allocate( ppovl(ngc,ngc) )
     call mkppovl2(alat,plat,qlat, &
          ngc,  ngvecc, &
          ngc,  ngvecc, &
          nbas, rmax, bas, &
          ppovl)
     if(smbasis()) then
        oo = oon
        deallocate(oon)
     else
        oo = 0d0
        do ipl1=1,nbloch
           oo(ipl1,ipl1) = 1d0
        enddo
        do ix=1,ngc
           do iy=1,ngc
              oo(nbloch+ix, nbloch+iy) = ppovl(ix,iy)
           enddo
        enddo
     endif

     allocate( oox(ngb,ngb) )
     oox = oo
     write(6,*)' --- goto eigen check1 --- '
     allocate(  vcoul0(ngb,ngb) )
     vcoul0 = vcoul(1:ngb,1:ngb)
     if(allochk) &
          write(*,*) 'allocate(hh(ngb,ngb),oo(ngb,ngb),oox,zz,eb,zzr)'
     allocate(hh(ngb,ngb),zz(ngb,ngb),eb(ngb),zzr(ngb))
     hh  = - vcoul0
     !          nmx = 15
     nmx = ngb
     call diagcv(oo,hh,zz,ngb, eb,nmx,1d99,nev)
     do ipl1=1,nev
        if(ipl1==11) write(6,*)' ... '
        if(ipl1>10 .AND. ipl1<nev-5) cycle
        write(6,'(i4,d23.16)')ipl1,-eb(ipl1)
     enddo
     write(6,"(' nev ngv q=',2i5,3f10.6)")nev,ngb,q

     !! -eb should be positive definite. However, we have one (or a few?) negative ones.
     !! I(kotani) think no problem to set eb=0 when -eb is negative.
     !! But this is a temporaly fix or better manner to calcuate coulomb matrix. strxq may be needed to be replaced
     do i=1,nev
        if(eb(i)>0 .AND. keeppositivecou) then
           write(6,"(a,d13.5)") &
                'KeepPositiveCou enforce : -eb<0 --> eb=0',eb(i)
           eb(i)=0d0
        endif
     enddo
     write(ifvcoud) ngb
     write(ifvcoud) q
     write(ifvcoud) -eb
     write(ifvcoud) zz
     write(6,*)
     write(6,'(" eig0 must be equal to the largest =", 2d24.16)') &
          sum(  dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1))  )
     write(6,'(" zz norm check=",d24.16)') &
          sum( dconjg(zz(1:ngb,1))*matmul(oox,zz(1:ngb,1)) )
     write(6,*)
     !          write(6,'(" --- vcoul(exact  no eee)=",d14.6," absq2=",d24.16)')
     !     &    fpi*voltot/(sum(tpiba**2*q(1:3)**2))
     !     &             , (sum(tpiba**2*q(1:3)**2))
     write(6,'(" --- vcoul(exact)=",d14.6," absq2=",d24.16)') &
          fpi*voltot/(sum(tpiba**2*q(1:3)**2)-eee) &
          , (sum(tpiba**2*q(1:3)**2)-eee)
     write(6,'(" --- vcoul(cal ) =",2d14.6)') &
          sum( dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1)) )*voltot
     ! cccccccccccccccccccccccccccccccccccccccc
     !          do igc=1,ngb
     !          qqx(1:3) = (q(1:3)+ matmul(qlat, ngvecc(1:3,igc)))
     !          write(6,'(" --- vcoul(exact) xxx =",d14.6," absq2=",d24.16)')
     !     &    fpi*voltot/(sum(tpiba**2*(qqx(1:3)**2)-eee))
     !     &             , (sum(tpiba**2*(qqx(1:3)**2)-eee))
     !          write(6,'(" --- vcoul(cal ) xxx =",2d14.6)')
     !     &    sum( dconjg(zz(1:ngb,igc))*matmul( vcoul0,zz(1:ngb,igc)) )*voltot
     !          enddo
     ! cccccccccccccccccccccccccccccccccccccccc
     deallocate(vcoul0)

     if( iqx-nqibz>=1 ) then
        if( wqt(iqx-nqibz)==0d0) then ! MIZUHO-IR
           !! --- To get the vector <Mixed basis| q=0> --------------
           if( .NOT. is_mix0vec()) then     !used original befor oct2006
              ! See switch.F ---> this is not used now.
              ifgb0vec_a =ifgb0vec1
              ifgb0vec_b =ifgb0vec
           else
              ! ismix0vec=1 is to avoid problem at BZ boundary when is_mix0vec()=0.
              ifgb0vec_a =ifgb0vec
              ifgb0vec_b =ifgb0vec1
           endif
           !1... Case1 to write ifgb0vec -------------------------------------
           write(6,*)' voltot=',voltot
           if(ngc==0) then
              continue
           else
              do igc=1,ngc
                 if( sum(abs( ngvecc(1:3,igc) ))==0 ) then
                    igc0=igc
                    exit
                 endif
              enddo
              write(6,*)' igc0=',igc0,ngvecc(1:3,igc0)
              zzr(nbloch+1:nbloch+ngc) = ppovl(1:ngc,igc0)
           endif
           allocate( gbvec(ngb), b0mat(nbloch) )
           !! ... get a vector <Product Basis| q+0>
           call mkb0( q, lxx,lx,nxx,nx, aa,bb,nr,nrx,rprodx, &
                alat,bas,nbas,nbloch, &
                b0mat)
           zzr(1:nbloch) = b0mat(1:nbloch)
           allocate(ooxi(ngb,ngb))
           ooxi=oox
           call matcinv(ngb,ooxi)
           gbvec = matmul(ooxi, zzr)
           deallocate(ooxi)
           dnorm = sqrt( sum(dconjg(gbvec)*zzr) )
           write(ifgb0vec_a,"(3d24.16,2i10,d24.16)") q, ngb,igc0,dnorm
           write(ifgb0vec_a,"(4d24.16)") (gbvec(i),zzr(i),i=1,ngb)
           deallocate( gbvec, b0mat)
           !1----------------------------------------------------
           !2... --- Case2 to write ifgb0vec c2 is problematic at BZ boundary...------
           dnorm  = 1d0
           zzr(:) = matmul (oox, zz(:,1))
           igc0 = 999999 !dummy now
           ! phasex ---just to clean. this is irrelevant
           phasex =1d0
           do i=1,ngb
              if(abs(zz(i,1)) > 1d-3) phasex = abs(zz(i,1))/zz(i,1)
           enddo
           do i=1,ngb
              zz(i,1)= phasex * zz(i,1)
              zzr(i) = phasex * zzr(i)
           enddo
           write (ifgb0vec_b,"(3d24.16,2i10,d24.16)") q, ngb,igc0,dnorm
           write (ifgb0vec_b,"(4d24.16)") (zz(i,1),zzr(i),i=1,ngb)
        endif
     endif ! MIZUHO-IR
     deallocate(hh,oo,zz,eb,oox,zzr)
     deallocate(ppovl)
     !2---------------------
     !        endif
     close(ifvcoud)
1001 enddo
  deallocate(ngvecc)
  call cputid(0)
  !      call MPI__Finalize
  if(imode==202) call rx0( ' OK! hvccfp0 imode=202 only for Q0P')
  if(imode==0) call rx0( ' OK! hvccfp0 imode=0')
  if(imode==3) call rx0( ' OK! hvccfp0 imode=3')
END PROGRAM hvccfp0


subroutine checkagree(a,b,char)
  real(8):: a(3),b(3)
  character*(*) :: char
  if(sum(abs(a-b))>1d-6) then
     write(6,*)' Error in checkagree:',char
     ! top2rx 2013.08.09 kino        stop ' Error in checkagree:'
     call rx( ' Error in checkagree:')
  endif
end subroutine checkagree

subroutine mkradmatch( p, nxdim, &
     rdmatch)
  !- make rdmatch
  !----------------------------------------------------
  !i  p(1,i): phi     at mt for i-th basis
  !i  p(2,i): dphi/dr at mt for i-th basis
  !o rdmatch(nxdim,nxdim)
  !-------
  !r    phinew_j(r) =sum_i phi_i(r)* rdmatch (i,j)
  !r     phinew_1(rmt)    =1      phinew_2(rmt)   =0
  !r   d phinew_1(rmt)/dr =0    d phinew_2(rmt)/dr=1
  !r for k >=3
  !r     phinew_k(rmt)    =0
  !r   d phinew_k(rmt)/dr =0
  !----------------------------------------------------
  implicit none
  integer(4):: nxdim,lbas,i,i1,i2,ix
  real(8):: p(1:2, 1:nxdim), rdmatch(1:nxdim,1:nxdim)
  real(8):: pd,p1,p1d,p2,p2d,s,t, eps=1d-3,delta
  !r                                       old     new
  !      write(6,"('mkradmatch: nxdim=',i4)") nxdim
  if(nxdim <=0) return
  ! top2rx 2013.08.09 kino      if(nxdim ==1) stop 'mkradmatch err nxdim==1'
  if(nxdim ==1) call rx( 'mkradmatch err nxdim==1')
  rdmatch=0d0
  !... pivot--- get better set of phi for augmentation
  do
     i1= nxdim
     i2= nxdim-1
     p1 = p(1, i1)
     p2 = p(1, i2)
     p1d= p(2, i1)
     p2d= p(2, i2)
     write(6,"('mkradmatch: i1 p1 p1d=',i3,2d13.6)") i1,p1,p1d
     write(6,"('mkradmatch: i2 p2 p2d=',i3,2d13.6)") i2,p2,p2d
     delta = p1*p2d-p2*p1d
     if(abs(delta) <eps*p1*p2) then
        if(i2==1) then
           write(6,"(' i1 i2=',2i5,2d13.6)") i1,i2,p1d/p1,p2d/p2
           ! top2rx 2013.08.09 kino            stop'mkradmatch: err poor linear dep'
           call rx( 'mkradmatch: err poor linear dep')
        endif
        i2=i2-1
     endif
     exit
  enddo
  !...
  call phimatch(1d0,0d0,  p1,p1d,p2,p2d, s,t)
  rdmatch(i1, 1)=  s
  rdmatch(i2, 1)=  t
  write(6,"('mkradmatch: 1 0    st=',2d13.5)") s,t
  call phimatch(0d0,1d0,  p1,p1d,p2,p2d, s,t)
  rdmatch(i1, 2)=  s
  rdmatch(i2, 2)=  t
  write(6,"('mkradmatch: 0 1    st=',2d13.5)") s,t

  ix=2
  do i= 1,nxdim
     if(i==i1 .OR. i==i2) cycle
     ix=ix+1
     !        write(6,"('mkradmatch: i p pd=',i3,2d13.5)") i,p(1,i),p(2,i)
     call phimatch(p(1,i),p(2,i),  p1,p1d,p2,p2d, s,t)
     rdmatch(i,  ix)=  1d0
     rdmatch(i1, ix)=  -s
     rdmatch(i2, ix)=  -t
     write(6,"('mkradmatch: ix st=',i3,2d13.5)") ix,s,t
  enddo
end subroutine mkradmatch

subroutine phimatch(p,pd, p1,p1d,p2,p2d, s,t)
  ! --- match for given p and pd
  !   phi = s phi1 + t phi2 !slope and value are at MT
  !     p  = s p1  + t p2
  !     pd = s pd1 + t pd2
  implicit none
  real(8):: matinv(2,2),p,pd,p1,p1d,p2,p2d,s,t,delta,ddd1,ddd2
  delta = p1*p2d-p2*p1d
  matinv(1,1) = 1/delta *  p2d
  matinv(1,2) = 1/delta * (-p2)
  matinv(2,1) = 1/delta * (-p1d)
  matinv(2,2) = 1/delta *  p1
  s = matinv(1,1) *p  + matinv(1,2) *pd
  t = matinv(2,1) *p  + matinv(2,2) *pd
  !... check
  ddd1 = abs(s*p1  + t*p2   -  p )
  ! top2rx 2013.08.09 kino      if(  ddd1 >1d-8 ) stop 'phimatch: ddd1 err'
  if(  ddd1 >1d-8 ) call rx( 'phimatch: ddd1 err')
  ddd2 = abs(s*p1d + t*p2d  -  pd)
  ! top2rx 2013.08.09 kino      if(  ddd2 >1d-8 ) stop 'phimatch: ddd2 err'
  if(  ddd2 >1d-8 ) call rx( 'phimatch: ddd2 err')
end subroutine phimatch

subroutine pmatorth(oo,oon,pmat,no,nn, pomat)
  ! get conversion matrix from old mixed basis(no) to augmented mixed basis(nn).
  ! pmatorth contains
  !   oo^{-1}_IJ
  implicit none
  integer(4):: no,nn,io,in,i
  complex(8):: pmat(no,nn),pomat(nn,no),oo(no,no),oon(nn,nn)
  complex(8),allocatable:: ooninv(:,:)
  real(8),allocatable:: eb(:)
  allocate(ooninv(nn,nn))
  ooninv = oon
  call matcinv(nn,ooninv) !generate ooninv
  !      pomat = matmul(ooninv, matmul(dconjg(transpose(pmat)),oo))
  pomat = transpose (matmul( oo, matmul(pmat,ooninv)))
  deallocate(ooninv)
end subroutine pmatorth
!      allocate(pp(nn,nn),ppin(nn,nn),eb(nn),zz(nn,nn),zze(nn,nn))
!      ppin = pp
!      call diagcvh(ppin,nn,eb,zz)
!      do i=1,nn
!        zze(:,i) =  zz(:,i)* sqrt(eb(i))
!      enddo
!      pomat = matmul(pmat, matmul(zze,dconjg(transpose(zz))))

subroutine diagcvh(hh,ngb,eb,zz)
  implicit none
  integer(4):: nmx,nev,i,ngb
  complex(8):: hh(ngb,ngb),oo(ngb,ngb),zz(ngb,ngb)
  real(8):: eb(ngb)
  nmx=ngb
  oo = 0d0
  do i=1,ngb
     oo(i,i) = 1d0
  enddo
  call diagcv(oo,hh,zz,ngb, eb,nmx,1d99,nev)
  write(6,*)' diagcvv: ngb,nev=',ngb,nev
  do i=1,nev
     write(6,'(i4,d23.16)')i, eb(i)
  enddo
end subroutine diagcvh
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine zgesvdnn2(no,nn, nnmx,epsmx, &
     pmat, &
     nnn)
  ! pmat(no,nn) ---> pmat(no,nnn)
  ! o input          pmat(no,nn)
  ! o output reduced pmat(no,nnn)
  implicit none
  integer(4):: lwork,info,nn,no,nnn,nnmx,i
  complex(8)::  pmat(no,nn),uu(no,no),vt(nn,nn)
  real(8):: ss(nn),epsmx
  real(8),allocatable:: rwork(:)
  complex(8),allocatable:: work(:),vtt(:,:),pmatx(:,:)
  !      write(6,*)' sumchk pmat=',sum(abs(pmat(1:no,1:nn)))
  lwork=4*no
  allocate(work(LWORK),rwork(5*no),pmatx(no,nn))
  pmatx =pmat
  call zgesvd('A','A',no,nn,pmat,no,SS,UU,no,VT,nn,work,lwork,rwork,info)
  nnn=-999
  do i=1,nn
     write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
     !         vtt(i,:)=ss(i)*vt(i,:)
     if(nnn==-999 .AND. ss(i)<epsmx) nnn = i-1
  enddo
  !      write(6,*) 'nnn=',nnn
  ! top2rx 2013.08.09 kino      if(nnn==0) stop 'strange: nnn=0'
  if(nnn==0) call rx( 'strange: nnn=0')
  if(nnn>nnmx) nnn=nnmx
  pmat=pmatx
  !      pmat(:,1:nnn) = uu(:,1:nnn)
  !      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
  !     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
  !      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk)))
  !     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
  !      deallocate(vtt)
end subroutine zgesvdnn2


!---------------------------------------------------------------------
subroutine mkb0( q, lxx,lx,nxx,nx, aa,bb, nrr,nrx,rprodx, &
     alat,bas,nbas,nbloch, &
     b0mat)
  !--make the matrix elementes < B_q | exp(iq r)>
  use m_lldata,only: ll
  implicit none
  integer(4) :: nlx,l,n,m,nr,ir,lm,ibl1,ibas,nrx,nbloch

  integer(4) :: nbas,lxx, lx(nbas), nxx, nx(0:lxx,nbas),nrr(nbas)
  real(8)    :: rprodx(nrx,nxx,0:lxx,nbas),aa(nbas),bb(nbas), &
       phi(0:lxx),psi(0:lxx), bas(3,nbas), &
       alat, &
       pi,fpi,tpiba,qg1(3),q(3),absqg,r2s,a,b

  complex(8) :: b0mat(nbloch),img=(0d0,1d0) ,phase

  integer(4),allocatable:: ibasbl(:), nbl(:), lbl(:), lmbl(:)
  real(8),allocatable :: ajr(:,:),rofi(:),rob0(:,:,:)
  real(8),allocatable::cy(:),yl(:)
  complex(8),allocatable :: pjyl(:,:)
  !$$$#ifdef COMMONLL
  !$$$      integer(4) ll(51**2)
  !$$$      common/llblock/ll
  !$$$#else
  !$$$      integer(4) ll
  !$$$#endif

  !-----
  write(6,*)'mkb0:'
  pi   = 4d0*datan(1d0)
  fpi  = 4*pi
  nlx  = (lxx+1)**2

  tpiba = 2*pi/alat
  qg1(1:3) = tpiba * q(1:3)
  absqg    = sqrt(sum(qg1(1:3)**2))

  allocate(ajr(1:nrx,0:lxx), pjyl(nlx,nbas),rofi(nrx), &
       ibasbl(nbloch), nbl(nbloch), lbl(nbloch), lmbl(nbloch), &
       cy(nlx),yl(nlx),rob0(nxx,0:lxx,nbas))

  call sylmnc(cy,lxx)
  call sylm( qg1/absqg,yl,lxx,r2s) !spherical factor Y( q+G )

  do ibas = 1,nbas
     a = aa(ibas)
     b = bb(ibas)
     nr= nrr(ibas)
     rofi(1)    = 0d0
     do ir      = 1, nr
        rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
        call besslggg(absqg**2*rofi(ir)**2,lx(ibas),phi,psi)
        do l  = 0,lx(ibas)
           ! ... bessel function
           ajr(ir,l) = phi(l)* rofi(ir) **(l +1 )
           ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
        enddo
     enddo

     ! ... Coefficients for j_l yl  on MT  in the expantion of of exp(i q r).
     phase = exp( img*sum(qg1(1:3)*bas(1:3,ibas))*alat  )
     do lm = 1,(lx(ibas)+1)**2
        l = ll(lm)
        pjyl(lm,ibas) = fpi *img**l *cy(lm)*yl(lm) *phase  *absqg**l
     enddo
     ! ... rob0
     do l = 0,lx(ibas)
        do n = 1,nx(l,ibas)
           call gintxx( ajr(1,l), rprodx(1,n,l,ibas), a,b,nr, &
                rob0(n,l,ibas) )
        enddo
     enddo
  enddo

  ! ... index (mx,nx,lx,ibas) order.
  ibl1 = 0
  do ibas= 1, nbas
     do l   = 0, lx(ibas) ! write(6,'(" l ibas nx =",3i5)') l,nx(l,ibas),ibas
        do n   = 1, nx(l,ibas)
           do m   = -l, l
              ibl1  = ibl1 + 1
              ibasbl(ibl1) = ibas
              nbl   (ibl1) = n
              lbl   (ibl1) = l
              lmbl  (ibl1) = l**2 + l+1 +m ! write(6,*)ibl1,n,l,m,lmbl(ibl1)
           enddo
        enddo
     enddo
  enddo
  ! ... pjyl * rob0
  do ibl1= 1, nbloch
     ibas= ibasbl(ibl1)
     n   = nbl  (ibl1)
     l   = lbl  (ibl1)
     lm  = lmbl (ibl1)
     b0mat(ibl1) = pjyl(lm,ibas) * rob0(n,l,ibas)
  enddo
  deallocate(ajr, pjyl,rofi, &
       ibasbl, nbl, lbl, lmbl, &
       cy,yl,rob0)
end subroutine mkb0