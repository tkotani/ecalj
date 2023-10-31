program hvccfp0   ! Coulomb matrix. <f_i | v| f_j>_q.
  ! output  VCCFP : the coulomb matrix vcoul(nblochpmx,nblochpmx) for all qibz.
  !    strx: structure constant for e=0 (means 1/|r-r'| )
  use m_xlgen,only:lgen
  use m_lattic,only:lctoff
  use m_genallcf_v3,only: Genallcf_v3,  alat,nbas=>natom, bas=>pos
  use m_hamindex,only:    Readhamindex,plat,qlat
  use m_read_bzdata,only: Read_bzdata, ginv,nqbz,qbz,nqibz,qibz,nq0i,wqt=>wt,q0i,nq0iadd
  use m_readqg,only:   readqg,readngmx
  use m_mpi,only: MPI__hx0fp0_rankdivider2,mpi__task,MPI__Initialize,MPI__Finalize,mpi__root, &
       MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,mpi__rank,mpi__size, &
       mpi__ranktab,MPI__consoleout,mpi__iend,mpi__iini,mpi__getrange
  use m_readgwinput,only: ReadGwinputKeys, keeppositivecou
  use m_lgunit,only: m_lgunit_init
  use m_vcoulq,only: vcoulq_4,mkjb_4,mkjp_4,genjh
  use m_pwmat,only: mkppovl2
  use m_hvccfp0_util,only: mkradmatch,pmatorth,mkb0,strxq
  implicit none
  integer :: ifvcfpout,ifhvccfp,is,  lmxcg,if1011,if3011, ifplane,ngpmx, ngcmx, nblochpmx, nbloch,&
       ibas,ic,lxx,nxx,nrx,l,n,k,isx,kdummy, nkdmx,nkqmx,lmax,nkdest,nkrest,ngp,ngc,nlxx,i,lnjcg,lnxcg, &
       nkd,nkq ,ibas1,ibas2,nlx1,nlx2, iqibz,ir,ig1,n1,n2, ngb,nev,nmx,iqx,ipl1,ipl2,igx1,igx2
  integer:: igc,igc0,ifgb0vec,ifgb0vec1,ix, iy, iqxini, iqxend,imode, ngc0, ifvcfporg,nqbz_in,nblochpmx_in
  integer:: ifprodmt,nl_r,lx_,nxx_r,nxdim,ibl1,nn,no,ngbnew, nmatch,ifpmatch,nmatch_q,ifpmatch_q,m,ifpomat,nbln,ibln,ngb_in,nnr,igc2
  integer:: nnmx ,ngcnn,ngbo,ifgb0vec_a,ifgb0vec_b , ifvcoud,idummy
  integer:: ifiwqfac,iqbz,iqbzx,nnn,ixyz,ifq0p,incwfin 
  integer::  nqnumc,ifiqgc 
  real(8) ::  q(3),p(3),voltot, tripl,alat0,epsx, tol,as,tpiba,qb0(3,3),vol0,rdist0,qdist0,radd,qadd, &
       a0,awald,alat1,tol1,r0,q0,awald0,qg(3),   absqg2,aaa,aaa12
  real(8):: eee,eees, q_org(3),screenfac
  complex(8):: pval,pslo,phasex
  complex(8):: phasep,img=(0d0,1d0)
  complex(8) :: xxx,trwv
  real(8)::absqq,qqx(3), epsmx,aaaa, sss1,sss2,dnorm, qqq(3),QpGcut_Cou,quu(3)
  integer,allocatable :: jcg(:),indxcg(:), lx(:),kmx(:),nblocha(:),nr(:),ificrb(:), &
       nx(:,:),ngvecp(:,:),ngvecc(:,:),ngvecci(:,:,:),iqibzx(:)
  integer,allocatable:: ngvecc0(:,:)
  integer,allocatable:: nx_r(:), ibl(:,:,:,:),imatcho(:),imatchn(:),imatcho_q(:),imatchn_q(:)
  real(8),allocatable:: prodmt(:,:,:,:),rdmatch(:,:,:,:)
  real(8),allocatable :: rmax(:), cg(:),rprodx(:,:,:,:),dlv(:,:),qlv(:,:),work(:),ngcn(:), &
       rojb(:,:,:), sgbb(:,:,:,:),aa(:),bb(:),rofit(:),phi(:),psi(:)
  real(8),allocatable:: wqfac(:),qbzwww(:,:)
  real(8),allocatable :: rkpr(:,:,:),rkmr(:,:,:),rofi(:,:)
  real(8),allocatable    :: eb(:)
  complex(8) ,allocatable :: vcoul(:,:),geig(:,:),strx(:,:,:,:), &
       sgpb(:,:,:,:),sgpp(:,:,:,:), fouvb(:,:,:,:),fouvp(:,:,:,:),vcoul0(:,:), &
       s(:,:),sd(:,:),rojp(:,:,:) , vcoulnn(:,:)
  complex(8),allocatable:: gbvec(:), ppovl(:,:), b0mat(:), hh1(:,:),oo1(:,:), vcoul_org(:,:),matp(:),matp2(:)
  complex(8),allocatable:: ppmt(:,:,:,:),pmat(:,:),pomat(:,:),oon(:,:), hh(:,:),oox(:,:),ooxi(:,:),oo(:,:),zz(:,:),zzr(:)
  logical :: checkeig, besseltest=.false. ,allochk=.false.,smbasis,debug=.false.,smbb, is_mix0vec,wvcc, cmdopt2
  character(20) :: xxt,outs=''
  character(3) :: charnum3
  character(10) :: i2char
  character(128):: vcoudfile,ixcc
  real(8),parameter::pi  = 4d0*datan(1d0), fpi = 4d0*pi
  call MPI__Initialize()
  call M_lgunit_init()
  if( mpi__root) write(6,"(' mode=0,3,202 (0 and 3 give the same results for given bas)' )")
  if(cmdopt2('--job=',outs)) then; read(outs,*) imode
  elseif( mpi__root ) then       ; read(5,*) imode;   endif
  call MPI__Broadcast(imode)
  write(ixcc,"('.mode=',i4.4)")imode
  call MPI__consoleout('hvccfp0'//trim(ixcc))
  call cputid (0)
  if(imode==202 )  then;   write(6,*)' hvccfp0: imode=',imode
  elseif(imode==0) then
  elseif(imode==3) then
  else; call rx( 'hvccfp0: now hvccfp0 support just normal mode=0 3 202 101'); endif
  incwfin= 0                !use ForX0 for core in GWIN
  call Genallcf_v3(incwfin) !in module m_genallcf_v3
  call Readhamindex()
  call Read_bzdata()
  call readngmx('QGcou',ngcmx)
  allocate(ngvecc(3,ngcmx))
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  call ReadGWinputKeys() !Readin dataset in GWinput
  allocate(lx(nbas),kmx(nbas),nblocha(nbas),nr(nbas),aa(nbas),bb(nbas),ificrb(nbas),rmax(nbas) )
  do ibas = 1,nbas !! Readin BASFP//atom. The product basis functions.
     ic = ibas !
     open(newunit=ificrb(ibas),file=trim('BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10))))
     read(ificrb(ibas),"(4i6,2d24.16)") lx(ibas), kmx(ibas), nblocha(ibas), nr(ibas),aa(ibas),bb(ibas)
     rmax(ibas) = bb(ibas)*(exp((nr(ibas)-1)*aa(ibas))-1d0)
  enddo
  lxx = maxval(lx)
  allocate( nx(0:lxx,nbas) )
  do ibas = 1,nbas
     read(ificrb(ibas),"(i5)") nx(0:lx(ibas),ibas)
  enddo
  nxx = maxval(nx)
  nrx = maxval(nr)
  allocate( rprodx(nrx,nxx,0:lxx,nbas) )
  do ibas = 1,nbas
     do l = 0, lx(ibas)
        do n = 1, nx(l,ibas)
           read(ificrb(ibas),"(3i5)"   ) k, kdummy,kdummy
           read(ificrb(ibas),"(d23.15)") (rprodx(i,n,l,ibas),i=1,nr(ibas))
        enddo
     enddo
     close(ificrb(ibas))
  enddo
  BesselTestBlock: if(besseltest) then
     write(6,*)' *** TEST case ***  rprodx is given by Bessel.'
     iqx  = 2 ! test G, corresponding <q+G|v|q+G> should be exact. e.g. ig1=1 and  ig1=35 for iqx=2 
     igx1 = 1
     igx2 = 35
     write(6,"(' iqx=',i3,' ig1 ig2=',2i3)") iqx,igx1,igx2
     write(6,"(a)") ' <q+G|v|q+G> for the corresponding iqx ig1 ig2 should be exact!'
     write(6,"(a)") ' See fort.196'
     write(6,"(a)") ' Errors will be from the radial function integrals !!!'
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
     allocate(rprodx(nrx,nxx,0:lxx,nbas),rofit(nrx),phi(0:lxx),psi(0:lxx))
     rofit(1) = 0d0
     do ir = 1, nrx;    rofit(ir) = bb(1)*( exp(aa(1)*(ir-1)) - 1d0);   enddo
     do n = 1, nxx
        if(n==1) ig1 = igx1
        if(n==2) ig1 = igx2
        qg(1:3) = tpiba * (qibz(1:3,iqx)+ matmul(qlat, ngvecci(1:3,ig1,iqx)))
        absqg2  = sum(qg(1:3)**2)
        do ir =1,nrx
           call bessl(absqg2*rofit(ir)**2,lxx,phi,psi)
           do ibas=1,nbas
              do l = 0, lx(ibas)
                 rprodx(ir,n,l,ibas) = phi(l)* rofit(ir) **(l +1 )
              enddo
           enddo
        enddo
     enddo
     do ibas=1,nbas
        do l = 0, lx(ibas) !orthogonalized rprodx.
           rprodx(1:nr(ibas),1,l,ibas)= rprodx(1:nr(ibas),1,l,ibas) + rprodx(1:nr(ibas),2,l,ibas)
           n = 1
           call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas),aa(ibas),bb(ibas),nr(ibas), aaa )
           aaa = 1d0/sqrt(aaa)
           rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
           if(nxx==1) cycle
           n1=1
           n2=2
           call gintxx(rprodx(1,n1,l,ibas),rprodx(1,n2,l,ibas),aa(ibas),bb(ibas),nr(ibas), aaa12 )
           rprodx(1:nr(ibas),n2,l,ibas) = rprodx(1:nr(ibas),n2,l,ibas) - aaa12*rprodx(1:nr(ibas),n1,l,ibas)
           n = 2
           call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas),aa(ibas),bb(ibas),nr(ibas), aaa )
           aaa = 1d0/sqrt(aaa)
           rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
        enddo
     enddo
  endif BesselTestBlock
  nbloch    = sum(nblocha)
  nblochpmx = nbloch + ngcmx
  !NOTE for CG coefficienets. <LM3|lm1 lm2> (confusing indexing, use subroutine clebsh)
  ! inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
  ! Injcg = indxcg(inxcg) to indxcg(inxcg)-1
  ! cg(inxcg)  : = <lm3|lm1 lm2>
  ! jcg(lnjcg) : = lm3
  lmxcg = lxx
  call scg_sizechk(lmxcg,lnjcg,lnxcg) !(lmax,c,cindx,js)
  write(6,*)'scg_sizechk= ',lnjcg,lnxcg
  !      if (lmxcg .le. 6) then;      lnjcg = 6500;   lnxcg = 1300
  !      elseif (lmxcg .le. 8) then;  lnjcg = 22700; lnxcg = 3400
  !      elseif (lmxcg .le. 10) then; lnjcg = 62200;lnxcg = 7400
  !      else; call rxi('setcg: cannot handle lmxcg=',lmxcg);   endif
  allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
  call scg(lmxcg,cg,indxcg,jcg)
  if(allochk) write(6,*)' end of scg: cg coefficients generated.'
  !! Apr2021 gomi. we still have a problem of posititve definiteness of the Coulomb matrix. Probably because of this routine   ! These default values ok?
  awald0 = 2d0   ! See p_lat_0
  tol    = 1d-14 ! It was 1d-9 before aug2019. Sakakibara had a problem when he treat slab model of KF for 18 atoms per cell. --> The lowest eigenvalue of coulomb matrix becomes negative.
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
  write(6,"(/' lattc:  as=',f6.3,'   tol=',1p,e9.2,'   lmax=',i2,'   awald=',0p,f7.4,'   v0=',f10.3/' alat1=',f9.5, &
       '   estimates:   nkd',i6,'   nkr',i6)") as,tol,lmax,awald,vol0,alat1,nkdest,nkrest
  call lgen(plat,r0+radd,nkd,nkdmx,dlv,work)
  write(6,"('  r0=',f9.4,'   rc=',f9.4,'   radd=',f9.4,'   nkd=', i7)") r0,r0*alat,radd,nkd
  call lgen(qb0,q0+qadd,nkq,nkqmx,qlv,work)
  write(6,"('  q0=',f9.4,'   qc=',f9.4,'   qadd=',f9.4,'   nkr=', i7)") q0,q0*tpiba,qadd,nkq
  deallocate(work)
  eee=screenfac()    
  if(imode==202) eee=0d0 !! for eps_lmf and epsPP_lmf mode, even the small eee=1d-4 can affect to dielectric function near q=0 when its values is large as one-hundred or more. Thus we set eee=0d0 to avoid this.
  write(6,"(' Coulomb is exp(sqrt(-eee)*r)/r. eee=',d13.6,d13.6)")eee
  ! bessel and hankel for the expansion of exp(-r/r_0)/r. bessel and hankel is renomarized so that its behaves as r^l and r^{-l-1} near r=0.
  !  rkpr means r^l*r for e=0 (r0c =infinity) case
  allocate(rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas),rofi(nrx,nbas))
  do ibas=1,nbas
     call genjh(eee,nr(ibas),aa(ibas),bb(ibas),lx(ibas), nrx,lxx,rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas))
  enddo
  allocate( rojb(nxx, 0:lxx, nbas), sgbb(nxx,  nxx,  0:lxx, nbas))
  do ibas = 1,nbas !--- onsite integrals <j(e=0)|B> and <B|v(onsite)|B>
     call mkjb_4( lxx, lx(ibas),nxx, nx(0:lxx,ibas), aa(ibas),bb(ibas), nr(ibas), nrx,rprodx(1,1,0,ibas), &
          rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas), rojb(1,0,ibas), sgbb(1,1,0,ibas))
  enddo
  nlxx= (lxx+1)**2
  allocate(ngvecc0(3,ngcmx))
  call readqg('QGcou',[0d0,0d0,0d0],  quu,ngc0, ngvecc0) !coulomb matrix for each q = qibz
  deallocate(ngvecc0)
  ngb = nbloch + ngc0
  if(allochk) write(*,*) 'allocate(vcoul)'
  allocate( vcoul(nblochpmx,nblochpmx) )
  vcoul  = 0d0
  write(6,'(a)') " Mix0vec.XXX is not empty only when the corresponding q is in Q0P with zero weight."
  iqxend = nqibz + nq0i+nq0iadd
  open(newunit=ifgb0vec, file="Mix0vec")
  open(newunit=ifgb0vec1,file="Mix0vec1")
  if(imode==202) then; iqxini= nqibz + 1
  else;                iqxini = 1 ;   endif
  if(imode==0) iqxini=1
  write(6,*)'iqxini iqxend=',iqxini,iqxend
  if(abs(sum(qibz(:,1)**2))/=0d0) call rx( 'hvccfp0: sum(q**2)==0d0')
  call MPI__getRange( mpi__iini, mpi__iend, iqxini, iqxend )
  mainforiqx: do 1001 iqx = mpi__iini, mpi__iend ! q=(0,0,0) is omitted!
     write(6,"('#### do 1001 start iqx=',5i5)")iqx,nqibz
     vcoudfile='Vcoud.'//i2char(iqx)  !this is closed at the end of do 1001
     open(newunit=ifvcoud,file=trim(vcoudfile),form='unformatted') !  !! Vcoud file, which contains E(\nu,I), given in PRB81,125102
     if(iqx > nqibz) then  !       iq = 1
        q  = q0i(:,iqx-nqibz)
     else                  !       iq = iqx
        q  = qibz(:,iqx)
     endif
     if(imode==202 .AND. abs(sum(q))<1d-8) cycle
     call readqg('QGcou',q,  quu,ngc, ngvecc ) !qq-->q ! q+G vector
     ngb = nbloch + ngc  !it was ngcnn(iq)
     write(6,'(" iqx q ngc =",i5,3f10.4,i5)') iqx,q,ngc
     allocate( strx(nlxx,nbas,nlxx,nbas)) !! strxq: structure factor.
     do ibas1 =1,nbas
        do ibas2 =1,nbas
           p = bas(:,ibas2)-bas(:,ibas1)
           phasep =exp(img*2*pi*sum(q*p))
           nlx1 = (lx(ibas1)+1)**2
           nlx2 = (lx(ibas2)+1)**2
           if(allochk) write(*,*) 'allocate( s(nlx1,nlx2))'
           allocate( s(nlx1,nlx2),sd(nlx1,nlx2)) !kino add sd----but sd is dummy
           call strxq(1,eee,q,p,nlx1,nlx2,nlx1,alat,voltot,awald,nkd,nkq,dlv,qlv,cg,indxcg,jcg, s,sd)
           strx(1:nlx1,ibas1,1:nlx2,ibas2) = fpi*s      !!! *phasep
           if(allochk) write(*,*)'deallocate( s )'
           deallocate( s,sd )
        enddo
     enddo
     if(allochk) write(*,*)'allocate(rojp,sgpb,fouvb)'
     allocate( rojp(ngc,      nlxx, nbas), sgpb(ngc, nxx, nlxx, nbas), fouvb(ngc, nxx, nlxx, nbas))
     do ibas = 1,nbas !  onsite integrals <j(e=0)|P^(q+G)_L> and <B|v(onsite)|B>
        call mkjp_4(q,ngc, ngvecc, alat, qlat, lxx, lx(ibas),nxx, nx(0:lxx,ibas), bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas), &
             nr(ibas), nrx, rprodx(1,1,0,ibas), eee, rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas), &
             rojp(1,1,ibas),  sgpb(1,1,1,ibas), fouvb(1,1,1,ibas))
     enddo
     if(allochk) write(6,*)' goto vcoulq_4'
     call vcoulq_4(q, nbloch, ngc, nbas, lx,lxx, nx,nxx, alat, qlat, voltot, ngvecc, strx, rojp,rojb, sgbb,sgpb, fouvb, nblochpmx, &
          bas,rmax, eee, aa,bb,nr,nrx,rkpr,rkmr,rofi, vcoul) !the Coulomb matrix
     if(allochk) write(6,*)' end of vcoulq_4'
     deallocate( strx, rojp,sgpb,fouvb)
     write(6,'(" vcoul trwi=",i6,2d22.14)') iqx,sum([(vcoul(i,i),i=1,nbloch)])
     write(6,'("### sum vcoul(1:ngb,      1:ngb) ",2d22.14,2x,d22.14)') sum(vcoul(1:ngb,1:ngb)), sum(abs(vcoul(1:ngb,1:ngb)))
     write(6,'("### sum vcoul(1:nbloch,1:nbloch) ",2d22.14,2x,d22.14)') &
          sum(vcoul(1:nbloch,1:nbloch)),sum(abs(vcoul(1:nbloch,1:nbloch)))
     write(6,*)
1101 continue
     ngbo=ngb
     if(debug) write(6,*) 'write out vcoul' !! == Write out VCCFP ==
     write(6,"(' ngc ngb/ngbo=',6i6)") ngc,ngb,ngbo
     if(allochk) write(*,*) 'allocate( ppovl(ngc,ngc))'
     allocate( oo(ngb,ngb) )
     allocate( ppovl(ngc,ngc) )
     call mkppovl2(alat,plat,qlat, ngc,  ngvecc, ngc,  ngvecc, nbas, rmax, bas, ppovl)
        oo = 0d0
        forall(ipl1=1:nbloch) oo(ipl1,ipl1) = 1d0
        do concurrent(ix=1:ngc,iy=1:ngc)
           oo(nbloch+ix, nbloch+iy) = ppovl(ix,iy)
        enddo
     allocate(oox,source=oo )
     write(6,*)' --- goto eigen check1 --- '
     allocate( vcoul0,source=vcoul(1:ngb,1:ngb) )
     if(allochk) write(*,*) 'allocate(hh(ngb,ngb),oo(ngb,ngb),oox,zz,eb,zzr)'
     allocate(hh(ngb,ngb),zz(ngb,ngb),eb(ngb),zzr(ngb))
     hh  = - vcoul0
     nmx = ngb
     call diagcv(oo,hh,zz,ngb, eb,nmx,1d99,nev) !! diagonalize the Coulomb matrix
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
           write(6,"(a,d13.5)")'KeepPositiveCou enforce : -eb<0 --> eb=0',eb(i)
           eb(i)=0d0
        endif
     enddo
     write(ifvcoud) ngb
     write(ifvcoud) q
     write(ifvcoud) -eb
     write(ifvcoud) zz
     write(6,*)
     write(6,'(" eig0 must be equal to the largest =", 2d24.16)') sum(  dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1))  )
     write(6,'(" zz norm check=",d24.16)')    sum( dconjg(zz(1:ngb,1))*matmul(oox,zz(1:ngb,1)) )
     write(6,*)
     write(6,'(" --- vcoul(exact)=",d14.6," absq2=",d24.16)') fpi*voltot/(sum(tpiba**2*q(1:3)**2)-eee) &
          , (sum(tpiba**2*q(1:3)**2)-eee)
     write(6,'(" --- vcoul(cal ) =",2d14.6)') sum( dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1)) )*voltot
     !          do igc=1,ngb
     !          qqx(1:3) = (q(1:3)+ matmul(qlat, ngvecc(1:3,igc)))
     !          write(6,'(" --- vcoul(exact) xxx =",d14.6," absq2=",d24.16)') fpi*voltot/(sum(tpiba**2*(qqx(1:3)**2)-eee))
     !     &             , (sum(tpiba**2*(qqx(1:3)**2)-eee))
     !          write(6,'(" --- vcoul(cal ) xxx =",2d14.6)') sum( dconjg(zz(1:ngb,igc))*matmul( vcoul0,zz(1:ngb,igc)) )*voltot
     !          enddo
     deallocate(vcoul0)
     if( iqx-nqibz>=1.and.wqt(iqx-nqibz)==0d0) then !! --- To get the vector <Mixed basis| q=0> --------------
        if( .NOT. is_mix0vec()) then     !used original befor oct2006               ! See switch.F ---> this is not used now.
           ifgb0vec_a =ifgb0vec1
           ifgb0vec_b =ifgb0vec
        else              ! ismix0vec=1 is to avoid problem at BZ boundary when is_mix0vec()=0.
           ifgb0vec_a =ifgb0vec
           ifgb0vec_b =ifgb0vec1
        endif
        write(6,*)' voltot=',voltot
        if(ngc/=0) then
           igc0=findloc([(sum(abs(ngvecc(1:3,igc) ))==0,igc=1,ngc)],value=.true.,dim=1)
           write(6,*)' igc0=',igc0,ngvecc(1:3,igc0)
           zzr(nbloch+1:nbloch+ngc) = ppovl(1:ngc,igc0)
        endif
        allocate( gbvec(ngb), b0mat(nbloch) )            !! ... get a vector <Product Basis| q+0>
        call mkb0( q, lxx,lx,nxx,nx, aa,bb,nr,nrx,rprodx, alat,bas,nbas,nbloch, b0mat)
        zzr(1:nbloch) = b0mat(1:nbloch)
        allocate(ooxi,source=oox)
        call matcinv(ngb,ooxi)
        gbvec = matmul(ooxi, zzr)
        deallocate(ooxi)
        dnorm = sqrt( sum(dconjg(gbvec)*zzr) )
        write(ifgb0vec_a,"(3d24.16,2i10,d24.16)") q, ngb,igc0,dnorm
        write(ifgb0vec_a,"(4d24.16)") (gbvec(i),zzr(i),i=1,ngb)
        deallocate( gbvec, b0mat)
        dnorm  = 1d0
        zzr(:) = matmul (oox, zz(:,1))
        igc0 = 999999 !dummy now            ! phasex ---just to clean. this is irrelevant
        phasex =1d0
        do i=1,ngb
           if(abs(zz(i,1)) > 1d-3) phasex = abs(zz(i,1))/zz(i,1)
        enddo
        zz(:,1)= phasex * zz(:,1)
        zzr(:) = phasex * zzr(:)
        write (ifgb0vec_b,"(3d24.16,2i10,d24.16)") q, ngb,igc0,dnorm
        write (ifgb0vec_b,"(4d24.16)") (zz(i,1),zzr(i),i=1,ngb)
     endif
     deallocate(hh,oo,zz,eb,oox,zzr)
     deallocate(ppovl)
     close(ifvcoud)
1001 enddo mainforiqx
  deallocate(ngvecc)
  call cputid(0)
  if(imode==202) call rx0( ' OK! hvccfp0 imode=202 only for Q0P')
  if(imode==0) call rx0( ' OK! hvccfp0 imode=0')
  if(imode==3) call rx0( ' OK! hvccfp0 imode=3')
END PROGRAM hvccfp0
