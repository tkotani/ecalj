      program hvccfp0
c- Coulomb matrix. <f_i | v| f_j>_q. ------------------------
c input files
c    HVCCIN      : some inputs by hbg0.
c    PLN         : plane wave expansion data by nbg0.
c    BASFP//atom : product basis by hbasfp0. ic=iatom should be kept!
c output
c    VCCFP : the coulomb matrix vcoul(nblochpmx,nblochpmx) for all qibz.
c-------------------------------------------------------------
c int
c    strx: structure constant for e=0 (means 1/|r-r'| )
c
      use m_readqg
      use keyvalue
      implicit none
      integer(4) :: ifvcfpout,iopen,ifhvccfp,is, nqbz, nbas,lmxcg,
     &  nband, ifplane,ngpmx, ngcmx, nblochpmx, nbloch,
     &  ibas,ic,lxx,nxx,nrx,l,n,k,isx,kdummy, iclose,
     &  nkdmx,nkqmx,lmax,nkdest,nkrest,ngp,ngc,nlxx,i,lnjcg,lnxcg,
     &  nkd,nkq ,ibas1,ibas2,nlx1,nlx2, nqibz,iqibz
      real(8) :: alat, plat(3,3),qlat(3,3),q(3),p(3),voltot,
     &  pi,fpi,tripl,alat0,epsx,
     &  tol,as,tpiba,qb0(3,3),vol0,rdist0,qdist0,radd,qadd,
     &  a0,awald,alat1,tol1,r0,q0,awald0,qg(3),   absqg2,aaa,aaa12
      integer(4),allocatable :: jcg(:),indxcg(:),
     &  lx(:),kmx(:),nblocha(:),nr(:),ificrb(:),
     &  nx(:,:),ngvecp(:,:),ngvecc(:,:),ngvecci(:,:,:),iqibzx(:)
      real(8),allocatable :: qbz(:,:),qibz(:,:),bas(:,:),rmax(:),
     &  cg(:),rprodx(:,:,:,:),dlv(:,:),qlv(:,:),work(:),ngcn(:),
     &  rojb(:,:,:), sgbb(:,:,:,:),aa(:),bb(:),rofit(:),phi(:),psi(:),
     &  wqt(:), q0i(:,:)
      complex(8) ,allocatable :: vcoul(:,:),geig(:,:),strx(:,:,:,:),
     &  sgpb(:,:,:,:),sgpp(:,:,:,:),
     &  fouvb(:,:,:,:),fouvp(:,:,:,:),vcoul0(:,:),
     &  s(:,:),sd(:,:),rojp(:,:,:) , vcoulnn(:,:)
      character*7,allocatable :: filename(:)
      character(20) :: xxt

c dimension work space
c      integer(4),parameter :: nsize= 1 000 000
c      real w(nsize)
c      common /w/ w

      complex(8):: phasep,img=(0d0,1d0)
      integer(4)::ir,ig1,n1,n2
c check
      complex(8),allocatable :: hh(:,:),oox(:,:),ooxi(:,:),oo(:,:),zz(:,:),zzr(:)
      real(8),allocatable    :: eb(:)

      complex(8),allocatable :: matp(:),matp2(:)
      complex(8) :: xxx,trwv
      integer(4) :: ngb,nev,nmx,iqx,ipl1,ipl2,nq0i,igx1,igx2
      logical checkeig
      logical:: besseltest=.false. !test
      real(8) :: sss1,sss2,dnorm
c

      complex(8),allocatable:: gbvec(:), ppovl(:,:), b0mat(:)

      integer(4) ::igc,igc0,ifgb0vec,ifgb0vec1,ix, iy

      integer(4) :: iqxini, iqxend,imode
      logical :: allochk=.false. !paralellx0=.true., 

      complex(8),allocatable:: hh1(:,:),oo1(:,:)
      integer(4):: bzcase, nqnumc,ifiqgc
c      character(5)  :: charnum5
c      integer(4),allocatable:: iqok(:)
      real(8):: qqq(3),QpGcut_Cou       ,qq(3)

      integer(4),allocatable:: ngvecc0(:,:)
      integer(4):: ngc0

      real(8):: ginv(3,3),quu(3),det

c---
      real(8),allocatable :: rkpr(:,:,:),rkmr(:,:,:),rofi(:,:)
      real(8):: eee,eees, q_org(3),screenfac
      integer(4):: ifvcfporg,nqbz_in,nblochpmx_in
      complex(8),allocatable:: vcoul_org(:,:)

      logical :: smbasis,debug=.false.,smbb
      integer(4)   :: ifprodmt,nl_r,lx_,nxx_r,nxdim,ibl1,nn,no,ngbnew,
     &  nmatch,ifpmatch,nmatch_q,ifpmatch_q,m,ifpomat,nbln,ibln,ngb_in,nnr,igc2
      character(3) :: charnum3
      character(5) :: charnum5
      character(11):: filenamep
      integer(4),allocatable:: nx_r(:), ibl(:,:,:,:)
     &  ,imatcho(:),imatchn(:),imatcho_q(:),imatchn_q(:)
      real(8),allocatable:: prodmt(:,:,:,:),rdmatch(:,:,:,:)
      complex(8),allocatable:: ppmt(:,:,:,:),pmat(:,:),pomat(:,:),oon(:,:)
      complex(8):: pval,pslo,phasex
      real(8)::absqq,qqx(3), epsmx,aaaa
      integer(4):: nnmx ,ngcnn,ngbo
cki      integer(4):: is_mix0vec ,ifgb0vec_a,ifgb0vec_b
      integer(4):: ifgb0vec_a,ifgb0vec_b , ifvcoud,idummy
      logical:: is_mix0vec,newaniso,wvcc
      character(128):: vcoudfile
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi

      write(6,"(' mode(now only =0)? +(pararell iqxini iqxend ?)')")
c      call readin5(imode,iqxini,iqxend)
      read(5,*) imode
      call headver('hvccfp0: start',imode)
      call cputid (0)
      if(imode==202 ) then
        write(6,*)' hvccfp0: imode=',imode
c      elseif(imode==101) then
c        write(6,*)' hvccfp0: imode=',imode
c        write(6,*)' remove_r0c is effective'
c        write(6,*)' Generate VCCFP = VCCFP.ORG - new_VCCFP'
c        ifvcfporg = iopen( "VCCFP.ORG",0,-1,0)
c      elseif(imode==102) then
c        write(6,*)' hvccfp0: imode=',imode
c        write(6,*)' remove_r0c is effective'
      elseif(imode==0) then
      else
        stop 'hvccfp0: now hvccfp0 support just normal mode=0 202 101'
      endif
c      if(iqxini< 2) paralellx0=.false.
c      if(paralellx0) then
c        write(6,"(' PARALELL.X0 mode: iqxini iqxend=',2i3)")
c     &  iqxini, iqxend
c      endif

C --- q, nqbz, alat, qlat, nbas, bas
      ifhvccfp = iopen('HVCCIN',0,-1,0)
      read (ifhvccfp) alat, plat,qlat, nqbz, nbas, nband
      if(allochk)
     & write(*,*) 'allocate(qbz(3,nqbz),bas(3,nbas),rmax(nbas))'
      allocate(qbz(3,nqbz),bas(3,nbas),rmax(nbas))
      read(ifhvccfp) qbz, bas,rmax
      read(ifhvccfp) nqibz
      if(allochk) 
     &  write(*,*)'allocate(qibz(3,nqibz),iqibzx(nqbz))'
      allocate(qibz(3,nqibz),iqibzx(nqbz))
      read(ifhvccfp) qibz(1:3,1:nqibz)
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      print *,' voltot=',voltot
      write(6,*)
      write(6,"(i4,3f13.6)")(i,qibz(1:3,i),i=1,nqibz)
c      print *," xxx=",alat, plat,qlat, nqbz, nbas, nband
c      print *," xxx=",qbz,bas,rmax
c      print *," xxx=",rmax
c      is = iclose ('HVCCIN')

c$$$c-----------------------------------------------------
c$$$      iqibzx =0
c$$$      do iqibz=1,nqibz
c$$$        do iq =1,nqbz
c$$$          if( sum(abs(qibz(:,iqibz)-qbz(:,iq)))<1d-8 ) then
c$$$            iqibzx(iq)=iqibz !iqibzx is the index for iqx in bz
c$$$            goto 119
c$$$          endif
c$$$        enddo
c$$$        stop " hvccfp: cannot find iqx"
c$$$  119 enddo
c$$$c --- Readin PLN. plane wave contributions 2000 May
c$$$      ifplane = iopen('PLN',0,-1,0)
c$$$      read (ifplane) ngpmx, ngcmx
c$$$c q+G vector
c$$$      if(allochk)
c$$$     & write(*,*) 'allocate( ngcn(nqbz), ngvecci(3,ngcmx,nqibz))'
c$$$      allocate( ngcn(nqbz), ngvecci(3,ngcmx,nqibz))
c$$$      do iq=1, nqbz
c$$$        read(ifplane) ngp, ngc
c$$$        if(allochk)write(*,*)'allocate( geig, ngvecp, ngvecc)'
c$$$        allocate( geig(ngp,nband), ngvecp(3,ngp), ngvecc(3,ngc))
c$$$        read(ifplane) ngvecp, ngvecc, geig
c$$$        if(iqibzx(iq) /=0) then
c$$$          iqibz = iqibzx(iq)
c$$$          ngcn(iqibz) = ngc
c$$$          ngvecci(1:3,1:ngc,iqibz) = ngvecc(1:3,1:ngc)
c$$$        endif
c$$$        if(allochk)write(*,*) 'deallocate( geig, ngvecp, ngvecc)'
c$$$        deallocate( geig, ngvecp, ngvecc)
c$$$      enddo
c-----------------------------------------------------

c q+G vector
c$$$      ifiqgc = 1302
c$$$      open(ifiqgc, file='QGcou',form='unformatted')
c$$$      read(ifiqgc ) nqnumc, ngcmx, QpGcut_Cou
c$$$      allocate( ngcn(nqbz), ngvecci(3,ngcmx,nqibz), ngvecc(3,ngcmx),iqok(nqibz))
c$$$      iqok=1
c$$$      do iq=1, nqnumc
c$$$        read (ifiqgc) qqq, ngc
c$$$        read (ifiqgc) ngvecc(1:3,1:ngc)
c$$$        do iqibz=1,nqibz
c$$$          if( sum(abs(qibz(:,iqibz)-qqq))<1d-8 ) then
c$$$            ngcn(iqibz) = ngc
c$$$            ngvecci(1:3,1:ngc,iqibz) = ngvecc(1:3,1:ngc)
c$$$            iqok(iqibz)=0
c$$$            exit
c$$$          endif
c$$$        enddo
c$$$        if(sum(iqok)==0) exit
c$$$      enddo
c$$$      if(sum(iqok)/=0) stop 'hvccfp0: iqok/=0;wrong QGcou?'
c$$$      deallocate(ngvecc,iqok)
c$$$      close(ifiqgc)


      call readngmx('QGcou',ngcmx)
      allocate(ngvecc(3,ngcmx))
c        allocate( ngcn(nqibz), ngvecci(3,ngcmx,nqibz),ngveccc0(3,ngcmx))
c        do iqibz = 1,nqibz
c          call readqg('QGcou',qibz(:,iqibz), ngcn(iqibz),ngvecci(1,1,iqibz))
c        enddo
c        call readqg('QGcou',(/0d0,0d0,0d0/), ngc0, ngvecc0(1,1))
c        call releaseqg('QGcou')



c --- Readin BASFP//atom. The product basis functions.
      if(allochk)
     &  write(*,*)'allocte(lx,kmx,nblocha,nr,aa,bb,filename,ificrb'
      allocate(lx(nbas),kmx(nbas),nblocha(nbas),
     &         nr(nbas),aa(nbas),bb(nbas),filename(nbas),
     &         ificrb(nbas) )

      do ibas = 1,nbas
        ic = ibas !
        filename(ibas)= 'BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10))
        ificrb(ibas)  = iopen ( filename(ibas),1,3,0)
        read(ificrb(ibas),"(4i6,2d24.16)") 
     &    lx(ibas), kmx(ibas), nblocha(ibas), nr(ibas),aa(ibas),bb(ibas)
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
cccccccccccccc
c      write(660+ibas,*)
c      write(660+ibas,'(" *** nlibas=",3i3)')  n,l,ibas
c      do i=1,nr(ibas)
c       write(660+ibas,'(2d16.8)') bb(ibas)*( exp(aa(ibas)*(i-1))- 1d0),
c     & rprodx(i,n,l,ibas)
c      enddo
cccccccccccccc
          enddo
        enddo
c      isx = iclose(filename(ibas))
      enddo


cccccccccccccccccccccccccccccccccccccccccccccc
cccc TEST cccccccccccccccc
c      open(117, file='xin')
c      do i=1,nr(1)
c        read(117,"(d24.16)") rprodx(i,1,0,1)
c      enddo
cccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccc
c TEST cccccccccccccccccccccccccccccccccccccc
      if(besseltest) then
        print *
        print *
        print *, ' *** TEST case ***  rprodx is given by Bessel.'
ccc test G, corresponding <q+G|v|q+G> should be exact. e.g. ig1=1 and  ig1=35 for iqx=2
ccc You can change these values for tests. ccccccccccccc
        iqx  = 2
        igx1 = 1
        igx2 = 35
c
        write(6,"(' iqx=',i3,' ig1 ig2=',2i3)") iqx,igx1,igx2
        write(6,"(a)")
     & ' <q+G|v|q+G> for the corresponding iqx ig1 ig2 should be exact!'
        write(6,"(a)") ' See fort.196'
        write(6,"(a)")
     & ' Errors will be from the radial function integrals !!!'
        write(6,"(a)") ' You can slso so similar test from hbasfp0.'
        write(6,"(a)") ' See test1 in basnfp0.'
c
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
        allocate(rprodx(nrx,nxx,0:lxx,nbas),rofit(nrx)
     &  ,phi(0:lxx),psi(0:lxx))
        rofit(1) = 0d0
        do ir   = 1, nrx
          rofit(ir) = bb(1)*( exp(aa(1)*(ir-1)) - 1d0)
        enddo
        do n = 1, nxx
          if(n==1) ig1 = igx1
          if(n==2) ig1 = igx2
          qg(1:3) = 
     &    tpiba * (qibz(1:3,iqx)+ matmul(qlat, ngvecci(1:3,ig1,iqx)))
          absqg2  = sum(qg(1:3)**2)
c
          do ir =1,nrx
            call bessl(absqg2*rofit(ir)**2,lxx,phi,psi)
            do ibas=1,nbas
              do l = 0, lx(ibas)
                rprodx(ir,n,l,ibas) = phi(l)* rofit(ir) **(l +1 )
              enddo
            enddo
          enddo
        enddo
c --- orthogonalized rprodx.
        do ibas=1,nbas
          do l = 0, lx(ibas)
            rprodx(1:nr(ibas),1,l,ibas)=
     &     rprodx(1:nr(ibas),1,l,ibas)
     &   + rprodx(1:nr(ibas),2,l,ibas)
            n = 1
            call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas)
     &    ,aa(ibas),bb(ibas),nr(ibas), aaa )
            aaa = 1d0/sqrt(aaa)
            rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
            if(nxx==1) cycle
            n1=1
            n2=2
            call gintxx(rprodx(1,n1,l,ibas),rprodx(1,n2,l,ibas)
     &    ,aa(ibas),bb(ibas),nr(ibas), aaa12 )
            rprodx(1:nr(ibas),n2,l,ibas) = rprodx(1:nr(ibas),n2,l,ibas)
     &      - aaa12*rprodx(1:nr(ibas),n1,l,ibas)
            n = 2
            call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas)
     &    ,aa(ibas),bb(ibas),nr(ibas), aaa )
            aaa = 1d0/sqrt(aaa)
            rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
          enddo
        enddo
      endif
cccc TEST end ccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccc


      nbloch    = sum(nblocha)
      nblochpmx = nbloch + ngcmx

c --- CG coefficienets. <LM3|lm1 lm2>
c inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
c Injcg = indxcg(inxcg) to indxcg(inxcg)-1
c cg(inxcg)  : = <lm3|lm1 lm2>
c jcg(lnjcg) : = lm3
      lmxcg = lxx

      call scg_sizechk(lmxcg,lnjcg,lnxcg) !(lmax,c,cindx,js)
      write(6,*)'scg_sizechk= ',lnjcg,lnxcg
c      if (lmxcg .le. 6) then
c        lnjcg = 6500
c        lnxcg = 1300
c      else if (lmxcg .le. 8) then
c        lnjcg = 22700
c        lnxcg = 3400
c      else if (lmxcg .le. 10) then
c        lnjcg = 62200
c        lnxcg = 7400
c      else
c        call rxi('setcg: cannot handle lmxcg=',lmxcg)
c      endif
c      if(allochk)
c     &  write(*,*) 'allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))'
      allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
      call scg(lmxcg,cg,indxcg,jcg)
      if(allochk) print *,' end of scg: cg coefficients generated.'


      call dinv33(qlat,0,ginv,det)

c --- Get real-space vectors and reciprocal-space vectors for Ewald sum.
C defaults values for ewald sum
c      call lattc(awald0,tol,alat,alat,plat0,gx,gy,gz,gam,plat,qlat,
c     .   lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx,w(owork))
c- taken from lattc.f

c default values ok?
      awald0 = 2d0   !See p_lat_0
      tol    = 1d-9
      nkdmx  = 800
      nkqmx  = 800
      lmax   = 2*lxx  !lxx or lmax=6 ???

      vol0= abs(tripl(plat,plat(1,2),plat(1,3)))
      as   = awald0
c      alat0= alat
      alat1= alat
c      if(alat1.le.0.5d0) alat1=alat
      tpiba=8.d0*datan(1.d0)/alat
      call cross_x(plat(1,2),plat(1,3),qb0)
      call cross_x(plat(1,3),plat(1,1),qb0(1,2))
      call cross_x(plat(1,1),plat(1,2),qb0(1,3))
      qb0(1:3,1:3) = qb0(1:3,1:3)/vol0

      rdist0=vol0**(1.d0/3.d0)
      qdist0=1.d0/rdist0
      radd=.7*rdist0
      qadd=.7*qdist0
      a0=as/rdist0
      awald=a0/alat
cccccccccccccccccccccccccccccccccccc
c takao
c      tol1= tol*alat**(lmax+1) *0.01
cccccccccccccccccccccccccccccccccccccc
      tol1= tol*alat**(lmax+1)
      if(allochk) write(*,*) 'allocate(dlv, qlv, work) '
      allocate(dlv(3,nkdmx), qlv(3,nkqmx), work(max0(nkdmx,nkqmx)) )
      call lctoff(a0,vol0,lmax,tol1,r0,q0)
      nkdest =4.18879*(r0+radd)**3/vol0+.5
      nkrest =4.18879*(q0+qadd)**3*vol0+.5
      write(6,340) as,tol,lmax,awald,vol0,alat1,nkdest,nkrest
  340 format(/' lattc:  as=',f6.3,'   tol=',1p,e8.2,'   lmax=',i2,
     .  '   awald=',0p,f7.4,'   v0=',f10.3/' alat1=',f9.5,
     .  '   estimates:   nkd',i6,'   nkr',i6)
      call lgen(plat,r0+radd,nkd,nkdmx,dlv,work)
      write(6,342) r0,r0*alat,radd,nkd
  342 format('  r0=',f9.4,'   rc=',f9.4,'   radd=',f9.4,'   nkd=', i7)
      call lgen(qb0,q0+qadd,nkq,nkqmx,qlv,work)
      write(6,341) q0,q0*tpiba,qadd,nkq
  341 format('  q0=',f9.4,'   qc=',f9.4,'   qadd=',f9.4,'   nkr=', i7)
      if(allochk) write(*,*) 'deallocate(work)'
      deallocate(work)

C... readin r0c
      if(newaniso()) then
        eee=screenfac() !takao feb2012
c      elseif(imode==101.or.imode==102) then
c        eee = eees()
      else
        eee=0d0
      endif

!! for eps_lmf and epsPP_lmf mode, 
!! even the small eee=1d-4 can affect to dielectric function near q=0 when its values is large as one-hundred or more.
!! Thus we set eee=0d0 to avoid this.
      if(imode==202) then !
         eee=0d0
      endif

      write(6,"(' Coulomb is exp(sqrt(-eee)*r)/r. eee=',d13.6,d13.6)")eee

C--- bessel and hankel for the expansion of exp(-r/r_0)/r.
c bessel and hankel is renomarized so that its behaves as r^l and r^{-l-1} near r=0.
c  rkpr means r^l*r for e=0 (r0c =infinity) case
      allocate(rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas),rofi(nrx,nbas))
      do ibas=1,nbas
        call genjh(eee,nr(ibas),aa(ibas),bb(ibas),lx(ibas), nrx,lxx,
     o     rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas))
      enddo

C--- onsite integrals <j(e=0)|B> and <B|v(onsite)|B>
cc      if(allochkw) write(*,*) ' allocate  rojb, sgbb '
c      allocate( rojb(nxx, 0:lxx, nbas), sgbb(nxx,  nxx,  0:lxx, nbas))
c      do ibas = 1,nbas
c        call mkjb( lxx, lx(ibas),nxx, nx(0:lxx,ibas),
c     i                 aa(ibas),bb(ibas), nr(ibas), nrx,
c     i                 rprodx(1,1,0,ibas),
c     o        rojb(1,0,ibas), sgbb(1,1,0,ibas))
c      enddo
      allocate( rojb(nxx, 0:lxx, nbas), sgbb(nxx,  nxx,  0:lxx, nbas))
      do ibas = 1,nbas
        call mkjb_4( lxx, lx(ibas),nxx, nx(0:lxx,ibas),  
     i                 aa(ibas),bb(ibas), nr(ibas), nrx,
     i                 rprodx(1,1,0,ibas),
     i     rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas),
     o        rojb(1,0,ibas), sgbb(1,1,0,ibas))
      enddo

c----------------
C--- coulomb matrix for each q = qibz
c----------------
      nlxx= (lxx+1)**2
c      ngb = nbloch + ngcn(1)
      allocate(ngvecc0(3,ngcmx))
      call readqg('QGcou',(/0d0,0d0,0d0/),ginv,  quu,ngc0, ngvecc0)
      deallocate(ngvecc0)
      ngb = nbloch + ngc0
      if(allochk) write(*,*) 'allocate( vcoul)'
      allocate( vcoul(nblochpmx,nblochpmx) )
c      if(imode==101) allocate( vcoul_org(nblochpmx,nblochpmx) )

      vcoul  = 0d0

C... q near zero
      write(6,*) '--- readin Q0P -------'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      if(allochk) write(*,*)'allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )'
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      write (6,"(d13.5,3x, 3d13.5)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      close(101)

      write(6,*) ' *** goto do iq nqibz nq0i=',nqibz,nq0i

C --- Check PARALELL.X0
c      INQUIRE (FILE = 'PARALELL.X0', EXIST = paralellx0)
      wvcc=.true.
      if(newaniso()) wvcc=.false.
      write(6,'(a)') " Mix0vec.XXX is not empty only when"
     &   //" the corresponding q is in Q0P with zero weight."
c      if(paralellx0) then
c        if(wvcc) ifvcfpout = iopen( "VCCFP." //xxt(iqxini,iqxend),0,-1,0)
c        ifgb0vec = iopen ( "Mix0vec."//xxt(iqxini,iqxend),1,3,0)
c        ifgb0vec1 = iopen ( "Mix0vec1."//xxt(iqxini,iqxend),1,3,0)
c      else
        iqxend = nqibz + nq0i
        if(wvcc) ifvcfpout = iopen('VCCFP',0,-1,0)
        ifgb0vec = iopen ( "Mix0vec",1,3,0)
        ifgb0vec1 = iopen ( "Mix0vec1",1,3,0)
c      endif

      if(imode==202) then
        iqxini= nqibz + 1
c      elseif(paralellx0) then
c     &  !skip
      elseif(bzcase()==1) then
!        iqxini = 2
        iqxini = 1 !oct2005
      else
        iqxini = 1
      endif
!!
      if(newaniso().and.imode==0) then
        iqxini=1
        iqxend=nqibz
      endif
      print *,'iqxini iqxend=',iqxini,iqxend
c qibz loop
c      epsx =  0.01d0
      if(bzcase()==1) then
        if(abs(sum(qibz(:,1)**2))/=0d0) stop 'hvccfp0: sum(q**2)==0d0'
      endif
      if(wvcc) write(ifvcfpout) nqbz, nblochpmx
c      if(imode==101) then
c        read(ifvcfporg) nqbz_in, nblochpmx_in
c        if(nqbz /= nqbz_in) stop 'nqbz /= nqbz_in VCCFP.ORG'
c        if(nblochpmx /= nblochpmx_in) 
c     &     stop 'nblochpmx /= nblochpmx_in VCCFP.ORG'
c      endif



C... Readin PRODMT into prodmt. oct2005
      smbb = smbasis()
      write(6,*) ' smooth mixed basis=',smbb
      if(smbasis()) then
        allocate( prodmt(2,nxx,0:lxx,nbas))
        allocate( nx_r(0:lxx))
        do ibas =1,nbas
          filenamep = 'PRODMT_'//charnum3(ibas)
          ifprodmt  = iopen(filenamep,0,-1,0)
          read(ifprodmt) nl_r
          if( 2*(nl_r-1) /= lxx ) then
            write(6,*) 2*(nl_r-1),lxx
            stop '2*nl_r-1 /= lxx '
          endif
          read(ifprodmt) nxx_r
          write(6,"(' nxx =',100i3)")nxx_r
          if(nxx_r>nxx) stop 'nxx_r>nxx'
          read(ifprodmt) nx_r(0:lxx)
          write(6,"(' nx_r=',100i3)") nx_r(0:lxx)
          lx_ = lx(ibas)
          if(sum(abs(nx(0:lx_,ibas)-nx_r(0:lx_))) /=0) then
            write(6,*)' debug: nx  =',nx(0:lx_,ibas)
            write(6,*)' debug: nx_r=',nx_r(0:lx_)
            stop 'nx /=nx_r'
          endif
          read(ifprodmt) prodmt(1:2, 1:nxx_r, 0:lxx, ibas)
          write(6,*)' sumcheck prodmt=',sum(abs(prodmt(:,:,:,ibas)))
          isx = iclose(filenamep)
        enddo


ccccccccccccccccccccccccccccccccccccccc
C... Check write for radial part of the product basis
        if(.false.) then
          do ibas= 1,1 !1,nbas
            do l   =  0,lx(ibas)
              open(1011,file='ProdOld_ibas'//charnum3(ibas)//'_l'//charnum3(l))
c         open(2011,file='ProdNew_ibas'//charnum3(ibas)//'_l'//charnum3(l))
              nxdim = nx(l,ibas)
              do ix=1,nxdim
                write(1011,"(' -- -- -- ',3i3,' --- ' )") ix,l,ibas
c         write(2011,"(' -- -- -- ',3i3,' --- ' )") ix,l,ibas
                do ir =1,nr(ibas)
                  write(1011,"(d13.5,2x,2d18.8)")   
     &   rofi(ir,ibas), rprodx(ir,ix,l,ibas) 
     &   , rprodx(ir,ix,l,ibas) /rofi(ir,ibas)
c         write(2011,"(d13.5,2x,2d18.8)")
c     &   rofi(ir,ibas), sum(rprodx(ir,1:nxdim,l,ibas)*rdmatch(1:nxdim,ix,l,ibas))
c     &   , sum(rprodx(ir,1:nxdim,l,ibas)*rdmatch(1:nxdim,ix,l,ibas))/rofi(ir,ibas)
                enddo
              enddo
              close(1011)
c         close(2011)
            enddo
          enddo
c       stop 'text end'
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        allocate( rdmatch(nxx,nxx,0:lxx,nbas) )
        do ibas= 1, nbas
          do l   = 0, lx(ibas)
            nxdim = nx(l,ibas)
            if(nxdim<=1)write(6,*)'hvccfp0:smbasis case error nxdim <=1'
!       pval  = prodmt(1, 1:nxdim, l,ibas)
!       pslo  = prodmt(2, 1:nxdim, l,ibas)
!       prod(r, inew) = \sum_iold rrmat(inew,iold) * prod(r,iold)
            write(6,"('goto mkradmatch ibas lnxdim =',3i4)")ibas,l,nxdim
            call mkradmatch(prodmt(1:2, 1:nxdim, l,ibas), nxdim,  
     o          rdmatch(1:nxdim,1:nxdim,l,ibas) )
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
          stop ' hvccfp0:smbasis mode  error ibl1/= nbloch'
        endif
! index (mx,nx,lx,ibas) ordering
ctttt
        nnr = 2 ! =2 new
        ! =0 equivalence with original mixed basis
        write(6,*)' sss:nbas lx=',nbas,lx(1:nbas)

        nbln=0
        do ibas= 1, nbas
          do l   = 0, lx(ibas)
            write(6,"('sss: nx=',3i4)") ibas,l,nx(l,ibas)
            if(nx(l,ibas)<=0) cycle
            if(nx(l,ibas)==1) stop 'nx(l,ibas) =1'
ccccccccccccccccd
ctttt
c         nnr = 2 ! =2 new
c               ! =0 equivalence with original mixed basis
c         if(l<=3) nnr=0
cccccccccccccccccc
            nbln = nbln + (2*l+1)*(nx(l,ibas)-nnr)
          enddo
        enddo
        allocate( pmat(nbloch+ngcmx, nbln+ngcmx) )
        pmat=0d0
        ibln = 0
        do ibas= 1, nbas
          do l   = 0, lx(ibas)
ccccccccccccccccccc
ctttt
c         nnr = 2 ! =2 new
c               ! =0 equivalence with original mixed basis
c         if(l<=3) nnr=0
cccccccccccccccccc
            do nn  = nnr+1, nx(l,ibas) !nn=1 and nn=2 corresponds to non-zero val sol
              do m   = -l, l
                ibln = ibln +1
                nxdim = nx(l,ibas)
                pmat( ibl(m,1:nxdim,l,ibas), ibln)
     &    =  rdmatch(1:nxdim, nn, l,ibas)
ctttt
c          pmat( ibl(m,nn,l,ibas), ibln)
c     &    =  1d0
ccccccccccccccccccccccccccccccccccccccccccccccccccc
              enddo
            enddo
          enddo
        enddo
C... Store matting matrix (imatchn,imatcho,pmatch)
        ifpomat   = iopen('POmat',0,-1,0)
c       write(6,*)'ttt= sumchk pmat(b)=',sum(abs(pmat(1:nbloch, 1:nbln)))
      endif
!! === open file Vcoud ===
!! This contains E(\nu,I), given in PRB81,125102

!! == main loop for iqx ==
      do 1001 iqx = iqxini, iqxend ! q in IBZ. avoid q=0 case for iqx=1
        write(*,"('#### do 1001 start iqx=',i5)")iqx
        vcoudfile='Vcoud.'//charnum5(iqx)  !this is closed at the end of do 1001
        ifvcoud = iopen(trim(vcoudfile),0,-1,0)
        if(iqx > nqibz) then  !       iq = 1
          q  = q0i(:,iqx-nqibz)
          qq = 0d0
        else                  !       iq = iqx
          q  = qibz(:,iqx)
          qq = q
        endif
cccccccccccccccccccccccc
c        if(imode==202) then !for iqx>nqibz
c           qq=q
c        endif
cccccccccccccccccccccccc

        if(.not. newaniso() ) then !this is for fe_epsPP_lmfh_chipm feb2012
          if(sum(q**2)<1d-12) q=(/1d-4,0d0,0d0/) !takao oct2006
        endif
!! ==== q+G vector ====
        call readqg('QGcou',qq,ginv,  quu,ngc, ngvecc )
        ngb = nbloch + ngc  !it was ngcnn(iq)
        write(6,'(" iqx q ngc =",i5,3f10.4,i5)') iqx,q,ngc

        if(newaniso()) then
          continue
        elseif(bzcase()==1.and.iqx==1) then
          goto 1101
        endif

c        ngc = ngcn(iq)
c        ngvecc(1:3,1:ngc) = ngvecci(1:3,1:ngc,iq)
c        write(6,*)' iq ngc=',iq, ngc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  q test
c       q=(/ 0.09d0,0.09d0,0.09d0/)
c      q = q+(/ 0.01d0,0.01d0,0.01d0/)
c        q=q/4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C--- strxq structure factor.
        if(allochk) print *,' goto strxq'
        if(allochk) write(*,*) 'allocate( strx(nlxx,nbas,nlxx,nbas))'
        allocate( strx(nlxx,nbas,nlxx,nbas))
        do ibas1 =1,nbas
          do ibas2 =1,nbas
            p = bas(:,ibas2)-bas(:,ibas1)
            phasep =exp(img*2*pi*sum(q*p))
            nlx1 = (lx(ibas1)+1)**2
            nlx2 = (lx(ibas2)+1)**2
            if(allochk) write(*,*) 'allocate( s(nlx1,nlx2))'
            allocate( s(nlx1,nlx2),sd(nlx1,nlx2)) !kino add sd----but sd is dummy
c          call strxq(1,0d0,q,p,nlx1,nlx2,nlx1,alat,voltot,
            call strxq(1,eee,q,p,nlx1,nlx2,nlx1,alat,voltot,
     i        awald,nkd,nkq,dlv,qlv, 
     i        cg,indxcg,jcg,  
     o        s,sd)
            strx(1:nlx1,ibas1,1:nlx2,ibas2) = fpi*s      !!! *phasep
            if(allochk) write(*,*)'deallocate( s )'
            deallocate( s,sd )
          enddo
        enddo
ccccccccccccccccccccccccccc
c        strx=0d0
cccccccccccccccccccccccccccc

C--- onsite integrals <j(e=0)|P^(q+G)_L> and <B|v(onsite)|B>
c$$$      if(.true.) then !==New version without sgpp and fouvp allocation June2004=====
        if(allochk) write(*,*)'allocate(rojp,sgpb,fouvb)'
        allocate(  rojp(ngc,      nlxx, nbas),
     &            sgpb(ngc, nxx, nlxx, nbas),  
     &           fouvb(ngc, nxx, nlxx, nbas))
c     &            sgpp(ngc, ngc, nlxx, nbas),
c     &           fouvp(ngc, ngc, nlxx, nbas) )
        do ibas = 1,nbas
          if(allochk) print *,' --- goto mkjp_4',ibas
          call mkjp_4(q,ngc, ngvecc, alat, qlat, 
     i      lxx, lx(ibas),nxx, nx(0:lxx,ibas),  
     i      bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas), 
     i      nr(ibas), nrx, rprodx(1,1,0,ibas),
     i   eee, rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas),
     o      rojp(1,1,ibas),  sgpb(1,1,1,ibas),
     o      fouvb(1,1,1,ibas))
c          call mkjp3(q,ngc, ngvecc, alat, qlat,
c     i      lxx, lx(ibas),nxx, nx(0:lxx,ibas),
c     i      bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas),
c     i      nr(ibas), nrx, rprodx(1,1,0,ibas),
c     o      rojp(1,1,ibas),  sgpb(1,1,1,ibas),
c     o      fouvb(1,1,1,ibas))
        enddo

C--- the Coulomb matrix
        if(allochk) print *,' goto vcoulq_4'
        call vcoulq_4(q, nbloch, ngc,  
     i                  nbas, lx,lxx, nx,nxx,
     i                  alat, qlat, voltot, ngvecc,
     i       strx, rojp,rojb, sgbb,sgpb, fouvb, !sgpp,fouvp, 
     i       nblochpmx, bas,rmax, 
     i    eee, aa,bb,nr,nrx,rkpr,rkmr,rofi,
     o          vcoul)
c       call vcoulq2(q, nbloch, ngc,
c     i                  nbas, lx,lxx, nx,nxx,
c     i                  alat, qlat, voltot, ngvecc,
c     i       strx, rojp,rojb, sgbb,sgpb, fouvb, !sgpp,fouvp,
c     i       nblochpmx, bas,rmax,
c     o          vcoul)
        if(allochk) print *,' end of vcoulq_4'
        deallocate( strx, rojp,sgpb,fouvb)

c$$$      else !===old version (allocation of sgpp and fouvp are required) ====
c$$$
c$$$       if(allochk) write(*,*) 'allocate(rojp,sgpb,sgpp,fouvb,fouvp)'
c$$$       allocate( rojp(ngc,      nlxx, nbas),
c$$$     &            sgpb(ngc, nxx, nlxx, nbas),
c$$$     &           fouvb(ngc, nxx, nlxx, nbas),
c$$$     &            sgpp(ngc, ngc, nlxx, nbas),
c$$$     &           fouvp(ngc, ngc, nlxx, nbas) )
c$$$       do ibas = 1,nbas
c$$$          print *,' xxx goto mkjp',ibas
c$$$          call mkjp2(q,ngc, ngvecc, alat, qlat,
c$$$     i      lxx, lx(ibas),nxx, nx(0:lxx,ibas),
c$$$     i      bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas),
c$$$     i      nr(ibas), nrx, rprodx(1,1,0,ibas),
c$$$     o      rojp(1,1,ibas),  sgpb(1,1,1,ibas),
c$$$     o      fouvb(1,1,1,ibas),
c$$$     o      sgpp(1,1,1,ibas),fouvp(1,1,1,ibas) )
c$$$       enddo
c$$$c--- the Coulomb matrix
c$$$       print *,' goto vcoulq'
c$$$       call vcoulq(q, nbloch, ngc,
c$$$     i                  nbas, lx,lxx, nx,nxx,
c$$$     i                  alat, qlat, voltot, ngvecc,
c$$$     i          strx, rojp,rojb, sgbb,sgpb,sgpp, fouvb,fouvp, nblochpmx,
c$$$     o          vcoul)
c$$$       if(allochk)
c$$$     &   write(*,*)'deallocate(strx, rojp,sgpb,sgpp, fouvb,fouvp)'
c$$$       deallocate( strx, rojp,sgpb,sgpp, fouvb,fouvp)
c$$$
c$$$      endif !=============================================================

c----check write
        trwv = 0d0
        do i = 1,nbloch
          trwv = trwv + vcoul(i,i)
        enddo
        write(6,'(" vcoul trwi=",i6,2d22.14)') iqx,trwv
        write(6,'("### sum vcoul(1:ngb,      1:ngb) ",2d22.14,2x,d22.14)')
     & sum(vcoul(1:ngb,1:ngb)), sum(abs(vcoul(1:ngb,1:ngb)))
        write(6,'("### sum vcoul(1:nbloch,1:nbloch) ",2d22.14,2x,d22.14)')
     & sum(vcoul(1:nbloch,1:nbloch)),sum(abs(vcoul(1:nbloch,1:nbloch)))
        write(6,*)
ccccccccccccccccccccccccccccccc
c      vcoul(:, nbloch+1:ngb)=0d0
c      vcoul(nbloch+1:ngb,:)=0d0
ccccccccccccccccccccccccccccccc

 1101   continue
        ngbo=ngb
C... Generate ppmt mattix oct2005 .......................
        if(smbasis()) then
          allocate( ppmt(2,(lxx+1)**2,nbas,ngc) )
          ppmt = 0d0
          call mkppmt(alat,plat,qlat, q, 
     i    ngc, ngvecc, 
     i    rmax, nbas,  bas, lx, lxx,
     o    ppmt) ! ppmt contains value and slove of e(i q+G r) at MT boundaries.
          ! ppmt(2,lmxaa,nbas)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c        write(6,*) 'lxx ppmtsum=',lxx, sum(abs(ppmt))
          write(6,*) 'nbln ngc',nbln,ngc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

C... Matching matrix pmtch. ppmt and prodmt
          pmat(:, nbln+1:nbln+ngc)=0d0
c        write(6,*) 'sss nbln ngc',nbln,ngc
          do igc=1,ngc
c          write(6,*) 'igc=',igc
            pmat(nbloch+igc, nbln+igc) = 1d0
            do ibas= 1, nbas
              do l  =  0, lx(ibas)
                do m  = -l, l
c           write(6,*) 'ibas l m=',ibas,l,m
                  pval= ppmt(1, l**2 + l+1 +m, ibas,igc)
                  pslo= ppmt(2, l**2 + l+1 +m, ibas,igc)
                  do n = 1,nx(l,ibas)
                    if(n==1.and.debug) write(6,"('ttt2: ')")
                    pmat(ibl(m,n,l,ibas), nbln+igc)
     &         =  rdmatch(n,1,l,ibas) * pval 
     &         +  rdmatch(n,2,l,ibas) * pslo
                    if(debug.and.abs(pmat(ibl(m,n,l,ibas), nbln+igc))/=0d0) 
     &        write(6,"('ttt2: i1 i2 pmat=',2i5,2d13.5)")
     &        ibl(m,n,l,ibas), nbln+igc, pmat(ibl(m,n,l,ibas), nbln+igc)
                  enddo
                enddo
              enddo
            enddo
          enddo
          deallocate(ppmt)
          nn = nbln  +ngc ! number for new smooth mixed basis.
          no = nbloch+ngc ! number for original size of mixed basis.
          if(debug) write(6,*) 'end of pmat'

C... oo(no,no). The original overlap matrix.
          allocate( pomat(nn,no) )
          allocate( ppovl(ngc,ngc),oo(no,no))
          call mkppovl2(alat,plat,qlat, 
     i          ngc,  ngvecc, 
     i          ngc,  ngvecc, 
     i          nbas, rmax, bas,  
     o          ppovl)
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


C... oon(nn,nn) is the overlap matrix with new basis
          allocate(oon(nn,nn))
          oon = matmul( dconjg(transpose(pmat(1:no,1:nn)))
     &               ,matmul(oo,pmat(1:no,1:nn)) )


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Reduction of pmat by SVD ... not meaningful
c      nnmx = 1000000
c      epsmx= 1D20
cc      print *,' sumchk pmat=',sum(abs(pmat(1:no,1:nn)))
c      call  zgesvdnn2(nn,nn, nnmx,epsmx,
c     o   oon, ! pmat is reduced to pmat(1:no,1:nnn) by SVD.
c     o   ngcnn)
cc      call  zgesvdnn2(no,ngc, nnmx,epsmx,
cc     i   pmat(1:no,nbln+1:nbln+ngc), ! pmat is reduced to pmat(1:no,1:nnn) by SVD.
cc     o   ngcnn)
cc      nn= nbln+ngcnn
c      write(6,*)' svd ngc ngcnn=',ngc, ngcnn
c      stop 'test end xxxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccc
          if(.false.) then
            open(3011,file='oontest'//charnum5(iqx))
            do ix=nbln+1,nn
              igc=ix-nbln
              qqx(1:3) = (q(1:3)+ matmul(qlat, ngvecc(1:3,igc)))
              absqq  =  sqrt(sum(qqx(1:3)**2))
c           absqg2x(ix) =sum( (2*pi/alat *q0i(1:3,nq0i))**2)
              do iy=nbln+ 1,nn
                igc2=iy-nbln
                if(ix==iy) then
                  write(3011,"('on : ',2i8,3i3,2x,3i3,f13.5,3x,2f20.10)")
     &       ix,iy, ngvecc(1:3,igc),ngvecc(1:3,igc2),absqq, oon(ix,iy)
                else
                  write(3011,"('off:', 2i8,3i3, 2f20.10)")ix,iy,
     &       ngvecc(1:3,igc)-ngvecc(1:3,igc2),  oon(ix,iy)
                endif
              enddo
            enddo
            close(3011)
          endif
ccccccccccccccccccccccccccccccc


C... Generat pomat
!    zmelt_new(K, ij) = \sum_I pomat(K,I)* zmelt(I, ij)
!    means <psi_i psi_j | K> where |K> denote new mixed basis.
!  See sxcf_fal2 and x0kf.
!   Be carefull its transpose procedure---it is a little confusing...
          call pmatorth(oo,oon, pmat(1:no,1:nn), no, nn,
     o    pomat)

ccccccccccccccccccccccccccccccccccccccccccccc
ctttt
c        pomat=0d0
c        do ix= 1,ngb
c          pomat(ix,ix)=1d0
c        enddo
c        do ix= 1,ngb
c        do iy= 1,ngb
c          if(pmat(ix,iy)/=0d0 )
c     &   write(6,"(' ttt: pmat=',2i3,2d13.6)")
c     &   ix,iy,pmat(ix,iy)
c        enddo
c        enddo
c        write(6,"(' ttt:sumchk=',2d13.6,2i4)")
c     &     sum(pomat(:,:)), no,nn
cccccccccccccccccccccccccccccccccccccccccccccc

          if( iqx <= nqibz ) deallocate(oon)
          deallocate(ppovl,oo)
C... Store matching matrix
          write(ifpomat) q,nn,no,iqx
          write(ifpomat) pomat
          deallocate(pomat)
        endif

        if(newaniso()) then
          continue
        elseif(bzcase()==1.and.iqx==1)then
          cycle
        endif

!! == Write out VCCFP ==
        if(debug) write(6,*) 'write out vcoul'
        if(smbasis()) then
          ngb= nn
          allocate(vcoulnn(ngb,ngb))
          vcoulnn= matmul(transpose(dconjg(pmat(1:no,1:nn)))
     &               ,matmul(vcoul(1:no,1:no),pmat(1:no,1:nn)))
          vcoul(1:ngb,1:ngb)= vcoulnn
          deallocate(vcoulnn)
        endif
        if(wvcc) then
          write(ifvcfpout) ngb
          write(ifvcfpout) vcoul(1:ngb,1:ngb),q
        endif  
        write(6,"(' ngc ngb/ngbo=',6i6)") ngc,ngb,ngbo

c Mix0vec ---------------------------------
!! diagonalize the Coulomb matrix
        if(.true.) then
c        if( iqx > nqibz .or. iqx==1) then !feb2012 add iqx==1 for newansio()=T
          if(allochk) write(*,*) 'allocate( ppovl(ngc,ngc))'
          allocate( oo(ngb,ngb) )
          allocate( ppovl(ngc,ngc) )
          call mkppovl2(alat,plat,qlat, 
     &          ngc,  ngvecc, 
     &          ngc,  ngvecc, 
     &          nbas, rmax, bas,  
     o          ppovl)
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
          print *,' --- goto eigen check1 --- '
          allocate(  vcoul0(ngb,ngb) )
          vcoul0 = vcoul(1:ngb,1:ngb)
          if(allochk)
     &    write(*,*) 'allocate(hh(ngb,ngb),oo(ngb,ngb),oox,zz,eb,zzr)'
          allocate(hh(ngb,ngb),zz(ngb,ngb),eb(ngb),zzr(ngb))
          hh  = - vcoul0
c          nmx = 15
          nmx = ngb
          call diagcv(oo,hh,zz,ngb, eb,nmx,1d99,nev)
          do ipl1=1,nev
            if(ipl1==11) write(6,*)' ... '
            if(ipl1>10.and.ipl1<nev-5) cycle
            write(6,'(i4,d23.16)')ipl1,-eb(ipl1)
          enddo
          write(6,"(' nev ngv q=',2i5,3f10.6)")nev,ngb,q

!! === save zz === apr2012takao
c          if( newaniso().and.iqx==1 ) then
c            if(sum(q**2)>1d-10) then
c               stop ' hvccfp0: sanity check. |q(iqx)| /= 0'
c            endif
          write(ifvcoud) ngb
          write(ifvcoud) q
          write(ifvcoud) -eb
          write(ifvcoud) zz
ccccccccccccccccccccccccccccccccccccccc

c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$        print *,' dddddddddddddddd q=',q
c$$$        do ix=1,ngb
c$$$        do iy=1,ngb
c$$$          aaaa=  sum(  dconjg(zz(1:ngb,ix))*matmul( oox,zz(1:ngb,iy))  )
c$$$           if(ix==iy .and.  abs(aaaa-1d0) >1d-8 ) then
c$$$            write(*,*)' dddd zcousum check',ix,iy,aaaa
c$$$           endif
c$$$           if(ix/=iy .and.  abs(aaaa) >1d-8 ) then
c$$$            write(*,*)' dddd zcousum check',ix,iy,aaaa
c$$$           endif
c$$$        enddo
c$$$        enddo
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



          print *
          write(6,'(" eig0 must be equal to the largest =", 2d24.16)')
     &       sum(  dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1))  )
          write(6,'(" zz norm check=",d24.16)')
     &    sum( dconjg(zz(1:ngb,1))*matmul(oox,zz(1:ngb,1)) )
          print *
c          write(6,'(" --- vcoul(exact  no eee)=",d14.6," absq2=",d24.16)')
c     &    fpi*voltot/(sum(tpiba**2*q(1:3)**2))
c     &             , (sum(tpiba**2*q(1:3)**2))
          write(6,'(" --- vcoul(exact)=",d14.6," absq2=",d24.16)')
     &    fpi*voltot/(sum(tpiba**2*q(1:3)**2)-eee)
     &             , (sum(tpiba**2*q(1:3)**2)-eee)
          write(6,'(" --- vcoul(cal ) =",2d14.6)') 
     &    sum( dconjg(zz(1:ngb,1))*matmul( vcoul0,zz(1:ngb,1)) )*voltot
ccccccccccccccccccccccccccccccccccccccccc
c          do igc=1,ngb
c          qqx(1:3) = (q(1:3)+ matmul(qlat, ngvecc(1:3,igc)))
c          write(6,'(" --- vcoul(exact) xxx =",d14.6," absq2=",d24.16)')
c     &    fpi*voltot/(sum(tpiba**2*(qqx(1:3)**2)-eee))
c     &             , (sum(tpiba**2*(qqx(1:3)**2)-eee))
c          write(6,'(" --- vcoul(cal ) xxx =",2d14.6)')
c     &    sum( dconjg(zz(1:ngb,igc))*matmul( vcoul0,zz(1:ngb,igc)) )*voltot
c          enddo
cccccccccccccccccccccccccccccccccccccccccc
          deallocate(vcoul0)

          if(wqt(iqx-nqibz)==0d0) then

C --- To get the vector <Mixed basis| q=0> --------------
cki            if(is_mix0vec()==0) then     !used original befor oct2006
            if(.not.is_mix0vec()) then     !used original befor oct2006
! See switch.F ---> this is not used now.
              ifgb0vec_a =ifgb0vec1
              ifgb0vec_b =ifgb0vec
cki            elseif(is_mix0vec()==1) then !oct2006 new case
            else
! ismix0vec=1 is to avoid problem at BZ boundary when is_mix0vec()=0.
              ifgb0vec_a =ifgb0vec
              ifgb0vec_b =ifgb0vec1
            endif
c1... Case1 to write ifgb0vec -------------------------------------
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
              print *,' igc0=',igc0,ngvecc(1:3,igc0)
              zzr(nbloch+1:nbloch+ngc) = ppovl(1:ngc,igc0)
            endif

            allocate( gbvec(ngb), b0mat(nbloch) )
            print *,' goto mkb0'

C ... get a vector <Product Basis| q+0>
            call mkb0( q, lxx,lx,nxx,nx, aa,bb,nr,nrx,rprodx,
     i        alat,bas,nbas,nbloch, 
     o        b0mat)
            zzr(1:nbloch) = b0mat(1:nbloch)
cccccccccccccccccccccccccccc
c           do igc=1,ngb
c             write(6,"('ssss: ',i5,2d14.6)") igc, zzr(igc)
c           enddo
cccccccccccccccccccccccccccc
            allocate(ooxi(ngb,ngb))
            ooxi=oox
            call matcinv(ngb,ooxi)
            gbvec = matmul(ooxi, zzr)

cccccccccccccc
c           do igc=1,ngb
c             write(6,"('ssss: ',i5,2d14.6)") igc, gbvec(igc)
c           enddo
ccccccccccccccc
            deallocate(ooxi)
            dnorm = sqrt( sum(dconjg(gbvec)*zzr) )
! remove /dnorm at 14June2008. See main/hx0fp0.
! dnorm corresponds to volume (or sum of MT volume if no IPW).
c            gbvec = gbvec /dnorm
c            zzr   = zzr   /dnorm
! Not dnorm=1 at 14June2008. See main/hx0fp0.
c            dnorm=1
            write(ifgb0vec_a,"(3d24.16,2i10,d24.16)") q, ngb,igc0,dnorm
            write(ifgb0vec_a,"(4d24.16)") (gbvec(i),zzr(i),i=1,ngb)
            deallocate( gbvec, b0mat)
c1----------------------------------------------------

c2... --- Case2 to write ifgb0vec c2 is problematic at BZ boundary...------
            dnorm  = 1d0
            zzr(:) = matmul (oox, zz(:,1))
            igc0 = 999999 !dummy now
c phasex ---just to clean. this is irrelevant
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
          if(allochk) !bugfix ---this was in inside or above if 7Feb2006
     &    write(*,*)'deallocate(hh,oo,zz,eb,oox,zzr)'
          deallocate(hh,oo,zz,eb,oox,zzr)
          deallocate(ppovl)
c2---------------------
        endif
        idummy=iclose(trim(vcoudfile))
 1001 continue

      deallocate(ngvecc)
      call cputid(0)
      if(imode==202) stop ' OK! hvccfp0 imode=202 only for Q0P'
c      if(imode==101) stop ' OK! hvccfp0 imode=101 overwirte VCCFP'
c      if(imode==102) stop ' OK! hvccfp0 imode=102 remove_r0c'
      if(imode==0) stop ' OK! hvccfp0 imode=0'
      end

      subroutine checkagree(a,b,char)
      real(8):: a(3),b(3)
      character*(*) :: char
      if(sum(abs(a-b))>1d-6) then
        print *,' Error in checkagree:',char
        stop ' Error in checkagree:'
      endif
      end

      subroutine mkradmatch( p, nxdim,  
     o          rdmatch)
C- make rdmatch
C----------------------------------------------------
Ci  p(1,i): phi     at mt for i-th basis
Ci  p(2,i): dphi/dr at mt for i-th basis
Co rdmatch(nxdim,nxdim)
C-------
Cr    phinew_j(r) =sum_i phi_i(r)* rdmatch (i,j)
Cr     phinew_1(rmt)    =1      phinew_2(rmt)   =0
Cr   d phinew_1(rmt)/dr =0    d phinew_2(rmt)/dr=1
Cr for k >=3
Cr     phinew_k(rmt)    =0
Cr   d phinew_k(rmt)/dr =0
C----------------------------------------------------
      implicit none
      integer(4):: nxdim,lbas,i,i1,i2,ix
      real(8):: p(1:2, 1:nxdim), rdmatch(1:nxdim,1:nxdim)
      real(8):: pd,p1,p1d,p2,p2d,s,t, eps=1d-3,delta
Cr                                       old     new
c      write(6,"('mkradmatch: nxdim=',i4)") nxdim
      if(nxdim <=0) return
      if(nxdim ==1) stop 'mkradmatch err nxdim==1'
      rdmatch=0d0
C... pivot--- get better set of phi for augmentation
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
            stop'mkradmatch: err poor linear dep'
          endif
          i2=i2-1
        endif
        exit
      enddo
C...
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
        if(i==i1.or.i==i2) cycle
        ix=ix+1
c        write(6,"('mkradmatch: i p pd=',i3,2d13.5)") i,p(1,i),p(2,i)
        call phimatch(p(1,i),p(2,i),  p1,p1d,p2,p2d, s,t)
        rdmatch(i,  ix)=  1d0
        rdmatch(i1, ix)=  -s
        rdmatch(i2, ix)=  -t
        write(6,"('mkradmatch: ix st=',i3,2d13.5)") ix,s,t
      enddo
      end

      subroutine phimatch(p,pd, p1,p1d,p2,p2d, s,t)
C --- match for given p and pd
c   phi = s phi1 + t phi2 !slope and value are at MT
c     p  = s p1  + t p2
c     pd = s pd1 + t pd2
      implicit none
      real(8):: matinv(2,2),p,pd,p1,p1d,p2,p2d,s,t,delta,ddd1,ddd2
      delta = p1*p2d-p2*p1d
      matinv(1,1) = 1/delta *  p2d
      matinv(1,2) = 1/delta * (-p2)
      matinv(2,1) = 1/delta * (-p1d)
      matinv(2,2) = 1/delta *  p1
      s = matinv(1,1) *p  + matinv(1,2) *pd
      t = matinv(2,1) *p  + matinv(2,2) *pd
C... check
      ddd1 = abs(s*p1  + t*p2   -  p )
      if(  ddd1 >1d-8 ) stop 'phimatch: ddd1 err'
      ddd2 = abs(s*p1d + t*p2d  -  pd)
      if(  ddd2 >1d-8 ) stop 'phimatch: ddd2 err'
      end

      subroutine pmatorth(oo,oon,pmat,no,nn, pomat)
C get conversion matrix from old mixed basis(no) to augmented mixed basis(nn).
C pmatorth contains
c   oo^{-1}_IJ
      implicit none
      integer(4):: no,nn,io,in,i
      complex(8):: pmat(no,nn),pomat(nn,no),oo(no,no),oon(nn,nn)
      complex(8),allocatable:: ooninv(:,:)
      real(8),allocatable:: eb(:)
      allocate(ooninv(nn,nn))
      ooninv = oon
      call matcinv(nn,ooninv) !generate ooninv
c      pomat = matmul(ooninv, matmul(dconjg(transpose(pmat)),oo))
      pomat = transpose (matmul( oo, matmul(pmat,ooninv)))
      deallocate(ooninv)
      end
c      allocate(pp(nn,nn),ppin(nn,nn),eb(nn),zz(nn,nn),zze(nn,nn))
c      ppin = pp
c      call diagcvh(ppin,nn,eb,zz)
c      do i=1,nn
c        zze(:,i) =  zz(:,i)* sqrt(eb(i))
c      enddo
c      pomat = matmul(pmat, matmul(zze,dconjg(transpose(zz))))

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
      print *,' diagcvv: ngb,nev=',ngb,nev
      do i=1,nev
        write(6,'(i4,d23.16)')i, eb(i)
      enddo
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zgesvdnn2(no,nn, nnmx,epsmx,
     i   pmat,
     o   nnn)
c pmat(no,nn) ---> pmat(no,nnn)
Cio input          pmat(no,nn)
Cio output reduced pmat(no,nnn)
      implicit none
      integer(4):: lwork,info,nn,no,nnn,nnmx,i
      complex(8)::  pmat(no,nn),uu(no,no),vt(nn,nn)
      real(8):: ss(nn),epsmx
      real(8),allocatable:: rwork(:)
      complex(8),allocatable:: work(:),vtt(:,:),pmatx(:,:)
c      print *,' sumchk pmat=',sum(abs(pmat(1:no,1:nn)))
      lwork=4*no
      allocate(work(LWORK),rwork(5*no),pmatx(no,nn))
      pmatx =pmat
      call zgesvd('A','A',no,nn,pmat,no,SS,UU,no,VT,nn,work,lwork,rwork,info)
      nnn=-999
      do i=1,nn
        write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
!         vtt(i,:)=ss(i)*vt(i,:)
        if(nnn==-999.and.ss(i)<epsmx) nnn = i-1
      enddo
c      write(6,*) 'nnn=',nnn
      if(nnn==0) stop 'strange: nnn=0'
      if(nnn>nnmx) nnn=nnmx
      pmat=pmatx
c      pmat(:,1:nnn) = uu(:,1:nnn)
!      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
!     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
!      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk)))
!     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
!      deallocate(vtt)
      end


c---------------------------------------------------------------------
      subroutine mkb0( q, lxx,lx,nxx,nx, aa,bb, nrr,nrx,rprodx,
     i        alat,bas,nbas,nbloch, 
     o        b0mat)
C--make the matrix elementes < B_q | exp(iq r)>
      implicit none
      integer(4) :: nlx,l,n,m,nr,ir,lm,ibl1,ibas,nrx,nbloch

      integer(4) :: nbas,lxx, lx(nbas), nxx, nx(0:lxx,nbas),nrr(nbas)
      real(8)    :: rprodx(nrx,nxx,0:lxx,nbas),aa(nbas),bb(nbas),
     &   phi(0:lxx),psi(0:lxx), bas(3,nbas),
     &   alat,
     &   pi,fpi,tpiba,qg1(3),q(3),absqg,r2s,a,b
c
      complex(8) :: b0mat(nbloch),img=(0d0,1d0) ,phase
c
      integer(4),allocatable:: ibasbl(:), nbl(:), lbl(:), lmbl(:)
      real(8),allocatable :: ajr(:,:),rofi(:),rob0(:,:,:)
      real(8),allocatable::cy(:),yl(:)
      complex(8),allocatable :: pjyl(:,:)
#ifdef COMMONLL
      integer(4) ll(51**2)
      common/llblock/ll
#else
      integer(4) ll
#endif

c-----
      write(6,*)'mkb0:'
      pi   = 4d0*datan(1d0)
      fpi  = 4*pi
      nlx  = (lxx+1)**2
c
      tpiba = 2*pi/alat
      qg1(1:3) = tpiba * q(1:3)
      absqg    = sqrt(sum(qg1(1:3)**2))
c
      allocate(ajr(1:nrx,0:lxx), pjyl(nlx,nbas),rofi(nrx),
     &  ibasbl(nbloch), nbl(nbloch), lbl(nbloch), lmbl(nbloch),
     &  cy(nlx),yl(nlx),rob0(nxx,0:lxx,nbas))
c
      call sylmnc(cy,lxx)
      call sylm( qg1/absqg,yl,lxx,r2s) !spherical factor Y( q+G )
c
      do ibas = 1,nbas
        a = aa(ibas)
        b = bb(ibas)
        nr= nrr(ibas)
        rofi(1)    = 0d0
        do ir      = 1, nr
          rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
          call bessl(absqg**2*rofi(ir)**2,lx(ibas),phi,psi)
          do l  = 0,lx(ibas)
c ... bessel function
            ajr(ir,l) = phi(l)* rofi(ir) **(l +1 )
            ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
          enddo
        enddo

c ... Coefficients for j_l yl  on MT  in the expantion of of exp(i q r).
        phase = exp( img*sum(qg1(1:3)*bas(1:3,ibas))*alat  )
        do lm = 1,(lx(ibas)+1)**2
          l = ll(lm)
          pjyl(lm,ibas) = fpi *img**l *cy(lm)*yl(lm) *phase  *absqg**l
        enddo
c ... rob0
        do l = 0,lx(ibas)
          do n = 1,nx(l,ibas)
            call gintxx( ajr(1,l), rprodx(1,n,l,ibas), a,b,nr,
     o                 rob0(n,l,ibas) )
          enddo
        enddo
      enddo

c ... index (mx,nx,lx,ibas) order.
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
c ... pjyl * rob0
      do ibl1= 1, nbloch
        ibas= ibasbl(ibl1)
        n   = nbl  (ibl1)
        l   = lbl  (ibl1)
        lm  = lmbl (ibl1)
        b0mat(ibl1) = pjyl(lm,ibas) * rob0(n,l,ibas)
      enddo
      deallocate(ajr, pjyl,rofi,
     &  ibasbl, nbl, lbl, lmbl,
     &  cy,yl,rob0)
      end
