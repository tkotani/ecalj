      program hbasfp0_v2
c-- Generates orthonormal optimal product basis and required radial integrals in each MT.
c input files
c  GWIN_V2   : input data for GW
c  LMTO   : fundamental data for crystal
c  PHICV  : radial functions Valence and Core
c
c output files
c  BASFP//ibas :: product basis for each ibas
c  PPBRD_V2_//ibas :: radial <ppb> integrals. Note indexing of ppbrd
c
c The main part of this routine is in the subroutine basnfp_v2
c-takao Apr 2002 -----
      use m_rgwinf_v3,only:rgwinf_v3, alat,nclass,natom,nspin,nl,nnv,nnc,nrx, cutbase,lcutmx,nindxc,
     &     nindxv,occv,unoccv,occc,unoccc,iclass
      use keyvalue,only: getkeyvalue
      implicit none
      real(8):: qbas(3,3),ginv(3,3)
      integer(4)::
     l     ifphiv(2),ifphic(2),
     l     iphiv(2),iphivd(2),iphic(2),iphi(2),iphidot(2),
     l     ifev(2),ifevf(2),ibas,ibas1,ic,icx,ifaln,ifinin,iflmto,ifphi,ii,ir,irad,isp,ix,lmx,lmx2,n,nbas,ncoremx,
     &     l,nn,icore,ifianf,nphi,nradmx,nsp,iopen,maxnn,iclose
      integer(4),allocatable:: mdim(:),lcutmxa(:)
      character(12) :: aaa
      integer(4),allocatable::nrofi(:), nocc(:,:),nunocc(:,:),nindx(:,:)
      logical :: ptest=.false. !See ptest in hvccfp0.f
      real(8),allocatable :: bb(:),zz(:), phic(:,:)
      integer(4) :: ndat
      integer(4),allocatable:: ncindx(:,:),
     &           lcindx(:,:),
     &           nrad(:),
     &           nindx_r(:,:),
     &           lindx_r(:,:),
     &           nc_max(:,:),ncore(:)
      real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:)
     &                     ,phitotr(:,:,:,:,:)
      character*11 :: ffaln
      integer(4)::incwfin,ret
      integer(4),allocatable:: idid(:),ibasf(:)
      logical :: checkdid ,anfexist
      integer(4):: iread, idummy
c--------------------------------------------------------------------------
      ifinin=-99999 !dummy
      write(6,'(a)') ' --- Input normal(=0); coremode(=3);'//
     & ' ptest(=4); Excore(=5); for core-valence Ex(=6);'//
     & ' val-val Ex(7);  normal+<rho_spin|B> (8); version(-9999) ?'
      call readin5(ix,iread,idummy)
      call headver('hbasfp0',ix)
      if(ix==3) then
        write(6,*)'  ### coremode; Product basis for SEXcore ### '
        incwfin = -2
      elseif(ix==0) then
        write(6,*)'  ### usual mode use occ and unocc for core ### '
        incwfin = 0
      elseif(ix==4) then
        write(6,*)
     & ' ### ptest mode. now special for Q0P. GWIN_V2 is neglected ### '
        write(6,*) '  See basnfp.f of ptest section.'
        incwfin = 0
      elseif(ix==5) then
        write(6,*)
     &  '  ### calculate core exchange energy ### ix==5'
        incwfin = 0
      elseif(ix==6) then
        write(6,*)
     &  '  ### calculate p-basis for core-valence Ex ix==6'
        write(6,*) ' occ=1:unocc=0 for all core'
        incwfin = -3
      elseif(ix==7) then
        write(6,*)
     &  '  ### calculate p-basis for val-val Ex ix==7'
        write(6,*) ' occ=0:unocc=0 for all core'
        incwfin = -4
      elseif(ix==8) then !May2005
        write(6,"('  ### usual mode use occ and unocc for core',
     &            ' and <rho_spin |B(I)> ### ')")
        incwfin = 0
      else
        print *,' hbasfp: Input must be 0 or 3 or 4'
        stop    ' hbasfp: Input must be 0 or 3 or 4'
      endif

      iflmto     = iopen('LMTO',1,0,0)
      if (iflmto <= 0) stop 'unit file for LMTO <= 0'
!!
      call rgwinf_v3 (iflmto,ifinin,incwfin) !these are inputs.
      nsp=nspin
      ! read data in m_rgwinf_v3
      ! Output are allocated and data are setted as above.
      print *,'end of rgwinf'

c-- readin lcutmxa ------------
      call getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifinin,status=ret)
      allocate(lcutmxa(1:natom))
      do
        read(ifinin,*,err=980) aaa
        if(aaa=='lcutmx(atom)') then
          read(ifinin,*) lcutmxa(1:natom)
c          write(6, '(" lcutmxa=",20i3)' ) lcutmxa(1:natom)
          goto 990
        endif
      enddo
  980 continue
      lcutmxa=lcutmx
  990 continue
      close(ifinin)
      write(6,"(' lcutmxa=',$)")
      write(6,'(20i3)') lcutmxa(1:natom)

      lmx        = 2*(nl-1)
      lmx2       = (lmx+1)**2
      nn         = maxnn (nindxv,nindxc,nl,nclass)
      nphi       = nrx*nl*nn*nclass

c -optimal orthonormal product basis
c> reindex nocc,nunocc,nindx
! For valence  from GWIN_V2
! occv   : occ    switch
! unoccv : unocc  switch
! nindexv: n index
!-----------------------------
! For core  from GWIN_V2
! occc   : occ  switch
! unoccc : unocc switch
! nindexc: n index
!-----------------------------
! For valence+core
! nocc
! nunocc
! nindx
      allocate( nocc(nl*nn,nclass),nunocc(nl*nn,nclass),
     &          nindx(nl,nclass) )
      call reindx  (occv,unoccv,nindxv,
     i              occc,unoccc,nindxc,
     d              nl,nn,nnv,nnc,nclass,
     o              nocc,nunocc,nindx)
      print *,' --- end of reindx ---'

c-----------
c read PHIVC  and reserve it to phitot
c----------
      ifphi  = iopen('PHIVC', 0,-1,0)     ! PHIV+PHIC augmentation wave and core
      read(ifphi) nbas, nradmx, ncoremx
      allocate(  ncindx(ncoremx,nbas),
     &           lcindx(ncoremx,nbas),
     &           nrad(nbas),
     &           nindx_r(1:nradmx,1:nbas),
     &           lindx_r(1:nradmx,1:nbas),
     &        aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) ,
     &        phitoto(nrx,0:nl-1,nn,nbas,nsp),
     &        phitotr(nrx,0:nl-1,nn,nbas,nsp),
     &        nc_max(0:nl-1,nbas),ncore(nbas) )

      read(ifphi) nrad(1:nbas)
      read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
c
      nc_max=0
      do ibas=1,nbas
        write(6,*)' --- read PHIVC of ibas=',ibas
        ic = ibas
        read(ifphi) ncore(ic), ncoremx                            !core
        read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
        read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
        if(ic/=icx) then
          print *, 'ic icx=',ic,icx
          stop 'hbasfp0: ic/=icx'
        endif
        read(ifphi) rr(1:nrofi(ic),ic)

        do isp = 1, nsp
          write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
          do icore = 1, ncore(ic)
            l =  lcindx(icore,ic)
            n =  ncindx(icore,ic)
c           read(ifphi) gcore(1:nrofi(ic),icore, ic,isp)   !core
            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
            phitotr(1:nrofi(ic),l,n, ic,isp)=              !core raw= core orthgonal
     &      phitoto(1:nrofi(ic),l,n, ic,isp)               !
c            print *,' sss1c=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
            if(n>nc_max(l,ic)) nc_max(l,ic)=n
          enddo

          do irad = 1, nrad(ic)
            l = lindx_r (irad,ic)
            n = nindx_r (irad,ic) + nc_max(l,ic)
c           read(ifphi) gx(1:nrofi(ic),l, n, ic,isp)
            read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
            read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
c            print *,' sss1=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
c            print *,' sss2=',sum(abs(phitotr(1:nrofi(ic),l,n, ic,isp)))
          enddo
        enddo
      enddo
c-----------


c-----check write
      ffaln ='PHIV.chk'
      ifaln = iopen(ffaln,1,-1,0)
      do ibas = 1,nbas
        ic = ibas
        do irad = 1, nrad(ic)
          l = lindx_r (irad,ic)
          n = nindx_r (irad,ic) + nc_max(l,ic)
          write(ifaln,"(a,5i5)")'------- ibas l n =',ibas,l,n
          do ir=1,nrofi(ic)
            write(ifaln,"(3d24.15)")rr(ir,ic), phitotr(ir,l,n,ic,1:nsp)
          enddo
c        ifaln = iclose(ffaln)

cccccccccccccccccccccccccc
c check write
c        if(ibas==1.and.l==2) then
c        write(2001,"(a,5i5)")'- - - --- ibas l n =',ibas,l,n
c        do ir=1,nrofi(ic)
c          write(2001,"(3d24.15)")rr(ir,ic), phitoto(ir,l,n,ic,1:nsp),
c     &                                      phitotr(ir,l,n,ic,1:nsp)
c        enddo
c        endif
ccccccccccccccccccccccccc

        enddo
      enddo
      ifaln = iclose(ffaln)

C ...  excore mode ---------
      if(ix==5 ) then
        call excore(nrx,nl,nnc,nclass,nsp,natom,
     &   phitotr(1:nrx,0:nl-1,1:nnc,1:nclass,1:nsp),   !core     
     &   nindxc,iclass,
     &    aa,bb,nrofi,rr)
        goto 998
      endif


C ... READIN ANFcond ------------------------------------
      inquire(file='ANFcond',EXIST=anfexist)
      print *,' ANFcond exist? =',anfexist
      if(anfexist) then
        allocate(ibasf(natom),idid(natom)) !,transaf(1:3,natom))
        ifianf = 211
        open(ifianf,file='ANFcond')
        read(ifianf,*)
        read(ifianf,*)
        read(ifianf,*)  !anfvec(1:3)
!       write(6,'(" ANFcond=",3f13.6)') anfvec
        do ibas=1,natom
          read(ifianf,*) ibas1, ibasf(ibas)
          if(ibas1/=ibas) stop ' readin ANFcond error:ibas1/=ibas'
          write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
        enddo
        close(ifianf)

C ... Check iclass =ibas ; CLASS file contains true classs information.
        do ibas=1,natom
          if(iclass(ibas)/=ibas) stop ' iclass(ibas)/=ibas: '
        enddo
        ii=0
        do ic=1,nclass
          ibas=ic
          if( checkdid(idid,ii, ibas) ) cycle
          phitotr(:,:,:,ibasf(ibas), :)=phitotr(:,:,:,ibas, :)
          ii=ii+2
          idid(ii-1)= ibas
          idid(ii  )= ibasf(ibas)
        enddo
        if( sum (idid(1:ii)) /= natom*(natom+1)/2) 
     &    stop 'hbasfp0:sum (idid(1:ii)) /= n(n+1)/2'
        print *,' end of readin ANFcond:'
      endif

ccccc
c      allocate( phitot(nrx,nl,nn,nclass, 2,nsp))
c      call dcopy(nphi,w(iphi(1)),1,    phitot(1,1,1,1,1,1),1)
c      call dcopy(nphi,w(iphidot(1)),1, phitot(1,1,1,1,2,1),1)
c      if(nsp==2) then
c        call dcopy(nphi,w(iphi(2)),1,    phitot(1,1,1,1,1,2),1)
c        call dcopy(nphi,w(iphidot(2)),1, phitot(1,1,1,1,2,2),1)
c      endif
c
c------------------------------
c      call genbasnfp  (
c     i              w(inindx),w(inocc),w(inunocc),nl,nn,nclass,
c     i              w(ia),w(inrofi),w(ir),nrx,
c     i              phitot,nsp
c     i   ,cutbase,lcutmxa)
c------------------------------
      do ic = 1,nclass
c        nr  = nrofi(ic)
c        ac  = aa(ic)
c        rws = rr(nr,ic)
c        xmt = dble(nr-1)*ac
c        b   = rws/( dexp(xmt) - 1.d0 )
        call basnfp_v2(nocc(1,ic),nunocc(1,ic),nindx(1,ic), ! Product Basis functions
     &    nl,nn,nrx, nrofi(ic),rr(1,ic),aa(ic),bb(ic),ic,
     &    phitoto,phitotr,nsp,nclass,
     i    cutbase, lcutmxa(ic),ix,iread,alat
     i  ,nc_max(0,ic) )
      end do
      if(ix==0) stop ' OK! hbasfp0 ix=0 normal mode '
      if(ix==3) stop ' OK! hbasfp0 ix=3 core mode '
      if(ix==4) stop ' OK! hbasfp0 ix=4 ptest mode  '
      if(ix==6) stop ' OK! hbasfp0 ix=6 Exx core-val mode  '
      if(ix==7) stop ' OK! hbasfp0 ix=7 Exx val-val mode  '
      if(ix==8) stop ' OK! hbasfp0 ix=8 normal(ix==0) + <B|spin den> '
 998  if(ix==5) stop ' OK! hbasfp0 ix=5 ex core mode  '
      end


      logical function checkdid (idid,ii, ibas)
      integer(4):: idid(ii),ix
      checkdid=.true.
      do ix=1,ii
        if(idid(ix)==ibas) return
      enddo
      checkdid=.false.
      end








