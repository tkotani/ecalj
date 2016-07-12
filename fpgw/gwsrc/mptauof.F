c ------------------------------------------------------------------------------
      subroutine getw0(llw,ii,ie,nq0i,dmlx,epinvq0i,wklm,wbz,lmxax,q0i,epinv, w0,llmat)
!! == Obtain effective screened Coulomnb interaction W-V at q=0 ==
!! Output
!!     w0(ii:ie), others are inputs.
!! INPUT
!! ii: start index of iw
!! iw: end index of iw
      implicit none
      integer,intent(in)::  ii,ie,nq0i,lmxax
      real(8),intent(in)::  epinvq0i(nq0i,nq0i),dmlx(nq0i,9), wklm((lmxax+1)**2),wbz,epinv(3,3,nq0i),q0i(3,nq0i)
      complex(8),intent(in):: llw(ii:ie,nq0i)
      complex(8),intent(out):: w0(ii:ie),llmat(3,3)
!!
      integer:: lm,lm1,lm2,nlklm,iw,ll,lm1x,lm2x
      complex(8):: llwyl(9,ii:ie),llw_invr(ii:ie,1:nq0i)
      real(8):: epinvq0i_m1(nq0i,nq0i)
      real(8)::  cg((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2)
      complex(8)::   cgllw((lmxax+1)**2,(lmxax+1)**2),sss,sss2
      real(8):: fpi,pi,r2s,emat(3,3),rrr(3),qnorm2(nq0i),ylm,rrr2(3)
      integer:: nlmm,lx,iq0i,iq0x, i
      real(8),allocatable:: cy(:),yl(:),yl2(:)
!!
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      nlklm=(lmxax+1)**2
      epinvq0i_m1=epinvq0i
      write(6,"(' getw0 start:',4i3,d13.6)") nq0i,lmxax,ii,ie,sum(abs(llw(ii:ie,1:nq0i)))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      nlm=(lmxax+1)**2
c      do lm=1,nlm
c        if(abs(wklm(lm))>1d-10) write(6,*)lm,ll(lm),wklm(lm)
c      enddo
c      stop 'test end vvvvvvvvvvvvvvvvvvvvvv'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call matinv(nq0i,epinvq0i_m1)
      llwyl=0d0
      do iw= ii,ie
        llw_invr(iw,1:nq0i) = matmul (epinvq0i_m1, llw(iw,1:nq0i))
        ! llw_invr: representation in a linear combination of invariat tensor of epsilon.
        lm=1
        llwyl(lm,iw) = sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
        do lm=5,9
          llwyl(lm,iw)= sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
        enddo
c   write(*,"('qqq iw llw',i3,3(2d12.4,x),3x,20(2d12.4,2x))")iw,llw(iw,1:nq0i),llwyl(1,iw),llwyl(5:9,iw)
      enddo

!! omega=0 3x3 matrix, 
      llmat=0d0
      do iq0i=1,nq0i
        llmat(1:3,1:3) = llmat(1:3,1:3) + llw_invr(ii,iq0i)*epinv(1:3,1:3,iq0i)
      enddo

      call rotcg(lmxax,(/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),1,cg) !no rotation

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      lx=2
c$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
c$$$      call sylmnc(cy,lx)
c$$$      do iw= ii,ie
c$$$        do iq0x=1,nq0i
c$$$        sss=0d0
c$$$        do iq0i=1,nq0i
c$$$            sss= sss  + llw_invr(iw,iq0i)*sum(q0i(:,iq0x)*matmul(epinv(:,:,iq0i),q0i(:,iq0x)))
c$$$     &                 /sum(q0i(:,iq0x)**2)
c$$$        enddo
c$$$          write(*,"(' ttt: epinv expansion=',2i4,2f10.5,2x,2d13.5)") iq0x,iw,sss,sss-llw(iw,iq0x)
c$$$        enddo
c$$$        write(*,*)
c$$$      enddo
c$$$      stop '---- ttt test1: reproduce llw at q0i -----------------------'

c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$!! === test for one r vector as for <ehat| epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat) ===
c$$$!! ===== generate YL for a test vector rrr (rrr is ehat above).=====
c$$$      lx=2
c$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
c$$$      call sylmnc(cy,lx)
c$$$      do i=1,3
c$$$         if(i==1) rrr =(/1d0,1d0,0d0/)
c$$$         if(i==2) rrr =(/1d0,-1d0,0d0/)
c$$$         if(i==3) rrr =(/0d0,0d0,1d0/)
c$$$         rrr = rrr/sqrt(sum(rrr**2))
c$$$c     write(*,"(' testttt: r=',3f10.5)") rrr
c$$$         call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
c$$$         do iw= ii,ie
c$$$            sss=0d0
c$$$            do iq0i=1,nq0i
c$$$               sss = sss + llw_invr(iw,iq0i)*sum(rrr*matmul(epinv(:,:,iq0i),rrr))
c$$$            enddo
c$$$c            if(abs(llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss)>1d-12) then
c$$$            write(*,"(' ttt: epinv expansion=',i3,2f10.5,2x,4d13.5)") iw,sss,
c$$$     &      llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss, llw(iw,i)-sss
c$$$c            endif
c$$$         enddo
c$$$      enddo
c$$$      stop '---- ttt testxxx: reproduce llw at q0i -----------------------'

c$$$      lx=lmxax
c$$$      if(allocated(cy)) deallocate(cy)
c$$$      if(allocated(yl)) deallocate(yl)
c$$$      allocate(cy((lx+1)**2),yl((lx+1)**2))
c$$$      call sylmnc(cy,lx)
c$$$      rrr=(/-0.36d0,0.20d0,0.4d0/)
c$$$      rrr=rrr/sqrt(sum(rrr**2))
c$$$      call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
c$$$      write(*,"(  ' test input: q=',3f10.5)") rrr

!! ---
      w0=0d0
      do iw =ii,ie
        lm2x=0
        cgllw=0d0
        do lm2=1,nlklm
          if(mod(ll(lm2),2)==1) cycle  !only even l
          lm2x=lm2x+1
          lm1x=0
          do lm1=1,nlklm
            if(mod(ll(lm1),2)==1) cycle !only even l
            lm1x=lm1x+1
            cgllw(lm1x,lm2x) = cg(lm2,1,lm1)*llwyl(1,iw) + sum(cg(lm2,5:9,lm1)*llwyl(5:9,iw))
          enddo
        enddo

!! === inversion of 1 = \sum_L1 llwyl(L1)Y_L1 * \sum_L2 K_L2 Y_L2===
!! Both sides are multipled by Y_L and spherically integraled.
!! warn: Klm should be not confused with wklm
!! Klm is defeined in Christoph's paper PRB 125102-8
        call matcinv(lm2x,cgllw(1:lm2x,1:lm2x))
        cgllw(1:lm2x,1:lm2x)= sqrt(fpi)*cgllw(1:lm2x,1:lm2x)
        !sqrt(fpi) comes from \int d\Omega Y_0 (right hand side of inversion eq).

!! === spherical integral of Klm ===
        lm1x=0
        do lm1=1, nlklm
          if(mod(ll(lm1),2)==1) cycle
          lm1x=lm1x+1
          w0(iw)= w0(iw)+ wklm(lm1)* fpi* cgllw(lm1x,1)/wbz      ! Klm=cgllw(:,1)
          if(lm1==1) w0(iw)= w0(iw) - wklm(lm1)*fpi**1.5d0/wbz   ! fpi**1.5 --> subtract Coulomb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          w0(iw)= wklm(1)*fpi**1.5*(1d0/llw(iw,1)-1d0)/wbz !test case of llw at q0i(:,iq0i=1).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$        lx=lmxax
c$$$        write(6,*)' lx=',lx
c$$$        if(allocated(cy)) deallocate(cy)
c$$$        if(allocated(yl)) deallocate(yl)
c$$$        allocate(cy((lx+1)**2),yl((lx+1)**2))
c$$$        call sylmnc(cy,lx)
c$$$        rrr=(/1d0,1d0,0d0/)
c$$$c        rrr=(/1d0,-1d0,0d0/)
c$$$c        rrr=(/0d0,0d0,1d0/)
c$$$        rrr=rrr/sqrt(sum(rrr**2))
c$$$        call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
c$$$c        write(*,"(  ' qqq test input: q=',3f10.5)") rrr
c$$$        sss=0d0
c$$$        lm2x=0
c$$$        do lm2=1,nlklm
c$$$           if(mod(ll(lm2),2)==1) cycle !only even l
c$$$           lm2x=lm2x+1
c$$$           sss= sss + cgllw(lm2x,1) *cy(lm2)*yl(lm2)
c$$$        enddo
c$$$        write(*,"(' qqq ep test:',i3,2f10.5,2x,2f10.5,2x,2f10.5)") iw,sss,
c$$$     &       1d0/llw(iw,1)
c$$$c     &      1d0/(llwyl(1,iw)*cy(1)*yl(1) + sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9)))
c$$$      stop ' --- test3 ttt vvvvvvvvvvvvvvvvv test end vvvvvvvvvvvvvvvvv ---'
      enddo
c$$$      do iw=ii,ie
c$$$        write(6,"('w0 & w0(1:nq0i)=',2d11.3,10(1x,2d11.3))")
c$$$     &  w0(iw), (wklm(1)*fpi**1.5*(1d0/llw(iw,iq0i)-1d0)/wbz,iq0i=1,nq0i)
c$$$      enddo
      end

      integer function ll(ilm)
      integer:: lla(169)
      data lla/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11,25*12 /
      ll=lla(ilm)
      end

      subroutine eibzgen(nqibz,symgg,ngrp,qibze,iqxini,iqxend,qbz,nqbz,timereversal,ginv,iprintx,
     o nwgt,igx,igxt,eibzsym,timerout)
      implicit none
!! === Obtain info for eibz symmetrization. See PRB81,125102-9 ===
!! For GaAs 4x4x4 (with timereversal), we have IBZxBZ=10x640 is reduced to be IBZxEBZ=286.
      integer:: nqibz,ngrp,iqxini,iqxend,iq,ig,igxx,neibzx,ieibz,ibz,nwgtsum,itimer,ntimer
      integer:: neibz(iqxini:iqxend),nwgt(nqbz,iqxini:iqxend),ik,nqbz,i,it
      logical:: timereversal
      real(8):: symgg(3,3,ngrp),qibze(3,iqxini:iqxend),qbz(3,nqbz),qeibz(3,nqbz,iqxini:iqxend),
     &  q(3),qx(3),qlat(3,3),qdiff(3),ginv(3,3),qxx(3),ddd,timer,qxi(3) ,
     &  sss(3,3),sumcheck1 !qout(3,nqbz,iqxini:iqxend,nqbz)
      integer:: igx(ngrp*2,nqbz,iqxini:iqxend),igxt(ngrp*2,nqbz,iqxini:iqxend),imm,im
      integer:: eibzsym(ngrp,-1:1,iqxini:iqxend),ibzxi(nqbz),ibzx,eallow(ngrp,-1:1),immxx
      logical:: iprintx,timerout

      integer,allocatable:: immx(:),mqbz(:)
      integer,allocatable:: nxx(:,:,:),ibznxx(:,:),nww(:,:),igx_(:,:),igxt_(:,:)

      real(8)::    epd=1d-5
      integer:: nkey(3),isig,kkk,kkk3(3),ik1,ik2,ik3,ix,ifi,pwmode,ifiese
      integer,allocatable:: ieord(:)
      integer,allocatable:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
      real(8),allocatable:: qbzl(:,:)
      logical :: inotime,nexist,debug=.false.
!!

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cnnn
      real(8):: alatx, plat(3,3)
      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw
      real(8)::ppin(3)
      real(8):: rlatp(3,3),xmx2(3),tolq=1d-8

      ifi = 6661
      open (ifi, file='LATTC')
      read(ifi,*) alatx
      read(ifi,*) plat(1:3,1)
      read(ifi,*) plat(1:3,2)
      read(ifi,*) plat(1:3,3)
      close(ifi)
      call minv33tp(plat,qlat)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c ginvtest
ccccccccccccccccccccc
c      print *,'ginv test'
c      ginv=ginv/2d0
c      write(6,"('ginv=',3f9.4)")ginv(:,:)
ccccccccccccccccccccc

      if(iprintx) write(6,"('eibzgen: TimeReversalSwitch ngrp= ',l,i3)") timereversal,ngrp
      ntimer=1
      if(timereversal) ntimer=-1
      sss(1,:)=(/1d0,0d0,0d0/)
      sss(2,:)=(/0d0,1d0,0d0/)
      sss(3,:)=(/0d0,0d0,1d0/)
      sumcheck1= sum(abs(symgg(:,:,1)-sss))
      if(sumcheck1>tolq) call rx( 'eibzgen: symgg(:,:,1) is not E')

!!=== Get key and nkey for each ix. (similar with those in readqg,readeigen) ===
      allocate(qbzl(3,nqbz),key(3,0:nqbz),ieord(nqbz))
      key(:,0)=0 !dummy
      do ibz=1,nqbz
        call rangedq(matmul(ginv,qbz(:,ibz)), qbzl(1,ibz))
cccccccccccccccccc
c            write(6,"('ibz=',i5,3f8.4)")ibz,qbzl(:,ibz)
cccccccccccccccccc
      enddo
      do ix =1,3
         call sortea(qbzl(ix,:),ieord,nqbz,isig)
         ik=0
         do i=1,nqbz
            kkk= ( qbzl(ix,ieord(i))+0.5*epd )/epd  !kkk is digitized by 1/epd
            if(i==1.or.key(ix,ik)<kkk) then
               ik=ik+1
               key(ix,ik) = kkk
c               write(6,*)ix, ik,i, key(ix,ik), qbzl(ix,ieord(i))
            elseif (key(ix,ik)>kkk) then
               write(6,*)ix, ik,i, key(ix,ik), qbzl(ix,ieord(i))
               call rx( 'iqindx: bug not sorted well')
            endif
         enddo
         nkey(ix)=ik
      enddo
      deallocate(ieord)
!!  key is reallocated. inverse mattping, iqkkk
      allocate( kk1(nkey(1)),kk2(nkey(2)),kk3(nkey(3)) )
      kk1(:) = key(1,1:nkey(1))
      kk2(:) = key(2,1:nkey(2))
      kk3(:) = key(3,1:nkey(3))
      deallocate(key)
      allocate( iqkkk(nkey(1),nkey(2),nkey(3)) )
      if(debug) write(6,*)' initqqq nqbz=',nqbz
      do i=1,nqbz
         kkk3= (qbzl(:,i)+0.5*epd)/epd !kkk is digitized by 1/epd
         call tabkk(kkk3(1), kk1,nkey(1), ik1)
         call tabkk(kkk3(2), kk2,nkey(2), ik2)
         call tabkk(kkk3(3), kk3,nkey(3), ik3)
         iqkkk(ik1,ik2,ik3)=i
         if(debug) write(6,"(' iqkkk ik1,ik2,ik3,=',i4,3i3,x,3i6)")i,ik1,ik2,ik3,kkk3(:)
      enddo
      if(debug) then
        write(6,*)'kk1=',kk1
        write(6,*)'kk2=',kk2
        write(6,*)'kk3=',kk3
      endif
      deallocate(qbzl)
!! -----------------------------------------------------------
ccccccccccccccccccc
      ifiese=9039
      inquire(unit=ifiese,opened=nexist)
      if(nexist) stop 'eibzgen: I assumed ifio=9039 is not used(too sloppy) pwmode'
      open(ifiese,file='PWMODE')
      read(ifiese,*)pwmode
      close(ifiese)
cccccccccccccccccc
      if(pwmode>0.and.pwmode<10) call shortn3_initialize(qlat,rlatp,xmx2)

!! ===  main iq loop ===
      eibzsym = 0
      do iq = iqxini,iqxend
         if(debug) write(6,"('iq =',i7,'  xxxxxxxxxxxxxxxxxxxxxxx')") iq
        q = qibze(:,iq) !q means k in eq.50 of PRB81,125102-9
!! Allowed operation to keep q,   eallow(ig,it) \in EIBZ(q)
        eallow=0d0
        do it=1,ntimer,-2
          timer=dble(it)
        do ig=1,ngrp
          qdiff= timer*matmul(symgg(:,:,ig),q) - q
          call rangedq(matmul(ginv,qdiff), qxx)
          if(debug) write(6,*)'ig qdiff',ig,qdiff,q,ginv,qxx
cccccccccccccccccccccccccccccccccccccccccccccc
          if(pwmode>0.and.pwmode<10) then
c     check wheter q is on 1st BZ boundary or not
             ppin=matmul(transpose(plat),q) 
             call shortn3(rlatp,xmx2,ppin,noutmx, nout,nlatout)
             if(ig/=1.and.nout>1) then
                do iout=1,nout
                   write(*,"(a,3i5,f10.4,3f8.4)")'rrrrn1 =',nlatout(:,iout),
     &                  sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
     &                  matmul(qlat(:,:),ppin+nlatout(:,iout))
                enddo
                print *
                cycle
             endif  
          endif
ccccccccccccccccccccccccccccccccccccccccccccc
          if(sum(abs(qxx))<tolq) then
            eallow(ig,it)=1
          endif  
        enddo
        enddo
        eibzsym(:,:,iq)=eallow ! === eibzsym is eallow -----

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        write(6,*)'xxx nqbz=',nqbz,ngrp,ntimer, 8*nqbz*ngrp/1000000000.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!! === Generate all nxx (integered q with epditit) generated from ibz. ===
        allocate(nxx(3,nqbz,ngrp*2),immx(nqbz)) !,nww(nbz,ngrp*2))
        allocate(igx_(ngrp*2,nqbz),igxt_(ngrp*2,nqbz))
        immxx=0
        do ibz=1,nqbz
          if(debug) write(6,"(' qbz=',i3,' ',3f9.4,' ',3f9.4)")ibz,qbz(:,ibz),matmul(ginv,qbz(:,ibz))
          imm=0
          do it=1,ntimer,-2 !ntimer is 1 or -1.
            timer = dble(it)
          do ig=1,ngrp
            if(eallow(ig,it)/=1) cycle !eallow=1 keeps q invariant. Allowed operation of EIBZ(q)
            qx = timer*matmul(symgg(:,:,ig),qbz(:,ibz))
            qxi= matmul(ginv,qx)
            call rangedq(qxi, qxx) ! qxx(i) takes from 0 to 1.
            do im=1,imm
              if(sum(abs(nxx(:,ibz,im) - (qxx(:)+.5*epd)/epd)) <2) then
c                 nww(ibz,im) = nww(ibz,im)+1
                 goto 3012      !<2 for safe
              endif   
            enddo
            imm=imm+1
            nxx(:,ibz,imm)= (qxx+.5*epd)/epd  !integer nxx(i) takes from 0 to 1/epd.
            igx_ (imm,ibz) = ig  !inefficient memory usage, I think.
            igxt_(imm,ibz) = it  !
            if(debug) write(6,"('imm ibz it ig nxx qxx qxi=',4i7,x,3i6,x, 3f8.4,x,3f8.4,' qx=',3f8.4)")
     &       imm,ibz,ig,it ,nxx(:,ibz,imm), qxx, qxi, qx
c            write(6,*)'qbz=',qbz(:,ibz)
c            write(6,*)'qx =',qx
c            write(6,*)'qxi=',qxi
!! inversion table:
!!             ibz<---from qbzrid
!! we can skip such ibz.
 3012       continue
          enddo
          enddo
          immx(ibz)=imm
          if(imm>immxx) immxx=imm
        enddo

        if(debug) write(6,*)'eibz aaaaaaa'
c        call cputid(0)

!! inequivalent points.
          allocate(mqbz(nqbz))
          mqbz=0

!! may2015 there zero-cleare is moved from next ix loop. (bug fix). correct?
          nwgt(:,iq)=0
          igx(:,:,iq) =999991
          igxt(:,:,iq)=999992
!!
          do ix=1,ntimer,-2 !try no time reversal for ix=1
             if(ix== 1) inotime=.true.
             if(ix==-1) inotime=.false.
             do ibz=1,nqbz
               if(debug) write(6,"('1111 ix,ibz',i3,i5,i5)")ix,ibz
               if(mqbz(ibz)==1) cycle
               do im=1,immx(ibz) 
                   if(inotime) then
                      if(igxt_(im,ibz) == -1) cycle
                   elseif(inotime) then
                      if(igxt_(im,ibz) ==  1) cycle
                   endif
                   kkk3 = nxx(:,ibz,im)
                   if(debug) write(6,*)'uuuuu ibz,im kkk3=',ibz,im,kkk3
                   call tabkk(kkk3(1), kk1,nkey(1), ik1)
                   call tabkk(kkk3(2), kk2,nkey(2), ik2)
                   call tabkk(kkk3(3), kk3,nkey(3), ik3)
                   if(ik1==-999999.or.ik2==-999999.or.ik3==-999999) cycle !! june2016
                   if(debug) write(6,*)'uuuuu ik=',ik1,ik2,ik3
                   ibzx = iqkkk(ik1,ik2,ik3)
                   if(mqbz(ibzx)==1) then
                     cycle
                   else  
                     mqbz(ibzx)=1
                     nwgt(ibz,iq)= nwgt(ibz,iq) + 1
                     igx (nwgt(ibz,iq),ibz,iq) = igx_(im,ibz)  
                     igxt(nwgt(ibz,iq),ibz,iq) = igxt_(im,ibz)  
                   endif   
               enddo
             enddo
             if(debug) then
             do ibz=1,nqbz
               write(6,"('2222 ix,ibz mqbz nqbz=',i2,i3,i5,i5)")ix,ibz,mqbz(ibz),nwgt(ibz,iq)
             enddo
             endif
             if(sum(nwgt(:,iq))==nqbz) then
                if(ix==1) eibzsym(:,-1,iq)=0
                goto 1202
             endif   
             write(6,*)' mptauo: we use time reversal'
          enddo
          call rx( 'eibzgen: bug. nwgt sum is not nqbz')
 1202     continue

C cccccccccccccccccccccccccccccccccc
C         do ibz=1,200 !nqbz
C           if(nwgt(ibz,iq)/=0) then
C             write(6,"('yyy0: ',i8,2x,25(i3,i2))") ibz,(igx(i,ibz,iq),igxt(i,ibz,iq),i=1,nwgt(ibz,iq))
C           endif
C         enddo
C cccccccccccccccccccccccccccccccccc

c          write(6,*)'pppppp nwgt(ibz,iq),nqbz=',iq,sum(nwgt(:,iq)),nqbz
c$$$!! === Search other ibz generated from ibz. ===
c$$$!!   ibznxx(ibz,imm) means "ibz rotated by igx(imm,ibz,iq),igxt(imm,ibz,iq)"
c$$$        allocate(mqbz(nqbz))
c$$$        do ibzx=1,nqbz
c$$$          do ibz=1,ibzx !search equivalent (ibzx,1) with first appeard (ibz,im) .
c$$$          do im=1,immx(ibz) 
c$$$             if(sum(abs(nxx(:,ibzx,1) - nxx(:,ibz,im))) <2) then
c$$$               nwgt(ibz,iq) = nwgt(ibz,iq) +  1
c$$$               igx (nwgt(ibz,iq),ibz,iq) = igx_(im,ibz)  
c$$$               igxt(nwgt(ibz,iq),ibz,iq) = igxt_(im,ibz)  
c$$$               goto 3312
c$$$             endif  
c$$$          enddo
c$$$          enddo 
c$$$ 3312     continue
c$$$        enddo
        deallocate(igx_,igxt_)
c        write(6,*)'eibz bbbbbb'
c        call cputid(0)
!
        neibzx=0
        do ibz=1,nqbz
           if(nwgt(ibz,iq) /=0) then
              neibzx=neibzx+1
              qeibz(:,neibzx,iq)=qbz(:,ibz)
              write(6,"('   ibz qeibz =',i7,3f10.5, ' nwgt=',i7)")ibz,qeibz(:,neibzx,iq),nwgt(ibz,iq)
           endif
        enddo   
        neibz(iq) = neibzx !number of eibz for iq
        if( iprintx ) then
          write(6,"('iq=',i8,' # of EIBZ: Full=',i8, 
     &      ' Used(TimeR 1 or -1)=',i3,'=',i2,'+',i2,' neibz= ',i7)")iq,sum(abs(eallow(:,:))),
     &      sum(eibzsym(:,:,iq)),sum(eibzsym(:,1,iq)),sum(eibzsym(:,-1,iq)),neibz(iq)
c          write(6,"('eibz: iq neibz nqbz= ',i3,3f11.5,3i7)") iq,q,neibz(iq),nqbz
        endif
        deallocate(nxx,immx,mqbz)
      enddo
!!
      nwgtsum = sum(nwgt(1:nqbz,iqxini:iqxend))
      timerout=timereversal
      if(minval(igxt)==1) timerout=.false.

      if(iprintx) then
        write(6,"(' nqbz,  sum(neibz(iq)), sum(ngwt)=sum(nqbz)= ',10i7)")
     &  nqbz,sum(neibz(iqxini:iqxend)), nwgtsum, nqbz*(iqxend-iqxini+1) 
      endif  
c      stop 'test stop -------- end of eibzgen ----------'
      end subroutine eibzgen

      subroutine tabkk(kkin, kktable,n, nout)
      integer, intent(in) :: n,kkin, kktable(n)
      integer, intent(out) :: nout
      integer:: i,mm,i1,i2
      i1=1
      i2=n
      if(kkin==kktable(1)) then
         nout=1
         return
      elseif(kkin==kktable(n)) then
         nout=n
         return
      endif
      do i=1,n
         mm=(i1+i2)/2
         if(kkin==kktable(mm)) then
            nout=mm
            return
         elseif(kkin>kktable(mm)) then
            i1=mm
         else
            i2=mm
         endif
      enddo
c$$$      do i=1,n                                                                                  
c$$$         if(kkin==kktable(i)) then                                                              
c$$$            nout=i                                                                              
c$$$            return                                                                              
c$$$         endif                                                                                  
c$$$      enddo                                                                                     
!jun2016
c      call rx( 'takk: error')
c      write(6,*) i1,i2,kkin
c      write(6,*) kktable(i1),kktable(i2)
      nout=-999999
      end subroutine tabkk

      subroutine mptauof(symops,ng,plat,nbas,bas,
     &             iclass,miat,tiat,invg,delta)
C- Mapping of each atom by points group operations------------------c
C Modded by okuda 1994 March.
C Simplified by kotani 1994 7/31.
Ci  Input
Ci     symops(1,ng),ng,plat,nbas,bas(3,nbas)
Ci     iclass(nbas); denote class for each atom
Co  Output
Co    miat(ibas  ,ig); ibas-th atom is mapped to miat-th atom, by the ig-th
Co    points group operation.  Origin is (0,0,0).
Co    tiat(k,ibas,ig);
Co    delta : shifting vector for non-symmorphic group.
Co            r' = matmul (am, r) + delta
Cr  Remarks
Cr
Cr (1) The ibas-th atom (position at bas(k,ibas) ) is mapped to
Cr
Cr    bas( k,miat(ibas,ig) )+ tiat(k,ibas,ig), k=1~3.
Cr
Cr (2) tiat= unit translation
Cr
C--------------------------------------------------------------------
      implicit none
      integer ng,nbas, miat(nbas,ng),iclass(nbas),invg(ng),
     &        nbmx, nsymx, ig,igd,i,j,ibas,mi,i1,i2,i3
      double precision SYMOPS(9,ng),plat(3,3),
     &                 tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas),
     &                 tr1,tr2,tr3,ep, dd1,dd2,dd3,t1,t2,t3
      integer  iprint
      integer ires(3, nbas, ng)
      integer(4):: ib1,ib2
      real(8) ::tran(3),delta(3,ng)
      data ep/1.0d-3/
c      data ep/1.0d-7/
c
      write(6,*)' mptauof: search miat tiat for wave function rotation'

      do 10 ig=1,ng
        do igd=1,ng
c seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ).le.ep.and.
     &        abs( symops(2,ig)-symops(4,igd) ).le.ep.and.
     &        abs( symops(3,ig)-symops(7,igd) ).le.ep.and.
     &        abs( symops(4,ig)-symops(2,igd) ).le.ep.and.
     &        abs( symops(5,ig)-symops(5,igd) ).le.ep.and.
     &        abs( symops(6,ig)-symops(8,igd) ).le.ep.and.
     &        abs( symops(7,ig)-symops(3,igd) ).le.ep.and.
     &        abs( symops(8,ig)-symops(6,igd) ).le.ep.and.
     &        abs( symops(9,ig)-symops(9,igd) ).le.ep  ) then
            invg(ig)=igd
            goto 16
          endif
        end do
 16     continue
c
        if(iprint() .ge.40) then
          write(6,*)' '
          write(6,*)' '
          write(6,*)' **** group ops no. ig (igd)= ', ig, invg(ig)
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
 1731     format (' ',3f9.4)
        endif

        do i=1,3
          do j=1,3
            am(i,j)=symops(i+3*(j-1),ig)
          end do
        end do
c
c trial shift vector tran
        do 20 ib1=1,nbas
        do 20 ib2=1,nbas
          tran =  bas(:,ib2)  - matmul(am,bas(:,ib1))
c
          do 30 ibas=1,nbas
            b1=am(1,1)*bas(1,ibas)
     &        +am(1,2)*bas(2,ibas)+am(1,3)*bas(3,ibas)
     &        +tran(1)
c     .        +( tr1*plat(1,1)+tr2*plat(1,2)+tr3*plat(1,3) )
            b2=am(2,1)*bas(1,ibas)
     &        +am(2,2)*bas(2,ibas)+am(2,3)*bas(3,ibas)
     &        +tran(2)
c     .        +( tr1*plat(2,1)+tr2*plat(2,2)+tr3*plat(2,3) )
            b3=am(3,1)*bas(1,ibas)
     &        +am(3,2)*bas(2,ibas)+am(3,3)*bas(3,ibas)
     &        +tran(3)
c     .        +( tr1*plat(3,1)+tr2*plat(3,2)+tr3*plat(3,3) )
c
            do 40 mi=1,nbas
              if( iclass(mi).ne.iclass(ibas) ) go to 40

              do 50 i1=-3,3
              do 50 i2=-3,3
              do 50 i3=-3,3
                dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )

                t1 = b1 - (bas(1,mi)+dd1)
                t2 = b2 - (bas(2,mi)+dd2)
                t3 = b3 - (bas(3,mi)+dd3)
                if(abs(t1).le.ep.and.abs(t2).le.ep.and.
     &             abs(t3).le.ep) go to 60
   50         continue
   40       continue
c seach failed, Not found mi and dd1. Try next (tr).
            goto 20

   60       continue
            miat(ibas,ig)  = mi
            tiat(1,ibas,ig)= dd1
            tiat(2,ibas,ig)= dd2
            tiat(3,ibas,ig)= dd3
            ires(1,ibas,ig)= i1
            ires(2,ibas,ig)= i2
            ires(3,ibas,ig)= i3
c
   30     continue
c When the do-30 loop has been completed, we get out of do-20 loop
          goto 21
   20   continue
        call rx( 'mptauo2: Can not find miat and tiat')
   21   continue
        delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000
c- have gotten the translation-> check write --------------------
        if(iprint().ge.40) then
c          write(6,4658)tr1,tr2,tr3
          write(6,4658)tran
 4658     format('  Obtained translation operation=',3d12.4)
          do 123  ibas=1,nbas
            write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig),
     &    tiat(2,ibas,ig), tiat(3,ibas,ig),
     &    ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
  150       format(' ibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
  123     continue
        endif
c---------------------------------------------------
   10 continue
      end



