      real(8) function xqconv(x)
!!! x --> xqconv : uniform to non uniform mesh converter
!! x is [0,1] --> xqcon = [0,1]
!! x can be -1 <= x =< 1
      use m_keyvalue,only: getkeyvalue
      real(8),intent(in)::x
      real(8),parameter:: pi=3.1415926535897932d0
      real(8),save:: adiv,bdiv
      logical,save:: oncew=.true.
      logical:: ggg
!! BZ division setting.
      if(oncew) then
         inquire(file='GWinput',exist=ggg)
         if(.not.ggg) then
            adiv=1d0
         else
            call getkeyvalue("GWinput","BZadiv",adiv,default=1d0)
         endif   
         write(6,"('BZadiv= ',f6.3)") adiv
         bdiv = (1d0 -adiv)/2d0
         oncew=.false.
      endif
!!   adiv+2b=1
      xqconv           = adiv*x + bdiv*(1-cos(pi*x))
      if(x>1d0) xqconv = adiv*x + bdiv*(1-cos(pi*(x-1d0))+2d0)
      if(x<0d0) xqconv = adiv*x + bdiv*(1-cos(pi*(x+1d0))-2d0)
c      print *,'x xqconv=',x,xqconv
      end function xqconv

      subroutine xconvv(xin,xout)
!!! xin is converted to xout (uniformmesh to non-uniform mesh).
      real(8)::xin(3),xout(3),xqconv
      integer:: i
      do i=1,3
        xout(i)=xqconv(xin(i))
      enddo
      end subroutine xconvv

!! ------------------------------------------------------------------------------------------------
      module m_get_bzdata1
!! In addition to these variables, this also write mtet file (search ifmtet) for mtet mode.
!! bzmesh, tet are included in this routine.
c------------------------------------------------------------
      implicit none  
!! all are output by getbzdata1
      integer,protected::             nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm 
      real(8),allocatable,protected:: qbz(:,:),wbz(:),qibz(:,:),wibz(:),qbzw(:,:) 
      real(8),protected::             dq_(3)!,qbasmc(3,3) !,dq_bzcase2(3),
      integer,allocatable,protected:: idtetf(:,:),ib1bz(:),idteti(:,:),irk(:,:),nstar(:),nstbz(:)
      real(8),allocatable,protected:: qbzm(:,:),qbzwm(:,:)
c------------------------------------------------------------
      contains

      subroutine getbzdata1(qlat,nnn, !n1q,n2q,n3q
     &    symops,ngrp,tetrai,tetraf,mtet,gammacellctrl) 
      use m_keyvalue,only: getkeyvalue
!! all arguments are inputs. getbzdata1 returns all variables in the module m_get_bzdata1
      implicit none
      logical :: tetrai,tetraf,multet,qbzreg
      integer :: ngrp,nnn(3), mtet(3),n1qm,n2qm,n3qm,nnnx(3), !n1q,n2q,n3q,
     &    itet,ix,im,ifiqmtet,iq,icase,n1qtet,n2qtet,n3qtet,IMC(0:1,0:1,0:1)
      real(8)  plat(3,3),qlat(3,3),ginv(3,3),qlattet(3,3),
     &    symops(3,3,ngrp),qc(3,0:3),qb(3,3),qb1(3,3)
      integer,allocatable:: IPQ(:,:,:),iw1(:),indexkw(:,:,:),indexk(:,:,:)
      real(8),allocatable:: qcm(:,:,:),wtet(:,:,:)  !,qibz_r(:,:) !,qbzz(:,:,:,:)
      logical:: skipgammacell
      integer:: gammacellctrl,iprint
      integer,allocatable::  idtetfm(:,:,:),ib1bzm(:),iwgt(:) !multipled tetrahedron.
      integer:: n,iqbz,intq(3),iqx,iqy,iqz,iqbzx,verbose !,n_index_qbzm
      real(8):: qout(3),qx(3),kvec(3,3),det33,  kkvec(3,0:3)
      real(8):: qbzshift(3)
      integer::  nmtet,nqbzwm,ntetfm,ifmtet       !,nqibz_r
      integer:: iqi,igrp,iccc,nnng(3),icasetet,nadd=0,iq4(4)
      real(8)::qmic(3,3),weight,x1,x2,x3,xqconv,x1i,x2i,x3i,wsum
      integer::kount,i1,i2,i3,n2,nnnv(3),ig,ir,j1,j2,j3,jj(3),k,m,i1x,i2x,i3x,ic,k1,k2,k3,kountw,ndx,
     &  ibtr(3,3),kcut(3,4,6),ngcell,i,ii,j,igamma
      real(8):: xvec(3),xvecc(3),xvecs(3),xvece(3),vv(3),xvv(3),xv(3),diff(3),diff2(3),swgt,v1(3),wfac
      real(8):: tolq=1d-8
      logical,allocatable::usediqig(:,:)
!! icase=1
c-------------------------------------------------------
      write(6,"('getbzdata1: n1n2n3=',3i5)") nnn(1:3)
      call minv33tp (qlat,plat)  !qlat --> plat 
      qbzshift=0d0
c      if(icase==2) qbzshift=0.5d0
c      if(icase>2) call rx( 'mkqg: wrong bzcase')
      if(gammacellctrl==2) then
        nadd=1  ! nddd=1 give end point of BZ (for gammacell case). See nqbz
        qbzshift=-1d0/2d0 ! gamma centerered.
              !this was (/n1q/2d0,n2q/2d0,n3q/2d0/) !this shift is so that gamma point is centerd.
      endif
      nnnv = nnn+nadd
c      wfac=dble(product(nnn))/product(nnnv)
      nqbz = product(nnnv) !(nnn(1)+nadd)*(nnn(2)+nadd)*(nnn(3)+nadd) 
      allocate(qbz(3,nqbz),wbz(nqbz))!,nstbz(nqbz))
      allocate(nstbz(nqbz)
     & ,indexk (0:nnnv(1)-1,0:nnnv(2)-1,0:nnnv(3)-1)
     & ,indexkw(0:nnn(1),0:nnn(2),0:nnn(3)))
      nstbz(1:nqbz)=0
      print *,'nnn =',nnn
      print *,'nnnv=',nnnv

c      call genqbz (qlat,nnn(1),nnn(2),nnn(3), !icase,
c     o     qbz,wbz,nstbz, nadd,qbzshift)
c      qmic(:,1)= qlat(:,1)/dble(nnn(1))
c      qmic(:,2)= qlat(:,2)/dble(nnn(2))
c      qmic(:,3)= qlat(:,3)/dble(nnn(3))
c      allocate(qbzz(3,nnnv(1),nnnv(2),nnnv(3)),qbzw(nnnv(1)+1,nnnv(2)+1,nnnv(3)+1))) 
c      weight = 1d0/dble(nnn1t*nnn2t*nnn3t)

!! Generage qbz
      wsum=0d0
      kount  = 0
!! In principle following do loop ordirng do not affect to the final result.
!! However, because of degeneracy (in the tetrahedron method), we have a slight effect.
      do i3 = 1,nnnv(3)
      do i2 = 1,nnnv(2)
      do i1 = 1,nnnv(1)
c      do i1 = 1,nnnv(1)
c      do i2 = 1,nnnv(2)
c      do i3 = 1,nnnv(3)
        xvec = (/I1-1,I2-1,I3-1/)/dble(nnn) + qbzshift
        call xconvv(xvec,           xvecc) !cell center
        kount=  kount + 1
        qbz(:,kount) = matmul(qlat,xvecc) !x1*qmic(:,1) + x2*qmic(:,2) + x3*qmic(:,3)
        if(i1<=nnn(1).and.i2<=nnn(2).and.i3<=nnn(3)) then
          call xconvv(xvec-0.5d0/nnn, xvecs) !cell start
          call xconvv(xvec+0.5d0/nnn, xvece) !cell end
c        print *,'i1i2i3=',i1,i2,i3
c        write(6,"('xvec=',3i5,3d13.5)")i1,i2,i3,xvec
          weight = product(xvece-xvecs)
c        qbzz(:,i1,i2,i3)= qbz(:,kount)
          wbz(kount) =  weight
          wsum = wsum + weight
        else
          wbz(kount) =  0d0
        endif
        indexk(i1-1,i2-1,i3-1)= kount
c         write(6,"(' iiiiiiii iq=',i8,' q=',3f17.13)")  kount, qbz(:,kount)
      enddo
      enddo
      enddo
      write(6,*) 'wsum=',wsum
      if( abs(wsum - 1d0) >tolq) call rx( 'getbzdata1: wrong wsum check')

!! Generate qbzw 
      kountw=0
      nqbzw = (nnn(1)+1)*(nnn(2)+1)*(nnn(3)+1)
      allocate(ib1bz(nqbzw), qbzw(3,nqbzw) )
      do i3 = 1,nnn(3)+1
      do i2 = 1,nnn(2)+1
      do i1 = 1,nnn(1)+1
c      do i1 = 1,nnn(1)+1 
c      do i2 = 1,nnn(2)+1
c      do i3 = 1,nnn(3)+1
        kountw = kountw+1
        indexkw(i1-1,i2-1,i3-1)= kountw
        xvec = (/I1-1,I2-1,I3-1/)/dble(nnn) + qbzshift
        call xconvv(xvec,           xvecc) !cell center
        qbzw(:,kountw) = matmul(qlat,xvecc) 
        if(nadd==0) ib1bz(kountw) = indexk(mod(i1-1,nnn(1)), mod(i2-1,nnn(2)), mod(i3-1,nnn(3)))
        if(nadd==1) ib1bz(kountw) = indexk(i1-1, i2-1, i3-1)
      enddo
      enddo
      enddo

c!x1*qmic(:,1) + x2*qmic(:,2) + x3*qmic(:,3)
c        write(6,"('xxxy=',3f9.4,2x,3f9.4,2x,3f9.4)")xvec,xvecc,weight

c
c        do iq=1,nqbz
c          write(6,"(' iq qbz nstbz=',i5,3f9.4,i5)")iq,qbz(:,iq),nstbz(iq)
c        enddo
c
c      if(icase/=1) then !get qibz_r, and nqibz_r
c        allocate(qibz_r(3,nqbz),ipq(nqbz),wibz(nqbz))
c        call bzmesh(1,plat,qbasmc,nnn(1),n2q,n3q,symops,ngrp,ipq,qibz_r
c     &      ,wibz,nqibz_r,nqbz,nadd,qbzshift) !Make q-points in IBZ.
c        deallocate(ipq,wibz)
c      else
c        nqibz_r=0
c      endif

c      nqibz_r=0
      allocate(qibz(3,nqbz),ipq(nnnv(1),nnnv(2),nnnv(3)),wibz(nqbz),iwgt(nqbz),usediqig(nqbz,ngrp))
      !nqibz instead of nqbz is enough (but nqibz is calculated from now on)

!! 
      WRITE(6,"(' getbzdata1 : ',5X,'plat',31X,'qlat')")
      DO K = 1, 3
        WRITE(6,"(3F10.5,5X,3F10.5)") (plat(M,K),M=1,3),(qlat(M,K),M=1,3)
      enddo
      ipq   = 0
      nqibz = 0
      wibz  = 0d0
      iwgt  = 0
      usediqig=.false.
!! Generate qibz. Gamma point first (by igamma loop).
      do 229 igamma=0,1 !igamma=0 is for gamma point.
      do 20 i3 = 1,nnnv(3)
      do 20 i2 = 1,nnnv(2)
      do 20 i1 = 1,nnnv(1)
c      do 20 i1 = 1,nnnv(1)
c      do 20 i2 = 1,nnnv(2)
c      do 20 i3 = 1,nnnv(3)
        xvec = (/I1-1,I2-1,I3-1/) / dble(nnn)
        call xconvv(xvec+qbzshift,           xvecc)
        if(sum(abs(xvecc))<tolq ) then
          if(igamma==1) cycle !skip Gamma since it is already in it.
        else
          if(igamma==0) cycle !get Gamma point only.
        endif

        vv = matmul(qlat, xvecc)
        if(i1<=nnn(1).and.i2<=nnn(2).and.i3<=nnn(3)) then
          call xconvv(xvec-0.5d0/nnn, xvecs)
          call xconvv(xvec+0.5d0/nnn, xvece)
          weight= 2d0* product(xvece-xvecs) 
        else
          weight=0d0
        endif
!!   Check vv can be mapped to q(iq,ig).
        do 22 iq = 1, nqibz
        do 22 ig = 1, ngrp
          if(usediqig(iq,ig)) cycle
          v1  = matmul( symops(:,:,ig), vv) 
          diff= matmul( qibz(:,iq)-v1, plat) ! [0,1] range because plat= transpose of qlat^-1 
c          print *,'fffffff ssss=',symops(:,:,ig)
          if(nadd==0) then !periodic case
            diff= diff - nint(diff)
          endif
          if(sum(abs(diff))<tolq) then !!we found qibz(:,iq) is mapped to vv.
            usediqig(iq,ig)=.true.
            iqx=iq
            goto 23
          endif
 22     continue
!!   New iqbz added
        nqibz = nqibz+1
        qibz(:,nqibz) = vv
        iqx = nqibz
 23     continue
        ipq(i1,i2,i3) = iqx
        iwgt(iqx)     = iwgt(iqx) + 1
        wibz(iqx)     = wibz(iqx) + weight
 20   continue
 229  continue
      SWGT = sum(wibz(1:nqibz))
c      if(gammacellctrl/=2) then
         if (dabs(swgt-2d0)> tolq) call rx( 'getbzdata: error in weights swgt')
c        write(6,"(/' getbzdata: sum swgt iwgt=',f12.8)") swgt
      if (sum(iwgt(1:nqibz))-product(nnnv)/=0) call rx('getbzdata: sum(iwgt(1:nqibz))-product(nnn)/=0')
c      endif
      write(6,"(/' getbzdata: nqibz sum.iwgt nqbz =',3i8)") nqibz,product(nnn)
      write(6,744) nqibz,product(nnnv)
  744 FORMAT(' getbzdata: ',I4,' irreducible k - points. from ',i5,' (',3i3,' )')
      write(6,"(13x,'Qx',8x,'Qy',8x,'Qz',6x,'Multiplicity    Weight')")
      do iq = 1, nqibz
        write(6,661) iq,qibz(1:3,iq),iwgt(iq),wibz(iq)
      enddo
  661 format(i5,2x,3f10.4,i10,f16.6)
      write(6,"(' nqibz= ',i8)")nqibz


!! Rotation checker
      do 30 i3 = 1,nnnv(3)
      do 30 i2 = 1,nnnv(2)
      do 30 i1 = 1,nnnv(1)
        xvec = (/I1-1,I2-1,I3-1/) / dble(nnn)
        call xconvv(xvec+qbzshift,           xvecc)
        vv = matmul(qlat, xvecc)
        do k  = 1,nqibz
        do ir = 1,ngrp
          diff = matmul(symops(:,:,ir),qibz(:,k)) - vv
          if(nadd==0) then
           call rangedq(matmul(diff,plat), diff2)
          else
            diff2=diff
          endif
          if(sum(abs(diff2))< tolq) then
             if(verbose()>50) then
             write(6,"(' i1i2i3= ',3i3,' q qibz k=',2x,3f7.3,2x,3f7.3,2i5)")
     &       i1,i2,i3, vv, qibz(:,k),k,k-ipq(i1,i2,i3) ! matmul(qp(:,NQx),rbas)*8, matmul(v, rbas)*8
             endif  
             goto 1022
           endif
        enddo
        enddo
        call rx( 'getbzdata1: not find irotk; it may require accurate symmetry.')
 1022   continue
 30   continue


c      call bzmesh(plat,qbasmc,nnn(1),nnn(2),nnn(3),symops,ngrp,ipq,qibz !icase,
c     &    ,wibz,nqibz,nqbz,nadd,qbzshift) !Make q-points in IBZ.

      write(6,"(/' nkstar: ----- nqbz nqibz=',2i8, ' -----')") nqbz,nqibz
      call minv33(qlat,ginv)
      allocate(nstar(nqibz), irk(nqibz,ngrp),iw1(nqbz))
      call nkstar  (qibz,qbz,symops,ginv, !ginv=plat
     d              nqibz,nqbz,ngrp,gammacellctrl,
     o              nstar,irk )

!! irk check. irk(iqi,igrp) specify
      iccc=0
      do iqi=1,nqibz
      do igrp=1,ngrp
         if(irk(iqi,igrp)/=0) then
c            write(6,*) 'iq igrp irk=',iqi,igrp,irk(iqi,igrp)
            iccc=iccc+1
         endif
      enddo
      enddo
      if(nqbz/=iccc) then
        write(6,*)'getbzdata1: total count=',iccc,nqbz
        call rx('getbzdata1: iccc check.error')
      endif
      ntetf=-1
      nteti=-1

!! Full bz divided into tetrahedron
      if(tetraf) then 
cccccccccccccccccccccccccc
c        goto 3331
cccccccccccccccccccccccccc
        skipgammacell=.false.
        if(gammacellctrl==1) skipgammacell=.true.
        write(6,"(/' tetfbzf ntetf =',2i8, ' -----')") ntetf
        call getkeyvalue("GWinput","ngcell",ngcell,default=1)
        ntetf = 6*nqbz
        allocate( idtetf(0:3,ntetf))
        ntetf = 0
        qb(:,1) = qlat(:,1)/nnn(1)
        qb(:,2) = qlat(:,2)/nnn(2)
        qb(:,3) = qlat(:,3)/nnn(3)
        CALL CCUTUP(QB,QB1,IBTR,KCUT) !This is from LMTO-3. Need to clean up ccutup
        DO 120  I3 = 1, nnn(3)
        DO 121  I2 = 1, nnn(2)
        DO 122  I1 = 1, nnn(1)
c          qb(:,1) = qbzw(:,i1+1)-qbzw(:,i1)
c          qb(:,2) = qbzw(:,i2+1)-qbzw(:,i2)
c          qb(:,3) = qbzw(:,i3+1)-qbzw(:,i3)
c          CALL CCUTUP(QB,QB1,IBTR,KCUT) !This is from LMTO-3. Need to clean up ccutup
!!   skip gamma cell option
          if(skipgammacell) then
           i1x = i1-nnn(1)
           i2x = i2-nnn(2)
           i3x = i3-nnn(3)
           if(i1<nnn(1)/2) i1x=i1
           if(i2<nnn(2)/2) i2x=i2
           if(i3<nnn(3)/2) i3x=i3
           ndx=0
           if(qbzreg()) ndx=1
           if( i1x<ngcell .and. -ngcell< i1x-ndx) then
           if( i2x<ngcell .and. -ngcell< i2x-ndx) then
           if( i3x<ngcell .and. -ngcell< i3x-ndx) then
             kount= indexkw(i1-1,i2-1,i3-1)
             write(6,"('qqqqx qbzw=',3i3, 3f9.5)") i1x,i2x,i3x, matmul(qbzw(1:3,kount),plat)
             cycle
           endif
           endif
           endif
          endif 
!!   set UP IDENTIFIERS AT 8 CORNERS OF MICROCELL ------
          IMC(0:1,0:1,0:1)  = indexkw(I1-1:I1,I2-1:I2,I3-1:I3)
c$$$          do K1 = 0, 1
c$$$            J1 = I1 -1 + K1
c$$$          do K2 = 0, 1
c$$$            J2 = I2 -1 + K2
c$$$          do K3 = 0, 1
c$$$            J3 = I3 -1 + K3
c$$$            IMC(K1,K2,K3)  = indexkw(J1,J2,J3)
c$$$          enddo
c$$$          enddo
c$$$          enddo
!!   LOOP OVER TETRAHEDRA 
        do 10 ITET = 1, 6 
          do IC = 1, 4
            K1 = KCUT(1,IC,ITET)
            K2 = KCUT(2,IC,ITET)
            K3 = KCUT(3,IC,ITET)
c            print *,'kkkkk',ic,itet,k1,k2,k3
            IQ4(IC) = IMC(K1,K2,K3)
          enddo  
          ntetf=ntetf+1
          idtetf(0:3,ntetf) = iq4(1:4)
   10   continue
c        write(6,"('ntetf=',3i3,2x,i6,2x,4i8)")i1,i2,i3,ntetf,idtetf(0:3,1)
 122   CONTINUE
 121   CONTINUE
 120   CONTINUE 
       write(6, "(1x,'getbzdata1 TETFBZF: ',2i8 )")  ntetf, 6*nqbz
c!! apr2016
c$$$ 3331  continue
c$$$        skipgammacell=.false.
c$$$        if(gammacellctrl==1) skipgammacell=.true.
c$$$        ntetf = 6*nqbz
c$$$        allocate( idtetf(0:3,ntetf))
c$$$        write(6,"(/' tetfbzf ntetf =',2i8, ' -----')") ntetf
c$$$c        call getkeyvalue("GWinput","ngcell",ngcell,default=1)
c$$$        call TETFBZF(qlat, nnn(1),nnn(2),nnn(3),qbz,nqbz,
c$$$     i      skipgammacell,qbzshift,nadd,
c$$$     o      idtetf,qbzw,ib1bz,ntetf)
      else
        ntetf=0
      endif

!! Teterahedron Irreducible case.
      if(tetrai) then  !ibz is divided into tetrahedron
        allocate( idteti(0:4,6*nqbz)) !0 is for number of equilvanent tetrahedron
        print *,' goto tetirr---' !this should work also for icase==2
        qb(:,1) = qlat(:,1)/nnn(1)
        qb(:,2) = qlat(:,2)/nnn(2)
        qb(:,3) = qlat(:,3)/nnn(3)
ccccccccccccccccccccccccc
c        call tetirr(qb,nnn(1),nnn(2),nnn(3),ipq,nqibz,
c     o              nteti,idteti)
c        goto 3333
ccccccccccccccccccccccccc
        CALL CCUTUP(qb,QB1,IBTR,KCUT)
        nteti = 0
!!   START LOOPING OVER MICROCELLS ---------
        DO 220  I3 = 1, nnn(3)
        DO 221  I2 = 1, nnn(2)
        DO 222  I1 = 1, nnn(1)
!!   SET UP IDENTIFIERS AT 8 CORNERS OF MICROCELL ------
          do K1 = 0, 1
            J1 = MOD(I1 + K1 -1,nnn(1)) + 1
          do K2 = 0, 1
            J2 = MOD(I2 + K2 -1,nnn(2)) + 1
          do K3 = 0, 1
            J3 = MOD(I3 + K3 -1,nnn(3)) + 1
            IMC(K1,K2,K3)  = IPQ(J1,J2,J3)
          enddo
          enddo
          enddo
!!   LOOP OVER TETRAHEDRA --------------
          do 110  ITET = 1, 6
            do IC = 1, 4
              K1 = KCUT(1,IC,ITET)
              K2 = KCUT(2,IC,ITET)
              K3 = KCUT(3,IC,ITET)
              IQ4(IC) = IMC(K1,K2,K3)
            enddo          
!!   order the identifiers -----
            do J=1,3
            do I=1,4-J
              IF(IQ4(I)>IQ4(I+1))THEN
                II=IQ4(I)
                IQ4(I)=IQ4(I+1)
                IQ4(I+1)=II
              ENDIF
            enddo
            enddo
            do n = 1, nteti
              if ( sum(abs(idteti(1:4,n)-iq4(1:4)))==0) then
                idteti(0,n) = idteti(0,n) + 1
                goto 101
              endif
            enddo
            nteti = nteti+1
            idteti(1:4,nteti) = iq4(1:4)
            idteti(0,nteti)=1
c            write(6,"(' nnnnn i1i2i3=',3i5,2x,10i3)") i1,i2,i3,idteti(:,nteti)
 101        CONTINUE
 110      CONTINUE
 222     CONTINUE
 221     CONTINUE
 220     CONTINUE
        write(6,"(1x,'getbzdata1 TETIRR: nteti=',i8,'  sum(idteti)= full tetrahedra=',2i8)")
     &    nteti, sum(idteti(0,1:nteti)),6*nqbz
c      stop 'xxxxxxxxxxxxxxxx test '
c 3333   continue
      else
        nteti=0
      endif
      deallocate(ipq)

c$$$!! for multipled tetrahedron. not used now...
c$$$      if(sum(abs(mtet))/=3) then
c$$$        icase=1
c$$$        print *, 'multitet mode: mtet=',mtet
c$$$        n1qm = mtet(1)*nnn(1)
c$$$        n2qm = mtet(2)*nnn(2)
c$$$        n3qm = mtet(3)*nnn(3)
c$$$        if(icase==2) then
c$$$          n1qm=n1qm+1
c$$$          n2qm=n2qm+1
c$$$          n3qm=n3qm+1
c$$$        endif
c$$$        nmtet   = mtet(1)*mtet(2)*mtet(3)
c$$$        nqbzm   = nmtet * nqbz
c$$$        nqbzwm  = (n1qm+1)* (n2qm+1)* (n3qm+1)
c$$$        ntetfm  = ntetf * nmtet
c$$$        allocate( 
c$$$     &       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
c$$$        ! Index for tetrahedron;    qbzmw(idtetfm) gives extended q vector.
c$$$     &       ib1bzm(nqbzwm), qbzm(3,nqbzm) )
c$$$        ! qbzm(1:3,ib1bz(iq)) gives q vector within the 1st bz.
c$$$! The datas idetetfm, qbzmw, ib1bzm, nmete, ntetf, nqbzm,
c$$$!     eigen(nband,nqbzm) are required for the multiply-devided tetrahedron method.
c$$$!
c$$$        allocate( qcm(1:3,0:3, nmtet), wtet(0:3, nmtet, ntetf), 
c$$$     &          indexkw(0:n1qm,0:n2qm,0:n3qm) )
c$$$        QB(:,1) = qlat(:,1)/(nnn(1)*mtet(1))
c$$$        QB(:,2) = qlat(:,2)/(nnn(2)*mtet(2))
c$$$        QB(:,3) = qlat(:,3)/(nnn(3)*mtet(3))
c$$$        call minv33(qb,qbi)
c$$$        call qwider(icase,qb, mtet(1)*nnn(1),mtet(2)*nnn(2),mtet(3)*nnn(3),
c$$$     &                       n1qm,n2qm,n3qm, ib1bzm, qbzwm, qbzm, indexkw )
c$$$        do itet=1, ntetf  !--- idtetfm(0:3, nmtet, ntetf), nmtet divition.
c$$$          do ix=0,3
c$$$            qc(:,ix)= qbzw(:, idtetf(ix,itet))
c$$$          enddo
c$$$ccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c        do ix = 1,3
c$$$c        kvec(1:3,ix) = qc(1:3,ix) - qc(1:3,0)
c$$$c        enddo
c$$$c        write(6,"('itet vol=',i5,d13.5)") itet,abs(det33(kvec(1:3,1:3))/6d0)
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$          call tetdevide(qc,qcm, wtet(:,:,itet),mtet(1),mtet(2),mtet(3))   !qc ---> qcm
c$$$          do im=1,nmtet ! Index for micro-devided tetrahedron
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c        do ix = 1,3
c$$$c        kvec(1:3,ix) = qcm(1:3,ix,im) - qcm(1:3,0,im)
c$$$c        enddo
c$$$c        write(6,"('itet im vol=',2i5,d13.5)") itet,im,abs(det33(kvec(1:3,1:3))/6d0)
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$            do ix=0,3     ! its four corners
c$$$              !qcm ---> idtetfm(ix,im,itet)
c$$$              nnnx(1:3) = matmul (qbi,qcm(1:3,ix,im))+ 1d-10
c$$$              idtetfm(ix,im,itet) = indexkw(nnnx(1),nnnx(2),nnnx(3))
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c                write(6,"(' nnnx=',3i3)") nnnx(1:3)
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$              ! qbzm is given as qbzm(1:3, ib1bz(idtetfm(ix,im,itet)) )
c$$$              if(abs(sum(nnnx - matmul (qbi,qcm(1:3,ix,im))))>1d-8) then
c$$$Cstop2rx 2013.08.09 kino                stop 'getbzdata1: nnn is not integer'
c$$$                call rx( 'getbzdata1: nnn is not integer')
c$$$              endif
c$$$c               write(6,"('n1 n2 n3=',3i5)") nnnx(1:3)
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cc       volt = 0d0
c$$$c       do itet = 1, ntetf
c$$$c       do im   = 1, nmtet
c$$$c          kkvec(1:3,0:3) = qbzwm (1:3, idtetfm(0:3,im,itet) )
c$$$c          do ix = 1,3
c$$$c          kkvec(1:3,ix) = kkvec(1:3,ix) - kkvec(1:3,0)
c$$$c          enddo
c$$$cc          volt = volt + abs(det33(kvec(1:3,1:3))/6d0)
c$$$ccccccccccccccccc
c$$$c        write(6,"('itet im vol=',2i5,d13.5)") itet,im,abs(det33(kkvec(1:3,1:3))/6d0)
c$$$ccccccccccccccc
c$$$c       enddo
c$$$c       enddo
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$c-- generate index_qbzm
c$$$c            n = 0
c$$$c            allocate(index_qbzm(1,1,1)) !dummy
c$$$c 1120       continue
c$$$c            n = n+1
c$$$c            print *,' =========================== n=',n
c$$$c            deallocate(index_qbzm)
c$$$c            allocate(index_qbzm(n,n,n))
c$$$c            index_qbzm = -9999
c$$$c            do iqbz = 1,nqbzm
c$$$c               call rangedq(matmul(ginv,qbzm(:,iqbz)), qout)
c$$$c               intq =  qout*n +1
c$$$c               call checkrange(intq(1),1,n) !sanity checks
c$$$c               call checkrange(intq(2),1,n)
c$$$c               call checkrange(intq(3),1,n)
c$$$c               if(index_qbzm(intq(1),intq(2),intq(3))/=-9999) then
c$$$c                 print *,' Go to test with new n'
c$$$c                 goto 1120
c$$$c               endif
c$$$c               index_qbzm(intq(1),intq(2),intq(3)) = iqbz
c$$$c            enddo
c$$$c            n_index_qbzm = n
c$$$c            do iqbz = 1,nqbzm
c$$$c              call rangedq(matmul(ginv,qbzm(:,iqbz)), qx)
c$$$c              intq  =  qx*n + 1
c$$$c              iqbzx = index_qbzm(intq(1),intq(2),intq(3))
c$$$c              if(verbose()>50) write(6,"(' iqbz qbz intq=',i5,3f8.3,3i5)")iqbzx,qbzm(:,iqbz),intq(1),intq(2),intq(3)
c$$$c            enddo
c$$$
c$$$c--- write Qmtet
c$$$c            ifiqmtet=501
c$$$c            open(ifiqmtet, file='Qmtet')
c$$$c            write(ifiqmtet,"(i10)") nqbzm
c$$$c            do iq=1,nqbzm
c$$$c              write(ifiqmtet,"(3d24.16)") qbzm(:,iq)
c$$$c            enddo
c$$$c            close(ifiqmtet)
c$$$c--- write mtet
c$$$        ifmtet=501
c$$$        open (ifmtet, file='mtet',form='unformatted')
c$$$        write(ifmtet) nmtet,nqbzwm,nqbzm,ntetfm !,n_index_qbzm
c$$$        write(ifmtet) idtetfm,ib1bzm,qbzm,qbzwm,wtet  !,index_qbzm
c$$$        close(ifmtet)
c$$$      endif
c$$$      print *,' end of getbzdata1'
      end subroutine getbzdata1

      subroutine nkstar (qibz,qbz,grp,ginv,
     d nqibz,nqbz,ngrp,gammacellctrl,
     o nstar,irotk)
c 92.02.22
c generates the no. stars of k
c i.e. the no. times k appears in the FBZ
c qibz  = k { IBZ
c qbz   = k { FBZ
c grp   = rotation matrices
c nqibz = no. k { IBZ
c nqbz  = no. k { FBZ
c ngrp  = no. rotation matrices
c nstar(k) = no. times k appears in the FBZ
c irotk(k{IBZ,R) = index to k{FBZ
      implicit none 
      integer:: nqibz,nqbz,ngrp,ir,k,kp,nsum,ivsum
      real(8):: qibz(3,nqibz),qbz(3,nqbz),grp(3,3,ngrp),ginv(3,3)
      integer:: nstar(nqibz),irotk(nqibz,ngrp)
      integer:: verbose,kout,iccc,gammacellctrl
      real(8):: diff(3),diff2(3),tolq=1d-8,gg(3,3)
      if(verbose()>104) then
        print *,' nkstar:'
        do kp = 1,nqbz
          write(6,"(' === kp=',i8,' qbz=',3f8.3)")kp,qbz(:,kp)
        enddo
      endif
      irotk=0
      nstar=0
      do kp = 1,nqbz
        do k  = 1,nqibz
          do ir = 1,ngrp
            if(verbose()>104) print *,' grp=',ir !;      print *, grp(:,ir)
            diff = matmul(grp(:,:,ir),qibz(:,k)) - qbz(:,kp)
            if(gammacellctrl/=2) then
              call rangedq(matmul(ginv,diff), diff2)
            else
              diff2=diff
            endif
            if(verbose()>104) write(6,"(' matmul(ginv,diff)=',3f8.3,' ',3f8.3)") diff, matmul(ginv,diff)
            if(sum(abs(diff2))< tolq) then
              irotk(k,ir)= kp
              kout=k
              nstar(k)   = nstar(k) + 1
              goto 1022
            endif
          enddo
        enddo
        call rx( 'nkstar: can not find irotk')
 1022   continue
c        write(6,"('   kp=',i8,' qbz=',3f18.14, 'ibz qibz',i8,3f18.14)")kp,qbz(:,kp),kout,qibz(:,kout)
      enddo
c      do k  = 1,nqibz
c          write(6,"(' k nstar=',3i6)") k,nstar(k)
c      enddo
c      print *,'sum nstart=',sum(nstar)
      iccc=0
      do k=1,nqibz
      do ir=1,ngrp
         if(irotk(k,ir)/=0) then
            iccc=iccc+1
         endif
      enddo
      enddo
      print *,'cccccc icc sum=',iccc,sum(irotk)
ccccccccccc

!! check that the sum of stars equal to the no. k{FBZ
      nsum = ivsum (nstar,nqibz)
      if (nsum .ne. nqbz) then
        print *,' nums nqbz=',nsum,nqbz
        do k  = 1,nqibz
          do ir = 1,ngrp
            if(irotk(k,ir)/=0) write(6,"(' k ir irotk=',3i6)") k,ir, irotk(k,ir)
          enddo
        enddo
        call rx( 'nkstar: wrong no. stars')
      endif
c$$$c write k { IBZ and no. stars to file KPNT
c$$$      ifkpnt     = ifile('KPNT')
c$$$      if (ifkpnt .gt. 0) then
c$$$        write (ifkpnt,*) 'irreducible k-points and no. stars'
c$$$        write (ifkpnt,*) 'k, k-vector, nstar '
c$$$        do       k = 1,nqibz
c$$$          write(ifkpnt,"(1x,i5,3f8.5,i3)")k,qibz(1,k),qibz(2,k),qibz(3,k),nstar(k)
c$$$        end do
c$$$      endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      stop 'test end'
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine nkstar
      end module

c--------------------------------------------------------------------------
      subroutine qwider(icase,qb,n1,n2,n3,n1w,n2w,n3w, ib1bz,qbzw,qbz,indexkw)
!! == Wider q point mesh. ==
      implicit none
      integer(4):: i1,i2,i3, kount,icase,n1w,n2w,n3w
      integer(4):: n1,n2,n3,ib1bz(*),indexkw(0:n1w,0:n2w,0:n3w)
      integer(4):: indexk(0:n1-1,0:n2-1,0:n3-1)
      real(8):: qb(3,3),qbzw(1:3,*),qbz(1:3,*),hf
      hf=0d0
c      if(icase==2) hf=0.5d0
      kount      = 0
      do      i1 = 1,n1
        do      i2 = 1,n2
          do      i3 = 1,n3
            kount    = kount + 1
            indexk(i1-1,i2-1,i3-1) = kount
            qbz(1:3,kount) = qb(1:3,1)*(i1-1+hf) +qb(1:3,2)*(i2-1+hf) +qb(1:3,3)*(i3-1+hf)
          end do
        end do
      end do
c
      kount      = 0
      do      i1 = 1,n1w+1
        do      i2 = 1,n2w+1
          do      i3 = 1,n3w+1
            kount    = kount + 1
            indexkw(i1-1,i2-1,i3-1) = kount
            qbzw(1:3,kount) =
     &    qb(1:3,1)*(i1-1+hf) + qb(1:3,2)*(i2-1+hf) + qb(1:3,3)*(i3-1+hf)
            ib1bz(kount) = indexk(mod(i1-1,n1), mod(i2-1,n2), mod(i3-1,n3))
          end do
        end do
      end do
      end
!!--------------------
      subroutine tetdevide(qc, qcm, wtet, mt1,mt2,mt3)   !qc ---> qcm
      integer(4):: mt1,mt2,mt3
      real(8):: qc(3,0:3), qcm(3,0:3,mt1*mt2*mt3),qq(3,0:9),wt(0:3,0:9)
     &         ,wtet(0:3,mt1*mt2*mt3)
      integer(4):: iq(0:3,8),itet,ic
!! == four tetrahedrons at ends ==
      iq(:,1) = (/0,7,8,9/)
      iq(:,2) = (/1,5,6,7/)
      iq(:,3) = (/2,4,6,9/)
      iq(:,4) = (/3,4,5,8/)
c octahedron into four tetrahedron.
      iq(:,5) = (/5,9,6,7/)
      iq(:,6) = (/5,9,7,8/)
      iq(:,7) = (/5,9,8,4/)
      iq(:,8) = (/5,9,4,6/)
      if(mt1==2.and.mt2==2.and.mt3==2) then
        qq(:,0:3)  =  qc(:,0:3)
        qq(:,3+1)  = (qc(:,2) +qc(:,3))/2d0
        qq(:,3+2)  = (qc(:,3) +qc(:,1))/2d0
        qq(:,3+3)  = (qc(:,1) +qc(:,2))/2d0
        qq(:,6+1)  = (qc(:,1) +qc(:,0))/2d0
        qq(:,6+2)  = (qc(:,3) +qc(:,0))/2d0
        qq(:,6+3)  = (qc(:,2) +qc(:,0))/2d0
c
        wt(:,0)  =  (/1d0,0d0,0d0,0d0/)
        wt(:,1)  =  (/0d0,1d0,0d0,0d0/)
        wt(:,2)  =  (/0d0,0d0,1d0,0d0/)
        wt(:,3)  =  (/0d0,0d0,0d0,1d0/)
        wt(:,3+1)  = (wt(:,2) +wt(:,3))/2d0
        wt(:,3+2)  = (wt(:,3) +wt(:,1))/2d0
        wt(:,3+3)  = (wt(:,1) +wt(:,2))/2d0
        wt(:,6+1)  = (wt(:,1) +wt(:,0))/2d0
        wt(:,6+2)  = (wt(:,3) +wt(:,0))/2d0
        wt(:,6+3)  = (wt(:,2) +wt(:,0))/2d0
        do itet=1,8
          do ic=0,3
            qcm (:,ic,itet) = qq(:,iq(ic,itet))
          enddo
          wtet(0,itet) =  sum( wt(0,iq(:,itet)) )/4d0
          wtet(1,itet) =  sum( wt(1,iq(:,itet)) )/4d0
          wtet(2,itet) =  sum( wt(2,iq(:,itet)) )/4d0
          wtet(3,itet) =  sum( wt(3,iq(:,itet)) )/4d0
c         write(6,"(' itet wtet=',i5,5f8.3)")itet,wtet(:,itet),sum(wtet(:,itet))
          if(abs(sum(wtet(:,itet))-1d0)>tolq) call rx( 'tetdevide: sumwtet/=1')
        enddo
      else
        call rx( ' tetdvide: only 2 2 2 has already implimented.')
      endif
      end
