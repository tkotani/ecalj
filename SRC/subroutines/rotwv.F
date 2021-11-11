!! See fpgw/gwsrc/readeigen.F. We will share the same code in fpgw and lm7k for maitenance.
!! Now rotwv rotwvip are moved to fpgw/gwsrc/readeigen.F
!! rotsig is essentially similar with rotwv.

      subroutine interparam(q,wgt)
C- Calculate interpolarion weights when values at corners of a cube is given.
C  The cube is devided into 6 piramids, whose center are the center of a cube.
Ci q(3) is within a cube, that is, q(i) \in (0,1) for i=1,3.
C     At first, we judge which piramids include the q vector.
C  Then we further devide the piramids into four tetrahedron.
C  We judge which tetrehedron includes the q vector.
C  The interpolation weight is calculated from the linear-interpolarion formular from its four corners.
C- cube
      implicit none
      integer:: ip,ix1,ix2,ix3
      real(8):: q(3),wgt(0:1,0:1,0:1),xx(3),x2d(2),wc
      real(8):: wgtx(0:1,0:1),maxx
      integer:: iz(1),ix
      xx = 2d0*q-(/1d0,1d0,1d0/)
c weight for a center.
      ix = maxloc(abs(xx),dim=1)
      !ix=iz(1)
      maxx = abs(xx(ix))
      if(xx(ix)>=0) ip=1
      if(xx(ix)<0)  ip=0

      print *,' xx=',xx
      print *,' ix=',ix
      print *,' maxx=',maxx

c
      wgt=1d0/8d0*(1-maxx) !center weight
      if(ix==1) then
        x2d=(/xx(2),xx(3)/)
        call prodx2d(x2d,wgtx,wc)
        wgt(ip,:,:) = wgt(ip,:,:) + wgtx + 1/4d0 * (maxx + wc)
        ! face-center weight --> 1/4d0*(maxx-1d0) + 1/4d0*wc
      elseif(ix==2)  then
        x2d=(/xx(1),xx(3)/)
        call prodx2d(x2d,wgtx,wc)
        wgt(:,ip,:) = wgt(:,ip,:) + wgtx + 1/4d0 * (maxx + wc)
      else
        x2d=(/xx(1),xx(2)/)
        call prodx2d(x2d,wgtx,wc)
        wgt(:,:,ip) = wgt(:,:,ip) + wgtx + 1/4d0 * (maxx + wc)
      endif
ccccccccccccccccccccccccccccccccccccccc
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3f8.3,a,3i3,f8.3)") '  q=', q,'  wgt=',ix1,ix2,ix3,wgt(ix1,ix2,ix3)
          enddo
        enddo
      enddo
cccccccccccccccccccccccccccccccccccccccc
      end
c
      subroutine prodx2d(x,wgt,wc)
      implicit none
      real(8):: x(2),wgt(0:1,0:1),wc,wpp,wpm
      real(8):: vpp(2),vpm(2)
      vpp=(/.5d0,.5d0/)
      vpm=(/.5d0,-.5d0/)
      wpp = sum(vpp*x)
      wpm = sum(vpm*x)
      print *,' wpp wpm=',wpp,wpm
      wgt=0d0
      if(wpp>=0) then
        wgt(1,1)=  wpp
        wc= -wpp
      else
        wgt(0,0)= -wpp
        wc=  wpp
      endif
      if(wpm>=0) then
        wgt(1,0)=  wpm
        wc= wc-wpm
      else
        wgt(0,1)= -wpm
        wc= wc+wpm
      endif
      end

!!== Obtain q parameter for given qin(3) ==
      subroutine qparam(qin,qsmesh,nk1,nk2,nk3,plat,  ikc,qpara,igshift)
!!  This find eight corners of a cube (each corner is one of the qsmesh). qin is within the cube.
Ci qsmesh(3,ik1,ik2,ik3); we assume it orders well; this routine should work even for regular mesh or off-gamma mesh.
Ci
Co ikc(1:3) ,qpara, igshift
Co
Cr not sophsticated well... At first, we need to clarify the definition and generality of qsmesh itself.
Cr
      implicit none
      integer::nk1,nk2,nk3, nkk(3),ik1,ik2,ik3,i,j,k,ikcc(3),
     &nx(0:1),ny(0:1),nz(0:1), mx(0:1), my(0:1), mz(0:1),ix1,ix2,ix3,
     &igshift(3,0:1,0:1,0:1), ikc(3,0:1,0:1,0:1)
      real(8):: qsmesh(3,nk1,nk2,nk3), qpara(3),qqq(3),qin(3),plat(3,3)
      real(8):: eps=1d-5, qmesh(3), qmeshm(3,0:1,0:1,0:1),qmin(3),qmax(3),qx(3),qqqm(3)
c      real(8),allocatable:: qmesh(:,:,:,:)

c Get projected q.-> ik1, ik2 ik3 for given q.
      qqq = matmul(qin,plat)
      nkk= (/nk1,nk2,nk3/)
      ikcc = qqq*nkk + 100*nkk !(1) ---> Takao suppose adding 100*nkk is enough to make ikcc positive.
      ikcc = ikcc    - 100*nkk !(2) Note that this conversion to integer, "(1) and (2)",
      !    do not equilvalent with ikc=qqq*nkk when qqq*nkk is negative.
      !
      qpara = nkk*qqq - ikcc
      ikcc= mod(ikcc+100*nkk,nkk)
      qqqm= (qpara + ikcc)/dble(nkk) !positive

      print *,'nkkk=',nkk
      write(6,"(a,3d13.5,3d13.5,3i5)")'  qin   =',qin
      write(6,"(a,3d13.5,3d13.5,3i5)")'  qqqm   =',qqqm
      write(6,"(a,3d13.5,3d13.5,3i5)")'  qpara =',qpara
      write(6,"(a,3i5)")              '  ikcc  =',ikcc

c qmesh is based on on qlat coordinate. Is this too time-consuming (maybe negligible)?
c         But maybe too general or
c      allocate( qmesh(3,nk1,nk2,nk3) )
      do ik1 = nk1,1,-1
        qmesh = matmul(qsmesh(:,ik1,1,1),plat)
        !print *,' aaa=',qmesh(1),qqqm(1)
        if( qmesh(1)-eps < qqqm(1) ) then
          nx(0)= ik1
          nx(1)= mod(ik1+1-1,nk1)+1
          mx(0) = 0
          mx(1) = (ik1+1- nx(1))/nk1
          go to 1011
        endif
      enddo
      stop ' no exit 1xx'
 1011 continue
      do ik2 = nk2,1,-1
        qmesh = matmul(qsmesh(:,1,ik2,1),plat)
        if( qmesh(2)- eps < qqqm(2) ) then
          ny(0)= ik2
          ny(1)= mod(ik2+1-1,nk2)+1
          my(0)= 0
          my(1)= (ik2+1- ny(1))/nk2
          goto 1012
        endif
      enddo
      stop 'no exit 2xx'
 1012 continue
      do ik3 = nk3,1,-1
        qmesh = matmul(qsmesh(:,1,1,ik3),plat)
        if( qmesh(3)- eps < qqqm(3)  ) then
          nz(0)= ik3
          nz(1)= mod(ik3+1-1,nk3)+1
          mz(0)= 0
          mz(1)= (ik3+1- nz(1))/nk3
          goto 1013
        endif
      enddo
      stop 'no exit 3xx'
 1013 continue

      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            ikc(1:3,ix1,ix2,ix3)     = (/nx(ix1),ny(ix2),nz(ix3)/)
            igshift(1:3,ix1,ix2,ix3) = (/mx(ix1),my(ix2),mz(ix3)/)
          enddo
        enddo
      enddo

ccccccccccccccccccccccccccccccccc
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3i3,3x,4i3)") 'ix1 ix2 ix3 ikc=',ix1,ix2,ix3,ikc(1:3,ix1,ix2,ix3)
          enddo
        enddo
      enddo
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3i3,3x,4i3)") '        igshift=',ix1,ix2,ix3,igshift(1:3,ix1,ix2,ix3)
          enddo
        enddo
      enddo
      do ix1=0,1
        do ix2=0,1
          do ix3=0,1
            write(6,"(a,3i3,3x,4i3)") '    ikc+nkk*igshift=',ix1,ix2,ix3,ikc(1:3,ix1,ix2,ix3) + nkk*igshift(1:3,ix1,ix2,ix3)
          enddo
        enddo
      enddo
ccccccccccccccccccccccccccccccccc

c check
      qmin  = matmul( qsmesh(:,nx(0),ny(0),nz(0)), plat )
      qx  = matmul( qsmesh(:,nx(1),ny(0),nz(0)), plat ) + igshift(1:3,1,0,0) !/dble(nk1)
      qmax(1)=qx(1)
      qx = matmul( qsmesh(:,nx(0),ny(1),nz(0)), plat )  + igshift(1:3,0,1,0) !/dble(nk2)
      qmax(2)=qx(2)
      qx = matmul( qsmesh(:,nx(0),ny(0),nz(1)), plat )  + igshift(1:3,0,0,1) !/dble(nk3)
      qmax(3)=qx(3)
c
      call regioncheck(qqqm,qmin,qmax)
      end
c
      subroutine regioncheck(qqqm,qmin,qmax)
      integer:: ix
      real(8):: qmin(3),qmax(3),qqqm(3),eps=1d-8
      do ix=1,3
        if( qmin(ix)-eps <= qqqm(ix) .and. qqqm(ix) <= qmax(ix)+eps ) then
        else
          goto 1010
        endif
      enddo
      return
 1010 continue
      write(6,*)' error regioncheck: ix  =',ix
      write(6,*)' error regioncheck: qqqm=',qqqm
      write(6,*)' error regioncheck: qmin=',qmin
      write(6,*)' error regioncheck: qmax=',qmax
      stop 'error in region check!'
      end


      subroutine rtoint(a,k,n)
      real(8):: a(n),eps=1d-12
      integer:: k(n),ix,n,i
      do i=1,n
        if(abs(a(i))<eps) then
          k(i)=0
        elseif(a(i)>0) then
          k(i) = a(i) + eps
        else
          k(i) = a(i) - eps
        endif
      enddo
      end
c
      subroutine rotsig(qin,qout,ndimh,napw_in,ldim,sig,sigout,ierr,iaf) !iaf is added @jun2015takao
c      subroutine rotsig(qin,qout,ndimh,napw_in,plat,qlat,ldim,sig,sigout,ierr)
      use m_hamindex,only: symops,invgx,miat,tiat,shtvg,qlat,plat,dlmm,ngrp,norbmto,
     &  ibastab,ltab,ktab,offl,offlrev,getikt,igv2,igv2rev,napwk,
     &  ngrp_original
      implicit none
!!=== sigma rotator. sig at qin to sig at qout. ===
!! obtain sigout for qout.
!! a little confusing since qin=symops(qout), and basis rotation. Need to clean up.
!! Both of q and qtarget should be in qq table(in m_hamindex) which is given by gen_hamindex 
!! (read from QGpsi).
!! Used idea is   <base(qin)|sigma(qin)|base(qin)> = <g(base)|sigma |g(base)>.   
!! where qin=g(qout).  qtarget=qin= g(q=qout)

      integer   :: ig,ndimh,napw_in,ibaso,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr,igg,ikt2,ikt,l,ibas,ig2,k,ix
      real(8):: qin(3),qout(3)
      real(8)   :: q(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3) !plat(3,3),qlat(3,3)
      complex(8):: phase,img=(0d0,1d0),img2pi

      complex(8):: sig(ndimh,ndimh),sigout(ndimh,ndimh)
      complex(8),allocatable:: sigx(:,:)
      integer :: ldim,debugmode,iaf,ngini,ngend
      character(300)::aaa
      img2pi=2*4d0*datan(1d0)*img
      ierr=1
      platt=transpose(plat) !this is inverse of qlat

!! find symops(3,3,ig),
      if(debugmode()>0) write(6,"('rotsig: qin qout=',3f9.4,x,3f9.4)") qin,qout
      qtarget= qin
      q      = qout  ! we find q

! we try to find qtrget = symops(igx) * q          (this means qin = symops(igc) qout).
! this definitions for q and qtarget are following what is used in rotwv.
!
     
      if(iaf==1) then !iaf mechanism is adde @jun2015takao for antiferro materials such as Nd2CuO4.
        ngini = 1
        ngend = ngrp_original
      elseif(iaf==2) then
        ngini = ngrp_original + 1
        ngend = ngrp
      else
        ngini = 1
        ngend = ngrp
      endif
ccccccccccccccccccccccccccccccccc
c      write(6,*)'ngini=',ngini,ngend
ccccccccccccccccccccccccccccccccc
!
      do igx=ngini,ngend
        if(debugmode()>0) print *, 'ddd=',matmul(platt,(qtarget-matmul(symops(:,:,igx),q)))
        call rangedq(   matmul(platt,(qtarget-matmul(symops(:,:,igx),q))), qx)
c        qx=matmul(platt,(qtarget-matmul(symops(:,:,igx),q)))
c        if(sum(abs(qx-anint(qx)))<1d-8) then
        if(sum(abs(qx))<1d-6) then
          igg=igx
c          write(6,"('ttttttt: qtarget=',3f9.3,2x,3f9.3)") qtarget,q
          if(debugmode()>0) then
            print *,'ttt: q      =',q
            print *,'ttt: qtarget=',qtarget
            print *,'ttt: matmul q =',matmul(symops(:,:,igx),q)
            print *,'ttt: rotsig: OK! igg=',igg
            print *
          endif
          goto 1012
        endif
      enddo
      write(aaa,"(a,3f7.3,2x,3f7.3)")' rotsig: qtarget is not a star of q',q,qtarget
      call rx(trim(aaa))
      print *
      return
 1012 continue

c ddd=ddd2 is not necessary. There is ambiguity. takao feb2011
c      ddd = - matmul(symops(:,:,invgx(igg)),shtvg(:,igg))
c      ddd2 = shtvg(:,invgx(igg))
ctakao think this is not necessary feb2012
c$$$      if(sum(abs(ddd-ddd2))>1d-8) then
c$$$         print *,'aaa: ddd ddd2=',ddd,ddd2
c$$$         print *,'rotsig: ddd/=ddd2'
c$$$         stop 'rotsig: ddd/=ddd2'
c$$$      endif
ccccccccccccccccccccccc
c      print *,'takao test 111 ddd=ddd2'
c      ddd=ddd2
ccccccccccccccccccccccc

      allocate(sigx(ndimh,ndimh))
      sigx=0d0

c      nlmto = ndimh-napw_in
      nlmto=ldim
      if(debugmode()>0) then
        print *,' tttt: invgx =',invgx(igg),shtvg(:,igg)
c        print *,' tttt: ',ddd !,ddd2
        print *,' tttt: ntorb napwin',norbmto,ndimh,napw_in,nlmto
      endif

c mto part ---
c          write(6,"(3i3,2x,2i5,3x,a)")
c     &    ib,l,k, offl(iorb,ib)+1,offl(iorb,ib)+2*l+1,trim(spid)
      if(nlmto/=0 )then
        ibaso=-999
        do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          !Index of Hamiltonian is devided into these blocks.
          ibas = ibastab(iorb)
          if(ibas/=ibaso) phase = exp( -img2pi*sum(qtarget*tiat(:,ibas,igg)) )
          ibaso=ibas
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+ 1
          iend1 = offl(iorb)+ 2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+ 1
          iend2 = offlrev(miat(ibas,igg),l,k)+ 2*l+1
          do ix=1,ndimh 
            sigx(ix,init1:iend1)= matmul(sig(ix,init2:iend2),dlmm(-l:l,-l:l,l,igg))*phase
          enddo
        enddo
      endif
c apw part ------------
      if(napw_in/=0) then
        write(*,*) ' Probably OK-->Remove this stop to use this branch.'
     &   //' But need to confirm two apw sections in this routines.(phase factors) '
     &   //' Idea of this routine: <i|\sigma|j>_qout= <g(i)|\sigma|g(i)>_qin, where qin=g(qout)'
        stop 'abort'
        ikt  = getikt(q)    !index for q
        ikt2 = getikt(qtarget) !index for qtarget
        print *,' rotsig ikt ikt2=',ikt,ikt2
        if(napw_in /= napwk(ikt) ) then
          call rx('rotsig: napw_in /= napw(ikt)')
        endif
        do ig = 1,napw_in
          qpg = q + matmul( qlat(:,:),igv2(:,ig,ikt))      !q+G
          qpgr = matmul(symops(:,:,igg),qpg)               !rotated q+G
          call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
          print *,ig,'nnn  ikt2=',nnn,ikt2
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2)
          phase= exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
          do ix=1,ndimh
            sigx(ix,nlmto+ig) = sig(ix,nlmto+ig2) * phase
          enddo
        enddo
      endif
      if(debugmode()>0) print *,' apw part end 111: ikt ikt2=',ikt,ikt2
c mto part ------
      if(nlmto/=0) then
        ibaso=-999
        do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          !Index of Hamiltonian is devided into these blocks.
          ibas = ibastab(iorb)
          if(ibas/=ibaso) phase = exp( img2pi*sum(qtarget*tiat(:,ibas,igg)) )
          ibaso=ibas
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+ 1
          iend1 = offl(iorb)+ 2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+ 1
          iend2 = offlrev(miat(ibas,igg),l,k)+ 2*l+1
          do ix=1,ndimh
          sigout(init1:iend1,ix)= phase * matmul(transpose(dlmm(-l:l,-l:l,l,igg)),sigx(init2:iend2,ix))
          enddo
        enddo
      endif
      if(debugmode()>0) print *,' end of 2nd mto part q=',q
c apw part ------------
      if(napw_in/=0) then
        ikt  = getikt(q)    !index for q
        ikt2 = getikt(qtarget) !index for qtarget
        if(debugmode()>0) print *,' rotsig 111 ikt ikt2=',ikt,ikt2
        if(napw_in /= napwk(ikt) ) then
          call rx('rotsig: napw_in /= napw(ikt)')
        endif
        do ig = 1,napw_in
          qpg = q + matmul( qlat(:,:),igv2(:,ig,ikt))      !q+G
          qpgr = matmul(symops(:,:,igg),qpg)               !rotated q+G
          call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
          if(debugmode()>0) print *,ig,'nnn  ikt2=',nnn,ikt2
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2)
          phase=exp(img2pi*sum(qpgr*shtvg(:,igg))) 
          do ix=1,ndimh
            sigout(nlmto+ig,ix) =   sigx(nlmto+ig2,ix) * phase
          enddo
        enddo
        if(debugmode()>0) print *,' apw part end 222: ikt ikt2=',ikt,ikt2
      endif
      ierr=0
      if(debugmode()>0) print *,' goto deallcate sigx'
      deallocate(sigx)
      if(debugmode()>0) print *,' end of rotsig'
      end subroutine rotsig


!!------------------------------------------------------
c$$$      subroutine conv2int(a,ii,n)
c$$$      implicit none
c$$$      integer::ii(n),i,n
c$$$      real(8):: a(n),eps=1d-8,aa
c$$$      do i=1,n
c$$$        ii(i)= nint(a(i))
c$$$        if(abs(ii(i)-a(i))>eps)then
c$$$           write(6,*)' ii a eps=',ii(i),a(i),eps
c$$$           call rx( 'conv2int: a is not close to integer')
c$$$        endif   
c$$$      enddo
c$$$      end
      subroutine conv2int(a,ii,n)
      implicit none
      real(8), intent(in) :: a(n)
      integer, intent(in) :: n
      integer, intent(out) :: ii(n)
      integer::i
      real(8):: eps=1d-5,aa !eps=1d-8 until aug2012
      do i=1,n
        ii(i)=nint(a(i))
        if(abs(ii(i)-a(i))>eps)then
           write(6,*)' ii a eps=',ii(i),a(i),eps
           call rx( 'conv2int: a is not close to integer')
        endif   
      enddo
      end subroutine conv2int

      subroutine rangedq(qin, qout) !renewed jun2015
      integer :: ix
      real(8):: qin(3),qout(3),tol=1d-12 
      qout= qin-nint(qin)
      do ix=1,3
        if(qout(ix)>0.5d0-tol) qout(ix)=-0.5d0 !this is needed to distinguish 0.5d0 and -0.5d0.
      enddo 
      end

