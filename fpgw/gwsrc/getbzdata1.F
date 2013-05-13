C  ... Get BZDATA case1
      module m_get_bzdata1
c In addition to these variables, this also write mtet file (search ifmtet) for mtet mode.
c------------------------------------------------------------
      implicit none  ! all are outputs
      integer(4):: nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm,nqibz_r
      real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:),wibz(:)
     &   ,qbzw(:,:) ,qibz_r(:,:)
      real(8):: qbasmc(3,3),dq_bzcase2(3)
      integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
     &   ,irk(:,:),nstar(:),nstbz(:)
      real(8),allocatable:: qbzm(:,:),qbzwm(:,:)

c------------------------------------------------------------
      contains
      subroutine getbzdata1(icase, plat,qlat,ginv,n1q,n2q,n3q
     &    ,symops,ngrp,tetrai,tetraf,mtet) !all are inputs
      implicit none
      logical :: tetrai,tetraf,multet
      integer(4):: ngrp,n1q,n2q,n3q, mtet(3),n1qm,n2qm,n3qm,nnn(3),
     &    itet,ix,im,ifiqmtet,iq,icase
      real(8):: plat(3,3),qlat(3,3),ginv(3,3),
     &    symops(3,3,ngrp),qc(3,0:3),qb(3,3),qbi(3,3)
      integer(4),allocatable:: IPQ(:),iw1(:),indexkw(:,:,:)
      real(8),allocatable:: qcm(:,:,:),wtet(:,:,:)
c
c---for multipuled tetrahedron.
      integer(4):: nmtet,nqbzwm,ntetfm,ifmtet
      integer(4),allocatable::
     &     idtetfm(:,:,:),ib1bzm(:)
      integer(4):: n,iqbz,intq(3),iqx,iqy,iqz,iqbzx,verbose !,n_index_qbzm
      real(8)::qout(3),qx(3),kvec(3,3),det33,  kkvec(3,0:3)
c          integer(4),allocatable:: index_qbzm(:,:,:)
      real(8):: hf
      integer:: iqi,igrp,iccc
c         logical:: testtr !invert qibz

      hf=0d0
      if(icase==2) hf=0.5d0
c-------------------------------------------------------
      write(6,*)'getbzdata1: in m_getbzdata1'
      write(6,"(' mkqg: ibzcase()=',i3,'-------------------' )") icase
      if(icase>2) stop 'mkqg: wrong bzcase'
      print *,' mkqg: CASE1:', n1q,n2q,n3q

      nqbz = n1q*n2q*n3q
      allocate(qbz(3,nqbz),wbz(nqbz),nstbz(nqbz))
      call genqbz (icase,qlat,n1q,n2q,n3q,
     o     qbz,wbz,nstbz )

ccccccccccccccccccccccccccccccccccccccccc
c        do iq=1,nqbz
c          write(6,"(' iq qbz nstbz=',i5,3f9.4,i5)")iq,qbz(:,iq),nstbz(iq)
c        enddo
c        stop 'xxxxxxxxxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccccc

      print *,' ---goto bzmesh--- '
      call pshprt(51)
      if(icase/=1) then !get qibz_r, and nqibz_r
        allocate(qibz_r(3,nqbz),ipq(nqbz),wibz(nqbz))
        call bzmesh(1,plat,qbasmc,n1q,n2q,n3q,symops,ngrp,ipq,qibz_r
     &      ,wibz,nqibz_r,nqbz) !Make q-points in IBZ.
        deallocate(ipq,wibz)
      else
        nqibz_r=0
      endif

      allocate(qibz(3,nqbz),ipq(nqbz),wibz(nqbz))
      call bzmesh(icase,plat,qbasmc,n1q,n2q,n3q,symops,ngrp,ipq,qibz
     &    ,wibz,nqibz,nqbz) !Make q-points in IBZ.

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c          if(testtr()) then
c           qibz_r=-qibz_r
c           qibz = -qibz
c          endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      print *,'  ---goto nkstar---'
      allocate(nstar(nqibz), irk(nqibz,ngrp),iw1(nqbz))
c          print *,'  ---ginv='
c          print *, ginv
      call nkstar  (qibz,qbz,symops,ginv, !ginv=plat
     w              iw1,
     d              nqibz,nqbz,ngrp,
     o              nstar,irk )

cccccccccccccccccccccccccccccccccccccccc
      iccc=0
      do iqi=1,nqibz
c         print *
      do igrp=1,ngrp
         if(irk(iqi,igrp)/=0) then
c            write(6,*) 'iq igrp irk=',iqi,igrp,irk(iqi,igrp)
            iccc=iccc+1
         endif
      enddo
      enddo
      write(6,*)'getbzdata1: total count=',iccc
cccccccccccccccccccccccccccccccccccccccccc

      ntetf=-1
      nteti=-1
      nqbzw=-1
      if(tetraf) then
        nqbzw = (n1q+1)*(n2q+1)*(n3q+1)
        ntetf = 6*nqbz
        allocate( idtetf(0:3,ntetf), ib1bz(nqbzw), qbzw(3,nqbzw) )
        print *,' goto tetfbzf---'
        call TETFBZF(icase,qlat, n1q,n2q,n3q,qbz,nqbz,
     o      idtetf,qbzw,ib1bz)
      endif
      if(tetrai) then
        allocate( idteti(0:4,6*nqbz))
        print *,' goto tetirr---' !this should work also for icase==2
        call tetirr(qbasmc,n1q,n2q,n3q,ipq,nqibz,
     o              nteti,idteti)
      endif
      deallocate(ipq,iw1)

c--- for bzcase==2
      dq_bzcase2 = -matmul(qlat(1:3,1:3),(/.5d0/n1q,.5d0/n2q,.5d0/n3q/) )


c---for multipled tetrahedron-----------------
      if(sum(abs(mtet))/=3) then
        print *, 'multitet mode: mtet=',mtet
        n1qm = mtet(1)*n1q
        n2qm = mtet(2)*n2q
        n3qm = mtet(3)*n3q
        if(icase==2) then
          n1qm=n1qm+1
          n2qm=n2qm+1
          n3qm=n3qm+1
        endif
        nmtet   = mtet(1)*mtet(2)*mtet(3)
        nqbzm   = nmtet * nqbz
        nqbzwm  = (n1qm+1)* (n2qm+1)* (n3qm+1)
        ntetfm  = ntetf * nmtet
        allocate( 
     &       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
        ! Index for tetrahedron;    qbzmw(idtetfm) gives extended q vector.
     &       ib1bzm(nqbzwm), qbzm(3,nqbzm) )
        ! qbzm(1:3,ib1bz(iq)) gives q vector within the 1st bz.
! The datas idetetfm, qbzmw, ib1bzm, nmete, ntetf, nqbzm,
!     eigen(nband,nqbzm) are required for the multiply-devided tetrahedron method.
!
        allocate( qcm(1:3,0:3, nmtet), wtet(0:3, nmtet, ntetf), 
     &          indexkw(0:n1qm,0:n2qm,0:n3qm) )
        QB(:,1) = qlat(:,1)/(n1q*mtet(1))
        QB(:,2) = qlat(:,2)/(n2q*mtet(2))
        QB(:,3) = qlat(:,3)/(n3q*mtet(3))
        call dinv33x(qb,qbi)
        qbi = transpose(qbi)
        call qwider(icase,qb, mtet(1)*n1q,mtet(2)*n2q,mtet(3)*n3q,
     &                       n1qm,n2qm,n3qm, ib1bzm, qbzwm, qbzm, indexkw )
        do itet=1, ntetf  !--- idtetfm(0:3, nmtet, ntetf), nmtet divition.
          do ix=0,3
            qc(:,ix)= qbzw(:, idtetf(ix,itet))
          enddo
ccccccccccccccccccccccccccccccccccccccccccccccc
c        do ix = 1,3
c        kvec(1:3,ix) = qc(1:3,ix) - qc(1:3,0)
c        enddo
c        write(6,"('itet vol=',i5,d13.5)") itet,abs(det33(kvec(1:3,1:3))/6d0)
cccccccccccccccccccccccccccccccccccccccccccccccc
          call tetdevide(qc,qcm, wtet(:,:,itet),mtet(1),mtet(2),mtet(3))   !qc ---> qcm
          do im=1,nmtet ! Index for micro-devided tetrahedron
cccccccccccccccccccccccccccccccccccccccccccccccc
c        do ix = 1,3
c        kvec(1:3,ix) = qcm(1:3,ix,im) - qcm(1:3,0,im)
c        enddo
c        write(6,"('itet im vol=',2i5,d13.5)") itet,im,abs(det33(kvec(1:3,1:3))/6d0)
cccccccccccccccccccccccccccccccccccccccccccccccc
            do ix=0,3     ! its four corners
              !qcm ---> idtetfm(ix,im,itet)
              nnn(1:3) = matmul (qbi,qcm(1:3,ix,im))+ 1d-10
              idtetfm(ix,im,itet) = indexkw(nnn(1),nnn(2),nnn(3))
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                write(6,"(' nnn=',3i3)") nnn(1:3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
              ! qbzm is given as qbzm(1:3, ib1bz(idtetfm(ix,im,itet)) )
              if(abs(sum(nnn - matmul (qbi,qcm(1:3,ix,im))))>1d-8) then
                stop 'getbzdata1: nnn is not integer'
              endif
c               write(6,"('n1 n2 n3=',3i5)") nnn(1:3)
            enddo
          enddo
        enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc       volt = 0d0
c       do itet = 1, ntetf
c       do im   = 1, nmtet
c          kkvec(1:3,0:3) = qbzwm (1:3, idtetfm(0:3,im,itet) )
c          do ix = 1,3
c          kkvec(1:3,ix) = kkvec(1:3,ix) - kkvec(1:3,0)
c          enddo
cc          volt = volt + abs(det33(kvec(1:3,1:3))/6d0)
ccccccccccccccccc
c        write(6,"('itet im vol=',2i5,d13.5)") itet,im,abs(det33(kkvec(1:3,1:3))/6d0)
ccccccccccccccc
c       enddo
c       enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





c-- generate index_qbzm
c            n = 0
c            allocate(index_qbzm(1,1,1)) !dummy
c 1120       continue
c            n = n+1
c            print *,' =========================== n=',n
c            deallocate(index_qbzm)
c            allocate(index_qbzm(n,n,n))
c            index_qbzm = -9999
c            do iqbz = 1,nqbzm
c               call rangedq(matmul(ginv,qbzm(:,iqbz)), qout)
c               intq =  qout*n +1
c               call checkrange(intq(1),1,n) !sanity checks
c               call checkrange(intq(2),1,n)
c               call checkrange(intq(3),1,n)
c               if(index_qbzm(intq(1),intq(2),intq(3))/=-9999) then
c                 print *,' Go to test with new n'
c                 goto 1120
c               endif
c               index_qbzm(intq(1),intq(2),intq(3)) = iqbz
c            enddo
c            n_index_qbzm = n
c            do iqbz = 1,nqbzm
c              call rangedq(matmul(ginv,qbzm(:,iqbz)), qx)
c              intq  =  qx*n + 1
c              iqbzx = index_qbzm(intq(1),intq(2),intq(3))
c              if(verbose()>50) write(6,"(' iqbz qbz intq=',i5,3f8.3,3i5)")iqbzx,qbzm(:,iqbz),intq(1),intq(2),intq(3)
c            enddo

c--- write Qmtet
c            ifiqmtet=501
c            open(ifiqmtet, file='Qmtet')
c            write(ifiqmtet,"(i10)") nqbzm
c            do iq=1,nqbzm
c              write(ifiqmtet,"(3d24.16)") qbzm(:,iq)
c            enddo
c            close(ifiqmtet)
c--- write mtet
        ifmtet=501
        open (ifmtet, file='mtet',form='unformatted')
        write(ifmtet) nmtet,nqbzwm,nqbzm,ntetfm !,n_index_qbzm
        write(ifmtet) idtetfm,ib1bzm,qbzm,qbzwm,wtet  !,index_qbzm
        close(ifmtet)
      endif
      end subroutine
      end module

c--------------------------------------------------------------------------
      subroutine qwider(icase,qb,n1,n2,n3,n1w,n2w,n3w, ib1bz,qbzw,qbz,indexkw)
c wider q point mesh.
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

      subroutine tetdevide(qc, qcm, wtet, mt1,mt2,mt3)   !qc ---> qcm
      integer(4):: mt1,mt2,mt3
      real(8):: qc(3,0:3), qcm(3,0:3,mt1*mt2*mt3),qq(3,0:9),wt(0:3,0:9)
     &         ,wtet(0:3,mt1*mt2*mt3)
      integer(4):: iq(0:3,8),itet,ic
c---
c four tetrahedrons at ends
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
          if(abs(sum(wtet(:,itet))-1d0)>1d-10) stop 'tetdevide: sumwtet/=1'
        enddo
      else
        stop ' tetdvide: only 2 2 2 has already implimented.'
      endif
      end
