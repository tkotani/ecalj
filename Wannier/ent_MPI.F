      subroutine zgesvdnn(ngb,zzz, SS,UU,VT)
C--- SVD zzz= matmul(UU,matmul(SS,VT)) ------------
c$$$C--- SVD of chipm test !SVD procedure is not required to calculate <m|chi|m>
c$$$      lwork=4*ngb
c$$$      allocate(work(LWORK),rwork(5*ngb),zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb),VTT(ngb,ngb),ss0(ngb))
c$$$      zzz= matmul(transpose(conjg(ovlpi)), matmul(zxq(1:ngb,1:ngb,iw), ovlpi) )
c$$$      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS0,UU,ngb,VT,ngb,work,lwork,rwork,info) 
c$$$      write(6,*)' -------SVD: Oinv*chipm*Oinv ----------------'
c$$$      aaax = 0d0
c$$$      do i=1,ngb
c$$$        addx= sum(svec(1:nbloch)*uu(1:nbloch,i)) *ss0(i)* sum(VT(i,1:nbloch)*svec(1:nbloch)) 
c$$$        write(6,"(' SVD OcO: eig_k <m|chi|m>_k=',i4,2x, d13.5,2x,2d14.6)")i,SS0(i),addx
c$$$        if(i<25) aaax= aaax+ addx
c$$$      enddo 
c$$$      aaax= mmnorm**2/aaax
c$$$      deallocate(work,rwork,zzz,uu,vt,vtt)
c$$$      deallocate(ovlpi)
      implicit none
      integer(4)::lwork,info,ngb,i
      complex(8):: zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb)
      real(8):: ss(ngb)
      real(8),allocatable:: rwork(:)
      complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
      lwork=4*ngb
      allocate(zw0bk(ngb,ngb))
      allocate(work(LWORK),rwork(5*ngb)) !,VTT(ngb,ngb))
      zw0bk = zzz
!      print *,' zgesvdnn: singular value decomp '
      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS,UU,ngb,VT,ngb,work,lwork,rwork,info) 
!      do i=1,ngb
!         write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
!         vtt(i,:)=ss(i)*vt(i,:)
!      enddo
!      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
!     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
!      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk))) 
!     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
!      deallocate(vtt)
      end

C===================================================================
      subroutine get_nwf_MPI(nwf)
      use rsmpi !RS
      implicit none
      integer(4):: nwf,nqbz,iko_ix,iko_fx,ifi
      integer(4):: iopen,iclose

      if (Is_IO_Root_RSMPI()) then
         ifi  = iopen('MLWU',0,0,0)
         read(ifi)nqbz,nwf,iko_ix,iko_fx
         ifi = iclose('MLWU')
      endif
      call MPI_Bcast(nwf,1,MPI_INTEGER,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(nwf)",ierror_rsmpi)
      end
C===================================================================
      subroutine choose_wanband_MPI(iwbnd,nwf,nqbze,nspin)
      use rsmpi !RS
      implicit none
      integer(4):: nwf,nqbze,nspin,iwbnd(nwf,nqbze,nspin)
      integer(4):: iopen,iclose
c local
      integer(4):: is,ifmlw,ifuu,nq0i,nwf2,nqbz,nqbz2,iko_ix,iko_fx,
     &             iko_ix2,iko_fx2,iqbz,iqbz2,iq0i,iq0i2,j1,j2,
     &             ikp,ib,iwf,nbnd
      integer(4),allocatable:: isort(:)
      real(8):: q(3)
      real(8),allocatable:: wbnd(:)
      complex(8),allocatable:: dnk(:,:,:,:),uum(:,:,:,:,:),cbwf(:,:)

      if (Is_IO_Root_RSMPI()) then
c --- Readin MLWU/D and UUq0U/D
      do is = 1,nspin

c fileopen
         if (is.eq.1) then
           ifmlw  = iopen('MLWU',0,0,0)
           ifuu   = iopen('UUq0U',0,0,0)
         else ! is
           ifmlw  = iopen('MLWD',0,0,0)
           ifuu   = iopen('UUq0D',0,0,0)
         endif ! is

c nqbz mesh-points
         read(ifmlw)nqbz,nwf2,iko_ix,iko_fx
         if (nwf2.ne.nwf)
     >      call RSMPI_Stop("choose_wanband: nwf error")
!         if (nqbze.ne.nqbz)
!     >      call RSMPI_Stop("choose_wanband: nqbz error")
         if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nspin))
         do iqbz = 1,nqbz
            read(ifmlw)iqbz2,q(1:3)
            if (iqbz2 .ne. iqbz)
     &       call RSMPI_Stop('choose_wanband: iqbz error')
            read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
         enddo ! iqbz

c shifted mesh points
         read(ifuu)
         read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
         if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nspin))
         if (nqbz2.ne.nqbz)
     >      call RSMPI_Stop("choose_wanband: nqbz2 error")
         if (iko_ix2.ne.iko_ix)
     >      call RSMPI_Stop("choose_wanband: iko_ix2 error")
         if (iko_fx2.ne.iko_fx)
     >      call RSMPI_Stop("choose_wanband: iko_fx2 error")
         do iqbz = 1,nqbz
         do iq0i =1,nq0i  
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz)
     >       call RSMPI_Stop('choose_wanband: iqbz error')
            if (iq0i2 .ne. iq0i)
     >       call RSMPI_Stop('choose_wanband: iq0i error')
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
         enddo ! iq0i
         enddo ! iqbz

c fileclose
         if (is.eq.1) then
            ifmlw  = iclose('MLWU')
            ifuu   = iclose('UUq0U')
         else ! is
            ifmlw  = iclose('MLWD')
            ifuu   = iclose('UUq0D')
         endif ! is

      enddo ! is

      allocate(cbwf(iko_ix:iko_fx,nwf),wbnd(iko_ix:iko_fx),
     &         isort(iko_fx-iko_ix+1))
      do ikp = 1,nqbze
         iqbz = mod(ikp,nqbz)
         if (iqbz .eq. 0) iqbz = nqbz
         iq0i = (ikp - iqbz)/nqbz
      do is = 1,nspin
         cbwf = 0d0
         if (iq0i.eq.0) then
             cbwf(:,:) = dnk(:,:,iqbz,is)
         else ! iq0i
c   <psi(k+q0,n) | psi(k+q0,m)^W> 
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^W>  
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^W>  
c
c psi^W : bloch fn. in the Wannier gauge
           do ib = iko_ix,iko_fx
           do iwf= 1,nwf
             cbwf(ib,iwf) =
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
           enddo ! iwf
           enddo ! ib
         endif ! iq0i

c choose bands
         wbnd = 0.0d0
         nbnd = iko_fx - iko_ix + 1
         do ib = iko_ix,iko_fx
         do iwf = 1,nwf
           wbnd(ib) = wbnd(ib) + abs(cbwf(ib,iwf))**2
         enddo ! iwf
         enddo ! ib   
         call sortbnd(wbnd(iko_ix:iko_fx),iko_fx-iko_ix+1,
     &                isort)
         do iwf = 1,nwf
            iwbnd(iwf,ikp,is) = isort(iwf) + iko_ix - 1
         enddo ! iwf
      enddo ! is
      enddo ! ikp
      deallocate(dnk,uum,cbwf,wbnd,isort)

      endif ! Is_IO_Root_RSMPI

      call MPI_Bcast(iwbnd,nwf*nqbze*nspin,MPI_INTEGER,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(iwbnd)",ierror_rsmpi)

      end
C===================================================================
      subroutine sortbnd(rin,n,isort)
      implicit none
      integer(4):: n,isort(n)
      integer(4):: i,j,itmp
      real(8):: rin(n),r(n),rtmp

      do j = 1,n
         isort(j) = j
      enddo ! j   

      r(:) = rin(:)
      do j = 2,n
        rtmp = r(j)
        do i = j-1,1,-1
           if (r(i).ge.rtmp) goto 999
           r(i+1) = r(i)
           isort(i+1) = isort(i)
        enddo ! i
        i = 0
 999    continue
        r(i+1) = rtmp
        isort(i+1) = j
      enddo ! j

      do j = 1,n-1
         if (rin(isort(j)).lt.rin(isort(j+1)))
     &   stop "hx0fp0: sortbnd error"
      enddo   

      end
C===================================================================
      subroutine sortr(a,n,isort)
      implicit real*8(a-h,o-z)
      real(8) :: eps=1.0d-6
      real (8) :: a(n),b(n)
      integer (4) :: isort(n)

      b = a
      do i = 1,n
         isort(i) = i
      enddo   

      do j = 2,n
         c = b(j)
         do i = j-1,1,-1
            if (b(i)<=c) goto 999
            b(i+1) = b(i)
            isort(i+1) = isort(i)
         enddo
         i = 0
 999     b(i+1) = c
         isort(i+1) = j
      enddo

      do i = 1,n-1
         if (b(i).gt.b(i+1)) stop 'sortr: sorting error!'
      enddo   
      do i = 1,n
         if (abs(b(i)-a(isort(i))).gt.eps) stop 'sortr: sorting error!'
      enddo   

      return
      end
c-----------------------------------------------------------------------
      subroutine  sortvec2(ndat,vec,dist,idat)

      implicit real*8(a-h,o-z)
      real (8) :: vec(3,ndat),vtmp(3,ndat),dist(ndat)
      integer (4) :: idat(ndat)

      vtmp = vec
      do i = 1,ndat
         dist(i) = dsqrt(sum(vtmp(:,i)**2))
         idat(i) = i
      enddo

      do j = 2,ndat
         d = dist(j)
         do i = j-1,1,-1
            if (dist(i)<=d) goto 999
            dist(i+1) = dist(i)
            idat(i+1) = idat(i)
         enddo
         i = 0
 999     dist(i+1) = d
         idat(i+1) = j
      enddo

      do i = 1,ndat
         vec(1:3,i) = vtmp(1:3,idat(i))
      enddo   

      do i = 1,ndat-1
         d1 = dsqrt(sum(vec(:,i)**2))
         d2 = dsqrt(sum(vec(:,i+1)**2))
         if (d1 .gt. d2) stop 'sortvec: sorting error!'
         if (abs(d1-dist(i)).gt.1.d-4)
     &   stop 'sortvec: sorting error in d!'
      enddo   

      return
      end
c-----------------------------------------------------------------------
      subroutine q2q0(q,plat,q0)

      implicit real*8(a-h,o-z)
      parameter (eps=1d-4)
      real(8) :: q(3),q0(3),plat(3,3)

      do ii = 1,3
         q0(ii) = sum(q(:)*plat(:,ii))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine q02q(q0,qbas,q)

      implicit real*8(a-h,o-z)
      real(8) :: q(3),q0(3),qbas(3,3)

      q(:) = 0d0
      do ii = 1,3
         q(:) = q(:) + qbas(:,ii)*q0(ii)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine q02q0g0(qin,qout,ng)
      implicit real*8(a-h,o-z)
c
c qin = qout + ng

      parameter (eps=1d-6)
      real(8) :: qin(3), qout(3)
      integer(4) :: ng(3)

      do i = 1,3
         q = qin(i) - eps
         ng(i) = nint(q)
c         g = anint(q)
         qout(i) = q - dble(ng(i)) + eps
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine read_syml(qbandx,nqbandx,nqband)
      implicit real*8(a-h,o-z)
      integer(4) :: nqbandx,nqband,ifsyml,nline,nlinemax
      integer(4),allocatable :: nqq(:)
      real(8) :: qbandx(3,nqbandx),qqx(3)
      real(8),allocatable :: qq1(:,:),qq2(:,:)

        nlinemax = 50
        allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax))
        ifsyml = 3001
        open(ifsyml,file='SYML')
        nline = 0
        do 
        nline = nline + 1
        read(ifsyml,*,err=601,end=601) 
     &    nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
        enddo
  601   continue
        close(ifsyml)
        nline = nline - 1
        write(6,"(/' Symmetry lines:'/' points',12x,'start',22x,'end')")
        do is=1,nline
          write(6,"(i6,2x,3f8.4,2x,3f8.4)") 
     &     nqq(is),(qq1(i,is),i=1,3),(qq2(i,is),i=1,3)
        enddo
        nqnumx = sum(nqq(1:nline))
        iqq = 0
        do is = 1,nline
         nk = nqq(is)
        do iq=1,nk
          xx = 0d0
          if(nk>1) xx=(iq-1d0)/(nk-1d0)
          qqx = xx*qq2(1:3,is)+(1d0-xx)*qq1(1:3,is) 
c          if(iqq>1 ) then
c            if(abs(sum(qqx-qq_rsband(:,iqq)))<1d-10) cycle
c          endif
          iqq = iqq + 1
          qbandx(1:3,iqq) = qqx
          write (6,"('  q=',3f7.3)") qbandx(1:3,iqq)
        enddo
        enddo
        nqband = iqq
        write (6,"(' Total number of q-points:',i5/)") nqband
        deallocate(nqq,qq1,qq2)
      return
      end
c-----------------------------------------------------------------------
      subroutine findk(q1,qbz,plat,qbas,nqbz,iq2)
      implicit real*8(a-h,o-z)
      real(8) :: q1(3),qbz(3,nqbz),plat(3,3),qbas(3,3),
     &       dqx(3),dqx0(3),dq0(3),dq(3),ddq(nqbz)
      integer(4) :: ndg(3),isort(nqbz)

      do i = 1,nqbz
        dqx = q1 - qbz(:,i)
        call q2q0(dqx,plat,dqx0)
        call q02q0g0(dqx0,dq0,ndg)
        call q02q(dq0,qbas,dq)
        ddq(i) = sum(dq*dq)
      enddo  
      call sortr(ddq,nqbz,isort)
      iq2 = isort(1)

      return
      end
c-----------------------------------------------------------------------
      subroutine findk4(qband,qbz,plat,qbas,nqband,nqbz,wqk4,iqk4)
      implicit real*8(a-h,o-z)
      real(8) :: qband(3,nqband),qbz(3,nqbz),plat(3,3),qbas(3,3),q1(3),
     &       wqk4(4,nqband),dqx(3),dqx0(3),dq0(3),dq(3),ddq(nqbz),
     &        amat(4,4),amat2(4,4),bvec(4),vol(3,3),qtmp(3,4),ipiv(4)
      integer(4) :: iqk4(4,nqband),ndg(3),isort(nqbz)
      logical :: lddq

      eps = 1.0d-4
      nn = min(nqbz,1000)
      iqk4 = -1
      wqk4 = 0.0d0
      do iq = 1,nqband
        q1 = qband(:,iq)
        nerr = 0
! sorting
        do i = 1,nqbz
          dqx = q1 - qbz(:,i)
          call q2q0(dqx,plat,dqx0)
          call q02q0g0(dqx0,dq0,ndg)
          call q02q(dq0,qbas,dq)
          ddq(i) = sum(dq*dq)
        enddo  
        call sortr(ddq,nqbz,isort)
        iqk4(1:4,iq) = isort(1:4)
!
        lddq = .false.
        if (ddq(isort(1)).lt.eps) lddq = .true.
        if (lddq) then
          wqk4(1,iq) = 1.0d0
! main part
        else ! ddq  
 200      continue 
          do ii=1,4
            call q2qg(qbz(:,iqk4(ii,iq)),q1,plat,qbas,qtmp(:,ii))
            do jj=1,3
              amat(ii,jj)=qtmp(jj,ii)
            enddo ! jj
            amat(ii,4) = 1.d0
          enddo ! ii
c-----------------------------------------------------------------------
c     Check if the quadruplet of points form a tetrahedron with non-zero
c     volume
c
          do ii=1,3
          do jj=1,3
            vol(ii,jj) = amat(ii+1,jj) - amat(1,jj)
          enddo
          enddo
      avol = vol(1,1)*vol(2,2)*vol(3,3) + vol(1,2)*vol(2,3)*vol(3,1) +
     >       vol(1,3)*vol(2,1)*vol(3,2) - vol(1,1)*vol(2,3)*vol(3,2) -
     >       vol(1,2)*vol(2,1)*vol(3,3) - vol(1,3)*vol(2,2)*vol(3,1)
      if (abs(avol).lt.eps*eps*eps) then
 300     continue
         nerr = nerr + 1
         if (nerr.eq.(nn-3)**3) then
            write(6,*) ' warning: nn parameter in findk4 too small!'
            write(6,"(3f10.5)")q1
            write(6,"(i5,6f10.5)")iqk4(1,iq),qbz(:,iqk4(1,iq)),qtmp(:,1)
            write(6,"(i5,6f10.5)")iqk4(2,iq),qbz(:,iqk4(2,iq)),qtmp(:,2)
            write(6,"(i5,6f10.5)")iqk4(3,iq),qbz(:,iqk4(3,iq)),qtmp(:,3)
            write(6,"(i5,6f10.5)")iqk4(4,iq),qbz(:,iqk4(4,iq)),qtmp(:,4)
            write(6,"('amat',4f12.6)")amat(1:4,1)
            write(6,"('amat',4f12.6)")amat(1:4,2)
            write(6,"('amat',4f12.6)")amat(1:4,3)
            write(6,"('amat',4f12.6)")amat(1:4,4)
            write(6,600) avol,q1,nn,nerr
 600                format(4f9.4,2i10)
            write(6,*) (isort(jj),jj=1,nn)
c
c     impossible to form a tetrahedron with non-zero volume from the
c     q-points listed in isort; at this point, the last quadruplet
c     of q-points has the 1st, (nn-2)-th, (nn-1)-th and nn-th nearest
c     q-points
c
c            goto 500
             write(6,*) 'ERROR: Cannot form tetrahedron'
             write(6,*) '  The program stops...'
             stop
         endif ! nerr
c
c     change the 2nd, 3rd and/or 4th q-points in the quadruplet and try
c     to find a tetrahedron again
c
         ii = int( nerr/(nn-3)**2 )
         jj = int( (nerr - ii*(nn-3)**2)/(nn-3) )
         kk = nerr - ii*(nn-3)**2 - jj*(nn-3)
         if (ii.gt.jj.or.jj.gt.kk) goto 300
         iqk4(2,iq) = isort(2 + ii)
         iqk4(3,iq) = isort(3 + jj)
         iqk4(4,iq) = isort(4 + kk)
         goto 200
      endif ! abs(avol)
c-----------------------------------------------------------------------
      do jj=1,4
        amat2=amat
        bvec=0d0
        bvec(jj)=1d0
c
c  call lapack subroutine
c
         call DGESV(4,1,amat2,4,ipiv,bvec,4,info)
         if (info.ne.0) then
            write(6,*) 'failure to determine wqk4. info =',info,
     >           q1(:),iqk4(1:4,iq),nerr
            write(6,*) 'ERROR: interpolation'
            write(6,*) '  The program stops...'
            stop
c            goto 500
         endif
         asum=0.d0
         do ii=1,3
            asum = asum + bvec(ii)*q1(ii)
         enddo
         wqk4(jj,iq) = asum + bvec(4)
      enddo     ! jj=1,4
 500    continue
!
        endif ! ddq  
      enddo ! iq

check
      do iq = 1,nqband
         wsum = -1.0d0
      do i4 = 1,4
         if (iqk4(i4,iq).lt.1.or.iqk4(i4,iq).gt.nqbz)
     &   stop 'findk4: iqk4 error'
         wsum = wsum + wqk4(i4,iq)
      enddo ! i4
         if (abs(wsum).gt.eps) stop 'findk4: wqk4 error'
      enddo ! iq

      return
      end
c-----------------------------------------------------------------------
      subroutine q2qg(q,q0,plat,qbas,qg)
      implicit none
c
      double precision q(3),q0(3),qg(3),diff(3),p,eps,
     &   plat(3,3),qbas(3,3),qq(3),qq0(3),qqg(3)
      integer i
c
      eps=1d-4
c
      call q2q0(q,plat,qq)
      call q2q0(q0,plat,qq0)
      do i=1,3
         qqg(i) = qq(i) - idnint(qq(i)-qq0(i)) 
      enddo

      diff(:) = qqg(:)-qq0(:)
      do i=1,3
        if (abs(diff(i)).gt.0.5+eps) stop 'q2qg error'
      enddo 

      call q02q(qqg,qbas,qg)

      return
      end
c-----------------------------------------------------------------------
      subroutine readuuq0(is,iko_ix,iko_fx,nqbz,nq0i,uuq0)
      use rsmpi !RS
      implicit real*8(a-h,o-z)
      complex(8) :: uuq0(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i)

      if (is.eq.1) then
        ifuu      = iopen('UUq0U',0,0,0)
      else
        ifuu      = iopen('UUq0D',0,0,0)
      endif

         read(ifuu)
         read(ifuu)nqbz2,nq0i_2,iko_ix2,iko_fx2

         if (nqbz2.ne.nqbz)
     >      call RSMPI_Stop("readuuq0: nqbz2 error")
         if (nq0i_2.ne.nq0i)
     >      call RSMPI_Stop("readuuq0: nq0i_2 error")
         if (iko_ix2.ne.iko_ix)
     >      call RSMPI_Stop("readuuq0: iko_ix2 error")
         if (iko_fx2.ne.iko_fx)
     >      call RSMPI_Stop("readuuq0: iko_fx2 error")
         do iqbz = 1,nqbz
         do iq0i =1,nq0i  
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz)
     >       call RSMPI_Stop('readuuq0: iqbz error')
            if (iq0i2 .ne. iq0i)
     >       call RSMPI_Stop('readuuq0: iq0i error')
            read(ifuu)
     &      ((uuq0(j1,j2,iqbz,iq0i),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
         enddo ! iq0i
         enddo ! iqbz

      if (is.eq.1) then
        ifu = iclose('UUq0U')
      else
        ifu = iclose('UUq0D')
      endif

      return
      end
c-----------------------------------------------------------------------
