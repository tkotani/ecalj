!! read BZDATA
      module m_read_bzdata
!! we set following data when you call read_BZDATA() 
      implicit none
      integer,allocatable :: ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,nqibz_r
      integer,allocatable :: n1,n2,n3,
     &    idtetf(:,:),ib1bz(:),idteti(:,:),
     &    nstar(:),irk(:,:),nstbz(:) !,index_qbz(:,:,:)
      real(8),allocatable:: qbas(:,:),ginv(:,:),qbasmc(:,:),dq_bzcase2(:)
     &    , qbz(:,:),wbz(:),qibz(:,:)
     &    , wibz(:),qbzw(:,:), qibz_r(:,:)
      logical:: done_read_bzdata=.false.
c-----------------------------------------------------------------
      contains
      subroutine read_BZDATA()
!! No input except BZDATA file.
!! After you call this, you can access Brillowin Zone datas above ----
      implicit none
      integer(4):: intq(3),iqbz,ifbz,n,verbose
      real(8) :: qout(3)
      ifbz = 6661
      write(6,*)' ### readin BZDATA ###'
      allocate(n1,n2,n3,ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,nqibz_r)
      allocate(qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2(3))
      open (ifbz, file='BZDATA')
      read(ifbz,*)  nqbz,nqibz, nqbzw, ntetf, nteti, ngrp, nqibz_r
      read(ifbz,*)  n1,n2,n3
      write(6,"(' read_bzdata',10i5)")nqbz,nqibz, nqbzw, ntetf, nteti,ngrp
      allocate(qbz(3,nqbz),wbz(nqbz))
      allocate(qibz(3,nqbz),wibz(nqbz),nstbz(nqbz))
      allocate(nstar(nqibz), irk(nqibz,ngrp))
      if(ntetf>0) then
        allocate( idtetf(0:3,ntetf), ib1bz(nqbzw), qbzw(3,nqbzw) )
      endif
      if(nteti>0) then
        allocate( idteti(0:4,6*nqbz))
      endif
      allocate( qibz_r(3,nqibz_r) )
      call rwbzdata(ifbz,1,
     &    ngrp,qbas,ginv,qbasmc,
     i    qbz, wbz,nstbz,           nqbz, 
     i    qibz,wibz, nstar,irk,  nqibz,
     i    idtetf, ntetf, qbzw,ib1bz, nqbzw,
     i    idteti, nteti,dq_bzcase2, qibz_r, nqibz_r )
      close(ifbz)

C ... Add index to specify qbz. index_qbz is used to find index in qbz.
C This algolism here can make index_qbz for any qbz(1:3,nqbz).
c With fbz2, you can get index iq for given q(1:3) so that
c            q(1:3)= qbz(1:3,iq) + some G vector.
c
c To find iqbz for given q,
c      call rangedq(matmul(ginv,q), qout)
c      intq =  qout*n +1
c      iqbz= index_qbz(intq(1),intq(2),intq(3))
c See fbz2.

c         n = 0
c         allocate(index_qbz(1,1,1)) !dummy
c 1120    continue
c         n = n+1
c         print *,' =========================== n=',n
c         deallocate(index_qbz)
c         allocate(index_qbz(n,n,n))
c         index_qbz = -9999
c         do iqbz = 1,nqbz
c           call rangedq(matmul(ginv,qbz(:,iqbz)), qout)
c           intq =  qout*n +1
cc           print *,' qbz=', qbz(:,iqbz) !,matmul(ginv,qbz(:,iqbz))
cc           print *,' qout='
cc           print *, qout
cc           print *, ' intq=', intq
c           if(verbose()>=100) write(6,"(' qbz=',3f10.5,'  index=',3i3,'   qout=',3f10.5)")
c     &       qbz(:,iqbz),intq,qout
c           call checkrange(intq(1),1,n) !sanity checks
c           call checkrange(intq(2),1,n)
c           call checkrange(intq(3),1,n)
c           if(index_qbz(intq(1),intq(2),intq(3))/=-9999) then
c              if(verbose()>=100) print *,'failed indexing with this n. Try to enlarged n'
c              goto 1120
c           endif
c           index_qbz(intq(1),intq(2),intq(3)) = iqbz
c         enddo
c         n_index_qbz=n
c         write(6,*) " O.K. index_qbz is generated; n=", n

      if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
        print *, 'sum (wibz)=', sum(wibz(1:nqibz))
Cstop2rx 2013.08.09 kino        stop 'read_BZDATA  sum (wibz) is not 2.'
        call rx( 'read_BZDATA  sum (wibz) is not 2.')
      endif
      done_read_bzdata=.true.
      end subroutine
      end module

      subroutine rwbzdata(ifbz,job,
     x  ngrp,qbas,ginv,qbasmc,
     x  qbz, wbz  ,nstbz, nqbz, 
     x  qibz,wibz, nstar,irk,  nqibz,
     x  idtetf, ntetf, qbzw,ib1bz, nqbzw,
     x  idteti, nteti,dq_bzcase2, qibz_r,nqibz_r )
C- ReadWrite  BZ mesh data reuired for GW ---
C----------------------------------------------------------------
      implicit none
      integer(4):: nqbz,n1q,n2q,n3q,ntetf,nteti,nqbzw,iqbz,ifbz
     & ,nqibz,iqibz,itet,ngrp,job,nqibz_r
      real(8)   :: plat(3,3),vol,qbasmc(3,3)
      integer(4),allocatable:: ipq(:),iw1(:)
      real(8):: qbz(3,nqbz),wbz(nqbz),qibz(3,nqibz),wibz(nqibz)
     &         ,qbzw(3,nqbzw),qbas(3,3),ginv(3,3),dq_bzcase2(3), qibz_r(3,nqibz_r)
      integer(4):: idtetf(0:3,ntetf),ib1bz(nqbzw),idteti(0:4,nteti)
     &       ,irk(nqibz,ngrp),nstar(nqibz),nstbz(nqbz),bzcase
      logical tetraf,tetrai

      if(job<0) write(ifbz,"(3d24.16)") qbas,ginv,qbasmc
      if(job>0) read (ifbz,"(3d24.16)") qbas,ginv,qbasmc

      do iqibz = 1,nqibz
        if(job<0) then
          write(ifbz,"(4d24.16,i9)") 
     &    qibz(1:3,iqibz),wibz(iqibz),nstar(iqibz)
          write(ifbz,"(100i8)") irk(iqibz,1:ngrp)
        else
          read(ifbz,"(4d24.16,i9)") 
     &    qibz(1:3,iqibz),wibz(iqibz),nstar(iqibz)
          read(ifbz,"(100i8)") irk(iqibz,1:ngrp)
        endif
      enddo

      if(job<0) then
        write(ifbz,"(i10)") nqibz_r
      else
        read (ifbz,"(i10)") nqibz_r
      endif
      do iqibz = 1,nqibz_r
        if(job<0) then
          write(ifbz,"(3d24.16)") qibz_r(1:3,iqibz)
        else
          read (ifbz,"(3d24.16)") qibz_r(1:3,iqibz)
        endif
      enddo

      do iqbz = 1,nqbz
        if(job<0) then
          write(ifbz,"(4d24.16,i10)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
        else
          read(ifbz, "(4d24.16,i10)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
! bug fix (kino); May 2004 it was write(ifbz,"(4d24.16)") read(ifbz,"(4d24.16)")
! But it was working in ifc
!   (the strange number in ifbz file was correctly recoverd to the original values!)
        endif
ccccccccccccccccccccccccccccccccccc
c         write(6, "(4d24.16,i5)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
ccccccccccccccccccccccccccccccccccc
      enddo
ccccccccccccccccccccccccccccccccccc
c      stop 'test end'
ccccccccccccccccccccccccccccccccccc
      if(ntetf>0) then
        if(job<0) then
          write(ifbz,"(4i10)") (idtetf(0:3,itet),itet=1,ntetf)
          write(ifbz,"(i9,3d24.16)") 
     &       (ib1bz(iqbz), qbzw(1:3,iqbz),iqbz=1,nqbzw)
        else
          read(ifbz,"(4i10)") (idtetf(0:3,itet),itet=1,ntetf)
          read(ifbz,"(i9,3d24.16)") 
     &       (ib1bz(iqbz), qbzw(1:3,iqbz),iqbz=1,nqbzw)
        endif
      endif
      if(nteti>0) then
        if(job<0) then
          write(ifbz,"(5i10)") (idteti(0:4,itet),itet=1,nteti)
        else
          read(ifbz,"(5i10)") (idteti(0:4,itet),itet=1,nteti)
        endif
      endif

c      if(bzcase()==2) then
      if(job<0) write(ifbz,"(3d24.16,' !dq_bzcase2')") dq_bzcase2
      if(job>0) read (ifbz,"(3d24.16)") dq_bzcase2
c      endif

      end

      subroutine checkrange(intq,n1,n2)
      if(intq<n1 .or. intq>n2) then
        print *,'checkrange: intq n1 n2= ',intq,n1,n2
Cstop2rx 2013.08.09 kino        stop    'checkrange: stop '
        call rx( 'checkrange: stop ')
      endif
      end


