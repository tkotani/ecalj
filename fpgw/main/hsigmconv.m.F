      program hsigmconv
c convert sigma to sigm
      implicit real*8 (a-h,o-z)
      real(8),allocatable ::  qqqx_m(:,:,:)
      complex(8), allocatable ::  sigma_m(:,:,:,:)
      integer(4) ::  n1,n2,n3, ifse_out
      ifse_out=3009

cccccccccccccccccccccccccccccccccccccccccccccccc
cswitch
c      go to 8888
cccccccccccccccccccccccccccccccccccccccccccccccccc

C ... read sigma ASCII format ---------------------------------------
! sigma to sigm
c      call headver('hsigmconv: convert sigma to sigm',0)
      open(UNIT=ifse_out, file='sigma')
      write(6,*) "=== Read sigma (ASCII) ========"
      call rwsigmaa2(ifse_out,sigma_m,qqqx_m, nspin,ntq,n1,n2,n3,nq,0)
      allocate(sigma_m(ntq,ntq,nq,nspin),qqqx_m(3,nq,nspin))
      rewind ifse_out
      call rwsigmaa2 (ifse_out,sigma_m,qqqx_m, nspin,ntq,n1,n2,n3,nq,1)
      close(ifse_out)
C ... write sigm Binary format ---------------------------------------
      open(UNIT=ifse_out, file='sigm',form='unformatted')
      write(6,*)
      write(6,*) "=== Write sigma to files (BINARY) ========"
      call rwsigma2 (-ifse_out,sigma_m,qqqx_m, nspin,ntq,n1,n2,n3,nq,1)
      close(ifse_out)
Cstop2rx 2013.08.09 kino 9999 stop ' OK! hsigmconv: convert from sigma to sigm '
 9999 call rx0( ' OK! hsigmconv: convert from sigma to sigm ')


 8888 continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! sigm to sigma
C ... read sigm binary format ---------------------------------------
      open(UNIT=ifse_out, file='sigm',form='unformatted')
      write(6,*) "=== Read sigma (BINARY) ========"
      call rwsigma2(ifse_out,sigma_m,qqqx_m, nspin,ntq,n1,n2,n3,nq,0)
      allocate(sigma_m(ntq,ntq,nq,nspin),qqqx_m(3,nq,nspin))
      rewind ifse_out
      call rwsigma2 (ifse_out,sigma_m,qqqx_m, nspin,ntq,n1,n2,n3,nq,1)
      close(ifse_out)
C ... write sigm ascii format ---------------------------------------
      write(6,*)
      write(6,*) "=== Write sigma to files (ASCII) ========"
      open(UNIT=ifse_out, file='sigma')
      call rwsigmaa2 (-ifse_out,sigma_m,qqqx_m, nspin,ntq,n1,n2,n3,nq,1)
      close(ifse_out)
Cstop2rx 2013.08.09 kino 9998 stop ' OK! hsigmconv: convert from sigm to sigma '
 9998 call rx0( ' OK! hsigmconv: convert from sigm to sigma ')
      end

C----------------------------------------------------------------------
      subroutine rwsigma2(ifs,sigma_m,qqqx_m,nspin,ntq,n1,n2,n3,nq,job)
      complex(8)::sigma_m(ntq,ntq,nq,nspin)
      real(8)::qqqx_m(3,nq,nspin),sumchk=0d0,qin(3)
      logical :: qlistchk,test
cccccccccccccccccccccccccdddddddddddddddddccccc
      logical :: qlistcase=.true.
ccccccccccccccccccccccccccccccccccccccdddddddddd
      ntm1=0;ntm2=0;ntm3=0
      nqx=nq

cccccccccccccccccccccccccc
      if(qlistcase) then
        do is=1,1
          nqx=0
          do ip=1,nq
            if(.not.qlistchk(qqqx_m(1:3,ip,is)) ) cycle
            nqx=nqx+1
          enddo
        enddo
        write(6,*)'nqx=',nqx
      endif
cccccccccccccccccccccccccccccccc

      if(ifs>0) print *, " rwsigma: Reading sigm (Binary)==="
      if(ifs<0) print *, " rwsigma  Writing sigm (Binary)==="
      if(ifs>0) read ( ifs) nspin,ntq,n1,n2,n3,nq, ntm1,ntm2,ntm3
      if(ifs<0) write(-ifs) nspin,ntq,n1,n2,n3,nqx, ntm1,ntm2,ntm3
      if(ifs>0) print *,' nspin ntq nq=', nspin,ntq,nq
      if(job==0) return
      do is=1,nspin
        do ip=1,nq
ccccccccccccccccccccccccccccccc
          if(qlistcase) then
            if(.not.qlistchk(qqqx_m(1:3,ip,is)) ) cycle
          endif
cccccccccccccccccccccccccc
          if(ifs>0) read (ifs)  qqqx_m(1:3,ip,is)
          if(ifs>0) read (ifs)  sigma_m(1:ntq,1:ntq,ip,is)
          if(ifs<0) write(-ifs) qqqx_m(1:3,ip,is)
          if(ifs<0) write(-ifs) sigma_m(1:ntq,1:ntq,ip,is)
          write(6,"('  === ',i5,i3,3f10.5,d24.16)")ip,is,qqqx_m(1:3,ip,is)
     &  ,sum(abs(sigma_m(1:ntq,1:ntq,ip,is)))
        enddo
      enddo
      print *,' === rwsigma:  sum check of sigma_m=',sum(abs(sigma_m))
      print *
      end
c
      subroutine rwsigmaa2(ifs,sigma_m,qqqx_m,nspin,ntq,n1,n2,n3,nq,job)
      complex(8)::sigma_m(ntq,ntq,nq,nspin)
      real(8)::qqqx_m(3,nq,nspin)
      integer(4):: job
      logical :: qlistchk
      ntm1=0;ntm2=0;ntm3=0
      if(ifs>0) print *, " rwsigmaa: Reading sigma (ASCII)==="
      if(ifs<0) print *, " rwsigmaa: Writing sigma (ASCII)==="
      if(ifs>0)read ( ifs,*)      nspin,ntq,n1,n2,n3,nq,ntm1,ntm2,ntm3
      if(ifs<0)write(-ifs,"(i10)")nspin,ntq,n1,n2,n3,nq,ntm1,ntm2,ntm3
      if(ifs>0) print *,' nspin ntq nq=', nspin,ntq,nq
      if(job==0) return
      do is=1,nspin
        do ip=1,nq
          if(ifs>0) read (ifs,*)  qqqx_m(1:3,ip,is)
          if(ifs<0) write(-ifs,"(3d24.16)") qqqx_m(1:3,ip,is)
          do j=1,ntq
            do i=1,ntq
              if(ifs>0) read (ifs,"(2d24.16)")  sigma_m(i,j,ip,is)
c          write(6,*) i,j,sigma_m(i,j,ip,is)
              if(ifs<0) write(-ifs,"(2d24.16,2i5)") sigma_m(i,j,ip,is),i,j
            enddo
          enddo
          write(6,"('  === ',i5,i3,3f10.5,d24.16)")ip,is,qqqx_m(1:3,ip,is)
     &  ,sum(abs(sigma_m(1:ntq,1:ntq,ip,is)))
        enddo
      enddo
      print *,' === rwsigmaa: sum check of sigma_m=',sum(abs(sigma_m))
      print *
      end


      logical function qlistchk(q)
      real(8),save:: qx(3,1000)
      real(8):: q(3)
      integer(4),save:: ifi,nx
      logical :: init=.true.
      if(init) then
        ifi=111
        open(ifi,file='qlist')
        ix=0
        do
          ix=ix+1
          read(ifi,*,end=101) qx(:,ix)
          write(6,*)qx(:,ix)
        enddo
 101    continue
        nx=ix-1
        print *,' nx=',nx
        init=.false.
      endif
      qlistchk=.true.
      do ix=1,nx
        if( sum(abs(q-qx(:,ix)))<1d-4 ) then
          write(6,'("    Find q in list q=",3f13.5)') q
          return
        endif
      enddo
c      write(6,'(" Not Find q in list q=",3f13.5)')q
      qlistchk=.false.
      end

