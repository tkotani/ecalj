c Information for plane-wave basis
      module m_QG
      implicit none
      integer ::  nqnum, ngpmx_qg, nnnn
      integer ::  nqnumc, ngcmx 
      double precision :: QpGcut_psi,QpGcut_Cou

      double precision,allocatable :: qqqa(:,:),qqqb(:,:)
      integer,allocatable :: ngp(:),ngvecp(:,:,:)
      integer,allocatable :: ngc(:),ngvecc(:,:,:)
      contains
cccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_QG()
      implicit none
      integer :: ikp
      integer :: ifiqg,ifiqgc

      write(6,*) '--- read_QG ---'
      ifiqg  = 101
      ifiqgc = 102

      open(ifiqg ,file='QGpsi',form='unformatted')
      open(ifiqgc,file='QGcou',form='unformatted')
      read(ifiqg  ) nqnum , ngpmx_qg, QpGcut_psi,nnnn
      read(ifiqgc ) nqnumc, ngcmx, QpGcut_Cou
      write(6,*) 'nqnum,nqnumc=',nqnum,nqnumc
      write(6,*) 'QpGcut_psi QpGcutCou nnnn='
     &   , QpGcut_psi,QpGcut_Cou ,nnnn
      if (nqnum.ne.nqnumc) then
        write(6,*) 'Error : nqnum!=nqnumc'
        write(6,*) 'nqnum,nqnumc=',nqnum,nqnumc
        stop 'Error : nqnum!=nqnumc'
      endif
      allocate(qqqa(3,nqnum),qqqb(3,nqnum))
      allocate(ngp(nqnum),ngc(nqnum))
      allocate( ngvecp(3,ngpmx_qg,nqnum),ngvecc(3,ngcmx,nqnum) )
      ngvecp(1:3,1:ngpmx_qg,1:nqnum)=0
      ngvecc(1:3,1:ngcmx,1:nqnum)=0
      do ikp = 1,nqnum
        read (ifiqg)  qqqa(1:3,ikp), ngp(ikp)
        read (ifiqgc) qqqb(1:3,ikp), ngc(ikp) 

c        write(6,"(i5,3f8.4,f10.5)")
c     &       ikp,qqqa(1:3,ikp),sum(abs(qqqa(1:3,ikp)-qqqb(1:3,ikp)))
        if (sum(abs(qqqa(1:3,ikp)-qqqb(1:3,ikp))).gt.1.0d-8) then
          stop 'qqqa!=qqqb'
        endif

        read (ifiqg ) ngvecp(1:3,1:ngp(ikp),ikp)
        read (ifiqgc) ngvecc(1:3,1:ngc(ikp),ikp)
      enddo
      close(ifiqg)
      close(ifiqgc)
      end subroutine read_QG
cccccccccccccccccccccccccccccccccccccccccccc
      end module m_QG

