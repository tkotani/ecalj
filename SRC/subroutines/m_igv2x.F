      module m_igv2x
      use m_struc_def,only: s_nv2
      public:: m_igv2xall_init, m_igv2x_setiq
      integer,protected,pointer,public :: igv2x(:,:)
      integer,protected,pointer,public :: napw,ndimh,ndimhx
      
      private
      integer,protected,private,target ::napw_z,ndimh_z,ndimhx_z
      type(s_nv2),allocatable,target,protected,private:: igv2xall(:)
      integer,allocatable,target,protected,private:: napwall(:),ndimhall(:),ndimhxall(:),igv2x_z(:,:)
      logical,private:: init=.true.
      contains
      
!!-----------------------------      
      subroutine m_Igv2x_setiq(iq) ! Get napw and so on for given qp
      integer:: iq
      napw  => napwall(iq) 
      ndimh => ndimhall(iq)
      ndimhx=> ndimhxall(iq)
      igv2x => igv2xall(iq)%nv2
      end
      
c$$$!!-----------------------------      
c$$$      subroutine m_Igv2x_set(qpin) ! Get napw and so on for given qp
c$$$      use m_qplist,only: qplist,iqini,iqend
c$$$      real(8)::qpin(3),qp(3)
c$$$      integer::iq,iqx
c$$$c      if(init) then
c$$$c         call m_igv2xall_init()
c$$$c         init=.false.
c$$$c      endif
c$$$      do iqx= iqini,iqend
c$$$         qp = qplist(:,iqx)
c$$$         if(sum((qp-qpin)**2)<1d-8) then
c$$$            iq=iqx
c$$$            exit
c$$$         endif
c$$$      enddo
c$$$      napw  => napwall(iq) 
c$$$      ndimh => ndimhall(iq)
c$$$      ndimhx=> ndimhxall(iq)
c$$$      igv2x => igv2xall(iq)%n
c$$$      end
!!-----------------------------      
      subroutine m_igv2xall_init()
      use m_qplist,only: qplist,iqini,iqend
      use m_MPItk,only: master_mpi,procid,master
      integer:: iq
      real(8):: qp(3)
      allocate(igv2xall(iqini:iqend),napwall(iqini:iqend),ndimhall(iqini:iqend),ndimhxall(iqini:iqend))
      do iq = iqini, iqend !This is a big iq loop
         qp = qplist(:,iq)
         call m_Igv2x_init(qp)   ! Get napw and so on for given qp
         allocate( igv2xall(iq)%nv2(3,napw_z) )
         igv2xall(iq)%nv2 = igv2x_z
         napwall(iq) = napw_z
         ndimhall(iq)= ndimh_z
         ndimhxall(iq)=ndimhx_z
c         write(6,*) 'xxxiiiiiiiiiqqqqqqqq',iqini,iqend,procid
cc         write(6,*) 'xxxiiiiiiiiiqqqqqqqq',iq,napw,ndimh,ndimhx
      enddo
      end subroutine
      
!!-----------------------------      
      subroutine m_igv2x_init(qp) !Set napw and igv2x_z for given qp
      use m_lmfinit,only: pwmode=>ham_pwmode,pwemin,pwemax,alat=>lat_alat,stdo,nspc
      use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
      use m_MPItk,only: master_mpi,procid,master
      use m_shortn3,only: shortn3_initialize,shortn3
      use m_lmfinit,only: ldham=>ham_ldham
      integer:: ifiese
      integer,allocatable ::  kv_iv(:,:) 
      real(8):: ppin(3),qp(3),qqq(3),pwgmin,pwgmax,dum
      logical:: debug,cmdopt0
      logical,save:: init=.true.
      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw,napwx
      debug    = cmdopt0('--debugbndfp')
c      print *,'mmmmmmm master_mpi=',master_mpi,procid,master
      if(pwmode>0.and.pwmode<10.and.init) then
         call shortn3_initialize(qlat)
         init=.false.
      endif
      
      if(pwmode>0.and.pwmode<10) then
         ppin=matmul(transpose(plat),qp) !basis on the qlat coordinate. qp in Cartesian.
         call shortn3(ppin,noutmx, nout,nlatout)
         if(debug) then
            do iout=1,nout
               write(*,"(a,3i5,f10.4,3f8.4)")'rrrrn1 =',nlatout(:,iout),
     &              sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
     &              matmul(qlat(:,:),ppin+nlatout(:,iout))
            enddo
         endif
      endif
      
      if(allocated(igv2x_z)) deallocate(igv2x_z)
      if (pwemax>0 .and. mod(pwmode,10)>0) then
         pwgmin = pwemin**.5d0
         pwgmax = pwemax**.5d0
         qqq = 0d0 
         if (mod(pwmode/10,10) == 1) qqq = qp !pwmode 1 in 10th digit means q-dependent nw
         call pshpr(1) !gvlst2 calls shortn3_initialize==>confusing(but initialization qlat is the same)
         call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,napwx,napwx,dum,dum,dum,dum) !get napw
         napw_z=napwx
         allocate(igv2x_z(3,napw_z), kv_iv(3,napw_z))
         call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napw_z,napw_z,kv_iv,dum,dum,igv2x_z)
         call poppr
         if (pwmode<10) then
            do iapw=1,napw_z
               igv2x_z(:,iapw)=igv2x_z(:,iapw)+nlatout(:,1)
            enddo
         endif
         deallocate(kv_iv)      !we only keep 
c         ndimh = ldim + napw    !ndimh is iq-dependent for pwmode >=10
c         if (ndimh .gt. ndham) then
c            call rx('bndfp:something wrong ndimh>ndham')
c         endif
      else
         napw_z=0 
         allocate(igv2x_z(1,1))   !dummy
      endif
      ndimh_z = ldham(1)+napw_z     ! ldim+napw
      if (mod(pwmode,10)==2) ndimh_z = napw_z !APW-only mode
      ndimhx_z = ndimh_z*nspc       !this is iq-dependent.


cccccccccccccccccccccc
c      if(allocated(igv2x)) deallocate(igv2x)
c      allocate(igv2x(3,napw_z))
c      igv2x=igv2x_z !xall(iq)%nv2
ccccccccccccccccccccc

      
      end
      end module
