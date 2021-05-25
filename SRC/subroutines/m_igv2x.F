      module m_igv2x
c      public:: m_igv2x_set,napw,igv2x
      integer,allocatable,protected:: igv2x(:,:)
      integer,protected::napw,ndimh,ndimhx
      
c      private
      contains
      subroutine m_igv2x_set(qp) !Set napw and igv2x for given qp
      use m_lmfinit,only: pwmode=>ham_pwmode,pwemin=>ham_pwemin,pwemax=>ham_pwemax,alat=>lat_alat,stdo,nspc
      use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
      use m_MPItk,only: master_mpi,procid,master
      use m_shortn3,only: shortn3_initialize,shortn3
      use m_suham,only: ldham=>ham_ldham
      integer:: ifiese
      integer,allocatable ::  kv_iv(:,:) 
      real(8):: ppin(3),qp(3),qqq(3),pwgmin,pwgmax,dum
      logical:: init=.true.,debug,cmdopt0
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
      
      if(allocated(igv2x)) deallocate(igv2x)
      if (pwemax>0 .and. mod(pwmode,10)>0) then
         pwgmin = pwemin**.5d0
         pwgmax = pwemax**.5d0
         qqq = 0d0 
         if (mod(pwmode/10,10) == 1) qqq = qp !pwmode 1 in 10th digit means q-dependent nw
         call pshpr(1) !gvlst2 calls shortn3_initialize==>confusing(but initialization qlat is the same)
         call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,napwx,napwx,dum,dum,dum,dum) !get napw
         napw=napwx
         allocate(igv2x(3,napw), kv_iv(3,napw))
         call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napw,napw,kv_iv,dum,dum,igv2x)
         call poppr
         if (pwmode<10) then
            do iapw=1,napw
               igv2x(:,iapw)=igv2x(:,iapw)+nlatout(:,1)
            enddo
         endif
         deallocate(kv_iv)      !we only keep 
c         ndimh = ldim + napw    !ndimh is iq-dependent for pwmode >=10
c         if (ndimh .gt. ndham) then
c            call rx('bndfp:something wrong ndimh>ndham')
c         endif
      else
         napw=0 
         allocate(igv2x(1,1))   !dummy
      endif
      ndimh = ldham(1)+napw     ! ldim+napw
      if (mod(pwmode,10)==2) ndimh = napw !APW-only mode
      ndimhx = ndimh*nspc       !this is iq-dependent.
      end
      end module
