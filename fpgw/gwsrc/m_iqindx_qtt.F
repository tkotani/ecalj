      module m_iqindx_qtt
!! To accelarate inverse mapping, q to 'integer index',
!! we prepare integer index of q, its resolution is given by epsd.
      use m_hamindex,only: qtt,nqtt
      implicit none

      public:: iqindx2_, init_iqindx_qtt
      private
      integer,allocatable,private:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
      integer,private:: nkey(3)
      real(8),private:: ginv(3,3)
      real(8),private:: epsd=1d-7 !key parameter to map to integer index

      contains

      subroutine iqindx2_(q,  iqindx,qu)
      intent(in)::        q
      intent(out)::           iqindx,qu
!> Find index as q=qq(:,iq) with modulo of premitive vector.===
!! ginv is the inverse of plat (premitive translation vector).
!! Use kk1,kk2,kk3,nkey(1:3),iqkkk to get iqindx.
      real(8) :: q(3)
      integer :: iqindx
      real(8) :: qu(3)
      integer:: i_out, iq,iqx ,kkk3(3),ik1,ik2,ik3
      real(8):: qx(3),qzz(3)
      logical::debug=.true.

ccccccccccccccccccccc
      debug=.false.
      if(abs(q(1)+0.1d0)+abs(q(2)+0.1d0)<1d-3) then
         debug=.true.
      endif
ccccccccccccccccccccc
      if(debug) write(6,"(' iqindx2_: q=',3f20.15)") q
      call rangedq(matmul(ginv,q), qzz)
      if(debug) write(6,"(' iqindx2_: q=',3f20.15)") qzz
!! we generate qzz integer index for qzz
      kkk3 = (qzz+0.5d0*epsd)/epsd 
      if(debug) write(6,*)'kkk3=',kkk3
      if(debug) write(6,*)'nkey=',nkey
      if(debug) write(6,*)'kk1=',kk1
      if(debug) write(6,*)'kk2=',kk2
      if(debug) write(6,*)'kk3=',kk3
      call tabkk(kkk3(1), kk1,nkey(1), ik1)
      call tabkk(kkk3(2), kk2,nkey(2), ik2)
      call tabkk(kkk3(3), kk3,nkey(3), ik3)
      if(debug) write(6,"(' 222222a q=',3i8,3f18.12)") kkk3,qzz
      if(debug) write(6,"(' 222222a ik1,ik2,ik3,q=',3i8,3f18.12)") ik1,ik2,ik3,q

      iqindx = iqkkk(ik1,ik2,ik3)
      if(debug) then
         do iqx=1,nqtt
         write(6,"(i5,3f13.5)")iqx,qtt(:,iqx)
         enddo
      endif
      qu =qtt(:,iqindx)
      if(debug) write(6,*) iqindx,qu
      end subroutine
!!----
      subroutine init_iqindx_qtt(ginv_)
      intent(in)::               ginv_
!! === mapping of qtt ===
!! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
!!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
!!       --->  ik1,ik2,ik3= tabkk(kkk,iqk,nkey) ---> iqkkk(ik1,ik2,ik3)
      real(8):: ginv_(3,3),qzz(3)
      real(8),allocatable:: qxx(:,:)
      integer:: isig,i,ix,kkk,kkk3(3),ik1,ik2,ik3,iq,ik
      integer,allocatable:: ieord(:)
      logical::debug=.false.
      ginv=ginv_
      allocate(ieord(nqtt))
      if(debug) write(6,"(a,2i5,20f9.4)")' iiiiii nqtt=',nqtt,size(qtt),ginv_(1:3,1:3),ginv(1:3,1:3)
      allocate(key(3,0:nqtt),qxx(3,nqtt))
      key=-99999
      do iq=1,nqtt
         call rangedq(matmul(ginv,qtt(:,iq)), qxx(:,iq))
c         write(6,"(a,i5,3f13.5,2x,3f13.5)") ' qqqttxx =',iq,qtt(:,iq),qxx(:,iq)
      enddo
c      write(6,*)'sssqqq=',sum(abs(qtt(1,1:nqtt))),      sum(abs(qtt(2,1:nqtt))),   sum(abs(qtt(3,1:nqtt)))

!! get key and nkey for each ix.
      key(:,0)=0 !dummy
      do ix =1,3
         call sortea(qxx(ix,:),ieord,nqtt,isig)
         ik=0
         do i=1,nqtt
            kkk=(qxx(ix,ieord(i))+0.5d0*epsd)/epsd  !kkk is digitized by 1/epsd
            if(i==1.or.key(ix,ik)<kkk) then
               ik=ik+1
               key(ix,ik) = kkk
            elseif (key(ix,ik)>kkk) then
               write(6,*)ix, ik,i, key(ix,ik), qxx(ix,ieord(i))
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
      iqkkk=-99999
c      print *,'eeee kk1=',kk1
c      print *,'eeee kk2=',kk2
c      print *,'eeee kk3=',kk3
c      write(6,*)' initqqq nqtt=',nqtt
      do i=1,nqtt
         kkk3= (qxx(:,i)+0.5d0*epsd)/epsd !kkk is digitized by 1/epsd
         call tabkk(kkk3(1), kk1,nkey(1), ik1)
         call tabkk(kkk3(2), kk2,nkey(2), ik2)
         call tabkk(kkk3(3), kk3,nkey(3), ik3)
         iqkkk(ik1,ik2,ik3)=i
c         print *,'111111 kkk3=',kkk3,qxx(:,i)
c         write(6,"(' 111111111a ik1,ik2,ik3,iqxx=',3i8,3f18.12)") kkk3,qxx(:,i)
c         write(6,"(' 111111111a ik1,ik2,ik3,iqkkk=',4i8,3f18.12)") ik1,ik2,ik3,i,qtt(:,i)
      enddo
      deallocate(qxx)
c      stop 'xxxxxxxxxxxxx'
      end subroutine

!! ---
      subroutine tabkk(kkin, kktable,n, nout)
      intent(in)::     kkin, kktable,n
      intent(out)::                     nout
      integer:: nout,n, kkin, kktable(n),i,mm,i1,i2
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
      write(6,*) 'xxxxx takk ', i1,i2,kkin
      write(6,*) 'xxxxx takk ',kktable(i1),kktable(i2)
      call rx( 'takk: error')
      end subroutine
      end module
