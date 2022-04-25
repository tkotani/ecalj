      subroutine zgesvdnn(ngb,zzz, SS,UU,VT)
! SVD zzz= matmul(UU,matmul(SS,VT)) ------------
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

