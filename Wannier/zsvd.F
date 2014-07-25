      subroutine zgesvdmn(ngb1,ngb2,zzz, SS,UU,VT)
      implicit none
      integer(4)::lwork,info,ngb1,ngb2,i
      complex(8):: zzz(ngb1,ngb2),UU(ngb1,ngb1),VT(ngb2,ngb2)
      real(8):: ss(ngb2)
      real(8),allocatable:: rwork(:)
      complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
      lwork=4*ngb1
      allocate(work(LWORK),rwork(5*ngb1)) 
      call zgesvd('A','A',ngb1,ngb2,zzz,ngb1,SS,UU,ngb1,VT,ngb2,work,lwork,rwork,info) 
      deallocate(work,rwork)
      end
!--------------------------------------------------------------------------
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
!--------------------------------------------------------------------------
