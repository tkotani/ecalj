!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
C---- get size
      subroutine getngbpomat(nqibze
     o  ,nnmx,nomx)
C- just to get the maximum size of ngb (mized basis size) from POmat
      implicit none
      integer(4):: ifpomat,nnmx,ikpo,nn_,noo,iqxxx,isx,
     &      nomx,nqibze,iopen,iclose
      real(8):: q_r(3)
      complex(8),allocatable:: pomat(:,:)
      open(newunit=ifpomat,file='POmat',form='unformatted')
      nnmx=0
      nomx=0
      do ikpo=1,nqibze
        read(ifpomat) q_r,nn_,noo,iqxxx !readin reduction matrix pomat
        write(6,"('smbasis: ikp q no nn=',i5,3f8.4,4i5)") ikpo,q_r, noo,nn_
        if(nn_>nnmx) nnmx=nn_
        if(noo >nomx) nomx=noo
        allocate( pomat(nn_,noo) )
        read(ifpomat) pomat
        deallocate(pomat)
      enddo
      close(ifpomat) 
      end

!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine qqsave (qi,nmax,ginv, qsave,imx)
C---accumulate q into qsave imx
c        if(allocated(qsave)) deallocate(qsave)
c        allocate( qsave(3,nmax))
c     imx=0
      integer:: nmax,i,imx
      real(8):: qi(3),qsave(3,nmax),qx(3),ginv(3,3)
      real(8):: tolq=1d-8
      do i = 1,imx
        if(sum(abs(qi-qsave(:,i)))<tolq) return
        call rangedq( matmul(ginv,qi-qsave(:,i)), qx)
        if(sum(abs(qx))< tolq) return
      enddo
      imx = imx+1
      if(imx>nmax) call rx( ' qqsave: imx>=nmax')
      qsave(:,imx) = qi
      end

!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine readd_iSigma_en(ifinin,iSigma_en)
      use m_keyvalue,only:getkeyvalue
      integer(4):: iSigma_en,ifinin
      call getkeyvalue("GWinput","iSigMode",iSigma_en )
      write(6,*)' iSigma_en=',iSigma_en
      end
      
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine getnemx(nbmx,ebmx,im,ipr)
C- Readin nbmx ebmx for hxofp0 hscfp0
      use m_keyvalue,only: getkeyvalue
      real(8)::ebmx
      integer (4)::nbmx,ret,im!,ifinin
      character(len=100):: recxxx=' '
      logical :: ipr !,readgwinput
c      if(readgwinput()) then
      if    (im==8) then
        call getkeyvalue("GWinput","nband_sigm",nbmx, default=99999 )
        call getkeyvalue("GWinput","emax_sigm", ebmx, default=1d10  )
      elseif(im==7) then
        call getkeyvalue("GWinput","nband_chi0",nbmx, default=99999 )
        call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
      endif
      if(ipr) write(6,"('  nbmx ebmx from GWinput=',i10,d13.6)") nbmx,ebmx
      return
      end

!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine getnemx8(nbmx,ebmx)
C- Readin nbmx ebmx for hscfp0
      use m_keyvalue,only: getkeyvalue
      integer(4):: ret
      real(8)::ebmx(2)
      integer (4)::nbmx(2)
c      character(len=100):: recxxx=' '
      call getkeyvalue("GWinput","nband_sigm",nbmx,1, default=(/9999999/),status=ret)
      write(6,*)' status 1=',ret
      call getkeyvalue("GWinput","emax_sigm", ebmx,1, default=(/1d10/),status=ret)
      write(6,*)' status 2=',ret
      return
      end
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine readin5(i0,i1,i2)
c     ! readin i0,i1,i2; these defaults value are 0 0 0 if these are not given.
      integer:: i0,i1,i2
      character(len=100):: recxxx
      character(len=106):: recxxx2
      read (5,"(a100)",end=1100) recxxx
 1100 continue
      recxxx2 = recxxx//' 0 0 0'
      read(recxxx2,*) i0, i1, i2
      end
!sssssssssssssssssssssssssssssssssssssss
      subroutine readin6(i0,i1,i2,i3)
c     ! readin i0,i1,i2; these defaults value are 0 0 0 if these are not given.
      integer:: i0,i1,i2,i3
      character(len=100):: recxxx
      character(len=106):: recxxx2
      read (5,"(a100)",end=1100) recxxx
 1100 continue
      recxxx2 = recxxx//' 0 0 0 0'
      read(recxxx2,*) i0, i1, i2,i3
      end
