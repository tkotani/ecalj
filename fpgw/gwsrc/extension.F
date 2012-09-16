C---- get size
      subroutine getngbpomat(nqibze
     o  ,nnmx,nomx)
C- just to get the maximum size of ngb (mized basis size) from POmat
      implicit none
      integer(4):: ifpomat,nnmx,ikpo,nn_,noo,iqxxx,isx,
     &      nomx,nqibze,iopen,iclose
      real(8):: q_r(3)
      complex(8),allocatable:: pomat(:,:)
      ifpomat = iopen('POmat',0,-1,0) !oct2005
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
      isx = iclose('POmat')
      end


      subroutine qqsave (qi,nmax,ginv, qsave,imx)
C---accumulate q into qsave imx
c        if(allocated(qsave)) deallocate(qsave)
c        allocate( qsave(3,nmax))
c        imx=0
      real(8):: qi(3),qsave(3,nmax),qx(3),ginv(3,3)
     . ,tolq=1d-5,tolq2=1d-4
      do i = 1,imx
        if(sum(abs(qi-qsave(:,i)))<tolq) return
c 27July2005
        call rangedq( matmul(ginv,qi-qsave(:,i)), qx)
        if(sum(abs(qx))< tolq) return
      enddo
      imx = imx+1
      if(imx>nmax) stop ' qqsave: imx>=nmax'
      qsave(:,imx) = qi
      end

      subroutine readd_iSigma_en(ifinin,iSigma_en)
      use keyvalue
      integer(4):: iSigma_en
c      logical :: readgwinput
c      if(readgwinput()) then
      call getkeyvalue("GWinput","iSigMode",iSigma_en )
c      else
c        ifinin=iopen('GWIN_V2',1,0,0)
c        read(ifinin,*)
c        read(ifinin,*)
c        read(ifinin,*) blank,blank,iSigma_en
cc        iSigma_en parameter which determines approximation for  self-energy,
cc        stored in GWIN_V2.
cc        iSigma_en==0 SE_nn'(ef)+imag intgrl:delta_nn'(SE_nn(e_n)-SE_nn(ef))
cc        iSigma_en==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
cc        iSigma_en==2 SE_nn'((e_n+e_n')/2)
cc        iSigma_en==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2
cc        iSigma_en==5 delta_nn' SE_nn(e_n)
cc        output file should contain hermitean part of SE for energies to be real
cc        (for example, hermitean conjunction of SE_nn'(e_n) means SE_n'n(e_n')^* )
c        is=iclose('GWIN_V2')
c      endif
      print *,' iSigma_en=',iSigma_en
      end
c--------------------------
      subroutine getnemx(nbmx,ebmx,im,ipr)
C- Readin nbmx ebmx for hxofp0 hscfp0
      use keyvalue
      real(8)::ebmx
      integer (4)::nbmx,ret
      character(len=100):: recxxx=' '
      logical :: ipr !,readgwinput
      ifinin = 8087
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
c      else
c        nbmx=999999
c        ebmx=1d10
c        open(ifinin,file='GWIN0')
c        do i=1,im; read(ifinin,*); enddo
c        read (ifinin,"(a100)",err=1001,end=1001) recxxx
c 1001   continue
c        read(recxxx,*,err=1012,end=1012) nbmx
c 1012   continue
c        read(recxxx,*,err=1013,end=1013) nbmx,ebmx
c 1013   continue
c        if(ipr) write(6,"('  nbmx ebmx from GWIN0=',i10,d13.6)") nbmx,ebmx
c        close(ifinin)
c      endif
      end

c--------------------------
      subroutine getnemx8(nbmx,ebmx)
C- Readin nbmx ebmx for hscfp0
      use keyvalue
      integer(4):: ret
      real(8)::ebmx(2)
      integer (4)::nbmx(2)
c      character(len=100):: recxxx=' '
      call getkeyvalue("GWinput","nband_sigm",nbmx,2, default=(/9999999,999999/),status=ret)
      if(ret<0) call getkeyvalue("GWinput","nband_sigm",nbmx(1), default=9999999)
      print *,' status 1=',ret
      call getkeyvalue("GWinput","emax_sigm", ebmx,2, default=(/1d10,1d10/),status=ret)
      print *,' status 2=',ret
      if(ret<0) call getkeyvalue("GWinput","emax_sigm",ebmx(1), default=1d10)
      print *,' status 3=',ret
      return
      end

c---------
      subroutine readin5(i0,i1,i2)
c! readin i0,i1,i2; these defaults value are 0 0 0 if these are not given.
      character(len=100):: recxxx
      character(len=106):: recxxx2
      read (5,"(a100)",end=1100) recxxx
 1100 continue
      recxxx2 = recxxx//' 0 0 0'
      read(recxxx2,*) i0, i1, i2
      end

      subroutine readin6(i0,i1,i2,i3)
c! readin i0,i1,i2; these defaults value are 0 0 0 if these are not given.
      character(len=100):: recxxx
      character(len=106):: recxxx2
      read (5,"(a100)",end=1100) recxxx
 1100 continue
      recxxx2 = recxxx//' 0 0 0 0'
      read(recxxx2,*) i0, i1, i2,i3
      end

c-------------
      character*8 function xt(num)
      integer(4) :: num
      if(num==0) then
        xt=''
        return
      endif
      xt = char(48+mod(num,10))
      if(num>9)     xt = char(48+mod(num/10,10))//xt
      if(num>99)    xt = char(48+mod(num/100,10))//xt
      if(num>999)   xt = char(48+mod(num/1000,10))//xt
      if(num>9999)  xt = char(48+mod(num/10000,10))//xt
      if(num>99999) xt = char(48+mod(num/100000,10))//xt
      if(num>999999) stop ' xt:can not produce'
      xt='.'//xt
      end
c-------------
      character*20 function xxt(num1,num2)
      integer(4) :: num,num1,num2
      num = num2
      xxt = char(48+mod(num,10))
      if(num>9)     xxt = char(48+mod(num/10,10))//xxt
      if(num>99)    xxt = char(48+mod(num/100,10))//xxt
      if(num>999)   xxt = char(48+mod(num/1000,10))//xxt
      if(num>9999)  xxt = char(48+mod(num/10000,10))//xxt
      if(num>99999) xxt = char(48+mod(num/100000,10))//xxt
      if(num>999999) stop ' xxt:can not produce'
      xxt ='to'//xxt
      num = num1
      xxt = char(48+mod(num,10))//xxt
      if(num>9)     xxt = char(48+mod(num/10,10))//xxt
      if(num>99)    xxt = char(48+mod(num/100,10))//xxt
      if(num>999)   xxt = char(48+mod(num/1000,10))//xxt
      if(num>9999)  xxt = char(48+mod(num/10000,10))//xxt
      if(num>99999) xxt = char(48+mod(num/100000,10))//xxt
      if(num>999999) stop ' xxt:can not produce'
      end
c
      character(5) function charnum5(num)
      integer(4) ::num
      charnum5 = char(48+mod(num/10000,10))//
     &           char(48+mod(num/1000,10))//
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
c
      character(3) function charnum3(num)
      integer(4) ::num
      charnum3 = 
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
c
      character(4) function charnum4(num)
      integer(4) ::num
      charnum4 = 
     &           char(48+mod(num/1000,10))//
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
c
      character*20 function xts(num1,num2)
      integer(4) :: num,num1,num2
c      character*(*) add
      num = num2
      xts = char(48+mod(num,10))
      if(num>9)     xts = char(48+mod(num/10,10))//xts
      if(num>99)    xts = char(48+mod(num/100,10))//xts
      if(num>999)   xts = char(48+mod(num/1000,10))//xts
      if(num>9999)  xts = char(48+mod(num/10000,10))//xts
      if(num>99999) xts = char(48+mod(num/100000,10))//xts
      if(num>999999) stop ' xts:can not produce'
      xts ='.L'//xts
      num = num1
      xts = char(48+mod(num,10))//xts
      if(num>9)     xts = char(48+mod(num/10,10))//xts
      if(num>99)    xts = char(48+mod(num/100,10))//xts
      if(num>999)   xts = char(48+mod(num/1000,10))//xts
      if(num>9999)  xts = char(48+mod(num/10000,10))//xts
      if(num>99999) xts = char(48+mod(num/100000,10))//xts
      if(num>999999) stop ' xts:can not produce'
      end
