      module m_ldau
      complex(8),allocatable,protected::  vorb(:,:,:,:)

      complex(8),allocatable,private::  dmato(:,:,:,:)
      logical,private:: init=.true.
      contains
!! LDA+U initialization
      subroutine m_ldau_init()
      use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu
      integer:: i,idmatu
      logical:: mmtargetx
      idmatu=nsp*nlibu*(lmaxu*2+1)**2
      if(init) then
      allocate( dmato(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
      allocate( vorb (-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
      init=.false.
      endif
      call sudmtu(dmato, vorb)  !read dmatu from dmatu.ext and generate vorbdmat from dmatu.ext
      inquire(file='mmtarget.aftest',exist=mmtargetx)
      if(mmtargetx) call vorbmodifyaftest()
      end
      
!! vorb LDA+U mixing
      subroutine m_ldau_vorbset(eks,dmatu)
      use m_chkdmu,only: Chkdmu
      use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu,nbas,stdo
      use m_MPItk,only: master_mpi
      complex(8):: dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      complex(8):: vorbav(-lmaxu:lmaxu,-lmaxu:lmaxu)
      real(8):: eks
      integer,parameter:: nx=1000
      logical:: mmtargetx
      real(8):: alpha,mmtarget, uhhist(nx),uhx,sss(1),fac,mmsite(nbas),mmhist(nx)
      integer:: nit,ibas,ifx,ncount=0,iblu,i
      call chkdmu(eks,dmatu, dmato, vorb)
      dmato=dmatu
      inquire(file='mmtarget.aftest',exist=mmtargetx)
      if(mmtargetx) call vorbmodifyaftest()
      end
      
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine vorbmodifyaftest()
      use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu,nbas,stdo
      use m_MPItk,only: master_mpi
      complex(8):: vorbav(-lmaxu:lmaxu,-lmaxu:lmaxu)
      integer,parameter:: nx=1000
      real(8):: alpha,mmtarget, uhhist(nx),uhx,sss(1),fac,mmsite(nbas),mmhist(0:nx)
      integer:: nit,ibas,ifx,ncount=0,iblu,i
      logical:: init=.true.
      character(3):: stat
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
!! Experimental block to add magnetic field to vorb to keep given size of magnetic moment
!! Need fixing. This is only for ibas=1 and ibas=2 are antiferro pair.
      if(master_mpi) then
             open(newunit=ifx,file='mmtarget.aftest')
             alpha  = .1d0
             read(ifx,*) mmtarget !,alpha
             close(ifx)
             open(newunit=ifx,file='mmagfield.aftest')
             nit=0
             mmsite=0d0
             mmhist=mmtarget !this gives we use uhx when mmsite are not read
             uhx=0d0
c             do
                read(ifx,*,end=1112,err=1112) uhx 
                read(ifx,*,end=1112,err=1112) (mmsite(ibas),ibas=1,nbas)
                nit=nit+1
                mmhist(nit)=(mmsite(1)-mmsite(2))/2d0
                uhhist(nit)=uhx
c             enddo
 1112        continue
             close(ifx)
             write(stdo,"('mmaftest: ',i5,f10.6,2x,12f10.6)")nit,uhx,(mmsite(ibas),ibas=1,nbas)
!     ! Generate new uhx based on the  history of uhx mmsites for given mm
!     ! test uh          
             if(nit>0) uhx= uhx + alpha*(mmhist(nit)-mmtarget)**2 - 2d0*(mmhist(nit)-mmtarget)
             write(stdo,"('mmhist0: MagF MMom ',i5, 2d13.4)" ) nit,uhx,mmhist(nit)
             sss(1)=uhx
             if(ncount==0) then
                open(newunit=ifx,file="mixmag.aftest")
                close(ifx,status='delete')
                ncount=1
             endif
             call mixmag(sss)
             uhx=sss(1)
             write(stdo,"('mmhist:  MagF MMom ',i5, 2d13.4)") nit,uhx,mmhist(nit)
c             stat='old'
c             if(init) stat='new'
c             open(newunit=ifx,file='mmagfield.aftest',position='append',status='new') !stat)
             open(newunit=ifx,file='mmagfield.aftest',status='unknown') !stat)
             write(ifx,"(d23.15,' !Magfield')") uhx 
             close(ifx)
             if(init) init=.false.
          endif
!! broadcast uhx       
c          call MPI_Bcast(uhx, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
          call mpibc1_real(uhx,1,'m_ldau_magfield') 
          vorbav=0d0
          do i=-lmaxu,lmaxu
             vorbav(i,i) =  1d0
          enddo
          if(nlibu==6) then
             do iblu =1,6
                fac=1d0
                if(iblu>3) fac=-1d0
                vorb(:,:,1,iblu) = vorb(:,:,1,iblu) - uhx*fac*vorbav(:,:) 
             enddo
          else
             do iblu =1,2       !6
                fac=1d0
                if(iblu==2) fac=-1d0 !>3) fac=-1d0
                vorb(:,:,1,iblu) = vorb(:,:,1,iblu) - uhx*fac*vorbav(:,:) 
             enddo
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end



!sssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine mixmag(sss)
c  subroutine pqmixa(nda,nmix,mmix,mxsav,beta,rms2,a,tj)
C- Mixing routine for sigma. Modified from pqmixa in subs/pqmix.f
C- Anderson mixing of a vector
Ci  mmix: number of iterates available to mix
Cio nmix: nmix > 0: number of iter to try and mix
Ci        nmix < 0: use mmix instead of nmix.
Co  nmix: (abs)  number of iter actually mixed.
Co        (sign) <0, intended that caller update nmix for next call.
C  MvS Feb 04 use sigin as input sigma if available (lsigin=T)
C             Add mixnit as parameter
      implicit none
c      logical lsigin
      integer,parameter:: nda=1
      integer:: nmix,mmix
      integer(4),parameter:: mxsav=10
      double precision rms2,tj(mxsav),beta
      integer im,imix,jmix,amix,onorm,okpvt,oa
      integer iprintxx,ifi,nitr,ndaf
      real(8)::sss(nda),sigin(nda)
      real(8):: tjmax
      real(8),allocatable::norm(:),a(:,:,:)
      integer(4),allocatable:: kpvt(:)
      integer(4)::ret
      character*20 fff
      logical fexist
      real(8):: acc
      integer(4):: ido,ifile_handle
      iprintxx = 30
      beta=.3d0
      allocate ( a(nda,0:mxsav+1,2) )
      fff="mixmag.aftest"
      INQUIRE (FILE =fff, EXIST = fexist)
      if(fexist)      write(6,*)'... reading file mixsigma'
      if(.not.fexist) write(6,*)'... No file mixsigma'
      ifi=ifile_handle()
      open(ifi,file=fff,form='unformatted')
      if(fexist) then
        read(ifi,err=903,end=903) nitr,ndaf
        if (ndaf .ne. nda) goto 903
        read(ifi,err=903,end=903) a
        goto 902
      endif
      goto 901
  903 continue
      print 368
  368 format(5x,'(warning) file mismatch ... mixing file not read')
  901 continue
      nitr = 0
 902  continue
      a(:,0,1) = sss   !output
      if(.not.fexist) a(:,0,2) = sss !input
C     if input sigma available, use it instead of file a(:,0,2)
c      if (lsigin) then
c        write(6,*)'... using input sigma read from sigm file'
c        a(:,0,2) = sigin  !input
c      endif
      write(6,*)'sum sss=',sum(abs(sss))
      imix=9
      mmix = min(max(nitr-1,0),imix)
      if (mmix > mxsav) mmix = mxsav
C     this information already printed out by amix
C     write(6,*)'mixing parameters for amix are fixed in mixsigma'
C     write(6,*)'   beta       =', beta
C     write(6,*)'   tjmax      =', tjmax
C     write(6,*)'   mmix mxsav =', mmix,mxsav
c     call getkeyvalue("GWinput","mixtj",acc,default=0d0,status=ret)
      acc=0d0
      if(acc/=0d0) then
        write(6,*)' readin mixtj from GWinput: mixtj=',acc
        tjmax=abs(acc)+1d-3
        if(mmix==1) then
          tj(1)=1d0
        else
          tj(1)= acc
          tj(2)= 1-acc
          mmix=2
        endif
        ido=2
      else
        tjmax=5d0
        ido=0
      endif
      allocate(norm(mxsav**2),kpvt(mxsav))
      imix = amix(nda,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax,
     .  norm,kpvt,a,tj,rms2)
      deallocate(norm, kpvt)
C ... Restore PQ array, updating new x
c      call dpscop(a,w(oa),nda,1+nda*(mxsav+2),1+nda*(mxsav+2),1d0)
c      call dcopy(nda*(mxsav+2)*2,w(oa),1,a,1)
c ...
      sss = a(:,0,2)
      rewind(ifi)
      write(ifi) nitr+1,nda
      write(ifi) a
      close(ifi)
      end
!sssssssssssssssssssssssssssssssssssssssssssssssssss
      end module
