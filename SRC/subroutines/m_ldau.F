      module m_ldau
      use m_ftox
      public:: m_ldau_init,m_ldau_vorbset
      complex(8),allocatable,protected,public::  vorb(:,:,:,:)
      real(8),protected,public:: eorb=0d0
!!
      private
      complex(8),allocatable,private::  dmato(:,:,:,:)
      logical,private:: init=.true.
      contains
!! LDA+U initialization
      subroutine m_ldau_init()
      use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu
c      use m_chkdmu,only: sudmtu
      integer:: i,idmatu
      logical:: mmtargetx
      call tcn('m_ldau_init')
      idmatu=nsp*nlibu*(lmaxu*2+1)**2
      if(init) then
      allocate( dmato(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
      allocate( vorb (-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
      init=.false.
      endif
      call sudmtu(dmato, vorb)  !read dmatu from dmatu.ext and generate vorbdmat from dmatu.ext
      inquire(file='mmtarget.aftest',exist=mmtargetx)
      if(mmtargetx) call vorbmodifyaftest()
      call tcx('m_ldau_init')
      end
      
!! vorb LDA+U mixing
      subroutine m_ldau_vorbset(eks,dmatu)
c      use m_chkdmu,only: Chkdmu
      use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu,nbas,stdo
      use m_MPItk,only: master_mpi
      complex(8):: dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      complex(8):: vorbav(-lmaxu:lmaxu,-lmaxu:lmaxu)
      real(8):: eks
      integer,parameter:: nx=1000
      logical:: mmtargetx
      real(8):: alpha,mmtarget, uhhist(nx),uhx,sss(1),fac,mmsite(nbas),mmhist(nx)
      integer:: nit,ibas,ifx,ncount=0,iblu,i
      call tcn('m_ldau_vorbset')
      call chkdmu(eks,dmatu, dmato, vorb)
      dmato=dmatu !dmato is kept in this module as previous dmat
      inquire(file='mmtarget.aftest',exist=mmtargetx)
      if(mmtargetx) call vorbmodifyaftest()
      call tcx('m_ldau_vorbset')
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
      use m_amix,only: amix
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
      integer im,imix,jmix,onorm,okpvt,oa !,amix
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
c      allocate(norm(mxsav**2),kpvt(mxsav))
      imix = amix(nda,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax,
     .  a,tj,rms2) !norm,kpvt,
c      deallocate(norm, kpvt)
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
      subroutine chkdmu(eks, dmatu,dmatuo,vorb)
      use m_struc_def,only: s_site,s_spec
      use m_lmfinit,only: stdl,nbas,nsp,nlibu,lmaxu,ssite=>v_ssite,sspec=>v_sspec,lldau,
     &     tolu=>mix_tolu,umix=>mix_umix,stdo
      use m_MPItk,only: master_mpi
      use m_mksym,only: g=>rv_a_osymgr,istab=>iv_a_oistab, ng =>lat_nsgrp
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
      implicit none
      intent(in)::       eks,       dmatuo
      intent(out)::                        vorb
! LDA+U potential vorb and total energy eorb.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu: number of U blocks
Ci   lmaxu :dimensioning parameter for U matrix
Ci   dmatu : dmatu produced in current iteration
Ci         : dmatu is passed in real harmonics
Ci   dmatuo: dmatu produced in prior iteration
Ci         : dmatuo is passed in real harmonics
Cixxx   idvsh=0 dmatu, dmatuo, vorb input/output in real harmonics
Ci   umix  :linear mixing parameter for density matrix
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   istab :site istab(i,ig) is transformed into site i by grp op ig
Co  vorb  :orbital dependent potential matrices
Co  eorb  : U contribution to LDA+U total energy
C ----------------------------------------------------------------------
      integer:: ierr
      include "mpif.h"
      integer:: idvsh=0,i_copy_size
      real(8):: eks, eorbxxx
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex dmatuo(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      integer l,idu(4),lmxa,ib,is,iblu,igetss,idmat,ivsiz,ifile_handle
      integer iprint,ipl,havesh
      double precision ddmat,uh(4),jh(4),eorbi,eterms(20),ddot,xx
      logical:: fexist,mmtargetx,eee
      real(8),allocatable:: uhall(:,:)
      real(8):: mmsite(nbas),uhxx,mmhist(10000),uhhist(10000),mmtarget,uhdiff
      real(8),save::uhxnew,uhx,alpha,alphax
      integer:: nn,ibas,ifx,key,i,nit
      integer,save::ncount=0
      complex(8):: dmatuav(-lmaxu:lmaxu,-lmaxu:lmaxu)
      real(8):: sss(1),sigin
      real(8):: fac
      if (nlibu .eq. 0) return
      havesh =0
      idvsh  =0  ! We assume real harmonics for i/o
      ipl = 1
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
C --- Symmetrize output dmatu (req. real harmonics); compare diff ---
      if(iprint()>=60) call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .' Unsymmetrized output dmats',dmatu)
      call symdmu(nlibu,dmatu, nbas,nsp, lmaxu, sspec, ssite, ng, g, istab, lldau, xx)
      ddmat = dsqrt(ddot(2*ivsiz,dmatu-dmatuo,1,dmatu-dmatuo,1)/(2*ivsiz))
      if(master_mpi)write(stdo,ftox)' chkdmu: LDA+U. RMSdiff of dmat from symmetrization =',ftod(xx,2)
C --- Compute U contribution to total energy; make vorb ---
      call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !mode=1, dmatu is from rh to sh 
      havesh = 1
      eorb = 0
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)
          lmxa=sspec(is)%lmxa
          idu=sspec(is)%idu
          uh=sspec(is)%uh
          jh=sspec(is)%jh
          do  l = 0, min(lmxa,3)
            if (idu(l+1) .ne. 0) then
              iblu = iblu+1
              eorbi = 999
              call ldau(100+idu(l+1),l,iblu,uh(l+1),jh(l+1),dmatu,nsp,lmaxu,vorb,eorbi)
              eorb = eorb + eorbi
            endif
          enddo
        endif
      enddo
C --- LDA total energy terms ---
      if(master_mpi)write(stdo,ftox)' eks =',ftof(eks),'e[U]=',ftof(eorb),'Etot(LDA+U)=',ftof(eks+eorb)
      if(master_mpi)write(stdl,ftox)' ldau EHK',ftof(eks),'U',ftof(eorb),'ELDA+U',ftof(eks+eorb)
      
C --- Restore dmatu, vorb to real harmonics 
      call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !-1, from sh to rh idvsh=0
      havesh = 0
      if(master_mpi)write(stdo,ftox)' LDA+U update density matrix ...'
      if(master_mpi)write(stdo,ftox)' RMS diff in densmat',ftod(ddmat) 
      dmatu = umix*dmatu+(1d0-umix)*dmatuo ! new*umix + old*(1-umix)
      call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !from rh to sh
      havesh = 1
      
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccc   
c$$$!! experimental block to keep magnetic moment for AF.
c$$$       inquire(file='mmtarget.aftest',exist=mmtargetx)
c$$$       if(mmtargetx) then
c$$$          if(master_mpi) then
c$$$             open(newunit=ifx,file='mmtarget.aftest')
c$$$             alpha  = .1d0
c$$$             read(ifx,*) mmtarget !,alpha
c$$$             close(ifx)
c$$$             open(newunit=ifx,file='uhval.aftest')
c$$$             nit=0
c$$$             mmsite=0d0
c$$$             uhx=0d0
c$$$             do
c$$$                read(ifx,*,end=1112,err=1112) uhx
c$$$                read(ifx,*,end=1112,err=1112) (mmsite(ibas),ibas=1,nbas)
c$$$                nit=nit+1
c$$$                mmhist(nit)=(mmsite(1)-mmsite(2))/2d0
c$$$                uhhist(nit)=uhx
c$$$             enddo
c$$$ 1112        continue
c$$$             close(ifx)
c$$$             write(6,"('uhval: ',i5,f10.6,2x,12f10.6)")nit,uhx,(mmsite(ibas),ibas=1,nbas)
c$$$!     ! Generate new uhx based on the  history of uhx mmsites for given mm
c$$$!     ! test uh          
c$$$             uhx= uhx + alpha*(mmhist(nit)-mmtarget)**2 - 2d0*(mmhist(nit)-mmtarget)
c$$$             write(6,"('mmhist0: UH mm',i5, 3d13.4)" ) nit,uhx,mmhist(nit)
c$$$             sss(1)=uhx
c$$$             call mixuh(sss)
c$$$             uhx=sss(1)
c$$$             write(6,"('mmhist:  UH mm ',i5, 3d13.4)") nit,uhx,mmhist(nit)
c$$$             if(ncount==0) then
c$$$                open(newunit=ifx,file="mixuh.aftest")
c$$$                close(ifx,status='delete')
c$$$                ncount=1
c$$$             endif
c$$$             open(newunit=ifx,file='uhval.aftest',position='append')
c$$$             write(ifx,"(d23.15,' !uhx')") uhx !, uhxnew
c$$$             close(ifx)
c$$$          endif
c$$$!! broadcast uhx       
c$$$          call MPI_Bcast(uhx, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
c$$$! Averaged dmatu up dn.
c$$$          dmatuav=0d0
c$$$          do i=-lmaxu,lmaxu
c$$$             dmatuav(i,i) =  1d0
c$$$          enddo
c$$$          do iblu =1,2
c$$$             fac=1d0
c$$$             if(iblu==2) fac=-1d0
c$$$c             dmatuav(:,:) =  .5d0*dmatu(:,:,1,iblu) +.5d0*dmatu(:,:,2,iblu)
c$$$             dmatu(:,:,1,iblu) =  fac*dmatuav(:,:)
c$$$             dmatu(:,:,2,iblu) = -fac*dmatuav(:,:)
c$$$          enddo
c$$$
c$$$       endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      vorb=1d99
      iblu = 0
      do  ib = 1, nbas
         if (lldau(ib) .ne. 0) then
            is = int(ssite(ib)%spec)
            lmxa=sspec(is)%lmxa
            idu=sspec(is)%idu
            uh=sspec(is)%uh
            jh=sspec(is)%jh
            do  l = 0, min(lmxa,3)
               if (idu(l+1) .ne. 0) then
                  iblu = iblu+1
                  call pshpr(iprint()-20)
c$$$  if(.not.mmtargetx) uhx = uh(l+1)
                  uhx = uh(l+1)
                  call ldau(idu(l+1),l,iblu,uhx,jh(l+1),dmatu,nsp, !uhxx,jh(l+1),dmatu,nsp, !
     .                 lmaxu,vorb,eorbxxx) !eorbxxx is dummy?
                  call poppr
               endif
            enddo
         endif
      enddo
C  ... Symmetrize vorb to check (symdmu requires real harmonics)
      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau) !vorb from sh to rh
      call symdmu(nlibu,vorb,nbas , nsp , lmaxu , sspec, ssite , ng , g , istab , lldau , xx )
      call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau) !vorb from rh to sh
      
!!=>  At this point, dmatu and vorb are in spherical harmonics
      if(Iprint()>20) write(6,ftox)'  RMS change in vorb from symmetrization =',ftod(xx)
      if(xx>.0001d0.and.iprint()>30) write(6,'(a)')' (warning) RMS change unexpectely large'
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite, ' Mixed dmats',dmatu)
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite, ' New vorb',vorb)
      if (master_mpi) then
         idmat = ifile_handle()
         open(idmat,file='dmats.'//trim(sname)) !havesh mush be 1
         call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' mixed dmats',dmatu)
         close(idmat)
      endif
!! write in real harmonics      
      call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !from sh to rh
      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)  !from sh to rh
      havesh=0                  !I recovered this 2022May8
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' Mixed dmats',dmatu)
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' New vorb',vorb)
      end subroutine chkdmu
!ssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine sudmtu(dmatu,vorb)
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
      use m_lmfinit,only: nbas,nsp,nlibu,lmaxu,lldau,ssite=>v_ssite,sspec=>v_sspec,stdo
      use m_mksym,only: g=>rv_a_osymgr,istab=>iv_a_oistab, ng =>lat_nsgrp
C- Initialize site density matrix and vorb  for LDA+U
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu : nlibu total number of U blocks
Ci   lmaxu :dimensioning parameter for U matrix
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:   *
Ci     Passed to: symdmu rotycs
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idu uh jh
Ci     Stored:   *
Ci     Passed to: symdmu rotycs
Ci   idvsh :0 dmatu and vorb returned in real harmonics
Ci         :1 dmatu and vorb returned in spherical harmonics
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   istab :site istab(i,ig) is transformed into site i by grp op ig
Co Outputs
Co   dmatu :density matrix for LDA+U orbitals
Co         :in real spherical harmonics basis
Co   vorb  :orbital dependent potential matrices
Co         :in real spherical harmonics basis
Cl Local variables
Cl   eorb  : U contribution to LDA+U total energy
Cr Remarks
Cr   Reads in diagonal occupation numbers from file occnum.ext or dmatu
Cr   given in order of m=-l,l, isp=1,2, and constructs initial vorb
Cu Updates
Cu   12 Nov 07 Generate dmatu and vorb in either real or spher. harmonics
Cu   07 May 07 Bug fix MPI mode, when reading occnum instead of dmats
Cu   31 Jan 06 Setup and printouts in spherical harmonics
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   27 Apr 05 Lambrecht first created
C-------------------------------------------------------------------
      implicit none
      integer:: idvsh=0,i_copy_size !nbas,nsp,nlibu,lmaxu,
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      logical rdstrn,parstr,mmtargetx,eee
      integer i,isp,ib,l,lmxa,m,m2,foccn,havesh,ivsiz,ipr,ifx
      double precision nocc(-3:3,2),iv(7)
      integer:: idu(4),igetss,is,idmat,fxst,iblu,nlm,nlmu,a2vec,nn,m1
      double precision uh(4),jh(4),eorb,xx !tmp(2,7,7)
      complex(8):: tmp(7,7),img=(0d0,1d0)
      real(8):: tempr(7,7),tempi(7,7)
      character str*80,spid*8,aaa*24
      complex(8) :: dmwk_zv(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      real(8):: uhx,uhxx
C ... MPI
      include "mpif.h"
      integer procid,master,mpipid,ierr
      logical mlog,occe,dexist,readtemp
      real(8)::sss
      character(128):: bbb
      call rxx(nsp.ne.2,'LDA+U must be spin-polarized!')
      procid = mpipid(1)
      master = 0
      call getpr(ipr)
!! When LDAU is dummy (usually just in order to print our dmats file).
      sss=0d0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)
          sss = sss + sum(abs(sspec(is)%uh))+sum(abs(sspec(is)%jh))
        endif
      enddo
      if(sss<1d-6) then
        dmatu = 0d0
        havesh = 1
        goto 1185
      endif
! Read in dmatu if file  dmats.ext  exists ---
      if(procid /= master) goto 1185
      inquire(file='dmats.'//trim(sname),exist=dexist)
      inquire(file='occnum.'//trim(sname),exist=occe) !if no dmats, try occnum. 
      if(dexist) then
         open(newunit=idmat,file='dmats.'//trim(sname)) 
 825     continue
         if (rdstrn(idmat,str,len(str),.false.)) then
            if (str(1:1) .eq. '#') goto 825
            i = 0
            if (parstr(str,'sharm ',len(str)-6,5,' ',i,m)) then
               havesh = 1
            else
               havesh = 0
            endif
         endif
         if(havesh ==1) bbb='spherical harmonics'
         if(havesh ==0) bbb='real harmonics'
         write(stdo,*)' sudmtu: reading density matrix from file dmats in '//trim(bbb)
         rewind idmat
         iblu = 0
         do  ib = 1, nbas
            if (lldau(ib) .ne. 0) then
               is = int(ssite(ib)%spec)
               lmxa=sspec(is)%lmxa
               idu=sspec(is)%idu
               do l = 0, min(lmxa,3)
                  if (idu(l+1) .ne. 0) then
                     iblu = iblu+1
                     nlm = 2*l+1
                     nlmu = 2*lmaxu+1
                     do  isp = 1, 2
                        readtemp=.false.
                        read(idmat,*,end=9888,err=9888)                     
                        do m1=1,nlm
                           read(idmat,*,end=9888,err=9888) (tempr(m1,m2),m2=1,nlm)
                        enddo
                        read(idmat,*,end=9888,err=9888)
                        do m1=1,nlm
                           read(idmat,*,end=9888,err=9888) (tempi(m1,m2),m2=1,nlm)
                        enddo
                        readtemp=.true.
 9888                   continue
c     readtemp = rdm(idmat,40,2*nlm**2,' ',tmp(1:nlm,1:nlm),nlm,nlm) .ne. 2
                        if(.not.readtemp) call rxi('sudmtu failed to read dmats for site',ib)
                        dmatu(-l:l,-l:l,isp,iblu)=tempr(1:nlm,1:nlm)+img*tempi(1:nlm,1:nlm)
                     enddo
                  endif
               enddo
            endif
         enddo
      elseif(occe) then         !if no occunum.* initial dmatu=0 
         write(stdo,*) ' sudmtu:  initial (diagonal) density-matrix from occ numbers'
         open(newunit=foccn,file='occnum.'//trim(sname))
         havesh = 1
 12      if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
         if (str(1:1) .eq. '#') goto 12
         if (str(1:1) .eq. '%') then
            i = 0
            if (parstr(str,'real ',len(str)-5,4,' ',i,m)) havesh = 0
         else
            backspace foccn
         endif
         iblu = 0
         do  ib = 1, nbas
            if (lldau(ib) .ne. 0) then
               is = int(ssite(ib)%spec)
               lmxa=sspec(is)%lmxa
               idu=sspec(is)%idu
               do l = 0,min(lmxa,3)
                  if (idu(l+1) .ne. 0) then
                     iblu = iblu+1
                     do  isp = 1, 2
 11                     continue
                        if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
C     Skip comment lines
                        if (str(1:1) .eq. '#') goto 11
                        i = 0
                        m = a2vec(str,len(str),i,4,', ',2,3,2*l+1,iv,nocc(-l,isp))
                        if (m .lt. 0) goto 99
                     enddo
                     do isp=1,2
                        write(stdo,ftox)' occ num: site',ib,'l',l,'isp',isp,' ',nocc(-l:l,isp)
                     enddo   
                     do isp = 1, 2
                        do m = -l, l
                           do m2 = -l, l
                              dmatu(m,m2,isp,iblu) = dcmplx(0d0,0d0)
                           enddo
                           dmatu(m,m,isp,iblu) = dcmplx(nocc(m,isp),0d0)
                        enddo
                     enddo
                  endif
               enddo
            endif
         enddo
         close(foccn)
      else
         dmatu=0d0
      endif
 1185 continue

C ... Initial printout
      call praldm(0,51,51,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' dmats read from disk',dmatu)
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call mpibc1(dmatu,2*ivsiz,4,mlog,'sudmtu','dmatu')
      call mpibc1(havesh,1,2,.false.,' ',' ')
C ... Density matrix in real or spherical harmonics (fixed by idvsh)
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = idvsh
      endif
C ... Symmetrize dmatu (symdmu requires real harmonics)
      dmwk_zv=dmatu
      if (havesh .eq. 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      call symdmu(nlibu,dmatu , nbas , nsp , lmaxu , sspec, ssite , ng , g , istab , lldau , xx )
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs(2*idvsh-1,dmwk_zv,nbas,nsp,lmaxu,sspec,ssite, lldau )
        havesh = idvsh
      endif
      if (ng .ne. 0) then
        write(stdo,ftox)' sudmtu:  RMS change in dmats'//
     .  ' from symmetrization',ftof(xx)
        if (xx .gt. .01d0) write(stdo,*)'(warning) RMS change unexpectely large'
        call daxpy ( ivsiz * 2 , - 1d0 , dmatu , 1 , dmwk_zv , 1 )
        if(ipr>=60) write(stdo,*)' change in dmat wrought by symmetrization'
        call praldm ( 0 , 60 , 60 , 0 , nbas , nsp , lmaxu , lldau , 
     .  sspec , ssite , ' ' , dmwk_zv )
      endif
C     Print dmats in specified harmonics
      dmwk_zv=dmatu
      if (havesh .ne. idvsh) then
        call rotycs ( 2 * idvsh - 1 , dmwk_zv , nbas , nsp , lmaxu, sspec , ssite , lldau )
      endif
      write(stdo,*)
      call praldm(0,30,30,idvsh,nbas,nsp,lmaxu,lldau,sspec , ssite , ' Symmetrized dmats' , dmwk_zv )
C     Print dmats in complementary harmonics
      i = 1-idvsh
      call rotycs(2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , sspec , ssite , lldau )
      write(stdo,*)' '
      call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,sspec,ssite, ' Symmetrized dmats' , dmwk_zv )
C ... Make Vorb (ldau requires spherical harmonics)
      if (havesh .ne. 1) then
        call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 1
      endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c$$$!! experimental block to keep magnetic moment for AF.
c$$$       inquire(file='mmtarget.aftest',exist=mmtargetx)
c$$$       inquire(file='mmagfield.aftest',exist=eee)
c$$$       if(mmtargetx.and. (procid==master)) then
c$$$          uhx=0d0
c$$$          open(newunit=ifx,file='mmagfield.aftest')
c$$$          do
c$$$             read(ifx,*,end=1112,err=1112) uhxx,aaa
c$$$             uhx=uhxx
c$$$             if(trim(aaa)=='!Magfield') uhx=uhxx
c$$$          enddo
c$$$ 1112     continue
c$$$          close(ifx)
c$$$          write(stdo,"('sudmtu: mmtarget mode. Readin Magfield from mmagfield.aftest=',f10.6)")uhx
c$$$          open(newunit=ifx,file='mmagfield.aftest')
c$$$          write(ifx,"(d23.15,1x,'!Magfield is read from previous mmagfield.aftest')") uhx
c$$$          close(ifx)
c$$$       endif
c$$$       call mpibc1_real(uhx,1,'sudmtu_uhx')
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      iblu = 0
      do  20  ib = 1, nbas
        if (lldau(ib) .eq. 0) goto 20
        is = int(ssite(ib)%spec)
        lmxa=sspec(is)%lmxa
        idu=sspec(is)%idu
        uh=sspec(is)%uh
        jh=sspec(is)%jh
        spid=sspec(is)%name
        i = min(lmxa,3)
        write(stdo,ftox)'Species '//spid//'mode',idu(1:i+1),'U',ftof(uh(1:i+1),2),'J',ftof(jh(1:i+1),2)
        do  22  l = 0, i
          if (idu(l+1) .ne. 0) then
            iblu = iblu+1
c            if(.not.mmtargetx) uhx=uh(l+1)
            uhx=uh(l+1)
            call ldau(idu(l+1),l,iblu,uhx,jh(l+1),dmatu,nsp,lmaxu,vorb,eorb)
          endif
   22   continue
   20 continue
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' Unsymmetrized vorb',vorb)
!==>     At this point, dmatu and vorb are in spherical harmonics
C ... Symmetrize vorb to check (symdmu requires real harmonics)
      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
      call symdmu (nlibu, vorb, nbas , nsp , lmaxu , sspec , ssite , ng , g , istab , lldau , xx )
C     EITHER: vorb =>  spherical harmonics OR dmatu => real harmonics
      if (idvsh .eq. 1) then
        call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 1
      else
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      if (ng .ne. 0) then
        write(stdo,ftox)' sudmtu:  RMS change in vorb from symmetrization = ',ftof(xx)
        if (xx .gt. .01d0) write(stdo,*)'          (warning) RMS change unexpectely large'
      endif
C     Print vorb in specified harmonics
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' Symmetrized vorb',vorb)
      i = 1-idvsh
      dmwk_zv=vorb
      call rotycs( 2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , sspec  , ssite , lldau )
      write(stdo,*) !vorb in complementary harmonics
      call praldm(0,30,30, i , nbas , nsp , lmaxu , lldau , sspec , ssite , ' Vorb' , dmwk_zv )
      eorb = 0d0
C --- Error exit ---
      return
   99 continue
      write(str,"('bad occnum file, site l= ',i5,i5)")ib,l
      call rx(str)
      end subroutine sudmtu
      
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine rotycs(mode,a,nbas,nsp,lmaxu,sspec,ssite,lldau)
      use m_struc_def  !Cgetarg
C- Rotate matrix a from real to spherical harmonics
C for LDA+U objects densmat and vorb
C-------------------------------------
Ci mode =1 from real to spherical
Ci      -1 from spherical to real
Ci a matrix to be transformed a(m,m,isp,iblu)  could be vorb or densmat
Ci nbas : number of sites
Ci nsp  : number of spins
Ci lmaxu: lmax for U
Ci sspec: species info
Ci ssite: sites info
Ci lldau  :lldau(ib)=0 => no U on this site otherwise
Ci        :U on site ib with dmat in dmats(*,lldau(ib))
Co a rotated in place
Cr Remarks
Cr order of cubic harmonics ls (l-1)s,ms...1s 0 1c mc... (l-1)c lc
Cr order of spherical harmonics -l:l
Cr Yl-m=(Ylmc-iYlms)/sqrt(2)  Ylm=(-1)**m*conjg(Yl-m)
Cu Updates
Cu   18 Jan 06 A. Chantis changed rotation matrices in accordance with
Cu             the definition of real harmonics used in the rest of
Cu             the code (Hund's rules satisfied as indicated by orb. moment)
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   30 Apr 05 Lambrecht first created
C----------------------------------------------------------------
      implicit none
      integer nbas,lldau(nbas),mode,lmaxu,nsp,i_copy_size
      complex(8),target:: a(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)
      integer:: ib,m,l,is,igetss,i,j,k,ll,isp,iblu
      complex(8):: b(2,2),c(2,2),add,bb(2,2)
      complex(8),parameter:: s2= 1/dsqrt(2d0), img=(0d0,1d0)
      complex(8),pointer:: rot(:,:)
      complex(8),save,target::rott(2,2,5,-1:1)
      logical,save:: init=.true.
      if(mode/=1.and.mode/=-1) call rx('ROTYCS: mode must be 1 or -1')
      if(init)then
         do m=1,5
            rott(1,1:2,m,-1) = [s2,          s2*(-1d0)**m] !mode -1: spherical to cubic basis
            rott(2,1:2,m,-1) = [img*s2, -img*s2*(-1d0)**m]
            rott(:,:,  m, 1) = transpose(dconjg(rott(:,:,m,-1))) !mode 1 : cubic to spherical
         enddo
         init=.false.
      endif
      iblu = 0
      do  ib = 1, nbas
         if (lldau(ib)==0) cycle
         is  = ssite(ib)%spec
         do  l = 0, min(sspec(is)%lmxa,3)
            if (sspec(is)%idu(l+1) ==0) cycle
            iblu = iblu+1
            do  isp = 1, 2
               do  m = 1, l
                  bb(1,:) = [a( m,m,isp,iblu), a( m,-m,isp,iblu)]
                  bb(2,:) = [a(-m,m,isp,iblu), a(-m,-m,isp,iblu)]
                  rot => rott(:,:,m,mode)
                  c = matmul(matmul(rot,bb),transpose(dconjg(rot))) !c=rot*b*rot^+
                  a(m,m,isp,iblu)   = c(1,1)
                  a(m,-m,isp,iblu)  = c(1,2)
                  a(-m,m,isp,iblu)  = c(2,1)
                  a(-m,-m,isp,iblu) = c(2,2)
               enddo
            enddo
         enddo
      enddo
      end subroutine rotycs
!---------------------------------------      
      end module m_ldau
