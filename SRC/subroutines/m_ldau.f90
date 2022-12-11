module m_ldau
  use m_mpitk,only:master_mpi
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
    !      use m_chkdmu,only: sudmtu
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
  end subroutine m_ldau_init

  !! vorb LDA+U mixing
  subroutine m_ldau_vorbset(eks,dmatu)
    !      use m_chkdmu,only: Chkdmu
    use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu,nbas,stdo
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
  end subroutine m_ldau_vorbset

  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine vorbmodifyaftest()
    use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu,nbas,stdo
    complex(8):: vorbav(-lmaxu:lmaxu,-lmaxu:lmaxu)
    integer,parameter:: nx=1000
    real(8):: alpha,mmtarget, uhhist(nx),uhx,sss(1),fac,mmsite(nbas),mmhist(0:nx)
    integer:: nit,ibas,ifx,ncount=0,iblu,i
    logical:: init=.true.
    character(3):: stat
    ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
       !             do
       read(ifx,*,end=1112,err=1112) uhx
       read(ifx,*,end=1112,err=1112) (mmsite(ibas),ibas=1,nbas)
       nit=nit+1
       mmhist(nit)=(mmsite(1)-mmsite(2))/2d0
       uhhist(nit)=uhx
       !             enddo
1112   continue
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
       !             stat='old'
       !             if(init) stat='new'
       !             open(newunit=ifx,file='mmagfield.aftest',position='append',status='new') !stat)
       open(newunit=ifx,file='mmagfield.aftest',status='unknown') !stat)
       write(ifx,"(d23.15,' !Magfield')") uhx
       close(ifx)
       if(init) init= .FALSE. 
    endif
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
    ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  end subroutine vorbmodifyaftest



  ! ssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine mixmag(sss)
    use m_amix,only: amix
    !  subroutine pqmixa(nda,nmix,mmix,mxsav,beta,rms2,a,tj)
    !- Mixing routine for sigma. Modified from pqmixa in subs/pqmix.f
    !- Anderson mixing of a vector
    !i  mmix: number of iterates available to mix
    ! o nmix: nmix > 0: number of iter to try and mix
    !i        nmix < 0: use mmix instead of nmix.
    !o  nmix: (abs)  number of iter actually mixed.
    !o        (sign) <0, intended that caller update nmix for next call.
    !  MvS Feb 04 use sigin as input sigma if available (lsigin=T)
    !             Add mixnit as parameter
    implicit none
    !      logical lsigin
    integer,parameter:: nda=1
    integer:: nmix,mmix
    integer(4),parameter:: mxsav=10
    double precision :: rms2,tj(mxsav),beta
    integer :: im,imix,jmix,onorm,okpvt,oa !,amix
    integer :: iprintxx,ifi,nitr,ndaf
    real(8)::sss(nda),sigin(nda)
    real(8):: tjmax
    real(8),allocatable::norm(:),a(:,:,:)
    integer(4),allocatable:: kpvt(:)
    integer(4)::ret
    character(20) :: fff
    logical :: fexist
    real(8):: acc
    integer(4):: ido,ifile_handle
    iprintxx = 30
    beta=.3d0
    allocate ( a(nda,0:mxsav+1,2) )
    fff="mixmag.aftest"
    INQUIRE (FILE =fff, EXIST = fexist)
    if(fexist)      write(6,*)'... reading file mixsigma'
    if( .NOT. fexist) write(6,*)'... No file mixsigma'
    ifi=ifile_handle()
    open(ifi,file=fff,form='unformatted')
    if(fexist) then
       read(ifi,err=903,end=903) nitr,ndaf
       if (ndaf /= nda) goto 903
       read(ifi,err=903,end=903) a
       goto 902
    endif
    goto 901
903 continue
    print 368
368 format(5x,'(warning) file mismatch ... mixing file not read')
901 continue
    nitr = 0
902 continue
    a(:,0,1) = sss   !output
    if( .NOT. fexist) a(:,0,2) = sss !input
    !     if input sigma available, use it instead of file a(:,0,2)
    !      if (lsigin) then
    !        write(6,*)'... using input sigma read from sigm file'
    !        a(:,0,2) = sigin  !input
    !      endif
    write(6,*)'sum sss=',sum(abs(sss))
    imix=9
    mmix = min(max(nitr-1,0),imix)
    if (mmix > mxsav) mmix = mxsav
    !     this information already printed out by amix
    !     write(6,*)'mixing parameters for amix are fixed in mixsigma'
    !     write(6,*)'   beta       =', beta
    !     write(6,*)'   tjmax      =', tjmax
    !     write(6,*)'   mmix mxsav =', mmix,mxsav
    !     call getkeyvalue("GWinput","mixtj",acc,default=0d0,status=ret)
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
    !      allocate(norm(mxsav**2),kpvt(mxsav))
    imix = amix(nda,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax, &
         a,tj,rms2) !norm,kpvt,
    !      deallocate(norm, kpvt)
    ! ... Restore PQ array, updating new x
    !      call dpscop(a,w(oa),nda,1+nda*(mxsav+2),1+nda*(mxsav+2),1d0)
    !      call dcopy(nda*(mxsav+2)*2,w(oa),1,a,1)
    ! ...
    sss = a(:,0,2)
    rewind(ifi)
    write(ifi) nitr+1,nda
    write(ifi) a
    close(ifi)
  end subroutine mixmag

  ! ssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine chkdmu(eks, dmatu,dmatuo,vorb)
    use m_lmfinit,only: stdl,nbas,nsp,nlibu,lmaxu,ispec,sspec=>v_sspec,lldau, &
         tolu=>mix_tolu,umix=>mix_umix,stdo,idu,uh,jh,ham_lsig
    use m_MPItk,only: master_mpi
    use m_mksym,only: g=>rv_a_osymgr,istab=>iv_a_oistab, ng =>lat_nsgrp
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    use m_ldauu,only: ldau
    implicit none
    intent(in)::       eks,       dmatuo
    intent(out)::                        vorb
    ! LDA+U potential vorb and total energy eorb.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlibu: number of U blocks
    !i   lmaxu :dimensioning parameter for U matrix
    !i   dmatu : dmatu produced in current iteration
    !i         : dmatu is passed in real harmonics
    !i   dmatuo: dmatu produced in prior iteration
    !i         : dmatuo is passed in real harmonics
    ! xxx   idvsh=0 dmatu, dmatuo, vorb input/output in real harmonics
    !i   umix  :linear mixing parameter for density matrix
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat in dmats(*,lldau(ib))
    !i   ng    :number of group operations
    !i   g     :point group operations
    !i   istab :site istab(i,ig) is transformed into site i by grp op ig
    !o  vorb  :orbital dependent potential matrices
    !o  eorb  : U contribution to LDA+U total energy
    ! ----------------------------------------------------------------------
    integer:: ierr
    include "mpif.h"
    integer:: idvsh=0
    real(8):: eks, eorbxxx
    double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    double complex dmatuo(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    integer :: l,lmxa,ib,is,iblu,igetss,idmat,ivsiz,ifile_handle !,idu(4)
    integer :: iprint,ipl,havesh
    double precision :: ddmat,eorbi,eterms(20),ddot,xx !,uh(4),jh(4)
    logical:: fexist,mmtargetx,eee
    real(8),allocatable:: uhall(:,:)
    real(8):: mmsite(nbas),uhxx,mmhist(10000),uhhist(10000),mmtarget,uhdiff,ddo
    real(8),save::uhxnew,uhx,alpha,alphax
    integer:: nn,ibas,ifx,key,i,nit
    integer,save::ncount=0
    complex(8):: dmatuav(-lmaxu:lmaxu,-lmaxu:lmaxu)
    real(8):: sss(1),sigin
    real(8):: fac,ssss
    if (nlibu == 0) return
    ssss=0d0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib) !ssite(ib)%spec
          ssss = ssss + sum(abs(uh(:,is)))+sum(abs(jh(:,is)))
       endif
    enddo
    havesh =0
    idvsh  =0  ! We assume real harmonics for i/o
    ipl = 1
    ivsiz = nsp*nlibu*(lmaxu*2+1)**2
    ! --- Symmetrize output dmatu (req. real harmonics); compare diff ---
    if(iprint()>=60) call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,' Unsymmetrized output dmats',dmatu)
    call symdmu(nlibu,dmatu, nbas,nsp, lmaxu, ng, g, istab, lldau, xx)
    if(master_mpi)write(stdo,ftox)
    if(master_mpi)write(stdo,ftox)'chkdmu: LDA+U. RMSdiff of dmat from symmetrization =',ftod(xx,2)
    ! --- Compute U contribution to total energy; make vorb ---
    call rotycs(1,dmatu,nbas,nsp,lmaxu,lldau) !mode=1, dmatu is from rh to sh
    havesh = 1
    eorb = 0
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib) 
          lmxa=sspec(is)%lmxa
          do  l = 0, min(lmxa,3)
             if (idu(l+1,is) /= 0) then
                iblu = iblu+1
                eorbi = 999
                call ldau(100+idu(l+1,is),l,iblu,uh(l+1,is),jh(l+1,is),dmatu,nsp,lmaxu,vorb,eorbi)
                eorb = eorb + eorbi
             endif
          enddo
       endif
    enddo
    ! --- LDA total energy terms ---
    if(master_mpi)write(stdo,ftox)'eks =',ftof(eks),'e[U]=',ftof(eorb),'Etot(LDA+U)=',ftof(eks+eorb)
    if(master_mpi)write(stdl,ftox)'ldau EHK',ftof(eks),'U',ftof(eorb),'ELDA+U',ftof(eks+eorb)
    ! --- Restore dmatu, vorb to real harmonics
    call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau) !-1, from sh to rh idvsh=0
    havesh = 0
    if(master_mpi)then
       ddmat = sum(abs(dmatu-dmatuo)**2)**.5
       write(stdo,ftox)'LDA+U update density matrix ... RMS diff in densmat',ftod(ddmat)
    endif   
    ddo = sum(abs(dmatuo)**2) !dmatuo=0d0 if no dmatu.* occnum.*
    if(ddo>1d-10.and.ssss>1d-10) dmatu = umix*dmatu+(1d0-umix)*dmatuo ! new*umix + old*(1-umix)
    call rotycs(1,dmatu,nbas,nsp,lmaxu,lldau) !from rh to sh
    havesh = 1
    vorb=1d99
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib)
          lmxa=sspec(is)%lmxa
          do  l = 0, min(lmxa,3)
             if (idu(l+1,is) /= 0) then
                iblu = iblu+1
                call pshpr(iprint()-20)
                uhx = uh(l+1,is)
                call ldau(idu(l+1,is),l,iblu,uhx,jh(l+1,is),dmatu,nsp, lmaxu,vorb,eorbxxx) 
                call poppr
             endif
          enddo
       endif
    enddo
    !  ... Symmetrize vorb to check (symdmu requires real harmonics)
    call rotycs(-1,vorb,nbas,nsp,lmaxu,lldau) !-1 means that vorb is converted from sh to rh
    call symdmu(nlibu,vorb,nbas , nsp , lmaxu , ng , g , istab , lldau , xx )
    call rotycs(1,vorb,nbas,nsp,lmaxu,lldau) !1 means that vorb is converted from rh to sh
    !!=>  At this point, dmatu and vorb are in spherical harmonics
    if(Iprint()>20) write(stdo,ftox)'RMS change in vorb from symmetrization =',ftod(xx)
    if(xx>.0001d0 .AND. iprint()>30) write(stdo,'(a)')'(warning) RMS change unexpectely large'
!    if(iprint()>0) write(6,ftox)'=== representation in spherical harmonics dmatu ==='
!    call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau, 'Mixed dmats',dmatu)
!    if(iprint()>0) write(6,ftox)'=== representation in spherical harmonics vorb ==='
!    call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau, 'New vorb',vorb)
    if (master_mpi) then
       idmat = ifile_handle()
       open(idmat,file='dmats.'//trim(sname)) !havesh mush be 1
       if(iprint()>0) write(stdo,*)' ... Writing density matrix dmats.*... '
       !write(idmat,ftox)'# === dmats in sph harmonics (read by lmfp-m_ldau_init-sudmtu)==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,' dmats !spherical harmonics',dmatu)
       write(idmat,ftox)'# --- We do not read following data. Just for check ---'
       write(idmat,ftox)'# vorb in sph harmonics (not readin) ==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau, 'New vorb',vorb)
    endif
    call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau) !dmatu is converted from sh to rh
    call rotycs(-1,vorb,nbas,nsp,lmaxu,lldau)  !vorb is converted from sh to rh
    havesh=0                  !I recovered this 2022May8
    if (master_mpi) then ! write in real harmonics
       write(idmat,ftox)'# dmatu in real harmonics (not readin) ==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,'Mixed dmats',dmatu)
       write(idmat,ftox)'# vorb  in real harmonics (not readin) ==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,'New vorb',vorb)
       close(idmat)
    endif
  end subroutine chkdmu
  ! sssssssssssssssssssssssssssssssssssssssssssssss
  subroutine sudmtu(dmatu,vorb)
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    use m_lmfinit,only: nbas,nsp,nlibu,lmaxu,lldau,ispec,sspec=>v_sspec,stdo,slabl,idu,uh,jh
    use m_mksym,only: g=>rv_a_osymgr,istab=>iv_a_oistab, ng =>lat_nsgrp
    !- Initialize site density matrix and vorb  for LDA+U
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlibu : nlibu total number of U blocks
    !i   lmaxu :dimensioning parameter for U matrix
    !i   ssite :struct for site-specific information; see routine usite
    !i     Elts read: spec
    !i     Stored:   *
    !i     Passed to: symdmu rotycs
    !i   sspec :struct for species-specific information; see routine uspec
    !i     Elts read: lmxa idu uh jh
    !i     Stored:   *
    !i     Passed to: symdmu rotycs
    !i   idvsh :0 dmatu and vorb returned in real harmonics
    !i         :1 dmatu and vorb returned in spherical harmonics
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat in dmats(*,lldau(ib))
    !i   ng    :number of group operations
    !i   g     :point group operations
    !i   istab :site istab(i,ig) is transformed into site i by grp op ig
    !o Outputs
    !o   dmatu :density matrix for LDA+U orbitals
    !o         :in real spherical harmonics basis
    !o   vorb  :orbital dependent potential matrices
    !o         :in real spherical harmonics basis
    !l Local variables
    !l   eorb  : U contribution to LDA+U total energy
    !r Remarks
    !r   Reads in diagonal occupation numbers from file occnum.ext or dmatu
    !r   given in order of m=-l,l, isp=1,2, and constructs initial vorb
    !u Updates
    !u   12 Nov 07 Generate dmatu and vorb in either real or spher. harmonics
    !u   07 May 07 Bug fix MPI mode, when reading occnum instead of dmats
    !u   31 Jan 06 Setup and printouts in spherical harmonics
    !u   09 Nov 05 (wrl) Convert dmat to complex form
    !u   27 Apr 05 Lambrecht first created
    !-------------------------------------------------------------------
    use m_ldauu,only: ldau
    implicit none
    integer:: idvsh=0
    complex(8):: dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),&
         Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    logical :: mmtargetx,eee
    integer :: i,isp,ib,l,lmxa,m,m2,foccn,havesh,ivsiz,ipr,ifx
    real(8):: nocc(-3:3,2),iv(7),eorb,xx
    integer:: igetss,is,idmat,fxst,iblu,nlm,nlmu,nn,m1 !idu(4),
    complex(8):: tmp(7,7),img=(0d0,1d0)
    real(8):: tempr(7,7),tempi(7,7)
    character str*80,spid*8,aaa*24
    complex(8) :: dmwk_zv(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    real(8):: uhx,uhxx
    ! ... MPI
    include "mpif.h"
    integer :: procid,master,mpipid,ierr
    logical :: mlog,occe,dexist,readtemp,cmdopt0
    real(8)::sss
    character(128):: bbb
    call rxx(nsp.ne.2,'LDA+U must be spin-polarized!')
    procid = mpipid(1)
    master = 0
    call getpr(ipr)
    !! When LDAU is dummy (usually just in order to print our dmats file).
    ! Read in dmatu if file  dmats.ext  exists ---
    if(procid /= master) goto 1185
    inquire(file='dmats.'//trim(sname),exist=dexist)
    inquire(file='occnum.'//trim(sname),exist=occe) !if no dmats, try occnum.
    if(dexist) then
       open(newunit=idmat,file='dmats.'//trim(sname))
825    continue
       read(idmat,*)str
       if(str(1:1) == '#') goto 825
       if(index(str,' sharm ')/=0) then
          havesh = 1 !spherical(complex) harmonics
       else
          havesh = 0 !real harmonics
       endif
       if(havesh ==1) bbb='spherical harmonics' !complex harmonics
       if(havesh ==0) bbb='real harmonics'      !real harmonics 
       if(master_mpi) write(stdo,*)' sudmtu: reading density matrix from file dmats in '//trim(bbb)
       rewind idmat
       iblu = 0
       do  ib = 1, nbas
          if (lldau(ib) /= 0) then
             is = ispec(ib) 
             lmxa=sspec(is)%lmxa
             do l = 0, min(lmxa,3)
                if (idu(l+1,is) /= 0) then
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
9888                  continue
                      if( .NOT. readtemp) call rxi('sudmtu failed to read dmats for site',ib)
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
12     continue
       read(foccn,*) str
       if (str(1:1) == '#') goto 12
       if (str(1:1) == '%') then
          if(index(str,' real ')/=0) havesh=0
       else
          backspace(foccn)
       endif
       iblu = 0
       do  ib = 1, nbas
          if (lldau(ib) /= 0) then
             is = ispec(ib)
             lmxa=sspec(is)%lmxa
             do l = 0,min(lmxa,3)
                if (idu(l+1,is) /= 0) then
                   iblu = iblu+1
                   do  isp = 1, 2
11                    continue
                      if (str(1:1) == '#') goto 11
                      read(foccn,*) nocc(-l:l,isp)
                      write(stdo,ftox)' occnum: site',ib,'l',l,'isp',isp,' ',ftof(nocc(-l:l,isp))
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
    ! ... Initial printout
    call praldm(0,51,51,havesh,nbas,nsp,lmaxu,lldau,' dmats read from disk',dmatu)
    ivsiz = nsp*nlibu*(lmaxu*2+1)**2
    call mpibc1(dmatu,2*ivsiz,4,mlog,'sudmtu','dmatu')
    call mpibc1(havesh,1,2,.false.,' ',' ')
    ! ... Density matrix in real or spherical harmonics (fixed by idvsh)
    if (havesh /= idvsh) then
       call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = idvsh
    endif
    ! ... Symmetrize dmatu (symdmu requires real harmonics)
    dmwk_zv=dmatu
    if (havesh == 1) then
       call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = 0
    endif
    call symdmu(nlibu,dmatu , nbas , nsp , lmaxu , ng , g , istab , lldau , xx )
    if (havesh /= idvsh) then
       call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu, lldau)
       call rotycs(2*idvsh-1,dmwk_zv,nbas,nsp,lmaxu, lldau )
       havesh = idvsh
    endif
    if (ng /= 0) then
       if(master_mpi)write(stdo,ftox)'sudmtu:  RMS change in dmats from symmetrization',ftof(xx)
       if (xx > .01d0) write(stdo,*)'(warning) RMS change unexpectely large'
       call daxpy ( ivsiz * 2 , - 1d0 , dmatu , 1 , dmwk_zv , 1 )
       if(ipr>=60) write(stdo,*)' change in dmat wrought by symmetrization'
       call praldm ( 0 , 60 , 60 , 0 , nbas , nsp , lmaxu , lldau ,' ' , dmwk_zv )
    endif
    !     Print dmats in specified harmonics
    dmwk_zv=dmatu
    if (havesh /= idvsh) then
       call rotycs ( 2 * idvsh - 1 , dmwk_zv , nbas , nsp , lmaxu, lldau )
    endif
    if(cmdopt0('--showdmat')) then
       if(master_mpi)write(stdo,*)
       call praldm(0,30,30,idvsh,nbas,nsp,lmaxu,lldau,' Symmetrized dmats' , dmwk_zv )
    endif   
    !     Print dmats in complementary harmonics
    i = 1-idvsh
    call rotycs(2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu, lldau )
    if(cmdopt0('--showdmat')) then
       if(master_mpi)write(stdo,*)
       call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau, ' Symmetrized dmats' , dmwk_zv )
    endif   
    ! ... Make Vorb (ldau requires spherical harmonics)
    if (havesh /= 1) then
       call rotycs(1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = 1
    endif
    if(master_mpi) write(stdo,*)
    iblu = 0
    do  20  ib = 1, nbas
       if (lldau(ib) == 0) goto 20
       is = ispec(ib) 
       lmxa=sspec(is)%lmxa
       spid=slabl(is) 
       i = min(lmxa,3)
       if(master_mpi) write(stdo,ftox)'Species '//spid//'mode',idu(1:i+1,is),'U',ftof(uh(1:i+1,is),2),'J',ftof(jh(1:i+1,is),2)
       do  22  l = 0, i
          if (idu(l+1,is) /= 0) then
             iblu = iblu+1
             uhx=uh(l+1,is)
             call ldau(idu(l+1,is),l,iblu,uhx,jh(l+1,is),dmatu,nsp,lmaxu,vorb,eorb)
          endif
22     enddo
20  enddo
    call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,' Unsymmetrized vorb',vorb)
    !!!! ==>     At this point, dmatu and vorb are in spherical harmonics (complex)
    
    ! ... Symmetrize vorb to check (symdmu requires real harmonics)
    call rotycs(-1,vorb,nbas,nsp,lmaxu,lldau)
    call symdmu (nlibu, vorb, nbas , nsp , lmaxu , ng , g , istab , lldau , xx )
    !     EITHER: vorb =>  spherical harmonics OR dmatu => real harmonics
    if (idvsh == 1) then
       call rotycs(1,vorb,nbas,nsp,lmaxu,lldau)
       havesh = 1
    else
       call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = 0
    endif
    if (master_mpi.and.ng /= 0) then
       write(stdo,ftox)' sudmtu:  RMS change in vorb from symmetrization = ',ftof(xx)
       if (xx > .01d0) write(stdo,*)'          (warning) RMS change unexpectely large'
       write(stdo,*)
    endif
    if(cmdopt0('--showdmat')) call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,' Symmetrized vorb',vorb)
    i = 1-idvsh
    dmwk_zv=vorb
    call rotycs( 2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , lldau )
    if(cmdopt0('--showdmat')) then
       write(stdo,*) ! vorb in complementary harmonics
       call praldm(0,30,30, i , nbas , nsp , lmaxu , lldau , ' Vorb' , dmwk_zv )
    endif   
    eorb = 0d0
    return
99  continue ! --- Error exit ---
    write(str,"('bad occnum file, site l= ',i5,i5)")ib,l
    call rx(trim(str))
  end subroutine sudmtu

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine rotycs(mode,a,nbas,nsp,lmaxu,lldau)
    use m_lmfinit,only:idu,ispec,sspec=>v_sspec
    !- Rotate matrix a from real to spherical harmonics
    ! for LDA+U objects densmat and vorb
    !-------------------------------------
    !i mode =1 from real a to spherical a
    !i      -1 from spherical a to real a
    !i a: matrix to be transformed a(m,m,isp,iblu)  could be vorb or dmat
    !i nbas : number of sites
    !i nsp  : number of spins
    !i lmaxu: lmax for U
    !i sspec: species info
    !i lldau  :lldau(ib)=0 => no U on this site otherwise
    !i        :U on site ib with dmat in dmats(*,lldau(ib))
    !o a rotated in place
    !r Remarks
    !r order of cubic harmonics ls (l-1)s,ms...1s 0 1c mc... (l-1)c lc
    !r order of spherical harmonics -l:l
    !r Yl-m=(Ylmc-iYlms)/sqrt(2)  Ylm=(-1)**m*conjg(Yl-m)
    !u Updates
    !u   18 Jan 06 A. Chantis changed rotation matrices in accordance with
    !u             the definition of real harmonics used in the rest of
    !u             the code (Hund's rules satisfied as indicated by orb. moment)
    !u   09 Nov 05 (wrl) Convert dmat to complex form
    !u   30 Apr 05 Lambrecht first created
    !----------------------------------------------------------------
    implicit none
    integer :: nbas,lldau(nbas),mode,lmaxu,nsp
    complex(8),target:: a(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nbas)
    integer:: ib,m,l,is,igetss,i,j,k,ll,isp,iblu
    complex(8):: b(2,2),c(2,2),add,bb(2,2)
    complex(8),parameter:: s2= 1/dsqrt(2d0), img=(0d0,1d0)
    complex(8),pointer:: rot(:,:)
    complex(8),save,target::rott(2,2,5,-1:1)
    logical,save:: init=.true.
    logical:: cmdopt0
    if(mode/=1 .AND. mode/=-1) call rx('ROTYCS: mode must be 1 or -1')
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
       is  = ispec(ib) !ssite(ib)%spec
       do  l = 0, min(sspec(is)%lmxa,3)
          if (idu(l+1,is) ==0) cycle
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
