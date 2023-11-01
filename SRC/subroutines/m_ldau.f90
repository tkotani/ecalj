!ldau density matrix and potential
module m_ldau 
  use m_ldau_util,only: chkdmu,sudmtu,mixmag
  use m_mpitk,only:master_mpi
  use m_ftox
  public:: m_ldau_init, m_ldau_vorbset
  complex(8),allocatable,protected,public::  vorb(:,:,:,:) !potential
  real(8),protected,public:: eorb=0d0 ! energy of U term
  private
  complex(8),allocatable,private::  dmato(:,:,:,:) !previous density matrix
  logical,private:: init=.true.
contains
  subroutine m_ldau_init() ! LDA+U initialization
    use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu
    integer:: i,idmatu
    logical:: mmtargetx
    call tcn('m_ldau_init')
    idmatu=nsp*nlibu*(lmaxu*2+1)**2
    if(init) then
       allocate( dmato(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
       allocate( vorb (-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
       init=.false.
    endif
    call sudmtu(dmato, vorb)  !Read dmatu from dmatu.ext and generate vorbdmat from dmatu.ext
    inquire(file='mmtarget.aftest',exist=mmtargetx)
    if(mmtargetx) call vorbmodifyaftest_experimental()
    call tcx('m_ldau_init')
  end subroutine m_ldau_init
  subroutine m_ldau_vorbset(eks,dmatu) ! Get vorb,eorb for given dmatu, and reserve dmato=previous dmatu.
    use m_lgunit,only:stdo
    use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu,nbas
    intent(in)::            eks,dmatu
    complex(8):: dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    complex(8):: vorbav(-lmaxu:lmaxu,-lmaxu:lmaxu)
    real(8):: eks
    integer,parameter:: nx=1000
    logical:: mmtargetx
    real(8):: alpha,mmtarget, uhhist(nx),uhx,sss(1),fac,mmsite(nbas),mmhist(nx)
    integer:: nit,ibas,ifx,ncount=0,iblu,i
    call tcn('m_ldau_vorbset')
    call chkdmu(eks,dmatu, dmato,vorb,eorb) !eorb is returned internally.
    dmato=dmatu !dmato is kept in this module as previous dmat
    inquire(file='mmtarget.aftest',exist=mmtargetx)
    if(mmtargetx) call vorbmodifyaftest_experimental()
    call tcx('m_ldau_vorbset')
  end subroutine m_ldau_vorbset
  subroutine vorbmodifyaftest_experimental()
    use m_lgunit,only:stdo
    use m_lmfinit,only: nlibu,nsp,lmaxu,lmaxu,nsp,nlibu,nbas
    complex(8):: vorbav(-lmaxu:lmaxu,-lmaxu:lmaxu)
    integer,parameter:: nx=1000
    real(8):: alpha,mmtarget, uhhist(nx),uhx,sss(1),fac,mmsite(nbas),mmhist(0:nx)
    integer:: nit,ibas,ifx,ncount=0,iblu,i
    logical:: init=.true.
    character(3):: stat
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
       read(ifx,*,end=1112,err=1112) uhx
       read(ifx,*,end=1112,err=1112) (mmsite(ibas),ibas=1,nbas)
       nit=nit+1
       mmhist(nit)=(mmsite(1)-mmsite(2))/2d0
       uhhist(nit)=uhx
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
    if(nlibu==6) then !experimental part for a specific system NiSe and so on.
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
  end subroutine vorbmodifyaftest_experimental
end module m_ldau

