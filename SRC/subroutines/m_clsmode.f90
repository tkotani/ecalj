module m_clsmode
  !! CLS: Core-level spectroscopy !We use CLSinput instead of --cls option.
  use m_lmfinit, only: lmet=>bz_lmet,nbas,nsp,ssite=>v_ssite,sspec=>v_sspec,nlmax,nspc,nl,lso,stdo
  use m_suham,only:   ndham=>ham_ndham
  use m_mkqp,only: nkp=>bz_nkp
  use m_MPItk,only: master_mpi
  use m_mkpot,only: ppnl_rv
  use m_ftox,only: ftox
  integer,parameter,private:: nsitmx = 256
  integer,private::  icls=0 , isite(nsitmx) , iclsl(nsitmx), iclsn(nsitmx),nsites
  complex(8),allocatable,private :: ausc_zv(:)
contains

  subroutine m_clsmode_init()
    logical:: cmdopt0
    character(10):: i2char
    integer::nevmx,ndhamx,i,ific
    character strn*120, clsopt*120
    character(512):: aaachar
    !! --- Options for core level specta (CLS) ---
    if (cmdopt0('--cls')) then
       !if (lmet/=2) call rx('For CLS restart with METAL=2')
       icls = 1
       !     clsopt = strn(6:)
       !     call suclst(nsitmx,nbas,nsp,ssite,sspec,clsopt,
       !     .     isite,iclsl,iclsn,nsites)
       !     open(newunit=ific,file='CLSinput') !2022apr28
       !     do i=1,nsites
       !     write(ific,ftox)isite(i),iclsl(i),iclsn(i),' ! site,core-l,core-n'
       !     enddo
       !     close(ific)
       open(newunit=ific,file='CLSinput') !2022apr28
       i=0
       do
          i=i+1
          read(ific,*,err=1018,end=1018) isite(i),iclsl(i),iclsn(i)
          !     write(stdo,ftox)'xxx clsmode for',isite(i),iclsl(i),iclsn(i),' ! site,core-l,core-n'
       enddo
1018   continue
       close(ific)
       nsites=i-1
       do i=1,nsites
          write(stdo,ftox)'clsmode for',isite(i),iclsl(i),iclsn(i),' ! site,core-l,core-n'
       enddo

       if (16*3*nlmax*ndham*nbas*nsp*nkp/1000000 > 24 .AND. master_mpi) then
          !     ! this path is go through by crn test.
          aaachar= &
               ' CLS: '//trim(i2char( 16*3*nlmax*ndham*nsites*nsp*nkp/1000000))// &
               ' Mb memory for aus: nlmax='//trim(i2char( nlmax))// &
               ' ndham='  //trim(i2char( ndham))// &
               ' nsistes='//trim(i2char( nsites))// &
               ' nsp='    //trim(i2char( nsp))// &
               ' nkp='    //trim(i2char( nkp))
          write(6,"(a)") aaachar
       endif
       allocate(ausc_zv(3*nlmax*ndham*nsites*nsp*nkp))
       ausc_zv=0d0
    else
       icls = 0
    endif
  end subroutine m_clsmode_init

  subroutine m_clsmode_set1(nmx,jsp,iq,qp,nev,t_zv)
    use m_igv2x,only: ndimhx
    integer:: nmx,jsp,iq,nev
    real(8)::qp(3)
    complex(8):: t_zv(1:ndimhx,1:nmx)
    call rxx(lso==1,'CLS not implemented in noncoll case')
    call makusq(0, nsites,isite, nev,jsp,iq,qp,t_zv, ausc_zv)
  end subroutine m_clsmode_set1

  subroutine m_clsmode_finalize(bz_ef,ndimh,ndhamx,nspx,nkp,dosw,evlall)
    integer:: ndimh,ndhamx,nspx,nkp
    real(8),intent(in) :: bz_ef,dosw(2)
    real(8):: eferm,evlall(ndhamx,nspx,nkp)
    eferm=bz_ef
    call vcdmel ( nl , ssite , sspec ,  nlmax , ndham ,&
         ndimh , nkp , nsp , nspc , eferm , evlall , ausc_zv , &
         nsites , isite , iclsl , iclsn ,dosw) !sbz,
    call rx0('done generating core level spectra')
  end subroutine m_clsmode_finalize
end module m_clsmode
