!>CLS: Core-level spectroscopy !We use CLSinput instead of --cls option.
module m_clsmode 
  use m_lmfinit, only: lmet=>bz_lmet,nbas,nsp,nlmax,nspc,lso
  use m_igv2x,only: nbandmx
  use m_lgunit,only:stdo
  use m_mkqp,only: nkp=>bz_nkp
  use m_MPItk,only: master_mpi
  use m_ftox,only: ftox
  integer,parameter,private:: nsitmx = 256
  integer,private::  icls=0 , isite(nsitmx) , iclsl(nsitmx), iclsn(nsitmx),nsites
  complex(8),allocatable,private :: ausc_zv(:),ausc(:)
contains
  subroutine m_clsmode_init()
    logical:: cmdopt0
    character(10):: i2char
    integer::i,ific
    character strn*120, clsopt*120
    character(512):: aaachar
    !! --- Options for core level specta (CLS) ---
    if (cmdopt0('--cls')) then
       icls = 1
       open(newunit=ific,file='CLSinput') !2022apr28
       i=0
       do
          i=i+1
          read(ific,*,err=1018,end=1018) isite(i),iclsl(i),iclsn(i)
       enddo
1018   continue
       close(ific)
       nsites=i-1
       do i=1,nsites
          write(stdo,ftox)'clsmode for',isite(i),iclsl(i),iclsn(i),' ! site,core-l,core-n'
       enddo
       if (16*3*nlmax*nbandmx*nbas*nsp*nkp/1000000 > 24 .AND. master_mpi) then ! this path is go through by crn test.
          aaachar= &
               ' CLS: '//trim(i2char( 16*3*nlmax*nbandmx*nsites*nsp*nkp/1000000))// &
               ' Mb memory for aus: nlmax='//trim(i2char( nlmax))// &
               ' nbandmx='  //trim(i2char( nbandmx))// &
               ' nsistes='//trim(i2char( nsites))// &
               ' nsp='    //trim(i2char( nsp))// &
               ' nkp='    //trim(i2char( nkp))
          write(6,"(a)") aaachar
       endif
       allocate(ausc_zv(3*nlmax*nbandmx*nsites*nsp*nkp))
       allocate(ausc(3*nlmax*nbandmx*nsites*nsp*nkp))
       ausc_zv=0d0
    else
       icls = 0
    endif
  end subroutine m_clsmode_init

  subroutine m_clsmode_set1(nmx,jsp,iq,qp,nev,t_zv)
    use m_igv2x,only: ndimhx
    use m_makusq,only: makusq
    integer:: nmx,jsp,iq,nev
    real(8)::qp(3)
    complex(8):: t_zv(1:ndimhx,1:nmx)
    call rxx(lso==1,'CLS not implemented in noncoll case')
    call makusq(nsites,isite, nev,jsp,iq,qp,t_zv, ausc)!ausc_zv is accumulating
    ausc_zv = ausc_zv + ausc
  end subroutine m_clsmode_set1

  subroutine m_clsmode_finalize(bz_ef,ndimh,nbandmx,nspx,nkp,dosw,evlall)
    use m_vcdmel,only:vcdmel
    integer:: ndimh,nbandmx,nspx,nkp
    real(8),intent(in) :: bz_ef,dosw(2)
    real(8):: eferm,evlall(nbandmx,nspx,nkp)
    call mpibc2_complex(ausc_zv,size(ausc_zv),'ausc_zv') !2023jan fixed for TestInstall/crn
    eferm=bz_ef
    if(master_mpi) call vcdmel(nlmax,ndimh,nbandmx,nspx,nkp,nsp,nspc,eferm,evlall,ausc_zv,nsites,isite,iclsl,iclsn,dosw)
  end subroutine m_clsmode_finalize
end module m_clsmode
