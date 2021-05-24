      module m_clsmode
!! CLS: Core-level spectroscopy !too confusing notation see suclst.
      use m_lmfinit, only: lmet=>bz_lmet,nbas,nsp,ssite=>v_ssite,sspec=>v_sspec,nlmax,nspc,nl
      use m_suham,only:   ndham=>ham_ndham
      use m_mkqp,only: nkp=>bz_nkp
      use m_MPItk,only: master_mpi
      use m_mkpot,only: ppnl_rv 

      
      integer,parameter,private:: nsitmx = 256
      integer,private::  icls=0 , isite(nsitmx) , iclsl(nsitmx), iclsn(nsitmx),nsites
      complex(8),allocatable,private :: ausc_zv(:)
      contains

      subroutine m_clsmode_init()
      logical:: cmdopt
      character(10):: i2char
      integer::nevmx,ndhamx
      character strn*120, clsopt*120
      character(512):: aaachar
!! --- Options for core level specta (CLS) ---
      if (cmdopt('--cls',5,0,strn)) then
         if (lmet/=2) call rx('For CLS restart with METAL=2')
         icls = 1
         clsopt = strn(6:)
         call suclst(nsitmx,nbas,nsp,ssite,sspec,clsopt,
     .     isite,iclsl,iclsn,nsites)
         if (16*3*nlmax*ndham*nbas*nsp*nkp/1000000 .gt. 24 .and. master_mpi) then
!! this path is go through by crn test.
            aaachar=
     &        ' CLS: '//trim(i2char( 16*3*nlmax*ndham*nsites*nsp*nkp/1000000))//
     &        ' Mb memory for aus: nlmax='//trim(i2char( nlmax))//
     &        ' ndham='  //trim(i2char( ndham))//
     &        ' nsistes='//trim(i2char( nsites))//
     &        ' nsp='    //trim(i2char( nsp))//
     &        ' nkp='    //trim(i2char( nkp))
            write(6,"(a)") aaachar
         endif
         allocate(ausc_zv(3*nlmax*ndham*nsites*nsp*nkp))
         ausc_zv=0.0d0
      else
         icls = 0
      endif
      end subroutine m_clsmode_init

      subroutine m_clsmode_set1(napw,ndimh,ndimhx,nmx,jsp,iq,qp,nev,igv2x,t_zv)
      integer:: ndimhx,nmx,jsp,iq,napw,ndimh,igv2x(3,napw),nev
      real(8)::qp(3)
      complex(8):: t_zv(1:ndimhx,1:nmx)      
      call makusq ( 0 , ssite , sspec ,  nbas , nsites 
     .              , isite , nlmax , ndham , ndimh , napw , igv2x , nev , 
     .              nsp , nspc , jsp , iq , qp , t_zv , ppnl_rv , ausc_zv   )
      end subroutine
      
      subroutine m_clsmode_finalize(bz_ef,ndimh,ndhamx,nspx,nkp,dosw,evlall)
      integer:: ndimh,ndhamx,nspx,nkp
      real(8),intent(in) :: bz_ef,dosw(2)
      real(8):: eferm,evlall(ndhamx,nspx,nkp)
      eferm=bz_ef
      call vcdmel ( nl , ssite , sspec ,  nlmax , ndham , !slat ,
     .        ndimh , nkp , nsp , nspc , eferm , evlall , ausc_zv , 
     .        nsites , isite , iclsl , iclsn ,dosw) !sbz,
            call rx0('done generating core level spectra')
      end subroutine
      end module m_clsmode
