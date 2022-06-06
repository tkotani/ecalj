!! == Main program for lmf part. Now we use MPI for all program
!  lmf-MPIK and lmfgw-MPIK
!     We use module-based programing.
!     In principle, all the data are generared and stored in some modules with 'protection'.
!     We can not modify data in a module by other modules (in prinicple, not everywhere yet).
program lmf
  use m_lmfinit,only:  M_Lmfinit_init,nlibu,plbnd,lrout,ctrl_lfrce,sspec=>v_sspec
  use m_writeham,only: M_writeham_init, M_writeham_write
  use m_ext,only:      M_Ext_init,sname
  use m_lattic,only:   M_Lattic_init
  use m_mkqp,only:     M_Mkqp_init,bz_nabc
  use m_rdfiln,only:   M_Rdfiln_init
  use m_MPItk,only:    M_MPItk_init, M_MPItk_finalize, nsize, master_mpi
  use m_hamindex, only:M_hamindex_init
  use m_hamindex0,only:M_hamindex0_init
  use m_supot,only:    M_supot_init
  use m_suham,only:    M_suham_init
  use m_qplist,only:   M_qplist_init, M_qplist_qspdivider, nkp
  use m_igv2x,only:    M_igv2xall_init
  use m_ldau,only:     M_ldau_init
  use m_gennlat,only:  M_gennlat_init
  use m_mksym,only:    M_mksym_init
  use m_lgunit,only:   M_lgunit_init, stdo,stdl
  use m_sugcut,only:sugcut
  implicit none
  integer:: k, iarg,jobgw,iprint,nit1,ifi,ifile_handle,nx,ny,nk1,nk2,nk3,i,j,ix
  logical:: fileexist,cmdopt0,cmdopt2, writeham,lbin,sigx
  character:: outs*20,aaa*512,sss*128
  character(8):: prgnam
  ! ===  Bootstrap building up of module variables by initialization routines m*_init().
  ! ===  In principle, all the fixed data required for calculation are pushed into modules by the initialzation.
  prgnam='LMF'
  if(cmdopt0('--jobgw')) prgnam='LMFGWD'
  call M_ext_init()         ! Get sname, e.g. trim(sname)=si of ctrl.si
  call M_MPItk_init(prgnam)
  call M_lgunit_init()
  aaa=''
  do iarg=1,iargc()
     call getarg(iarg,sss)
     aaa=aaa//' '//trim(sss)
  enddo
  aaa='===START '//trim(prgnam)//' with  '//trim(aaa)//' ==='
  if(master_mpi) call show_programinfo(stdo)
  if(master_mpi) write(stdo,"(a)") trim(aaa)
  if(master_mpi) write(stdl,"(a)") trim(aaa)
  if(master_mpi) write(stdo,*) 'mpisize=',nsize
  if(master_mpi) write(stdl,*) 'mpisize=',nsize
  if(cmdopt0('--help')) then  !help and quit
     call M_lmfinit_init(prgnam) ! show help and quit for --input
     call Rx0('end of help mode')
  endif

  !! --writepdos mode and exit
  if( Cmdopt0('--writepdos') ) then
     !! See job_pdos. New pdos mode (use --mkprocar and --fullmesh together).
     !! We use all k points (--fullmesh), instead of using crystal symmetry. See job_pdos
     if(master_mpi) write(6,*) '... Doing writepdos mode. Wait a while ...'
     if(master_mpi) write(6,*) '... See job_pdos to know how to call --writepdos mode'
     call Writepdos(trim(sname))
     call Rx0('done: end of --writepdos mode.')
  endif

  !! --wdsawada mode and exit
  if( Cmdopt0('--wdsawada') ) then !! Sawada's simple mode and exit
     if(master_mpi)write(6,*) '... write Dos from tetraf.dat and eigenf.dat. '
     if(master_mpi)write(6,*) '...  eigenf.dat is for qplistf.dat '
     call Writedossawada()
     call Rx0('done: end of --wdsawada mode.')
  endif

  !! Read ctrl.* into m_rdfiln
  ! m_rdfiln is too complicated to maintain. We will use simple reader.
  ! 'math operation' and '-v substitution' may be convenient, but probably we will separate them away.
  !  'math operation' allows expression such as ALAT=7.88*1.1 in ctrl file
  call M_rdfiln_init() ! Read recrd of ctrl file into m_lmfinit
  ! Set all the initial conditions in the module m_lmfinit. All variables except v_sspec and v_ssite are protected.
  call M_lmfinit_init(prgnam) ! Computational settings from ctrl file is stored in m_lmfinit
  call M_lattic_init()      ! lattice setup (for ewald sum)
  call M_mksym_init(prgnam) !symmetry go into m_lattic and m_mksym
  if(trim(prgnam)=='LMF') call M_mkqp_init() ! data of BZ go into m_mkqp
  !!======================================================================
  !! Main program. lmf-MPIK/lmfgw-MPIK. 'call lmfp'---------------------
  !! --rs=3 mode is removed. (--rs=3 meand fixed density Harris-foukner MD).
  !! Sep2020 comment " Shorten site positions" removed. (we are useing shortn3 mainly now)
  ! Get jobgw for lmfgw mode. Quit for job=0
  jobgw=-1
  if(prgnam =='LMFGWD') then
     if( cmdopt2('--job=',outs) ) then
        read(outs,*) jobgw
     elseif(master_mpi) then
        write(stdo,*)
        write(stdo,*) ' === lmfgw-MPIK: Choose one of following jobs: ==='
        write(stdo,*) '  0 : init mode creates HAMindex0'
        write(stdo,*) '  1 : GW setup for eigenfunctions: HAMindex gwa.* gwb.*'
        write(stdo,*) '  jobgw 0 or 1?'
        read (5,*) jobgw
     endif
     call Mpibc1_int(jobgw,1,'lmfp_jobgw')
  endif
  if(jobgw==0) then ! Index for hamiltonian gen_hamindex
     if(master_mpi) call M_hamindex0_init()
     call Rx0(' OK! '//'lmfgw mode=0 generated HAMindex0')
  endif
  !! Array allocated in supot rhoat smrho.
  call M_supot_init() ! get G vectors for charge
  call Sugcut(1)
  call M_suham_init()   ! Get estimated dimension of Hamiltonian (probably simplified in future).
  if(nlibu>0) call M_ldau_init() !! LDA+U initialization
  if(cmdopt0('--quit=dmat')) call Rx0('--quit=dmat')
  call M_qplist_init(plbnd,jobgw==1) ! Get q point list at which we do band calculations
  call M_qplist_qspdivider()  !generate iqini:iqend,isini,isend  for each rank
  call M_igv2xall_init(1,nkp) !G vectors for qplist. (1,nkp) is needed for gwb.head

  !! Madelung mode here may need to be recovered if necessary.
  !!  allocate(madelung(nbas**2)); call madmat(madelung) !Monopole Madelung matrix (kept for future).
  !! Shear mode here is currently commented out.=>probably shear mode should be outside of fortran.

  call M_hamindex_init(jobgw) !bugfix. moved here 2022apr04
  writeham= cmdopt0('--writeham') !m_writeham_* is after m_hamindex_init 2022apr22
  if(writeham .AND. master_mpi) call M_gennlat_init(bz_nabc) !for interpolation of Hamiltonian
  if(writeham .AND. master_mpi) call M_writeham_init()
  if(writeham .AND. master_mpi) call M_writeham_write()
  if( cmdopt0('--quit=ham') ) call Rx0('quit = ham')

  !! (we need check. a simple approximaiton to determine VBM and CBM. Need fixing if necessary).
  if(cmdopt0('--vbmonly')) then !Get VBM and CBM relative to vaccum
     if(master_mpi) call Vbmmode()
     call Rx0('--vbmonly mode done')
  endif
  if(cmdopt0('--getq')) then ! Current version is not for spin dependent, with many restrictions.
     if(master_mpi) call Getqmode()
     call Rx0('--getq mode done')
  endif
  call Praugm(sspec,0) !print out properties of species
  !! Main rouitne.
  call Lmfp(jobgw==1)
  call Rx0("OK! end of "//trim(prgnam)//" ======================")
end program Lmf

include "show_programinfo.fpp" !this is preprocessed from show_programinfo.f90 by Makefile
