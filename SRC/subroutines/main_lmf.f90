! Main program for lmf part. Now we use MPI for all program based on the paper
! [1] TK. H.Kino H.Akai, JPSJ84,034702(2015): We will write newer paper (2023plan)
!   Cite ecalj/Document/PAPERandPRESENTATION/KotaniKinoAkai2015FormulationPMT.pdf
!
! We use a module-based programing. In principle, all the variables are generared and stored in modules with 'protection'.
! This assure that we can not modify data in a module by other modules.
! Bootstrap sequence of module initialzation. The variables in modules are proteted except m_density. Use variables with 'use only'.
module m_lmf
contains
  subroutine lmf(commin) bind(C)
    use m_args,only:     argall
    use m_ext,only:      m_ext_init,sname
    use m_MPItk,only:    m_MPItk_init, nsize, master_mpi
    use m_lgunit,only:   m_lgunit_init, stdo,stdl
    use m_cmdpath,only:  setcmdpath
    use m_lmfinit,only:  m_lmfinit_init,nlibu,plbnd,nbas
    use m_lattic,only:   m_lattic_init,readatompos,  qlat=>lat_qlat
    use m_mksym,only:    m_mksym_init
    use m_mkqp,only:     m_mkqp_init,bz_nabc
    use m_lattic,only:   Setopos
    use m_hamindex0,only:m_hamindex0_init
    use m_supot,only:    m_supot_init
    use m_sugcut,only:   sugcut
    use m_suham,only:    m_suham_init
    use m_ldau,only:     m_ldau_init
    use m_qplist,only:   m_qplist_init, m_qplist_qspdivider, nkp
    use m_igv2x,only:    m_igv2xall_init
    use m_hamindexW, only:m_hamindexW_init
    use m_hamindex, only:readhamindex
    use m_gennlat,only:  m_gennlat_init
    use m_writeham,only: m_writeham_init, m_writeham_write
    use m_lmfp,only:     lmfp  !this is main part lmfp-->bndfp
    use m_writeband,only: writepdos,writedossawada
    use m_ftox
    implicit none
    integer,optional:: commin
    integer:: iarg,iprint,jobgw=-1,ierr,ifi
    logical:: cmdopt0,cmdopt2, writeham,sigx
    character:: outs*20,aaa*512,sss*128
    character(32):: prgnam
    integer:: comm
    include "mpif.h" 
    comm = MPI_COMM_WORLD
    if(present(commin)) comm= commin  
    if(cmdopt2('--jobgw=',outs))then
       prgnam='LMFGWD' !GW set up mode
       read(outs,*) jobgw
       if(jobgw/=0.and.jobgw/=1) call rx0(' Set --jobgw=0 or 1')
    else
       prgnam='LMF'
    endif
    call m_MPItk_init(comm)  ! MPI info
    call m_ext_init()    ! Get sname, e.g. trim(sname)=si of ctrl.si
    call m_lgunit_init() ! Set file handle of stdo(console) and stdl(log)    !print *, 'len_trim(argall)=',trim(argall),len_trim(argall),master_mpi
    aaa='===START '//trim(prgnam)//' with  '//trim(argall)//' ==='
    if(master_mpi) call show_programinfo(stdo)
    if(master_mpi) write(stdo,"(a)") trim(aaa)
    if(master_mpi) write(stdl,"(a)") trim(aaa)
    if(master_mpi) write(stdo,"(a,g0)")'mpisize=',nsize
    if(master_mpi) write(stdl,"(a,g0)")'mpisize=',nsize
    !  call setcmdpath()         ! Set self-command path (this is for call system at m_lmfinit)
    WritePdosmode: if( cmdopt0('--writepdos') ) then  ! See job_pdos. New pdos mode (use --mkprocar and --fullmesh together).
       ! We use all k points (--fullmesh), instead of using crystal symmetry. See job_pdos
       if(master_mpi) write(stdo,*) '... Doing writepdos mode. Wait a while ...'
       if(master_mpi) write(stdo,*) '... See job_pdos to know how to call --writepdos mode'
       call writepdos(trim(sname))
       call rx0('done: end of --writepdos mode.')
    endif WritePdosmode
    WriteDOSsawadamode:  if( cmdopt0('--wdsawada') ) then !! Sawada's simple mode and exit
       if(master_mpi)write(stdo,*) '... write Dos from tetraf.dat and eigenf.dat. '
       if(master_mpi)write(stdo,*) '...  eigenf.dat is for qplistf.dat '
       call writedossawada()
       call rx0('done: end of --wdsawada mode.')
    endif WriteDOSsawadamode
    if(master_mpi) then
       open(newunit=ifi,file='save.'//trim(sname),position='append')
       write(ifi,"(a)")'Start '//trim(prgnam)//trim(argall)
       close(ifi)
    endif   
    if(master_mpi) call ConvertCtrl2CtrlpByPython() !convert ctrl file to ctrlp.
    if(cmdopt0('--quit=ctrlp')) call rx0('--quit=ctrlp')
    call MPI_BARRIER( comm, ierr)
    call m_lmfinit_init(prgnam,comm)! Read ctrlp into module m_lmfinit.
    call m_lattic_init()       ! lattice setup (for ewald sum)
    call m_mksym_init()  !symmetry go into m_lattic and m_mksym
    if(trim(prgnam)=='LMF') call m_mkqp_init() ! data of BZ go into m_mkqp
    ! Sep2020:  Shorten site positions" removed. (we are useing shortn3 mainly now)
    ReadingPos: block !read atomic positions from AtomPos if it exists. Overwrite pos in m_lattic.
      logical:: irpos
      call ReadAtomPos(irpos)
      if(irpos) write(stdo,*) 'Readin AtomPos.'//trim(sname)//' !!!!'
    endblock ReadingPos
    if(jobgw==0) then ! Index for hamiltonian gen_hamindex ! Get jobgw for lmfgw mode. Quit for job=0
       if(master_mpi) call m_hamindex0_init()
       call rx0(' OK! '//'lmfgw mode=0 generated HAMindex0')
    endif
    call m_supot_init() ! get G vectors for charge ! Array allocated in supot rhoat smrho.
    call sugcut(1)
    call m_suham_init()   ! Get estimated dimension of Hamiltonian (probably simplified in future).
    if(nlibu>0) call m_ldau_init() !! LDA+U initialization
    if(cmdopt0('--quit=dmat')) call rx0('--quit=dmat')
    call m_qplist_init(plbnd,jobgw==1) ! Get q point list at which we do band calculations
    call m_qplist_qspdivider()  !generate iqini:iqend,isini,isend  for each rank
    call m_igv2xall_init(1,nkp) !G vectors for qplist. (1,nkp) is needed for gwb.head
    !Madelung mode here may need to be recovered if necessary.
    !   allocate(madelung(nbas**2)); call madmat(madelung) !Monopole Madelung matrix (kept for future).
    !Shear mode here is currently commented out.=>probably shear mode should be outside of fortran.
    inquire(file='sigm.'//trim(sname),exist=sigx)
    if(jobgw>=0.or.sigx) call m_hamindexW_init() !write HAMindex for GW part, or for reading sigm
    if(sigx) call readhamindex() !for reading sigm
    writeham= cmdopt0('--writeham') !m_writeham_* is after m_hamindex_init 2022apr22
    if(writeham .AND. master_mpi) call m_gennlat_init(bz_nabc) !for interpolation of Hamiltonian
    if(writeham .AND. master_mpi) call m_writeham_init()
    if(writeham .AND. master_mpi) call m_writeham_write()
    if( cmdopt0('--quit=ham') ) call rx0('quit = ham')
    if(cmdopt0('--vbmonly')) then !Get VBM and CBM relative to vaccum ! (a simple approximaiton to determine VBM and CBM. Need fixing if necessary).
       if(master_mpi) call Vbmmode()
       call rx0('--vbmonly mode done')
    endif
    if(cmdopt0('--getq')) then ! Current version is not for spin dependent, with many restrictions.
       if(master_mpi) call Getqmode()
       call rx0('--getq mode done')
    endif
    MainRoutine: block
      call lmfp(jobgw==1) 
    endblock MainRoutine
    if(master_mpi) write(stdo,"(a)")"OK! end of "//trim(prgnam)//" ======================"
  end subroutine lmf
end module m_lmf

subroutine ConvertCtrl2CtrlpByPython()
  use m_MPItk,only: master_mpi
  use m_args,only: argall
  use m_ext,only :sname        ! sname contains extension. foobar of ctrl.foobar
  use m_cmdpath,only:cmdpath
  use m_lgunit,only: stdo,stdl
  implicit none
  character(512):: cmdl
  logical:: fileexist
  integer::ifi
  inquire(file='ctrl.'//trim(sname),exist=fileexist)
  if(.NOT.fileexist) call rx("No ctrl file found!! ctrl."//trim(sname))
  cmdl=trim(cmdpath)//'ctrl2ctrlp.py '//trim(argall)//'<ctrl.'//trim(sname)//' >ctrlp.'//trim(sname)
  write(stdo,"(a)")'cmdl for python='//trim(cmdl)
  call system(cmdl) !See  results ctrlp.* given by ctrl2ctrl.py 
end subroutine ConvertCtrl2CtrlpByPython

!include "../exec/show_programinfo.fpp" !this is for 'call show_programinfo' ! preprocessed from show_programinfo.f90 by Makefile
