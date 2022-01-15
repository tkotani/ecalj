!! == Main program for lmf part. Now we use MPI for all program
!  lmf-MPIK and lmfgw-MPIK
!     We use module-based programing.
!     In principle, all the data are generared and stored in some modules with 'protection'.
!     We can not modify data in a module by other modules (in prinicple, not everywhere yet).
      program lmf
      use m_lmfinit,only:m_Lmfinit_init,nlibu,plbnd,irs1,irs1x10,lrout,ctrl_lfrce,sspec=>v_sspec
      use m_writeham,only: m_writeham_init, m_writeham_write
      use m_ext,only:      m_Ext_init,sname
      use m_lattic,only:   m_Lattic_init, Setopos
      use m_mkqp,only:     m_Mkqp_init
      use m_rdfiln,only:   m_Rdfiln_init
      use m_MPItk,only:    m_MPItk_init,m_MPItk_finalize,nsize,master_mpi
      use m_lmfgw,only:    m_lmfgw_init
      use m_hamindex,only: m_hamindex_init
      use m_supot,only:    m_supot_init
      use m_suham,only:    m_suham_init
      use m_qplist,only:   m_qplist_init,m_qplist_qspdivider,nkp
      use m_igv2x,only:    m_igv2xall_init
      use m_ldau,only:     m_ldau_init
      use m_iors,only:  Iors
      use m_lmaux,only: lmaux
      implicit none
      integer:: k,lgunit,stdo,stdl
      logical:: fileexist,cmdopt,llmfgw=.false.,cmdopt0,cmdopt2
      logical:: writeham,lbin
      character:: outs*20,aaa*512,sss*128
      character prgnam*8
      character*8 charext
      integer:: iarg,jobgw,iprint,nit1,ifi,irs1x,ifile_handle
#if LMFGWD
      data prgnam /'LMFGWD'/
#else
      data prgnam /'LMF'/
#endif
      call M_ext_init()        ! Get sname, e.g. trim(sname)=si of ctrl.si
      call M_MPItk_init(prgnam)
      stdo = lgunit(1)          !standard output stdo=6 usually
      stdl = lgunit(2)
      aaa=''
      do iarg=1,iargc()
         call getarg(iarg,sss)
         aaa=aaa//' '//trim(sss)
      enddo
      aaa= '===START '//trim(prgnam)//'  '//trim(aaa)//' ==='
      if(master_mpi) write(stdo,"(a)") trim(aaa)
      if(master_mpi) write(stdl,"(a)") trim(aaa)
      if(master_mpi) write(stdo,*) 'mpisize=',nsize
      if(master_mpi) write(stdl,*) 'mpisize=',nsize
      if(cmdopt2('--pr=',outs)) then ! Set iprint() from cmd line ; the verbosity of console out put.
        read(outs,*) k
        call Setpr(k)
      endif
      if(.not.master_mpi) call setpr(-100) !iprint() is negative except master
      if (cmdopt0('--help')) then  !help and quit
         call M_lmfinit_init(prgnam) ! show help and quit for --input
         call Rx0('end of help mode')
      endif
!! New pdos mode (use --mkprocar and --fullmesh together). and exit ===
!! We use all k points (--fullmesh), instead of using crystal symmetry. See job_pdos
      if( Cmdopt0('--writepdos') ) then
         write(6,*) '... Doing writepdos mode. Wait a while ...'
         write(6,*) '... See job_pdos to know how to call --writepdos mode'
         call Writepdos(trim(sname))
         call Rx0('done: end of --writepdos mode.')
      endif
!! Sawada's simple mode and exit
      if( Cmdopt0('--wdsawada') ) then
         write(6,*) '... write Dos from tetraf.dat and eigenf.dat. '
         write(6,*) '...  eigenf.dat is for qplistf.dat '
         call Writedossawada()
         call Rx0('done: end of --wdsawada mode.')
      endif
      
!! Read ctrl.* into m_rdfiln
!! m_rdfiln is too complicated to maintain. We will use simple reader.
!! 'math operation' and '-v substitution' may be convenient, but probably we will separate them away.
!!    it allow expression like ALAT=7.88*1.1
      call M_rdfiln_init()
!! Read recrd of ctrl file into m_lmfinit
!!  Set all the initial conditions in the module m_lmfinit.
!!  All variables except v_sspec and v_ssite are protected.
      call M_lmfinit_init(prgnam) ! Computational settings go into m_lmfinit
!! Lattice and Symmetry
      call Latsym(prgnam) !---> lattice(for ewald sum) and symmetry go into m_lattic and m_mksym
!! q points info     
      if(trim(prgnam) == 'LMF') then
         call M_mkqp_init() ! data of BZ go into m_mkqp
      endif
!!======================================================================      
!! Main program. lmf-MPIK/lmfgw-MPIK. 'call lmfp'---------------------
!! --rs=3 mode is removed. (--rs=3 meand fixed density Harris-foukner MD).
!! Sep2020 comment " Shorten site positions" removed. (we are useing shortn3 mainly now)
!!Get jobgw for lmfgw mode. Quit for job=0
      if(prgnam =='LMFGWD') llmfgw=.true.
      jobgw=0
      if(llmfgw) then
      if( cmdopt2('--job=',outs) ) then 
        read(outs,*) jobgw
        elseif(master_mpi) then
         write(stdo,*)
         write(stdo,*) '=== lmfgw-MPIK: Choose one of following jobs: ==='
         write(stdo,*) '   0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima'
         write(stdo,*) '   1 : GW setup mode'
         write(stdo,*) ' job 0 or 1?'
         read (5,*) jobgw
        endif
        call mpibc1_int(jobgw,1,'lmfp_jobgw')
      endif
!! Index for hamiltonian gen_hamindex
      call m_hamindex_init(jobgw==0) 
      if(llmfgw.and.jobgw==0) call m_lmfgw_init() 
      if(llmfgw.and.jobgw==0) call rx0(' OK! '//'lmfgw mode=0 generated LATTC CLASS NLAindx ldima')
!! Array allocated in supot rhoat smrho. We have to store all fixed quantites in modules.
      call m_supot_init() ! get G vectors for charge
      call sugcut(1) 
      call m_suham_init() ! Get estimated dimension of Hamiltonian (probably simplified in future).
      writeham= cmdopt0('--writeham')
      if(writeham.and.master_mpi) call m_writeham_init()
      if(writeham.and.master_mpi) call m_writeham_write()
      if( cmdopt0('--quit=ham') ) call rx0('quit = ham')
!! LDA+U initialization
      if(nlibu>0) call m_ldau_init()
      if(cmdopt0('--quit=dmat')) call rx0('--quit=dmat')
      write(6,*) 'end of reading rst or atm'
!! Get q point list at which we do band calculations
      call m_qplist_init(plbnd,llmfgw)
      call m_qplist_qspdivider()  !generate iqini:iqend,isini,isend  for each rank
!! G vectors for qplist
      call m_igv2xall_init(1,nkp)  !(iqini,iqend)   !(1,nkp) is needed for gwb.head
!! madelung mode may need to be recovered.      
!! allocate(rv_a_omad(nbas**2)); call madmat (rv_a_omad ) !Monopole Madelung matrix (kept for future).
c      
!!=== NOTE: shear mode is currently commented out.
!!   (probably shear mode should be recovered outside of fortran).
c
!! Initial density Read atm.* or rst.*  rhoat smrho in supot
!!  --- Read restart file or overlap free atom densities ---
!!    irs1x tells what to read and whether to invoke smshft.
!!     0    read from atom file  atm;  1  read from binary  rst ;  2 read from ascii  rsta
!!    +10 -> invoke smshft(1) after file read.
      write(6,*) 'goto end of reading rst or atm',irs1x
      irs1x = irs1
      if (irs1x == 0) then ! read overlap free-atom densities
         call rdovfa()     ! superposition of spherical atomic densities by lmfa.
         nit1 = 0
      elseif( mod(irs1x,4)==1 .or. mod(irs1x,4)==2 ) then 
         lbin = (irs1x==1) ! lbin=T:rst for --rs=1, or lbin=F:rsta for --rs=2,
         k = -1
         if(master_mpi) then
            ifi=ifile_handle() 
            if(lbin)       open(ifi,file='rst.'//trim(sname),form='unformatted')
            if(.not. lbin) open(ifi,file='rsta.'//trim(sname)) 
         endif
         call mpibc1_int(ifi,1,'lmfp_ifi')
         k = iors(nit1 , lbin, ifi=ifi )! read rst file. sspec ssite maybe modified 
         if (master_mpi) close(ifi) 
         call mpibc1_int(k,1,'lmfp_k')
         call Setopos() !set pos read in iors into rv_a_opos
         if(k<0) then !initial potential from atm file if rst can not read
            irs1x = 0
            call rdovfa() 
            nit1 = 0
         elseif (irs1x10) then  
            call smshft(1)
         endif
      endif
      write(stdo,*) 'end of reading rst or atm'
!! (we need check. a simple approximaiton to determine VBM and CBM. Need fixing if necessary).
      if(cmdopt0('--vbmonly')) then !Get VBM and CBM relative to vaccum
         if(master_mpi) call vbmmode()
         call rx0('--vbmonly mode done')
      endif
      if(cmdopt0('--getq')) then ! Current version is not for spin dependent, with many restrictions.
         if(master_mpi) call getqmode()
         call rx0('--getq mode done')
      endif
!!      
      if(iprint()>=30) call praugm(sspec,0) !print out properties of species
      if(lrout == 0 .and. ctrl_lfrce /= 0) then
        write(stdo,"(a)") 'Error: output density required when forces sought.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call rx('incompatible input. see the end of console output')
      endif
      if(lrout == 0 .and. cmdopt0('--etot')) then
        write(stdo,"(a)") 'Error: output density required with --etot switch.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call rx('incompatible input. see the end of console output')
      endif
!! Main rouitne. After module initializations, call lmfp  -----------------------
      call Lmfp(llmfgw,irs1x)
!! finalization      
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
      call m_MPItk_finalize()
      if(master_mpi) write(6,"(a)") "OK! end of "//trim(prgnam)//" ======================"
      call Exit(0)
      end program Lmf
