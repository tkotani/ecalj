!! For lmchk. Single core mode only. No MPI modes,although our makefile use mpif90 because of simplicity.
      program lmchk
c      use m_work,only: iv_a_oipqx
      use m_mksym,only: iv_a_oipc,iv_a_oics,iv_a_onrc,rv_a_oclabl,iv_a_onrcp
      use m_clsprm,only: clsprp,clsprm
      use m_struc_def
c      use m_globalvariables
      use m_ext,only: getext,sname
      use m_lmfinit, only: lso,
     &     slabl_,v_sspec,v_ssite,!v_sham,!,v_slat,v_slat2 v_spot,v_sarry,v_sarry2,v_smix,
     &     v_ssite2,sstrnmix,sstrnsymg, !v_sbz,,v_sctrlv_sctrl2,
     &     ctrl_nl,ctrl_nspin,ctrl_nclass,ctrl_nl,ctrl_nspec,ctrl_nbas,ctrl_quit,rdctrl2
      use m_clsprm,only: rv_a_opnu,
     &     rv_a_obxc , rv_a_ogrrme ,  rv_a_opmpol , rv_a_opnu , rv_a_opp , rv_a_opprel,
     &     rv_a_oqnu , rv_a_oqpp , rv_a_oqt , rv_a_orhrmx ,
     &     rv_a_osop , rv_a_ovdif , rv_a_oves , rv_a_ovintr , rv_a_ovrmax
      use m_supot,only:rv_a_omad
      implicit none
      character prgnam*8, vrsion(2)*6, ext*100
      logical lgors
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     .     igets,a2vec,NULLI,mxspec
      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)
      data prgnam /'LMCHK'/ auxmod /1/
      parameter (mxrecs=10000,recln0=511)
      character*8 alabl
      character*(mxrecs*recln0):: recrd
      real(8),parameter::    NULLR =-99999
      real(8)::pwmode,pwemin,pwemax
      integer::ic,ibas
      real(8):: bz_w,dval       !test
#if MPI | MPIK
      include "mpif.h"
      integer fext
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
c      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog
      integer  nsize, id
!! =================================================================
      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      if(id==0) write(*,*) 'mpi size=',nsize
#endif
!!
      if(nsize/=1) call rx('lmchk allows only -np 1!')
C     --- Version ---
      call show_programinfo(6)
      stdo = lgunit(1)
c      globalvariables%stdo = dble(stdo)
c      globalvariables%l_stdo = globalvariables%l_stdo +1
      if (cmdopt('--h',3,0,outs)) call lmhelp(prgnam)
      call headl2(prgnam,0,stdo)
      i = dble(stdo)
      call finits(2,0,0,i)
c      call getext()
c      ext=adjustl(sname)
     
      call poseof(fopn('LOG'))
!! Get extension of ctrl file      
c      call finits(2,0,0,i)
c      call getext(master==procid)
c      print *,' trim sname=',trim(sname)
c      ext=adjustl(sname)
c      i = fextg(ext)
c      call word(ext,1,i,j)
c      if (ext(i:i) .eq. '.') i=i+1
c      if (ext(j:j) .eq. '.') j=j-1
c      if (j .ge. i) call ptenv('EXT='//ext(i:j))
!!
      
C     ... Abort with error message if ctrl file is missing (swtmp = .true.)
      swtmp = .false.
c     if (cmdopt('--input',6,0,strn)) then
c     if (nproc .gt. 0) call rx('--input not allowed with MPI')
c     else
      if (fxst('CTRL') .ne. 1) then
         call awrit0(' '//prgnam//'%a:%9pmissing ctrl file',' ',80,
     .        i1mach(2))
         swtmp = .true.
      endif
C     ... File logical units
      i = fadd('TMP',-1,4)
      i = fadd('BAND',-1,4)
#if LMDOS
      i = fadd('MOMS',-1,4)
      i = fadd('CLS',-1,4)
#endif
C     --- Set the top-level verbosity if specified from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
         i = 4
         i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
         if (i .eq. 1) call setpr(k)
      endif
C     --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
C     ki      allocate( recrd( 0:mxrecs*recln0-1 ) )
c     if (procid .eq. master) then
      if (.not.cmdopt('--input',7,0,strn)) then
         nfilin = fopna('CTRL',-1,1)
c     takao
         call findctrlstart(nfilin)
         alabl = '#{}% ct '
         if (cmdopt('--show',6,0,strn)) alabl = '#{}% ctp'
         call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs)
c     To show variables.
c     call shosyv(0,-1,0,6)
         i = 60
         if (cmdopt('--show',6,0,strn)) then
            i = 1
            call info0(i,0,0,' ------------------------ End '//
     .           'of input file ----------------------')
         endif
         call info2(i,0,1,' '//prgnam//'%a : %i lines read from'//
     .        ' input file',nrecs,0)
         if (cmdopt('--showp',7,0,strn)) call cexit(0,1)
      endif
C     write(stdo,*)' === TEST out recrd === '
C     j = 0
C     do i = 1, nrecs
C     write(stdo,"(i2,1x,256a1)") procid,recrd(j:j+recln0-1)
C     j = j + recln0
C     enddo
C     write(stdo,*)' --- end of rdfile --- '
C     call rx('done')
      call rdctrl2(recrd,recln0,nrecs,prgnam)
c$$$c     .     slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite, !v_sstr,
c$$$c     .     v_sarry,sstrnmix,sstrnsymg) 
c$$$C     --- Lattice setup ---
c$$$c      call setcg(v_slat,8,12)
c$$$c      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
c$$$      call lattic()!v_slat,ctrl_nbas)
c$$$      do ibas=1,ctrl_nbas
c$$$        v_ssite(ibas)%pos = rv_a_opos(3*(ibas-1)+1:3*(ibas-1)+3)
c$$$c        call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos
c$$$c     .   , i_copy_size , i_spackv + 1 - 1 , slat%rv_a_opos )
c$$$      enddo
c$$$C     --- Generate symmetry operations; split species into classes  ---
c$$$      strn = 'find'
c$$$      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg) 
c$$$      print *,' strn sym=',trim(strn)
c$$$c      call upacks('strn symg',i,j)
c$$$c      if (j .ge. i)strn = sstrn(i:j+1)
c$$$      if (cmdopt('--nosym',7,0,outs)) strn = ' '
c$$$      lc = 20
c$$$cccccc
c$$$c      if (  iand(1,int(ctrl_lqp)) == 0 ) lc = lc+2 !iand(2**n,i) means 2**n ditit of i.
c$$$c      if ( .not. iand(1,int(ctrl_lqp)) .ne.0 ) lc = lc+2
c$$$      if( ctrl_noinv==0 ) lc = lc+2 !no inversion
c$$$ccccccc
c$$$  call mksym(lc,slabl_,strn,v_sctrl,v_slat,v_ssite,v_sarry)
      
c      v_sarry%nclasp=ctrl_nclass
C     --- Allocate permanent class arrays, maps and other initialization ---
      call clsprm(1,v_sspec)!, v_slat) !,v_sbz) !,v_sstr)!,sstrn)v_sctrl,v_sham,v_spot,,v_sarry
C     --- Read available class parameters from file ---
c#if LMCHK
      call aiocls(.false.,0,!v_sarry,!v_spot, !v_ctrl v_sham,
     .     v_sspec,1,0)!v_slat
      call aiocls(.false.,17,!v_sarry,!v_spot, !v_sham,
     .     v_sspec,1,0)!,v_slat
      if (lc .ge. 20) call clsprp(1,v_sspec)!, v_slat)!,v_sbz)! ,v_sstr)!,sstrn)v_sham,v_spot,,v_sarry
c     #endif
      if (cmdopt('--sfill',7,0,strn)) call rx('no --sfill option in lm7K')
C     ... Set default values for species data
      call defspc(v_sspec)
C     ... quit after SHOW
      if ( int(ctrl_quit) .eq. 1 ) then 
         call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
         call rx0(prgnam)
      endif
c#if LMCHK & FP
!! takao did not check --fp mode recently
      if (cmdopt('--fp',4,0,strn)) then
         call fpchk(v_sspec,v_ssite)
         call cexit(0,1)
      endif
c#endif
c#if LMCHK
      call lmaux(prgnam,v_sspec, !v_sham,,v_slat v_spot,
     .     v_ssite,slabl_,auxmod,0) !wksize),sstrn, v_sstr,v_sbz,v_sctrl,v_sarry,
c#endif
c#if LMDOS
c      call asados(prgnam,v_sbz,v_sctrl,v_sham,v_slat, v_sspec,v_ssite,v_sarry)
c#endif
 1000 continue
      call rx0(prgnam)
      end program lmchk
