      program lmfa
      use m_struc_def
c      use m_rdctrl2_func, only: rdctrl2
      use m_lmfinit,only: globalvariables
c      use m_globalvariables
      use m_lmfinit, only: !rdctrl2,,v_sctrl,v_sctrl2 v_sham,
     &     v_sspec,!,v_ssite,  !v_sarry2,v_sbz ,v_slat,v_slat2 ,v_spot slabl_,v_smix,
     &     rdctrl2 !v_ssite2,sstrnmix,sstrnsymg,
      use m_work,only: rv_a_opos
      implicit none
      integer procid, master, mpipid, nproc
c      integer,parameter:: lstrn=10000
c      character(lstrn):: sstrnmix,sstrnsymg='' !,jobid !header,
      character prgnam*8, vrsion(2)*6, ext*100
c      logical lgors
c      integer mxspec
c      type(s_bz):: v_sbz
c      type(s_ctrl):: v_sctrl
c      type(s_lat):: v_slat
c      type(s_ham):: v_sham
c      type(s_pot):: v_spot
c      type(s_mix):: v_smix
c      type(s_array):: v_sarry
c      type(s_spec),allocatable:: v_sspec(:)
c      type(s_site),allocatable:: v_ssite(:)
c      parameter (mxspec=256)
c      character*8 slabl(mxspec)
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     *igets,a2vec,NULLI
      double precision dglob !
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)
      data prgnam /'LMFA'/
      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      parameter (mxrecs=10000,recln0=511)
      character*8 alabl
c      character:: recrd_(0:mxrecs*recln0-1)
      character*(recln0):: recrd(mxrecs)
c      equivalence ( recrd, recrd_ )
      real(8),parameter::    NULLR =-99999
      integer:: ibas

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
      if(nsize/=1) call rx('lmfav7.F allows only -np 1!')
      call show_programinfo(6)
      stdo = lgunit(1)
c     i = dglob('stdo',dble(stdo),1)
      globalvariables%stdo = dble(stdo)
     
      master = 0
      procid = mpipid(1)
      nproc  = mpipid(0)

      if (cmdopt('--h',3,0,outs)) call lmhelp(prgnam)!,vn(2),vn2(2),0)

      call headl2(prgnam,0,stdo)
      i = dble(stdo)
      call finits(2,0,0,i)

      call poseof(fopn('LOG'))
c$$$      i = fextg(ext)
c$$$      call word(ext,1,i,j)
c$$$      if (ext(i:i) .eq. '.') i=i+1
c$$$      if (ext(j:j) .eq. '.') j=j-1
c$$$      if (j .ge. i) call ptenv('EXT='//ext(i:j))

C     ... Abort with error message if ctrl file is missing (swtmp = .true.)
      swtmp = .false.
      if (cmdopt('--input',6,0,strn)) then
        if (nproc .gt. 0) call rx('--input not allowed with MPI')
      else
        if (fxst('CTRL') .ne. 1) then
          call awrit0(' '//prgnam//'%a:%9pmissing ctrl file',' ',80,
     *    i1mach(2))
          swtmp = .true.
        endif
      endif
      call mpibc1(swtmp,1,1,.false.,'lmf','error')
      if (swtmp) call cexit(-1,1)
C --- Set the top-level verbosity if spec'd from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i .eq. 1) call setpr(k)
      endif
C --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
      if (procid .eq. master) then
        if (.not.cmdopt('--input',7,0,strn)) then
          nfilin = fopna('CTRL',-1,1)
          call findctrlstart(nfilin)
          alabl = '#{}% ct '
          if (cmdopt('--show',6,0,strn)) alabl = '#{}% ctp'
          call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs)
          i = 60
          if (cmdopt('--show',6,0,strn)) then
            i = 1
            call info0(i,0,0,' ------------------------ End '//
     *      'of input file ----------------------')
          endif
          call info2(i,0,1,' '//prgnam//'%a : %i lines read from'//
     *    ' input file',nrecs,0)
          if (cmdopt('--showp',7,0,strn)) call cexit(0,1)
        endif
      endif
C --- Read recrd ---
      call rdctrl2(recrd,recln0,nrecs,prgnam)
c$$$c     *slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite, !v_sstr,
c$$$c     *v_sarry,sstrnmix,sstrnsymg) 
c$$$      print *,'end of rdctrl2 in imfav7'
c$$$c      call setcg(v_slat,8,12)
c$$$c      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
c$$$      call lattic()!v_slat,v_sctrl%nbas)
c$$$      do ibas=1,v_sctrl%nbas
c$$$        v_ssite(ibas)%pos = rv_a_opos(3*(ibas-1)+1:3*(ibas-1)+3)
c$$$c        call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos
c$$$c     .   , i_copy_size , i_spackv + 1 - 1 , slat%rv_a_opos )
c$$$      enddo
c$$$C --- Generate symmetry operations; split species into classes  ---
c$$$      strn = 'find'
c$$$      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
c$$$      if (cmdopt('--nosym',7,0,outs)) strn = ' '
c$$$      lc = 20
c$$$c      if ( .not. iand(1,int(v_sctrl%lqp)) .ne.0 ) lc = lc+2
c$$$      if ( .not. iand(1,int(v_sctrl%noinv)) .ne.0 ) lc = lc+2
c$$$      print *,'goto mksym:',trim(strn)
c$$$      call mksym(lc,slabl_,strn,v_sctrl,v_slat,v_ssite)!,v_sarry)
c$$$      print *,'end of mksym x'
c$$$C ... Set default values for species data
c$$$      print *,'goto defspc'
c$$$      call defspc(v_sspec)
c$$$      print *,'end of defspc'
c$$$      if ( int(v_sctrl%quit) .eq. 1 ) then
c$$$        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
c$$$        call rx0(prgnam)
c$$$      endif
c      print *,'goto freeat uuuuuuuuu'
      call freeat(v_sspec)!,v_spot)!v_sctrl,v_sham,
      call rx0(prgnam)
      end

