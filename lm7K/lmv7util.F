!     ! For lmdos and lmchk. Single core mode only. No MPI modes.
      program lmfutil
      use m_struc_def
      use m_rdctrl2_func, only: rdctrl2
      use m_globalvariables
      use m_rdctrl,only:lso
      implicit none
C     ... Controls for IO
c      integer lstrn
c      parameter (lstrn=20000)
      character prgnam*8, vrsion(2)*6, ext*100
C     ... Structure for strings
c      character*(lstrn) sstrn
C     ... For structures
      logical lgors
      integer mxspec
c      type(s_sstrn)::v_sstrn !may2015
      integer,parameter:: lstrn=10000
      character(lstrn):: sstrnmix,sstrnsymg !,jobid !header,

      type(s_bz):: v_sbz
      type(s_ctrl):: v_sctrl
      type(s_lat):: v_slat
      type(s_ham):: v_sham
      type(s_pot):: v_spot
      type(s_mix):: v_smix
      type(s_array):: v_sarry
      type(s_move):: v_smove
      type(s_str):: v_sstr
c      type(s_tb):: v_stb
      type(s_spec),allocatable:: v_sspec(:) 
      type(s_site),allocatable:: v_ssite(:) 
      parameter (mxspec=256)
      character*8 slabl(mxspec)
C     ... miscellaneous local variables
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     .     igets,a2vec,NULLI
      double precision dglob !vn(2),vn2(2),
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)
C     ... Program-dependent name and help
#if  LMDOS
      data prgnam /'LMDOS'/
#endif
#if LMCHK
      data prgnam /'LMCHK'/ auxmod /1/
#endif
      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      parameter (mxrecs=10000,recln0=511)
      character*8 alabl
      character*(mxrecs*recln0):: recrd
      real(8),parameter::    NULLR =-99999
      real(8)::pwmode,pwemin,pwemax
      integer::ic
      real(8):: bz_w,dval       !test

C     --- Version ---
      call show_programinfo(6)
c$$$      vrsion(1) = 'LM'
c$$$      vrsion(2) = ' '
c$$$      vn(1) = 7.00d0
c$$$      vn2(1) = 0
c$$$      vn(2) = 0.00d0
c$$$      vn2(2) = 0
c$$$      vn(1) = vn(1) + vn2(1)/1d5
c$$$      vn(2) = vn(2) + vn2(2)/1d5
      stdo = lgunit(1)
      globalvariables%stdo = dble(stdo)
      globalvariables%l_stdo = globalvariables%l_stdo +1
      i = dble(stdo)
C     --- Help ---
      swtmp = .false.
      if (swtmp .or. cmdopt('--h',3,0,outs)) call lmhelp(prgnam)!,vn(2),vn2(2),0)
c$$$      if (cmdopt('--version',9,0,outs)) then
c$$$         write(*,'(f4.2)') vn(2)
c$$$         call cexit(0,1)
c$$$      endif
      call headl2(prgnam,0,stdo)
      call finits(2,0,0,i)
      call poseof(fopn('LOG'))
      i = fextg(ext)
      call word(ext,1,i,j)
      if (ext(i:i) .eq. '.') i=i+1
      if (ext(j:j) .eq. '.') j=j-1
      if (j .ge. i) call ptenv('EXT='//ext(i:j))
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
      call rdctrl2(recrd,recln0,nrecs,prgnam,F, !,vrsion,vn,vn2
     .     slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     .     v_sarry,v_smove,sstrnmix,sstrnsymg) !,v_stb
C     --- Lattice setup ---
      call setcg(v_slat,8,12)
      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
C     --- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg) 
      print *,' strn sym=',trim(strn)
c      call upacks('strn symg',i,j)
c      if (j .ge. i)strn = sstrn(i:j+1)
      if (cmdopt('--nosym',7,0,outs)) strn = ' '
      lc = 20
      if ( v_sctrl%lqp==0 ) lc = lc+2 
      call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
C     --- Allocate permanent class arrays, maps and other initialization ---
      call clsprm(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec, v_slat,v_sbz,v_sstr)!,sstrn)
C     --- Read available class parameters from file ---
#if LMCHK
      call aiocls(.false.,0,v_sarry,v_sctrl,v_sham,v_spot,
     .     v_sspec,v_slat,1,0)
      call aiocls(.false.,17,v_sarry,v_sctrl,v_sham,v_spot,
     .     v_sspec,v_slat,1,0)
      if (lc .ge. 20)
     .     call clsprp(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .     v_slat,v_sbz,v_sstr)!,sstrn)
#endif
      if (cmdopt('--sfill',7,0,strn)) then
         call rx('no --sfill option in lm7K')
      endif
C     ... Set default values for species data
      call defspc(v_sspec)
C     ... quit after SHOW
      if ( int(v_sctrl%quit) .eq. 1 ) then 
         call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
         call rx0(prgnam)
      endif
#if LMCHK & FP
      if (cmdopt('--fp',4,0,strn)) then
         call fpchk(v_sspec,v_ssite)
         call cexit(0,1)
      endif
#endif
#if LMCHK
      call lmaux(prgnam,v_sctrl,v_sham,v_spot,v_slat,v_sspec,
     .     v_ssite,v_sstr,v_sarry,v_sbz,slabl,auxmod,0) !wksize),sstrn
#endif
#if LMDOS
      call asados(prgnam,v_sbz,v_sctrl,v_sham,v_slat, v_sspec,v_ssite,v_sarry)
#endif
 1000 continue
      call rx0(prgnam)
      end program lmfutil


      subroutine lmhelp(prgnam)!,vn,vn2,wksize)
C     - Help printout
C     ----------------------------------------------------------------------
C     i Inputs
C     i   prgnam:name of main program
C     i   vn    :program main version
C     i   vn2   :program minor version
C     i   wksize:worksize
C     o Outputs
C     o   message written to stdout
C     r Remarks
C     u Updates
C     u   11 Apr 03
C     ----------------------------------------------------------------------
      implicit none
C     ... Passed parameters
      character prgnam*8
c      double precision vn,vn2
c      integer wksize
C     ... Local parameters
      integer i1,i2
      character ch*1,outs*1000

      call locase(prgnam)
      call info0(0,0,0,' usage:  '//prgnam//
     .     '%a [--OPTION] [-var-assign] [ext]')

      print 343
      print 344
 343  format(/' --h'/' --help',t17,'Print this message, and quit'
     .     /' --input',t17,
     .     'List categories, tokens, and data program expects, and quit'
     .     /' --show',t17,
     .     'Print control file after parsing by preprocessor,'/t17,
     .     'and echo input data as read from the control file'
     .     /' --showp',t17,
     .     'Same as --show, but quit after input parsed'
     .     /' --iactiv',t17,'(--no-iactiv) ',
     .     'Turn on (off) interactive mode'/
     .     t17,'This switch overrides input file setting',
     .     /' --pr#1[,#2...]',t17,
     .     'Set the verbosity (stack) to values #1,#2, ...'
     .     /' --time=#1[,#2]',t17,
     .     'Print timing info to # levels (#1=summary; #2=on-the-fly)'/
     .     /' -vnam=expr',t17,
     .     'Define numerical variable "nam"; set to result of ''expr'''
     .     /' -cnam=strn',t17,
     .     'Define character variable "nam"; set to ''strn'''
     .     )


 344  format(
     .     /' --rpos=filnam',t17,
     .     'After reading input file, read site positions from "filnam"'/
     .     ' --fixpos[:tol=#]',
     .     ' Adjust positions slightly, rendering them as'/
     .     t17,'exactly as possible consistent with the symmetry group')

      if (.true.) then
         outs = '%N '//prgnam//'%a-specific options:'
         call strip(outs,i1,i2)
         call info0(0,0,0,outs(1:i2))

         if (prgnam .eq. 'lmdos') then
            call info0(0,0,0,
     .           '%N%1f--dos:options  modifies number and kinds of dos '//
     .           'generated;%N%16fsee documentation (doc/lmto.html)')
         endif

         if (prgnam .eq. 'lmchk') then
            call info0(0,0,0,
     .           '%N%1f--shell[:v][:e][:r=#][:sites:site-list]'//
     .           '[:pairs:pair-list]...'//
     .           '%N%8f...[:tab[=#]][:disp=fnam][:nn][:fn=fnam]'//
     .           '%N --mino[:dxmx=#][:xtol=#][:maxit=#][:style=#]:list'//
     .           '%N --findes'//
     .           '%N --wpos=fnam'//'%N --angles[=#]'//'%N --terse')
         endif

      endif

c$$$      if (nint(vn2) .ne. 0) then
c$$$         ch = char(nint(vn2)+ichar('a')-1)
c$$$      else
c$$$         ch = ' '
c$$$      endif
c$$$      call info2(0,0,0,'%N version %,2d'//ch//'   worksize %i K',
c$$$     .     vn-vn2/1d5,wksize/1000)
      call fexit(0,0,' ',0)
      end subroutine lmhelp 


      subroutine findctrlstart(nfilin)
      use m_globalvariables
      character(len=9):: strn, ccc
      do
         read(nfilin,"(a)",err=1010,end=1010) strn
         if(strn == 'ctrlstart') return
      enddo
 1010 continue
      rewind(nfilin)
C     getarg       end
      end subroutine findctrlstart 

      subroutine clsprm(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,sstr) !,sstrn)
      use m_struc_def           !Cgetarg
C     - Allocate memory for, and try to read, class-related data
C     ----------------------------------------------------------------------
C     i    mode :1s digit
C     i         : 1 allocates ASA arrays (see spot,sham,sarray)
C     i         : 2 allocates TBE arrays:  oifrlx
C     i   sctrl :struct for program flow parameters; see routine uctrl
C     i     Elts read: nbas nclass nl nspin nspec loptc nmap
C     i     Stored:    *
C     i     Passed to: lgors ioqpp aiocls clsprp
C     i   sarray:structure containing offsets to various arrays
C     i     Elts read: nclasp oics oipc oclabl opos
C     i     Stored:    ogroup ormax ohave oifrlx ozos
C     i     Passed to: aiocls clsprp
C     i   sham  :struct for parameters defining hamiltonian; see routine uham
C     i     Elts read: *
C     i     Stored:    osop oivso
C     i     Passed to: aiocls clsprp
C     i   spot  :struct for information about the potential; see routine upot
C     i     Elts read: ovrmax orhrmx opnu oqnu opp oves
C     i     Stored:    opnu oqnu oqc oqt oves ovdif ovrmax orhrmx opp orhos
C     i                osop ogrrme ovintr opmpol oqpp
C     i     Passed to: ioqpp aiocls clsprp
C     i   sspec :struct for species-specific information; see routine uspec
C     i     Elts read: rmt idmod z lmxa
C     i     Stored:    *
C     i     Passed to: aiocls clsprp
C     i   slat  :struct for lattice information; see routine ulat
C     i     Elts read: plat
C     i     Stored:    *
C     i     Passed to: aiocls clsprp
C     i   sbz   :struct for the Brillouin Zone; see routine ubz
C     i     Elts read: ndos
C     i     Stored:    *
C     i     Passed to: clsprp
C     i
C     i   sstr  :struct for parameters for screened strux; see routine ustr
C     i     Elts read: n map
C     i     Stored:    *
C     i     Passed to: *
C     i   sstrn :struct for global strings
C     i     Elts read: map
C     i     Stored:    *
C     i     Passed to: clsprp
C     u Updates
C     u   08 Nov 07 (J. Xu) LDA+U implementation; qpp is complex
C     u   30 Sep 04 Reads/writes relativistic ppar's
C     u   02 Apr 04 SO array redimensioned; make and pack pot->nrhos
C     u   18 Mar 03 Allocate space for relativistic potential parameters
C     u   07 Feb 03 SO array redimensioned
C     u   19 Feb 02 entry clsprp copies class data to doubly padded layers
C     u   28 Apr 98 allocate space for radial matrix elements of grad
C     ----------------------------------------------------------------------
C     implicit none
      integer mode
c      character*(*) sstrn
      type(s_ctrl)::sctrl
      type(s_array)::sarray
      type(s_ham)::sham
      type(s_pot)::spot
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_bz)::sbz
      type(s_str)::sstr
C     Local variables
      logical lgors,cmdopt
      character*80 strn,clabl*8
      integer:: ogroup , oifrlx , oivso , oipc , opos , ozos , n0 
      parameter (n0=10)
      integer nbas,nclass,nspec,nl,nsp,nlspc,igets,nclasp,i1,i2,nmap,ic,
     .     is,ib,iclbsj,idmod(n0),iprint,lmx,nclspp,nqpp,lnsph,isw,neula,
     .     nbf,nrhos,nlibu
      double precision plat(9),z,xx,dval
      integer,allocatable:: ics(:)

      nbas=sctrl%nbas
      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec
      nlibu=sham%nlibu
      nclasp=sarray%nclasp
      nclspp = max(2*nclasp-nclass,nspec)
      nlspc = nl*nsp*max(nclspp,nspec)

C     ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then
         allocate(spot%rv_a_opnu(nlspc))
         spot%rv_a_opnu(:)=0.0d0
         allocate(spot%rv_a_oqnu(3*nlspc))
         spot%rv_a_oqnu(:)=0.0d0
         allocate(spot%rv_a_opp(6*nlspc))
         spot%rv_a_opp(:)=0.0d0
         allocate(spot%rv_a_opprel(32*nl*nl*max(nclspp,nspec)))
         spot%rv_a_opprel(:)=0.0d0
         allocate(spot%rv_a_oqt(nclasp))
         spot%rv_a_oqt(:)=0.0d0
         allocate(sarray%iv_a_oipq(nclasp))
         sarray%iv_a_oipq(:)=0
         allocate(spot%rv_a_oves(nclspp))
         spot%rv_a_oves(:)=0.0d0
         allocate(spot%rv_a_ovdif(nclspp))
         spot%rv_a_ovdif(:)=0.0d0
         allocate(spot%rv_a_orhrmx(nclspp))
         spot%rv_a_orhrmx(:)=0.0d0
         allocate(spot%rv_a_ovrmax(2*nclspp))
         spot%rv_a_ovrmax(:)=0.0d0
         allocate( ics(nclspp) )
         call icopy ( nclspp , sarray%iv_a_oics , 1 , ics , 1 ) 
         print *,' qqqqq ics nclspp=',ics,nclspp
         allocate(sarray%rv_a_ormax(abs(nclspp)))
         call dcopy ( nclspp , sspec ( ics ( 1:nclspp ) ) %rmt , 1 , sarray%rv_a_ormax, 1 ) 
C     Even in collinear case, bxc might be written to
         allocate(spot%rv_a_obxc(abs(-3*nclasp)))
         spot%rv_a_obxc(:)=0.0d0
         spot%nrhos=0 
c         if ( sctrl%lncol /=0 ) then 
         if ( lso /=0 ) then 
            neula=sham%neula
            nbf=sham%nbf
            nrhos = max(nl,neula,nbf)
            if (allocated(spot%rv_a_obxc)) deallocate(spot%rv_a_obxc)
            allocate(spot%rv_a_osop(nlspc*nsp*9))
            spot%rv_a_osop(:)=0.0d0
            allocate(spot%rv_a_obxc(3*nclasp))
            spot%rv_a_obxc(:)=0.0d0
            spot%nrhos=nrhos 
         endif
         if (allocated(spot%rv_a_ogrrme)) deallocate(spot%rv_a_ogrrme)
c$$$         if ( int ( sctrl%loptc ) .gt. 0 ) then 
c$$$            allocate(spot%rv_a_ogrrme(abs(-16*nl*nsp*nclasp)))
c$$$            if (-16*nl*nsp*nclasp<0) spot%rv_a_ogrrme(:)=0.0d0
c$$$         endif 
         if ( iand(3,int(sctrl%lsx)) .ne.0 .or. iand(-1,int(sctrl%lscr)) .ne.0 ) then 
            allocate(spot%rv_a_ovintr(abs(-nclasp*(nl*nsp)**2)))
            spot%rv_a_ovintr(:)= 0d0
         endif
         lnsph = isw ( .false.)
C     ... Allocate memory for, attempt to read qpp from file
         if (lnsph .ne. 0 .or. nlibu .ne. 0) then
            allocate(spot%rv_a_opmpol(abs(-(2*nl-1)*nl**2*3*nsp*nclasp)))
            spot%rv_a_opmpol(:)= 0d0
            i1 = nl**2
            nqpp = (i1*(i1+1))/2
            allocate(spot%rv_a_oqpp(2*nqpp*4*nsp*nbas))
            call dvset ( spot%rv_a_oqpp , 1 , 1 , - 1d0 ) 
            call ioqpp(.false.,sctrl,spot)
         endif
      endif
      return

      entry clsprp(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,
     .     sstr) !,sstrn)
C     - Copy class data to doubly padded case
      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec
      nclasp=sarray%nclasp

      if (nclasp .eq. nclass) return
      nclspp = max(2*nclasp-nclass,nspec)
C     ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then
         i1 = nclasp-nclass
         call dpscop ( spot%rv_a_opnu , spot%rv_a_opnu , nl * nsp * i1 , nl * nsp 
     .        * nclass + 1 , nl * nsp * nclasp + 1 , 1d0 ) 
         call dpscop ( spot%rv_a_oqnu , spot%rv_a_oqnu , 3 * nl * nsp * i1 , 3 * 
     .        nl * nsp * nclass + 1 , 3 * nl * nsp * nclasp + 1 , 1d0 ) 
         call dpscop ( spot%rv_a_oves , spot%rv_a_oves , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
         call dpscop ( spot%rv_a_ovdif , spot%rv_a_ovdif , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
         call dpscop ( spot%rv_a_ovrmax , spot%rv_a_ovrmax , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
         call dpscop ( spot%rv_a_orhrmx , spot%rv_a_orhrmx , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
      else
         call rxi('clsprp: bad mode',mode)
      endif
      end subroutine clsprm 

C================================================================
      subroutine ioqpp(lio,sctrl,spot)
      use m_struc_def           !Cgetarg
C     - File i/o for phi-phi, phi-dot, dot-dot products
C     ----------------------------------------------------------------------
C     i Inputs:
C     i   lio: true for write, false for read
C     i          <0 write
C     u Updates
C     u   08 Nov 07 (J. Xu) qpp is complex
C     ----------------------------------------------------------------------
C     implicit none
C     ... Passed parameters
      logical lio
C     getarg       double precision sctrl(1),spot(1)
      type(s_ctrl)::sctrl
      type(s_pot)::spot
C     ... Local parameters
      integer:: i1 , i2 , nbas , nl , nsp , igets , ifi , fopna , rdm 
     .     , ipr 
      call getpr(ipr)
      nbas=sctrl%nbas
      nl=sctrl%nl
      nsp=sctrl%nspin
      i1 = nl**2*(nl**2+1)
      i2 = 4*nsp*nbas
      ifi = fopna('qpp',-1,4+8)
      if (lio) then
         call ywrm ( 1 , 'lmasa' , 1 , ifi , ' ' , spot%rv_a_oqpp , 1 , i1 
     .        , i1 , i2 ) 
         if (ipr .ge. 30) print *, 'IOQPP:  wrote qpp to disk'
      else
         if (ipr .ge. 30) print *, 'IOQPP:  reading qpp from disk ...'
         call pshpr(0)
         if ( rdm ( ifi , 2 , i1 * i2 , ' ' , spot%rv_a_oqpp , i1 , i2 ) .lt. 
     .        0 ) then 
            if (ipr .ge. 0) print *,'IOQPP:  (warning) failed to qpp file'
            call dvset ( spot%rv_a_oqpp , 1 , 1 , - 1d0 ) 
         endif
         call poppr
      endif
      call fclose(ifi)
      end subroutine ioqpp 


C================================================================
c     take subs/spinav.F here
      subroutine spinav(mode,nclass,nl,nsp,pnu,qnu)
C     - Averages up+down spin moments + pp's for all classes
C     ----------------------------------------------------------------------
C     i Inputs
C     i   mode  :0 average spins
C     i         :1 do not average, but exchange spins
C     i   nclass:number of inequivalent classes
C     i   nl    :(global maximum l) + 1
C     i   nsp   :2 for spin-polarized case, otherwise 1
C     i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
C     i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
C     i   qnu   :energy-weighted moments of the sphere charges
C     o Outputs :moments are spin-averaged
C     i   pnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
C     i   qnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
C     o   nsp   :set to 1 on output (mode=0)
C     l Local variables
C     l         :
C     r Remarks
C     r   
C     u Updates
C     u   10 Jan 06 Added mode
C     ----------------------------------------------------------------------
C     implicit none
      integer mode,nclass,nl,nsp,ic
      double precision pnu(nl,nsp,nclass),qnu(3,nl,nsp,nclass)

      if (nsp .eq. 1) return

      if (mode .gt. 1) then
         call rx('spinav: bad mode')
      elseif (mode .eq. 1) then
         do  ic = 1, nclass
            call dswap(nl,pnu(1,2,ic),1,pnu(1,1,ic),1)
            call dswap(3*nl,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
         enddo
         return
      endif

      do  ic = 1, nclass
         call daxpy(nl,1d0,pnu(1,2,ic),1,pnu(1,1,ic),1)
         call dscal(nl,.5d0,pnu(1,1,ic),1)
         call daxpy(3*nl,1d0,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
      enddo

      do  ic = 2, nclass
         call dcopy(nl,pnu(1,1,ic),1,pnu(1,ic,1),1)
         call dcopy(3*nl,qnu(1,1,1,ic),1,qnu(1,1,ic,1),1)
      enddo
      nsp = 1
      end
c     #define ONE_ATOM_UNIT 1
      subroutine aiocls(lio,mode,sarray,sctrl,sham,spot,sspec,slat,ic1,ic2)
      use m_struc_def           !Cgetarg
      use m_globalvariables
C     - File I/O atomic data for classes ic1..ic2
C     ----------------------------------------------------------------------
C     i Inputs
C     i   lio    F for read, T for write
C     i   mode   ones digit
C     i          1 use data from first class corresponding to same species,
C     i            if data from own class is missing
C     i          2 like 1, but attempt to read data from disk anyway
C     i          3 make no attempt to read from disk but copy info from
C     i            first class corresponding to same species, if data missing
C     i          4 Add 4 if to use default P,Q when not otherwise supplied
C     i          10s digit
C     i          1 assemble background rho
C     i   sarray:structure containing offsets to various arrays
C     i     Elts read: nclasp oclabl ohave oics onrc
C     i     Stored:
C     i     Passed to:
C     i   sctrl :struct for program flow parameters; see routine uctrl
C     i     Elts read: nl nspin zbak nbas nclass
C     i     Stored:    zbak
C     i     Passed to: lgors
C     i   sham  :struct for parameters defining hamiltonian; see routine uham
C     i          (not used now)
C     i   spot  :struct for information about the potential; see routine upot
C     i     Elts read: opp osop ogrrme opnu oqnu orhrmx ovrmax ovintr opmpol
C     i                oves
C     i     Stored:
C     i     Passed to:
C     i   sspec :struct for species-specific information; see routine uspec
C     i     Elts read: lmxa idmod p q
C     i     Stored:
C     i     Passed to:
C     i   slat  :struct for lattice information; see routine ulat
C     i     Elts read: vol
C     i     Stored:
C     i     Passed to:
C     i   ic1,ic2: range of classes to read data
C     r Remarks
C     r   Right now, read always takes data from file if available
C     u Updates
C     u   09 Nov 07 Corrected sign of default moment (paioc2)
C     u   29 Sep 04 Reads/writes relativistic ppar's
C     u   26 Apr 03 Added MPI calls
C     u   07 Feb 03 adjusted for redimensioned sop
C     u   30 May 02 Assign better default P
C     u   28 Apr 98 I/O of radial matrix elements of grad
C     u   28 Sep 00 Added setting default P,Q
C     ----------------------------------------------------------------------
C     implicit none
C     Passed parameters
      logical lio
      integer mode,ic1,ic2
      type(s_array)::sarray
      type(s_ctrl)::sctrl
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_ham)::sham
      type(s_pot)::spot
C     Local variables
      character*8 clabl,alabel, outs1*20, outs2*20, outs3*20, outs*80
      logical sw,aiomom,aiopar,aiopot,aiova,lpot,lcor,scat,lgors,
     .     aiogen,aiosop,aiorme,aiocor,aiomp,lrell,lgen
      logical havepq,havepp,haveso,haveop,haveva,havemp,
     .     readpq,readpp,readso,readop,readva,readmp
      integer::  onrc , n0 , nglob !olpp , olppr , olsop , olgrme , olva , olmp , opot  , ocor ,
      parameter(n0=10)
      integer ic,is,lmx,k,nl,nsp,ifi,jfi,fopn,lmxx,nspx,nrx,nr,
     .     idmod(n0),isw,bitand,i2,nclasp,iclbsj,icmap,jc,nbas,nclass,
     .     nclspp,iprint,lgunit,mode0,mode00,lrel
c     integer mpipid,procid
      double precision rhrmx,vrmax(2),ves,z,rmxx,ax,qc,dq,vrmxx(2),
     .     sumec,sumtc,sumev,thrpv,ekin,utot,rhoeps,etot,a,rmax,dval,
     .     zbak(2),dgets,pdf(n0,2),qdf(n0,2),pnuloc(100),qnuloc(100)
      real(8):: olpp(100),olppr(3000),olsop(500),olgrme(500),olva(100),olmp(3000),opot(3000),ocor(3000)
      
      nl=sctrl%nl
      nsp=sctrl%nspin
      nclass=sctrl%nclass
      nclasp=sarray%nclasp
      lrel = globalvariables%lrel
      mode0 = mod(mode,10)
      mode00 = mod(mode0,4)
      i2 = ic2
      if (i2 .eq. 0) i2 = nclasp
      do  10  ic = ic1, i2
         is=sarray%iv_a_oics(ic)
         icmap = iclbsj ( is , sarray%iv_a_oics , - nclasp , 1 ) 
         if (icmap .eq. ic .or. mode00 .eq. 0) icmap = 0
         lmx=sspec(is)%lmxa
         i_copy_size=size(sspec(is)%idmod) 
         call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 
         i_copy_size=size(sspec(is)%p) 
         call dcopy(i_copy_size,sspec(is)%p,1,pdf,1) 
         i_copy_size=size(sspec(is)%q) 
         call dcopy(i_copy_size,sspec(is)%q,1,qdf,1) 
         call dpzero(pnuloc,100)
         call dpzero(qnuloc,100)

C     --- Open the atom file ---
         outs1 = ' '
         outs2 = ' '
         outs3 = ' '
         call r8tos8 ( dval ( sarray%rv_a_oclabl , ic ) , clabl ) 
         if (mode00 .ne. 3) then
            ifi = fopn(clabl)
         endif

C     --- Copy what is passed through to holding arrays ---
         k = nl*nsp
         call pvaioc ( sarray%iv_a_oipq , 1 , ic , icmap , havepq , jc ) 
         readpq = .not. havepq .or. mode00 .eq. 2 .and. jc .ne. ic
         readpq = readpq .and. mode00 .lt. 3
         if (havepq) then
            call dpscop ( spot%rv_a_opnu , pnuloc , k , 1 + ( jc - 1 ) * k , 1  , 1d0 ) 
            call dpscop ( spot%rv_a_oqnu , qnuloc , 3 * k , 1 + ( jc - 1 ) * 3 
     .           * k , 1 , 1d0 ) 
            call dpscop ( spot%rv_a_orhrmx , rhrmx , 1 , jc , 1 , 1d0 ) 
            call dpscop ( spot%rv_a_ovrmax , vrmax , 2 , 2 * jc - 1 , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a pq,',outs1,len(outs1),0)
         endif
         k = 6*nl*nsp
         call pvaioc ( sarray%iv_a_oipq , 2 , ic , icmap , havepp , jc ) 
         readpp = .not. havepp .or. mode00 .eq. 2 .and. jc .ne. ic
         readpp = readpp .and. mode00 .lt. 3
         if (havepp) then
            call dpscop ( spot%rv_a_opp , olpp , k , 1 + ( jc - 1 ) * k  , 1 , 1d0 ) 
            if (lrel .eq. 2) then
               k = 4*nl*2*nl*2*2
               call dpscop ( spot%rv_a_opprel , olppr , k , 1 + ( jc - 1 ) 
     .              * k , 1 , 1d0 ) 
            endif
            call dpscop ( spot%rv_a_oves , ves , 1 , jc , 1 , 1d0 ) 

            if (ic.ne.jc) call awrit0('%a pp,',outs1,len(outs1),0)
         endif
         k = (nl*nsp)**2
         call pvaioc ( sarray%iv_a_oipq , 8 , ic , icmap , haveva , jc ) 
         haveva = haveva .and. allocated(spot%rv_a_ovintr)
         readva = .not. haveva .or. mode00 .eq. 2 .and. jc .ne. ic
         readva = readva .and. allocated(spot%rv_a_ovintr) .and. mode00 .lt. 3
         if (haveva) then
            call dpscop ( spot%rv_a_ovintr , olva , k , 1 + ( jc - 1 ) * k , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a va,',outs1,len(outs1),0)
         endif
         k = nl**2*(2*nl-1)*3*nsp
         call pvaioc ( sarray%iv_a_oipq , 16 , ic , icmap , havemp , jc ) 
         havemp = havemp .and. allocated(spot%rv_a_opmpol)
         readmp = .not. havemp .or. mode00 .eq. 2 .and. jc .ne. ic
         readmp = readmp .and. allocated(spot%rv_a_opmpol) .and. mode00 .lt. 3
         if (havemp) then
            call dpscop ( spot%rv_a_opmpol , olmp , k , 1 + ( jc - 1 ) * k 
     .           , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a mp,',outs1,len(outs1),0)
         endif
         k = nl*nsp*nsp*9
         call pvaioc ( sarray%iv_a_oipq , 4 , ic , icmap , haveso , jc ) 
         haveso = haveso .and. allocated(spot%rv_a_osop)
         readso = .not. haveso .or. mode00 .eq. 2 .and. jc .ne. ic
         readso = readso .and. allocated(spot%rv_a_osop) .and. mode00 .lt. 3
         if (haveso) then
            call dpscop ( spot%rv_a_osop , olsop , k , 1 + ( jc - 1 ) * 
     .           k , 1 , 1d0 ) 

            if (ic.ne.jc) call awrit0('%a sop,',outs1,len(outs1),0)
         endif
         k = 16*nl*nsp
         call pvaioc ( sarray%iv_a_oipq , 32 , ic , icmap , haveop , jc ) 
         haveop = haveop .and. allocated(spot%rv_a_ogrrme)
         readop = .not. haveop .or. mode00 .eq. 2 .and. jc .ne. ic
         readop = readop .and. allocated(spot%rv_a_ogrrme) .and. mode00 .lt. 3
         if (haveop) then
            call dpscop ( spot%rv_a_ogrrme , olgrme , k , 1 + ( jc - 1 )* k , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a opp,',outs1,len(outs1),0)
         endif

C     --- File WRITE ---
         if (lio) then
            lgen = .false.
            lpot = .false.
            lcor = .false.
C     ... Pick up GEN and POT, if available, to save again
            if (scat(iabs(ifi),'GEN:',':',.true.)) then
               lgen = aiogen(alabel,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,dq,
     .              vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
            endif
            if (scat(iabs(ifi),'POT:',':',.true.)) then
               read(ifi,102) nr,nsp,a,rmax
 102           format(2i5,2f12.5)
               lpot = aiopot(nr,nsp,a,rmax,-99d0,opot,ifi)
            endif
            lcor = aiocor(nr,nsp,a,rmxx,ocor,sumec,sumtc,ifi)

            rewind ifi
            jfi = -ifi
            if (lgen) sw = aiogen(clabl,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,
     .           dq,vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,jfi)
            if (havepq) sw = aiomom(clabl,pnuloc,qnuloc,idmod,
     .           nl,lmx,nsp,rhrmx,vrmax,jfi)
            if (havepp) sw = aiopar(clabl,lrel,olpp,olppr,ves,nl,
     .           lmx,nsp,jfi)
            if (haveva) sw = aiova(clabl,olva,nl,lmx,nsp,jfi)
            if (havemp) sw = aiomp(clabl,olmp,nl,2*nl-2,nsp,jfi)
            if (haveso) sw = aiosop(clabl,olsop,nl,lmx,nsp,jfi)
            if (haveop) sw = aiorme(clabl,olgrme,nl,nsp,jfi)
            if (lpot)   sw = aiopot(nr,nsp,a,rmax,-99d0,opot,jfi)
            if (lcor) lcor = aiocor(nr,nsp,a,rmxx,ocor,sumec,sumtc,jfi)

C     --- File READ ---
         else

C     ... Copy whatever is available on disk to holding arrays
            if (readpq .or. mode0 .ge. 4) then
               if (readpq) rewind ifi
               if (readpq) readpq = aiomom(clabl,pnuloc,qnuloc,idmod,
     .              nl,lmx,nsp,rhrmx,vrmax,ifi)
C     Couldn't read from atom file ; take default values
               if (readpq) call awrit0('%a pq,',outs2,len(outs2),0)
               if (mode0 .ge. 4 .and. .not. (readpq .or. havepq)) then
C     call dmcpy(pdf,n0,1,pnuloc,nl,1,nl,nsp)
                  call paioc2(nsp,nl,n0,pdf,qdf,pnuloc,qnuloc)
                  call awrit0('%a pq,',outs3,len(outs2),0)
                  call dvset(vrmax,1,2,-.7d0)
                  rhrmx = .1d0
                  readpq = .true.
               endif
            endif
            if (readpp) then
               readpp = aiopar(clabl,lrel,olpp,olppr,ves,nl,lmx,nsp,
     .              ifi)
               if (readpp) call awrit0('%a pp,',outs2,len(outs2),0)
            endif
            if (readso .and. allocated(spot%rv_a_osop) ) then
               readso = aiosop(clabl,olsop,nl,lmx,nsp,ifi)
               if (readso) call awrit0('%a so,',outs2,len(outs2),0)
            endif
            if (readop .and. allocated(spot%rv_a_ogrrme)) then
               readop = aiorme(clabl,olgrme,nl,nsp,ifi)
               if (readop) call awrit0('%a op,',outs2,len(outs2),0)
            endif
            if (readva .and. allocated(spot%rv_a_ovintr) ) then
               readva = aiova(clabl,olva,nl,lmx,nsp,ifi)
               if (readva) call awrit0('%a va,',outs2,len(outs2),0)
            endif
            if (readmp .and. allocated(spot%rv_a_opmpol) ) then
               readmp = aiomp(clabl,olmp,nl,2*nl-2,nsp,ifi)
               if (readmp) call awrit0('%a mp,',outs2,len(outs2),0)
            endif

C     ... Update what parameters are available
            sarray%iv_a_oipq(ic) = isw(haveop.or.readop)*32+
     .           isw(havemp.or.readmp)*16+
     .           isw(haveva.or.readva)*8 +
     .           isw(haveso.or.readso)*4 +
     .           isw(havepp.or.readpp)*2 +
     .           isw(havepq.or.readpq)*1 +
     .           sarray%iv_a_oipq(ic) - bitand(sarray%iv_a_oipq(ic),63)
            k = nl*nsp
            if (havepq .or. readpq) then
               call dpscop ( pnuloc , spot%rv_a_opnu , k , 1 , 1 + ( ic - 1 ) * k, 1d0 ) 
               call dpscop ( qnuloc , spot%rv_a_oqnu , 3 * k , 1 , 1 + ( ic - 1 ) 
     .              * 3 * k , 1d0 ) 
               call dpscop ( rhrmx , spot%rv_a_orhrmx , 1 , 1 , ic , 1d0 ) 
               call dpscop ( vrmax , spot%rv_a_ovrmax , 2 , 1 , 2 * ic - 1 , 1d0 ) 
            endif
            k = 6*nl*nsp
            if (havepp .or. readpp) then
               call dpscop ( olpp , spot%rv_a_opp , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 
               if (lrel .eq. 2) then
                  k = 4*nl*2*nl*2*2
                  call dpscop ( olppr , spot%rv_a_opprel , k , 1 , 1 + ( ic - 
     .                 1 ) * k , 1d0 ) 
               endif
               call dpscop ( ves , spot%rv_a_oves , 1 , 1 , ic , 1d0 ) 
            endif
            k = (nl*nsp)**2
            if ( ( haveva .or. readva ) .and. allocated(spot%rv_a_ovintr ) ) call 
     .           dpscop ( olva  , spot%rv_a_ovintr , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 

            k = nl**2*(2*nl-1)*3*nsp
            if ( ( havemp .or. readmp ) .and. allocated(spot%rv_a_opmpol ) ) call 
     .           dpscop ( olmp  , spot%rv_a_opmpol , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 
            k = nl*nsp*nsp*9
            if ( ( haveso .or. readso ) .and. allocated(spot%rv_a_osop ) ) 
     .           call dpscop ( olsop  , spot%rv_a_osop , k , 1 , 1 + ( ic - 1) * k , 1d0 ) 

            k = 16*nl*nsp
            if ( ( haveop .or. readop ) .and. allocated(spot%rv_a_ogrrme ) ) 
     .           call dpscop ( olgrme  , spot%rv_a_ogrrme , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 
         endif

         if (mode00 .lt. 3) call fclr(clabl,ifi)
         if (iprint() .gt. 40) then
            outs = ' '
            if (outs1 .ne. ' ') then
               call awrit1('%x '//clabl//'%a: copied '//outs1//'%a%b from '//
     .              'class %i',outs,len(outs),0,jc)
            endif
            if (outs2 .ne. ' ') then
               call awrit0('%x '//clabl//'%a: read '//outs2//'%a%b from '//
     .              'disk',outs,len(outs),0)
            endif
            if (outs3 .ne. ' ') then
               if (outs2 .eq.  ' ') then
                  call awrit0('%x '//clabl//'%a: use defaults for: '//outs3//
     .                 '%a%b',outs,len(outs),0)
               else
                  call awrit0('%a; use defaults for: '//outs3//'%a%b',
     .                 outs,len(outs),0)
               endif
            endif
            if (outs .eq. ' ')
     .           call awrit0(' '//clabl//'%a: nothing read',outs,len(outs),0)
            call awrit0(' aiocls class'//outs,' ',-len(outs),lgunit(1))
         endif
 10   continue
      end subroutine aiocls 

      subroutine pvaioc(initc,mask,ic0,icmap,lhave,ic)
C     - Find whether data avaiable either in class or mapped class
C     implicit none
      logical lhave
      integer initc(1),mask,ic0,icmap
      integer ic
      ic = ic0
      lhave = mod(initc(ic)/mask,2) .eq. 1
      if (.not. lhave .and. icmap .ne. 0) then
         lhave = mod(initc(icmap)/mask,2) .eq. 1
         ic = icmap
      endif
      end subroutine pvaioc 

      subroutine paioc2(nsp,nl,n0,pat,qat,pnu,qnu)
C     - Widget to copy pat,qat to pnu,qnu
C     implicit none
      integer n0,nl,nsp
      double precision pat(n0,2),qat(n0,2),pnu(nl,nsp),qnu(3,nl,nsp)
      integer i,il
      do  i = 1, nsp
         do  il = 1, nl
C     pnu(il,i) = int(pat(il,i)) + .5d0
            pnu(il,i) = pat(il,i)
            qnu(1,il,i) = qat(il,1)/nsp
            if (nsp .eq. 2) then
               if (pat(il,i) .eq. 0) pnu(il,i) = pnu(il,1)
               qnu(1,il,i) = qat(il,1)/nsp + qat(il,2)/2*dble(3-2*i)
            endif
            qnu(2,il,i) = 0d0
            qnu(3,il,i) = 0d0
         enddo
      enddo
      end subroutine paioc2 


