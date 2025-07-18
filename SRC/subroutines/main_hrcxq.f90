!> Calculate Im(chi0) and do Hilbert transformation.
module m_hrcxq
  contains
subroutine hrcxq() bind(C)
  !  Output: rxcq.iq files.
  !   After set up a kind of enviromental variables, by calling module functions,
  !   we read tetrahedron weight via 'call X0kf_v4hz_init_read(iq,is)'.
  !   Then we calculate Im(chi0) by x0kf_v4hz.
  !  Module coding rule:
  !    (1) Read files and readin data are stored in modules. All data in modules are protected.
  !    (2) To set Enviromental variables before main loop (for exaple, do 1001 in this code),
  !        we call module funcitons one by one. It is a bootstrap sequence of calling modules.
  !    (3) During the main loop, a few of module variables are rewritten by module functions
  !        (tetrahedron weight, matrix elements ...). Be careful, and clarify it.
  !    (4) Do now write long fortran program. One MPI loop and one OpenMP loop.
  use m_ReadEfermi,only: Readefermi,ef
  use m_readqg,only: Readngmx2,ngpmx,ngcmx
  use m_hamindex,only: Readhamindex
  use m_readeigen,only: Init_readeigen,Init_readeigen2,Readeval
  use m_read_bzdata,only:Read_bzdata,nqbz,nqibz,n1,n2,n3,ginv,dq_,qbz,wbz,qibz,wibz, ntetf,idtetf,ib1bz, qbzw,nqbzw,q0i,nq0i,nq0iadd !for tetrahedron   !     &     idteti, nstar,irk,nstbz
  use m_genallcf_v3,only: Genallcf_v3,natom,nspin,nl,nn,nlnmx,nctot,alat,il,in,im,nlnm,plat,pos,ecore,&
       nband
  use m_rdpp,only: Rdpp,nxx,lx,nx,mdimx,nbloch,cgr,ppbrd ,nblochpmx,mrecl,nprecx ! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
  use m_zmel,only: Mptauof_zmel !Set data for "call get_zmelt" zmelt= matrix element <phi |phi MPB>.
  use m_itq,only: Setitq !set itq,ntq,nband,ngcmx,ngpmx to m_itq
  use m_freq,only: Getfreq2,frhis,freq_r,freq_i,nwhis,nw_i,nw,npm,niw ! Frequency !output of getfreq
  use m_tetwt,only: Tetdeallocate,Gettetwt,whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  use m_w0w0i,only:       W0w0i,    w0,w0i,llmat   !      use m_readq0p,only:     Readq0p,  wqt,q0i,nq0i ,nq0iadd,ixyz
  use m_readgwinput,only: ReadGwinputKeys,egauss,ecut,ecuts,mtet,ebmx,nbmx,imbas
  use m_qbze,only:    Setqbze,nqbze,nqibze,qbze,qibze
!  use m_readhbe,only: Readhbe,nband !,nprecb,mrecb,mrece,nqbzt,nband,mrecg !  use m_eibz,only:    Seteibz,nwgt,neibz,igx,igxt,eibzsym
  use m_x0kf,only: x0kf_zxq,deallocatezxq,deallocatezxqi
  use m_llw,only: WVRllwR,WVIllwI,w4pmode,MPI__sendllw
  use m_mpi,only: MPI__Initialize,MPI__root,MPI__rank,MPI__size,MPI__consoleout,comm, &
                & MPI__SplitXq, MPI__Setnpr_col, comm_b, comm_k, mpi__root_k, mpi__root_q,ipr
  use m_lgunit,only: m_lgunit_init,stdo
  use m_ftox
  use m_readVcoud,only: Readvcoud,ngb
  use m_gpu,only: gpu_init
!  use m_dpsion,only: dpsion5
  implicit none
  real(8),parameter:: pi = 4d0*datan(1d0),fourpi = 4d0*pi,sqfourpi= sqrt(fourpi)
  integer:: iq,kx,ixc,iqxini,iqxend,is,iw,ifwd,ngrpx,verbose,npr,nmbas,ifif
  integer:: i_red_npm,i_red_nwhis,ierr,ircxq,npmx
  integer:: ipart,iwhis,igb1,imb,igb2,isf
  real(8):: ua=1d0,qp(3),quu(3),hartree,rydberg,schi=-9999,q00(3)
  logical :: debug=.false. ,realomega,imagomega !,nolfco=.false.,crpa=.false.
  logical :: hx0,iprintx=.false.,chipm=.false.,localfieldcorrectionllw   !eibzmode,eibz4x0,
  logical:: cmdopt2,emptyrun,cmdopt0
  character(10) :: i2char
  character(20):: outs=''
  real(8),allocatable :: symope(:,:),ekxx1(:,:),ekxx2(:,:)
  complex(8),allocatable:: zxq(:,:,:),zxqi(:,:,:),zzr(:,:)!,rcxq(:,:,:,:)
  logical,allocatable::   mpi__Qtask(:)
  integer,allocatable::   mpi__Qrank(:)
  integer :: n_kpara = 1, n_bpara = 1, npr_col, worker_inQtask
  call MPI__Initialize()
  call gpu_init(comm) 
  call M_lgunit_init()
  emptyrun=cmdopt0('--emptyrun')
  call MPI__consoleout('hrcxq')
  call cputid (0)
  if(verbose()>=100) debug= .TRUE. 
  hartree= 2d0*rydberg()
  call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
  call Read_BZDATA(hx0)      !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
  call Readefermi() !Readin EFERMI
!  call Readhbe()    !Read dimensions
  call ReadGWinputKeys() !Readin dataset in GWinput   !      call Readq0p()    !Readin Offset Gamma
  call Readngmx2()  !Get ngpmx and ngcmx in m_readqg
  call Setqbze()    ! extented BZ points list
  !  write(stdo,*)' ngcmx ngpmx=',ngcmx,ngpmx !ngcmx: max of PWs for W,ngpmx: max of PWs for phi
  !! Get space-group transformation information. See header of mptaouof.
  !! But we only use symops=E in hx0fp0 mode. c.f. hsfp0.sc
  call Mptauof_zmel(symops=reshape([1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0],[3,3]),ng=1)
  !! Rdpp gives ppbrd: radial integrals and cgr = rotated cg coeffecients. --> call Rdpp(ngrpx,symope) is moved to Mptauof_zmel \in m_zmel
  call Setitq()         ! Set itq in m_zmel
  call Readhamindex()
  call Init_readeigen() ! Initialization of readEigen !readin m_hamindex
  call Init_readeigen2()
  realomega = .true.
  imagomega = .true.
  call Getfreq2(.false.,realomega,imagomega,ua,iprintx) ! Getfreq gives frhis,freq_r,freq_i, nwhis,nw,npm
  if(MPI__root) call writewvfreq() 
  ! nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB
  iqxini = 1
  iqxend = nqibz + nq0i + nq0iadd ! [iqxini:iqxend] range of q points.

  if(cmdopt2('--nk=', outs)) read(outs,*) n_kpara
  n_bpara = max(mpi__size/(n_kpara*(iqxend - iqxini + 1)), 1) !Default setting of parallelization. k-parallel is 1.
  if(cmdopt2('--nb=', outs)) read(outs,*) n_bpara
  worker_inQtask = n_bpara * n_kpara
  if(ipr) write(stdo,'(1X,A,3I5)') 'MPI: worker_inQtask:(n_bpara,n_kpara)', worker_inQtask, n_bpara, n_kpara
  call MPI__SplitXq(n_bpara, n_kpara)
  ! allocate( mpi__Qrank(iqxini:iqxend), source=[(mod(iq-1,mpi__size)           ,iq=iqxini,iqxend)])
  ! allocate( mpi__Qtask(iqxini:iqxend), source=[(mod(iq-1,mpi__size)==mpi__rank,iq=iqxini,iqxend)])
  allocate( mpi__Qrank(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)*worker_inQtask           ,iq=iqxini,iqxend)])
  allocate( mpi__Qtask(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)==mpi__rank/worker_inQtask,iq=iqxini,iqxend)])
  if(ipr) write(stdo,ftox)'mpi_rank',mpi__rank,'mpi__Qtask=',mpi__Qtask
  if(ipr) write(stdo,ftox) 'mpi_qrank', mpi__qrank
  call flush(stdo)
  if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
  MainLoopToObtainZxq: do 1001 iq = iqxini,iqxend
    if( .NOT. MPI__Qtask(iq) ) cycle
    if(ipr) write(stdo,*)'mpi_rank in IQ loop:', iq, mpi__rank
    call cputid (0)
    qp = qibze(:,iq)
    if(ipr) write(stdo,ftox)'do 1001: iq q=',iq,ftof(qp,4),' of nq=',iqxend !4 means four digits below decimal point (optional).
    call Readvcoud(qp, iq,NoVcou=chipm) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix
    npr=ngb
    call MPI__Setnpr_col(npr, npr_col) ! set the npr_col : split of npr(column) for MPI color_b
    call x0kf_zxq(realomega,imagomega,qp,iq,npr,schi, crpa=.false.,chipm=.false.,nolfco=.false.) !Get zxq,zxqi in m_x0kf
    if(debug) print *,'sumchk zxq=',sum(zxq),sum(abs(zxq)),' zxqi=',sum(zxqi),sum(abs(zxqi))
    if(mpi__root_k) then
      call WVRllwR(qp,iq,npr,npr_col) !WV=W-v in RandomPhaseApproximation along realaxis. Write big files WVR.  --emptyrun in it
      call deallocatezxq()
      call WVIllwI(qp,iq,npr,npr_col) !WV=W-v along imagaxis Write WVI --emptyrun in it 
      call deallocatezxqi()
    endif
    call mpi_barrier(comm_k, ierr)
    call mpi_barrier(comm_b, ierr)
1001 enddo MainLoopToObtainZxq
   GetEffectiveWVatGammaCell: block !Get W-v(q=0): Divergent part and non-analytic constant part of W(0) calculated from llw
    ! we have wing elemments: llw, llwi LLWR, LLWI
    call MPI_barrier(comm,ierr)
    call MPI__sendllw(iqxend,MPI__Qrank) ! Send all LLW data to mpi_root.
    ! Get effective W0,W0i, and L(omega=0) matrix. Modify WVR WVI with w0 and w0. Files WVI and WVR are modified.
    if(MPI__rank==0) call W0w0i(nw_i,nw,nq0i,niw,q0i) 
  endblock GetEffectiveWVatGammaCell
  if(ipr) write(stdo,ftox) '--- end of hrcxq --- irank=',MPI__rank
  call cputid(0)
  call rx0( ' OK! hrcxq WV generated')
  
  contains
  subroutine writewvfreq() !writeonly
     open(newunit=ifwd, file='__WV.d')
     write(ifwd,"(1x,10i14)") nprecx, mrecl, nblochpmx, nw+1,niw, nqibz + nq0i-1, nw_i
     close(ifwd)
     open(newunit=ifif,file='freq_r') ! Write number of frequency points nwp and frequensies
     write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
     do iw= nw_i,-1
        write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw !negative frequecncies for x0
     enddo
     do iw= 0,nw
        write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw    !positive frequecncies for x0
     enddo
     close(ifif)
   end subroutine writewvfreq
end subroutine hrcxq
end module m_hrcxq
