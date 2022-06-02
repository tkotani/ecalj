program hx0zmel
  !  Calculate Im(chi0), writint to rcxq.* files.
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
  use m_hamindex,only:   Readhamindex
  use m_readeigen,only: Init_readeigen,Init_readeigen2,Readeval
  use m_read_bzdata,only: Read_bzdata, nqbz,nqibz,n1,n2,n3,ginv, &
       dq_,qbz,wbz,qibz,wibz, &
       ntetf,idtetf,ib1bz, qbzw,nqbzw,nq0i ,nq0iadd !for tetrahedron
  use m_genallcf_v3,only: Genallcf_v3, nclass,natom,nspin,nl,nn, &
       nlmto,nlnmx, nctot,  alat, clabl,iclass, il, in, im, nlnm, plat, pos, ecore
  !! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
  use m_rdpp,only: Rdpp, nxx,lx,nx,mdimx,nbloch,cgr,ppbrd ,nblochpmx,mrecl,nprecx
  !! Set data for "call get_zmelt" zmelt= matrix element <phi |phi MPB>.
  use m_zmel,only:  Mptauof_zmel, Setppovlz 
  use m_itq,only: Setitq !set itq,ntq,nband,ngcmx,ngpmx to m_itq
  !! Frequency
  use m_freq,only: Getfreq2, frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,niw !output of getfreq
  !! Antiferro
  !     use m_anf,only: anfcond,
  !     & laf,ibasf !,ldima,pos,natom
  !! Tetwt
  use m_tetwt,only: Tetdeallocate, Gettetwt, whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  use m_w0w0i,only:       W0w0i,     w0,w0i,llmat
  !      use m_readq0p,only:     Readq0p,   wqt,q0i,nq0i ,nq0iadd,ixyz
  use m_readVcoud,only:   Readvcoud, vcousq,zcousq,ngb,ngc
  use m_readgwinput,only: ReadGwinputKeys, &
       egauss,ecut,ecuts,nbcut,nbcut2,mtet,ebmx,nbmx,imbas
  use m_qbze,only:    Setqbze, nqbze,nqibze,qbze,qibze
  use m_readhbe,only: Readhbe, nband !, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_eibz,only:    Seteibz, nwgt,neibz,igx,igxt,eibzsym
  use m_x0kf,only:    X0kf_v4hz, X0kf_v4hz_symmetrize, X0kf_v4hz_init, &
       X0kf_v4hz_init_write,X0kf_v4hz_init_read,X0kf_zmel, kc,ncount
  use m_llw,only:     WVRllwR,WVIllwI,w4pmode,MPI__sendllw
  !! MPI
  use m_mpi,only: MPI__hx0fp0_rankdivider2Q, MPI__Qtask, &
       MPI__Initialize, MPI__Finalize,MPI__root, &
       MPI__Broadcast, MPI__rank,MPI__size, MPI__consoleout,MPI__barrier
  use m_lgunit,only: m_lgunit_init
  !! ------------------------------------------------------------------------
  implicit none
  integer:: MPI__Ss,MPI__Se
  real(8),parameter:: pi = 4d0*datan(1d0),fourpi = 4d0*pi, sqfourpi= sqrt(fourpi)
  integer:: iq,isf,kx,ixc,iqxini,iqxend,is,iw,ifwd,ngrpx,verbose,nmbas1,nmbas2,nmbas_in,ifif
  real(8):: ua=1d0, qp(3), quu(3), hartree, rydberg, schi=-9999
  real(8),allocatable :: symope(:,:), ekxx1(:,:),ekxx2(:,:)
  complex(8),allocatable:: zxq(:,:,:),zxqi(:,:,:),zzr(:,:), rcxq(:,:,:,:)
  logical :: debug=.false. , realomega, imagomega, nolfco=.false.
  logical :: hx0, eibzmode, eibz4x0,iprintx=.false.,chipm=.false., localfieldcorrectionllw
  integer:: i_red_npm, i_red_nwhis,  i_red_nmbas2,ierr,ircxq,npmx,kold,isold
  character(8) :: charext
  character(20):: outs=''
  integer:: ipart,icount,k
  logical:: cmdopt2
  !-------------------------------------------------------------------------
  call MPI__Initialize()
  call M_lgunit_init()
  call MPI__consoleout('hx0zmel')
  call cputid (0)
  if(verbose()>=100) debug= .TRUE. 
  hartree= 2d0*rydberg()
  call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
  call Read_BZDATA(hx0)      !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
  if(nclass /= natom) call Rx( ' hx0fp0_sc: nclass /= natom ') !WE ASSUME iclass(iatom)= iatom
  call Readefermi() !Readin EFERMI
  call Readhbe()    !Read dimensions
  call ReadGWinputKeys() !Readin dataset in GWinput
  !      call Readq0p()    !Readin Offset Gamma
  call Readngmx2()  !Get ngpmx and ngcmx in m_readqg
  call Setqbze()    ! extented BZ points list
  !  write(6,*)' ngcmx ngpmx=',ngcmx,ngpmx !ngcmx: max of PWs for W, ngpmx: max of PWs for phi
  !! Get space-group transformation information. See header of mptaouof.
  !! But we only use symops=E in hx0fp0 mode. c.f. hsfp0.sc
  ngrpx = 1 !no space-group symmetry operation in hx0fp0. ng=1 means E only.
  allocate(symope(3,3))
  symope(1:3,1) = [1d0,0d0,0d0]
  symope(1:3,2) = [0d0,1d0,0d0]
  symope(1:3,3) = [0d0,0d0,1d0]
  call Mptauof_zmel(symope,ngrpx)
  !! Rdpp gives ppbrd: radial integrals and cgr = rotated cg coeffecients.
  !!       --> call Rdpp(ngrpx,symope) is moved to Mptauof_zmel \in m_zmel
  call Setitq()             ! Set itq in m_zmel
  call Readhamindex()
  call Init_readeigen() ! Initialization of readEigen !readin m_hamindex
  call Init_readeigen2()
  !! Getfreq gives frhis,freq_r,freq_i, nwhis,nw,npm
  realomega = .true.
  imagomega = .true.
  call Getfreq2(.false.,realomega,imagomega,ua,iprintx)!tetra,
  !! We first accumulate Imaginary parts. Then do K-K transformation to get real part.
  ! nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB
  iqxini = 1
  iqxend = nqibz + nq0i + nq0iadd ! [iqxini:iqxend] range of q points.
  eibzmode = eibz4x0()                ! EIBZ mode
  call Seteibz(iqxini,iqxend,iprintx) ! EIBZ mode
  !!    call Setw4pmode() !W4phonon. !still developing...
  !! Rank divider
  call MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
  MPI__Ss = 1
  MPI__Se = nspin
  !! Obtain zmel
  if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
  do 2501 iq = iqxini,iqxend
     if( .NOT. MPI__Qtask(iq) ) cycle
     kold=-999
     isold=-999
     qp = qibze(:,iq)
     write(6,"('do 1501: iq q=',i5,3f9.4)")iq,qp
     call Readvcoud(qp,iq,NoVcou=.false.) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix
     if(iq > nqibz .AND. ( .NOT. localfieldcorrectionllw())  ) then
        nolfco =.true.
        nmbas_in = 1
     else          ! We usually use localfieldcorrectionllw()=T
        nolfco = .false.
        nmbas_in = ngb
     endif
     call Setppovlz(qp,matz=.not.eibz4x0())
     do 2003 is = MPI__Ss,MPI__Se !is=1,nspin. rcxq is acuumulated for spins
        isf = is
        call X0kf_v4hz_init_read(iq,is) !readin icount data (index sets and tetrahedron weight)
        do 2001 icount = 1,ncount
           k = kc(icount)
           if(k/=kold .OR. is/=isold) then
              call x0kf_zmel(qp,iq,k, is,isf)
              kold=k
              isold=is
           endif
2001    enddo
2003 enddo
2501 enddo
  write(6,*) '--- end of hx0zmel --- irank=',MPI__rank
  call cputid(0)
  call rx0( ' OK! hx0zmel')
end program hx0zmel
