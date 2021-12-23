      program hrcxq
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
!      
      use m_ReadEfermi,only: Readefermi,ef
      use m_readqg,only: Readngmx2,ngpmx,ngcmx
      use m_readeigen,only: Init_readeigen,Init_readeigen2,Readeval
      use m_read_bzdata,only: Read_bzdata, !<--- 'call read_bzdata' sets up following data.
     &   ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,qbas,ginv,
     &   dq_,qbz,wbz,qibz,wibz,
     &     ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
c     &     idteti, nstar,irk,nstbz
      use m_genallcf_v3,only: Genallcf_v3,
     &     nclass,natom,nspin,nl,nn,ngrp,
     &     nlmto,nlnmx, nctot,!niw, !nw_input=>nw,
     &     alat, symgrp,clabl,iclass, !diw,dw,delta,deltaw,
     &     invg, il, in, im, nlnm, plat, pos, ecore, symgg 
!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: Rdpp,      !"call rdpp" generate following data.
     &     nxx,lx,nx,mdimx,nbloch,cgr,ppbrd ,nblochpmx,mrecl,nprecx
!! Set data for "call get_zmelt" zmelt= matrix element <phi |phi MPB>.
      use m_zmel,only:  !these data set are stored in this module, and used when 
     &      Mptauof_zmel, Setppovlz !Ppbafp_v2_zmel,
c     &        ppovlz,  shtvg, miat,tiat  !ppbir,itq, ntq, nband,ngcmx,ngpmx,
      use m_itq,only: Setitq !set itq,ntq,nband,ngcmx,ngpmx to m_itq
!! Frequency
      use m_freq,only: Getfreq2,
     &     frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,niw !output of getfreq
!! Antiferro
c     use m_anf,only: anfcond,
c     & laf,ibasf !,ldima,pos,natom
!! Tetwt
      use m_tetwt,only: Tetdeallocate, Gettetwt, !followings are output of 'L871:call gettetwt')
     &                                   whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
      use m_w0w0i,only:       W0w0i,     w0,w0i,llmat
      use m_readq0p,only:     Readq0p,   wqt,q0i,nq0i ,nq0iadd,ixyz
      use m_readVcoud,only:   Readvcoud, vcousq,zcousq,ngb,ngc
      use m_readgwinput,only: ReadGwinputKeys,
     &                                   egauss,ecut,ecuts,nbcut,nbcut2,mtet,ebmx,nbmx,imbas
      use m_qbze,only:    Setqbze, nqbze,nqibze,qbze,qibze
      use m_readhbe,only: Readhbe, nband !, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      use m_eibz,only:    Seteibz, nwgt,neibz,igx,igxt,eibzsym
      use m_x0kf,only:    X0kf_v4hz, X0kf_v4hz_symmetrize, X0kf_v4hz_init,x0kf_v4hz_init_write,x0kf_v4hz_init_read,x0kf_zmel
      use m_llw,only:     WVRllwR,WVIllwI,w4pmode,MPI__sendllw
!! MPI
      use m_mpi,only: MPI__hx0fp0_rankdivider2Q, MPI__Qtask, 
     &     MPI__Initialize, MPI__Finalize,MPI__root,
     &     MPI__Broadcast, MPI__rank,MPI__size, MPI__consoleout,MPI__barrier  
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
      integer:: i_red_npm, i_red_nwhis,  i_red_nmbas2,ierr,ircxq,npmx
      character*8 :: charext
      character(20):: outs=''
      integer:: ipart
      logical:: cmdopt2
!-------------------------------------------------------------------------
      call MPI__Initialize()
      call MPI__consoleout('hx0fp0_sc')
      call cputid (0)
      if(verbose()>=100) debug=.true.
      hartree= 2d0*rydberg()
      call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
      call Read_BZDATA(hx0)      !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
      if(nclass /= natom) call Rx( ' hx0fp0_sc: nclass /= natom ') !WE ASSUME iclass(iatom)= iatom
      call Readefermi() !Readin EFERMI
      call Readhbe()    !Read dimensions 
      call ReadGWinputKeys() !Readin dataset in GWinput
      call Readq0p()    !Readin Offset Gamma 
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
      call Setitq()         ! Set itq in m_zmel
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
      
!! Obtain rcxq -------------------           
      if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
      do 1001 iq = iqxini,iqxend
        if( .not. MPI__Qtask(iq) ) cycle
        call cputid (0)
        qp = qibze(:,iq)
        write(6,"('do 1001: iq q=',i5,3f9.4)")iq,qp
        call Readvcoud(qp,iq,NoVcou=.false.) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix
        if(iq > nqibz .and. (.not.localfieldcorrectionllw())  ) then 
           nolfco =.true.
           nmbas_in = 1 
        else          ! We usually use localfieldcorrectionllw()=T
           nolfco = .false.
           nmbas_in = ngb
        endif
        nmbas1 = nmbas_in !We (will) use nmbas1 and nmbas2 for block division of matrices.
        nmbas2 = nmbas_in
!! We set ppovlz for calling get_zmelt (get matrix elements) \in m_zmel \in subroutine x0kf_v4hz
        call Setppovlz(qp,matz=.not.eibz4x0())
        allocate( rcxq(nmbas1,nmbas2,nwhis,npm)) 
        rcxq = 0d0
        do 1003 is = MPI__Ss,MPI__Se !is=1,nspin. rcxq is acuumulated for spins
           write(6,"(' ### ',2i4,' out of nqibz+n0qi+nq0iadd nsp=',2i4,' ### ')")iq,is,nqibz+nq0i+nq0iadd,nspin
           isf = is
           call X0kf_v4hz_init_read(iq,is) !readin icount data (index sets and tetrahedron weight)
           call x0kf_v4hz(qp, is,isf, iq, nmbas_in, eibzmode=eibzmode, nwgt=nwgt(:,iq),rcxq=rcxq)
                !  rcxq is accumulating 
 1003   continue
!! Symmetrize and convert to Enu basis (diagonalized basis for the Coulomb matrix).
!!   That is, we get dconjg(tranpsoce(zcousq))*rcxq*zcousq for eibzmode
        if(eibzmode)  then
           call x0kf_v4hz_symmetrize( qp, iq, 
     i       nolfco, zzr, nmbas_in, chipm, eibzmode=eibzmode, eibzsym=eibzsym(:,:,iq),
     o       rcxq=rcxq)              !  crystal symmetry of rcxq is recovered for EIBZ mode.
        endif
!! only output in this program        
        open(newunit=ircxq,file='rcxq.'//trim(charext(iq)),form='unformatted')
        write(ircxq) nmbas1,nmbas2
        write(ircxq) rcxq
        close(ircxq)
        deallocate(rcxq)
 1001 continue 
      write(6,*) '--- end of hrcxq --- irank=',MPI__rank
      call cputid(0)
      call rx0( ' OK! hrcxq')
      end program hrcxq
