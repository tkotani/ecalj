      program hhilbert
!!  Hilbert transformation Img(chi0) to W-V matrix
!!    dpsion5: calculate real part by the Hilbert transformation from the Im part
!!      LLW: wing part of dielectric function. See Eq.(40) in PRB81 125102
!!    W0W0i: Rewrite W-V at q=0. We have to define effective W(q=0) to avoid divergece.
!!
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
      use m_hamindex,only:   Readhamindex
      use m_readeigen,only: Init_readeigen,Init_readeigen2,Readeval
      use m_read_bzdata,only: Read_bzdata, !<--- 'call read_bzdata' sets up following data.
     &   ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,qbas,ginv,
     &   dq_,qbz,wbz,qibz,wibz,
     &     ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
c     &     idteti, nstar,irk,nstbz
      use m_genallcf_v3,only: Genallcf_v3,
     &     nclass,natom,nspin,nl,nn,
     &     nlmto,nlnmx, nctot, alat, clabl,iclass, 
     &     il, in, im, nlnm, plat, pos, ecore
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
c      use m_eibz,only:    Seteibz, nwgt,neibz,igx,igxt,eibzsym
      use m_x0kf,only:    X0kf_v4hz, X0kf_v4hz_symmetrize, X0kf_v4hz_init,x0kf_v4hz_init_write,x0kf_v4hz_init_read
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
      logical :: hx0, eibzmode, crpa, eibz4x0,iprintx=.false.,chipm=.false., localfieldcorrectionllw
      integer:: i_red_npm, i_red_nwhis,  i_red_nmbas2,ierr,ircxq,npmx
      character*8 :: charext
      character(20):: outs=''
      integer:: ipart
      logical:: cmdopt2
!-------------------------------------------------------------------------
      call MPI__Initialize()
      call MPI__consoleout('hhilbert')
      if(MPI__root) iprintx=.true.
      call cputid (0)
      if(verbose()>=100) debug=.true.
      write(6,*) ' --- hhilbert ----------------'
      hartree= 2d0*rydberg()
      call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
      call Read_BZDATA(hx0)      !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
      if(nclass /= natom) call Rx( ' hhilbert: nclass /= natom ') ! WE ASSUME iclass(iatom)= iatom
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
!!   --> call Rdpp(ngrpx,symope) is moved to Mptauof_zmel \in m_zmel
      call Readhamindex()
      call Init_readeigen() ! Initialization of readEigen !readin m_hamindex
c      call Init_readeigen2() 
!! Getfreq gives frhis,freq_r,freq_i, nwhis,nw,npm 
      realomega = .true.
      imagomega = .true.
      call Getfreq2(.false.,realomega,imagomega,ua,iprintx)!tetra,
!! We first accumulate Imaginary parts. Then do K-K transformation to get real part.
! nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB 
      iqxini = 1 
      iqxend = nqibz + nq0i + nq0iadd ! [iqxini:iqxend] range of q points.
!! Rank divider
      call MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
      if(sum(qibze(:,1)**2)>1d-10) call rx(' hhilbert: sanity check. |q(iqx)| /= 0')
!!-Hilbert transformation -----------
      do 1201 iq = iqxini,iqxend
        if( .not. MPI__Qtask(iq) ) cycle
        qp = qibze(:,iq)
!! Read Coulomb matrix 
        call Readvcoud(qp,iq,NoVcou=.false.) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix 
!! Read rcxq        
        open(newunit=ircxq,file='rcxq.'//trim(charext(iq)),form='unformatted')
        read(ircxq) nmbas1,nmbas2
        if(allocated(rcxq)) deallocate(rcxq)
        allocate( rcxq(nmbas1,nmbas2,nwhis,npm))
        read(ircxq) rcxq
        close(ircxq)
!! Hilbert transform . Genrerate Real part from Imaginary part. ======
        if(realomega) allocate( zxq(nmbas1,nmbas2,nw_i:nw) )
        if(imagomega) allocate( zxqi(nmbas1,nmbas2,niw)    )
        write(6,'("goto dpsion5: nwhis nw_i niw nw_w nmbas1 nmbas2=",6i5)')
     &       nwhis,nw_i,nw,niw,nmbas1,nmbas2
        call dpsion5(           ! frhis,nwhis, freq_r, nw, freq_i,niw,
     i   realomega, imagomega, 
     i   rcxq, nmbas1,nmbas2,   ! rcxq is alterd---used as work npm,nw_i, 
     o   zxq, zxqi,
     i   chipm, schi,is,  ecut,ecuts)
        deallocate(rcxq)
        if(debug) print *,'sumchk zxq=',sum(zxq),sum(zxqi),sum(abs(zxq)),sum(abs(zxqi))
!! ===  RealOmega === W-V: WVR and WVI. Wing elemments: llw, llwi LLWR, LLWI
        if (realomega) then
           call WVRllwR(qp,iq,zxq,nmbas1,nmbas2)
           deallocate(zxq)
        endif 
!! === ImagOmega ===
        if (imagomega) then
           call WVIllwI(qp,iq,zxqi,nmbas1,nmbas2)
           deallocate(zxqi)
        endif 
!! === ImagOmega end ===
 1201 continue
!! == Divergent part and non-analytic constant part of W(0) ==
      call MPI__barrier()
      call MPI__sendllw(iqxend) !!! Send all LLW data to mpi_root.
!! Get effective W0,W0i, and L(omega=0) matrix. Modify WVR WVI
!!  With w0 and w0i, we modify W0W0i. Files WVI and WVR are modified. jun2020
      if(MPI__rank==0) call W0w0i(nw_i,nw,nq0i,niw,q0i)
      write(6,*) '--- end of hhilbert --- irank=',MPI__rank
      call cputid(0)
      call rx0( ' OK! hhilbert')
      end program hhilbert
