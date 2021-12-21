      program hx0init
!!  Calculate W-V for QSGW mode. Cleaned up jun2020. 
!!  We calculate dielectric chi0 by the follwoing three steps.
!!    gettetwt: tetrahedron weights
!!    x0kf_v4hz: Accumlate Im part of the Lindhard function. Im(chi0)
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
      complex(8),allocatable:: zxq(:,:,:),zxqi(:,:,:), rcxq(:,:,:,:)
      logical :: debug=.false. , realomega, imagomega, nolfco=.false.
      logical :: hx0, eibzmode, crpa, eibz4x0,iprintx=.false.,chipm=.false., localfieldcorrectionllw
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
      write(6,*) ' --- hx0fp0_sc Choose modes below ----------------'
      write(6,*) '  ixc= 11,10011' 
      if(cmdopt2('--job=',outs)) then
        read(outs,*) ixc
      else
        if( MPI__root ) then
        read(5,*) ixc 
        endif
        call MPI__Broadcast(ixc)
      endif
      crpa = .false.
      if(ixc==11) then
         write(6,*) " OK ixc=11 normal mode "
      elseif(ixc==10011) then
         write(6,*) " OK ixc=10011 crpa mode "
         crpa=.true.
      else
         write(6,*)'we only allow ixc==11 or 10011. Given ixc=',ixc
         call Rx( 'error: give allowed arg for hx0fp0_sc.')
      endif
      if(MPI__root) iprintx=.true.
      hartree= 2d0*rydberg()
      call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
      call Read_BZDATA(hx0)      !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
!!!!  CAUTION: WE ASSUME iclass(iatom)= iatom,  nclass = natom.  !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) call Rx( ' hx0fp0_sc: nclass /= natom ')
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
      if(realomega .and. mpi__root) then  ! Write freq_r (real omega). Read from hsfp0.
         open(newunit=ifif,file='freq_r') ! Write number of frequency points nwp and frequensies 
         write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
         do iw= nw_i,-1
            write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw !negative frequecncies for x0
         enddo
         do iw= 0,nw
            write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw    !positive frequecncies for x0
         enddo
         close(ifif)
      endif
!! We first accumulate Imaginary parts. Then do K-K transformation to get real part.
! nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB 
      iqxini = 1 
      iqxend = nqibz + nq0i + nq0iadd ! [iqxini:iqxend] range of q points.
      if(MPI__root) then  ! I think it is not so meaningful to give a subroutine for writing WV.d.
        open(newunit=ifwd, file='WV.d', action='write')
        write (ifwd,"(1x,10i14)") nprecx, mrecl, nblochpmx, nw+1,niw, nqibz + nq0i-1, nw_i
        close(ifwd)
      endif
      allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz)) ! For eigenvalues.
      eibzmode = eibz4x0()                ! EIBZ mode
      call Seteibz(iqxini,iqxend,iprintx) ! EIBZ mode
!!    call Setw4pmode() !W4phonon. !still developing...
!! Rank divider
      call MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
      MPI__Ss = 1
      MPI__Se = nspin
      
!! == Calculate x0(q,iw) and W == main loop 1001 for iq. 
!! NOTE:o iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
!! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
!! (or need to modify do 1001 loop).
!! ---------------------------------------------------------------
!! === do 1001 loop over iq ============================================
!     ! ---------------------------------------------------------------
!! We have to do sum for iq,is,k,it,itp,jpm,iw,igb1,igb2      
!!   iq (q vector IBZ), is (spin)> k (k vector BZ)> it,itp (band)
!!      > jpm,iw (omega) >igb1,igb2 (MPB index)
!!   (usually, we only use jpm=1 only--- This meand no negative omega needed).
!! I think, iq,igb1,igb2,(it,itp) are suitable for decomposition (computation, and memory distribution).
!!     !note The pair (it,itp) gives very limited range of allowed iw.
!!
      
      if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
!! do 1101 is for whw and index for x0kf_v4hz
      do 1101 iq = iqxini,iqxend
        if( .not. MPI__Qtask(iq) ) cycle
        call cputid (0)
        qp = qibze(:,iq)
        write(6,"('do 1001: iq q=',i5,3f9.4)")iq,qp
        do 1103 is = MPI__Ss,MPI__Se !is=1,nspin
           write(6,"(' ### ',2i4,' out of nqibz+n0qi+nq0iadd nsp=',2i4,' ### ')")iq,is,nqibz+nq0i+nq0iadd,nspin
           if(debug) write(6,*)' niw nw=',niw,nw
           isf = is
           do kx = 1, nqbz
              ekxx1(1:nband,kx) = readeval(qbz(:,kx),    is ) ! read eigenvalue
              ekxx2(1:nband,kx) = readeval(qp+qbz(:,kx), isf) !
           enddo
           call gettetwt(qp,iq,is,isf,nwgt(:,iq),ekxx1,ekxx2,nband=nband,eibzmode=eibzmode) ! Tetrahedron weight for x0kf_v4hz
           ierr=x0kf_v4hz_init(0,qp,is,isf,iq,nmbas_in, eibzmode=eibzmode, nwgt=nwgt(:,iq),crpa=crpa)
           ierr=x0kf_v4hz_init(1,qp,is,isf,iq,nmbas_in, eibzmode=eibzmode, nwgt=nwgt(:,iq),crpa=crpa)
           call X0kf_v4hz_init_write(iq,is)
           call tetdeallocate()
 1103   continue
 1101 continue
      write(6,*) '--- end of hx0init --- irank=',MPI__rank
      call cputid(0)
      if(ixc==11     ) call rx0( ' OK! hx0init ixc=11 ')
      if(ixc==10011  ) call rx0( ' OK! hx0init ixc=10011')
      end program hx0init
