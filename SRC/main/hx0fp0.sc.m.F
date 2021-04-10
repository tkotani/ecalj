      program hx0fp0_sc
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
      use m_readqg,only: Readngmx2,Readqg0,ngpmx,ngcmx
      use m_readeigen,only: Init_readeigen,Init_readeigen2,Readeval
      use m_read_bzdata,only: Read_bzdata, !<--- 'call read_bzdata' sets up following data.
     &   ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,qbas,ginv,
     &   dq_,qbz,wbz,qibz,wibz,
     &     ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
c     &     idteti, nstar,irk,nstbz
      use m_genallcf_v3,only: Genallcf_v3,
     &     nclass,natom,nspin,nl,nn,ngrp,
     &     nlmto,nlnmx, nctot,!niw, !nw_input=>nw,
     &     alat, deltaw,symgrp,clabl,iclass, !diw,dw,delta,
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
      use m_x0kf,only:    X0kf_v4hz, X0kf_v4hz_symmetrize, X0kf_v4hz_init
      use m_llw,only:     WVRllwR,WVIllwI,MPI__sendllw
      use m_w0w0i,only:   Finalizew4p,w0w0i
!! MPI
      use m_mpi,only: MPI__hx0fp0_rankdivider2Q,MPI__hx0fp0_rankdivider2S,
     &     MPI__Qtask,MPI__InitializeQSPBM,MPI__Finalize,MPI__root,
     &     MPI__Broadcast,MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__rank,MPI__size,
     &     MPI__Qranktab,MPI__consoleout,MPI__Ss,MPI__Se, MPI__allreducesumS,
     &     MPI__barrier, MPI__rankQ,MPI__rootQ,MPI__rootS
!! ------------------------------------------------------------------------      
      implicit none
      real(8),parameter:: pi = 4d0*datan(1d0),fourpi = 4d0*pi, sqfourpi= sqrt(fourpi)
      integer:: iq,isf,kx,ixc,iqxini,iqxend,is,iw,ifwd,ngrpx,verbose,nmbas1,nmbas2,nmbas_in,ifif
      real(8):: ua=1d0, q(3), quu(3), hartree, rydberg, schi=-9999
      real(8),allocatable :: symope(:,:), ekxx1(:,:),ekxx2(:,:)
      complex(8),allocatable:: zxq(:,:,:),zxqi(:,:,:),zzr(:,:), rcxq(:,:,:,:)
      logical :: debug=.false. , realomega, imagomega, nolfco=.false.
      logical :: hx0, eibzmode, crpa, eibz4x0,iprintx=.false.,chipm=.false., localfieldcorrectionllw
      integer:: i_reduction_npm, i_reduction_nwhis,  i_reduction_nmbas2
!-------------------------------------------------------------------------
      call MPI__InitializeQSPBM()
      call MPI__consoleout('hx0fp0_sc')
      call cputid (0)
      if(verbose()>=100) debug=.true.
      write(6,*) ' --- hx0fp0_sc Choose modes below ----------------'
      write(6,*) '  ixc= 11,10011' !,or 1011 '
      if( MPI__root ) then
        read(5,*) ixc 
      endif
      call MPI__Broadcast(ixc)
      if(MPI__root) iprintx=.true.
      crpa = .false.
      if(ixc==11) then
         write(6,*) " OK ixc=11 normal mode "
      elseif(ixc==10011) then
         write(6,*) " OK ixc=10011 crpa mode "
         crpa=.true.
c$$$      elseif(ixc==1011) then
c$$$         write(6,*) 'OK ixc=1011 Add W0W0I part at q=0'
      else
         write(6,*)'we only allow ixc==11 or 10011. Given ixc=',ixc
         call Rx( 'error: give allowed arg for hx0fp0_sc.')
      endif
      allocate(zzr(1,1))       !dummy !zzr is required for chi^+- mode for hx0fp0
      hartree= 2d0*rydberg()
      call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
      call Read_BZDATA(hx0)      !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
!!!!  CAUTION: WE ASSUME iclass(iatom)= iatom,  nclass = natom.  !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) call Rx( ' hx0fp0_sc: nclass /= natom ')
      call Readefermi() !Readin EFERMI
      call Readhbe()    !Read dimensions 
      call ReadGWinputKeys() !Readin dataset in GWinput
c      write(6,"(' nbcut nbcut2=',2i5)") nbcut,nbcut2
c      if(ngrp/= ngrp2)  call rx( 'ngrp inconsistent: BZDATA and GWIN_V2')
c      if(nlmto/=nlmtot) call rx( ' hx0fp0_sc: nlmto/=nlmtot in hbe.d')
c      if(nqbz /=nqbzt ) call rx( ' hx0fp0_sc: nqbz /=nqbzt  in hbe.d')
!! --- Readin Offset Gamma --------
      call Readq0p()
c      write(6,"(' ### nqibz nq0i nq0iadd=', 3i5)")nqibz,nq0i,nq0iadd
!! --- Read Qpoints
c      call Readngmx('QGpsi',ngpmx), call Readngmx('QGcou',ngcmx)
      call Readngmx2()  !Get ngpmx and ngcmx in m_readqg
      call Setqbze()    ! extented BZ points list
c      write(6,*)' ngcmx ngpmx=',ngcmx,ngpmx !ngcmx: max of PWs for W, ngpmx: max of PWs for phi
!! Get space-group transformation information. See header of mptaouof.
!! But we only use symops=E for hx0fp0 mode. c.f hsfp0.sc 
      ngrpx = 1 !no space-group symmetry operation in hx0fp0. ng=1 means E only.
      allocate(symope(3,3))
      symope(1:3,1) = [1d0,0d0,0d0]
      symope(1:3,2) = [0d0,1d0,0d0]
      symope(1:3,3) = [0d0,0d0,1d0]
      call Mptauof_zmel(symope,ngrpx) !Since we only use E, 
c      if(verbose()>=40) write (*,*)' hsfp0.sc.m.F: end of mptauof_x'
!! Rdpp gives ppbrd: radial integrals and cgr = rotated cg coeffecients.
!!       --> call Rdpp(ngrpx,symope) is moved to Mptauof_zmel \in m_zmel
!! Set itq in m_zmel
      call Setitq()
!! ... initialization of readEigen !readin m_hamindex
      call Init_readeigen() !nband,mrece)!EVU EVD are stored in m_readeigen
      call Init_readeigen2() !mrecb,nlmto,mrecg)
!! Getfreq gives frhis,freq_r,freq_i, nwhis,nw,npm 
      realomega = .true.
      imagomega = .true.
      call Getfreq2(.false.,realomega,imagomega,ua,iprintx)!tetra,
!! Write freq_r (real omega). Read from hsfp0
      if(realomega .and. mpi__root) then
         open(newunit=ifif,file='freq_r') !write number of frequency points nwp and frequensies in 'freq_r' file
         write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
         do iw= nw_i,-1
            write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw 
         enddo
         do iw= 0,nw
            write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw 
         enddo
         close(ifif)
      endif
!! We first accumulate Imaginary parts. Then do K-K transformation to get real part.
c      noccxv = maxocc2 (nspin,ef, nband, qbze,nqbze) 
c        !max no. of occupied valence states
c      if(noccxv>nband) call Rx( 'hx0fp0_sc: all the bands filled! too large Ef')
c     noccx  = noccxv + nctot ! We usually assume nctot=0
c     nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB 
      iqxini = 1 
      iqxend = nqibz + nq0i + nq0iadd ! [iqxini:iqxend] range of q points.
!!      
      if(MPI__root) then  ! I think it is not so meaningful to give a subroutine for writing WV.d.
        open(newunit=ifwd, file='WV.d', action='write')
        write (ifwd,"(1x,10i14)") nprecx, mrecl, nblochpmx, nw+1,niw, nqibz + nq0i-1, nw_i
        close(ifwd)
      endif
!! Array for eigenvalues.
      allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))
!! EIBZ mode
      eibzmode = eibz4x0()
      call Seteibz(iqxini,iqxend,iprintx)
c$$$!! jump for w0w0i-only mode (test mode)
c$$$      if(ixc==1011)then !ixc==11 is a debug mode to test contrib. at \Gamma point.
c$$$         goto 1191
c$$$      endif
!! W4phonon. !still developing...
!!      call Setw4pmode()
!! Rank divider
      call MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
      call MPI__hx0fp0_rankdivider2S(nspin)
      write(6,'("irank=",i5," allocated(MPI__qtask)=",L5)')MPI__rank,allocated(MPI__qtask)
      do iq = iqxini,iqxend
        if(MPI__qtask(iq)) write(6,'("irank iq=",i5,i5)') MPI__rank,iq
      enddo
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
      do 1001 iq = iqxini,iqxend
        if( .not. MPI__Qtask(iq) ) cycle
        call cputid (0)
        q = qibze(:,iq)
        write(6,"('do 1001: iq q=',i5,3f9.4)")iq,q
!! Read Coulomb matrix 
        call Readvcoud(q,iq,NoVcou=.false.) ! Readin vcousq,zcousq ngb ngc for the Coulomb matrix for given q
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
        call Setppovlz(q,matz=.not.eibz4x0())
        allocate( rcxq(nmbas1,nmbas2,nwhis,npm)) 
        rcxq = 0d0
        do 1003 is = MPI__Ss,MPI__Se !is=1,nspin
           write(6,"(' ### ',2i4,' out of nqibz+n0qi+nq0iadd nsp=',2i4,' ### ')") 
     &          iq, is, nqibz + nq0i+nq0iadd, nspin
           if(debug) write(6,*)' niw nw=',niw,nw
           isf = is
           do kx = 1, nqbz
              call readeval(qbz(:,kx),   is,  ekxx1(1:nband, kx) ) !read eigenvalue
              call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
           enddo
           call gettetwt(q,iq,is,isf,nwgt(:,iq),ekxx1,ekxx2,nband,eibzmode) ! Tetrahedron weight.
           call x0kf_v4hz_init(0, q, is, isf, iq,  nmbas_in, eibzmode, nwgt(:,iq),crpa)
           call x0kf_v4hz_init(1, q, is, isf, iq,  nmbas_in, eibzmode, nwgt(:,iq),crpa)
           call x0kf_v4hz(        q, is, isf, iq,  nmbas_in, eibzmode, nwgt(:,iq),crpa, 
     o          rcxq)           !  rcxq is the accumulating variable for spins 
           call tetdeallocate()
 1003   continue
!! Symmetrize and convert to Enu basis (diagonalized basis for the Coulomb matrix).
!!   That is, we get dconjg(tranpsoce(zcousq))*rcxq*zcousq for eibzmode
        if(eibzmode)  then
           call x0kf_v4hz_symmetrize( q, iq, 
     i       nolfco, zzr, nmbas_in, chipm, eibzmode, eibzsym(:,:,iq),
     o       rcxq)              !  crystal symmetry of rcxq is recovered for EIBZ mode.
        endif
!! Reduction rcxq in the Spin-axis
        write(6,*) 'MPI__AllreduceSumS start'
        do i_reduction_npm=1,npm
        do i_reduction_nwhis=1,nwhis
        do i_reduction_nmbas2=1,nmbas2
          call MPI__AllreduceSumS( rcxq(1,i_reduction_nmbas2,i_reduction_nwhis,i_reduction_npm), nmbas1)
        enddo
        enddo
        enddo
        write(6,*) 'MPI__AllreduceSumS end'
!! We only need rcxq hereafter. Correct? 
        if(.not.MPI__rootS) exit
!! Hilbert transform . Genrerate Real part from Imaginary part. ======
        if(realomega) allocate( zxq(nmbas1,nmbas2,nw_i:nw) )
        if(imagomega) allocate( zxqi(nmbas1,nmbas2,niw)    )
        write(6,'("goto dpsion5: nwhis nw_i niw nw_w nmbas1 nmbas2=",6i5)') nwhis,nw_i,nw,niw,nmbas1,nmbas2
        call dpsion5(           !frhis,nwhis, freq_r, nw, freq_i,niw,
     i   realomega, imagomega, 
     i   rcxq, nmbas1,nmbas2, ! rcxq is alterd---used as work npm,nw_i, 
     o   zxq, zxqi,
     i   chipm, schi,is,  ecut,ecuts)
        if(allocated(rcxq) ) deallocate(rcxq)
!! ===  RealOmega === W-V: WVR and WVI. Wing elemments: llw, llwi LLWR, LLWI
        if(debug) print *,'sumchk zxq=',sum(zxq),sum(zxqi),sum(abs(zxq)),sum(abs(zxqi))
        if (realomega) then
           call WVRllwR(q,iq,zxq,nmbas1,nmbas2)
           deallocate(zxq)
        endif 
!! === ImagOmega ===
        if (imagomega) then
           call WVIllwI(q,iq,zxqi,nmbas1,nmbas2)
           deallocate(zxqi)
        endif 
!! === ImagOmega end ===
 1001 continue ! =end of loop over q point =================================

      
!! == W(0) divergent part and W(0) non-analytic constant part.==. Hereafter can be in another fortran program.
      call MPI__barrier()
      call MPI__sendllw(iqxend) !!! Send all LLW data to mpi_root.
      if(MPI__rank == 0 ) then  ! MIZUHO-IR only on root node
c$$$!! ix=1011 is a special mode to overwrite llw and llwI for test purpose
c$$$        if(ixc==1011) then  
c$$$          open(newunit=ifw0w0i,file='W0W0I',form='unformatted')
c$$$          read(ifw0w0i) nw_ixx,nwxx,niw_,nq0ix
c$$$          write(6,*)'w0w0i: n=',nw_ixx,nwxx,niw_,nq0ix
c$$$          if(nq0i/=nq0ix)  call rx('nq0i/=nq0ix')
c$$$          if(nw_i/=nw_ixx) call rx(nw_i/=nw_ixx)
c$$$          allocate( llw(nw_i:nw,1:nq0i), llwI(1:niw_,1:nq0i))
c$$$          if(nw/=nwxx) call rx(nw/=nwxx)
c$$$          read(ifw0w0i) llw(nw_i:nw,1:nq0i)
c$$$          read(ifw0w0i) llwI(1:niw_,1:nq0i)
c$$$          close(ifw0w0i)
c$$$        endif  
!! A file W0W0I is generated by call w0w0i, but usually unused anywhere.
!! get w0 and w0i (diagonal element at Gamma point)
!! This return w0 and w0i. (llw and llwi are input)
!! Outputs w0,w0i,llmat. See use m_w0w0i at the begining of this routine.
!!    With w0 and w0i, wow0i modify WVI and WVR. jun2020
        call W0w0i(nw_i,nw,nq0i,niw,q0i)  !Get effective W0,W0i, and L(omega=0) matrix. Overwrite WVR WVI
        call FinalizeW4p()  
      endif
      write(6,*) '--- end of hx0fp0_sc --- irank=',MPI__rank
      call cputid(0)
c      call flush(6)
c      call MPI__Finalize
      if(ixc==11     ) call rx0( ' OK! hx0fp0_sc ixc=11 Sergey F. mode')
      if(ixc==10011  ) call rx0( ' OK! hx0fp0_sc ixc=11 Sergey F. mode')
c$$$      if(ixc==1011) call rx0( ' OK! hx0fp0_sc ixc=1011 W0W0Ionly mode')
      end program hx0fp0_sc 


