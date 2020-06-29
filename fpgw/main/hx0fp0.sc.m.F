      program hx0fp0_sc
!!  Calculate W-V for QSGW mode.
!! We calculate chi0 by the follwoing three steps.
!!  tetwt5: tetrahedron weights
!!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
!!  dpsion5: calculate real part by the Hilbert transformation from the Im part
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
     &     invg, il, in, im, nlnm, 
     &     plat, pos, ecore, symgg 
!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: Rdpp,      !"call rdpp" generate following data.
     &     nxx,lx,nx,mdimx,nbloch,cgr,ppbrd 
!! Set data for "call get_zmelt".
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
     &     whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
!! w0 and w0i (head part at Gamma point)
      use m_w0w0i,only: W0w0i,
     &     w0,w0i,llmat
!! MPI
      use m_mpi,only: MPI__hx0fp0_rankdivider2Q,MPI__hx0fp0_rankdivider2S,
     &     MPI__Qtask,MPI__InitializeQSPBM,MPI__Finalize,MPI__root,
     &     MPI__Broadcast,MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__rank,MPI__size,
     &     MPI__Qranktab,MPI__consoleout,MPI__Ss,MPI__Se, MPI__allreducesumS,
     &     MPI__barrier, MPI__rankQ,MPI__rootQ,MPI__rootS
!! q0p
      use m_readq0p,only: Readq0p,
     &     wqt,q0i,nq0i ,nq0iadd,ixyz
!! vcou
      use m_readVcoud,only: Readvcoud,
     &     vcousq,zcousq
!!
      use m_readgwinput,only: ReadGwinputKeys,
     &     egauss,ecut,ecuts,nbcut,nbcut2,mtet,ebmx,nbmx,imbas
      use m_qbze,only: Setqbze,
     &     nqbze,nqibze,qbze,qibze
      use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      use m_eibz,only:Seteibz, nwgt,neibz,ngrpt,igx,igxt,eibzsym
      use m_x0kf,only: X0kf_v4hz, X0kf_v4hz_symmetrize
!! ------------------------------------------------------------------------      
      implicit none
      integer::maxocc2,iclose, ixc,iqxini,iqxend, !iopen,
c     &     ifhbe,  nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &     i,nq0ix,ngrpmx,mxx,ini,ix,ngrpx !nqbze,nqibze,
     &     ,nblochpmx,ndummy1,ndummy2,ifcphi,is,nwp, !ifvcfpout,,mdimx,nbloch
     &     ifepscond,ifvxcpout,ifgb0vec
     &     ,nw0,iw,ifinin,iw0,noccxv,noccx
     &     ,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &     ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &     ,kx,isf,kqxx,kp,job,nwmax !,ifev1,ifev2 !,nhwtot
     &     ,ihis,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &     ,  ifchipmlog ,   nw_w,nwmin  !,ngpmx,ngcmx
      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm, dwry,dwh,omg2, q(3),  qgbin(3),qx(3) 
      real(8):: ua=1d0          ! this is a dummy.
      integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2), ndble=8, nword
      integer,allocatable ::  iqib(:),ifppb(:) !,lx(:) ngvecc(:,:),ngveccB(:,:),
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &     zw(:,:),zw0(:,:), zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: eqt(:), 
     &     aaa(:,:),symope(:,:),
     &     ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:)!, qbze(:,:)!,qibze(:,:) !,ecore(:,:)
c     &  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: rcxq(:,:,:,:)
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      real(8) ::qbzx(3)
      logical :: debug=.false.
c      integer,allocatable:: ibasf(:)
      logical :: realomega, imagomega 
      complex(8),allocatable:: zzr(:,:)
      complex(8) :: epxxx,vcmean
      character*9 fileps
      character*15 filepsnolfc
c      logical :: paralellx0=.true. !, hist
      character(5) :: charnum5
      character(20):: xxt
      real(8) :: Emin, Emax      ,emax2,emin2
c     integer :: iSigma_en  !sf..21May02  !iSigma_en is integer
                                !parameter stored in GWIN_V2
                                !which determines approximation for  self-energy.
                                !Self-energy should be made hermitian for energies to be real
cxxx  !iSigma_en==0 SE_nn'(ef)+img integral:delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
cxxx  !iSigma_en==1 SE_nn'(ef)+delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
                                !iSigma_en==2 [SE_nn'((e_n+e_n')/2)+h.c.]/2
                                !iSigma_en==3 [(SE_nn'(e_n)+SE_nn'(e_n'))/2+h.c.]/2
      real(8) :: omg2max,omg1max,wemax
      logical::imagonly=.false. , noq0p !,readgwinput
      integer::nwin, incwfin, verbose,ifpomat,nnmx,ikpo,nn_,noo,iqxxx,nomx
c      real(8)::efin
      logical :: nolfco=.false.
      integer:: isp1,isp2, ngc!,mrecg ! bzcase,
      real(8)::  quu(3),deltaq(3),qqq(3)=0d0 !
      complex(8),allocatable:: wgt(:,:,:)
      logical :: qbzreg 
!    logical ::smbasis !smbasis will be implemented in m_zmel.f which generates <phi|phi M>
      real(8):: q_r(3)
      complex(8),allocatable:: pomat(:,:)
      logical   :: timereversal,onceww
      integer :: jpm,ncc
      real(8) :: frr !, sciss
      integer :: idummy,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig,nmbas1x !ifepstinv,
      complex(8),allocatable:: epstinv(:,:),epstilde(:,:) !,zcousqrsum(:,:,:),zcousqr(:,:),eemat(:,:),zcousq0(:,:)
      real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq
!! Eq.(40) in PRB81 125102
c     complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:), 
     &     w_k(:),w_ks(:),w_kI(:), w_ksI(:) 
      complex(8),allocatable:: llw(:,:), llwI(:,:),aaamat(:,:)
      integer:: lxklm,nlxklm,ifidmlx,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
      complex(8):: vc1vc2
      real(8),allocatable:: aik(:,:,:,:)
      integer,allocatable:: aiktimer(:,:)
      integer:: nmbas_in , iqxendx,imb2 !iqqv,l2nl, 
      logical:: eibz4x0,tiii,iprintx=.false.,chipm=.false.,iqinit,localfieldcorrectionllw
      real(8)::hartree,rydberg,pi
c      integer :: iqeibz
      complex(8):: epslfc, axxx(10)
      integer:: src,dest
      integer:: ifw0w0i
      logical :: symmetrize,eibzmode
      real(8):: schi=-9999 !dummy
      integer:: i_reduction_npm, i_reduction_nwhis,  i_reduction_nmbas2
      logical:: crpa
      integer:: ibasx,ificlass,ifile_handle,ifiq0p
      logical:: readw0w0itest=.false. !tetra ,
      integer::nw_ixx,nwxx

      logical:: w4pmode
      complex(8),allocatable:: wmu(:,:),wmuk(:,:)
      integer:: ifw4p,ngbq0,igb
      real(8):: qv(3,3)

      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

      logical:: hx0
      integer:: imbas1,imbas2,ipm,niw_,ifif
!-------------------------------------------------------------------------
      call MPI__InitializeQSPBM()
      call MPI__consoleout('hx0fp0_sc')
      call cputid (0)
      if(verbose()>=100) debug=.true.
      write(6,*) ' --- hx0fp0_sc Choose omodes below ----------------'
      write(6,*) '  ixc= 11,10011,or 1011 '
      write(6,*) ' --- Put number above ! -----------------'
      if( MPI__root ) then
         read(5,*) ixc !c     call readin5(ixc,iqxini,iqxend)
      endif
      call MPI__Broadcast(ixc)
      if(MPI__root) iprintx=.true.
      crpa = .false.
      if(ixc==11) then
         write(6,*) " OK ixc=11 normal mode "
      elseif(ixc==10011) then
         write(6,*) " OK ixc=10011 crpa mode "
         crpa=.true.
      elseif(ixc==1011) then
         write(6,*) 'OK ixc=1011 Add W0W0I part at q=0'
      else
         write(6,*)'we only allow ixc==11. given ixc=',ixc
         call Rx( 'error: give allowed arg for hx0fp0_sc.')
      endif

      allocate( zzr(1,1))       !dummy
      hartree= 2d0*rydberg()
      pi     = 4d0*datan(1d0)
      fourpi = 4d0*pi
      sqfourpi= sqrt(fourpi)
      
!!  newaniso2=T is fixed now
      call Read_BZDATA(hx0) !Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
c      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp
!! === Readin by genallcf ===
!! See "use m_genallcf_v3" at the begining of this routine
      call Genallcf_v3(incwfx=0) ! incwfin= 0 takes 'ForX0 for core' in GWinput
!!!!  CAUTION: WE ASSUME iclass(iatom)= iatom,  nclass = natom.  !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) call Rx( ' hx0fp0_sc: nclass /= natom ')
      tpioa=2d0*pi/alat
      call Readefermi() !Readin EFERMI
      call Readhbe()    !Read dimensions 
      call ReadGWinputKeys() !Readin dataset in GWinput
c      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and GWIN_V2')
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
      if(verbose()>=40) write (*,*)' hsfp0.sc.m.F: end of mptauof_x'
!! Rdpp gives ppbrd: radial integrals and cgr = rotated cg coeffecients.
c      call Rdpp(nl, ngrpx, nn, nclass, nspin, symope,qbas)
c      call Rdpp(ngrpx,symope) !moved to m_zmel
!! Set itq in m_zmel
      call Setitq()
!!
!! ... initialization of readEigen !readin m_hamindex
      call Init_readeigen() !nband,mrece)!EVU EVD are stored in m_readeigen
      call Init_readeigen2() !mrecb,nlmto,mrecg)

!! Getfreq gives frhis,freq_r,freq_i, nwhis,nw,npm 
      realomega = .true.
      imagomega = .true.
      call Getfreq2(.false.,realomega,imagomega,ua,iprintx)!tetra,
!! Write freq_r
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
      nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB
      iqxini = 1 
      iqxend = nqibz + nq0i + nq0iadd ! [iqxini:iqxend] range of q points.
      nwp = nw+1              ! energy mesh along real zxis 
      nprecx = ndble          ! We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()!record size
      nspinmx = nspin
      allocate( zw(nblochpmx,nblochpmx) )
      if(MPI__root) then  ! I think it is not so meaningful to give a subroutine for writing WV.d.
        open(newunit=ifwd, file='WV.d', action='write')
        write (ifwd,"(1x,10i14)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
        close(ifwd)
      endif
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2
      ! multitet=T ==> micro tetrahedron method (divided-tetrahedron). Not used so much now...

!! Array for eigenvalues.
      allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))

!! EIBZ mode
      eibzmode = eibz4x0()
      call Seteibz(iqxini,iqxend,iprintx)
!!
      allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
      llw = 1d99
      llwI= 1d99
      if(ixc==1011)then !ixc==11 is a debug mode to test contrib. at \Gamma point.
         goto 1191
      endif
      
!! w4phonon. all nodes have wmu array. !still developing mode
      w4pmode=.false.
      if(sum(ixyz)/=0) w4pmode=.true.
      if(w4pmode) then
        allocate( wmuk(2:nblochpmx,3))
        wmuk=1d99
      endif
!! Get ngbq0 (for q=0) and broadcast for w4p
      if( MPI__root.and. w4pmode ) then
        q = (/0d0,0d0,0d0/)
        call readqg0('QGcou', q,  quu,ngc)
        ngbq0 = nbloch+ngc
      endif
      call MPI__Broadcast(ngbq0)

!! Rank divider
      call MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
      call MPI__hx0fp0_rankdivider2S(nspinmx)
      write(6,'("irank=",i5," allocated(MPI__qtask)=",L5)')MPI__rank,allocated(MPI__qtask)
      do iq = iqxini,iqxend
        if(MPI__qtask(iq)) write(6,'("irank iq=",i5,i5)') MPI__rank,iq
      enddo
      
!! == Calculate x0(q,iw) and W == main loop 1001 for iq. 
!! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
!! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
!! (or need to modify do 1001 loop).
!! ---------------------------------------------------------------
!! === do 1001 loop over iq ============================================
!! ---------------------------------------------------------------

!! main loop ========================================================
      if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
      do 1001 iq = iqxini,iqxend
        if( .not. MPI__Qtask(iq) ) cycle
        call cputid (0)
        q = qibze(:,iq)
        write(6,"('do 1001: iq q=',i5,3f9.4)")iq,q
        if(MPI__rootS) then
!! open WVR (W-V along real axis), WVI (along imag axis) for given q (iq).
         open(newunit=ifrcwi,file='WVI.'//charnum5(iq),form='unformatted',access='direct',recl=mrecl)
         open(newunit=ifrcw, file='WVR.'//charnum5(iq),form='unformatted',access='direct',recl=mrecl)
        endif
        call Readqg0('QGcou',q,   quu,ngc) ! ngc: the number of IPW for the interaction matrix (in QGcou),
                                           ! quu: equivalent to q.
        ngb = nbloch + ngc    ! number of MPB = PB part + IPW part
                              ! ngc is shown at the end of console output of qg4gw.
        write(6,*)'nbloch ngb ngc=',nbloch,ngb,ngc
!! Read Coulomb matrix 
        call Readvcoud(q,iq,ngb) ! Readin vcousq,zcousq for the Coulomb matrix
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
        call Setppovlz(q,ngc,ngb,matz=.not.eibz4x0())
!! ----------------------
!! === loop over spin=== 
!! ----------------------
        allocate( rcxq(nmbas1,nmbas2,nwhis,npm)) 
        rcxq = 0d0
        do 1003 is = MPI__Ss,MPI__Se !is=1,nspinmx
           write(6,"(' ### ',2i4,' out of nqibz+n0qi+nq0iadd nsp=',2i4,' ### ')") 
     &          iq, is, nqibz + nq0i+nq0iadd, nspin
           if(debug) write(6,*)' niw nw=',niw,nw
           isf = is
           do kx = 1, nqbz
              call readeval(qbz(:,kx),   is,  ekxx1(1:nband, kx) ) !read eigenvalue
              call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
           enddo
           call gettetwt(q,iq,is,isf,nwgt(:,iq),ekxx1,ekxx2,nband,eibzmode) ! Tetrahedron weight. 
           call x0kf_v4hz(
     i          q, is,isf,  iq,ngb,ngc,
     o          rcxq,           !  rcxq is the accumulating variable for spins 
     i          nmbas_in, 
     i          eibzmode, 
     i          nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq),eibzsym(:,:,iq), crpa) 
           call tetdeallocate()
 1003   continue
!! Symmetrize and convert to Enu basis.
!!   That is, we get dconjg(tranpsoce(zcousq))*rcxq*zcousq for eibzmode
        if(eibzmode)  then
          call x0kf_v4hz_symmetrize(!npm,!ncc,   
     i     q,iq,ngb,ngc, 
     o     rcxq,                ! rcxq is the accumulating variable for spins 
     i     nolfco,zzr,nmbas_in, 
     i     chipm,eibzmode,      
     i     eibzsym(:,:,iq)) 
        endif
!! Reduction rcxq in the S-axis
        write(6,*) 'MPI__AllreduceSumS start'
        do i_reduction_npm=1,npm
        do i_reduction_nwhis=1,nwhis
        do i_reduction_nmbas2=1,nmbas2
           call MPI__AllreduceSumS(
     .          rcxq(1,i_reduction_nmbas2,i_reduction_nwhis,i_reduction_npm), nmbas1)
        enddo
        enddo
        enddo
        write(6,*) 'MPI__AllreduceSumS end'
!! Hilbert transform.  Genrerate Real part from Imaginary part. ======
        if(allocated(zxq) ) deallocate(zxq,zxqi)
        allocate( zxq(nmbas1,nmbas2,nw_i:nw), zxqi(nmbas1,nmbas2,niw) )
        write(6,'("goto dpsion5: nwhis nw_i niw nw_w nmbas1 nmbas2=",6i5)') nwhis,nw_i,nw,niw,nmbas1,nmbas2
        call dpsion5(           !frhis,nwhis, freq_r, nw, freq_i,niw,
     i   realomega, imagomega, 
     i   rcxq, nmbas1,nmbas2, ! rcxq is alterd---used as work npm,nw_i, 
     o   zxq, zxqi,
     i   chipm, schi,is,  ecut,ecuts)
        if(allocated(rcxq) ) deallocate(rcxq)
        write(6,*)' --- end of dpsion5 ----',sum(abs(zxq)),sum(abs(zxqi))
!!        
        allocate( zw0(ngb,ngb) ) 
!! ===  RealOmega ===
        if (realomega) then
          if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
          nwmax = nw
          nwmin = nw_i
!! prepare for iq0.
          iq0 = iq - nqibz
          allocate(epstilde(ngb,ngb))
          allocate(epstinv(ngb,ngb))
          write(6, *)" === trace check for W-V === nwmin nwmax=",nwmin,nwmax
          do 1015 iw  = nwmin,nwmax
            frr= dsign(freq_r(abs(iw)),dble(iw))
            imode = 1
            if(iq<=nqibz) then  !for mmmw
              if(iq==1) then
                ix=1
                zw0(:,1)=0d0
                zw0(1,:)=0d0
              else
                ix=0
              endif
!!  Eqs.(37),(38) in PRB81 125102 (Friedlich)
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                enddo
              enddo
              epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
              call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))

!! w4p writing eps
              if(iw==0.and.w4pmode) then 
                !static epstinv is saved. For q=0 epstilde (mu=1 skipped). For q/=0 full matrix inversion.
                             !(ix=1 is set for q=0) 
                ifw4p = ifile_handle()
                open(ifw4p,file='W4PHONON.'//charnum5(iq),form='unformatted')
                write(ifw4p) iq,q,ngb,ix !ix=0, or ix=1 for q=0 (iq=1)
                write(ifw4p) epstinv(ix+1:ngb,ix+1:ngb) 
                close(ifw4p)
              endif  

c$$$  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$  cmmm direct inversion vs. block inversion
c$$$  if(iq>nqibz) then
c$$$  c direct inversion
c$$$  ix=0
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  do igb1=1+ix,ngb
c$$$  do igb2=1+ix,ngb
c$$$  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
c$$$  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
c$$$  enddo
c$$$  enddo
c$$$  c              write(6,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$  write(6,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,1d0/epstinv(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$  c             write(6,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(1,2:10:3),zw0(1,63:70:3)
c$$$  c block inversion
c$$$  ix=1
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  absq=sqrt(sum(q**2*tpioa**2))
c$$$  sk(  1:ngb)= zxq(1,1:ngb,iw)
c$$$  sks( 1:ngb)= zxq(1:ngb,1,iw)
c$$$  w_k(1) =0d0
c$$$  w_ks(1)=0d0
c$$$  w_k( 2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(vcousq(1)*sk(2:ngb)*vcousq(2:ngb),epstinv(2:ngb,2:ngb))
c$$$  w_ks(2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(epstinv(2:ngb,2:ngb),vcousq(1)*sks(2:ngb)*vcousq(2:ngb))
c$$$  llw(iw,iq0)=
c$$$  &             1d0
c$$$  &            -vcousq(1)*sk(1)*vcousq(1) ! sk(1,1,iw)=sks(1,1,iw)=H of Eq.(40).
c$$$  &            -vcousq(1)*vcousq(1)* sum( vcousq(2:ngb)*sk(2:ngb) * matmul(epstinv(2:ngb,2:ngb),sks(2:ngb)*vcousq(2:ngb)))
c$$$  write(6,"('mmmmzwp99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0), !(1d0/llw(iw,iq0)-1d0)*vcousq(1)**2,
c$$$  c     &                  w_k(2:10:3)/llw(iw,iq0), w_k(63:70:3)/llw(iw,iq0)
c$$$  &                  w_ks(2:10:3)/llw(iw,iq0), w_ks(63:70:3)/llw(iw,iq0)
c$$$  write(6,"('mmmmzwp99x ')")
c$$$  endif
c$$$  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              do igb1=1+ix,ngb
                do igb2=1+ix,ngb
                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                enddo
              enddo
c$$$                     if(iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
              zw(1:ngb,1:ngb) = zw0
              if (MPI__rootS)then
                write(ifrcw, rec= iw-nw_i+1 ) zw !  WP = vsc-v
              endif
              call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
            endif
           
            if(iq>nqibz) then
!! Full inversion to calculalte eps with LFC.
              vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
              if(localfieldcorrectionllw()) then
                ix=0
                do igb1=ix+1,ngb
                  do igb2=ix+1,ngb
                    if(igb1==1.and.igb2==1) then
                      epstilde(igb1,igb2)= 1d0 - vcou1*zxq(1,1,iw)
                      cycle
                    endif
                    epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                    if(igb1==igb2) then
                      epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                    endif   
                  enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                if(iq0<=nq0i) llw(iw,iq0)= 1d0/epstinv(1,1)
!! Wing elements calculation july2016
!! We need check nqb is the same as that of q=0
                if(ixyz(iq0)/=0.and.iw==0) then
                   if(ngb/=ngbq0) then
                      write(6,*)q,iq0,ngb,ngbq0
                      call rx('hx0p0_sc: ngb/=ngbq0')
                   endif
                   wmuk(2:ngb,ixyz(iq0))=epstinv(1,2:ngb)/epstinv(1,1)
                   ! this is dot(q(:)*w_mu(:,igb)). See PRB125102(2016) eq.(36)
                endif
              else
                if(iq0<=nq0i) llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw) 
              endif 
              if(iq0<=nq0i) write(6,"('iq iw_R omg(iw) eps(wFC) eps(woLFC) ',
     &         2i5,x,10(d13.6,2x,d13.6,x,d13.6,2x,d13.6,x,d13.6))")
     &         iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
            endif
 1015     continue              !iw
          if( allocated(zzr) ) deallocate(zzr)
        endif 
!! === RealOmega end ===

        
!! === ImagOmega ===
        if (imagomega) then
          write(6,*)' goto imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
          imode=1
          do 1016 iw  = 1,niw
c                  if( newaniso2 .and. iq<=nqibz ) then
            if( iq<=nqibz ) then
!!  Eqs.(37),(38) in PRB81 125102
              if(iq==1) then
                ix=1
                zw0(:,1)=0d0
                zw0(1,:)=0d0
              else
                ix=0
              endif
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                enddo
              enddo
              epstinv=epstilde
              call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                enddo
              enddo
              zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
              if (MPI__rootS) then
                write(ifrcwi, rec= iw)  zw !  WP = vsc-v
              endif
              call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
            endif

            if(iq>nqibz) then
!! Full inversion to calculalte eps with LFC.
              vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
              if(localfieldcorrectionllw()) then
                ix=0
                do igb1=ix+1,ngb
                  do igb2=ix+1,ngb
                    if(igb1==1.and.igb2==1) then
                      epstilde(igb1,igb2)= 1d0 - vcou1*zxqi(1,1,iw)
                      cycle
                    endif
                    epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                    if(igb1==igb2) then
                      epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                    endif   
                  enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                if(iq0<=nq0i) llwI(iw,iq0)= 1d0/epstinv(1,1)
              else
                if(iq0<=nq0i) llwI(iw,iq0)=  1d0 -vcou1*zxqi(1,1,iw) 
              endif  
              if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)")
     &         iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
            endif

 1016     continue
          deallocate(epstinv)
          if(allocated(epstilde)) deallocate(epstilde)
        endif 
!! === ImagOmega end ===
        if(allocated(vcoul)) deallocate(vcoul)
        if(allocated(zw0)) deallocate(zw0)
        if(allocated(zxq )) deallocate(zxq)
        if(allocated(zxqi)) deallocate(zxqi)
        if (MPI__rootS) then
          close(ifrcwi)
          close(ifrcw)
        endif
 1001 continue
c============end of loop over q point =================================
c=======================================================================
      call MPI__barrier()


      
      
!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
      if(MPI__size/=1) then
        do iq=nqibz+1,iqxend
          iq0 = iq - nqibz
          if(MPI__Qranktab(iq)/=0) then !jan2012
            if(MPI__Qranktab(iq) == MPI__rankQ) then
              dest=0
              if(iq0<=nq0i) then
                 call MPI__DbleCOMPLEXsendQ(llw(nw_i,iq0),(nw-nw_i+1),dest)
                 call MPI__DbleCOMPLEXsendQ(llwI(1,iq0),niw,dest)
              endif
              if(ixyz(iq0)/=0) then
                call MPI__DbleCOMPLEXsendQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,dest)
              endif
            elseif(MPI__rootQ) then
              src=MPI__Qranktab(iq)
              if(iq0<=nq0i) then
                 call MPI__DbleCOMPLEXrecvQ(llw(nw_i,iq0),(nw-nw_i+1),src)
                 call MPI__DbleCOMPLEXrecvQ(llwI(1,iq0),niw,src)
              endif
              if(ixyz(iq0)/=0) then
                call MPI__DbleCOMPLEXrecvQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,src)
              endif
            endif
          endif
        enddo  
      endif

c commentout block inversion
c$$$!! Add LFC (local field correction) to llw and llwI
c$$$         if(newaniso2 .and. MPI__rank == 0 ) then ! only on root node
c$$$            iq=1 !for q=0
c$$$            vcoudfile='Vcoud.'//charnum5(iq)
c$$$            ifvcoud = iopen(trim(vcoudfile),0,-1,0)
c$$$            read(ifvcoud) ngb0
c$$$            read(ifvcoud) qvv
c$$$            if(sum(abs(qvv))>1d-10) then
c$$$               write(6,*)'qvv =',qvv
c$$$               stop 'hx0fp0: qvv/=0 hvcc is not consistent'
c$$$            endif
c$$$            if(allocated(zcousq0)) deallocate( zcousq0,vcousq0 )
c$$$            allocate( zcousq0(ngb0,ngb0),vcousq0(ngb0))
c$$$            read(ifvcoud) vcousq0
c$$$            read(ifvcoud) zcousq0
c$$$            idummy=iclose(trim(vcoudfile))
c$$$            vcousq=sqrt(vcousq)
c$$$            allocate(epstinv(ngb0,ngb0),w_k(ngb0),w_ks(ngb0),w_kI(ngb0),w_ksI(ngb0),eemat(ngb0,ngb0))
c$$$
c$$$            do iq0=1,nq0i
c$$$              iq = iq0 + nqibz
c$$$              q = qibze(:,iq)
c$$$
c$$$              vcoudfile='Vcoud.'//charnum5(iq)
c$$$              ifvcoud = iopen(trim(vcoudfile),0,-1,0)
c$$$              read(ifvcoud) ngb
c$$$              read(ifvcoud) qvv
c$$$              if(sum(abs(qvv-q))>1d-10) then
c$$$               write(6,*)'qvv =',qvv
c$$$               stop 'hx0fp0: qvv/=0 hvcc is not consistent'
c$$$              endif
c$$$              if(allocated(zcousq)) deallocate(zcousq)
c$$$              if(allocated(vcousq)) deallocate(vcousq)
c$$$              allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
c$$$              read(ifvcoud) vcousq
c$$$              read(ifvcoud) zcousq
c$$$              idummy=iclose(trim(vcoudfile))
c$$$              vcousq=sqrt(vcousq)
c$$$
c$$$              ifepstinv = iopen('EPS0inv',0,0,0)
c$$$              read(ifepstinv) ngb
c$$$
c$$$               ngc=ngb-nbloch
c$$$               if(allocated(ppovlz)) deallocate(ppovlz)
c$$$               if(allocated(ppovl)) deallocate(ppovl)
c$$$               allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
c$$$               call readppovl0(q,ngc,ppovl) !q was qq
c$$$               ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
c$$$               ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
c$$$
c$$$!  eemat: Z\mu_i(\bfk=0)^* <i|j> Z\nu_j(\bfk) 
c$$$               eemat =matmul(transpose(dconjg(zcousq0)),matmul(ppovlz,zcousq))
c$$$               vcou1  = fourpi/sum(q**2*tpioa**2) ! test-->vcousq(1)**2 !fourpi/sum(q**2*tpioa**2-eee)
c$$$               vcou1sq = vcou1**.5
c$$$               write(6,*)
c$$$
c$$$              do iw=nwmin,nwmax
c$$$                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$                epstinv(2:ngb,2:ngb) = matmul( transpose(dconjg(eemat(2:ngb,2:ngb))),
c$$$     &                                matmul(epstinv(2:ngb,2:ngb),eemat(2:ngb,2:ngb)) )
c$$$                if(iw/=iwx) then
c$$$                write(6,*)'iw iwx=',iw,iwx
c$$$                stop 'hx0fp0_sc: iw/=iwx'
c$$$                endif
c$$$                w_k(2:ngb) = vcou1sq*matmul( epstinv(2:ngb,2:ngb), sk(2:ngb,iw,iq0)*vcousq(2:ngb))
c$$$                epslfc = -vcou1sq*sum( sks(2:ngb,iw,iq0) * w_k(2:ngb) *vcousq(2:ngb) )
c$$$                llw(iw,iq0) = llw(iw,iq0)  + epslfc
c$$$                write(6,"('eps(on real) iq iw',2i4,2f9.3,2x,2f9.3)") iq0,iw, llw(iw,iq0)-epslfc,llw(iw,iq0)
c$$$              enddo
c$$$              do iw=1,niw
c$$$                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$                if(iw/=iwx) then
c$$$                 write(6,*)'iw iwx=',iw,iwx
c$$$                 stop 'hx0fp0_sc: iw/=iwx'
c$$$                endif
c$$$                w_kI(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), skI(2:ngb,iw,iq0)*vcousq(2:ngb))
c$$$                epslfc=- vcou1sq*sum( sksI(2:ngb,iw,iq0)* w_kI(2:ngb)*vcousq(2:ngb) )
c$$$                llwI(iw,iq0)= llwI(iw,iq0)+epslfc 
c$$$                write(6,"('eps(on img ) iq iw',2i4,2f9.3,2x,2f9.3)")iq0,iw, llwI(iw,iq0)-epslfc,llwI(iw,iq0)
c$$$              enddo
c$$$              ifepstinv = iclose('EPS0inv')
c$$$           enddo
c$$$         endif


      
!! == W(0) divergent part and W(0) non-analytic constant part.==
 1191 continue

      if(MPI__rank == 0 ) then  ! MIZUHO-IR only on root node
!! ix=1011 is a special mode to overwrite llw and llwI for test purpose
!! A file W0W0I is  generated by call w0w0i, but usually unused at anywhere.
        if(ixc==1011) then  
          open(newunit=ifw0w0i,file='W0W0I',form='unformatted')
          read(ifw0w0i) nw_ixx,nwxx,niw_,nq0ix
          write(6,*)'w0w0i: n=',nw_ixx,nwxx,niw_,nq0ix
          if(nq0i/=nq0ix)  call rx('nq0i/=nq0ix')
          if(nw_i/=nw_ixx) call rx(nw_i/=nw_ixx)
          if(nw/=nwxx) call rx(nw/=nwxx)
          read(ifw0w0i) llw(nw_i:nw,1:nq0i)
          read(ifw0w0i) llwI(1:niw_,1:nq0i)
c            read(ifw0w0i) w0(nw_i:nw)
c            read(ifw0w0i) w0i(1:niw)
          close(ifw0w0i)
        endif  
!! get w0 and w0i (diagonal element at Gamma point)
!! This return w0 and w0i. (llw and llwi are input)
!! Outputs w0,w0i,llmat. See use m_w0w0i at the begining of this routine.
        call w0w0i(llw,llwI,nw_i,nw,nq0i,niw,q0i)  !all inputs. get effective W0,W0i, and L(omega=0) matrix.
!! Finalize w4phonon
!! wmuk(ix)= matmul(wmu,qv) ==> wmu= matmul(wmuk,qvinv)
        if(w4pmode) then
          do i=1,3
            qv(:,i)= tpioa*q0i(:,ixyz(i))
            qv(:,i)= qv(:,i)/sqrt(sum(qv(:,i)**2))
          enddo
          call matinv(3,qv)
          allocate( wmu(2:ngbq0,3) )
          do igb=2,ngbq0
            wmu(igb,:) =matmul(wmuk(igb,:),qv)
          enddo
          ifw4p = ifile_handle()
          open(ifw4p,file='W4PHONON.HeadWing',form='unformatted')
          write(ifw4p) llmat(1:3,1:3),ngbq0 !for q~0
          write(ifw4p) wmu(2:ngbq0,1:3) !for q~0
          close(ifw4p)
          deallocate(wmu,wmuk)
        endif

!! Read WVR and WVI at Gamma point, and give correct W(0) (averaged in the Gamma cell, where
!! Gamma cell) is the micro cell of BZ including Gamma point).
!! === w0,w0i are stored to zw for q=0 ===
!! === w_ks*wk are stored to zw for iq >nqibz ===
! We assume iq=1 is for rank=0
        do iq = 1,1             !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
          q = qibze(:,iq)
          do ircw=1,2
            if    (ircw==1) then
              nini=nw_i
              nend=nw
              open(newunit=ifrcwx,  file='WVR.'//charnum5(iq), form='unformatted',
     &          status='old',access='direct',recl=mrecl)
            elseif(ircw==2) then;  nini=1;      nend=niw;
              open(newunit=ifrcwx,  file='WVI.'//charnum5(iq), form='unformatted',
     &          status='old',access='direct',recl=mrecl)
            endif
            do iw=nini,nend
              read(ifrcwx, rec= iw-nini+1 ) zw !(1:ngb,1:ngb)
              if( iq==1 ) then
                if(ircw==1) zw(1,1) = w0(iw)
                if(ircw==2) zw(1,1) = w0i(iw)
              endif
              write(ifrcwx,rec=iw-nini+1) zw !(1:ngb,1:ngb)
           enddo
           close(ifrcwx)
          enddo
        end do
      endif
      call cputid(0)
      write(6,*) '--- end of hx0fp0_sc --- irank=',MPI__rank
      call flush(6)
      call MPI__Finalize
      if(ixc==11) call rx0( ' OK! hx0fp0_sc ixc=11 Sergey F. mode')
      if(ixc==1011) call rx0( ' OK! hx0fp0_sc ixc=1011 W0W0Ionly')
      end program hx0fp0_sc 


C===================================================================
      subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
      implicit none
      integer:: nblochpmx,nbloch,ngb,iw,i,iq
      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
      real(8):: freqq
      character*(*)::tagname
      trwv=0d0
      do i = 1,nbloch
        trwv = trwv + zw(i,i)
      enddo
      trwv2 = 0d0
      do i = 1,ngb
         trwv2 = trwv2 + zw(i,i)
      enddo                     !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
      write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
c     do i = 1,ngb
c     write(6,'("iii i=",i4,a,f10.4,2i5,4d22.14)')i,tagname,freqq,iq,iw,zw(i,i)
c     enddo
      end


