      module m_bzintegration
      use m_addrbl,only: swtk
      use m_subzi,only: rv_a_owtkb !rewritten! need to protect them somehow
      contains
      subroutine m_bzintegration_init(evlall, eferm,sev,sumqv,vnow)!,goto99)
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq
      use m_lmfinit, only: lso,nsp,ham_scaledsigma,nlibu,lmaxu,bz_w,
     & lmet=>bz_lmet,stdo,nbas,epsovl=>ham_oveps,nspc,bz_n,bz_fsmommethod,qbg=>zbak,fsmom=>bz_fsmom,ndos=>bz_ndos
      use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
      use m_lmfinit,only: ham_ldham , nchan=>pot_nlma, nvl=>pot_nlml
      use m_qplist,only: m_qplist_init, qplist,nkp,xdatt,labeli,labele,dqsyml,etolc,etolv,
     &     nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml
      use m_mkpot,only: qval
c     o      vrmt, sv_p_osig, sv_p_otau, sv_p_oppi 
c     o  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
c     o     , vval_rv , fes1_rv 
      use m_subzi, only: lwtkb,lswtk
      use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx,nspx=>ham_nspx 
      use m_ext,only: sname
      implicit none
      logical:: lfill=.false.,ltet
      logical:: debug,cmdopt0 !goto99,
      integer:: ierr,ifimag,i,ifi,ifile_handle,iobzwt
      real(8):: dosrng,evlall(*),sev,sumqv(3,*),eferm,vnow,ef0,bz_ef
      real(8),parameter::    NULLR =-99999
c     goto99=.false.
      sev=0d0
      ltet = ntet>0
c      qbg = ctrl_zbak(1) !homogenious background charge
      debug    = cmdopt0('--debugbndfp')
      if (master_mpi) ierr=unlink('MagField') !delete
!! --- BZ integration for fermi level, band sum and qp weights ---
c      if(lmet>=0) then
         dosrng = 8
         if (bz_n .lt. 0) dosrng = 16
         if( debug) print *, 'tttttttttt 99999999 6 call bzwtsf'
         if(bz_fsmommethod == 1) then !takao dec2010
!     ! vnow june22013 !vnow !(in Ry) contains magnetic field
!     ! For eigenvalus, add  -vnow/2 for isp=1, and +vnow/2 for isp=2.
            call bzwtsf2 ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval-qbg, !note qval is output
     .           fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .            , dosrng , rv_a_owtkp , evlall ,  lswtk , swtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill,vnow)
         else
!     ! vnow june22013
            call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval-qbg , 
     .           fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .           , dosrng , rv_a_owtkp , evlall , lswtk , swtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill, vnow)
         endif
!     ! june2013 magfield is added
         if(fsmom/=NULLR .and. master_mpi) then
            open(newunit=ifimag,file='MagField',status='unknown')
            write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow
            close(ifimag)
         endif
!     !         Store val charge & magnetic moment in sumqv(1..2)
c         if (lmet .ne. 4) then
            sumqv(1,1) = sumqv(1,2)
            sumqv(2,1) = sumqv(2,2)
c         endif
c         if (lmet .ne. 4) then
c            ef0 = eferm
c            bz_ef=ef0
c         endif
         if (lmet > 0) then
            if (master_mpi) then
               ifi=ifile_handle()
               open(ifi,file='wkp.'//trim(sname),form='unformatted') 
               i = iobzwt ( 0 , ndhamx , nkp , nspx , eferm , rv_a_owtkb , -ifi )
               close(ifi) !call fclr('wkp',ifi)
            endif
         endif
!     !
c         write(6,*) 'ssssssssss sev=',sev
c      endif
      end subroutine
      end module
