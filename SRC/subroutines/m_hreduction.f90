module m_hreduction
! use m_cmdopt_registry, only: c0_gs, c0_mlo_ortho, c0_mlo_orthonorm !the options were removed 2026-10-02
contains
  !> How many of the lowest PMT eigenstates at this k are NOT model-like: the
  !  count of leading states whose weight sum_j |<Psi^PMT_i|Psi^MTO_j>|^2 in the
  !  model subspace is below 1/2 (semicore local orbitals, O 2s bands, ...).
  !  m_HamPMT takes the MINIMUM of this over all k (and spins, and MPI ranks)
  !  and hands it to Hreduction as nskip_auto, so the projector drops the same
  !  number of states at every k. Deciding it per k (the rule until 2026-09-18)
  !  flipped between 0 and 1 across k for Cu's d-only model, where the s band is
  !  band 1 at some k and not at others, and put kinks into the MLO bands.
  subroutine Hreduction_nskip(ndimPMT,hamm,ovlm,ndimMTO,ix, nskipin, evlout)
    use m_zhev,only:zhev_tk4
    use m_nvfortran, only: findloc
    use m_lmfinit,only:oveps
    implicit none
    integer,intent(in):: ndimPMT,ndimMTO,ix(ndimMTO)
    complex(8),intent(in):: hamm(ndimPMT,ndimPMT),ovlm(ndimPMT,ndimPMT)
    integer,intent(out):: nskipin
    real(8),intent(out),optional:: evlout(:)   ! the lowest size(evlout) PMT eigenvalues (for the gap check)
    integer:: i,j,nev,nmx
    real(8):: evlmto(ndimMTO),evl(ndimPMT),wsum(ndimPMT)
    complex(8):: evecmto(ndimMTO,ndimMTO),evecpmt(ndimPMT,ndimPMT),ovlmx(ndimPMT,ndimPMT),hammx(ndimPMT,ndimPMT)
    complex(8):: sv(ndimPMT,ndimMTO)
    real(8),parameter:: epscore=0.5d0
    ovlmx=ovlm; hammx=hamm; nmx=ndimMTO
    call zhev_tk4(ndimMTO,hammx(ix(1:ndimMTO),ix(1:ndimMTO)),ovlmx(ix(1:ndimMTO),ix(1:ndimMTO)), nmx,nev, evlmto, evecmto, oveps)
    ovlmx=ovlm; hammx=hamm; nmx=ndimPMT
    call zhev_tk4(ndimPMT,hammx,ovlmx, nmx,nev, evl,evecpmt, oveps)
    sv = matmul(ovlm(:,ix(1:ndimMTO)), evecmto)          ! S |Psi_MTO_j>
    wsum = 0d0
    do i=1,nev
       do j=1,ndimMTO
          wsum(i) = wsum(i) + abs(sum(dconjg(evecpmt(:,i))*sv(:,j)))**2
       enddo
    enddo
    nskipin = findloc(wsum(1:nev) > epscore, value=.true., dim=1) - 1
    if (nskipin < 0) nskipin = nev
    if (present(evlout)) then
       evlout = 1d99
       evlout(1:min(size(evlout),nev)) = evl(1:min(size(evlout),nev))
    endif
  end subroutine Hreduction_nskip

  subroutine Hreduction(mlomethod,iprx,ndimPMT,hamm,ovlm,ndimMTO,ix,fff1, hammout,ovlmout, qp, cmlo,nev, zMLO, nskip_auto, rnorm, lowdin) !> Reduce H(ndimPMT) to H(ndimMTO)
    ! cmlo= <Psi^MPT i|F^MLO j>
   use m_zhev,only:zhev_tk4
   use m_nvfortran, only: findloc
   use m_readqplist,only: eferm, ecbot
!   use m_HamPMT,only: GramSchmidt!,epsovl
   use m_lgunit,only:stdo
   use m_lmfinit,only:oveps
   use m_keyvalue,only: getkeyvalue
   use m_GWinput, only: gwinput_init, gwinput_loaded, &
                        tg_mlo_w => mlo_w, &   ! tg_mlo_nskip => mlo_nskip retired 2026-09-18
                        tg_mlo_emax => mlo_emax, &
                        tg_mlo_delta => mlo_delta, tg_mlo_wfrz => mlo_wfrz, &
                        tg_mlo_down => mlo_down !, tg_mlo_low => mlo_low, tg_mlo_wlow => mlo_wlow, &
                        !tg_mlo_pcut => mlo_pcut, tg_mlo_pw => mlo_pw   ! hidden cuts, commented out 2026-09-17
   implicit none
   integer::i,j,ndimPMT,ndimMTO,nx,nmx,ix(ndimMTO),nev,nxx,jj,ndimPMTx,nvpmt,mlomethod,nskip,nskipin
   real(8)::beta,emu,val,wgt(ndimPMT),evlmto(ndimMTO),evl(ndimPMT),evlx(ndimPMT),qp(3),eww,eadd
   complex(8):: evecmto(ndimMTO,ndimMTO),evecpmt(ndimPMT,ndimPMT)
   complex(8):: ovlmx(ndimPMT,ndimPMT),hammx(ndimPMT,ndimPMT),fac(ndimPMT,ndimMTO),ddd(ndimMTO,ndimMTO)
   complex(8):: hamm(ndimPMT,ndimPMT),ovlm(ndimPMT,ndimPMT)
   complex(8):: hammout(ndimMTO,ndimMTO),ovlmout(ndimMTO,ndimMTO)
   complex(8),optional,intent(out):: cmlo(ndimPMT,ndimMTO) !<Psi^PMT_i|F^MLO_k>, PMT-eigenstate-basis coefficients
   complex(8),optional,intent(out):: zMLO(ndimPMT,ndimMTO) !|F_MLO_k> = sum_m |chi^PMT_m> zMLO(m,k)
   integer,optional,intent(in):: nskip_auto  ! number of semicore-LO states to drop (k-independent, from m_HamPMT)
   logical,optional,intent(in):: lowdin ! Loewdin orthonormalization after the rnorm normalization (the standard MLO, 2026-10-02)
   real(8),optional,intent(in):: rnorm(ndimMTO) ! square integral of the real-space MLOs (HamRsMLO); F^MLO_i/sqrt(rnorm_i): one constant
                                                ! for all k, so the MLO in real space has norm 1 (2026-10-02)
   complex(8):: cmlo_loc(ndimPMT,ndimMTO) !internal working array
!   complex(8),optional:: zcplz(ndimPMT,ndimMTO)
   complex(8),allocatable :: Amat(:,:)
   real(8):: fff1,fff !epsovl=1d-8 epsovlm=0d0 ,
   logical:: iprx
   ! choosed MTO Hamiltonian. Get evlmto,evecmto
   ovlmx= ovlm
   hammx= hamm
   nmx = ndimMTO !  write(stdo,*)'Start Hreduction: 111'
   call zhev_tk4(ndimMTO,hamm(ix(1:ndimMTO),ix(1:ndimMTO)),ovlm(ix(1:ndimMTO),ix(1:ndimMTO)), nmx,nev, evlmto, evecmto, oveps)
   if(nev/=ndimMTO) call rx('Hreduction: nev/=ndimMTO We didnot get eigenfuncitons of ndimMTO. Linear dependency problem?')
   ! PMT Hamiltonian
   ovlm= ovlmx
   hamm= hammx
   nmx = ndimPMT !  write(stdo,*)'Start Hreduction: 222'
   call zhev_tk4(ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), nmx,nev, evl,evecpmt, oveps) !PMT
   ovlm=ovlmx
   ndimPMTx=nev !obtained. oveps may reduce ndimPMT to be ndimPMTx
   fac = (0d0,0d0)
   FacPMTMTO: block !fac(i,j) = <Psi_PMT_i|S|Psi_MTO_j> (rows i>nev stay 0), the matrix element Amat is made from (2026-09-27 18:55: two products)
     use m_blas,only: zmm_h, m_op_C
     complex(8):: sx(ndimPMT,ndimMTO), sv(ndimPMT,ndimMTO)
     integer:: istat
     sx = ovlmx(:,ix(1:ndimMTO))
     istat = zmm_h(sx, evecmto, sv, m=ndimPMT, n=ndimMTO, k=ndimMTO)                              ! S |Psi_MTO_j>
     istat = zmm_h(evecpmt, sv, fac, m=nev, n=ndimMTO, k=ndimPMT, opA=m_op_C, ldc=ndimPMT)        ! <Psi_PMT_i| S |Psi_MTO_j>
   endblock FacPMTMTO
   ModifyMatrixElements :block
      use m_nvfortran,only : findloc
      use m_ftox
      integer:: ie,nidxevlmto,nidxevl,ibx,jx,idxevlmto(ndimMTO),idxevl(ndimPMT),jbx,nval,nnn,imx,nbx,ii
      real(8):: eee,fffx,ecut,xxx,rydberg,facww,sss,fff,epscore,emax,alpha,emin,ww(ndimPMTx),dex,ddd,ewuse,efrz,ewfrz,dwin,wfrz,down !,ewcutf
      !real(8):: elow, ewlow, thlow, pcut, pw, thp, pwt(ndimPMTx)   ! hidden cuts (mlo_low / mlo_pcut), commented out
      !logical:: lowcut, charcut
      real(8),allocatable::mulfac(:,:),mulfacw(:,:)
      complex(8):: imag=(0d0,1d0)
      ! Normalization check.  sum_i |<Psi_PMT_i|Psi_MTO_j>|^2 = 1 requires the PMT
      ! eigenvectors to be complete; zhev_tk4 drops near-linearly-dependent
      ! directions (oveps), so a model that spans (nearly) the whole MTO block --
      ! e.g. mlo_lm + mlo_lm2 listing both radial sets -- loses a little weight
      ! there.  Report the worst case and renormalize when it is small; abort
      ! only when the loss is large enough to change the model.
      NormalizationCheck: block
        real(8):: dev, devmax
        integer:: jworst
        devmax=0d0; jworst=0
        do j=1,ndimMTO
          dev = sum(abs(fac(:,j))**2)-1d0
          if(abs(dev)>abs(devmax)) then; devmax=dev; jworst=j; endif
        enddo
        if(abs(devmax)>1d-2) then
          write(stdo,"(a,i5,a,f10.6)")' Hreduction: PMT completeness loss too large: band',jworst,' dev=',devmax
          call rxi('Hreduction: normalization error band index=',jworst)
        endif
        if(abs(devmax)>1d-4) then
          if(iprx) write(stdo,"(a,i5,a,es10.2,a)") &
            ' Hreduction: PMT completeness loss, worst band',jworst,' dev=',devmax,' -> renormalized'
          do j=1,ndimMTO
            dev = sum(abs(fac(:,j))**2)
            if(dev>1d-8) fac(:,j)=fac(:,j)/sqrt(dev)
          enddo
        endif
      endblock NormalizationCheck
      if(iprx) then
        do j=1,ndimMTO !Amat is corrected matrix element of <psi_PMT|psi_MTO>
          do i=1,ndimPMTx
            if(abs(fac(i,j))**2>.1) write(stdo,ftox)'fac matrix ',j,i,ftof(abs(fac(i,j))**2)
          enddo
        enddo
      endif
      
      ! Determine nskip: the lowest nskip PMT eigenstates (semicore local orbitals) are
      ! removed from the projector. The count comes from the basis (number of semicore
      ! LO functions, m_HamPMT%nsemicore) so it is the same at every k. The old rule
      ! -- count the lowest states whose weight in the model subspace is < 0.5 -- is
      ! k-dependent: for Cu's d-only model the s band is band 1 at some k and not at
      ! others, and the projector jumped between them (kinks in the MLO bands,
      ! Samples/MLOsamples/BackUp_notes/mlo_nskip_cu_problem.md). Kept as fallback
      ! when the caller cannot supply the count.
      epscore=0.5d0
      if (present(nskip_auto)) then
         nskipin = nskip_auto                       ! k-independent: min over k of the per-k count (m_HamPMT)
      else
         nskipin = findloc( sum(abs(fac(:,:))**2,dim=2) > epscore, value=.true.,dim=1)-1 ! per-k fallback (callers without a pre-pass)
      endif
      call gwinput_init()
      if (.not. gwinput_loaded) call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
      nskip = nskipin
      ! mlo_nskip (manual override) is retired 2026-09-18: the automatic rule is the definition now.
      !nskip = tg_mlo_nskip
      !if (nskip == -huge(0)) nskip = nskipin   ! sentinel ⇒ key absent
      write(stdo,ftox) 'nnnnn nskip',nskip !,ftof(sum(abs(fac(:,:))**2,dim=2))
      do i = 1, nskip   ! a dropped state that IS model-like would be a mistake; say so
         if (sum(abs(fac(i,:))**2) > epscore) write(stdo,ftox) &
              ' Hreduction: WARNING skipped PMT state',i,'has weight',ftof(sum(abs(fac(i,:))**2)),'in the model at q=',ftof(qp)
      enddo

      ! === Usage ===
      ! Simple version ---> Set mlo_method 0 with mlo_emax for Semiconductor (\lesssim VBM) or Al2O3_Cr (7eV or higher)
      !                     Not necessary for NiO, Ru2O3.  
      ! Generally speaking, we need emax for localized bands, while we need mlomethod2 for semiconductors (broad sp bands, smooth cutoff).
      ! 1. Only localized bands, I think no switch needed.
      ! 2. For semiconductors, set emax = Efermi (or even -9999) around ( ---> then mlomethod0 is close to mlomethod2).
      ! 3. For Al2O3_Cr (sp and d bandd), we need to set emax a little about the localized bands.
      ! --------------------
      !
      ! For sp bands smooth cutoff.
      !  Semiconductors: Si, GaAs,
      !    + We have to set mlo_emax 0 or something. For Al2O3_Cr, we need to set mlo_emax as 15 eV or so.
      !    (For Al2O3_Cr, we found mlomethod 1 with emax= 7 eV works well).
      !    
      !  NiO, Ru2O3 
      !    mlo_method 1 auto emax, or mlo_method 0 auto emax. auto emax 
      ! 
      !  Extract 3d or 4f bands
      !     + mlomethod 0 works.
      !     mlomethod 2 works for 4f extraction. No emax
      !
      !  We have to set mlo_emax, up to which we have to include i for fac=<Psi_PMT(i)|Psi_MTO(j)> for semiconductors or broad band included.
      !

      ! When P = \sum_i \sum_j |Psi^PMT_i> <Psi^PMT_i|Psi^MTO_j> <Psi^MTO_j|, we have |F^MTO_k>=  P| F^MTO_k>. That is P is identical operator.
      ! Instead of <Psi^PMT_i|Psi^MTO_j>, we use Amat which is a modified version.
      if (gwinput_loaded) then
         ! mlo_* keys are all in eV (like mlo_emax). Convert once, here.
         ! mlo_w absent: method 4 takes 2.0 eV (optimized for it), the older
         ! methods take 0.2 Ry so the shipped samples reproduce exactly.
         if (tg_mlo_w == huge(0d0)) then
            eww = merge(2.0d0, 0.2d0*rydberg(), mlomethod == 4)/rydberg()
         else
            eww = tg_mlo_w/rydberg()
         endif
         dwin  = tg_mlo_delta /rydberg()
         wfrz  = tg_mlo_wfrz /rydberg()
         down  = tg_mlo_down /rydberg()
         ! Two hidden cuts were tried here on 2026-09-17 and left commented out
         ! (Samples/MLOsamples/BackUp_notes/mlo_nskip_cu_problem.md):
         !  - lower cut  theta_j *= sigma((EF+mlo_low - eps)/mlo_wlow)  -- Cu d model: 98 -> 97 meV, no help
         !  - character cut theta *= sigma((p_i - mlo_pcut)/mlo_pw), p_i = sum_j |<Psi_PMT_i|Psi_MTO_j>|^2 -- not tried
         !charcut = tg_mlo_pcut /= huge(0d0)
         !pcut = 0d0; pw = tg_mlo_pw
         !if (charcut) then
         !   pcut = tg_mlo_pcut
         !   if (iprx) write(stdo,ftox) ' Hreduction: mlo_pcut (hidden) character cut at p=',ftof(pcut),' width',ftof(pw)
         !endif
         !lowcut = tg_mlo_low /= huge(0d0)
         !elow   = 0d0; ewlow = eww
         !if (lowcut) then
         !   elow = eferm + tg_mlo_low/rydberg()
         !   if (tg_mlo_wlow /= huge(0d0)) ewlow = tg_mlo_wlow/rydberg()
         !   if (iprx) write(stdo,ftox) ' Hreduction: mlo_low (hidden) lower cut at EF',ftof(tg_mlo_low),&
         !        'eV, width',ftof(ewlow*rydberg()),'eV'
         !endif
      else
         call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
!         call getkeyvalue("GWinput","mlo_w",eww,default=0.2d0) !smoothing cutoff
      endif
      emax = evl(ndimMTO+nskip) - eferm   ! emax is the max of evl at ndimMTO+nskip. This is mainly useful for localized bands range.
!      emax = evlmto(ndimMTO) - eferm
      if (gwinput_loaded) then
         eee = tg_mlo_emax
         if (eee == huge(0d0)) eee = emax*rydberg()  ! sentinel ⇒ key absent
      else
         call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
!         call getkeyvalue("GWinput","mlo_emax",eee,default=emax*rydberg())  !eV relative to Ef.
      endif
      emax=eee/rydberg()+eferm

      !Amat is a modification of fac(ndimPMTx,ndimMTO), which is <Psi_PMT_i |Psi_MTO j>. Psi are eigenfunctions.
      allocate(Amat(ndimPMTx,ndimMTO),source=(0d0,0d0))!this is to avoid bug in ifort18.0.5
      !pwt = 1d0
      !if (charcut) then
      !   do i=1,ndimPMTx
      !      pwt(i) = fermidist((pcut - sum(abs(fac(i,1:ndimMTO))**2))/pw)   ! ~1 for model-like states
      !   enddo
      !endif
      mloloop : do j=1,ndimMTO
        ewuse = eww
        efrz  = -1d99 ! hard-freeze edge; active only for mlomethod=3
        ewfrz = 1d0
        if(mlomethod==0) then          ! Determine ecut for j to determine maxmum i index for PMT.
          ecut = max(emax, evlmto(j))  !  emax(relative to ef) is the rigid limit for localized MTOs
        elseif(mlomethod==1) then
          ecut = emax
        elseif(mlomethod==2) then
          ecut = evlmto(j)
        elseif(mlomethod==9) then
          ! Diagnostic (hidden): no energy cut at all, theta = 1 for every PMT state.
          ! Then Amat = fac and P = sum_i |Psi_PMT_i><Psi_PMT_i| is the identity on the
          ! PMT space, so |F_MLO> = |F_MTO> and the MLO bands are the MTO-only bands.
          ecut = 1d99
        elseif(mlomethod==4) then
          ! For target (1): a minimal model that reproduces the energy region
          ! around EF. method 0 with the hand-set emax replaced by an automatic
          ! floor, the same rule for every material:
          !   ecut_j = max( ecbot + mlo_delta , eps^MTO_j )
          ! mlo_emax is ignored; mlo_delta (Ry) is how far above the band edge the
          ! model is required to be accurate.
          !
          ! Two design points, both measured on the MLOsamples set:
          !  - ONE sigmoid, as in method 0, so eps^MTO_j enters the CUT POSITION.
          !    method 3 combines two sigmoids with max(), and since eps_frz <=
          !    ecut_j always holds, eps^MTO_j survives only in the tail and only
          !    while wfrz < eww. That is the weaker way to use it: at the same
          !    floor Fe gives 0.017 eV here vs 0.031 eV for method 3.
          !  - The floor is k-INDEPENDENT. Tracking the lowest unoccupied state
          !    at each k makes the window mean something different at every k and
          !    measures much worse on metals (Fe 0.075 eV).
          ! The floor is referenced to the global conduction edge ecbot, not to
          ! EF. For metals and narrow gaps the two coincide, but a wide gap needs
          ! the floor ABOVE the CBM or the CBM is left unconstrained: Al2O3:Cr
          ! (CBM at EF+6.2 eV) gives a gap error of -128 meV from EF+dwin and
          ! +7 meV from ecbot+dwin (m* ratio 0.58 -> 1.08). That is what the
          ! hand-set mlo_emax = 7 eV in the sample was doing.
          ecut = max(ecbot + dwin, evlmto(j))
        elseif(mlomethod==3) then
          ! Two-stage form: a freeze edge at the local conduction edge plus the
          ! method-0 style per-orbital cut, combined by max().
          !   theta_j = max( sigma((eps-efrz)/wfrz), sigma((eps-ecut_j)/eww) )
          !   efrz    = evl(nocc+1) + mlo_delta        (local band edge + dwin)
          !   ecut_j  = max( efrz, evl^MTO_j + mlo_down )
          ! Kept for comparison. Two measured drawbacks vs method 4:
          !  - efrz <= ecut_j always, so with wfrz = eww the first sigmoid wins
          !    everywhere and evl^MTO_j drops out entirely; it survives only in
          !    the tail, and only while wfrz < eww.
          !  - efrz follows the lowest unoccupied state AT EACH k, so the window
          !    means something different at every k (Fe: 0.075 eV vs 0.017 eV for
          !    a k-independent floor at the same height).
          block
            integer:: nocc
            nocc  = count(evl(nskip+1:ndimPMTx) < eferm)
            efrz  = evl(min(nskip+nocc+1,ndimPMTx)) + dwin
            ecut  = max(efrz, evlmto(j) + down)
            ewfrz = wfrz
            ewuse = eww
          endblock
        endif
        pmtloop: do i=nskip+1,ndimPMTx
          !thlow = 1d0
          !if (lowcut) thlow = fermidist((elow-evl(i))/ewlow)   ! rises through 1/2 at elow
          !Amat(i,j)= fac(i,j) * thlow * pwt(i) * max( fermidist((evl(i)-efrz)/ewfrz), fermidist((evl(i)-ecut)/ewuse) )
          Amat(i,j)= fac(i,j) * max( fermidist((evl(i)-efrz)/ewfrz), fermidist((evl(i)-ecut)/ewuse) )
          ! v6.3 two-stage theta-bar: near-unity inside the frozen window (narrow 0.05Ry edge -> band-edge
          ! curvature undistorted), the usual broad eww tail outside (rank + graded completeness).
        enddo pmtloop
      enddo mloloop
      Amat(1:nskip,:)=0d0
! do we need GramSchmidt orthogonalizaition? We expect lower is enphasized more for mode0 and for mode2.
      ! GramSchmidt of Amat (--gs): removed 2026-10-02 (an experimental option)
      ! if(c0_gs) call GramSchmidt(ndimPMTx,ndimMTO,Amat) !Amat= ¥bar{<Psi_PMT_i |Psi_MTO j>}
      
      ! cmlo(i,k) = \sum_j ¥bar{<Psi_PMT_i |Psi_MTO j>} <Psi_MTO j|F_MTO k> = <Psi_PMT_i|F_MLO k>
      ! |F^MLO_k> = P |F_MTO_k> = (\sum_{i,j} |Psi_PMT_i> ¥bar{<Psi_PMT_i |Psi_MTO j>} <Psi_MTO j|) |F_MTO k>=  |Psi_PMT> C * zMTO 
      ! Here P is a projector-like (probably not a projector because of GramSchmidt).
      nx = ndimPMTx
      cmlo_loc(ndimPMTx+1:ndimPMT,1:ndimMTO)=0d0
      cmlo_loc(1:ndimPMTx,1:ndimMTO) = matmul(Amat(1:ndimPMTx,1:ndimMTO),&
           matmul(transpose(dconjg(evecmto(:,:))),ovlmx(ix(1:ndimMTO),ix(1:ndimMTO)))) ! where <Psi_MTO j|MTO_k> = (evecmto*) @ ovlmx

      ! Per-orbital diagonal normalization at each k, |F^MLO_i> -> |F^MLO_i>/sqrt(<F^MLO_i|F^MLO_i>) (--mlo_diagnorm,
      ! --mlo_feb4): removed 2026-10-01 23:46 (user). The factor depends on k, so it changes the real-space orbital and the
      ! interpolated bands (C 0.74 eV, Al 0.10 eV off the mesh). The MLOs are normalized in real space instead, by one
      ! constant per orbital where the matrix elements of v and W are made (m_mlo_wfs, set_rnorm).
      ! MLODiagonalNormalize: if (c0_mlo_diagnorm .or. c0_mlo_feb4) then
      !    do i = 1, ndimMTO
      !       ddd = sum(dconjg(cmlo_loc(1:nx,i))*cmlo_loc(1:nx,i)) !<F^MLO|F^MLO>
      !       cmlo_loc(1:nx,i) = cmlo_loc(1:nx,i)/sqrt(ddd)
      !    enddo
      ! endif MLODiagonalNormalize
      ! Real-space normalization (2026-10-02, user): the MLO in real space, F_i0 = (1/N_k) sum_k F_i(k), has the square integral
      ! rnorm_i = (1/N_k) sum_k O_ii(k) = O_ii(R=0) of HamRsMLO (m_HamPMT, where the model is defined). Dividing by sqrt(rnorm_i),
      ! the same at every k, changes neither the shape of the orbital nor any band (H, O -> D H D, D O D), and makes (ii|ii) of v
      ! and W the U of a normalized orbital (the raw MLOs: rnorm 0.38 for Ni d, 0.18 for SrVO3 t2g).
      if(present(rnorm)) then
        forall(i=1:ndimMTO) cmlo_loc(1:nx,i) = cmlo_loc(1:nx,i)/sqrt(rnorm(i))
      endif
      ! Loewdin orthonormalization, the standard MLO (2026-10-02 user: "Loewdin orthogonalization is the standard of the MLO
      ! model; the functions may oscillate, they are localized, the bands do not change"). With the MLOs normalized above,
      !   |F~_j(k)> = sum_i |F_i(k)> X_ij(k),  X = O(k)^-1/2,  O = C^+ C                                       (Loewdin)
      ! the projected Wannier functions of the MLO subspace: orthonormal at every k (O~ = 1), the same subspace and bands, and
      ! they keep the symmetry and orbital labels (X commutes with the rotations within the subspace). Loewdin is not invariant
      ! under a rescaling of the input orbitals, so it is taken after the real-space normalization (rnorm): the model of
      ! m_HamPMT (LowdinModel) applies the same X = D (D O D)^-1/2 to its raw H and O, D = diag rnorm^-1/2.
      ! The old --mlo_ortho/--mlo_orthonorm (removed earlier on 2026-10-02) did the same without the normalization.
      if(present(lowdin)) then
        if(lowdin) then
          LowdinOrthonormal: block
            use m_lapack, only: zhev => zhev_h
            complex(8) :: oo(ndimMTO,ndimMTO), zo(ndimMTO,ndimMTO)
            real(8) :: eo(ndimMTO)
            integer :: istat
            oo = matmul(dconjg(transpose(cmlo_loc(1:nx,:))), cmlo_loc(1:nx,:))
            istat = zhev(oo, n=ndimMTO, evl=eo)
            if(minval(eo) <= 0d0) call rx('Hreduction: the MLO overlap is not positive definite (Loewdin)')
            forall(i=1:ndimMTO) zo(:,i) = oo(:,i)/sqrt(eo(i))
            cmlo_loc(1:nx,:) = matmul(cmlo_loc(1:nx,:), matmul(zo, dconjg(transpose(oo))))
          endblock LowdinOrthonormal
        endif
      endif


      ! |F^MLO j'>= |F^PMT_i'> z^PMT_i'i cmlo(i,j)
      do i=1,ndimMTO
        do j=1,ndimMTO
          hammout(i,j)= sum( dconjg(cmlo_loc(1:nx,i))*evl(1:nx)*cmlo_loc(1:nx,j)) !|F^MLO_i> = |Psi^PMT_j> cmlo(j,i)
          ovlmout(i,j)= sum( dconjg(cmlo_loc(1:nx,i))*cmlo_loc(1:nx,j) ) !<F^MLO|F^MLO>
        enddo
      enddo
      if(present(cmlo)) cmlo = cmlo_loc
      if(present(zMLO)) zMLO = matmul(evecpmt(:,1:nx), cmlo_loc(1:nx,:)) !PMT-basis-function coefficients

!       GetZCZ: if(present(zcplz)) then ! This is for k mesh of GWinput. zcplz(ndimPMT, ndimMLO), ndimMLO=ndimMTO
!         diagonalizeMTO: block         !  |F^MLO j> = |Psi^PMT i> zcplz(i,j)  !MLO eigenfunctions constructed from the PMT basis.
!           real(8):: evlx(ndimMTO),oveps=0d0
!           complex(8):: evecmlo(ndimMTO,ndimMTO),hh(ndimMTO,ndimMTO),oo(ndimMTO,ndimMTO)
!           hh=hammout
!           oo=ovlmout
!           call zhev_tk4(ndimMTO,hh,oo,ndimMTO,nev, evlx,evecmlo, oveps)
!           ! nmx=ndimMTO.
!           ! CAUTION for zhev_tk4: If nmx=0, only eigenvalues returned. Diangonalize (hamm- evl ovlm) z=0
!           ! evecmlo= z^MLO_jj' convert to eigenfunction Psi
! !          zcplz = matmul(matmul(evecpmt(1:ndimPMT,1:nx), cmlo(1:nx,1:ndimMTO)), evecmlo) !zPMT* C * zMTO
!           zcplz = cmlo(1:nx,1:ndimMTO), evecmlo) !  |F^MLO j> = |Psi^PMT i> cpmo(i,j)
!                                             !MLO eigenfunctions constructed from the PMT basis.
!           !   do i=1,ndimMTO
!           !     write(stdo,ftox)'eigen111',i,ftof(qp,3),'  ',ftof(evl(i+nskip)),' ',ftof(evlx(i))
!           !   enddo
!         endblock diagonalizeMTO
!       endif GetZCZ
       
     endblock ModifyMatrixElements
     return
   end subroutine Hreduction
   real(8) function fermidist(x)
     real(8),intent(in) :: x
     if(x>100d0) then
       fermidist=0d0
     elseif(x<-100d0) then
       fermidist=1d0
     else
       fermidist=1d0/(exp(x)+1)
     endif
   end function fermidist
   subroutine GramSchmidt(nv,n,zmel)
     integer:: igb=1,it,itt,n,nv
     complex(8):: ov(n),vec(nv),dnorm2(nv),zmel(nv,n)
     real(8):: dnorm
     do it = 1,n
       vec(:)= zmel(:,it)
       do itt = 1,it-1
         ov(itt) = sum( dconjg(zmel(:,itt))*vec(:))
       enddo
       vec = vec - matmul(zmel(:,1:it-1),ov(1:it-1))
       zmel(:,it) = vec/sum(dconjg(vec)*vec)**.5d0
     enddo
   end subroutine GramSchmidt
end module m_hreduction
