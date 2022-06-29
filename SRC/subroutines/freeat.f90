module m_freeat
  use m_lgunit,only: stdo,stdl
contains
  subroutine freeat()
    use m_ext,only:sname
    use m_lmfinit,only: smalit,ctrl_lxcf,ham_seref,nsp,nspec, sspec=>v_sspec,&
         idmod,slabl,vmtz,eref,rs3,eh3,nmcore,coreh,coreq,rcfa,pnux=>pnu,pzx=>pz,qnu
    use m_ftox
    !- For each species, makes free atom self-consistent
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   sspec :struct for species-specific information; see routine uspec
    !i     Elts read: rsmfa rfoca coreq z rmt a nr p q idmod lmxa pz eref
    !i                rs3 eh3 vmtz
    !i     Stored:    name coreh z rmt a nr norp ntorb
    !o Outputs via iofa
    !l Local variables
    !l   ccof  :coefficient to fit of core tail to smoothed Hankel
    !l   ceh   :energy of core tail to smoothed Hankel
    !l   sumtc :core kinetic energy
    !r Remarks
    !u Updates
    !u   01 Feb 06 Enables renormalized free atom density
    !u   01 Jul 05 Skips spheres with Z=0 and R=0
    !u   21 Jun 04 Added fit of sm. Hankel tails to local orbitals
    !u   18 Sep 03 (ATP) Enabled partial core occupation
    !u   06 Sep 03 Constrain rsm in fit to FA wave function
    !u   18 Mar 03 Altered sign of magnetic moment to conform to std
    !u   19 Apr 02 Redesigned freats call to avoid the use of structures
    !u   22 Dec 01 Adjustments to accomodate changes in phidxasz
    !u   22 Mar 01 Added printout of reference energy
    !u   10 Jun 00 spin polarized
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ifi,iprint,is,nglob,nr,nrmt,nrmx, &
         n0,nkap0,nxi,nxi0,nrmix,lxcfun,igets,lmxa,kcor,lcor
    character(8) :: spid,chole(8)
    parameter ( nrmx=1501, nxi0=10, n0=10, nkap0=3)
    double precision :: qc,ccof,ceh,z,rmt,rfoca,qcor(2),a,sumec, &
         sumtc,seref,dgets,dgetss,etot
    double precision :: hfc(nxi0,2),exi(nxi0),hfct(nxi0,2)
    double precision :: v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),rofi(nrmx*2)
    double precision :: pnu(n0,2),pz(n0,2),qat(n0,2)!,rcfa(2)
    double precision :: rtab(n0,2),etab(n0,2),rsmfa
!    double precision :: rs3,eh3!,vmtz
    !integer :: idmod(n0)
    integer:: iofa, i_dum,ifile_handle
    integer:: ifives,ifiwv
    character strn*120
    logical :: cmdopt
    ifi = ifile_handle()
    open(ifi,file='atm.'//trim(sname))
    rewind ifi
    hfct = 0d0
    open(newunit=ifives,file='vesintatm.'//trim(sname)//'.chk')
    ifiwv = ifile_handle()
    open(ifiwv,file='veswavatm.'//trim(sname)//'.chk')
    do  is = 1, nspec ! Takao found that Li requires -0.5 to have positive smooth rho.
       if(sspec(is)%z<3.5) then !At least for Li, fitting is not good (negative smooth rho).
          exi(1) = -0.5
          exi(2) = -1
          exi(3) = -2
          exi(4) = -4
          exi(5) = -6
          exi(6) = -9
          exi(7) = -15
          nxi = 7
       else !this is original setting. For CrN. This is necessary(maybe long range cutoff is required).
          exi(1) = -1
          exi(2) = -2
          exi(3) = -4
          exi(4) = -6
          exi(5) = -9
          exi(6) = -15
          nxi = 6
       endif
       nrmix=smalit
       lxcfun = int(ctrl_lxcf)
       spid = slabl(is) !sspec(is)%name
       !       rsmfa= sspec(is)%rsmfa
       rfoca= sspec(is)%rfoca
       qcor = coreq(:,is)
       chole= coreh(is)
       call gtpcor(sspec,is,kcor,lcor,qcor)
       z   = sspec(is)%z
       rmt = sspec(is)%rmt
       rsmfa=.5d0*rmt            ! moved to here 2022-6-27
       a   = sspec(is)%a
       nrmt= sspec(is)%nr
       if (z == 0 .AND. rmt == 0) cycle !floating orbital
       pnu(:,1)=  pnux(1:n0,1,is) ! sspec(is)%p
       if(nsp==2) pnu(:,2)= pnu(:,1)
       qat(1:n0,1:nsp)=  qnu(1:n0,1:nsp,is) !sspec(is)%q
!       idmod= sspec(is)%idmod
       lmxa = sspec(is)%lmxa
       pz(:,1) =  pzx(1:n0,1,is) ! sspec(is)%pz
       if(nsp==2) pz(:,2)= pz(:,1)
!       eref = eref(is)
       !rs3=sspec(is)%rs3
       !eh3=sspec(is)%eh3
       !vmtz=sspec(is)%vmtz
       !rcfa=sspec(is)%rcfa
       print *,'goto freats'
       call freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,kcor,lcor,qcor, &
            nrmix,1,lxcfun,z,rmt,a,nrmt,pnu,pz,qat,rs3(is),eh3(is),vmtz(is),rcfa(:,is), &
            idmod(:,is),lmxa,eref(is),rtab,etab,hfc,hfct,nr,rofi,rho,rhoc,qc,ccof, &
            ceh,sumec,sumtc,v,etot,nmcore(is),ifives,ifiwv)
       print *,'end of freats: spid nmcore=',spid,nmcore(is)
       if (iprint()>40) write(stdo,"(a)")' write free atom data for species  '//trim(spid)
       if (nsp == 2 .AND. nr > nrmt) then
          call dcopy(nrmt,rho(1+nr),1,rho(1+nrmt),1)
          call dcopy(nrmt,rhoc(1+nr),1,rhoc(1+nrmt),1)
          call dcopy(nrmt,v(1+nr),1,v(1+nrmt),1)
       endif
       i_dum = iofa(spid,nxi0,nxi,exi,hfc,hfct,rsmfa,z,rmt, &
            a,nrmt,qc,ccof,ceh,sumtc,rho,rhoc,v,-ifi)
    enddo
    close(ifi)
    if(iprint() > 30)  then
       seref = ham_seref
       write(stdo,"(a,f35.12)") 'Sum of reference energies: ',seref
    endif
    close(ifives)
    close(ifiwv)
  end subroutine freeat

  subroutine freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,kcor,lcor,qcor, &
       nrmix,lwf,lxcfun,z,rmt,a,nrmt,pnu,pz,qat,rs3,eh3,vmtz,rcfa, &
       idmod,lmxa,eref,rtab,etab,hfc,hfct,nr,rofi,rho,rhoc,qc,ccof,ceh, &
       sec,stc,v,etot,nmcore,ifives,ifiwv)
    use m_lmfinit,only: nsp,lrel
    use m_ftox

    use m_ext,only:sname
    use m_getqvc
    !- Makes one free atom self-consistent, fits rho tails to smoothed Hankels
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   spid  :species label (for printout)
    !i   is    :species index
    !i   nxi0  :nxi0: leading dimension of hfc,hfct
    !i   nxi   :number of hankel functions used in fitting of tails
    !i   exi   :hankel energies used in fitting of tails
    !i   rfoca :smoothing radius for hankel fit to core
    !i   rsmfa :smoothing radius for hankel fit to valence
    !i   kcor  :(partial core occupation) p.q.n for occupation
    !i   lcor  :(partial core occupation) l quantum for occupation
    !i   qcor  :(partial core occupation) core charge and moment
    !i   nrmix :nrmix(1) = maximum number of interations in sphere
    !i         :           before giving up on self-consistency
    !i         :xxx nrmix(2) = no prior iterations Anderson mixing in
    !i         :           self-consistency cycle.
    !i   lwf   :1 print information about wave functions
    !i   lxcfun:selects local exchange-correlation functional
    !i   z     :nuclear charge
    !i   rmt   :augmentation radius, in a.u.
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nrmt  :number of mesh points from origin to rmt
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   pz    :boundary conditions for local orbitals
    !i   qat   :valence charges for each l channel
    !i   rs3   :minimum allowed smoothing radius in attaching Hankel tails
    !i         :to local orbitals
    !i   eh3   :Hankel energy when attaching Hankel tails to high-lying
    !i         :local orbitals
    !i   vmtz  :parameter used in attaching Hankel tails to local orbitals
    !i         :It is used as a constant shift to Hankel energies for the
    !i         :fitting of local orbitals to Hankel tails. Thus vmtz
    !i         :is an estimate for the potential at the MT radius.
    !i   idmod :0,1 or 2, specifing how the enu is set for an l-channel
    !i   lmxa  :augmentation l-cutoff
    !i   eref  :reference energy (used for printout)
    !o Outputs
    !o   rtab  :smoothing radius for optimized wave function
    !o   etab  :energy for optimized wave function
    !o   hfc   :fit coeffs for valence density,
    !o   hfct  :contains fit coeffs for full density (not calc. now)
    !o   nr    :number of radial mesh points for spherical rho
    !o   rofi  :rofi(1..nr)=radial mesh for points
    !o         :rofi(nr+1..2*nr)=radial mesh weights
    !o   rho   :free-atom valence density
    !o   rhoc  :free-atom core density
    !o   qc    :Sphere core charge
    !o   ccof  :coefficient to fit of core tail to smoothed Hankel
    !o   ceh   :energy of core tail to smoothed Hankel
    !o   sec   :sum of core eigenvalues
    !o   stc   :core kinetic energy
    !o   v     :spherical potential
    !l Local variables
    !l   itab  :itab(l+1)=1  a wave function was optimzed for this l
    !l   pnul  :EITHER : pnu for valence state, OR
    !l         :local orbital if DEEPER than valence (pz<pnu)
    !r Remarks
    !u Updates
    !u   01 Feb 06 Enables renormalized free atom density
    !u   19 Apr 02 Redesigned input to avoid the use of structures
    !u   10 Apr 02 Redimensionsed etab,rtab to accomodate larger lmax
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: nrmx,nrmt,is,nxi0,nxi,nrmix,lwf,lxcfun,n0,kcor,lcor
    parameter (nrmx=1501,n0=10)
    character(8) :: spid
    double precision :: rsmfa,rfoca,qc,ccof,ceh,sec,stc,z,rmt,a,eref, &
         v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),hfc(nxi0,*),hfct(nxi0,*), &
         exi(*),rtab(n0,2),etab(n0,2),rofi(nrmx*2),rs3,eh3,vmtz,qcor(2)
    ! ... Local parameters
    logical :: cmdopt
    integer :: ncmx,nvmx
    parameter (ncmx=50, nvmx=20)
    integer :: idmod(n0)
    !     integer idmoz(n0)
    character str*8,strn*32
    double precision :: rmax,b,etot,dq, &
         ec(ncmx),ev(nvmx),sumev,vrmax(2),exrmax(2),ekin,utot,rhoeps, &
         amgm,rhrmx,qvt,qtot,qct,qvin,qcin,r,wt0,wt1,qtt,qtin,pnul, &
         pzl,pnu(n0,2),qat(n0,2),pl(n0,2),rhoin(nrmx*2), &
         rhot(nrmx*2),ql(3,n0,2),pz(n0,2),rcfa(2) ,qatbk(n0,2)
    !     double precision qz(n0,2)
    integer :: i,ifi,ipr,iprint,isp,isw,l,lfrz, &
         lgrad,lmxa,lplfa,nitmax,nmix,nr,lplawv,irchan(n0)
    real(8) ,allocatable :: g_rv(:)
    real(8) ,allocatable :: psi_rv(:)
    integer :: itab(n0,2)
    integer ::iwdummy
    integer:: ipl, ipz,iplx,iplv,ix, nmcore, ifives,ifiwv,iplz,miplzipl
    real(8):: qcc,qzz, qvalv,qvaltot,qval,vsum,qelectron,qvalz
    real(8),allocatable:: plplus(:,:),qlplus(:,:)
    ipr    = iprint()
    lfrz   = 0
    lgrad  = lxcfun/10
    if(ipr>29) then
       print *
       do i = 1, nsp
          do l = 0, lmxa
             write(stdo,"(a,2i3,2f10.3)")'ttt: pnu qat=',i,l,pnu(l+1,i),qat(l+1,i)
          enddo
       enddo
    endif

    ! --- Get species data ----
    call dpzero(pl,2*n0)
    call dpzero(ql,2*3*n0)
    allocate(plplus(0:lmxa,nsp),qlplus(0:lmxa,nsp))
    !! plplus, qlplus are higher ones (when PZ exist).
    plplus=0d0
    qlplus=0d0
    !      print *,'NOTE: when we have two valence: P and PZ, We assume eigen(PZ) is deeper than eigen(P).'
    do  101  i = 1, nsp
       do  10  l = 0, lmxa
          pnul = pnu(l+1,i)
          pzl  = mod(pz(l+1,1),10d0)
          if(nsp==1) qvaltot = qat(l+1,1)
          if(nsp==2) qvaltot = qat(l+1,1)/2d0 - qat(l+1,2)/2d0*dble(2*i-3)
          qval=qvaltot
          !!   Make pz the valence state:
          !! NOTE: Because of historical reason,  PZ channel is treated as P; P channel is treated by plplus and qlplus.
          !! plplus and qlplus are for P channel when deeper PZ exits.
          !!       jun2012 feb2011
          if (pzl /= 0) then
             if(int(pnul-1) == int(pzl)) then !pz is deeper than pnul
                plplus(l,i)= int(pnul) + .5d0
                if(nsp==2 .AND. pnul==0) plplus(l+1,i) = plplus(l+1,1)
                ! plplus and qlplus are for P channel
                if(nsp==1) qval = min(qvaltot, 2d0*(2d0*l+1d0)) ! We assume eigen(PZ) is deeper then eigen(P).
                if(nsp==2) qval = min(qvaltot, 2d0*l+1d0)
                qvalv=0d0
                if(qvaltot-qval>0d0) qvalv = qvaltot-qval
                qlplus(l,i)= qvalv ! For P channel. Preserve pl and ql to plplus,and qlplus.
                pnul = pzl
             elseif (int(pnul+1) == int(pzl)) then !pz is higher than pnul. plplus and qlplus are for
                plplus(l,i) = pzl
                if(nsp==2 .AND. pnul==0) plplus(l+1,i) = plplus(l+1,1)
                if(nsp==1) qval = min(qvaltot, 2d0*(2d0*l+1d0))
                if(nsp==2) qval = min(qvaltot, 2d0*l+1d0)
                qvalz = 0d0
                if(qvaltot-qval>0d0) qvalz = qvaltot-qval
                qlplus(l,i)= qvalz
             else
                call fexit3(-1,111,' Exit -1 freeat, l=%i:  '// &
                     'sc PZ=%d incompatible with valence (you may need to add P simultaneously.) P=%;3d',l,pzl,pnul)
             endif
          endif
          !! pl,ql are for PZ when PZ= nonzero exist... somehow confusing.
          pl(l+1,i) = int(pnul) + .5d0
          ql(1,l+1,i) = qval
          if (nsp == 2 .AND. pnul == 0) pl(l+1,i) = pl(l+1,1)
          ql(2,l+1,i) = 0d0
          ql(3,l+1,i) = 0d0

          if(iprint()>60) then
             if (pzl/=0 .AND. int(pnul-1)==int(pzl)) then !pz is deeper than pnul
                if(nsp==1)write(stdo,"(a,i2,3f8.3)")'=== Charge for l: Qtot Qpz Qv=',l,qvaltot,ql(1,l+1,i),qlplus(l,i)
                if(nsp==2)write(stdo,"(a,i2,3f8.3)")'=== Charge for l: Qtot Qpz Qv=',l,qvaltot,ql(1,l+1,i),qlplus(l,i)
             elseif (pzl/=0) then !pz is higher than pnul
                if(nsp==1)write(stdo,"(a,i2,3f8.3)")'=== Charge for l: Qtot Qpz Qv=',l,qvaltot,qlplus(l,i),ql(1,l+1,i)
                if(nsp==2)write(stdo,"(a,i2,3f8.3)")'=== Charge for l: Qtot Qpz Qv=',l,qvaltot,qlplus(l,i),ql(1,l+1,i)
             else
                if(nsp==1)write(stdo,"(a,i2,3f8.3)")'=== Charge for l: Qtot=Qv=    ',l,qvaltot,ql(1,l+1,i)
                if(nsp==2)write(stdo,"(a,i2,3f8.3)")'=== Charge for l: Qtot=Qv=    ',l,qvaltot,ql(1,l+1,i)
             endif
          endif
10     enddo
101 enddo
    !!
    print *
    write(stdo,"(a)")'conf:------------------------------------------------------'
    write(stdo,"(a,a)")'conf:SPEC_ATOM= '//trim(spid),' : --- Table for atomic configuration ---'
    !      write(stdo,"(a)")'conf When int(P)z .ne. int(P), Qval: Q for MTOcore(PZ)+MTO(P)'
    write(stdo,"(a)")'conf:  isp  l  int(P) int(P)z    Qval     Qcore   CoreConf'
    do i = 1, nsp
       do l = 0, lmxa
          pzl  = mod(pz(l+1,1),10d0)
          iplz=999
          if (pzl/=0 .AND. int(pnu(l+1,i)-1)==int(pzl)) then !PZ exists and is deeper than P.
             ipl   = int(plplus(l,i)) !P
             iplz  = int(pl(l+1,i)) !PZ  deeper
          else                   !PZ exist and higher than P. Or No PZ exist.
             ipl  = int(pl(l+1,i)) !P  deeper
             iplz = int(plplus(l,i)) !PZ
          endif
          miplzipl = ipl
          if(iplz/=0) miplzipl = min(iplz,ipl)
          qcc= (miplzipl-l-1)*(4*l+2)/nsp
          !!  plplus, qlplus for P;   ql for Pz
          write(stdo,"('conf: ', i4,i3,3x,i5,i3,6x,f8.3,1x,f8.3,' => ',10(i1,','))") &
               i,l,ipl,iplz,qlplus(l,i)+ql(1,l+1,i),qcc, (ix,ix=l+1,miplzipl-1)
       enddo
       !      write(stdo,"(a)")'conf:-----------------------------------------------------'
    Enddo
    write(stdo,"('usedQ=',10f10.3)") (sum(qlplus(l,1:nsp)+ql(1,l+1,1:nsp)),l=0,lmxa)

    ! ov 2010 QvalCheck
    do i = 1, nsp
       do l = 0, lmxa
          if( .NOT. cmdopt('--skip_qvalcheck',16,0,strn)) then
             if(ql(1,l+1,i)<-1d-10) then
                call rx('conf:negative Qval. Check SPEC_ATOM_Q & MMOM or Use --skip_qvalcheck')
             endif
          endif
       enddo
    enddo

    call getqvc(nsp, n0, lmxa, z, pl, ql, 0, 0, kcor, lcor, qcor, qc, qtot, amgm )
    qtot=qtot+sum(qlplus(:,:)) !takao feb2011
    amgm=amgm+sum(qlplus(:,1)-qlplus(:,nsp)) !takao feb2011
    if(ipr>=20) write(stdo,ftox)'conf: Species ',spid,'Z=',ftof(z,2), &
         'Qc=',ftof(qc,3),'R=',ftof(rmt),'Q=',ftof(qtot),'nsp=',nsp,'mom=',ftof(amgm)

    ! --- Set up radial mesh for free atom ---
    rmax = 50d0
    if (z < 10) rmax = 25
    if (z <=  6) rmax = 20
    call pshpr(0)
    call rmesh(z,rmt,lrel,lgrad,nrmx,a,nrmt)
    call poppr
    b = rmt/(dexp(a*nrmt-a)-1d0)
    nr = 1d0+dlog(1d0+rmax/b)/a
    if (mod(nr,2) == 0) nr = nr-1
    rmax = b*(dexp(a*(nr-1))-1d0)
    call info5(21,0,0,'conf:   rmt=%,6;6d  rmax=%,6;6d'// &
         '  a=%d  nr=%i  nr(rmax)=%i',rmt,rmax,a,nrmt,nr)

    ! --- Make atom self-consistent ---
    nitmax = nrmix
    !      nmix = nrmix(2)
    nmix = -30
    !     call pshpr(min(iprint(),40))
    ec(1) = 0
    print *,'goto atomc xxx'
    qelectron = z+qtot     !total number of electrons. 22mar2013
    call atomsc(.false.,n0,nsp,lmxa,z,0d0,kcor,lcor,qcor,rmax,a,nr, &
         rofi,ec,ev,pl,ql,idmod,v,0d0,rhoin,rho,rhoc,nmix,qc,sec,stc, &
         sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,dq,exrmax,'gue', &
         nitmax,lfrz,plplus,qlplus,nmcore,qelectron,vsum)
    print *,'end of atomsc xxxxx'
    print *,'vsum=',vsum,is

    if(ifives>0) then
       print *,'ifives=',ifives
       write(ifives,"(f23.15,i5,a)") vsum,is, ' !spatical integral of electrostatic potential'
       print *,' write end of ifives'
    endif
    if(ifiwv>0) write(ifiwv,"(f23.15,' ! total charge')") sum(qlplus(0:lmxa,1:nsp)+ql(1,1:lmxa+1,1:nsp))

    dq=dq+sum(qlplus(:,:)) !takao feb2011
    if (ipr>=20)write(stdo,ftox)'sumev=',ftof(sumev),'etot=',ftof(etot), &
         'eref=',ftof(eref),'etot-eref=',ftof(etot-eref)

    !! june2012takao comment out followings. because confusing.
    !      call dcopy(lmxa+1,ql,3,qat,1)
    !      call dcopy(lmxa+1,ql(1,1,nsp),3,qat(1,2),1)
    !      call awrit4('fa  Pl %n:-1d  Ql %n:-1d',
    !     .' ',80,stdl,lmxa+1,pl,lmxa+1,qat)
    !      if (nsp .eq. 2)
    !     .call awrit4('fa  Pl2%n:-1d  Ql2%n:-1d',
    !     .' ',80,stdl,lmxa+1,pl(1,nsp),lmxa+1,qat(1,nsp))

    if (dabs(dq) > 1d-5 .AND. iprint() >= 10) &
         write(stdo,ftox)' freeat (warning) atom not neutral, Q=',ftof(dq)
    ! .. Subtract core from density to make valence density
    do  isp = 1, nsp
       do  i = 1, nr
          rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)
          rho(i+(isp-1)*nr) = rho(i+(isp-1)*nr)-rhoc(i+(isp-1)*nr)
       enddo
    enddo

    ! --- Renormalize atom density or potential ---
    call ivset(irchan,1,n0,0)
    !! takao jun2012: rnatm is not tested ---> this is related to SPEC_ATOM_RCFA on.
    call rnatm(pl,qat,n0,irchan,lmxa,z,a,b,rofi,ev,nr,rcfa,nsp,v,rho,plplus,qlplus)
    !     call prrmsh('starting total rho',rofi,rhot,nr,nr,nsp)
    do  isp = 1, nsp
       do  i = 1, nr
          rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)+rhoc(i+(isp-1)*nr)
       enddo
    enddo
    !     call prrmsh('ending total rho',rofi,rhot,nr,nr,nsp)

    ! --- Print info about free-atom wavefunctions ---
    if (lwf /= 0) then
       if (ipr > 30) then
          allocate(g_rv(nr*2))
          allocate(psi_rv(nr*(lmxa+1)*nsp))
          lplawv = 0
          if (ipr >= 50) lplawv = 1
          call pratfs ( spid , lplawv , z , a , nr , rmax , nrmt , lmxa &
               , pl , nsp , v , rofi , g_rv , psi_rv, plplus,qlplus,ifiwv)
          deallocate(psi_rv)
          deallocate(g_rv)
       endif

       ! --- Optimise smooth-Hankel basis ---
       call dvset(rtab,1,n0,-1d0)
       call dpzero(etab,n0)
       i = 1
       if ( .NOT. cmdopt('--noopt',7,0,strn)) then
          if (cmdopt('--norscnst',10,0,strn)) i = 0
          if (z > 0) then
             !c     stop 'xxxxxxxxxxxxxxx end cccccxxx000'
             !             write(stdo,*) '111 ',i,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,nsp
             !             write(stdo,*) '222 ',sum(v(1:nr)),sum(rofi(1:nr)),spid
             ! cccccccccccccccccccccccccc
             call optfab(i,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,nsp,v, &
                  rofi,spid,itab,rtab,etab)
             !   ... Fit value and slope of local orbitals
             call ftfalo(i,z,a,nr,rmax,nrmt,rmt,lmxa,pnu,pz,rs3,eh3,vmtz, &
                  nsp,v,rofi,spid)
          endif
       endif
    endif

    ! --- Print charges within/outside MT sphere ---
    if (z > 0) then
       qvt = 0d0
       qct = 0d0
       qvin = 0d0
       qcin = 0d0
       do  301  isp = 1, nsp
          do  30  i = 1, nr
             r = rofi(i)
             wt0 = rofi(i+nr)
             !       wt0 = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
             wt1 = wt0
             if (i > nrmt) wt1 = 0d0
             if (i == 1 .OR. i == nrmt) wt1 = a*(r+b)/3d0
             if (i == 1 .OR. i == nr)   wt0 = a*(r+b)/3d0
             qvin = qvin + wt1*rho(i+(isp-1)*nr)
             qcin = qcin + wt1*rhoc(i+(isp-1)*nr)
             qvt = qvt + wt0*rho(i+(isp-1)*nr)
             qct = qct + wt0*rhoc(i+(isp-1)*nr)
30        enddo
301    enddo
       qtt = qvt + qct
       qtin = qvin + qcin
       if (ipr >= 40) write (stdo,550) qvin,qvt-qvin,qvt,qcin,qct-qcin, &
            qct,qtin,qtt-qtin,qtt
550    format(/' Charges:     inside',7x,'outside',7x,'sum' &
            /' valence',3f13.6/' core   ',3f13.6/' total  ',3f13.6)

       write (stdl,710) z,rmax,qc,qct-qcin,dq,etot
710    format('fa Z',f6.1,'   rm',f7.2,'  qc',f6.2,'  qspl',f8.5, &
            '  dq',f8.5,'  Etot',f15.6)
    else
       qvt = 0
    endif

    ! --- Attach smooth Hankel tails to valence density ---
    !     lplfa = nglob('lplfa')
    lplfa = 0
    if (qvt > 1d-6) then
       if (lplfa == 1) then
          write(stdo,344)
344       format(/' write plot file with valence density..')
          if (is < 10) write (str,'(''pl'',i1)') is
          if (is >= 10) write (str,'(''pl'',i2)') is
          !          ifi = fopna(str,-1,0)
          open(newunit=ifi,file=trim(str)//'.'//trim(sname))
          write (ifi,490) spid,rmt,rsmfa,nxi
490       format('# fit to fa density: ',a/ &
               '# rmt=',f7.3,'   rsm=',f7.3,'   nxi=',i2)
          close(ifi)
       endif
       call tailsm(0,nr,nrmt,nsp,a,b,rmt,rsmfa,nxi0,nxi,exi,rofi, &
            rho,rhot,hfc,hfct)
       !       call prrmsh('rho-fa',rofi,rho,nr,nr,1)
    else
       call dpzero(hfc, nxi0*nsp)
       call dpzero(hfct, nxi0*nsp)
    endif
    ! --- Fit analytical expression to tail of core density ---
    call fctail(nr,nrmt,a,b,rfoca,rofi,rhoc,ccof,ceh)
  end subroutine freats


  subroutine pratfs(spid,lplawv,z,a,nr,rmax,nrmt,lmaxa,pl,nsp,v, &
       rofi,g,psi,plplus,qlplus,ifiwv)
    use m_ext,only:sname
    !- Prints out core and valence energy levels of free-atom
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   spid  :species label
    !i   z     :nuclear charge
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rmax  :free-atom radius, in a.u.
    !i   nrmt  :mesh for MT radius
    !i   lmaxa :muffin-tin cutoff
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
    !i         :pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   v     :spherical potential (atomsr.f)
    !i   rofi  :radial mesh points
    !i   g     :normalized wave function times r (work array)
    !o Outputs
    !o   psi   :normalized wave functions for each l
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: n0,nr,lmaxa,nrmt,lplawv,nsp
    parameter (n0=10)
    double precision :: z,a,rmax,pl(n0,nsp),v(nr,nsp),rofi(nr), &
         g(2*nr),psi(nr,0:lmaxa,nsp)
    ! ... Local parameters
    integer :: nglob,isp,l,konfig,nn,nre,i,konf, &
         konfg(0:8),ifi,lmaxc
    double precision :: ev(0:20),pi,b,tol,eb1,eb2,dl,val,slo,sum,pzero, &
         pmax,ctp,ecor,rmt
    character(1) :: lsym(0:n0-1), cc
    character:: spid*8
    character(15)::  str
    data lsym /'s','p','d','f','g','5','6','7','8','9'/
    integer::iz,ifiwv
    !! feb2011 "plplus,qlplus" mechanism is a fixing. Previous version did not allow SPEC_ATOM_Q setting for l with PZ.
    !! Now, the "plplus,qlplus" mechanism allows to set Q (valence charge, not including semicore charge).
    !! Our current version assumes MTOcore(specified by PZ) is below MTO(specified by P).
    real(8):: plplus(0:lmaxa,nsp),qlplus(0:lmaxa,nsp)
    real(8):: sumr !mar2013
    pi   = 4d0*datan(1d0)
    if (lmaxa > n0-1) call rx('pratfs:  lmax too large')
    b = rmax/(dexp(a*nr-a)-1d0)
    rmt = b*(dexp(a*nrmt-a)-1d0)
    tol = 1d-8
    write(stdo,580)
580 format(/' Free-atom wavefunctions:')

    !      allocate(qinrmt(0:lmxa),

    if(ifiwv>0) write(ifiwv,"(i2,i3,' !nsp,lmaxa')")nsp,lmaxa
    do  80  isp = 1, nsp
       ! --- Valence states ---
       if (isp == 1) write(stdo,401)
       if (isp == 2) write(stdo,'(/'' spin 2:'')')
       eb1 = -50d0
       eb2 =  50d0
       do  201  l = 0, lmaxa
          do 20 iz= 0,1
             if(iz==1 .AND. plplus(l,isp)<= 0d0) cycle !feb2011
             konfig = pl(l+1,isp)
             if(iz==1) konfig=plplus(l,isp)
             dl = dtan(pi*(0.5d0-pl(l+1,isp)))
             nn = konfig-l-1
             ev(l) = -0.5d0
             val = rmax
             slo = dl+1
             if (rmax > 9.99d0) then
                val = 1d-30
                slo = -val
             endif
             call rseq(eb1,eb2,ev(l),tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b, &
                  rofi,nr,nre)
             call gintsl(g,g,a,b,nr,rofi,sum)
             call gintsl(g,g,a,b,nrmt,rofi,pmax)
             sum = sum - pmax
             sumr=pmax

             call ppratf(ev(l),z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
             cc = ' '
             if (dabs(ctp-rmax) < 1d-3) cc = '*'
             write(stdo,400) konfig,lsym(l),ev(l),pzero,pmax,ctp,cc,sum
400          format(i4,a1,f14.5,2x,3f12.3,a,f12.6)
401          format(' valence:',6x,'eval',7x,'node at',6x,'max at',7x, &
                  'c.t.p.   rho(r>rmt)')
             if(ifiwv>0) write(ifiwv,"(i2,i3,i3,d23.15,d23.15,' !isp,l,eval,last is norm within MT')") isp,l,konfig,ev(l),sumr
             !   ... Copy valence wavefunction to psi
             do  24  i = 1, nr
                psi(i,l,isp) = g(i)
24           enddo
20        enddo
201    enddo

       ! --- Core states ---
       write(stdo,403)
       eb1 = -2.5d0*z*z-5d0
       eb2 = 50d0
       call config(pl,lmaxa,z,konfg,lmaxc)

       do  4011  konf = 1, 8
          do  40  l = 0, min(konf-1,lmaxc)
             konfig = konfg(l)
             if (konf >= konfig) goto 40
             nn = konf-l-1
             ecor = -50d0
             val = 1d-30
             slo = -val
             call rseq(eb1,eb2,ecor,tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,rofi, &
                  nr,nre)
             call gintsl(g,g,a,b,nr,rofi,sum)
             call gintsl(g,g,a,b,nrmt,rofi,pmax)
             sum = sum - pmax
             call ppratf(ecor,z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
             write(stdo,400) konf,lsym(l),ecor,pzero,pmax,ctp,' ',sum
403          format(/' core:        ecore',7x,'node at',6x,'max at', &
                  7x,'c.t.p.   rho(r>rmt)')
40        enddo
4011   enddo

80  enddo
    if(ifiwv>0) write(ifiwv,"('9 9 9 9 9 --------! this is terminator for an atom section')")

    ! --- Write file with valence wavefunctions
    if (lplawv == 1) then
       write (str,'(''wf_'',a)') spid
       write (stdo,344) str
344    format(/' Write valence wavefunctions to plot file: ',a)
       !        ifi = fopna(str,-1,0)
       open(newunit=ifi,file=trim(str)//'.'//trim(sname))
       write (ifi,490) spid,rmax,rmt,nr,lmaxa,nr,1+nsp*(lmaxa+1)
490    format('# Free-atom wavefunctions (divided by r) for species ', &
            a/'# rmax=',f7.3,'   rmt=',f7.3,'   nr=',i5,'   lmax=',i3/ &
            '% rows ',i5,' cols ',i3)
       do  50  i=1,nr
          write (ifi,495) rofi(i),((psi(i,l,isp),l=0,lmaxa),isp=1,nsp)
495       format(f9.5,1p,16d14.5)
50     enddo
       close(ifi)
       !        call fclr(str,ifi)
    endif

!!! print out qbyl
    !      open('qinrmt.'//trim(sname))
    !      do  ib = 1, nbas
    !        ispec=ssite(ib)%spec
    !        write(ifqbyl,"(i3,10f12.6)") sspec(ispec)%lmxa, (sum(qbyl_rv(il,1:nsp,ib)),il=1,lmaxa)
    !      enddo
    !      close(ifqbyl)


  end subroutine pratfs


  subroutine ppratf(e,z,nr,nre,rofi,a,b,v,g,pzero,pmax,ctp)

    !- Find outermost node and maximum of wavefct
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   e     :wave function eigenvalue
    !i   z     :nuclear charge
    !i   nr    :number of radial mesh points
    !i   nre   :last point for which wf is calculated
    !i   rofi  :radial mesh points
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   v     :spherical potential (atomsr.f)
    !i   g     :normalized wave function times r
    !o Outputs
    !o   pzero :outermost node
    !o   pmax  :outermost maximum
    !o   ctp   :classical turning point
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: nr,nre
    double precision :: a,b,ctp,e,pmax,pzero,z,rofi(nr),v(nr),g(nr)
    ! ... Local parameters
    integer :: i,ir
    double precision :: g1,g2,rho1,rho2,rho3,x

    ! ... Find the classical turning point
    do  20 i = nr-1, 5, -1
       ir = i
       if (e > v(i)-2d0*z/rofi(i)) goto 21
20  enddo
21  g1 = e-v(ir) + 2d0*z/rofi(ir)
    g2 = e-v(ir+1) + 2d0*z/rofi(ir+1)
    ctp = rofi(nr)
    if (g1*g2 < 0d0) ctp = (rofi(ir)*g2-rofi(ir+1)*g1)/(g2-g1)

    ! ... Find the outermost node
    do  10  i = nre-1, 5, -1
       ir = i
       if (g(i)*g(i+1) < 0d0) goto 11
10  enddo
11  continue
    pzero = 0d0
    g1 = g(ir)
    g2 = g(ir+1)
    if (ir > 5) pzero = (rofi(ir)*g2-rofi(ir+1)*g1)/(g2-g1)

    ! ... Find the outermost maximum
    do  30  i = nre-2, 5, -1
       ir = i
       rho1 = g(i)*g(i)
       rho2 = g(i+1)*g(i+1)
       rho3 = g(i+2)*g(i+2)
       if (rho1 < rho2) goto 31
30  enddo
31  pmax = 0
    if (ir > 5) then
       x = -0.5d0*(rho3-rho1)/(rho1+rho3-2*rho2)
       pmax = b*(dexp(a*(ir+x))-1d0)
    endif
  end subroutine ppratf


  subroutine optfab(isw,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,nsp,v,rofi, &
       spid,itab,rtab,etab)
    use m_ftox
    !- Optimise a minimal smooth-Hankel basis for the free atom.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   isw   :1 constrain rsm to be <= rmt
    !i   z     :nuclear charge
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rmax  :muffin-tin radius, in a.u.
    !i   nrmt  :number of points between 0..rmt
    !i   rmt   :muffin-tin radius, in a.u.
    !i   lmxa  :muffin-tin l-cutoff
    !i   pl    :boundary conditions.  If Dl = log. deriv. at rmax,,
    !i         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   ql    :sphere moments
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   v     :spherical potential (atomsr.f)
    !i   rofi  :radial mesh points
    !i   spid  :species label
    !o Outputs
    !o   itab  :itab(l+1)=1  optimized wave function was found for this l
    !o   rtab  :smoothing radius for optimized wave function
    !o   etab  :energy for optimized wave function
    !l Local variables
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    integer :: lmxa,nr,nrmt,nsp,n0,isw
    double precision :: a,rmax,rmt,z
    parameter (n0=10)
    double precision :: rofi(1),v(nr,nsp),pl(n0,nsp),ql(3,n0,nsp)
    character spid*8
    ! ... Local parameters
    logical :: cmdopt
    character strn*80
    integer :: itab(n0,2)
    double precision :: rtab(n0,2),etab(n0,2)
    integer:: ipr , iprint , irep , isp , istife , istifr , jpr , &
         konfig , l ,  lplawv , nn , nrep
    real(8) ,allocatable :: g_rv(:)
    real(8) ,allocatable :: gp_rv(:)
    real(8) ,allocatable :: h_rv(:)
    real(8) ,allocatable :: psi_rv(:)
    double precision :: b,deh,deh0,dphi,dphip,drsm,drsm0,e1,e2,e3,eadd, &
         eaddx,eh,elim1,elim2,enew,enu,eval,p,phi,phip,pnu,qvl,radd, &
         raddx,rlim1,rlim2,rnew,rsm,stife,stifr,sume1,sume2,qrmt,rmtmx
    ipr = iprint()
    !      if (z .lt. 0.99d0) lrel = 0 !takao think lrel here was not used may2021
    if (lmxa > 8) call rx('optfab:  lmax too large')
    b = rmax/(dexp(a*nr-a)-1d0)
    do  80  isp = 1, nsp
       if(ipr>=20) write(stdo,ftox) &
            ' Optimise free-atom basis for species '//spid, 'Rmt=',ftof(rmt)
       allocate(h_rv(nr))
       allocate(g_rv(2*nr))
       allocate(gp_rv(2*nr*4))

       ! --- Parameters for minimisation ---
       drsm0 = 0.1d0
       ! takao makes "safer setting"
       ! akao
       !      rlim1 = 0.9d0
       !        rlim1 = 0.3d0 !original
       rlim1 = 0.5*rmt !original
       !      rlim1 = 0.8d0

       rlim2 = 2*rmt
       ! akao
       !      rlim2 = rmt+1d-6

       raddx = 0.2d0

       deh0  = 0.05d0
       ! akao
       elim1 = -2.5d0
       elim2 = -0.10d0
       !     elim2 = -0.20d0
       eaddx = 0.099d0
       jpr=0
       if (ipr >= 50) jpr=1

       ! --- Loop over bound valence states ---
       if (ipr >= 20) write (stdo,261)
       sume1 = 0d0
       sume2 = 0d0
       do  10  l = 0, lmxa
          itab(l+1,isp) = 0
          konfig = pl(l+1,isp)
          nn = konfig-l-1
          qvl = ql(1,l+1,isp)
          !   ... Get exact fa wavefunction, eigval, pnu at rmt
          call popta3 ( 0 , l , z , nn , rmt , nr , nrmt , rofi , v ( 1 &
               , isp ) , a , b , eval , pnu , g_rv )

          if (eval > 0d0) goto 10
          sume1 = sume1 + qvl*eval
          !   ... Potential parameters at MT sphere
          call popta4 ( l , z , rmt , nrmt , rofi , v ( 1 , isp ) , g_rv &
               , gp_rv , a , b , pnu , enu , p , phi , dphi , phip , dphip &
               )

          !          rsm = rmt
          rsm = rmt*.5
          eh = -1
          if (jpr > 0) write (stdo,340)
340       format('  L   parin    aux      E1       E2       E3', &
               '       stiff    Eout     parout')
          do  12  irep = 1, 50
             nrep = irep
             !     ... Get center energy
             call popta1 ( rsm , eh , l , z , rmt , nr , nrmt , rofi , h_rv &
                  , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip , dphip &
                  , e2 , qrmt )
             !            print *,' ttt center e=',irep,eh,e2
             !     ... Vary rsm
             !            drsm = drsm0
             !            call popta1 ( rsm + drsm , eh , l , z , rmt , nr , nrmt , rofi
             !     .      , h_rv , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip
             !     .      , dphip , e3 , qrmt )

             !            call popta1 ( rsm - drsm , eh , l , z , rmt , nr , nrmt , rofi
             !     .      , h_rv , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip
             !     .      , dphip , e1 , qrmt )

             !            call popta2(l,rsm,eh,drsm,e1,e2,e3,rlim1,rlim2,raddx,rnew,
             !     .      stifr,jpr)

             !     ... Vary eh
             deh = deh0
             !         if (eh+deh.gt.-0.01d0) deh=-eh-0.01d0
             call popta1 ( rsm , eh + deh , l , z , rmt , nr , nrmt , rofi &
                  , h_rv , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip &
                  , dphip , e3 , qrmt )
             call popta1 ( rsm , eh - deh , l , z , rmt , nr , nrmt , rofi &
                  , h_rv , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip &
                  , dphip , e1 , qrmt )
             call popta2(l,eh,rsm,deh,e1,e2,e3,elim1,elim2,eaddx,enew, &
                  stife,jpr)
             !            radd = rnew-rsm
             radd=0d0
             eadd = enew-eh
             !            rsm = rnew
             eh = enew
             if (dabs(radd) < 5d-3 .AND. dabs(eadd) < 5d-3) goto 90
12        enddo
90        continue

          !   ... End of iteration loop

          sume2 = sume2 + qvl*e2
          if (ipr >= 20) &
               write (stdo,260) l,nrep,rsm,eh,stifr,stife,e2,eval,pnu,qvl
260       format(i2,i4,2f8.3,1x,2f9.1,1x,2f10.5,f8.2,f7.2)
261       format(' l  it    Rsm      Eh     stiffR   stiffE', &
               '      Eval      Exact     Pnu    Ql')
          istifr = stifr+0.5d0
          istife = stife+0.5d0
          write (stdl,710) l,nrep,rsm,eh,istifr,istife,e2,eval,pnu,qvl
710       format('fa op',i2,i4,2f7.3,'  stf',2i6,'  ev',2f9.5, &
               '  pq',2f6.2)

          !   ... Possibly constrain rsm
          ! akao
          ! akao dec15 2010
          rmtmx= rmt*.5d0+1d-10
          if(rmtmx<.5d0) rmtmx=.5d0
          !          if(l>=2) rmtmx=rmt

          if (mod(isw,10) == 1 .AND. rsm > rmtmx) then
             if (ipr >= 20) &
                  write(stdo, &
                  '('' ...rsm exceeded rmtmx.. repeat with rsm= rmtmx ='',f8.5)') rmtmx
             !            rsm = rmtmx
             rsm = 0.5d0*rmtmx

             sume2 = sume2 - qvl*e2

             do  112  irep = 1,50
                nrep = irep
                print *,' ttt2 center e=',e2
                !     ... Get center energy
                call popta1 ( rsm , eh , l , z , rmt , nr , nrmt , rofi , h_rv &
                     , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip , dphip &
                     , e2 , qrmt )

                !     ... Vary eh
                deh = deh0
                !         if (eh+deh.gt.-0.01d0) deh=-eh-0.01d0
                call popta1 ( rsm , eh + deh , l , z , rmt , nr , nrmt , rofi &
                     , h_rv , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip &
                     , dphip , e3 , qrmt )

                call popta1 ( rsm , eh - deh , l , z , rmt , nr , nrmt , rofi &
                     , h_rv , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip &
                     , dphip , e1 , qrmt )

                call popta2(l,eh,rsm,deh,e1,e2,e3,elim1,elim2,eaddx,enew, &
                     stife,jpr)

                eadd = enew-eh
                eh = enew
                if (dabs(eadd) < 5d-3) goto 190
112          enddo
190          continue
             !   ... End of iteration loop

             sume2 = sume2 + qvl*e2
             if (ipr >= 20) &
                  write (stdo,260) l,nrep,rsm,eh,stifr,stife,e2,eval,pnu,qvl
             istife = stife+0.5d0
             write (stdl,710) l,nrep,rsm,eh,istifr,istife,e2,eval,pnu,qvl
          endif

          itab(l+1,isp) = 1
          rtab(l+1,isp) = rsm
          etab(l+1,isp) = eh
10     enddo


       if (ipr >= 20) write (stdo,320) sume1,sume2,sume2-sume1
320    format(' eigenvalue sum:  exact',f10.5,'    opt basis',f10.5, &
            '    error',f8.5)
       write (stdl,720) sume1,sume2,sume2-sume1
720    format('fa op sumev',f11.5,'   opt basis',f11.5,'   err',f9.5)

       ! i
       if (allocated(gp_rv)) deallocate(gp_rv)
       if (allocated(g_rv)) deallocate(g_rv)
       if (allocated(h_rv)) deallocate(h_rv)


80  enddo

    ! --- Make plot file ---
    !     lplawv=nglob('lplawv')
    lplawv = 0
    if (cmdopt('--plotwf',8,0,strn)) lplawv = 1
    if (lplawv == 1) then
       if (nsp == 2) call rx('optfab is not spinpol yet')
       allocate(psi_rv(nr*lmxa))

       call popta5 ( lmxa , rtab , etab , itab , z , pl , rmax , rmt &
            , nr , nrmt , rofi , psi_rv , v , g_rv , a , b , spid )

       if (allocated(psi_rv)) deallocate(psi_rv)

    endif

    ! i      call rlse (oh)

    ! i#error ERROR, try to release name= opsi ,but list=null at linenumber= 1031 list= (None)


  end subroutine optfab


  subroutine ftfalo(icst,z,a,nr,rmax,nrmt,rmt,lmxa,pnu,pz,rs3,eh3, &
       vmtz,nsp,v,rofi,spid)

    !- Fit value and slope of local orbitals to smoothed Hankel
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   icst  :1 constrain rsm to be <= rmt
    !i   z     :nuclear charge
    !i   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rmax  :muffin-tin radius, in a.u.
    !i   nrmt  :number of points between 0..rmt
    !i   rmt   :muffin-tin radius, in a.u.
    !i   lmxa  :muffin-tin l-cutoff
    !i   pl    :boundary conditions for valence wavefunctions.
    !i   pz    :boundary conditions for local orbital. pz=0 -> no loc. orb.
    !i         :10s digit controls how local orbital included in hamiltonian
    !i         :10s digit nonzero -> smooth Hankel tail is attached.
    !i   rs3   :minimum allowed smoothing radius in attaching Hankel tails
    !i         :to local orbitals
    !i   eh3   :Hankel energy when attaching Hankel tails to high-lying
    !i         :local orbitals
    !i   vmtz  :parameter used in attaching Hankel tails to local orbitals
    !i         :It is used as a constant shift to Hankel energies for the
    !i         :fitting of local orbitals to Hankel tails. Thus vmtz
    !i         :is an estimate for the potential at the MT radius.
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   v     :spherical potential (atomsr.f)
    !i   rofi  :radial mesh points
    !i   spid  :species label
    !o Outputs
    !l Local variables
    !r Remarks
    !u Updates
    !u   16 Jun 04 First created
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: lmxa,nr,nrmt,nsp,n0,icst
    double precision :: a,rmax,rmt,z,rs3,eh3,vmtz
    parameter (n0=10)
    double precision :: rofi(1),v(nr,nsp),pz(n0,nsp),pnu(n0,nsp)
    ! ... Local parameters
    logical :: cmdopt
    character spid*8, strn*80, flg(2)*1
    integer:: ipr , iprint , i , konfig , l , info , nn &
         , lplawv , loclo=-999 , nfit , isw
    real(8) ,allocatable :: g_rv(:)
    real(8) ,allocatable :: gp_rv(:)
    real(8) ,allocatable :: h_rv(:)

    real(8) ,allocatable :: psi_rv(:)

    double precision :: b,dasum,dphi,dphip,e2,eh,eval,p,phi,phip, &
         pnul,rsm,rsmin,rsmax,ekin
    !     emin and emax are the maximum allowed ranges in Hankel energies
    !     for the fitting of local orbitals to Hankel tails.
    double precision :: emin,emax,tphi
    !     For plotting wave functions
    integer :: itab(n0,2)
    double precision :: rtab(n0,2),etab(n0,2),pl(n0,nsp),qrmt
    data flg/'*',' '/
    ipr = iprint()
    if (lmxa > 8) call rx('ftfalo:  lmax too large')
    b = rmax/(dexp(a*nr-a)-1d0)
    nfit = 0

    if (dasum(lmxa+1,pz,1) == 0) return

    do  80  i = 1, nsp
       allocate(h_rv(nr))

       allocate(g_rv(2*nr))

       allocate(gp_rv(2*nr*4))


       ! --- Loop over local orbitals ---
       !      sume1 = 0d0
       !      sume2 = 0d0
       do  10  l = 0, lmxa

          itab(l+1,i) = 0
          pnul = pnu(l+1,i)
          pl(l+1,i) = pnu(l+1,i)
          konfig = mod(pz(l+1,1),10d0)

          !       Skip all but local orbitals with tails attached
          if (pz(l+1,1) < 10) goto 10

          !       Case local orbital deeper than valence
          if (int(pnul-1) == int(mod(pz(l+1,1),10d0))) then
             loclo = 1
             !         Not needed, actually, since overwritten by popta3
             !         pnul = mod(pz(l+1,1),10d0)

             !       Case local orbital higher than the valence state
          elseif (int(pnul+1) == int(mod(pz(l+1,1),10d0))) then
             pnul = mod(pz(l+1,1),10d0)
             loclo = 0

             !       Local orbital neither one: error
          else
             call fexit3(-1,111,' Exit -1 freeat, l=%i:  sc '// &
                  'PZ=%d incompatible with valence P=%;3d',l,pz(l+1,1),pnul)
          endif

          !       Skip high-lying local orbitals unless specifically sought
          if (loclo == 0 .AND. .NOT. cmdopt('--getallloc',11,0,strn)) &
               goto 10

          nfit = nfit + 1
          if (nfit == 1) then
             call info2(20,1,0, &
                  ' Fit local orbitals to sm hankels, species '//spid// &
                  '%a, rmt=%;7g',rmt,0)
             if (ipr >= 20) write (stdo,261)
          endif

          !   ... Get exact fa wavefunction, eigval, pnu_l at rmt
          if (loclo == 1) then
             nn = konfig-l-1
             call popta3 ( 0 , l , z , nn , rmt , nr , nrmt , rofi , v ( 1 &
                  , i ) , a , b , eval , pnul , g_rv )

             !       Finish if in future, need w.f. at r>rmt
             !        else
             !          call popta3(1,l,z,nn,rmt,nr,nrmt,rofi,v(1,i),a,b,eval,
             !     .      pnul,w(og))
          endif
          pl(l+1,i) = pnul

          !   ... Potential parameters at MT sphere
          call popta4 ( l , z , rmt , nrmt , rofi , v ( 1 , i ) , g_rv &
               , gp_rv , a , b , pnul , eval , p , phi , dphi , phip , dphip &
               )


          !   ... Set conditions on envelope functions ... For now
          rsmin = rs3
          rsmax = 5
          if (icst == 1) rsmax = rmt
          !       Use r->infty value for energy
          eh = min(-.02d0,eval)

          !   ... Match Hankel to phi,dphi
          !        rsm = rsmin
          !        emax = -.02d0
          !        emin = -5d0
          !        call mtchre(100,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,phi,
          !     .    dphi,rsm,eh,ekin,info)

          !   ... Match slope and K.E. of Hankel to phi,dphi
          tphi = eval - (v(nrmt,i)-2*z/rmt)
          rsm = 0
          eh = min(eval-vmtz,-.02d0)
          emax = -.02d0
          emin = -10d0
          !       if (ipr .ge. 20) call pshpr(max(ipr,50))
          call mtchre(003,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,tphi, &
               dphi,rsm,eh,ekin,info)
          !       if (ipr .ge. 20) call poppr
          !       Match failed ... turn up verbosity and repeat for info
          if (info == -1) then
             call info2(0,2,1, &
                  ' *** ftfalo (fatal) cannot fit smooth Hankel to w.f.'// &
                  ' class '//spid// &
                  '%N ... possibly reduce RS3 (current value = %,1d)',rs3,0)
             call pshpr(max(ipr,110))
             call mtchr2(1,l,emin,emax,(emin+emax)/2, &
                  rmt,phi,dphi,rsmin,eh,ekin,i)
             call poppr
             !         call pshpr(max(ipr,110))
             !         call mtchre(103,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,tphi,
             !    .      dphi,rsm,eh,ekin,info)
             call fexit2(-1,111, &
                  ' Exit -1 : ftfalo : failed to match log der=%,1;3d'// &
                  ' to envelope, l=%i',dphi/phi,l)
          endif

          !  ... Get energy of this wave function
          call popta1 ( rsm , eh , l , z , rmt , nr , nrmt , rofi , h_rv &
               , v ( 1 , i ) , a , b , eval , p , phi , dphi , phip , dphip &
               , e2 , qrmt )


          if (ipr >= 20) &
               write (stdo,260) l,rsm,eh,qrmt,e2,eval,pnul,tphi,ekin, &
               flg(2-isw(dabs(ekin-tphi) > 1d-5))

260       format(i2,2f8.3,3f10.5,f9.3,2f10.5,a1,f10.5)
261       format(' l    Rsm     Eh     Q(r>rmt)   Eval', &
               '      Exact      Pnu     K.E.    fit K.E.')

          itab(l+1,i) = 1
          rtab(l+1,i) = rsm
          etab(l+1,i) = eh

10     enddo
       ! akao moved deallocation to here
       deallocate(gp_rv)
       deallocate(g_rv)
       deallocate(h_rv)
80  enddo

    ! --- Make plot file ---
    !     lplawv=nglob('lplawv')
    lplawv = 0
    if (cmdopt('--plotwf',8,0,strn)) lplawv = 1
    if (lplawv == 1) then
       if (nsp == 2) call rx('optfab is not spinpol yet')
       allocate(psi_rv(nr*lmxa))

       call popta5 ( lmxa , rtab , etab , itab , z , pl , rmax , rmt &
            , nr , nrmt , rofi , psi_rv , v , g_rv , a , b , spid )

       deallocate(psi_rv)

    endif


  end subroutine ftfalo


  subroutine popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,v,a,b,enu,p, &
       phi,dphi,phip,dphip,eval,qrmt)
    use m_hansr,only:hansmd
    !- Calculate expectation value for smooth Hankel
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   rsm   :smoothing radius of basis function
    !i   eh    :energy of basis function
    !i   l     :l quantum number
    !i   z     :nuclear charge
    !i   rmt   :muffin-tin radius
    !i   nr    :number of radial mesh points
    !i   nrmt  :number points to muffin-tin radius
    !i   rofi  :radial mesh points
    !i   h     :work array
    !i   v     :spherical potential (atomsr.f)
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   enu   :enu's for making charge density
    !i   p     :<gp**2> (potential parameter)
    !o   phi   :wave function at rmt
    !o   dphi  :radial derivative of of phi at rmt
    !o   phip  :energy derivative of phi
    !o   dphip :radial derivative of dphi
    !o Outputs
    !o   eval  :expectation value
    !o   qrmt  :fraction of (wave function)^2 for r>rmt
    !r Remarks
    !u Updates
    !u   24 Sep 04 return qrmt
    !u   16 Jun 04 Adapted to new hansmd, mtchae
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: l,nr,nrmt
    double precision :: a,b,dphi,dphip,eh,enu,eval,p,phi,phip,rmt,rsm,z, &
         rofi(nr),h(nr),v(nr),qrmt
    ! ... Local parameters
    integer :: i
    double precision :: alfa,beta,det(1),drdi,hlap, &
         hum,hum1,hum2,r,sum,sum1,sum2,tum2,vum2,wt
    !     double precision xi(0:20)

    double precision :: hs(0:l),dhs(0:l),ddhs(0:l)

    !     pi = 4d0*datan(1d0)
    !     asm = 1d0/rsm
    !     lp1 = l+1

    ! ... Integrals over smooth Hankel on mesh
    !     gfac = (asm*asm/pi)**1.5d0 * dexp(eh*rsm*rsm/4d0)
    !     ta2 = 2d0*asm*asm
    tum2 = 0d0
    sum2 = 0d0
    vum2 = 0d0

    do  10  i = nrmt, nr
       r = rofi(i)
       !   ... Make r*h and r Laplacian h, including L^2
       call hansmd(2,r,eh,rsm,l,hs,dhs,ddhs,det,det,det)
       h(i) = hs(l)*r
       hlap = ddhs(l)*r
       !C      Old : r*h and r Laplacian h, including L^2
       !C      h = r*radial part of sm. Hankel
       !       call hansmr(r,eh,asm,xi,l)
       !       h(i) = xi(l)*(r**lp1)
       !C      radial part of Gaussian
       !       gl = gfac * dexp(-asm*asm*r*r) * ta2**l * (r**lp1)
       !C      r * (nabla_r - l(l+1)/r^2) h_l
       !       hlap = -4d0*pi*gl - eh*h(i)

       !  ...  Accumulate <h h>, <h v h>, <h -nabla h>
       wt = 2*(mod(i+1,2)+1)/3d0
       if (i == nrmt .OR. i == nr) wt = 1d0/3d0
       drdi = a*(r+b)
       sum2 = sum2 + wt*drdi*h(i)*h(i)
       vum2 = vum2 + wt*drdi*h(i)*h(i)*(v(i)-2d0*z/r)
       tum2 = tum2 + wt*drdi*h(i)*(-hlap)
10  enddo
    hum2 = tum2+vum2

    ! --- BC's: match phi,phidot to envelope at RMT ---
    call mtchae(0,rsm,eh,l,rmt,phi,dphi,phip,dphip,alfa,beta)
    !C    OLD matching
    !C    Match value, slope fl,dfl to linear combination of phi,phidot
    !     call hansmr(rmt,eh,asm,xi,l+1)
    !C    Value and radial derivative of h (JMP 39, 3393, Eq. 4.7)
    !     fl = xi(l)*rmt**l
    !     flp1 = xi(l+1)*rmt**(l+1)
    !     dfl = l*fl/rmt-flp1
    !C    Match fl,dfl to linear combination of phi,phidot
    !C    Use  phi=phi(R); phip=phidot(R) dphi=phi'(R); dphip=phidot'(R)
    !C    (phi  phip ) (alpha)   (fl )    (alpha)    1  (dphip -phip) (fl )
    !C    (          ) (     ) = (   ) -> (     ) = --- (           ) (   )
    !C    (dphi dphip) (beta )   (dfl)    (beta )   det (-dphi  phi ) (dfl)
    !     det = phi*dphip-dphi*phip
    !     alfa = (fl*dphip-dfl*phip)/det
    !     beta = (dfl*phi-fl*dphi)/det

    !     O = alpha^2 <phi | phi> + beta^2 <phidot | phidot>
    sum1 = alfa*alfa + beta*beta*p
    hum1 = alfa*alfa*enu + alfa*beta + beta*beta*enu*p

    sum = sum1+sum2
    hum = hum1+hum2
    eval = hum/sum

    qrmt = sum2/sum
  end subroutine popta1


  subroutine popta2(l,x0,y0,dx,e1,e2,e3,xmin,xmax,xshx,xnew,stiff, &
       jpr)

    !- Find minimum from three values
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   l     :angular momentum
    !i   x0    :starting value
    !i   y0    :used for printout
    !i   dx    :excursion in x for numerical differentiation
    !i   e1    :function value at x0-dx
    !i   e2    :function value at x0
    !i   e3    :function value at x0+dx
    !i   xmin  :boundary: estimated minimum must be >= xmin
    !i   xmax  :boundary: estimated minimum must be <= xmax
    !i   xshx  :maximum step size
    !i   jpr   :printout verbosity
    !o Outputs
    !o   xnew  :new estimate for the minimum
    !o   stiff :estimated curvature
    !r Remarks
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: jpr,l
    double precision :: dx,e1,e2,e3,stiff,x0,xmax,xmin,xnew,xshx,y0
    ! ... Local parameters
    integer :: ie0
    double precision :: a,aa,b,c,ee1,ee2,ee3,een,enew,xadd

    c = e2
    b = (e3-e1)/(2*dx)
    a = (e1+e3-2*e2)/(2*dx*dx)
    if (a <= 0d0) then
       xadd = -xshx
       enew = e1
       if (e3 < e1) xadd = xshx
       if (e3 < e1) enew = e3
    else
       xadd = -b/(2*a)
       enew = a*xadd*xadd + b*xadd + c
    endif
    aa = 2*1d3*a

    if (xadd > xshx)  xadd = xshx
    if (xadd < -xshx) xadd = -xshx
    xnew = x0+xadd
    if (xnew > xmax) xnew = xmax
    if (xnew < xmin) xnew = xmin

    ie0 = e2
    ee1 = 1d3*(e1-ie0)
    ee2 = 1d3*(e2-ie0)
    ee3 = 1d3*(e3-ie0)
    een = 1d3*(enew-ie0)
    stiff = aa

    if (jpr > 0) write (stdo,810)l,x0,y0,ee1,ee2,ee3,aa,een,xnew
810 format(i3,f8.3,f8.3,f10.3,2f9.3,f9.1,f10.3,f8.3,a)

  end subroutine popta2

  subroutine popta3(mode,l,z,nn,rmt,nr,nrmt,rofi,v,a,b,evl,pnu,g)

    !- Get exact fa wavefunction, eigval, pnu at Rmt.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 boundary condition is val,slo = 0 at nr
    !i         :1 boundary condition is that w.f. satisfy pnu at nrmt
    !i         :  (under development)
    !i   l     :angular momentum
    !i   z     :nuclear charge
    !i   rmax  :sphere radius
    !i   rmt   :muffin-tin radius, in a.u.
    !i   nr    :number of radial mesh points
    !i   nrmt  :number of radial mesh points to rmt
    !i   rofi  :radial mesh points
    !i   v     :spherical potential (atomsr.f)
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    ! o Inputs/Outputs
    ! o  nn    :number of nodes (input mode 0; output mode 1)
    ! o  pnu   :boundary condition at rmt (output mode 0; input mode 1)
    !o Outputs
    !o   g     :normalized wave function times r
    !o   evl   :eigenvalue
    !r Remarks
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: mode,l,nn,nr,nrmt
    double precision :: a,b,evl,pnu,rmt,z,rofi(nr),v(nr),g(nr*2)
    ! ... Local parameters
    integer :: nre,konfig,nri,nn2
    double precision :: d0l,p0l,dphi,drdi,du,eb1,eb2,g1,g2,g3,g4,g5,pi, &
         slo,slou,sum,tol,val,valu,dnu
    pi = 4d0*datan(1d0)

    eb1 = -50 !ccctakao 30
    eb2 = 20
    tol = 1d-10
    val = 1d-30
    slo = -val
    evl = -0.5d0
    nri = nr
    if (mode == 1) then
       konfig = pnu
       nn = konfig-l-1
       dnu = dtan(pi*(0.5d0-pnu))
       val = rmt
       slo = dnu+1d0
       nri = nrmt
    endif
    call rseq(eb1,eb2,evl,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nri,nre)
    if (mode == 1) then
       !       integration becomes rather strange for r>>rmt.
       !       Need to truncate radius.
       call rsq1(nri,evl,l,z,v,nr,g,val,slo,nn2,a,b,rofi,nr)
       !       call prrmsh('g',rofi,g,nr,nr,1)
       call rx('not finished mode 1')
    endif
    g1 = g(nrmt-2)
    g2 = g(nrmt-1)
    g3 = g(nrmt)
    g4 = g(nrmt+1)
    g5 = g(nrmt+2)
    drdi = a*(rmt+b)
    valu = g3
    slou = (-2*g5+16*g4-16*g2+2*g1)/(24d0*drdi)
    du   = rmt*slou/valu
    dphi = du-1
    pnu  = nn+l+1 + (0.5d0-datan(dphi)/pi)

    ! ... Don't set too low..
    d0l = l
    p0l = nn+l+1 + 0.5d0-datan(d0l)/pi
    p0l = nn+l+1 + 0.1d0
    if (pnu < p0l) then
       write (stdo,145) l,pnu,p0l
145    format(' l=',i1,'  increase Pnu=',f8.3,'  to ',f8.3)
       pnu = p0l
    endif

  end subroutine popta3

  subroutine popta4(l,z,rmt,nrmt,rofi,v,g,gp,a,b,pnu,enu,p,phi,dphi, &
       phip,dphip)

    !- Potential parameters at MT sphere
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   l     :angular momentum
    !i   z     :nuclear charge
    !i   rmt   :muffin-tin radius, in a.u.
    !i   nrmt  :number of radial mesh points to rmt
    !i   rofi  :radial mesh points
    !i   v     :spherical potential (atomsr.f)
    !i   g     :normalized wave function times r
    !i   gp    :energy derivative(s) of g
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   enu   :enu's for making charge density
    !o Outputs
    !o   phi   :wave function at rmt
    !o   dphi  :radial derivative of of phi at rmt
    !o   phip  :energy derivative of phi
    !o   dphip :radial derivative of dphi
    !o   p     :<gp**2> (potential parameter)
    !r Remarks
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: l,nrmt
    double precision :: a,b,dphi,dphip,enu,p,phi,phip,pnu,rmt,z
    double precision :: rofi(nrmt),v(nrmt),g(nrmt),gp(nrmt,4)
    ! ... Local parameters
    integer :: konfig,nn,nre
    double precision :: dnu,eb1,eb2,pi,slo(5),sum,tol,val(5)
    pi = 4d0*datan(1d0)
    eb1 = -50 !ccctakao -30
    eb2 = 20
    tol = 1d-10

    konfig = pnu
    nn = konfig-l-1
    dnu = dtan(pi*(0.5d0-pnu))
    val(1) = rmt
    slo(1) = dnu+1d0
    enu=-0.5d0

    call rseq(eb1,eb2,enu,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nrmt, &
         nre)
    val(1) = val(1)/dsqrt(sum)
    slo(1) = slo(1)/dsqrt(sum)

    !      call phidot(z,l,v,enu,a,b,rofi,nrmt,g,val,slo,tol,nn,gp,phi,dphi,
    !     .  phip,dphip,p)

    call phidx(1,z,l,v,0d0,0d0,rofi,nrmt,2,tol,enu,val,slo,nn,g,gp, &
         phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
    !     dphip = (slo(2)-phip)/rmt


    !|    write(stdo,200) l,enu,p,phi,dphi,phip,dphip
    ! 200 format(' PP',i2,'  e',f10.5,'  p',f10.5,'  bc',4f10.5)

  end subroutine popta4


  subroutine popta5(lmax,rtab,etab,itab,z,pl,rmax,rmt,nr,nrmt, &
       rofi,psi,v,g,a,b,spid)
    use m_ext,only:sname
    use m_hansr,only:hansmr
    !- Write wave functions to plot file
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   lmax  :maximum l for a given site
    !i   rtab  :smoothing radii for wavefunction, each l
    !i   etab  :smoothed hankel energies for wavefunction, each l
    !i   itab  :1 if a wave function calculated, 0 if not
    !i   z     :nuclear charge
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
    !i         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   rmax  :muffin-tin radius, in a.u.
    !i   rmt   :muffin-tin radius
    !i   nr    :number of radial mesh points
    !i   nrmt  :number points to muffin-tin radius
    !i   rofi  :radial mesh points
    !i   psi   :wave function tabulated on the rofi mesh
    !i   v     :spherical potential (atomsr.f)
    !i   g     :normalized wave function times r
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   spid
    !o Outputs
    !    wave functions written to disk
    !l Local variables
    !l         :
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: itab(0:*),lmax,nr,nrmt,n0
    parameter (n0=10)
    double precision :: a,b,rmax,rmt,z,rtab(0:*),etab(0:*), &
         rofi(nr),psi(nr,0:*),g(nr,2),v(nr),pl(0:n0-1)
    character spid*8
    ! ... Local parameters
    integer :: i,ifi,konfig,l,lp1,m,n,nn,nre
    integer :: ltab(n0)
    double precision :: asm,dfl,drdi,eb1,eb2,eh,evl,fac,fl,flp1,r,rsm, &
         slo,sum1,sum2,tol,val,wt,xi(0:20)
    character str*32

    eb1 = -20
    eb2 = 20
    tol = 1d-8
    n = 0

    do  10  l = 0, lmax
       if (itab(l) == 0) goto 10
       n = n+1
       ltab(n) = l
       lp1 = l+1
       rsm = rtab(l)
       eh = etab(l)
       asm = 1d0/rsm
       konfig = pl(l)
       nn = konfig-l-1

       ! ...   Smooth hankel fct outside rmt
       sum2 = 0d0
       do  12  i = nrmt, nr
          r = rofi(i)
          call hansmr(r,eh,asm,xi,l)
          psi(i,n) = xi(l)*(r**lp1)
          wt = 2*(mod(i+1,2)+1)/3d0
          if (i == nrmt .OR. i == nr) wt=1d0/3d0
          drdi = a*(r+b)
          sum2 = sum2 + wt*drdi*psi(i,n)**2
12     enddo

       ! ...   Attach numerical solution inside MT sphere
       call hansmr(rmt,eh,asm,xi,l+1)
       fl   = xi(l)*rmt**l
       flp1 = xi(l+1)*rmt**(l+1)
       dfl  = l*fl/rmt-flp1
       val = rmt*fl
       slo = rmt*dfl+fl
       evl = -0.5d0
       call rseq(eb1,eb2,evl,tol,z,l,nn,val,slo,v, &
            g,sum1,a,b,rofi,nrmt,nre)
       fac = val/(g(nrmt,1)*dsqrt(sum1+sum2))
       do  14  i = 1, nrmt
          psi(i,n) = fac*g(i,1)
14     enddo
       fac = 1d0/dsqrt(sum1+sum2)
       do  16  i = nrmt+1,nr
          psi(i,n) = psi(i,n)*fac
16     enddo
10  enddo

    ! ... Write the plot file
    write (str,'(''wfa_'',a)') spid
    write (stdo,344) str
344 format(/' Write fit wavefunctions to plot file: ',a)
    !      ifi = fopna(str,-1,0)
    open(newunit=ifi,file=trim(str)//'.'//trim(sname))
    write (ifi,490) spid,rmax,rmt,(ltab(i),i=1,n)
490 format('# Free-atom opt basis (divided by r) for species ', &
         a/'# rmax=',f7.3,'   rmt=',f7.3,'   l=',8i2)
    write (ifi,'(''% rows '',i5,'' cols '',i3)') nr,n+1
    do  30  i = 1, nr
       write (ifi,495) rofi(i),(psi(i,m),m=1,n)
495    format(f9.5,1p,8d14.5)
30  enddo
    close(ifi)
  end subroutine popta5


  subroutine fctail(nr,nrmt,a,b,rsm,rofi,rhoc,c,eh)
    use m_lmfinit,only: nsp


    !- Fit one Hankel to tail of core density.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nr    :number of radial mesh points
    !i   nrmt  :number points to muffin-tin radius
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   rsm   :smoothing radius
    !i   rofi  :radial mesh points
    !i   rhoc  :core density
    !o Outputs
    !o   c     :coefficient to fit of rhoc(spin+)+rhoc(spin-)
    !o   eh    :energy
    !l Local variables
    !l   rmt   :muffin-tin radius
    !r Remarks
    !b Bugs
    !b   Should this be fit to smoothed function??
    !u Updates
    !u   19 Apr 02 Make rmt a local variable.
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: nr,nrmt
    double precision :: a,b,c,eh,rmt,rsm,rofi(nr),rhoc(nr,2)
    ! ... Local parameters
    integer :: i,ipr
    double precision :: ak1,akap,fit,q,q0,r,s,v0,wt
    character sout*80
    call getpr(ipr)
    rmt = rofi(nrmt)
    q0 = 0d0
    do  10  i = nrmt, nr
       r = rofi(i)
       wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
       if (i == nrmt .OR. i == nr) wt = a*(r+b)/3d0
       q0 = q0 + wt*(rhoc(i,1)+rhoc(i,nsp))/(3-nsp)
10  enddo
    v0 = (rhoc(nrmt,1)+rhoc(nrmt,nsp))/(3-nsp)/(rmt*rmt)
    write(stdo,"('conf: Core rhoc(rmt)=',f9.6,' spillout=',f9.6)")v0,q0
    !      write(stdo,"(a)")'conf:-----------------------------------------------------'
    !      sout = ' '
    !      call awrit3('%?#(n>=30)#%N## coretail: q=%;3g, spill out=%;3g.',
    !     .sout,len(sout),0,ipr,v0,q0)

    !      if (ipr .ge. 20) write (stdo,339) v0,q0
    !  339 format(/' coretail:  rho(rmt)=',f12.8,'   charge=',f12.8)
    if (dabs(q0) < 1d-6) then
       c = 0d0
       eh = -1d0
       return
    endif

    ! ... Parameters of hankel fct
    s = dsqrt(rmt**4 * v0**2 + 4*rmt*q0*v0)
    ak1 = (rmt*rmt*v0+s)/(2d0*q0)
    !     ak2 = (rmt*rmt*v0-s)/(2d0*q0)
    !|      write(stdo,975) ak1,ak2
    !|  975 format('ak1,ak2=',2f14.8)
    akap = ak1
    c = rmt*v0*dexp(akap*rmt)
    eh = -akap*akap

    if (ipr >= 20) then
       write(stdo,"(' Fit with Hankel e=',f10.6,' coeff=',f10.6)")eh,c
       !        call awrit2('%a  Fit with Hankel e=%;5g  coeff=%;5g',sout,
       !     .  len(sout),-stdo,eh,c)
    endif

    ! ... Test
    if (ipr > 30) then
       write (stdo,501)
       q = 0d0
       do  20  i = nrmt, nr
          r = rofi(i)
          wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
          if (i == nrmt .OR. i == nr) wt = a*(r+b)/3d0
          fit = c*dexp(-akap*r)*r
          q = q+wt*fit
          if ((rhoc(i,1)+rhoc(i,nsp))/(3-nsp) < 1d-8) goto 90
          if (mod(i,5) == 0 .OR. i == nrmt) &
               write (stdo,500) r,(rhoc(i,1)+rhoc(i,nsp))/(3-nsp),fit
500       format(f12.6,2f14.8)
501       format(6x,'r',12x,'rhoc',10x,'fit')
20     enddo
90     continue
       !|      v=c*dexp(-akap*rmt)/rmt
       !|      write(stdo,885) q,q0,v,v0
       !|  885 format('q,q0,v,v0=',4f14.8)
    endif

    ! ... look at smoothed core..
    !|      rg=0.4
    !|      qc=36
    !|      sum0=-c*dexp(eh*rg*rg/4d0)/eh
    !|      cg=qc-sum0
    !|      write(stdo,888) qc,sum0,cg
    !|  888 format(' qcore=',f10.4,'  sum0=',f12.6,'   cg=',f12.6)
    !|      ag=1d0/rg
    !|      fac=4d0*pi*(ag*ag/pi)**1.5d0
    !|      q=0d0
    !|      do i=1,nr
    !|        r=rofi(i)
    !|        wt=2*(mod(i+1,2)+1)*a*(r+b)/3d0
    !|        if (i.eq.1 .or. i.eq.nr) wt=a*(r+b)/3d0
    !|        call hansmr(r,eh,ag,xi,1)
    !|        fit=c*xi(0)*r*r + cg*fac*dexp(-ag*ag*r*r)*r*r
    !|        if (rhoc(i).gt.1d-10) write(49,490) r,rhoc(i),fit
    !|  490   format(f12.6,2f16.8)
    !|        q=q+wt*fit
    !|      enddo
    !|      write(stdo,965) q
    !|  965 format(' integral over smoothed core:',f10.5)

  end subroutine fctail


end module m_freeat
