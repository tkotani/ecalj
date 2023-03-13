module m_freeat !free-standing spherical atom calculaitons for initial contdition
  use m_lgunit,only: stdo,stdl
  use m_rseq,only: rseq
  public freeat,freats
  private
contains
  subroutine freeat() !For all species, we make free atom self-consistent, and get density to files.
    use m_ext,only:sname
    use m_lmfinit,only: smalit,lxcf,ham_seref,nsp,nspec, sspec=>v_sspec,& 
         idmod,slabl,vmtz,eref,rs3,eh3,nmcore,coreh,coreq,pnux=>pnusp,pzx=>pzsp,qnu
    use m_ftox
    !Inputs  are module variables of m_lmfinit
    !Outputs are via iofa, atmpnu are pnu (logarismic derivatives of atoms).
    ! Memo
    !   ccof  :coefficient to fit of core tail to smoothed Hankel
    !   ceh   :energy of core tail to smoothed Hankel
    !   sumtc :core kinetic energy
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ifi,iprint,is,nglob,nr,nrmt,nrmx, n0,nkap0,nxi,nxi0,nrmix,igets,lmxa,kcor,lcor,iofa, i_dum,ifives,ifiwv
    character(8) :: spid,chole(8)
    parameter ( nrmx=1501, nxi0=10, n0=10, nkap0=3)
    real(8) :: qc,ccof,ceh,z,rmt,rfoca,qcor(2),a,sumec, sumtc,seref,dgets,dgetss,etot
    real(8) :: hfc(nxi0,2),exi(nxi0),hfct(nxi0,2), v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),rofi(nrmx*2)
    real(8) :: pnu(n0,2),pz(n0,2),qat(n0,2), rtab(n0,2),etab(n0,2),rsmfa
    character strn*120
    logical :: cmdopt
    open(newunit=ifi,file='atm.'//trim(sname))
    open(newunit=ifives,file='vesintatm.'//trim(sname)//'.chk')
    open(newunit=ifiwv,file='veswavatm.'//trim(sname)//'.chk')
    hfct = 0d0
    do  is = 1, nspec ! Takao found that Li requires -0.5 to have positive smooth rho.
       if(sspec(is)%z<3.5) then !At least for Li, fitting is not good (negative smooth rho).
          exi(1:7) = [real(8):: -0.5, -1,-2,-4,-6,-9,-15]
          nxi = 7
       else !this is original setting. For CrN. This is necessary(maybe long range cutoff is required).
          exi(1:6) = [real(8)::       -1,-2,-4,-6,-9,-15]
          nxi = 6
       endif
       nrmix= smalit
       spid = slabl(is) 
       rfoca= sspec(is)%rfoca
       qcor = coreq(:,is)
       chole= coreh(is)
       call gtpcor(is,kcor,lcor,qcor)
       z   = sspec(is)%z
       rmt = sspec(is)%rmt
       a   = sspec(is)%a
       nrmt= sspec(is)%nr
       rsmfa=.5d0*rmt            ! moved to here 2022-6-27
       if (z == 0 .AND. rmt == 0) cycle !floating orbital
       pnu(:,1)=  pnux(1:n0,1,is) 
       if(nsp==2) pnu(:,2)= pnu(:,1)
       qat(1:n0,1:nsp)=  qnu(1:n0,1:nsp,is) 
       lmxa = sspec(is)%lmxa
       pz(:,1) =  pzx(1:n0,1,is) 
       if(nsp==2) pz(:,2)= pz(:,1) !       write(6,ftox)'xxx isp pz=',is,ftof(pz(1:lmxa+1,1),6)
       write(stdo,"(a)")'freats:'
       call freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,kcor,lcor,qcor, &
            nrmix,1,lxcf,z,rmt,a,nrmt,pnu,pz,qat,rs3(is),eh3(is),vmtz(is),& !rcfa(:,is), &
            idmod(:,is),lmxa,eref(is),rtab,etab,hfc,hfct,nr,rofi,rho,rhoc,qc,ccof, &
            ceh,sumec,sumtc,v,etot,nmcore(is),ifives,ifiwv)
       print *,'end of freats: spid nmcore=',spid,nmcore(is)
       if (iprint()>40) write(stdo,"(a)")' write free atom data for species  '//trim(spid)
       if (nsp == 2 .AND. nr > nrmt) then
          call dcopy(nrmt,rho(1+nr),1,rho(1+nrmt),1)
          call dcopy(nrmt,rhoc(1+nr),1,rhoc(1+nrmt),1)
          call dcopy(nrmt,v(1+nr),1,v(1+nrmt),1)
       endif
       i_dum = iofa(spid,nxi0,nxi,exi,hfc,hfct,rsmfa,z,rmt,a,nrmt,qc,ccof,ceh,sumtc,rho,rhoc,v,ifi,'write')
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
       nrmix,lwf,lxcf,z,rmt,a,nrmt,pnu,pz,qat,rs3,eh3,vmtz, & !,rcfa removed 2023feb. rnatm meaninful?
       idmod,lmxa,eref,rtab,etab,hfc,hfct,nr,rofi,rho,rhoc,qc,ccof,ceh, &
       sec,stc,v,etot,nmcore,ifives,ifiwv)
    use m_lmfinit,only: nsp,lrel
    use m_ftox
    use m_ext,only:sname
    use m_getqvc,only: getqvc
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
    !i   lxcf:selects local exchange-correlation functional
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
    integer :: nrmx,nrmt,is,nxi0,nxi,nrmix,lwf,lxcf,n0,kcor,lcor
    parameter (nrmx=1501,n0=10)
    character(8) :: spid
    real(8) :: rsmfa,rfoca,qc,ccof,ceh,sec,stc,z,rmt,a,eref, v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),hfc(nxi0,*),hfct(nxi0,*), &
         exi(*),rtab(n0,2),etab(n0,2),rofi(nrmx*2),rs3,eh3,vmtz,qcor(2)
    logical :: cmdopt
    integer :: ncmx,nvmx
    parameter (ncmx=50, nvmx=20)
    integer :: idmod(n0)
    character str*8,strn*32
    real(8) :: rmax,b,etot,dq, ec(ncmx),ev(nvmx),sumev,vrmax(2),exrmax(2),ekin,utot,rhoeps, &
         amgm,rhrmx,qvt,qtot,qct,qvin,qcin,r,wt0,wt1,qtt,qtin,pnul, pzl,pnu(n0,2),qat(n0,2),pl(n0,2),rhoin(nrmx*2), &
         rhot(nrmx*2),ql(3,n0,2),pz(n0,2) ,qatbk(n0,2) !,rcfa(2)
    integer :: i,ifi,ipr,iprint,isp,isw,l, lgrad,lmxa,lplfa,nitmax,nmix,nr,lplawv !,irchan(n0)
    real(8) ,allocatable :: g_rv(:)
    real(8) ,allocatable :: psi_rv(:)
    integer :: itab(n0,2)
    integer ::iwdummy
    integer:: ipl, ipz,iplx,iplv,ix, nmcore, ifives,ifiwv,iplz,miplzipl
    real(8):: qcc,qzz, qvalv,qvaltot,qval,vsum,qelectron,qvalz
    real(8),allocatable:: plplus(:,:),qlplus(:,:)
    logical::lfrz=.false.
    ipr    = iprint()
    lgrad  = lxcf/10
    if(ipr>49) then
       print *
       do i = 1, nsp
          do l = 0, lmxa
             write(stdo,"(a,2i3,2f10.3)")'ttt: pnu qat=',i,l,pnu(l+1,i),qat(l+1,i)
          enddo
       enddo
    endif
    ! --- Get species data ----
    pl=0d0
    ql=0d0
    allocate(plplus(0:lmxa,nsp),qlplus(0:lmxa,nsp),source=0d0)! plplus, qlplus are higher ones (when PZ exist).
    ! print *,'NOTE: when we have two valence: P and PZ, We assume eigen(PZ) is deeper than eigen(P).'
    nsploop: do  101  i = 1, nsp
       lloop: do  10  l = 0, lmxa
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
                call rx3('freeat: sc incompatible with valence (you may need to add P simultaneously) l pz pnuzl=',l,pzl,pnul)
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
10     enddo lloop
101 enddo nsploop
    write(stdo,*)
    write(stdo,"(a)")'conf:------------------------------------------------------'
    write(stdo,"(a,a)")'conf:SPEC_ATOM= '//trim(spid),' : --- Table for atomic configuration ---'
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
          !  plplus, qlplus for P;   ql for Pz
          write(stdo,"('conf: ', i4,i3,3x,i5,i3,6x,f8.3,1x,f8.3,' => ',10(i1,','))") &
               i,l,ipl,iplz,qlplus(l,i)+ql(1,l+1,i),qcc, (ix,ix=l+1,miplzipl-1)
       enddo
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
    if(ipr>=20) write(stdo,ftox)'conf: Species ',spid,'Z=',ftof(z,2),'Qc=',ftof(qc,3),'R=',ftof(rmt),'Q=',ftof(qtot),&
         'nsp=',nsp,'mom=',ftof(amgm)
    rmax = 50d0
    if (z < 10) rmax = 25
    if (z <=  6) rmax = 20
    call pshpr(0)
    call rmesh(z,rmt,lrel,lgrad,nrmx,a,nrmt) !Set up radial mesh for free atom ---
    call poppr
    b = rmt/(dexp(a*nrmt-a)-1d0)
    nr = 1d0+dlog(1d0+rmax/b)/a
    if (mod(nr,2) == 0) nr = nr-1
    rmax = b*(dexp(a*(nr-1))-1d0)
    if(iprint()>20)write(stdo,ftox)'conf: rmt rmax a=',ftof(rmt),ftof(rmax),ftof(a),'nrmt nr=',nrmt,nr
    ! --- Make atom self-consistent ---
    nitmax = nrmix
    nmix = -30
    ec(1) = 0
    print *,'goto atomc xxx'
    qelectron = z+qtot     !total number of electrons. 22mar2013
    call atomsc(.false.,n0,nsp,lmxa,z,0d0,kcor,lcor,qcor,rmax,a,nr, rofi,ec,ev,pl,ql,idmod,v,0d0,rhoin,rho,rhoc,nmix,qc,sec,stc, &
         sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,dq,exrmax,'gue', nitmax,lfrz,plplus,qlplus,nmcore,qelectron,vsum)
    print *,'end of atomsc xxxxx'
    print *,'vsum=',vsum,is
    if(ifives>0) then
       print *,'ifives=',ifives
       write(ifives,"(f23.15,i5,a)") vsum,is, ' !spatical integral of electrostatic potential'
       print *,' write end of ifives'
    endif
    if(ifiwv>0) write(ifiwv,"(f23.15,' ! total charge')") sum(qlplus(0:lmxa,1:nsp)+ql(1,1:lmxa+1,1:nsp))
    dq=dq+sum(qlplus(:,:)) !takao feb2011
    if (ipr>=20)write(stdo,ftox)'sumev=',ftof(sumev),'etot=',ftof(etot),'eref=',ftof(eref),'etot-eref=',ftof(etot-eref)
    !! june2012takao comment out followings. because confusing.
    !      call dcopy(lmxa+1,ql,3,qat,1)
    !      call dcopy(lmxa+1,ql(1,1,nsp),3,qat(1,2),1)
    !      call awrit4('fa  Pl %n:-1d  Ql %n:-1d',
    !     .' ',80,stdl,lmxa+1,pl,lmxa+1,qat)
    !      if (nsp .eq. 2)
    !     .call awrit4('fa  Pl2%n:-1d  Ql2%n:-1d',
    !     .' ',80,stdl,lmxa+1,pl(1,nsp),lmxa+1,qat(1,nsp))
    if (dabs(dq) > 1d-5 .AND. iprint() >= 10) write(stdo,ftox)' freeat (warning) atom not neutral, Q=',ftof(dq)
    ! .. Subtract core from density to make valence density
    do  isp = 1, nsp
       do  i = 1, nr
          rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)
          rho(i+(isp-1)*nr) = rho(i+(isp-1)*nr)-rhoc(i+(isp-1)*nr)
       enddo
    enddo
    ! --- Renormalize atom density or potential ---
    !irchan=0 !call ivset(irchan,1,n0,0)
    !! takao jun2012: rnatm is not tested ---> this is related to SPEC_ATOM_RCFA on.
    !call rnatm(pl,qat,n0,irchan,lmxa,z,a,b,rofi,ev,nr,rcfa,nsp,v,rho,plplus,qlplus)
    !     call prrmsh('starting total rho',rofi,rhot,nr,nr,nsp)
    do  isp = 1, nsp
       do  i = 1, nr
          rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)+rhoc(i+(isp-1)*nr)
       enddo
    enddo
    allocate(g_rv(nr*2))
    allocate(psi_rv(nr*(lmxa+1)*nsp))
    lplawv = 0
    if (ipr >= 50) lplawv = 1 !Print info about free-atom wavefunctions ---
    call pratfs (is, spid , lplawv , z , a , nr , rmax , nrmt , lmxa, pl,pz, nsp , v , rofi , g_rv , psi_rv, plplus,qlplus,ifiwv)
    deallocate(psi_rv,g_rv)
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
       if(ipr>=40)write(stdo,"(/' Charges:     inside',7x,'outside',7x,'sum'/' valence',3f13.6/' core   ',3f13.6/' total  ', &
            3f13.6)")                           qvin,       qvt-qvin,    qvt,   qcin,qct-qcin,qct,   qtin,qtt-qtin,  qtt
       write (stdl,"('fa Z',f6.1,'   rm',f7.2,'  qc',f6.2,'  qspl',f8.5,'  dq',f8.5,'  Etot',f15.6)")z,rmax,qc,qct-qcin,dq,etot
    else
       qvt = 0
    endif
    lplfa = 0
    if (qvt > 1d-6) then
       if (lplfa == 1) then
          write(stdo,"(/' write plot file with valence density..')")
          if (is < 10) write (str,'(''pl'',i1)') is
          if (is >= 10) write (str,'(''pl'',i2)') is
          open(newunit=ifi,file=trim(str)//'.'//trim(sname))
          write (ifi,"('# fit to fa density: ',a/'# rmt=',f7.3,'   rsm=',f7.3,'   nxi=',i2)") spid,rmt,rsmfa,nxi
          close(ifi)
       endif
       call tailsm(0,nr,nrmt,nsp,a,b,rmt,rsmfa,nxi0,nxi,exi,rofi, rho,rhot,hfc,hfct) !Attach smooth Hankel tails to valence density ---
    else
       call dpzero(hfc, nxi0*nsp)
       call dpzero(hfct, nxi0*nsp)
    endif
    call fctail(nr,nrmt,a,b,rfoca,rofi,rhoc,ccof,ceh) !Fit analytical expression to tail of core density ---
  end subroutine freats
  subroutine pratfs(is,spid,lplawv,z,a,nr,rmax,nrmt,lmaxa,pl,pz,nsp,v,rofi,g,psi,plplus,qlplus,ifiwv)!Prints out core and valence energy levels of free-atom
    use m_ext,only:sname
    use m_ftox
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
    implicit none
    integer :: n0,nr,lmaxa,nrmt,lplawv,nsp
    parameter (n0=10)
    real(8) :: z,a,rmax,pl(n0,nsp),v(nr,nsp),rofi(nr), pz(n0,nsp),g(2*nr),psi(nr,0:lmaxa,nsp)
    integer :: nglob,isp,l,konfig,nn,nre,i,konf,  konfg(0:8),ifi,lmaxc
    real(8) :: ev(0:20),pi,b,tol,eb1,eb2,dl,val,slo,sum,pzero, pmax,ctp,ecor,rmt
    character(1) :: lsym(0:n0-1), cc
    character:: spid*8
    character(15)::  str
    data lsym /'s','p','d','f','g','5','6','7','8','9'/
    integer::iz,ifiwv,ifipnu,is,ipz
    real(8):: pnu
    !! feb2011   "plplus,qlplus" mechanism is a fixing. Previous version did not allow SPEC_ATOM_Q setting for l with PZ.
    !! Now, the "plplus,qlplus" mechanism allows to set Q (valence charge, not including semicore charge).
    !! Our current version assumes MTOcore(specified by PZ) is below MTO(specified by P).
    real(8):: plplus(0:lmaxa,nsp),qlplus(0:lmaxa,nsp),pfree
    real(8):: sumr !mar2013
    character(8):: charext
    pi   = 4d0*datan(1d0)
    if (lmaxa > n0-1) call rx('pratfs:  lmax too large')
    b = rmax/(dexp(a*nr-a)-1d0)
    rmt = b*(dexp(a*nrmt-a)-1d0)
    tol = 1d-8
    write(stdo,"(/' Free-atom wavefunctions:')")
    if(ifiwv>0) write(ifiwv,"(i2,i3,' !nsp,lmaxa')")nsp,lmaxa
    open(newunit=ifipnu,file='atmpnu.'//trim(charext(is))//'.'//trim(sname))
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
             ipz=0
             if(konfig==mod(int(pz(l+1,isp)),10)) ipz=1 !ipz=1 for local orbital
             !write(6,*) 'konfig=',konfig,l,iz,'pz=',pz(l+1,isp)
             dl = dtan(pi*(0.5d0-pl(l+1,isp)))
             nn = konfig-l-1
             ev(l) = -0.5d0
             val = rmax
             slo = dl+1
             if (rmax > 9.99d0) then
                val = 1d-30
                slo = -val
             endif
             call rseq(eb1,eb2,ev(l),tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,rofi,nr,nre)
             call gintsl(g,g,a,b,nr,rofi,sum)
             call gintsl(g,g,a,b,nrmt,rofi,pmax)
             block
               real(8)::g1,g2,g3,g4,g5,drdi,valu,slou,du,dphi
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
               pnu  = nn+l+1  + (0.5d0-datan(dphi)/pi)
               pfree = nn+l+1 + 0.5d0 - datan(dble(l))/pi
               if(pnu<pfree) pnu=pfree
             endblock
             sum = sum - pmax
             sumr=pmax
             call ppratf(ev(l),z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
             cc = ' '
             if (dabs(ctp-rmax) < 1d-3) cc = '*'
             write(stdo,400) konfig,lsym(l),ev(l),pzero,pmax,ctp,cc,sum,pnu,ipz
             write(ifipnu,"(f23.16,i2,2x,i2,i2,i4,a1,f14.5,x,i2)") pnu,ipz,l,isp,konfig,lsym(l),ev(l)
             !   write default pnu setting to atmpnu file.
400          format(i4,a1,f14.5,2x,3f12.3,a,f12.6,f12.3,x,i2)
401          format(' valence:',6x,'eval',7x,'node at',6x,'max at',7x,'c.t.p.   rho(r>rmt)       pnu')
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
             call rseq(eb1,eb2,ecor,tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,rofi,nr,nre)
             call gintsl(g,g,a,b,nr,rofi,sum)
             call gintsl(g,g,a,b,nrmt,rofi,pmax)
             sum = sum - pmax
             call ppratf(ecor,z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
             write(stdo,400) konf,lsym(l),ecor,pzero,pmax,ctp,' ',sum 
403          format(/' core:        ecore',7x,'node at',6x,'max at',7x,'c.t.p.   rho(r>rmt)')
40        enddo
4011   enddo
80  enddo
    close(ifipnu)
    if(ifiwv>0) write(ifiwv,"('9 9 9 9 9 --------! this is terminator for an atom section')")
    ! --- Write file with valence wavefunctions
    if (lplawv == 1) then
       write (str,'(''wf_'',a)') spid
       write (stdo,344) str
344    format(/' Write valence wavefunctions to plot file: ',a)
       open(newunit=ifi,file=trim(str)//'.'//trim(sname))
       write (ifi,490) spid,rmax,rmt,nr,lmaxa,nr,1+nsp*(lmaxa+1)
490    format('# Free-atom wavefunctions (divided by r) for species ',a/'# rmax=',f7.3,'   rmt=',f7.3,'   nr=',i5,'   lmax=',i3/ &
            '% rows ',i5,' cols ',i3)
       do  50  i=1,nr
          write (ifi,495) rofi(i),((psi(i,l,isp),l=0,lmaxa),isp=1,nsp)
495       format(f9.5,1p,16d14.5)
50     enddo
       close(ifi)
    endif
!!! print out qbyl
    !      open('qinrmt.'//trim(sname))
    !      do  ib = 1, nbas
    !        ispec(ib)
    !        write(ifqbyl,"(i3,10f12.6)") sspec(ispec)%lmxa, (sum(qbyl_rv(il,1:nsp,ib)),il=1,lmaxa)
    !      enddo
    !      close(ifqbyl)
  end subroutine pratfs
  subroutine ppratf(e,z,nr,nre,rofi,a,b,v,g,pzero,pmax,ctp) ! Find outermost node and maximum of wavefct
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
    implicit none
    integer :: nr,nre
    real(8) :: a,b,ctp,e,pmax,pzero,z,rofi(nr),v(nr),g(nr)
    integer :: i,ir
    real(8) :: g1,g2,rho1,rho2,rho3,x
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
  subroutine fctail(nr,nrmt,a,b,rsm,rofi,rhoc,c,eh)!Fit one Hankel to tail of core density.
    use m_lmfinit,only: nsp
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
    implicit none
    integer :: nr,nrmt,i,ipr
    real(8) :: a,b,c,eh,rmt,rsm,rofi(nr),rhoc(nr,2),ak1,akap,fit,q,q0,r,s,v0,wt
    character sout*80
    call getpr(ipr)
    rmt = rofi(nrmt)
    q0 = 0d0
    do  i = nrmt, nr
       r = rofi(i)
       wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
       if (i == nrmt .OR. i == nr) wt = a*(r+b)/3d0
       q0 = q0 + wt*(rhoc(i,1)+rhoc(i,nsp))/(3-nsp)
    enddo
    v0 = (rhoc(nrmt,1)+rhoc(nrmt,nsp))/(3-nsp)/(rmt*rmt)
    write(stdo,"('conf: Core rhoc(rmt)=',f9.6,' spillout=',f9.6)")v0,q0
    if (dabs(q0) < 1d-6) then
       c = 0d0
       eh = -1d0
       return
    endif
    s = dsqrt(rmt**4 * v0**2 + 4*rmt*q0*v0) !Parameters of hankel fct
    ak1 = (rmt*rmt*v0+s)/(2d0*q0)
    akap = ak1
    c = rmt*v0*dexp(akap*rmt)
    eh = -akap*akap
    if(ipr >= 20) write(stdo,"(' Fit with Hankel e=',g0,' coeff=',g0)")eh,c
    if(ipr > 30) then
       write (stdo,"(6x,'r',12x,'rhoc',10x,'fit')")
       q = 0d0
       do  i = nrmt, nr
          r = rofi(i)
          wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
          if (i == nrmt .OR. i == nr) wt = a*(r+b)/3d0
          fit = c*dexp(-akap*r)*r
          q = q+wt*fit
          if ((rhoc(i,1)+rhoc(i,nsp))/(3-nsp) < 1d-8) exit
          if (mod(i,5) == 0 .OR. i == nrmt) write (stdo,"(f12.6,2f14.8)") r,(rhoc(i,1)+rhoc(i,nsp))/(3-nsp),fit
       enddo
    endif
  end subroutine fctail
  subroutine atomsc(  & !- Makes an atomic sphere self-consistent and get atomic charges
       lgdd,nl,nsp,lmax,z,rhozbk,kcor,lcor,qcor,rmax,a, nr,rofi,ec,ev,pnu,qnu,idmod,v,dv,rhoin,rho,rhoc,nrmix,qc,sumec, &
       sumtc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,qtot,exrmax, job,niter,lfrz,plplus,qlplus,nmcore,qelectron,vsum)
    use m_lmfinit,only: lrel
    use m_getqvc
    use m_lgunit,only: stdo
    use m_amix,only: amix
    ! ----------------------------------------------------------------
    !i Inputs
    !i   lgdd  :T  add q2 phi phidd into the density
    !i         :F  add q2 <phi phidd> phi phi into the density
    !i         :NB: both produce the same integrated density.
    !i         :Which one should be used depends on the context.
    !i   nl    :leading dimension of pnu,qnu.
    !i         :Also, total charge cutoff; see Remarks
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmax  :maximum l for this site (but see Remarks)
    !i   z     :nuclear charge
    !i   rhozbk:constant nuclear background density (jellium) added to z
    !i   kcor  :(partial core occupation) p.q.n for occupation
    !i   lcor  :(partial core occupation) l quantum for occupation
    !i   qcor  :(partial core occupation) core charge and moment
    !i   rmax  :potential, density calculated to rmax
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :                 -//-
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
    !i         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
    !i   qnu   :energy moments of charge (see Remarks)
    !i   idmod :0,1 or 2, specifing how the enu is set for an l-channel
    ! i   pz   :pnu for second panel (npan=2)
    ! i   qz   :qnu for second panel (npan=2)
    ! i   idmoz:idmod for second panel; specifies how semicore state
    ! i        :is treated
    ! i   npan :>1 for two-panel calculations
    !i   v     :spherical potential (job='pot', else input v not used)
    !i   dv    :constant potential shift added to sphere.
    !i         :Not used now.
    !i   rhoin:input density (job='rho', else input rhoin not used)
    !i         :used internally as a work array
    !i   nrmix :number of prior densities to mix to accelerate convergence
    !i         :to self-consistency (Anderson mixing); see Remarks.
    !i         :nrmix<0, linear mixing, with mixing beta |nrmix/100|
    !i   job   :(see Remarks)
    !i         :job='pot': start with potential
    !i         :job='rho': start with rhoin.
    !i         :job='gue': start with guessed charge density.
    !i   niter :number of iterations to attempt convergence to
    !i          self-consistency (see Remarks)
    !i   lfrz  :0 for soft core, >0 for frozen core
    ! o Input/Outputs
    ! o  ec    :core eigenvalues.  On input, these are guessed values.
    ! o        :if input ec(1) = 0, atomsc makes an internal choice fo ec.
    ! o  ev    :valence eigenvalues  On input, these are guessed values.
    ! o        :if input ec(1) = 0, atomsc makes an internal choice for ev.
    !o Outputs
    !o   rofi  :dimensioned (nr,2).
    !o         :(*,1) radial mesh points
    !o         :(*,2) weights
    !o   rho   :spherical charge density = 4 pi r**2 rhotrue
    !o   rhoc  :core charge density (unchanged if lfrz=1)
    !o   qc:   :core electronic charge
    !o   sumec :core single particle energy
    !o   sumtc :core kinetic energy (unchanged if lfrz=1)
    !o   ekin  :total kinetic energy
    !o   utot  :electrostatic energy
    !o   rhoeps:exchange-correlation energy
    !o   etot  :sphere total energy
    !o   amgm  :difference between spin up and spin down charge
    !o   rhrmx :true density at rmax
    !o   vrmax :true potential at rmax
    !o   qtot  :net charge in sphere
    !o   exrmax:exchange-correlation energy at rmax
    !r Remarks
    !r   Boundary conditions pnu, moments qnu, and the electrostatic
    !r   potential at rmax are all that is required to uniquely determine
    !r   a self-consistent spherical charge and potential.  'job'
    !r   determines how the starting potential is made, but the final
    !r   potential  and density should be independent of the initial choice.
    !r
    !r   atomsc uses the boundary condition that the potential at rmax
    !r   is zero, regardless of the total net charge inside the sphere.
    !r   See subroutine madpot for discussion of how this choice affects
    !r   the Madelung energy.
    !r
    !r   Sphere total energy is sum of K.E., Hartree energy, XC energy:
    !r      etot = ekin + utot + rhoeps
    !r   The kinetic energy is computed via double-counting terms
    !r     ekin = sumev + sumec + dsumec - rhovh - rhomu
    !b Bugs
    !b   Total energy terms need to be cleaned up and simplified.
    ! ----------------------------------------------------------------
    use m_ftox
    use m_vxcatom,only: vxc0sp
    implicit none
    logical :: lfrz,lgdd
    character job*3
    integer :: nr,nsp,nl,nrmix,niter,kcor,lcor,ncmx,nvmx,lmax, idmod(0:nl-1)
    parameter (ncmx=50, nvmx=20)
    double precision :: ec(ncmx),ev(nvmx),rofi(nr,2),v(nr,nsp),rho(nr,nsp),rhoc(nr,nsp),rhoin(nr,nsp),pnu(nl,2), &
         qnu(3,nl,nsp),z,rmax,a,qc,vrmax(nsp),exrmax(2),dv,rhrmx, rhozbk,qcor(2),amgm,ekin,etot,qtot,rhoeps,sumec,sumev,sumtc,utot
    logical :: last,cmdopt,ltmp
    integer :: iprint,ir,isp,iter,jmix,nglob,ipr1,l,ii
    character strn*10
    real(8):: b,ddot,decay,dl,dold,drdi,drho,dsumec,ea,fac,pi,rho0t,rhomu,rhovh,rmsdel,ro,ssum,tolrsq,vhrmax,vnucl,vrhoc,vsum, &
         zvnucl,rvh(2),rho0(2),reps(2),rmu(2),sec(2),stc(2),sev(2),dasum
    integer:: nmix 
    real(8) ,allocatable :: rho_rv(:)
    real(8):: plplus(0:lmax,nsp),qlplus(0:lmax,nsp),qelectron
    double precision :: norm(10,10),awk(10,2),beta,beta1
    parameter (beta = 0.3d0)
    double precision :: tolch,tl
    parameter (tolch=5d-5,tolrsq=1d-12)! tolch is tolerance for change in the charge density, tolv for rseq
    integer:: nmcore
    print *,'atomsc nmcore=',nmcore
    if (lmax >= nl) call rx('atomsc:  lmax too large')
    pi = 4d0*datan(1d0)
    b = rmax/(dexp(a*(nr-1)) - 1)
    vhrmax = 0d0
    nmix = min(max(nrmix,0),10)
    sec = 0d0
    stc = 0d0
    allocate(rho_rv(nr*nsp*2*(nmix+2)))
    ! --- Core charge, radial mesh points and weights ---
    if (kcor /= 0) then
       if (qcor(1) /= 0 .OR. qcor(2) /= 0) then
          write(stdo,ftox)'Add core hole:  kcor=',kcor,' lcor=',lcor,' qcor=',ftof(qcor),'amom=',ftof(qcor(2))
       endif
    endif
    call getqvc(nsp,nl,lmax,z,pnu,qnu,0,0,kcor,lcor,qcor,qc,qtot,amgm)
    ! --- Guesses for core and valence eigenvalues ---   ! takao probably not needed here
    if (ec(1) == 0) call getqvc(nsp,nl,lmax,z,pnu,qnu,ncmx,nvmx,kcor,lcor,qcor,qc,qtot,amgm,ec,ev) !only for core number check since we set ev and ev below.
    ! ill be replace by simple code.
    call radmsh(rmax,a,nr,rofi)
    call radwgt(rmax,a,nr,rofi(1,2))
    !! initialize ev here(overide getqvc now).feb2011
    ev = -0.5d0
    ec = -5.0d0
    if (job == 'pot') then !Initial charge density ---
       call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rhoin,rhoc,kcor, &
            lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tolrsq,nsp,lfrz,000,plplus,qlplus,nmcore)
       if (niter == 0) then
          rho=rhoin
          if(allocated(rho_rv)) deallocate(rho_rv)
          return
       endif
    else if (job == 'gue') then
       decay = 1d0+z/10d0
       decay = dmin1(decay,5d0)
       decay = 5
       rhoin(:,1) = dexp(-decay*rofi(:,1))*rofi(:,1)**2
       fac = z/(sum(rhoin(:,1)*a*(rofi(:,1)+b))*nsp)
       rhoin(:,1)=rhoin(:,1)*fac
       if(nsp==2) rhoin(:,2)=rhoin(:,1)
    elseif (job /= 'rho') then
       call rx('atomsc: job not pot|rho|gue')
    endif
    rho_rv(1+nr*nsp*(nmix+2):nr*nsp+nr*nsp*(nmix+2))= reshape(rhoin(1:nr,1:nsp),[nr*nsp])
    drho = 100d0
    last = .false.
    if (iprint() >= 41) write(stdo,341)
    jmix = 0
    dold = 1
    beta1 = beta
    vrhoc=-1d99
    if (nrmix < 0) beta1 = dble(-nrmix)/100
    do  35  iter = 1, niter !Start self-consistency loop ---
       tl = tolrsq
       call addzbk(rofi,nr,nsp,rhoin,rhozbk,-1d0)
       if(abs(rmax-rofi(nr,1))>1d-6) call rx('atomsr.F:something wrong. abs(rmax-rofi(nr,1))>1d-3')
       vhrmax=2d0*(qelectron-z)/rmax
       call poiss0(z,a,b,rofi,rhoin,nr,vhrmax,v,rvh,vsum,nsp) !  Hartree potential
       vsum = vsum + 4d0*pi*(z-qelectron)*rmax**2
       call addzbk(rofi,nr,nsp,rhoin,rhozbk,1d0)
       vnucl = v(1,1)  
       if (last .AND. iprint() >= 50) call pshpr(80)
       call vxc0sp(a,b,rofi,rhoin,nr,v,rho0,reps,rmu,nsp,exrmax(1))! Exchange-correlation potential
       if (last .AND. iprint() >= 50) call poppr
       fac = 4*pi*rofi(nr,1)**2
       rhrmx = rhoin(nr,1)/fac
       if (nsp == 2) then !       Get rhrmx, exrmax
          exrmax(2) = exrmax(1)
          rhrmx = rhrmx + rhoin(nr,2)/fac
       endif
       ipr1 = 0
       if (last .AND. iprint()>= 40) ipr1 = 1
       if (last .AND. iprint() > 40) ipr1 = 2
       if( .NOT. last) call pshpr(15) !low print()
       call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rho,rhoc, &
            kcor,lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tl,nsp,lfrz,ipr1,plplus,qlplus,nmcore)
       if( .NOT. last) call poppr !set back to original print()
       drho = 0d0
       ssum = 0d0
       vrhoc = 0d0
       rho0t = 0d0
       do  40  isp = 1, nsp
          rho0t = rho0t + rho0(isp)
          ssum = ssum + ddot(nr,rofi(1,2),1,rho(1,isp),1)
          do  42  ir = 1, nr
             drdi = a*(rofi(ir,1) + b)
             drho = drho + rofi(ir,2)/drdi*dabs(rho(ir,isp)-rhoin(ir,isp))
42        enddo
          do  41  ir = 2, nr
             vrhoc = vrhoc + rofi(ir,2)*(v(ir,isp)-2*z/rofi(ir,1))*rhoc(ir,isp)
41        enddo
40     enddo
       rho_rv(1:nr*nsp) = reshape(rho(1:nr,1:nsp),[nr*nsp])
       jmix = amix(nr*nsp , min ( jmix , nmix ) , nmix , 0 , beta1, iprint ( ) - 70 , .9d0 ,  rho_rv , awk , rmsdel )
       rhoin(1:nr,1:nsp) = reshape(rho_rv(1+nr*nsp*(nmix+2):nr*nsp+nr*nsp*(nmix+2)),[nr,nsp])
       if (last) goto 90
       if (iprint() >= 41 .AND.(drho < tolch .OR. iter == niter-1 .OR. iter == 1)) &
            write(stdo,340) iter,ssum,drho,vnucl,rho0t,vsum,beta1
340    format(i5,f12.6,1p,e12.3,0p,f14.4,e14.4,f14.4,f7.2)
341    format(/'  iter     qint',9x,'drho',10x,'vh0',10x,'rho0',10x,'vsum',5x,'beta')
       last = (drho .lt. tolch .or. iter .eq. niter-1)
       jmix = jmix+1
       beta1 = min(max((1-drho/dold)/beta1,beta1-.2d0,beta), 1d0,beta1+.2d0)!Beta for next iteration
       if (nmix > 0 .AND. drho < 1) beta1 = 1
       if (nrmix < 0) beta1 = dble(-nrmix)/100
       dold = drho
35  enddo
90  continue
    if (allocated(rho_rv)) deallocate(rho_rv)
    if (iprint() >= 30) write(stdo,'(1x)')
    ! --- Collect terms for total energy ---
    rhoeps = 0d0
    rhomu  = 0d0
    sumev  = 0d0
    if ( .NOT. lfrz) then
       sumec = 0d0
       sumtc = 0d0
    endif
    rhovh  = 0d0
    do  isp = 1, nsp
       if(nsp==2.AND.iprint()>30)write(stdo,"(' Spin',i2,':',/' vrmax=',f12.5,'    sumev= ',f12.5,'    sumec=',f12.5,/' rhovh=',&
            f12.5,'    rhoeps=',f12.5,'    rhomu=',f12.5)") isp,v(nr,isp)-2*z/rmax, sev(isp),sec(isp),rvh(isp),reps(isp),rmu(isp)
       rhoeps = rhoeps + reps(isp)
       rhomu = rhomu + rmu(isp)
       sumev = sumev + sev(isp)
       if ( .NOT. lfrz) then
          sumec = sumec + sec(isp)
          sumtc = sumtc + stc(isp)
       endif
       rhovh = rhovh + rvh(isp)
    enddo
    zvnucl = -z*vnucl
    utot = .5d0*(rhovh + zvnucl)
    dsumec = vrhoc - (sumec-sumtc)! Correction to core eigenvalues if sumec not obtained from this V
    ekin = sumev + sumec + dsumec - rhovh - rhomu
    etot = ekin + utot + rhoeps
    if (iprint() >= 40) write(stdo,139) sumev,sumec,vnucl,rhovh,zvnucl,utot,rhomu,rhoeps,dsumec,ekin,sumtc,etot
139 format(/' sumev=',f13.6,'    sumec =',f13.6,'   vnucl =',f13.6 /' rhovh=',f13.6,'    zvnucl=',f13.6,'   utot  =',f13.6 &
         /' rhomu=',f13.6,'    rhoeps=',f13.6,'   dsumec=',f13.6     /' ekin= ',f13.6,'    tcore =',f13.6,'   etot  =',f13.6)
    vrmax(1) = -2*z/rmax + sum(v(nr,1:nsp))/nsp
    vrmax(2) = 0d0
    if(nsp==2) vrmax(2) = v(nr,1)-v(nr,2)
  end subroutine atomsc
  subroutine addzbk(rofi,nr,nsp,rho,rhozbk,scale)
    implicit none
    integer :: nr,nsp,ir,isp
    double precision :: rofi(nr),rho(nr,*),rhozbk,scale,s
    if(rhozbk == 0) return
    s = 16d0*datan(1d0)*scale*rhozbk
    do isp = 1, nsp
       rho(:,isp) = rho(:,isp) + s*rofi(:)**2
    enddo
  end subroutine addzbk
  subroutine newrho(z,lrel,lgdd,nl,nlr,lmax,a,b,nr,rofi,v,rho,rhoc, &
       kcor,lcor,qcor,pnu,qnu,sumec,sumtc,sumev,ec,ev,tol,nsp,lfrz,ipr,plplus,qlplus,nmcore)
    use m_lmfinit,only: c=>cc
    !! ev is dummy now. feb2010 takao
    !- Makes spherical charge density for a spherical potential.
    !  ---------------------------------------------------
    !i Inputs:
    !i   z     :nuclear charge
    !i   lrel  :0 for nonrelativistic, 1 for relativistic
    !i   lgdd  :T q2 is coefficient to phidot**2 + phi*phidotdot
    !i         :F q2 is coefficient to phidot**2 - p phi**2
    !i         :Both produce the same integrated density; see Remarks.
    !i         :lgdd=F follows Stuttgart conventions.
    !i   nl    :(global maximum l) + 1
    !i   nlr   :second dimension of rho:
    !r         :1 if spherical part of rho is to be generated
    !i         :nl if generated rho is to be decomposed by l
    !i         :In the latter case, the core is not calculated
    !i   lmax  :maximum l for a given site
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   v     :spherical potential (atomsr.f)
    !i   kcor  :(partial core occupation) p.q.n for occupation
    !i   lcor  :(partial core occupation) l quantum for occupation
    !i   qcor  :(partial core occupation) core charge and moment
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
    !i         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
    !i   qnu   :energy moments of charge (see Remarks)
    !i   tol   :precision to which wave functions are integrated.
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   g     :normalized wave function times r
    !i   gp    :energy derivatives of g
    !i   lfrz  :T, do not make core rho
    !i   ipr   :0 no printout
    !i         :1 summary printout of core
    !i         :2 detailed printout of core
    !o Outputs:
    !o   rho   :spherical charge density times 4*pi*r*r
    !o   rhoc  :core density times 4*pi*r*r
    !o   sumec :sum of core eigenvalues
    !o   sumtc :core kinetic energy = sumec - v*rhoc
    !o   sumev :sum of valence eigenvalues
    !o   ec    :core eigenvalues
    !o   ev    :valence eigenvalues !not output now
    !r Remarks:
    !r   rho is determined by boundary conditions pnu and moments qnu.
    !r   For rmax>10 sets phidot, phidotdot are made zero.
    !r
    !r   Recall that val,slo correspond to u = r*phi, so
    !r     rmax * slo / val = D_u = D_phi + 1.
    !r     For val=rmax  slo = D + 1
    !r
    !r   Switch lgdd concerns the contribution of <phidot|phidot> to the
    !r   sphere charge.  The zeroth moment may be defined as the amount of
    !r   charge in channel l, i.e.
    !r     q^0 = (amount of <phi|phi>) +  p* (amount of <phidot|phidot>)
    !r   where
    !r     p = <phidot|phidot>
    !r   Then q^0 is not the amount of phi*phi inside the sphere, but rather
    !r   but instead, it is
    !r     (q_0 - p q_2)
    !r   since q2 is the amount of <phidot|phidot>.  The charge density is
    !r     rho= (q_0 - p q_2) phi*phi + 2 q_1 phi*phidot + q_2 phidot*phidot
    !r   This is the Stuttgart convention (lgdd=F)
    !r
    !r   Methfessel convention:  to avoid explicit dependence of rho on
    !r   p, he approximated p*phi*phi with -phi*phidotdot (they have the
    !r   same integrated charge).  Then
    !r     rho= q_0 phi*phi + 2 q_1 phi*phidot +
    !r          q_2 (phidot*phidot + phi*phidotdot)
    !  ---------------------------------------------------
    use m_rhocor,only:rhocor
    use m_phidx,only: phidx
    use m_lmfinit,only:cc
    implicit none
    logical :: lgdd,lfrz
    integer :: nl,nlr,lmax,nr,nsp,lrel,kcor,lcor,ipr,iz
    double precision :: z,a,b,tol,sumev(nsp),sumec(nsp),sumtc(nsp),ec(*),ev(*),qcor(2),v(nr,nsp),rofi(nr),&
         rho(nr,nlr,nsp),rhoc(nr,nsp), qnu(3,nl,nsp),pnu(nl,nsp)
    logical :: free
    integer :: konfig(0:10),l,isp,ir,ival,nn,nre,jr,lmaxc,lr,k,nrmx
    parameter (nrmx=1501)
    real(8):: rocrit,eb1,eb2,q0,q1,q2,rmax,eval,dl,val(5),slo(5),sum,ro,phi,dphi,phip,dphip,p,fllp1,r,tmc,gfac,q00,&
         g(2*nrmx),gp(2*nrmx*4)
    real(8):: plplus(0:lmax,nsp),qlplus(0:lmax,nsp)
    real(8),parameter:: pi = 4d0*datan(1d0)
    integer::nmcore !jun2012
    if (nr > nrmx) call rxi(' newrho: increase nrx, need',nr)
    lr = 1
    rocrit = 0.002d0/4
    eb1 = -50d0
    eb2 =  50d0
    rmax = rofi(nr)
    free = rmax>9.99d0
    call config(pnu,lmax,z,konfig,lmaxc)
    if (kcor > 0) then
       lmaxc = max(lmaxc,lcor)
       konfig(lcor) = max(konfig(lcor),kcor+1)
    endif
    ! --- Calculate core density ---
    if (nlr == 1) then
       if ( .NOT. lfrz) then
          rhoc=0d0
          call rhocor(0,z,lmaxc,nsp,konfig,a,b,nr,rofi,v,g,kcor,lcor,qcor,tol,ec,sumec,sumtc,rhoc,nmcore=nmcore,ipr=ipr)
       endif
       call dcopy(nr*nsp,rhoc,1,rho,1)
    endif
    eval=-0.5d0 !initial condition. the same as getqvc ! --- Loop over valence states ---
    sumev= 0d0
    do  202  isp = 1, nsp
       do  201  l = 0, lmax
          izloop: do  20  iz=0,1  !takao feb2011
             if (nlr == nl) lr = l+1
             q0 = max(qnu(1,l+1,isp),0d0)
             q1 = qnu(2,l+1,isp)
             q2 = qnu(3,l+1,isp)
             if (q0 < 1d-6) goto 20
             nn = int(pnu(l+1,isp)) - l - 1
             val(1) = rmax
             dl = dtan(pi*(0.5d0 - pnu(l+1,isp)))
             slo(1) = dl + 1
             !! Setting when iz=1 override setting when iz=1
             if(iz==1)then
                if(plplus(l,isp)>0d0 .AND. qlplus(l,isp)>0d0) then ! note this is a case of Qv(lower principle quantum number).
                   !  Search "=== Charge for l" in freeat.F, and NOTE: above it.   print *,'vvvvvv', plplus(l,isp)
                   q0=qlplus(l,isp)
                   q1=0d0
                   q2=0d0
                   nn=int(plplus(l,isp)) - l - 1
                   dl = dtan(pi*(0.5d0 - plplus(l,isp)))
                   slo(1) = dl + 1
                else
                   cycle
                endif
             endif
             if (free) val(1) = 1d-30
             if (free) slo(1) = -val(1)
             call rseq(eb1,eb2,eval,tol,z,l,nn,val(1),slo(1),v(1,isp),g, sum,a,b,rofi,nr,nre) 
             sumev(isp) = sumev(isp) + eval*q0 + q1
             ro = g(nr)**2
             if(.NOT.free.AND.ro<rocrit) write(*,"(' NEWRHO (warning): PHP,PHPP set to zero,l,nn,nre,rho=',3i5,2f8.4)") l,nn,nre,ro
             if (free .OR. ro < rocrit) then
                gp=0d0 !call dpzero(gp,8*nr)
                p = 0
             else
                val(1) = val(1)/dsqrt(sum)
                slo(1) = slo(1)/dsqrt(sum)
                call phidx(1,z,l,v(1,isp),0d0,0d0,rofi,nr,2,tol,eval,val,slo,nn,g,gp,phi,dphi,phip,dphip,p,0d0,[0d0],0d0,[0d0])
             endif
             fllp1 = l*(l+1)
             !  ...  Case add q2 phi phidd rho
             if (lgdd) then
                k = 2*nr
                do  ir = 2, nre
                   jr = ir + nr
                   r = rofi(ir)
                   tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
                   gfac = 1d0 + fllp1/(tmc*r)**2
                   rho(ir,lr,isp) =  rho(ir,lr,isp) +  q0*(gfac*g(ir)**2 + g(jr)**2) + &
                        2*q1*(gfac*g(ir)*gp(ir)+g(jr)*gp(jr)) + q2*(gfac*(gp(ir)**2 + g(ir)*gp(ir+k)) + gp(jr)**2 + g(jr)*gp(jr+k))
                enddo
             else!Case add -p q2 phi phi into rho
                q00 = q0-p*q2
                do ir = 2, nre
                   jr = ir + nr
                   r = rofi(ir)
                   tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
                   gfac = 1d0 + fllp1/(tmc*r)**2
                   rho(ir,lr,isp) = rho(ir,lr,isp) &
                        + q00*(gfac*g(ir)**2+g(jr)**2) + 2*q1*(gfac*g(ir)*gp(ir)+g(jr)*gp(jr)) + q2*(gfac*gp(ir)**2+gp(jr)**2)
                enddo
             endif
20        enddo izloop
201    enddo
202 enddo
    !      call tcx('newrho')
  end subroutine newrho
  subroutine gintsl(g1,g2,a,b,nr,rofi, ssum) ! Integrate inner product of two wave equations, large component only
    implicit none ! g1,g2  : wave function, mesh rofi(i) = b [e^(a(i-1)) -1],i=1,nr
    integer :: nr,ir
    real(8) :: a,b,g1(nr),g2(nr),rofi(nr),ssum
    ssum = (4d0*sum([((rofi(ir)+b)*(g1(ir)*g2(ir)),ir=2,nr-1,2)]) & !2,4,... nr-1
         +  2d0*sum([((rofi(ir)+b)*(g1(ir)*g2(ir)),ir=3,nr-2,2)]) & !3,5,7,...
         +            (rofi(nr)+b)*(g1(nr)*g2(nr)) )*a/3d0
  end subroutine gintsl
subroutine tailsm(lrhot,nr,nrmt,nsp,a,b,rmt,rsm,nxi0,nxi,exi,rofi, rho,rhot,hfc,hfct)  !- Fit tails of rho to smoothed Hankel functions
  use m_lgunit,only:stdo
  use m_hansr,only :hansr,hansmr
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lrhot :0 make fit for rho only; 1 make fit for rho and rhot
  !i   nr    :number of radial mesh points
  !i   nrmt  :number of points between 0..rmt
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   rmt   :muffin-tin radius, in a.u.
  !i   rsm   :smoothing radius for smoothed hankel fit
  !i   nxi0  :leading dimension of hfc, nfct
  !i   nxi   :number of energies to include in fit
  !i   exi   :energies to include in fit
  !i   rofi  :radial mesh points
  !i   rho   :spherical valence charge density times 4*pi*r*r
  !i   rhot  :total charge density (used if lrhot=1)
  !o Outputs
  !o   hfc   :coefficients to h.f. fits to rho
  !o   hfct  :coefficients to h.f. fits to rhot (if lrhot=1)
  !l Local variables
  !l   rsq   :work array holding rofi**2
  !r Remarks
  !r   A fit is make to tails of the valence charge density for r>rmt
  !r   using smoothed hankel functions of smoothing radius rsm.
  !r   Fit is constrained to agree with integrated charge.
  !u Updates
  !u   19 Apr 02 Move rsq to a local array.  Altered argument list.
  ! ----------------------------------------------------------------------
  implicit none
  integer :: lrhot,nsp,nr,nrmt,nxi0,nxi
  double precision :: a,b,rmt,rsm
  double precision :: hfc(nxi0,nsp),exi(nxi),rho(*),rofi(nr),hfct(nxi0,nsp),rhot(*)
  logical :: lzero
  integer:: lx0(20) , isp , ie , ir , ipr
  integer,allocatable :: idx_rv(:)
  real(8) ,allocatable :: xi_rv(:)
  double precision :: sr,x0(0:2),xi(0:2),fpi,qout,qcst(20),qcst0(20), err,rsq(nr)
  print *,'tailsm: init'
  fpi = 16*datan(1d0)
  call getpr(ipr)
  ! --- Tabulate smoothed Hankels on a radial mesh ---
  rsq = rofi**2
  lx0=0d0
  allocate(xi_rv(nr*nxi))
  allocate(idx_rv(nr*2))
  !     Patch to properly handle case when rsm=0
  lzero = rsq(1) .eq. 0 .and. rsm .eq. 0
  if (lzero) rsq(1) = rsq(2)/100
  call hansr ( rsm , 0 , 0 , nxi , lx0 , exi , rsq , nr , nr , idx_rv, 0 , xi_rv )
  if (lzero) rsq(1) = 0
  deallocate(idx_rv)
  ! --- Fit smoothed hankels to rho ---
  isploop: do  20  isp = 1, nsp
     if (isp == 1 .AND. ipr >= 30) write(stdo,333) nxi,rmt,rsm
333  format(/' tailsm: fit tails to',i2,' smoothed hankels, rmt=', f8.5,', rsm=',f8.5)
     if (isp == 2 .AND. ipr >= 20) write(stdo,'(/'' tailsm: spin 2 ...'')')
     !   ... Fitting constraints for smoothed Hankels
     do  ie = 1, nxi
        if (rsm < 1d-9) then
           sr = dsqrt(-exi(ie))*rmt
           qcst(ie) = -dsqrt(fpi)*(sr+1)*dexp(-sr)/exi(ie)
        else
           call hansmr(rmt,0d0,1/rsm,x0,1)
           call hansmr(rmt,exi(ie),1/rsm,xi,1)
           qcst(ie) = dsqrt(fpi)/exi(ie)*(-dexp(rsm**2/4*exi(ie)) &
                - rmt**3*(xi(1)-dexp(rsm**2/4*exi(ie))*x0(1)))
        endif
        qcst0(ie) = qcst(ie)
     enddo
     !   ... Fit for this spin
     call hnsmft ( rofi , rho ( 1 + ( isp - 1 ) * nr ) , nr , qout &
          , a , b , nrmt , exi , qcst , xi_rv , hfc ( 1 , isp ) , nxi , err )
     if (isp==1 .AND. ipr>=20)  write(stdo,"(' tailsm:  fit tails to ',i8,' functions with' &
          //' rsm=',d13.5,' rms error=',d13.5)") nxi,rsm,err
     !   ... Fit a second time for the full density in rhot
     if (lrhot /= 0) then
        do  ie=1,nxi
           qcst(ie)=qcst0(ie)
        enddo
        call hnsmft ( rofi , rhot ( 1 + ( isp - 1 ) * nr ) , nr , qout &
             , a , b , nrmt , exi , qcst , xi_rv , hfct ( 1 , isp ) , nxi , err )
     endif
20 enddo isploop
  if (allocated(xi_rv)) deallocate(xi_rv)
end subroutine tailsm
end module m_freeat
