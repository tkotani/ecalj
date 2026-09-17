module m_heftet
  use m_ftox
  use m_lgunit,only: stdo
  public heftet
contains
  subroutine heftet() bind(C)! Calculates the Fermi energy by tetrahedron method.
    use m_GWinput, only: gwinput_init, gwinput_loaded, tg_tetrakbt => tetrakbt
    use m_read_bzdata,only: read_bzdata, idteti,qbz,qibz,dq_,nqibz,ntetf,nteti,ginv,nqbz
    use m_genallcf_v3,only: genallcf_v3,natom,nspin,nl,nn,nnv,nnc,nlnmx, nctot,niw,nspx
    use m_genallcf_v3,only: alat, deltaw,esmr, plat, pos,z,ecore, konf,nlnx,valn=>qval
    use m_hamindex,only:   Readhamindex,qtt,nqtt
    use m_readeigen,only: init_readeigen,readeval
    use m_tetrakbt,only: tetrakbt_init, kbt
    use m_keyvalue,only: getkeyvalue
    use m_genallcf_v3,only: nprecb,mrecb,mrece,nqbzt,nband,mrecg !  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg  
    use m_mpi,only: MPI__Initialize
    use m_lgunit,only: m_lgunit_init
    use m_bzints,only: bzints2x
    use m_cmdopt_registry, only: c2_job
    implicit none
    integer :: mxclass,ngnmax, ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd, &
         nprecx,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx, ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq, &
         i,k, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,iimdim, ifreqx,ifreqw,iwx,iexpa,mxkp,nqibzxx,ntet,nene,iqi, ix,iw, &
         nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx, nq2,ntq2,iqx,itx
    integer:: nkp,itt,ifief ,imode,iftote,it, nptdos_kbt = 100001, ifief_kbt, nbandx, nptdos=101,  nwin,bzcase=1,ifi
    integer,allocatable ::idtetx(:,:),idtet(:,:),ipq(:),iene(:,:,:),ibzx(:)   ! fermi_kbt
    real(8) ::tpia,vol,voltot,rs,alpha, qfermi,efx,efnew,edummy,efz,qm,xsex,egex, dscdw1,dscdw2,dscdw,zfac,efxx2,tripl,&
         qqex(1:3), eex,exsp,eee, exwgt,qq(3),        elda, eqp01, eqp02, efermi_kbt, e11, e22
    real(8):: elo,ehi,e1,e2,efermi,dum,dosef,efxx,rydberg,dum111(1,1,1)  , tol=1d-8,toql=1d-8, volwgt, ddq(3),bandgap,tolq=1d-8
    real(8),allocatable:: dos(:), qz(:,:), eband(:,:,:),eband2(:,:,:), ene(:), dos_kbt(:), dosef_kbt
    real(8),parameter:: pi= 4d0*datan(1d0)
    logical :: metal,qbzreg,usetetrakbt
    call MPI__Initialize()
    call M_lgunit_init()
    if (c2_job >= 0) then
      imode = c2_job
    else
      write(stdo,*) 'mode=(1-5)?'
      read(5,*) imode
    endif
    write(stdo,*) '--- heftet: calculation mode =',imode
    if(imode<1 .OR. imode>5) call rx( 'mode out of range(1-4)')
    call read_BZDATA() ;          write(stdo,ftox)' heftet: nqibz, ntetf,nteti,nqbz ', nqibz,ntetf,nteti,nqbz
    call genallcf_v3(incwfx=-1) !; if(nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
    tpia   = 2d0*pi/alat
    voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3))) !cell volume
    ! ef for rs for the empty-sphere test case 
    !$$$      legas = .false.; INQUIRE (FILE = 'LEGAS', EXIST = legas)
    !$$$      if(legas) then !!! test for electron gas case.
    !$$$        write(stdo,*)' find LEGAS. legas =',legas
    !$$$        iflegas = 2101
    !$$$        open (iflegas,file='LEGAS')
    !$$$        read(iflegas,*)rs !rs parameter for homogas
    !$$$        close(iflegas)
    !$$$        alpha = (9*pi/4d0)**(1d0/3d0)
    !$$$        qfermi = alpha/rs
    !$$$        efx  = qfermi**2
    !$$$        valn = efx**1.5d0*voltot/3d0/pi**2
    !$$$        write (stdo,*)'  #### egas test mode  legas=T #### given rs =',rs
    !$$$        write (stdo,*)' egas  Exact Fermi momentum  qf  =', qfermi
    !$$$        write (stdo,*)' egas  Exact Fermi energy    Ef  =', efx
    !$$$      endif
    write(stdo, *) ' --- computational conditions --- '
    write(stdo,'("    deltaw  =",f13.6)') deltaw
    write(stdo,'("    esmr    =",f13.6)') esmr
    write(stdo,'("    alat voltot =",2f13.6)') alat, voltot
    if(imode==1 .OR. imode==5) then
      call Readhamindex()
      call init_readeigen() !initialization of readEigen
      nbandx= nband
      nkp   = nqibz
      allocate( eband2(nbandx,nspx,nqbz) )
      ddq = 0d0
      if( .NOT. qbzreg()) ddq= -dq_
      write(stdo,"(' nspx nqbz, q is replaced by q+dq_; dq_ =',2i5,3f13.5)") nspx,nqbz,dq_
      do is  = 1,nspx !Readin eband
        ix = 0
        do iqi = 1,nqbz
          ix = ix+1
          eband2(:,is,ix)= readeval(qbz(:,iqi)+ddq, is)
        enddo
      enddo
      call efrang3(nspx,nqbz,nbandx,valn,eband2,  e1,e2,elo,ehi,bandgap)
      if(qbzreg() .AND. bzcase==1 ) then
        allocate( eband(nbandx,nspx,nqibz))
        do is  = 1,nspx
          do iqi = 1,nqibz
            eband(:,is,iqi)=readeval(qibz(:,iqi),is)
            write(stdo,"(' -------- qibz=',2i5,3f9.4)")iqi,is,qibz(:,iqi)
          enddo
        enddo
      endif
    elseif(imode==2 .OR. imode==3 .OR. imode==4) then
      open(newunit=iftote,file='TOTE.UP')
      read(iftote,*) nq,ntq, efxx
      if(nq/= nqibz) call rx( ' heftet: nq TOTE/= nqibz')
      nbandx = ntq
      allocate( eband(nbandx,nspx,nqibz), qz(3,nqibz) )
      qz=qibz 
      do  is  = 1,nspx !Readin eband
        if(is==2) then
          open(newunit=iftote,file='TOTE.DN')
          read(iftote,*) nq2,ntq2, efxx2
          if(nq2/=nq)   call rx( 'nq dirrerent TOTE.UP TOTE.DN')
          if(ntq2/=ntq) call rx( 'ntq dirrerent TOTE.UP TOTE.DN')
          if(efxx/=efxx)call rx( 'efxx dirrerent TOTE.UP TOTE.DN')
        endif
        do  iqi = 1,nqibz
          do   it = 1,nbandx
            read(iftote,*) qq(1:3), itx, iqx, elda, eqp01, eqp02, zfac
            if(it /= itx ) call rx( ' heftet: it /=itx ')
            if(sum(abs(qq(1:3)-qz(1:3,iqi)))>tolq ) then
              write(stdo,*) 'iqi it=',iqi,it
              write(stdo,*) qq(1:3)
              write(stdo,*) qz(1:3,iqi)
              call rx( ' heftet: q in TOTE /=qz (iqi) ')
            endif
            if(imode==2) eband(it,is,iqi) = elda /rydberg() + efxx
            if(imode==3) eband(it,is,iqi) = eqp01/rydberg() + efxx
            if(imode==4) eband(it,is,iqi) = eqp02/rydberg() + efxx
          enddo
        enddo
        close(iftote)
      enddo
      nkp =nqibz
      call efrang3(nspx,nkp,nbandx,valn,eband,  e1,e2,elo,ehi,bandgap)
    endif
    write(stdo,"(' end of efrang3: e1 e2 nteti=',2d15.7,i6)") e1,e2,nteti
    write(stdo,"(' elo ehi = ',2d15.7,i6)") elo,ehi
    if(e1/=e2.and. (.NOT.(qbzreg().AND.bzcase==1)) ) &
         call rx("Fermi energy finder for off-regular mesh (Chi_RegQbz off)or bzcase==2 is not implimented yet...")
    if(imode==1) open(newunit=ifief,file='EFERMI')
    if(imode==5) open(newunit=ifief,file='EFERMI')
    if(imode==2) open(newunit=ifief,file='EFERMI.check')
    if(imode==3) open(newunit=ifief,file='EFERMI.QP')
    if(imode==4) open(newunit=ifief,file='EFERMI.QPz=1')
    if(e1 == e2) then ! --- Find band gap  e1=e2 is at the middle of gap---
      efermi     = e1
      write(stdo,*)'heftet: only filled or empty bands encountered.'
      write(ifief,"(2d23.15,a)") efermi,bandgap," ! efermi, bandgap are obtained by heftet"
      write(ifief,"(a)")'heftet: only filled or empty bands encountered.'
    else ! --- Find the Fermi energy by the Tetrahedron method ---
      write(stdo,*)' goto bzints2 for dos'
      allocate(dos(nptdos), dos_kbt(nptdos_kbt))
      volwgt = (3.d0 - nspin) / ntetf ! ntetf was =6*n1*n2*n3
      write (stdo,"(1x,'heftet : ---- Tetrahedron Integration ---- ')")
      do itt = 1, 20
        call bzints2x(volwgt,eband,dum111,nkp,nbandx,nbandx,nspx,e1,e2,dos,nptdos,efermi,1,nteti,idteti)
        if (dos(1) > valn) then ! Quick fix to enlarge window [e1 e2] for very few number of k points for metal (then efrang3 may not work). 2024-8-17
          e1=e1 -0.1d0
          cycle
        endif
        if (dos(nptdos) < valn) then
          e2=e2 +0.1d0
          cycle
        endif
        call fermi2( valn,dos,nptdos,e1,e2,   efermi,e1,e2,dosef) ! write(*,"(i3,7x,5(d13.6,1x))")itt, efermi, e1, e2, e2-e1, dosef
        if(e2 - e1 < tol) goto 444
      enddo
      call rx( 'heftet:bug in fermi level finder or tol too small')
444   continue
      deallocate(dos)
      bandgap=0d0
      write(ifief,"(2d23.15,a)") efermi,bandgap," ! efermi bandgap are obtained by heftet"
    endif
    ! Fermi energy at finite temperature (EFERMI_kbt). Written for metals AND
    ! insulators (the insulator branch above used to skip it, so hgw died on an
    ! empty EFERMI_kbt for every fresh insulator with tetrakbt=true, 2026-09-17).
    call gwinput_init()
    if (gwinput_loaded) then
      usetetrakbt = tg_tetrakbt
    else
      call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
    endif
    ! Used by m_tetwt when tetrakbt=true. (imode==1 is the gwsc GW flow.)
    if ((imode==1 .OR. imode==5) .AND. usetetrakbt) then
      call tetrakbt_init() !! read kbt
      FermiKbtTetraNOS: block ! EF^T from tetrahedron NOS (x) thermal kernel (2026-06-11).
        ! Replaces the discrete k-sum as the EFERMI_kbt source: T->0 reduces exactly to
        ! the tetrahedron EF (no staircase error at kbt below the mesh level spacing).
        ! For an insulator efermi is mid-gap; the NOS window then holds the gap and
        ! fermi_kbt_tetra returns the middle of the degenerate solution interval.
        integer,parameter:: nptk=2401
        real(8):: ek1, ek2, dosk(nptk), efkbt_disc
        volwgt = (3.d0 - nspin) / ntetf
        ek1 = efermi - 60d0*kbt - 1d-9
        ek2 = efermi + 60d0*kbt + 1d-9
        call bzints2x(volwgt,eband,dum111,nkp,nbandx,nbandx,nspx,ek1,ek2,dosk,nptk,efermi,1,nteti,idteti)
        call fermi_kbt_tetra(valn, dosk, nptk, ek1, ek2, kbt, efermi, efermi_kbt)
        call fermi_kbt_discrete(valn,eband2,nbandx,nspx,nqbz,nspin,kbt,elo,ehi, efkbt_disc) ! diagnostic
        write(stdo,ftox)' heftet: EFERMI_kbt(tetraNOSconv)=',ftof(efermi_kbt), &
             ' [diagnostic discrete k-sum=',ftof(efkbt_disc),'] kbt=',ftof(kbt)
      endblock FermiKbtTetraNOS
      open(newunit=ifief_kbt,file='EFERMI_kbt')
      write(ifief_kbt,"(2d23.15,a)") efermi_kbt, kbt,' ! efermi_kbt (finite-T, tetra-NOS conv) by heftet'
      close(ifief_kbt)
    endif
    close(ifief)
    if((imode/=1 .OR. imode/=5) .and. (any(eband(nbandx,:,:)<efermi)) ) &
         call rx(' heftet: WARNING! eband(maxband) is less than efermi: enlarge nbandx to get efermi by tetrahedron?')
    if(imode==1) then
      write(stdo,"(' Tet EFERMI gap = ',2f24.15)") efermi,bandgap
      call rx0( ' OK! heftet mode=1 EFERMI generated ')
    elseif(imode==2) then
      write(stdo,"(' Tet EFERMI.check gap= ',2f24.15)") efermi,bandgap
      call rx0( ' OK! heftet mode=2 EFERMI.check generated ')
    elseif(imode==3) then
      write(stdo,"(' Tet EFERMI.QP gap   = ',2f24.15)") efermi,bandgap
      call rx0( ' OK! heftet mode=3 EFERMI.QP generated ')
    elseif(imode==4) then
      write(stdo,"(' Tet EFERMI.QPz=1 gap= ',2f24.15)") efermi,bandgap
      call rx0( ' OK! heftet mode=4 EFERMI.QPz=1 generated ')
    elseif(imode==5) then
      write(stdo,"(' Tet EFERMI EFERMI_kbt = ',2f24.15)") efermi,efermi_kbt
      call rx0( ' OK! heftet mode=5 EFERMI and EFERMI_kbt ')
    endif
  end subroutine heftet

  subroutine efrang3(nsp,nkp,nband,zval,eband, e1,e2,elo,ehi,bandgap) !Find range of Fermi energy.
    implicit none
    intent(in) ::    nsp,nkp,nband,zval,eband
    intent(out)::                              e1,e2,elo,ehi,bandgap
    ! e1=e2 is just middle of HOMO and LUMO for insulator
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp,nkp,nbmax,nband,eband
    !i   zval no. of valence electrons
    !o Outputs
    !o   e1,e2: e1 < ef < e2
    !o   elo, ehi:  lowest and highest band found
    !r Remarks
    !r    For an even no. of electrons ef is between the bottom of the
    !r    zval/2+1'th band and the top of the zval/2 'th band. If the
    !r    former is higher that the latter we have an insulator and e1=e2.
    !r    For an odd no. of electrons ef is between the bottom and top
    !r    of the (zval+1)/2 'th band.
    !r    For spin pol case: bottom of the zval+1'th band < ef <
    !r                       top of the zval'th band.
    !r                       If the bottom is higher that the top then
    !r                       we have an insulator and e1=e2.
    !r    If nsp is zero, the upper and lower band limits are returned
    ! ----------------------------------------------------------------------
    integer :: nsp,nkp,nband, i,j,ikp,isp,iba,nval,nbbot,nbtop
    real(8):: zval,e1,e2,eband(nband,nsp,nkp), ebbot(nband*nsp),ebtop(nband*nsp),elo,ehi,em, e,d1mach, bandgap
    logical ::  zvalisinteger
    elo = 1d99
    ehi = -elo
    nval = nint(zval)
    zvalisinteger =  nval == zval
    if(.not.zvalisinteger) write(stdo,ftox) 'efrang3, comment:  zval is not an integer'
    if(nsp==1) then ! --- nsp = 1 ---
      if ( zvalisinteger ) then
        if (nval == 0) then
          nbtop = 1
          nbbot = 1
        else if ( (nval/2) * 2 == nval ) then  ! even integer
          nbtop = nval / 2
          nbbot = nbtop + 1
        else
          nbtop = nval / 2 + 1    ! odd integer or others
          nbbot = nbtop
        endif
      else    ! zval is not integer
        nbbot = int(zval/2) + 1 !nval / 2  +1
        nbtop = nbbot
      endif
      e1 = eband(nbbot,1,1) ! e1 is for bottom of conduction.
      e2 = eband(nbtop,1,1) ! e2 is for top of valence.  for insulator. ! e1<e2 for metal
      do ikp = 2 , nkp
        if( eband(1,1,ikp) < elo ) elo = eband(1,1,ikp)
        if( eband(nband,1,ikp) > ehi ) ehi = eband(nband,1,ikp)
        if( eband(nbbot,1,ikp) < e1 ) e1 = eband(nbbot,1,ikp)
        if( eband(nbtop,1,ikp) > e2 ) e2 = eband(nbtop,1,ikp)
      enddo
      if( e1 > e2 ) then ! if e1>e2 . This is insulator case.
        bandgap = e1-e2
        em = 0.5d0*(e1+e2)
        e1 = em
        e2 = em
      else
        bandgap=0d0
      endif
    else
      if (zvalisinteger) then
        do  isp = 1, nsp
          do  iba = 1, nband
            ebbot(iba+(isp-1)*nband) = eband(iba,isp,1)
            ebtop(iba+(isp-1)*nband) = eband(iba,isp,1)
          enddo
        enddo
      else
        do   isp = 1, nsp
          do   iba = 1, nband
            ebbot(iba+(isp-1)*nband) = minval(eband(iba,isp,1:nkp))
            ebtop(iba+(isp-1)*nband) = maxval(eband(iba,isp,1:nkp))
          enddo
        enddo
      endif
      do   ikp = 1, nkp
        do   isp = 1, nsp
          do   iba = 1, nband
            if( eband(1,isp,ikp)     < elo ) elo = eband(1,isp,ikp)
            if( eband(nband,isp,ikp) > ehi ) ehi = eband(nband,isp,ikp)
            if( eband(iba,isp,ikp) < ebbot(iba+(isp-1)*nband) ) ebbot(iba+(isp-1)*nband) =  eband(iba,isp,ikp)
            if( eband(iba,isp,ikp) > ebtop(iba+(isp-1)*nband) ) ebtop(iba+(isp-1)*nband) =  eband(iba,isp,ikp)
          enddo
        enddo
      enddo
      do  i = 1, nband*nsp - 1
        do  j = 1, nband*nsp - i
          if( ebbot(j) > ebbot(j+1) ) then
            e = ebbot(j)
            ebbot(j) = ebbot(j+1)
            ebbot(j+1) = e
          endif
          if( ebtop(j) > ebtop(j+1) ) then
            e = ebtop(j)
            ebtop(j) = ebtop(j+1)
            ebtop(j+1) = e
          endif
        enddo
      enddo
      if (zvalisinteger) then
        e1 = ebbot(nval+1)
        if( nval+1 > nband*nsp ) e1 = ebbot(nval)
        e2 = ebtop(nval)
      else
        e1 = ebbot(nval+1)
        e2 = ebtop(nval+1)
      endif
      if( e1 > e2 ) then !insulator case
        bandgap = e1-e2
        em = 0.5d0*(e1+e2)
        e1 = em
        e2 = em
      else
        bandgap=0d0
      endif
    endif
  end subroutine efrang3

  subroutine fermi_kbt_discrete(valn,eband,nb,nsp,nq,nspin,kbt,elo,ehi, efkbt)
    ! Finite-temperature Fermi level by a discrete k-sum (no DOS/tetrahedron):
    !   N(EF) = (3-nspin)/nq * sum_{s=1..nsp, iq=1..nq, n=1..nb} f(e-EF;kbt) = valn,
    ! solved by bisection (N is monotonic in EF). eband is the FULL-BZ eigenvalue
    ! array (equal weight 1/nq). Consistent with the finite-T chi0 (tetrakbt).
    implicit none
    intent(in)::  valn,eband,nb,nsp,nq,nspin,kbt,elo,ehi
    intent(out):: efkbt
    integer:: nb,nsp,nq,nspin, is,iq,n,it
    real(8):: valn,eband(nb,nsp,nq),kbt,elo,ehi,efkbt, a,b,c,fc,sf,x
    sf = dble(3-nspin)/dble(nq)            ! spin degeneracy * k-weight
    a = elo - 30d0*kbt ; b = ehi + 30d0*kbt
    do it = 1,300
      c = 0.5d0*(a+b)
      fc = 0d0
      do is=1,nsp ; do iq=1,nq ; do n=1,nb
        x = (eband(n,is,iq)-c)/kbt
        if    (x >  40d0) then ; continue
        elseif(x < -40d0) then ; fc = fc + 1d0
        else                   ; fc = fc + 1d0/(exp(x)+1d0)
        endif
      enddo ; enddo ; enddo
      fc = fc*sf
      if(fc > valn) then ; b = c ; else ; a = c ; endif
      if(b-a < 1d-12) exit
    enddo
    efkbt = 0.5d0*(a+b)
  end subroutine fermi_kbt_discrete
  subroutine fermi_kbt_tetra(valn, nos, npt, ea, eb, kbt, ef0, efkbt)
    ! Finite-temperature Fermi level from the TETRAHEDRON NOS (2026-06-11):
    !   N_T(E) = int de DOS(e) f_T(e-E) = int de NOS(e) (-f'_T(e-E))
    !          = sum_i ker_i * NOS(E + 2*kbt*t_i)        [B'-type GL20 kernel]
    ! solved by bisection. T->0 reduces EXACTLY to the tetrahedron EF (same NOS),
    ! removing the discrete-staircase O(level spacing) error of fermi_kbt_discrete
    ! at kbt below the coarse-mesh level spacing (e.g. 30K on 4^3/6^3 meshes).
    ! nos(1:npt): tetrahedron NOS on the linear grid [ea,eb] (from bzints2x job=1).
    use m_fpiint,only: gausq
    implicit none
    intent(in)::  valn, nos, npt, ea, eb, kbt, ef0
    intent(out):: efkbt
    integer:: npt, i, it
    integer,parameter:: NE=20
    real(8):: valn, nos(npt), ea, eb, kbt, ef0, efkbt
    real(8):: tg(NE), wg(NE), ker(NE), a,b,c, fc, de, eup, elo
    call gausq(NE, -6d0, 6d0, tg, wg, 0, 0)
    ker = wg*0.5d0/cosh(tg)**2
    ker = ker/sum(ker)
    de  = (eb-ea)/(npt-1)
    ! N_T(E) is flat (= valn) across a gap wider than the kernel's reach, so the
    ! solution is an interval there: bisect to its upper edge (first E with
    ! N_T > valn) and to its lower edge (last E with N_T < valn) and take the
    ! middle. For a metal both edges coincide.
    eup = edge(.true.)
    elo = edge(.false.)
    efkbt = 0.5d0*(eup+elo)
  contains
    real(8) function edge(upper)
      logical,intent(in):: upper
      real(8),parameter:: eps = 1d-10
      a = ef0 - 45d0*kbt
      b = ef0 + 45d0*kbt
      do it = 1,200
        c  = 0.5d0*(a+b)
        fc = 0d0
        do i=1,NE
          fc = fc + ker(i)*nosi(c + 2d0*kbt*tg(i))
        enddo
        if(upper) then
          if(fc > valn*(1d0+eps)) then; b = c; else; a = c; endif
        else
          if(fc < valn*(1d0-eps)) then; a = c; else; b = c; endif
        endif
        if(b-a < 1d-12) exit
      enddo
      edge = 0.5d0*(a+b)
    end function edge
    real(8) function nosi(e) ! clamped linear interpolation of nos on [ea,eb]
      real(8),intent(in):: e
      real(8):: p
      integer:: j
      if(e <= ea) then;     nosi = nos(1)
      elseif(e >= eb) then; nosi = nos(npt)
      else
        p = (e-ea)/de;  j = int(p)+1;  p = p-(j-1)
        nosi = (1d0-p)*nos(j) + p*nos(j+1)
      endif
    end function nosi
  end subroutine fermi_kbt_tetra
end module m_heftet
