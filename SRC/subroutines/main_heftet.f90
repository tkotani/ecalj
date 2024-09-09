module m_heftet
  contains
subroutine heftet() bind(C)! Calculates the Fermi energy by tetrahedron method. 
  use m_read_bzdata,only: read_bzdata, idteti,qbz,qibz,dq_,nqibz,ntetf,nteti,ginv,nqbz
  use m_genallcf_v3,only: genallcf_v3, nclass,natom,nspin,nl,nn,nnv,nnc, nlmto,nlnmx, nctot,niw
  use m_genallcf_v3,only: alat, delta,deltaw,esmr,clabl,iclass, plat, pos,z,ecore, konf,nlnx,valn=>qval
  use m_hamindex,only:   Readhamindex,qtt,nqtt
  use m_readeigen,only: init_readeigen,readeval
  use m_tetrakbt,only: tetrakbt_init, kbt
  use m_keyvalue,only: getkeyvalue
  use m_genallcf_v3,only: nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
!  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg  
  use m_mpi,only: MPI__Initialize
  use m_lgunit,only: m_lgunit_init
  use m_ftox
  use m_lgunit,only: stdo
  use m_bzints,only: bzints2x,efrang3
  implicit none
  integer :: mxclass,ngnmax, ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd, &
       nprecx,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx, ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq, &
       i,k, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,iimdim, ifreqx,ifreqw,iwx,iexpa,mxkp,nqibzxx,ntet,nene,iqi, ix,iw, &
       nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx, nq2,ntq2,     ndble=8,ifev(2)
  real(8) ::tpia,vol,voltot,rs,alpha, qfermi,efx,efnew,edummy,efz,qm,xsex,egex, &
       zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,efxx2,        lowesteb
  logical :: lqall
  integer,allocatable :: itq(:)
  real(8),allocatable    :: q(:,:)
  integer,allocatable :: ngvecpB(:,:,:),ngveccB(:,:,:), ngvecp(:,:), ngvecc(:,:),ngpn(:),ngcni(:),iqib(:), &
       kount(:,:), nx(:,:),nblocha(:),lx(:),ngveccBr(:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), wqt(:), wgt0(:,:),q0i(:,:), &
       ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:), ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:,:), &
       ppb(:),pdb(:),dpb(:),ddb(:), eq(:,:), eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
  complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)
  logical :: exchange, legas
  real(8):: qreal(3), tripl!ntot,nocctotg2,
  logical ::nocore
  integer(4) ::ib,igp,iii,ivsumxxx,isx
  real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
  integer,allocatable:: ieord(:),itex1(:,:,:)
  real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,qq(3),        elda, eqp01, eqp02
  integer:: itmx,ipex,itpex,itex,nnex,isig,iex,ifexspx,ifexspxx , itx, iqx
  character(12) :: filenameex
  logical :: exspwrite=.false.
  integer:: nptdos=101 !,iflegas
  logical :: metal      ,qbzreg
  integer:: nkp,itt,ifief ,imode,iftote,it
  real(8) :: elo,ehi,e1,e2,efermi,dum,dosef,efxx,rydberg,dum111(1,1,1)  , tol=1d-10,toql=1d-8, volwgt
  real(8),allocatable:: dos(:)
  integer:: nwin,bzcase=1,ifi
  real(8)::ddq(3),bandgap,tolq=1d-8   ! space group infermation
  integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
  real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)   ! tetra
  real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:), eband(:,:,:),eband2(:,:,:), ene(:)
  integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:),iene(:,:,:),ibzx(:)   ! fermi_kbt
  integer:: nptdos_kbt = 100001, ifief_kbt, nbandx
  real(8), allocatable:: dos_kbt(:), dosef_kbt
  real(8) :: efermi_kbt, e11, e22
  real(8),parameter:: pi= 4d0*datan(1d0)
  logical:: usetetrakbt, cmdopt2
  character(20):: outs=''
  call MPI__Initialize()
  call M_lgunit_init()
  if(cmdopt2('--job=',outs)) then
     read(outs,*) imode
  else
     write(6,*) 'mode=(1-5)?'
     read(5,*) imode
  endif
  write(6,*) '--- heftet: calculation mode =',imode
  if(imode<1 .OR. imode>5) call rx( 'mode out of range(1-4)')
  call read_BZDATA()
  write(6,ftox)' heftet: nqibz, ntetf,nteti,nqbz ', nqibz,ntetf,nteti,nqbz
  call genallcf_v3(incwfx=-1) !in module m_genallcf_v3
!  call Readhbe()   !      write(6,"(' nqbz nqibz ngrp=',3i8)") nqbz,nqibz,ngrp  ! if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
  if(nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
  tpia = 2d0*pi/alat
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  ! ef is taken as rs for the empty-sphere test case of legas=T case -------------
  !$$$      legas = .false.
  !$$$      INQUIRE (FILE = 'LEGAS', EXIST = legas)
  !$$$      if(legas) then !!! test for electron gas case.
  !$$$        write(6,*)' find LEGAS. legas =',legas
  !$$$        iflegas = 2101
  !$$$        open (iflegas,file='LEGAS')
  !$$$        read(iflegas,*)rs
  !$$$        close(iflegas)
  !$$$        alpha = (9*pi/4d0)**(1d0/3d0)
  !$$$        qfermi = alpha/rs
  !$$$        efx  = qfermi**2
  !$$$        valn = efx**1.5d0*voltot/3d0/pi**2
  !$$$        write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
  !$$$        write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
  !$$$        write (6,*)' egas  Exact Fermi energy    Ef  =', efx
  !$$$      endif
  write(6, *) ' --- computational conditions --- '
  write(6,'("    deltaw  =",f13.6)') deltaw
  write(6,'("    esmr    =",f13.6)') esmr
  write(6,'("    alat voltot =",2f13.6)') alat, voltot
!  open(newunit=ifi,file='efermi.lmf')
!  read(ifi,*)
!  read(ifi,*)
!  read(ifi,*)
!  read(ifi,*) valn !Get number of valence electron valn
!  close(ifi)
  if(imode==1 .OR. imode==5) then
     call Readhamindex()
     call init_readeigen() !initialization of readEigen
     nbandx= nband
     nkp   = nqibz
     allocate( eband2(nbandx,nspin,nqbz) )
     ddq = 0d0
     if( .NOT. qbzreg()) ddq= -dq_
     write(6,"(' nsp nqbz, q is replaced by q+dq_; dq_ =',2i5,3f13.5)") nspin,nqbz,dq_
     do is  = 1,nspin !Readin eband
        ix = 0
        do iqi = 1,nqbz
           ix = ix+1
           eband2(:,is,ix)= readeval(qbz(:,iqi)+ddq, is)
        enddo
     enddo
     call efrang3(nspin,nqbz,nbandx,valn,eband2,  e1,e2,elo,ehi,bandgap)
     if(qbzreg() .AND. bzcase==1 ) then
        allocate( eband(nbandx,nspin,nqibz))
        do is  = 1,nspin
           do iqi = 1,nqibz
              eband(:,is,iqi)=readeval(qibz(:,iqi),is)
              write(6,"(' -------- qibz=',2i5,3f9.4)")iqi,is,qibz(:,iqi)
           enddo
        enddo
     endif
  elseif(imode==2 .OR. imode==3 .OR. imode==4) then
     open(newunit=iftote,file='TOTE.UP')
     read(iftote,*) nq,ntq, efxx
     if(nq/= nqibz) call rx( ' heftet: nq TOTE/= nqibz')
     nbandx = ntq
     allocate( eband(nbandx,nspin,nqibz), qz(3,nqibz) )
     qz=qibz 
     do  is  = 1,nspin !Readin eband
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
                 write(6,*) 'iqi it=',iqi,it
                 write(6,*) qq(1:3)
                 write(6,*) qz(1:3,iqi)
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
     call efrang3(nspin,nkp,nbandx,valn,eband,  e1,e2,elo,ehi,bandgap)
  endif
  write(6,"(' end of efrang3: e1 e2 nteti=',2d15.7,i6)") e1,e2,nteti
  write(6,"(' elo ehi = ',2d15.7,i6)") elo,ehi
  write(6,*)'e1=',e1
  write(6,*)'e2=',e2
  if(e1 /= e2) then
     if( .NOT. ( qbzreg() .AND. bzcase==1 ) ) then
        write(6,*)"Fermi energy finder for off-regular mesh (Chi_RegQbz off)or bzcase==2 is not implimented yet..."
        call rx("Ef finder for off-regular mesh. NotYet implimented")
     endif
  endif
  if(imode==1) open(newunit=ifief,file='EFERMI')
  if(imode==5) open(newunit=ifief,file='EFERMI')
  if(imode==2) open(newunit=ifief,file='EFERMI.check')
  if(imode==3) open(newunit=ifief,file='EFERMI.QP')
  if(imode==4) open(newunit=ifief,file='EFERMI.QPz=1')
  if(e1 == e2) then ! --- Find band gap  e1=e2 is at the middle of gap---
     efermi     = e1
     write(6,*)'heftet: only filled or empty bands encountered.'
     write(ifief,"(2d23.15,a)") efermi,bandgap," ! efermi, bandgap are obtained by heftet"
     write(ifief,"(a)")'heftet: only filled or empty bands encountered.'
  else ! --- Find the Fermi energy by the Tetrahedron method ---
     write(6,*)' goto bzints2 for dos'
     allocate(dos(nptdos), dos_kbt(nptdos_kbt))
     volwgt = (3.d0 - nspin) / ntetf ! ntetf was =6*n1*n2*n3
     write (6,"(1x,'heftet : ---- Tetrahedron Integration ---- ')")
     do itt = 1, 20
        call bzints2x(volwgt,eband,dum111,nkp,nbandx,nbandx,nspin,e1,e2,dos,nptdos,efermi,1,nteti,idteti)
        if (dos(1)      > valn) then ! Quick fix to enlarge window [e1 e2] for very few number of k points for metal (then efrang3 may not work). 2024-8-17
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
     write(*,*) 'heftet:bug in fermi level finder or tol too small'
     call rx( 'heftet:bug in fermi level finder or tol too small')
444  continue
     ! Fermi energy at finite temperature (EFERMI_kbt); Okumura (Feb.2020)
     call getkeyvalue("GWinput","tetrakbt",usetetrakbt,default=.false.)
     if (imode==5 .AND. usetetrakbt) then
        call tetrakbt_init() !! read kbt
        e11 = efermi - 10*kbt
        e22 = efermi + 10*kbt
        nptdos_kbt = 10001
        call bzints2x(volwgt,eband,dum111,nkp,nbandx,nbandx, nspin,e11,e22,dos_kbt, nptdos_kbt ,efermi,1,nteti,idteti)
        call fermi_kbt(valn,dos_kbt,nptdos_kbt,e11,e22, kbt, efermi, efermi_kbt) !Fermi level at finite tempearture.
        deallocate(dos_kbt)
        open(newunit=ifief_kbt,file='EFERMI_kbt')
        write(ifief_kbt,"(2d23.15,a)") efermi_kbt, kbt,' ! This efermi kbt are obtained by heftet:'
        close(ifief_kbt)
     endif
     deallocate(dos)
     bandgap=0d0
     write(ifief,"(2d23.15,a)") efermi,bandgap," ! efermi bandgap are obtained by heftet"
  endif
  close(ifief)
  ebandchek: if(imode/=1 .OR. imode/=5) then
     if(any(eband(nbandx,:,:)<efermi)) then
        write(6,*)' heftet: WARNING! eband(maxband) is less than efermi: enlarge nbandx to get efermi by tetrahedron?'
        goto 666
     endif   
     write(6,*) ' check OK! eband(nbandx,:,:) is greater than efermi.'
666  continue
  endif ebandchek
  if(imode==1) then
     write(6,"(' Tet EFERMI gap = ',2f24.15)") efermi,bandgap
     call rx0( ' OK! heftet mode=1 EFERMI generated ')
  elseif(imode==2) then
     write(6,"(' Tet EFERMI.check gap= ',2f24.15)") efermi,bandgap
     call rx0( ' OK! heftet mode=2 EFERMI.check generated ')
  elseif(imode==3) then
     write(6,"(' Tet EFERMI.QP gap   = ',2f24.15)") efermi,bandgap
     call rx0( ' OK! heftet mode=3 EFERMI.QP generated ')
  elseif(imode==4) then
     write(6,"(' Tet EFERMI.QPz=1 gap= ',2f24.15)") efermi,bandgap
     call rx0( ' OK! heftet mode=4 EFERMI.QPz=1 generated ')
  elseif(imode==5) then
     write(6,"(' Tet EFERMI EFERMI_kbt = ',2f24.15)") efermi,efermi_kbt
     call rx0( ' OK! heftet mode=5 EFERMI and EFERMI_kbt ')
  endif
end subroutine heftet
end module m_heftet
