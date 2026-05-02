!> Read values from GWinput
module m_readgwinput 
  use m_struct_from_lmf,only: nspin; use m_core_state,only: nctot
  implicit none
  real(8),protected:: egauss,ecut,ecuts,ebmx ,ebmx_sig,ua_
  integer,protected:: nbmx,nbmx_sig !,nbcutlow_sig !nbcut,nbcut2,
  ! iSigmode: declared but no external user; commented out 2026-05-02 (dead)
  !integer,protected:: iSigmode
  integer,protected:: mtet(3),nmbas
  integer,protected,allocatable:: imbas(:)
  logical,protected:: keeppositivecou
  logical,protected,public:: corehole
  real(8),allocatable:: wcorehole(:,:)
contains
  ! SetIsigmode: only callsite was inside this module; with iSigmode dead,
  ! this routine has no purpose. Commented out 2026-05-02.
  !subroutine SetIsigmode(ism)
  !  integer :: ism
  !  isigmode=ism
  !end subroutine SetIsigmode
  subroutine ReadGWinputKeys() ! Read values from GWinput
    use m_struct_from_lmf,only: natom,nband
    use m_GWinput, only: gwinput_init, gwinput_loaded, &
                         tg_ecut_p     => ecut_p, &
                         tg_ecuts_p    => ecuts_p, &
                         tg_nband_chi0 => nband_chi0, &
                         tg_emax_chi0  => emax_chi0, &
                         tg_multitet   => multitet, &
                         tg_nband_sigm => nband_sigm, &
                         tg_emax_sigm  => emax_sigm, &
                         tg_gauss_img  => gauss_img, &
                         tg_KPC        => KeepPositiveCou, &
                         tg_MagAtom    => MagAtom
    logical:: cmdopt0
    integer :: ifcorehole,it,nctot,nspin
    call gwinput_init()
    if (.not. gwinput_loaded) call rx('m_readgwinput: GWinput.toml is required.')
    ! ---- TOML path: copy from m_GWinput module variables ----
    ecut            = tg_ecut_p
    ecuts           = tg_ecuts_p
    nbmx            = tg_nband_chi0
    ebmx            = tg_emax_chi0
    mtet            = tg_multitet
    nbmx_sig        = tg_nband_sigm
    ebmx_sig        = tg_emax_sigm
    ua_             = tg_gauss_img
    KeepPositiveCou = tg_KPC
    if (allocated(tg_MagAtom)) then
       nmbas = size(tg_MagAtom)
       allocate(imbas(nmbas))
       imbas = tg_MagAtom
    else
       nmbas = 0
       allocate(imbas(0))
    endif
    if (nmbas>0) write(6,"('Readin MagAtom (TOML) nmbas =',i3,' imbas= ',10i3)") nmbas,imbas(1:nmbas)

!   ---- Legacy GWinput reader (disabled 2026-05-02 -- TOML only) ----
!   use m_keyvalue,only: Getkeyvalue
!   integer:: istat
!   integer,allocatable:: imbasd(:)
!   call Getkeyvalue("GWinput","ecut_p" ,ecut, default=1d10 )
!   call Getkeyvalue("GWinput","ecuts_p",ecuts,default=1d10 )
!   call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband)
!   call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
!   call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
!   call getkeyvalue("GWinput","nband_sigm",nbmx_sig, default=9999999)
!   call getkeyvalue("GWinput","emax_sigm", ebmx_sig, default=1d10)
!   call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
!   call getkeyvalue("GWinput","KeepPositiveCou",KeepPositiveCou,default=.true.)
!   nmbas=natom
!   allocate(imbas(nmbas),imbasd(nmbas))
!   imbasd = -9999
!   istat=-9999
!   call getkeyvalue("GWinput","MagAtom", imbas,nmbas,status=istat,default=imbasd)
!   nmbas = istat

    corehole = cmdopt0('--corehole')
    if(corehole) then
       open(newunit=ifcorehole,file='CoreHole')
       allocate(wcorehole(nctot,nspin))
       do it=1,nctot
          read(ifcorehole,*) wcorehole(it,1:nspin)
       enddo
       close(ifcorehole)
       write(*,*) 'corehole mode: end of reading CoreHole'
    endif

  end subroutine ReadGWinputKeys
end module m_readgwinput

