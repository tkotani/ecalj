!! vbmmode !Get VBM and CBM relative to vaccum (a simple approximaiton. need fixing.).
subroutine vbmmode()
  !! write vbmcbm.*
  use m_lmfinit,only: sspec=>v_sspec,ispec,nbas,ctrl_nspec,vol
  use m_ext,only:sname
  implicit none
  character(120):: vbmlll
  real(8):: rydberg=13.6058d0,esold
  integer:: ifvesatm,ifves,ifvesintloc,ib,ifvbm,ifvesintatm,ispe
  real(8):: vessm,sumvesloc,sumvesatm,vref,vesloc,evbm,ecbm
  real(8),allocatable::vesatm(:)
  logical:: lfill=.false.,ixx
  !     ! === VBM CBM section. Mar2013takao === just output
  print *,'--vbmonly mode'
  !      vol = lat_vol
  open(newunit=ifvesintatm,file='vesintatm.'//trim(sname),status='old',err=9898)
  !! vesintatm is given by lmfa. electrostatic potential integrals.
  allocate(vesatm(ctrl_nspec))
  do ispe=1,ctrl_nspec
     if(sspec(ispe)%z==0d0 .AND. sspec(ispe)%rmt==0d0 ) cycle
     read(ifvesintatm,*,err=9898,end=9898) vesatm(ispe)
  enddo
  close(ifvesintatm)
  open(newunit=ifves,file='vessm.'//trim(sname),status='old')
  read(ifves,*) vessm       !smooth part of electrostatic pot. given in smves.F. (es on MT average).
  close(ifves)
  open(newunit=ifvesintloc,file='vesintloc.'//trim(sname),status='old',err=9898)
  sumvesloc=0d0
  sumvesatm=0d0
  do ib=1,nbas
     !ispec=ssite(ib)%spec
     if(sspec(ispec(ib))%z==0d0 .AND. sspec(ispec(ib))%rmt==0d0 )  cycle
     read(ifvesintloc,*) vesloc
     sumvesloc = sumvesloc + vesloc
     sumvesatm = sumvesatm + vesatm(ispec(ib))
  enddo
  close(ifvesintloc)
  deallocate(vesatm)
  write(*,*)
  !     print *,'vessm=',vessm
  vref= sumvesatm/vol - vessm -sumvesloc/vol
  write(*,"('### VBM: Add vref to eigval to estimate eV relative to vaccum. vref(eV)=',f12.6)") vref*Rydberg
  write(*,"('### VBM: Mean estatic pot by superposition of atoms(eV)=',f12.6)")sumvesatm/vol*rydberg
  !      if(cmdopt0('--vbmonly')) then
  open(newunit=ifvbm,file='vbmcbm.'//trim(sname))
  read(ifvbm,"(a)") vbmlll
  print *,'readin vbmbm--> ',vbmlll
  read(vbmlll(9:),*)evbm
  read(vbmlll(30:),*)ecbm
  read(vbmlll(59:),*)esold
  write(*,*)evbm,ecbm,esold
  evbm= evbm +sumvesatm/vol*rydberg-esold
  ecbm= ecbm +sumvesatm/vol*rydberg-esold
  !      elseif(lfill) then
  !         evbm=(evtop+vref)*rydberg
  !         ecbm=(ecbot+vref)*rydberg
  !      else
  !         evbm=(eferm+vref)*rydberg
  !         ecbm=(eferm+vref)*rydberg
  !      endif
  open(newunit=ifvbm,file='vbmcbm.'//trim(sname))
  write(*,    "('### VBM: VBM= ',f12.6,' eV CBM= ',f12.6,' eV estatic_pot=',f12.6,' eV')") &
       evbm, ecbm, sumvesatm/vol*rydberg
  write(ifvbm,"('### VBM: VBM= ',f12.6,' eV CBM= ',f12.6,' eV estatic_pot=',f12.6,' eV')") &
       evbm, ecbm, sumvesatm/vol*rydberg
  close(ifvbm)
9898 continue
  print *,' if vesintatm.'//trim(sname)//' ---> VBM: is not shown. Need to repeat new lmfa for VBM:'
end subroutine vbmmode
