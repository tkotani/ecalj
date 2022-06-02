module m_mkehkf
  use m_ftox
contains
  subroutine m_mkehkf_etot1(sev, etot)
    use m_struc_def
    use m_lmfinit,only: nsp,stdo,stdl
    use m_mkpot,only: utot,rhoexc,xcore,valvef,amom !readonly
    use m_mkrout,only: sumec,sumtc,sumt0 !readonly
    !- Make Harris energy or Kohn-Sham energy
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :1 Harris Foulkes energy etot=eh
    !i         :2 Kohn-Sham energy      etot=eks
    !i   sev   :band sum
    ! o Inputs/Outputs
    ! o  sumtv :valence kinetic energy :
    ! o        :output for mode=1
    ! o        :input for mode=2
    !o Outputs
    !o   etot  :Harris energy (mode = 1)
    !o         :Hohnberg-Kohn-Sham energy (mode = 2)
    implicit none
    integer :: mode
    integer:: i_copy_size
    real(8):: sev , etot ,ekin,ttcor, eh,eks,sumtv
    integer :: ipr,ipl
    call tcn('m_mkehkf_etot1')
    call getpr(ipr)
    ipl = 1
    sumtv = sev - valvef
    ttcor = sumec - xcore + sumt0
    ekin  = sumtv + ttcor
    etot  = ekin + utot + rhoexc
    eh    = etot
    if (ipr >= 30) then
       write(stdo,660) sev,valvef,sumtv,sumec,xcore,ttcor,rhoexc,utot, &
            eh
660    format(/' Harris energy:' &
            /' sumev= ',f15.6,'  val*vef=',f15.6,'   sumtv=',f15.6 &
            /' sumec= ',f15.6,'  cor*vef=',f15.6,'   ttcor=',f15.6 &
            /' rhoeps=',f15.6,'     utot=',f15.6,'    ehar=',f15.6)
    elseif (ipr >= 10) then
       write(stdo,ftox)' ekin=',ftof(ekin),'rho*v=',ftof(utot+rhoexc),'sumev=',ftof(sev), &
            'ehf=',ftof(eh)
    endif
    if (ipr > 0 .AND. ipl > 0) write (stdl,720) utot,ekin,sev,etot
720 format('fp Har: U',f15.6,'  T',f15.6,'  sev',f14.6,'  EH ',f14.6)
    call tcx('m_mkehkf_etot1')
  end subroutine m_mkehkf_etot1
  !!
  subroutine m_mkehkf_etot2(sev,sumtv, etot)
    use m_struc_def
    use m_lmfinit,only: nsp,stdo,stdl
    use m_mkpot,only: utot,rhoexc,xcore,valvef,amom !readonly
    use m_mkrout,only: sumec,sumtc,sumt0 !readonly
    !- Make Harris energy or Kohn-Sham energy
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :1 Harris Foulkes energy etot=eh
    !i         :2 Kohn-Sham energy      etot=eks
    !i   sev   :band sum
    ! o Inputs/Outputs
    ! o  sumtv :valence kinetic energy :
    ! o        :output for mode=1
    ! o        :input for mode=2
    !o Outputs
    !o   etot  :Harris energy (mode = 1)
    !o         :Hohnberg-Kohn-Sham energy (mode = 2)
    implicit none
    integer :: mode
    integer:: i_copy_size
    real(8):: sev , etot ,ekin,ttcor, eh,eks,sumtv
    integer :: ipr,ipl
    call tcn('m_mkehkf_etot2')
    call getpr(ipr)
    ipl = 1
    ekin = sumtv + sumtc
    etot = ekin + utot + rhoexc
    eks  = etot
    if (ipr >= 30) then
       if ( .FALSE. ) then !rhosig /= -99 .AND. rhosig /= 0) then
       else
          write (stdo,410) sumtv,sumtc,ekin,rhoexc,utot,eks
410       format(/' Kohn-Sham energy:' &
               /' sumtv= ',f15.6,'  sumtc=',f17.6,'   ekin=',f16.6 &
               /' rhoep= ',f15.6,'   utot=',f17.6,'   ehks=',f16.6)
       endif
       if (nsp == 2) write (stdo,412) amom
412    format(' mag. mom=',f13.6) !,'  (bands)',f16.6,'  (output rho)')

    elseif (ipr >= 10) then
       write(stdo,ftox)' ekin=',ftof(ekin),'rho*v=',ftof(utot+rhoexc),'ehf=',ftof(eh), &
            'ehks=',ftof(eks) !%,6;6d',' ',80,stdo,ekin,,eh,eks)
    endif
    if (ipr > 0 .AND. ipl > 0) write (stdl,721) utot,ekin,rhoexc,eks
721 format('fp KS: U',f16.7,' T',f16.7,' Exc',f15.7,' EKS',f15.7)
    call tcx('m_mkehkf_etot2')
  end subroutine m_mkehkf_etot2
end module m_mkehkf
