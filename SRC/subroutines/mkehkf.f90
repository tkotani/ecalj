module m_mkehkf
  use m_ftox
contains
  subroutine m_mkehkf_etot1(sev, etot)
    use m_lmfinit,only: nsp,stdo,stdl
    use m_mkpot,only: utot,rhoexc,xcore,valvef,amom 
    use m_mkrout,only: sumec,sumt0 
    !- Make Harris energy 
    !i   sev   :band sum
    !o   etot  :Harris energy 
    implicit none
    real(8):: sev , etot ,ekin, eh,eks,sumtv,ekincore
    integer :: ipr
    call tcn('m_mkehkf_etot1')
    call getpr(ipr)
    sumtv = sev - valvef         !kineric energy for valence
    ekincore = sumec-xcore + sumt0  !kinetic energy (sumec-xcore is sum for lfoc=1 atoms, sumt0 is sum for lfoc=0 atoms)
    ! This is in the spilits of Harris energy.
    ! Test for Zrt shows this is advantatious than ttcor=sumt0 for lfoc=0
    ekin  = sumtv + ekincore
    etot  = ekin + utot + rhoexc
    eh    = etot
    if (ipr > 0) then
       write(stdo,660) sev,valvef,sumtv, ekincore,rhoexc,utot,eh !sumec,xcore,
660    format(/' m_mkehkf_etot1: Harris energy: (B.1) in JPSJ84,034702' &
            /' Eb(band sum)= ',f15.6,' Vin*nin=',f15.6,' Ek=Eb-Vin*nin=',f15.6 &
            /' Ek(core)=',f15.6,' Exc=',f15.6,' Ees=',f15.6,' Eharris=',f15.6)
       write (stdl,720) utot,ekin,sev,etot
720    format('fp Har: U',f15.6,'  T',f15.6,'  sev',f14.6,'  EH ',f14.6)
    endif
    call tcx('m_mkehkf_etot1')
  end subroutine m_mkehkf_etot1
  subroutine m_mkehkf_etot2(sumtv, etot) !Make Kohn-Sham energy
    use m_lmfinit,only: nsp,stdo,stdl,nbas,ispec
    use m_mkpot,only: utot,rhoexc,xcore,valvef,amom 
    use m_mkrout,only: sumtc
    !i  sumtv :valence kinetic energy 
    !o   etot :Hohnberg-Kohn-Sham energy
    implicit none
    real(8):: etot ,ekin,ttcor, eks,sumtv
    integer :: ipr
    call tcn('m_mkehkf_etot2')
    call getpr(ipr)
    ekin = sumtv + sumtc !sumtc is sum of the kinetic energies for cores
    etot = ekin + utot + rhoexc
    eks  = etot
    if (ipr > 0) then
       write (stdo,410) sumtv,sumtc,ekin,rhoexc,utot,eks
410    format(/' m_mkehkf_etot2: Kohn-Sham energy:  Ek = Eband-Vin*nout' &
            /' Ek= ',f15.6,' Ekcore=',f17.6,' Ektot    =',f16.6 &
            /' Exc=',f15.6,' Ees   =',f17.6,' EKohnSham=',f16.6)
       if (nsp == 2) write(stdo,"(' Magnetic moment=',f13.6)") amom
       write (stdl,721) utot,ekin,rhoexc,eks
721    format('fp KS: U',f16.7,' T',f16.7,' Exc',f15.7,' EKS',f15.7)
    endif
    call tcx('m_mkehkf_etot2')
  end subroutine m_mkehkf_etot2
end module m_mkehkf
