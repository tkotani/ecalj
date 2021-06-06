      subroutine mkehkf(mode,sev,sumtv,etot)
      use m_struc_def
      use m_lmfinit,only: nsp,stdo,stdl
      use m_energy_mkpot,only: utot,rhoexc,xcore,valvef,amom !readonly
      use m_mkrout,only: sumec,sumtc,sumt0 !readonly
C- Make Harris energy or Kohn-Sham energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 Harris Foulkes energy etot=eh
Ci         :2 Kohn-Sham energy      etot=eks
Ci   sev   :band sum
Cio Inputs/Outputs
Cio  sumtv :valence kinetic energy :
Cio        :output for mode=1
Cio        :input for mode=2
Co Outputs
Co   etot  :Harris energy (mode = 1)
Co         :Hohnberg-Kohn-Sham energy (mode = 2)
      implicit none
      integer mode
      integer:: i_copy_size
      real(8):: sev , etot , sumtv , ekin,ttcor, eh,eks
      integer ipr,ipl
      call getpr(ipr)
      ipl = 1
      if (mode .eq. 1) then
        sumtv = sev - valvef
        ttcor = sumec - xcore + sumt0
        ekin  = sumtv + ttcor
        etot  = ekin + utot + rhoexc
        eh    = etot
      elseif (mode .eq. 2) then
        ekin = sumtv + sumtc
        etot = ekin + utot + rhoexc
        eks  = etot
      endif
      if (mode .eq. 0) return
C --- Printout ---
      if (mode .eq. 1) then
        if (ipr .ge. 30) then
          write(stdo,660) sev,valvef,sumtv,sumec,xcore,ttcor,rhoexc,utot,
     .    eh
  660     format(/' Harris energy:'
     .    /' sumev= ',f15.6,'  val*vef=',f15.6,'   sumtv=',f15.6
     .    /' sumec= ',f15.6,'  cor*vef=',f15.6,'   ttcor=',f15.6
     .    /' rhoeps=',f15.6,'     utot=',f15.6,'    ehar=',f15.6)
        elseif (ipr .ge. 10) then
          call awrit4('%N ekin=%,6;6d  rho*v=%,6;6d  sumev=%,6;6d'//
     .    ' ehf =%,6;6d',' ',80,stdo,ekin,utot+rhoexc,sev,eh)
        endif
        if (ipr.gt.0 .and. ipl.gt.0) write (stdl,720) utot,ekin,sev,etot
  720   format('fp Har: U',f15.6,'  T',f15.6,'  sev',f14.6,'  EH ',f14.6)
      elseif (mode .eq. 2) then
         if (ipr .ge. 30) then
          if (.false.) then !rhosig .ne. -99 .and. rhosig .ne. 0) then
          else
            write (stdo,410) sumtv,sumtc,ekin,rhoexc,utot,eks
  410       format(/' Kohn-Sham energy:'
     .      /' sumtv= ',f15.6,'  sumtc=',f17.6,'   ekin=',f16.6
     .      /' rhoep= ',f15.6,'   utot=',f17.6,'   ehks=',f16.6)
          endif
          if (nsp .eq. 2) write (stdo,412) amom
  412     format(' mag. mom=',f13.6) !,'  (bands)',f16.6,'  (output rho)')

        elseif (ipr .ge. 10) then
          call awrit4('%N ekin=%,6;6d  rho*v=%,6;6d ehf =%,6;6d'//
     .    '  ehks =%,6;6d',' ',80,stdo,ekin,utot+rhoexc,eh,eks)
        endif
        if (ipr.gt.0 .and. ipl.gt.0) write (stdl,721) utot,ekin,rhoexc,eks
  721   format('fp KS: U',f16.7,' T',f16.7,' Exc',f15.7,' EKS',f15.7)
      endif
      end subroutine mkehkf


