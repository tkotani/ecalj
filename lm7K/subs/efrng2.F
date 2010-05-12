      logical function efrng2(nsp,nkp,nbmax,nband,zval,eband,ebbot,
     .  ebtop,elo,ehi,e1,e2)
C- Find range of Fermi energy.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nbmax :leading dimension of eband
Ci   nband :number of bands
Ci   zval  :no. of valence electrons
Ci   eband :energy bands
Co Outputs
Co   e1,e2: e1 < ef < e2
Co   elo, ehi:  lowest and highest band found
Co   efrng2:flags whether metal or insulator.
Co         :F metal (highest occ band crossed lowest unoccupied one)
Co         :T insulator (highest occ band did not cross lowest unocc)
Cr Remarks
Cr    For an even no. of electrons ef is above the bottom of the
Cr    zval/2+1'th band and below the top of the zval/2 'th band. If the
Cr    former is higher that the latter we have an insulator, with
Cr    these two numbers estimates for the conduction band minimum and
Cr    valence band maximum, respectively.
Cr    For an odd no. of electrons ef is between the bottom and top
Cr    of the (zval+1)/2 'th band.
Cr
Cr    For spin pol case:
Cr      bottom of the zval+1'th band < ef < top of the zval'th band.
Cr      If the bottom is higher that the top then we have an insulator,
Cr      with the bottom an estimate for the conduction band minimum and
Cr      to and estimate for the valence band maximum.
Cr      and e1=e2.
Cu Updates
Cu   08 Jul 08 Extend to case where number of bands can be q dependent
Cu             eband=99999 => not calculated: ignore
Cu   01 Apr 03 Set e1=e2 if every band is filled
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nsp,nkp,nbmax,nband
      double precision zval,e1,e2,eband(nbmax,nsp,nkp),ebbot(nband,nsp),
     .                 ebtop(nband,nsp),elo,ehi
C Local parameters
      double precision xx,d1mach,enull
      integer ikp,isp,iba,nval,nbbot,nbtop,nfound
      parameter (enull=99999d0)

      elo = enull
      ehi = -enull
      nval = zval + 1d-7

C --- Find bottom and top of each band ---
      call dvset(ebbot,1,nband*nsp,elo)
      call dvset(ebtop,1,nband*nsp,ehi)
      do  ikp = 1, nkp
        do  isp = 1, nsp
          do  iba = 1, nband
            if (eband(iba,isp,ikp) .ne. enull) then
              ebbot(iba,isp) = min(ebbot(iba,isp),eband(iba,isp,ikp))
              ebtop(iba,isp) = max(ebtop(iba,isp),eband(iba,isp,ikp))
            endif
          enddo
        enddo
      enddo
C     Set all -enull to enull to float to top when sorted
      do  isp = 1, nsp
        do  iba = 1, nband
          if (ebtop(iba,isp) .eq. -enull) ebtop(iba,isp) = enull
        enddo
      enddo
C     Sort bands irrespective of spin
      call dshell(nband*nsp,ebbot)
      call dshell(nband*nsp,ebtop)
      nfound = nband*nsp
   10 continue
      if (ebtop(nfound,1).eq.enull .or. ebbot(nfound,1).eq.enull) then
        nfound = nfound-1
        if (nfound .eq. 0) call rx('efrng2: no bands')
        goto 10
      endif

C --- Find limits ---
      nbtop = (nval+2-nsp)/(3-nsp)
      if (zval .gt. nval) nbtop = nbtop+1
      nbbot = nval/(3-nsp) + 1
      if (nbtop .gt. nfound) nbtop = nfound
      if (nbbot .gt. nfound) nbbot = nfound
      elo = ebbot(1,1)
      ehi = ebtop(nfound,1)
      if (elo .eq. enull) call rx('efrng2: no bands')
      e1  = ebbot(nbbot,1)
      e2  = ebtop(nbtop,1)
      efrng2 = .false.
C     if (e1 .gt. e2) then
      if (e1-e2 .gt. d1mach(3)) then
        xx = e1
        e1 = e2
        e2 = xx
        efrng2 = .true.
      endif
C     Every band is filled ??
C     if (nbbot .eq. nbtop) e1 = e2

C --- Printout ---
C      if (iprint() .ge. 50) then
C        print '(1x)'
C        call awrit3(' efrng2:  emin=%d  emax=%d'//
C     .    '  de=%d',' ',80,i1mach(2),e1,e2,e2-e1)
C      endif
      end

