      subroutine tbtote(sctrl,slat,sarray,ssite,nl,nsp,idxdn,nclass,
     .  nrclas,dclabl,pv,force,nbas,iax,npr,iam,npm,V0,qnu,sumev,
     .  alpha,entrpy,emad,ecorr,r,stni,uh,jh,mmom,emag,e,eri,ppdip,
     .  thrpv,tpvq,f,fe,fmax,fnou,fnom,eatm,erep,etot,efree,emg,amgm)
C- Calculates tight-binding total energy and pair potential terms
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl,nsp,idxdn,nclass,nrclas,dclabl
Ci   iax, neighbor lists; npr, see tbham
Ci   iam,npm:  see remarks
Ci   V0, parameters for pair potential, see rdtbh
Ci   qnu(1),qnu(2): initial charge and site energies (1st and 2nd moms)
Ci   sumev, sum of occupied eigenvalues
Ci   alpha, see Phys Rev B, 53, 15381 (1996)
Ci   entrpy, entropy term (actually kTS)
Ci   emad, electrostatic energy (if U=T)
Ci   ecorr, correction to UL=T total energy
Ci   r,   s,p, and d Mulliken charges
Ci   mmom,emag, local moments and magnetic energies (Stoner model)
Ci   e,   band energies for each atom
Ci   eri, work array for atom decomposition of pair energy
Ci   ppdip: atomic contributions to dipole moment
Ci   thrpv,tpvq: 3pV from bands and multipoles respectively
Ci   f,   forces from bands; fe, forces from electrostatix
Ci   fnou: force from overlap dependence of monopoles (Hubbard part)
Ci   fnom: force from overlap dependence of monopoles (Madelung part)
Ci   ppdip, accumulated atomic dipole moments (from tbesel) 
Ci          RELAX=0 in ctrl is misused 
Ci          to permit the dipole to be accumulated over a subset
Ci          of the sites.
Co Outputs
Co   f,    force on each atom, has pair and e'static contributions added
Co   fmax, maximum force
Co   eatm, ref energy for each atom, unchanged if non-zero on entry
Co   erep, total energy contribution from pair potential
Co   etot, total energy
Co   efree, free energy: \Phi, eq (5), Phys Rev B, 53, 15381 (1996)
Co   emg,  magnetic energy (Stoner model)
Co   amgm, magnetic moment
Co   thrpv, 3pV: pV is the "internal virial" ie, (1/3)\sum r.f
Co         on output thrpv has pair and electrostatic terms added
Co   sumev is overwritten with Tr[rho H_0] if TRH is set
Cr Remarks
Cr   iam(1,kk) and iam(2,kk) are the two classes for the kkth ME pair
Cr   k = iam(3,kk) is a pointer to tabme(*,k) which holds the
Cr   matrix elements <iam(1,kk) | H | iam(2,kk)>
Cr   npm(0,i) is number matrix elts tabulated for ith class
Cr   npm(1,i) is number matrix elts in all sites preceding ith class
Cu Updates
Cu   04 Jun 08 (ATP) multipole contribution to 3PV
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nsp,nclass,nbas,niax
      parameter (niax=10)
      integer idxdn(nl,nclass),nrclas(1),iax(niax,2),npr(0:1,2),
     .  iam(3,2),npm(0:1,1)
      double precision sctrl(1),slat(1),sarray(1),ssite(1)
      double precision alat,sumev,alpha,entrpy,emad,thrpv,tpvq,erep,
     .                 etot,efree,amgm,emg,ecorr
      double precision plat(3,3),V0(9,2),qnu(3,nl,nsp,nclass),
     .  r(nl,2,nbas),e(nbas,nsp),eri(nbas),f(3,1),fe(3,1),fnou(3,1),
     .  fnom(3,1),stni(nclass),uh(4,1),jh(4,1),dclabl(nclass),
     .  eatm(nclass,nsp),mmom(nbas),emag(nbas),ppdip(3)
      double precision getavJ
      character*8 clabli,clablj
      logical pv,force
C Heap:
      integer w(1)
      common /w/ w
C Local variables
      integer i,k,n,iat,ic,l,ifc,ipr,iatmx,isp,lgunit,ltb,lstnr(3),
     .  nsgrp,oipc,opos,i1mach,ifrlx(3)
      double precision wk(0:3),erepi,derepi,eref,dv(3),dt(3),trh0(2),
     .  q(2),tpvpr,dp,fabs,fmax,d1mach,dsum,qtot(2),pos(3),
     .  pdip(3),ptot(3),dipole,dip1,dip2,q0,Dq,J,tpv
      logical bittst,lov,trh,TBU,UL,charge,stoner,MOL
      character*80 outs

      call upack('ctrl ltb lstonr nbas',sctrl,ltb,lstnr,nbas,0,0)
      call upack('lat alat plat nsgrp',slat,alat,plat,nsgrp,0,0)
      call upack('array oipc opos',sarray,oipc,opos,0,0,0)

      lov    = bittst(ltb,1)
      trh    = bittst(ltb,2**10)
      charge = bittst(ltb,2**11)
      UL     = bittst(ltb,2**15)
      TBU    = bittst(ltb,2**13)
C     Uav = (.not. bittst(ltb,2**14))
      MOL    = bittst(ltb,2**18)
      stoner = lstnr(1) .ne. 0
      if (UL .or. TBU) then
        emad = ecorr
        trh = .true.
      endif
      amgm = 0d0
      call getpr(ipr)

C --- Get reference energy ---
      eref = 0d0
      do  ic = 1, nclass
        if (dabs(eatm(ic,1))   .lt. 10*d1mach(3) .and.
     .      dabs(eatm(ic,nsp)) .lt. 10*d1mach(3)) then
          eatm(ic,1) = 0d0
          if (nsp .eq. 2) eatm(ic,2) = 0d0
          do  isp = 1, nsp
            do  l = 1, nl
              if (idxdn(l,ic) .gt. 1) goto 1
              eatm(ic,isp) = eatm(ic,isp)
     .                     + qnu(1,l,isp,ic)*qnu(2,l,isp,ic)
    1         continue 
            enddo
          enddo
        endif
        eref = eref + eatm(ic,1)*nrclas(ic)
        if (nsp .eq. 2) eref = eref + eatm(ic,2)*nrclas(ic)
      enddo

C --- Loop over all atoms ---
      fmax = 0d0
      erep = 0d0
      tpvpr = 0d0
      call dpzero(eri,nbas)
      if (ipr .gt. 20 .and. ipr .le. 30 .and. force) write (*,102)
      do  iat = 1, nbas
C --- Get pair potential contributions ---
        ic = w(oipc+iat-1)
        n = npr(1,iat)
        if (force) call dpzero(dv,3)
        do  i = 2, npr(0,iat)
          call dlmn(nbas,plat,w(opos),iax(1,i+n),wk)
          call meptr(w(oipc+iax(1,i+n)-1),w(oipc+iax(2,i+n)-1),
     .               iam,npm,k)
          if (k .eq. 0) goto 2
          call makv0(wk(0)*alat,v0(1,k),erepi,derepi,dp)
          if (force) call daxpy(3,derepi,wk(1),1,dv,1)
          erep = erep + erepi
          eri(iat) = eri(iat) + erepi/2
          tpvpr = tpvpr - dp
    2     continue 
        enddo
C --- Add to and print forces ---
        if (force) then
          do  ifc = 1, 3
            dt(ifc) = f(ifc,iat) + fe(ifc,iat) + fnou(ifc,iat)
     .              + fnom(ifc,iat) + dv(ifc)
          enddo
          fabs = dsqrt(dt(1)**2 + dt(2)**2 + dt(3)**2)
          if (fabs .gt. fmax) then
            fmax = fabs
            iatmx = iat
          endif
          call r8tos8(dclabl(ic),clabli)
          call dpscop(w(opos),pos,3,3*iat-2,1,1d0)
          if (ipr .gt. 20 .and. ipr .le. 30)
     .      print 101, iat,clabli,(pos(ifc),ifc=1,3),dt
          if (ipr .gt. 30) then
            if (UL .or. TBU) then
              if (lov) then
                print 103,iat,clabli,(pos(ifc),ifc=1,3),
     .                  (f(ifc,iat),ifc=1,3),(fe(ifc,iat),ifc=1,3),
     .                  (fnou(ifc,iat),ifc=1,3),(fnom(ifc,iat),ifc=1,3),
     .                  dv,dt
              else
                print 104,iat,clabli,(pos(ifc),ifc=1,3),
     .                  (f(ifc,iat),ifc=1,3),(fe(ifc,iat),ifc=1,3),
     .                  dv,dt
              endif
            else
              print 100, iat,clabli,(pos(ifc),ifc=1,3),
     .                   (f(ifc,iat),ifc=1,3),dv,dt
            endif
          endif
          call dcopy(3,dt,1,f(1,iat),1)
        endif
      enddo
      if (force .and. fmax .gt. 1d-6 .and. ipr .gt. 20) then
        call r8tos8(dclabl(w(oipc+iatmx-1)),clablj)
        call awrit2('%x= %;6,6d on atom %i ('//
     .    clablj//'%a)',outs,80,0,fmax,iatmx)
        call awrit0('%N Maximum force'//outs,' ',-80,lgunit(1))
      endif
      print *, ' '

  100 format(/' Forces on atom ',i4,4x,'Species: ',a4/
     .  '  Coordinates:',3f14.8/
     .  '  From bands :',3f14.8/'  From pairs :',3f14.8/
     .  '  Total      :',3f14.8)
  103 format(/' Forces on atom ',i4,4x,'Species: ',a4/
     .  '  Coordinates        :',3f14.8/
     .  '  From bands         :',3f14.8/
     .  '  From e''stx         :',3f14.8/
     .  '  From overlap (U)   :',3f14.8/
     .  '  From overlap (Mad) :',3f14.8/
     .  '  From pairs         :',3f14.8/
     .  '  Total              :',3f14.8)
  104 format(/' Forces on atom ',i4,4x,'Species: ',a4/
     .  '  Coordinates        :',3f14.8/
     .  '  From bands         :',3f14.8/
     .  '  From e''stx         :',3f14.8/
     .  '  From pairs         :',3f14.8/
     .  '  Total              :',3f14.8)
  102 format(/'   Site',16x,'pos',30x,'force')
  101 format(i4,' ',a,3f10.5,3x,3f10.5)

C --- Print atom- or class-specific charges and energies ---
      emg = 0d0
      if (trh .or. charge .or.  UL .or. TBU .or. stoner) then
        trh0(1) = 0d0
        trh0(2) = 0d0
        q(1) = 0d0
        q(2) = 0d0
        call dcopy(3,0d0,0,pdip,1)
        do  iat = 1, nbas
          call dpscop(w(opos),pos,3,3*iat-2,1,1d0)
          ic = w(oipc+iat-1)
          if (trh .or. charge .or. UL .or. TBU) then
            qtot(1) = dsum(nl,r(1,1,iat),1)
            q(1) = q(1) + qtot(1)
            qtot(2) = 0d0
            if (nsp .eq. 2) then
              qtot(2) = dsum(nl,r(1,2,iat),1)
              q(2) = q(2) + qtot(2)
            endif
C  --- Dipole moment ---
            call upack2('site relax',ssite,iat,ifrlx)
            q0 = dsum(nl,qnu(1,1,1,ic),3)
            if (nsp .eq. 2) then
              q0 = q0 + dsum(nl,qnu(1,1,2,ic),3)
            endif
            Dq = (qtot(1)+qtot(2)) - q0
            if (ifrlx(1) .eq. 1) then
              call daxpy(3,Dq,pos,1,pdip,1)
            endif
          endif
          call r8tos8(dclabl(ic),clabli)
          if ((nsp .eq. 2 .and. ipr .gt. 20)
     .      .or. ipr .gt. 30) write (*,310) iat,clabli
C  --- Mulliken charges ---
          if ( (trh .or. charge .or. UL .or. TBU)
     .      .and. ((nsp .eq. 2 .and. ipr .gt. 20)
     .            .or. ipr .gt. 30) ) then
            if (nsp .eq. 2) then
              if (nl .eq. 3)
     .        write (*,400) nl,((r(l,1,iat)+r(l,2,iat)), l = 1, nl),
     .           dsum(nl,r(1,1,iat),1)+dsum(nl,r(1,2,iat),1)
              if (nl .eq. 2)
     .        write (*,401) nl,((r(l,1,iat)+r(l,2,iat)), l = 1, nl),
     .           dsum(nl,r(1,1,iat),1)+dsum(nl,r(1,2,iat),1)
              if (nl .eq. 1)
     .        write (*,402) nl,((r(l,1,iat)+r(l,2,iat)), l = 1, nl),
     .           dsum(nl,r(1,1,iat),1)+dsum(nl,r(1,2,iat),1)
              write (*,410)    ((r(l,1,iat)-r(l,2,iat)), l = 1, nl)
              if (nl .gt. 2 .and. (UL .or. TBU)) then
                call awrit2('  d-band magnetic energy: I=%d, E_X=%d',
     .            ' ',120,i1mach(2),jh(3,iat),
     .            -0.25d0*(r(3,1,iat)-r(3,2,iat))**2*stni(ic))
              endif
            else 
              if (nl .eq. 3) 
     .        write (*,450) nl,(r(l,1,iat), l = 1, nl),
     .                      dsum(nl,r(1,1,iat),1)
              if (nl .eq. 2) 
     .        write (*,451) nl,(r(l,1,iat), l = 1, nl),
     .                      dsum(nl,r(1,1,iat),1)
              if (nl .eq. 1) 
     .        write (*,452) nl,(r(l,1,iat), l = 1, nl),
     .                      dsum(nl,r(1,1,iat),1)
            endif
          endif
          erepi = eri(iat)
C --- Band energies, pair energies ---
          trh0(1) = trh0(1) + e(iat,1)
          e(iat,1) = e(iat,1) - eatm(ic,1)
          if (nsp .eq. 2) then
            trh0(2) = trh0(2) + e(iat,2)
            e(iat,2) = e(iat,2) - eatm(ic,2)
          endif
          if (ipr .gt. 20) then
            write (*,560) pos, Dq
          endif
C --- Stoner model ---
          if (stoner .and. ipr .ge. 20) then
            amgm = amgm + mmom(iat)
            emg  = emg  + emag(iat)
            write (*,575) mmom(iat),emag(iat)
          endif
C  --- Magnetic energy and moment ---
          if (nsp .eq. 2 .and. .not. stoner) then
            amgm = q(1) - q(2)
            if (UL .or. TBU) then
              J = getavJ(nl,jh(1,iat),idxdn,ic)
              emg = emg - 0.25d0 * J * mmom(iat)**2
            endif
          endif
        enddo
      endif

C --- dipole moment ---
      if ((charge .or.  UL .or. TBU) .and. MOL) then
        call dscal(3,-2.541748d0,pdip,1)
        call dscal(3,-2.541748d0,ppdip,1)
        do  i = 1, 3
          ptot(i) = pdip(i) + ppdip(i)
        enddo
        dipole = dsqrt(ptot(1)**2+ptot(2)**2+ptot(3)**2)
        call dscal(3,1d0/dipole,ptot,1)
        dip1 = dsqrt(pdip(1)**2+pdip(2)**2+pdip(3)**2)
        if (dip1 .gt. d1mach(3)) then
          call dscal(3,1d0/dip1,pdip,1)
        endif
        dip2 = dsqrt(ppdip(1)**2+ppdip(2)**2+ppdip(3)**2)
        if (dip2 .gt. d1mach(3)) then
          call dscal(3,1d0/dip2,ppdip,1)
        endif
        print *, ' '
        print *, ' Molecular dipole moment in Debye ... (unit vector)'
        call awrit2('  From point charges:  %d (%3:1d)',' ',128,
     .              i1mach(2),dip1,pdip)
        if (dip2 .gt. d1mach(3)) then
          call awrit2('  From atomic dipoles: %d (%3:1d)',' ',128,
     .                i1mach(2),dip2,ppdip)
          call awrit2('  Total moment:        %d (%3:1d)',' ',128,
     .                i1mach(2),dipole,ptot)
        endif
      endif

  310 format(/' Atom ',i4,'   Species ',a4)
  400 format('  Charges: NL=',i1,': ',3(1x,f10.6),' (Total: ',f10.6,')')
  401 format('  Charges: NL=',i1,': ',2(1x,f10.6),' (Total: ',f10.6,')')
  402 format('  Charges: NL=',i1,': ',1(1x,f10.6),' (Total: ',f10.6,')')
  410 format(9x,'Moment: ',4(1x,f10.6))
  450 format('  Charges: NL=',i1,': ',3(1x,f10.6),' (Total: ',f10.6,')')
  451 format('  Charges: NL=',i1,': ',2(1x,f10.6),' (Total: ',f10.6,')')
  452 format('  Charges: NL=',i1,': ',1(1x,f10.6),' (Total: ',f10.6,')')
  560 format('  POS=',3f10.6,'  Dq/e=',f10.6)
  575 format('  Stoner Model:  MMOM=',f14.8,'   EMAG=',f14.8)

C --- Print total charges ---
      if (ipr .gt. 10) print *
      if (ipr .ge. 20) then
        if (nsp .eq. 2) then
          if (trh) then
            write (*,600) q(1)+q(2),amgm,trh0(1),trh0(2),trh0(1)+trh0(2)
          elseif (charge) then
            write (*,600) q(1)+q(2),amgm
          endif
        else
          if (trh) then
            if (stoner) then
              write (*,675) q(1),amgm,trh0(1)
            else
              write (*,650) q(1),trh0(1)
            endif
          elseif (charge .or. UL .or. TBU) then
            if (stoner) then
              write (*,675) q(1),amgm
            else
              write (*,650) q(1)
            endif
          endif
        endif
      endif

C --- Total energy ---
      erep = erep/2
      tpvpr = tpvpr/2
      if (UL .or. TBU) then
        if (nsp .eq. 1) then
          etot = trh0(1) + emad - eref + erep
        else
          etot = trh0(1) + trh0(2) + emad - eref + erep
        endif
      else
        etot = sumev + emad - eref + erep
      endif
      if (stoner) etot = etot + emg
      efree = etot - (1d0 - alpha)*entrpy
      if (ipr .ge. 10) then
        if (UL .or. TBU) then
          if (nsp .eq. 2) then
            write(*,700) sumev,emad,erep,eref,emg,amgm
          else
            write(*,700) sumev,emad,erep,eref
          endif
        else
          if (stoner) then
            write(*,750) sumev,erep,eref,emg
          else
            write(*,750) sumev,erep,eref
          endif
        endif
        write(*,800) etot, efree
      endif
      if (trh) then
        sumev = trh0(1) + trh0(2)
      endif

C --- Pressure ---
      if (pv) then
        tpv = thrpv + tpvpr + tpvq
      else
        tpv = tpvpr
      endif
      if (ipr .ge. 10) then
        if (pv) then
          write(*,900) tpvpr,thrpv,tpvq,tpv
        else
          write(*,900) tpvpr
        endif
      print *
      endif
      thrpv = tpv

  600 format( ' Tr[rho]         total   :  ',f16.8,
     .       /'                moment   :  ',f16.8,:,
     .       /' Tr[rho][H_0]       up   :  ',f16.8,
     .       /'                  down   :  ',f16.8,
     .       /'                 total   :  ',f16.8)
  650 format( ' Tr[rho]                 :  ',f16.8,:,
     .       /' Tr[rho][H_0]            :  ',f16.8)
  675 format( ' Tr[rho]                 :  ',f16.8,:,
     .       /' Stoner magnetic moment  :  ',f16.8,:,
     .       /' Tr[rho][H_0]            :  ',f16.8)
  700 format( ' band structure energy   :  ',f16.8,
     .       /' E_2                     :  ',f16.8,:,
     .       /' pair potential energy   :  ',f16.8,:,
     .       /' reference energy        :  ',f16.8,:,
     .       /' Stoner magnetic energy  :  ',f16.8,:,
     .       /' Magnetic moment         :  ',f16.8)
  750 format( ' band structure energy   :  ',f16.8,
     .       /' pair potential energy   :  ',f16.8,:,
     .       /' reference energy        :  ',f16.8,:,
     .       /' Stoner magnetic energy  :  ',f16.8)
  800 format( ' total energy            :  ',f16.8,:,
     .       /' free energy             :  ',f16.8)
  900 format( ' 3PV              pair   :  ',f16.8,:,
     .       /'                 bands   :  ',f16.8,
     .       /'                 charges :  ',f16.8,
     .       /'                 total   :  ',f16.8)

      end
