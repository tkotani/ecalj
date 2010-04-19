      subroutine tbfrce(mode,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp,
     .  lmx,ipc,nrclas,idxdn,indxsh,iwk,ldim,nev,zval,ikp,bk,wt,wts,
     .  norder,width,metal,efermi,nsite,iax,npr,xyzfrz,eband,z,h,h0,dh,
     .  dhcf,vso,hso,s,ds,dscf,ip,charge,rl,rn,trh,pv,force,lso,lov,
     .  cryf,ocryf,xm0,wk,wk2,wkcf,sumev,entrpy,f,thrpv,e,rho,rholm,
     .  rhoc,rhon,drhosl)
C- TB band forces, energies, charges, and pressure for one k-point
C-----------------------------------------------------------------------
Ci Inputs
Ci   In mode 0 (sampling or insulators only) the Fermi energy or no. of
Ci   states is input; in mode 1 the k-point weights are input in wts;
Ci   lmol: molecule (cluster) branch
Ci   lgamma: gamma point only (hamiltonian real); must be TRUE if lmol
Ci   plat,nbas,nl,lmx,ipc,nrclas,indxsh
Ci   nspc = 2 for coupled spins (empirical S-O)
Ci   nsp1 = 2 if S-O or TB+U (dimensions eband and BZ integration)
Ci   nsp  = 2 if TB+U or TB-L spin polarised, isp is then current spin
Ci   iwk, number of orbitals for each atom
Ci   ldim, dimension of "lower" block of hamiltonian matrix
Ci   nev, number of eigenvectors found from diagno
Ci   zval, total no. of valence electrons
Ci   ikp,bk,wt: k-point index, k-point, and k-point degeneracy weight
Ci   wts, band weights for all k-points
Ci   norder, width, BZ sampling parameters (metals only); metal,efermi;
Ci   nsite, total number of neighbors in all clusters;
Ci   iax, neighbor lists; npr, see tbham;
Ci   eband, bands this k; z, eigenvectors this k;
Ci   h, h0, dh(x,y,z,r), real space hamiltonian and derivatives;
Ci   dhcf(x,y,z,r), real space crystal field derivatives
Ci   vso,hso: table of spin-orbit parameters and the hamiltonian
Ci   s, ds(x,y,z,r), real space overlap and derivatives;
Ci   dscf(x,y,z,r), real space overlap crystal field derivatives
Ci   if ip >= 1 print timings
Ci   Switches: charge, return s,p,d charges on each atom;
Ci             rl, return on-site c*_RL c_RL' in rhoc (see Outputs)
Ci             rn, return on-site c*_RL S c_RL' in rhon (see Outputs)
Ci             trh, return band energy on each atom;
Ci             pv, calculate 3PV;
Ci             force, calculate forces;
Ci             lso, include spin-orbit interactions
Ci             lov, include overlap matrix (non-orthogonal TB)
Ci   cryf, true if crystal field terms in Hamiltonian
Ci   ocryf, true if overlap crystal field terms
Ci   xm0,wk,wk2,wkcf: work arrays
Co Outputs
Co   sumev, sum occupied levels if lmol (i.e., a molecule)
Co   entrpy, entropy term (actually TS)
Co   f(3,nbas) forces on each atom in the basis accumulated for this k;
Co   thrpv is 3PV for this k;
Co   e(nbas,nsp1) band energy for each atom for this k
Co   rho(nl,2,nbas) s,p,d charges for each atom 
Co                                       accumulated for this k, spin
Co   rholm(nl,2,nbas) {lm} Mulliken charges for each atom 
Co                                       accumulated for this k, spin
Co                         non spin pol for now
Co   rhoc(nl**2,nl**2,nbas) s,p,d eigenvector products  c*_RL c_RL'
Co                                       accumulated for this k and spin
Co                                       and summed over spin
Co   rhon(nl**2,nl**2,nbas,nsp=2) occupation numbers for each atom 
Co                                       accumulated for this k and spin
Co   drhosl(ldim,ldim,3) c*_RL dS/dR c_RL' + cc
Co                                       accumulated for this k
Cr Remarks
Cr   Forces, rho, rho[ln] are symmetrized after k-point sum, see symfor
Cr   symr and symrtb
Cu Updates
Cu   04 Jun 08 (ATP) new molecule (cluster) mode
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer mode,nbas,nl,nspc,nsp,isp,nsp1,ldim,nev,ikp,norder,
     .        nsite,ip
      integer lmx(1),ipc(1),nrclas(1),indxsh(1),iwk(nbas),iax(1),
     .  npr(1),idxdn(0:nl-1,1)
      double precision zval,wt,width,efermi,thrpv,sumev,entrpy
      double precision plat(1),bk(3),wts(ldim,nsp,1),eband(ldim),
     .  z(ldim,ldim,2),h(nl**4*nsite*nspc**2*nsp),
     .  h0(nl**4*nsite*nspc**2*nsp),
     .  dh(nl**4*nsite*nspc**2,4),
     .  dhcf(nl**4*nsite*nspc**2,4),vso(1),hso(1),
     .  s(nl**4*nsite*nspc**2*nsp),
     .  ds(nl**4*nsite*nspc**2,4),dscf(nl**4*nsite*nspc**2,4),
     .  wk(ldim,ldim,2),wk2(ldim,ldim,2),wkcf(nbas),f(3,nbas),
     .  e(nbas,nsp1),rho(0:nl-1,2,nbas),rholm(nl**2,2,nbas),
     .  xm0(ldim*nspc),drhosl(ldim,ldim,3),
     .  rhon(0:nl**2-1,0:nl**2-1,nbas,2),
     .  rhoc(0:nl**2-1,0:nl**2-1,nbas)
      logical lmol,lgamma,metal,charge,rl,rn,trh,pv,force,lso,lov,cryf,
     .  ocryf,xyzfrz(3)
C Heap:
      integer w(1)
      common /w/ w
C Local variables
      integer oicl,oiwk,odhk,odsk
      integer nstate,nst0,ntry,itry,n,n2,l,m,i,j,ib,id,ibas,indx,iprint,
     .        lsp,i1mach,is,k,ilm
      double precision dr,zv0,e1,e2,dosef,cv,ddot,dsum
      character*80 outs
      logical Fls,T,cmdopt
      data  Fls /.false./ T /.true./

      call tcn('tbfrce')
      n = ldim
      n2 = n*n
      lsp = ldim / nspc
      call rxx(lgamma.and.nspc.eq.2,
     .  'TBFRCE not set up for S-O and gamma point only')
      call rxx((.not.trh).and.lov.and.(rn.or.rl),
     .  'TBFRCE: for s-c TB with overlap restart with TRH=T in ctrl')
      call defi(oiwk, nbas)
      if (lgamma .and. iprint() .ge. 30) then
        if ((iprint() .gt. 30 .and. ldim .lt. 10)
     .       .or. iprint() .gt. 60) then
          print *, 'TBFRCE: real eigenvectors ..'
          do  i = 1, ldim
            write(*,200) (z(i,j,1),j=1,ldim)
          enddo
  200     format (1028f8.4)
        endif
      endif

C --- Set k-point weights if insulator or mode 0 ---
      nstate = nev
      if (.not. metal .and. .not. lmol) then
        nstate = (zval + 1d-4) / (3 - nspc)
        call dcopy(nstate,wt,0,wts(1,isp,ikp),1)
        entrpy = 0d0
      elseif (mode .eq. 0) then
        if (lmol) then
          ntry = 100
          e1 = eband(1)
          e2 = eband(ldim)
          nst0 = (zval + 1d-4) / (3 - nspc)
          efermi = eband(nst0)
          do  itry = 1, ntry
            call pshprt(0)
            call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,
     .                  T,sumev,wts(1,isp,ikp),zv0,entrpy,dosef,cv)
            call popprt
            if (zv0 .gt. zval) e2 = efermi
            if (zv0 .le. zval) e1 = efermi
C           ... could use qtol here:
            if (dabs(zval - zv0) .lt. 1d-12) goto 1
            efermi = 0.5d0*(e1 + e2)
          enddo
          if (iprint() .gt.10) then
            print *, 'TBFRCE: ***warning*** cannot find HOMO/LUMO'
          endif
    1     continue 
          if (iprint() .ge. 30) then
            call awrit1(' TBFRCE: locate molecule ''Fermi'' energy '//
     .                  '... %i tries,',' ',128,i1mach(2),itry)
            if (norder .ge. 0) then
              call awrit6(' N=%i, W=%d, E_F=%d, sumev=%d, entropy term:'
     .                    //' %d, %d electrons',' ',256,i1mach(2),
     .                    norder,width,efermi,sumev,entrpy,zv0)
            else
              call awrit5(' T=%dK, E_F=%d, sumev=%d, TS=%d,'
     .                    //' %d electrons',' ',256,i1mach(2),
     .                    0.1579d6*width,efermi,sumev,entrpy,zv0)
            endif
            if (iprint() .gt. 30) then
              call splwts(1,nstate,nstate,1,wt,eband,norder,width,
     .                    efermi,T,sumev,wts(1,isp,ikp),zv0,entrpy,
     .                    dosef,cv)
            endif
          endif
        else
          call pshprt(0)
          call splwts(1,nstate,nstate,1,wt,eband,norder,width,efermi,
     .                metal,sumev,wts(1,isp,ikp),zv0,entrpy,dosef,cv)
          call popprt
        endif
      endif

C --- Get z_RL dS_RLR'L'/dR z_R'L' ---
      call tcn('make density matrix')
      if (lov .and. (rl .or. rn) .and. force
     .    .and. .not. cmdopt('--nodr',6,0,outs)) then
        call tcn('make drhosl')
        do  k = 1, 3
          call pshprt(0)
          call tbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,
     .         nsite,iax,npr,ds(1,k),vso,hso,Fls,ldim,wk2,w(oiwk))
          call popprt
          do  i = 1, ldim
            do  j = 1, ldim
              do  is = 1, nstate
                if (lgamma) then
                 drhosl(i,j,k) = drhosl(i,j,k)
     .          + z(i,is,1) * wk2(i,j,1) * z(j,is,1) * wts(is,isp,ikp)
                else
                 drhosl(i,j,k) = drhosl(i,j,k)
     .              + (
     .                  z(i,is,1) *
     .              ( wk2(i,j,1) * z(j,is,1) - wk2(i,j,2) * z(j,is,2) )
     .               +  z(i,is,2) *
     .              ( wk2(i,j,2) * z(j,is,1) + wk2(i,j,1) * z(j,is,2) )
     .                 ) * wts(is,isp,ikp)
                endif
              enddo
            enddo
          enddo
        enddo
        if (iprint() .gt. 80) then
          print *,   ' TBFRCE: drhosl:'
          do  k = 1, 3
            if (k. eq. 1) print *, '         x :'
            if (k. eq. 2) print *, '         y :'
            if (k. eq. 3) print *, '         z :'
            do  i = 1, ldim
              write(*,400) (drhosl(i,j,k),j=1,ldim)
  400         format (1028f10.6)
            enddo
          enddo
        endif
        call tcx('make drhosl')
      endif
      if (rl .or. rn) then
C --- Accumulate z*_RL z_RL' in rhoc
        call tcn('mkrhol')
        call mkrhol(lgamma,ldim,nl,nsp,isp,nstate,nbas,ipc,wts,ikp,
     .              z,z,idxdn,rhoc)
        call tcx('mkrhol')
      endif
C --- Get Mulliken s, p, and d charges summed over bands ---
      if (.not. (trh .or. charge)) goto 3

C --- Get S_R'L'RL z_RL ---
      if (lov) then
        call tcn('make SZ')
        call pshprt(0)
        call tbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,
     .              nsite,iax,npr,s,vso,hso,Fls,ldim,wk2,w(oiwk))
        call popprt
        if (lgamma) then
          call dmpy(wk2,n,1,z,n,1,wk,n,1,n,nstate,n)
        else
          call zmpy(wk2,n,1,n2,z,n,1,n2,wk,n,1,n2,n,nstate,n)
        endif
        call tcx('make SZ')
      endif

      if (rn) then
C --- Accumulate z*_RL S z_RL' in rhon (this spin) ---
        call tcn('make ZSZ')
        if (lov) then
          call mkrhol(lgamma,ldim,nl,nsp,isp,nstate,nbas,ipc,wts,ikp,
     .      z,wk,idxdn,rhon(0,0,1,isp))
        else
          call mkrhol(lgamma,ldim,nl,nsp,isp,nstate,nbas,ipc,wts,ikp,
     .      z,z,idxdn,rhon(0,0,1,isp))
        endif
        if (iprint() .gt. 60) then
          do  ibas = 1, nbas
            call awrit3('isp=%i ikp=%i ibas=%i; z S z :',' ',
     .        120,i1mach(2),isp,ikp,ibas)
            do  i = 0, nl**2-1
              write (*,'(5x,9f10.6)') (rhon(i,j,ibas,isp),j=0,nl**2-1)
            enddo
          enddo
        endif
        call tcx('make ZSZ')
      else
C --- Collect Mulliken charges ---
        call tcn('Mulliken charges')
        do ib = 1, nstate
          call dpzero(xm0,nbas*nspc*nl**2)
C --- Make orbital-decomposed weights ---
          do  i = 1, ldim
            if (lov) then
              if (lgamma) then
                xm0(i) = z(i,ib,1)*wk(i,ib,1)
              else
                xm0(i) = z(i,ib,1)*wk(i,ib,1) + z(i,ib,2)*wk(i,ib,2)
              endif
            else
              if (lgamma) then
                xm0(i) = z(i,ib,1)**2
              else
                xm0(i) = z(i,ib,1)**2 + z(i,ib,2)**2
              endif
            endif
          enddo
          id = 0
          do  ibas = 1, nbas
            ilm = 0
            do  l = 0, nl-1
              do  m = -l, l
                ilm = ilm + 1
                id = id + 1
                j = indxsh(id)
C ...           Exclude "higher" orbitals ---
                if (j .gt. lsp) goto 2
                dr = xm0(j)*wts(ib,isp,ikp)
                rholm(ilm,isp,ibas) = rholm(ilm,isp,ibas) + dr
                rho(l,isp,ibas) = rho(l,isp,ibas) + dr
                if (nspc .eq. 2) then
                  dr = xm0(j+lsp)*wts(ib,isp,ikp)
                  rho(l,2,ibas) = rho(l,2,ibas) + dr
                endif
    2           continue
              enddo
            enddo
          enddo
        enddo
        call tcx('Mulliken charges')
      endif

C --- Get band energy sum for each atom, use H_in, not H ---
      if (.not. trh) goto 3
      call tcn('Tr[rho H]')
      call defdc(odhk, n2)
      call pshprt(0)
      call tbloch(lgamma,bk,nl,nspc,nsp,isp,nbas,plat,lmx,ipc,indxsh,
     .            nsite,iax,npr,h0,vso,hso,Fls,ldim,w(odhk),w(oiwk))
      call popprt
      if (lgamma) then
        call dmpy(w(odhk),n,1,z,n,1,wk,n,1,n,nstate,n)
      else
        call zmpy(w(odhk),n,1,n2,z,n,1,n2,wk,n,1,n2,n,nstate,n)
      endif
      call rlse(odhk)
      do  ib = 1, nstate
        indx = 1
        do  ibas = 1, nbas
          if (lgamma) then
            dr = ddot(iwk(ibas),z(indx,ib,1),1,wk(indx,ib,1),1)
          else
            dr = ddot(iwk(ibas),z(indx,ib,1),1,wk(indx,ib,1),1)
     .         + ddot(iwk(ibas),z(indx,ib,2),1,wk(indx,ib,2),1)
          endif
          e(ibas,isp) = e(ibas,isp) + dr*wts(ib,isp,ikp)
          if (nspc .eq. 2) then
            if (lgamma) then
            dr = ddot(iwk(ibas),z(indx+lsp,ib,1),1,wk(indx+lsp,ib,1),1)
            else
            dr = ddot(iwk(ibas),z(indx+lsp,ib,1),1,wk(indx+lsp,ib,1),1)
     .         + ddot(iwk(ibas),z(indx+lsp,ib,2),1,wk(indx+lsp,ib,2),1)
            endif
            e(ibas,2) = e(ibas,2) + dr*wts(ib,isp,ikp)
          endif
          indx = indx + iwk(ibas)
        enddo
      enddo
      call tcx('Tr[rho H]')

C --- Get 3PV ---
    3 continue
      call tcx('make density matrix')
      if (pv) then
        call tcn('3PV')
        if (lgamma) then
          call defdr(odhk, n2)
        else
          call defdc(odhk, n2)
        endif
        call pshprt(0)
        call tbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,
     .       indxsh,nsite,iax,npr,dh(1,4),vso,hso,Fls,ldim,
     .       w(odhk),w(oiwk))
        call popprt
        if (lgamma) then
          call dmpy(w(odhk),n,1,z,n,1,wk,n,1,n,nstate,n)
        else
          call zmpy(w(odhk),n,1,n2,z,n,1,n2,wk,n,1,n2,n,nstate,n)
        endif
C --- Printout ---
        if ((ikp .eq. 1 .and. iprint() .gt. 70) .or. (iprint()
     .    .ge. 110)) then
          call defi(oicl, nbas*nspc)
          call icopy(nbas,ipc,1,w(oicl),1)
          if (nspc .eq. 2) call icopy(nbas,ipc,1,w(oicl+nbas),1)
          print *,'radial component'
          call prmsk2(lmx,lmx,w(oicl),2,ldim,ldim,w(odhk),'?','dH',1,
     .      ldim,indxsh)
          call rlse(oicl)
        endif
        call rlse(odhk)
C --- Overlap contribution ---
        if (lov) then
          if (lgamma) then
            call defdr(odsk, n2)
          else
            call defdc(odsk, n2)
          endif
          call pshprt(0)
          call tbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,
     .      nsite,iax,npr,ds(1,4),vso,hso,Fls,ldim,w(odsk),w(oiwk))
          call popprt
          if (lgamma) then
            call dmpy(w(odsk),n,1,z,n,1,wk2,n,1,n,nstate,n)
          else
            call zmpy(w(odsk),n,1,n2,z,n,1,n2,wk2,n,1,n2,n,nstate,n)
          endif
C --- Printout ---
          if ((ikp .eq. 1 .and. iprint() .gt. 70) .or. (iprint()
     .      .ge. 110)) then
            call defi(oicl, nbas*nspc)
            call icopy(nbas,ipc,1,w(oicl),1)
            if (nspc .eq. 2) call icopy(nbas,ipc,1,w(oicl+nbas),1)
            print *,'radial component'
            call prmsk2(lmx,lmx,w(oicl),2,ldim,ldim,w(odsk),'?','dS',1,
     .        ldim,indxsh)
            call rlse(oicl)
          endif
          call rlse(odsk)
        endif
C --- Sum over bands ---
        do  ib = 1, nstate
          if (lov) then
            call daxpy(n,-1d0*eband(ib),wk2(1,ib,1),1,wk(1,ib,1),1)
            if (.not. lgamma) then
              call daxpy(n,-1d0*eband(ib),wk2(1,ib,2),1,wk(1,ib,2),1)
            endif
          endif
          if (lgamma) then
            dr = ddot(n,z(1,ib,1),1,wk(1,ib,1),1)
          else
            dr = ddot(n,z(1,ib,1),1,wk(1,ib,1),1)
     .         + ddot(n,z(1,ib,2),1,wk(1,ib,2),1)
          endif
          if (cryf) then
            call dpzero(wkcf,nbas)
            call cffor(nbas,nl,nspc,nsite,ldim,ib,lmx,ipc,iax,npr,
     .        indxsh,dhcf(1,4),z,w(oiwk),pv,wkcf)
            dr = dr + dsum(nbas,wkcf,1)
          endif
          if (ocryf) then
            call dpzero(wkcf,nbas)
            call cffor(nbas,nl,nspc,nsite,ldim,ib,lmx,ipc,iax,npr,
     .        indxsh,dscf(1,4),z,w(oiwk),pv,wkcf)
            dr = dr - eband(ib)*dsum(nbas,wkcf,1)
          endif
          thrpv = thrpv - dr*wts(ib,isp,ikp)
        enddo
        call tcx('3PV')
      endif
      
C --- Loop over x, y, z : get Hellmann-Feynman forces ---
      if (.not. force) goto 5
      call tcn('HF forces')
      do  i = 1, 3
        if (xyzfrz(i)) goto 4
        if (lgamma) then
          call defdr(odhk, n2)
        else
          call defdc(odhk, n2)
        endif
        call pshprt(0)
        call tbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,
     .    nsite,iax,npr,dh(1,i),vso,hso,Fls,ldim,w(odhk),w(oiwk))
        call popprt
        if (lgamma) then
          call dmpy(w(odhk),n,1,z,n,1,wk,n,1,n,nstate,n)
        else
          call zmpy(w(odhk),n,1,n2,z,n,1,n2,wk,n,1,n2,n,nstate,n)
        endif
C --- Printout ---
        if ((ikp .eq. 1 .and. iprint() .gt. 70) .or. (iprint()
     .    .ge. 110) .and. .not. lgamma) then
          if (i .eq. 1) then
            call defi(oicl, nbas*nspc)
            call icopy(nbas,ipc,1,w(oicl),1)
            if (nspc .eq. 2) call icopy(nbas,ipc,1,w(oicl+nbas),1)
            print *,'x-component'
          endif
          if (i .eq. 2) print *,'y-component'
          if (i .eq. 3) print *,'z-component'
          call prmsk2(lmx,lmx,w(oicl),2,ldim,ldim,w(odhk),'?','dH',1,
     .      ldim,indxsh)
          if (i .eq. 3) call rlse(oicl)
        endif
        call rlse(odhk)
C --- Overlap contribution ---
        if (lov) then
          if (lgamma) then
            call defdr(odsk, n2)
          else
            call defdc(odsk, n2)
          endif
          call pshprt(0)
          call tbloch(lgamma,bk,nl,nspc,1,1,nbas,plat,lmx,ipc,indxsh,
     .      nsite,iax,npr,ds(1,i),vso,hso,Fls,ldim,w(odsk),w(oiwk))
          call popprt
          if (lgamma) then
            call dmpy(w(odsk),n,1,z,n,1,wk2,n,1,n,nstate,n)
          else
            call zmpy(w(odsk),n,1,n2,z,n,1,n2,wk2,n,1,n2,n,nstate,n)
          endif
C --- Printout ---
          if ((ikp .eq. 1 .and. iprint() .gt. 70) .or. (iprint()
     .      .ge. 110) .and. .not. lgamma) then
            if (i .eq. 1) then
              call defi(oicl, nbas*nspc)
              call icopy(nbas,ipc,1,w(oicl),1)
              if (nspc .eq. 2) call icopy(nbas,ipc,1,w(oicl+nbas),1)
              print *,'x-component'
            endif
            if (i .eq. 2) print *,'y-component'
            if (i .eq. 3) print *,'z-component'
            call prmsk2(lmx,lmx,w(oicl),2,ldim,ldim,w(odsk),'?','dS',1,
     .        ldim,indxsh)
            if (i .eq. 3) call rlse(oicl)
          endif
          call rlse(odsk)
        endif
C --- Sum over bands ---
        do  ib = 1, nstate
          indx = 1
          if (lov) then
            call daxpy(n,-1d0*eband(ib),wk2(1,ib,1),1,wk(1,ib,1),1)
            if (.not. lgamma) then
              call daxpy(n,-1d0*eband(ib),wk2(1,ib,2),1,wk(1,ib,2),1)
            endif
          endif
          if (cryf .or. ocryf) then
            call dpzero(wkcf,nbas)
            if (ocryf) then
              call cffor(nbas,nl,nspc,nsite,ldim,ib,lmx,ipc,iax,npr,
     .          indxsh,dscf(1,i),z,w(oiwk),Fls,wkcf)
              call dscal(nbas,-1d0*eband(ib),wkcf,1)
            endif
            if (cryf) call cffor(nbas,nl,nspc,nsite,ldim,ib,lmx,ipc,
     .        iax,npr,indxsh,dhcf(1,i),z,w(oiwk),Fls,wkcf)
            call dscal(nbas,-1d0*wts(ib,isp,ikp),wkcf,1)
            if (iprint() .ge. 110) then
              if (i .eq. 1) print *,'Crystal field forces: x-components'
              if (i .eq. 2) print *,'Crystal field forces: y-components'
              if (i .eq. 3) print *,'Crystal field forces: z-components'
              write(*,500) ib,(wkcf(ibas),ibas = 1, nbas)
  500         format(2x,'Band=',i4,(/4x,9f8.4:))
            endif
          endif
          do  ibas = 1, nbas
            if (cryf .or. ocryf) f(i,ibas) = f(i,ibas) + wkcf(ibas)
            if (lgamma) then
              dr = ddot(iwk(ibas),z(indx,ib,1),1,wk(indx,ib,1),1)
            else
              dr = ddot(iwk(ibas),z(indx,ib,1),1,wk(indx,ib,1),1)
     .           + ddot(iwk(ibas),z(indx,ib,2),1,wk(indx,ib,2),1)
            endif
            f(i,ibas) = f(i,ibas) - 2d0*dr*wts(ib,isp,ikp)
            if (nspc .eq. 2) then
             dr = ddot(iwk(ibas),z(indx+lsp,ib,1),1,wk(indx+lsp,ib,1),1)
     .          + ddot(iwk(ibas),z(indx+lsp,ib,2),1,wk(indx+lsp,ib,2),1)
             f(i,ibas) = f(i,ibas) - 2d0*dr*wts(ib,isp,ikp)
            endif
            indx = indx + iwk(ibas)
          enddo
        enddo
    4   continue 
      enddo
      call tcx('HF forces')

    5 call rlse(oiwk)

      call tcx('tbfrce')
      end
      subroutine mkrhol(lg,ldim,nl,nsp,isp,nstate,nbas,ipc,wts,ikp,
     .                  zl,zr,idxdn,rhol)
C- Make zl*_RL zr_RL'
C ----------------------------------------------------------------------
Ci Inputs: 
Ci   lg: gamma point only: z real
Ci   eigenvectors and k-point weights from BZ integration
Co Outputs:
Co   rhol
Cr Remarks
Cr   Assembles rho_lm (see notes) from the eigenvectors zl (left) and
Cr   zr (right). In a non-orthogonal basis then on entry zr could
Cr   be Sz. Otherwise zl and zr can point to the same address.
Cr   If zl and zr are complex then mkrhol makes (zl* zr + zl zr*)/2,
Cr   in other words, the real part of zl times zr.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      logical lg
      integer ikp,ldim,nl,nsp,isp,nstate,nbas,ipc(1),idxdn(0:nl-1,1)
      double precision wts(ldim,nsp,1),zl(ldim,ldim,2),zr(ldim,ldim,2),
     .                 rhol(0:nl**2-1,0:nl**2-1,nbas)
C Local Variables
      integer ib,is,ic,ptr,ptr0,m1,m2,m3,m4
      double precision r

      ptr = 1
      do  ib = 1, nbas
        ic = ipc(ib)
        ptr0 = ptr
C ... s block
        if (idxdn(0,ic) .eq. 1) then
          do  is = 1, nstate
            if (lg) then
              r =  zl(ptr0,is,1)*zr(ptr0,is,1) * wts(is,isp,ikp)
            else
              r = (zl(ptr0,is,1)*zr(ptr0,is,1)
     .          +  zl(ptr0,is,2)*zr(ptr0,is,2)) * wts(is,isp,ikp)
            endif
            rhol(0,0,ib) = rhol(0,0,ib) + r
          enddo
          ptr0 = ptr0 + 1
        endif
C ... p block
        if (nl .gt. 1 .and. idxdn(1,ic) .eq. 1) then
          do  is = 1, nstate
            do  m1 = 0, 2
              do  m2 = 0, 2
                if (lg) then
                  r = zl(ptr0+m1,is,1)*zr(ptr0+m2,is,1)
     .              * wts(is,isp,ikp)
                else
                  r = (zl(ptr0+m1,is,1)*zr(ptr0+m2,is,1)
     .              +  zl(ptr0+m1,is,2)*zr(ptr0+m2,is,2))
     .              * wts(is,isp,ikp)
                endif
                rhol(m1+1,m2+1,ib) = rhol(m1+1,m2+1,ib) + r
              enddo
            enddo
          enddo
          ptr0 = ptr0 + 3
        endif
C ... d block
        if (nl .gt. 2 .and. idxdn(2,ic) .eq. 1) then
          do  is = 1, nstate
            do  m1 = 0, 4
              do  m2 = 0, 4
                if (lg) then
                  r = zl(ptr0+m1,is,1)*zr(ptr0+m2,is,1)
     .              * wts(is,isp,ikp)
                else
                  r = (zl(ptr0+m1,is,1)*zr(ptr0+m2,is,1)
     .              +  zl(ptr0+m1,is,2)*zr(ptr0+m2,is,2))
     .              * wts(is,isp,ikp)
                endif
                rhol(m1+4,m2+4,ib) = rhol(m1+4,m2+4,ib) + r
              enddo
            enddo
          enddo
          ptr0 = ptr0 + 5
        endif
C ... sp block
        if (nl .gt. 1 .and. idxdn(0,ic) .eq. 1 .and. idxdn(1,ic) .eq. 1)
     .    then
          do  is = 1, nstate
            do  m1 = 1, 3
              if (lg) then
                r = zl(ptr,is,1)*zr(ptr+m1,is,1) * wts(is,isp,ikp)
              else
                r = (zl(ptr,is,1)*zr(ptr+m1,is,1)
     .            +  zl(ptr,is,2)*zr(ptr+m1,is,2)) * wts(is,isp,ikp)
              endif
              rhol(0,m1,ib) = rhol(0,m1,ib) + r
              rhol(m1,0,ib) = rhol(m1,0,ib) + r
            enddo
          enddo
        endif
C ... sd block
        if (nl .gt. 2 .and. idxdn(0,ic) .eq. 1 .and. idxdn(2,ic) .eq. 1)
     .    then
          do  is = 1, nstate
            do  m1 = 4, 8
              m2 = m1
              if (idxdn(1,ic) .ne. 1) m2 = m1 - 3
              if (lg) then
                r = zl(ptr,is,1)*zr(ptr+m2,is,1) * wts(is,isp,ikp)
              else
                r = (zl(ptr,is,1)*zr(ptr+m2,is,1)
     .            +  zl(ptr,is,2)*zr(ptr+m2,is,2)) * wts(is,isp,ikp)
              endif
              rhol(0,m1,ib) = rhol(0,m1,ib) + r
              rhol(m1,0,ib) = rhol(m1,0,ib) + r
            enddo
          enddo
        endif
C ... pd block
        if (nl .gt. 2 .and. idxdn(1,ic) .eq. 1 .and. idxdn(2,ic) .eq. 1)
     .    then
          do  is = 1, nstate
            do  m1 = 1, 3
              do  m2 = 4, 8
                if (idxdn(0,ic) .eq. 1) then
                  m3 = m1
                  m4 = m2
                else
                  m3 = m1 - 1
                  m4 = m2 - 1
                endif
                if (lg) then
                  r = zl(ptr+m3,is,1)*zr(ptr+m4,is,1)
     .              * wts(is,isp,ikp)
                else
                  r = (zl(ptr+m3,is,1)*zr(ptr+m4,is,1)
     .              +  zl(ptr+m3,is,2)*zr(ptr+m4,is,2))
     .              * wts(is,isp,ikp)
                endif
                rhol(m1,m2,ib) = rhol(m1,m2,ib) + r
                rhol(m2,m1,ib) = rhol(m2,m1,ib) + r
              enddo
            enddo
          enddo
        endif
        ptr = ptr0
      enddo
      call rxx(ptr .ne. ldim+1,'Bug in mkrhol')
      end
