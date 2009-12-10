      subroutine mkwttb(nfileb,nbas,nl,nsp,nclass,ipc,nrclas,idxdn,
     .  indxsh,nev,ldim,mull,imfil,nchan,nlm,natm,ifi,ng,istab,g,gd,lov,
     .  zll,sll,pdos,bc,bcos,pwk,bwk,eband,doswt)
C- Calculate Mulliken DOS weights for tight-binding at 1 k-point
C ----------------------------------------------------------------------
Ci Inputs
Ci   nfileb : file handle for BAND
Ci   nbas,nl,nsp,nclass,ipc,nrclas,idxdn,indxsh
Ci   nev : number of eigenvectors found in secmtb
Ci   ldim : dimension of "lower" block of hamiltonian matrix
Ci   mull : Mulliken DOS switch, determines type of DOS:
Ci           0, Partial DOS,     resolved by {class, L, spin}
Ci           1, Partial DOS,     resolved by {atom, L, M, and spin}
Ci          10, Bond charge DOS, resolved by class pairs
Ci          11, Bond charge DOS, resolved by atom pairs
Ci          20, Bond charge DOS, resolved by {class, L} pairs
Ci          21, Bond charge DOS, resolved by {atom, L} pairs
Ci          30, Bond charge DOS, resolved by {class, L, M, spin} pairs
Ci          31, Bond charge DOS, resolved by {atom, L, M, spin} pairs
Ci   imfil : 0, Accumulate all DOS channels (mull < 20)
Ci           1, only accumulate channels from file MULL (mull < 20)
Ci   nchan : number of DOS channels, see remarks
Ci   nlm,natm : dimensions for DOS work arrays:
Ci          mull =  0,        nlm  = nl*nsp
Ci          mull =  1,        nlm  = nl*nl*nsp
Ci          mull = 10 or 11,  nlm  = 1
Ci          mull = 20 or 21,  nlm  = nl
Ci          mull = 30 or 31,  nlm  = nl*nl*nsp
Ci          mod(mull,10) = 0, natm = nclass
Ci          mod(mull,10) = 1, natm = nbas
Ci   ifi : file handle for MULL, see remarks
Ci   ng : number of group operations
Ci   istab : atom istab(i,ig) is transformed into atom i by operation ig
Ci   g : 3x3 rotation matrices (Cartesian)
Ci   gd : transformation table for d-orbitals, see symtbd
Ci   lov : true if have overlap matrix (non-orthogonal TB)
Ci   zll : eigenvectors this k
Ci   sll : Bloch-transformed overlap matrix this k
Ci   pdos : work array for partial DOS
Ci   bc : work array for bond charge DOS
Ci   bcos : work array for on-site bond charge DOS
Ci   pwk,bwk : work arrays used in symmetrization
Ci   eband : bands for this k, written to file BAND
Co Outputs
Co   doswt : Mulliken DOS for this k, must be initialized to zero,
Co           written to file BAND
Cr Remarks
Cr   For (mull=0,1; imfil=0) all the partial DOS channels are generated
Cr     automatically.  All of the angular momentum channels for a given
Cr     atom appear consecutively in the DOS file (see below for specific
Cr     ordering).  The atom ordering follows the class indices for
Cr     (mull=0) and the basis indices for (mull=1).  In the case of
Cr     (nsp=2) all of the spin up channels for a given atom or class
Cr     appear consecutively followed by all of the spin down channels
Cr     for the same atom or class.
Cr
Cr  For (mull=0,1; imfil=1) the partial DOS channels are specified in
Cr    the file MULL using the class indices (mull=0) or the basis
Cr    indices (mull=1) and the angular momentum indices below.  The
Cr    format of the lines in this free format file is:
Cr
Cr      <ang mom ind 1> <atom ind 1>
Cr      ...
Cr
Cr    In the case of (nsp=2) the choice of spin up or down is specified
Cr    in the file MULL by choosing the appropriate lm index (see below).
Cr
Cr  For (mull=10,11; imfil=0) the bond charge is automatically
Cr    generated for all atom pairs.  The ordering of the channels is
Cr    done according to the class indices for (mull=10) and the basis
Cr    indices for (mull=11).  The ordering is as follows:
Cr
Cr         (mull=11)            (mull=10)
Cr
Cr      Atom 1   Atom 2      Atom 1   Atom 2    
Cr      ------   ------      ------   ------
Cr         1        1           1        1    On-site contribution
Cr         1        2           1        1    All others (if existent)
Cr         1        3           1        2
Cr       ...      ...         ...      ...       
Cr         2        2           2        2    On-site contribution
Cr         2        3           2        2    All others (if existent)
Cr         2        4           2        3
Cr       ...      ...         ...      ...
Cr
Cr  For (mull=10,11, imfil=1) the bond charge DOS channels are specified
Cr    in the file MULL using the class indices (mull=10) or the basis
Cr    indices (mull=11).  The format of the lines in this free format
Cr    file is:
Cr
Cr      <atom ind 1> <atom ind 2> <ionsit>
Cr      ...
Cr
Cr     The switch (ionsit) must always be present and must be equal to
Cr     zero unless [mod(mull,10)=0] and the two atom class indices are
Cr     the same.  In this case (ionsit=1) refers to the on-site
Cr     contribution and (ionsit=0) refers to all other contributions
Cr     (if they exist).
Cr
Cr   For (mull=20,21) and (mull=30,31) the atoms are indexed according
Cr     to the class index (mull=20,30) or the basis index (mull=21,31).
Cr     Ths DOS channels are specified in the file MULL using these atom
Cr     indices and the angular momentum indices indicated below.  The
Cr     format of the lines in this file is (free format):
Cr
Cr       <lm ind 1> <atom ind 1> <lm ind 2> <atom ind 2> <ionsit>
Cr       ...
Cr
Cr     The switch (ionsit) is the same as for (mull=10,11; imfil=1)
Cr
Cr     In the case of (nsp=2) the spin combination is specified in the
Cr     file MULL by choosing the appropriate lm indices (see below).
Cr
Cr   The angular momentum channels for each atom or class for
Cr     (mull=0,20,21) are indexed as follows:
Cr       lm = 1 : s
Cr       lm = 2 : p (summed over m = -1, 0, 1)
Cr       lm = 3 : d (summed over m = -2, -1, 0, 1, 2)
Cr       ...
Cr
Cr     In the case of (nsp=2) all of the spin up lm indices appear
Cr     consecutively followed by all of the spin down lm indices
Cr     (i.e. if nl=2 then lm=4 corresponds to spin down p).
Cr
Cr   The angular momentum channels for each atom or class for
Cr     (mull=1,30,31) are indexed as follows:
Cr       lm = 1     : s
Cr       lm = 2,3,4 : p_x, p_y, p_z
Cr       lm = 5-9   : d_xy, d_yz, d_zx, d_(x^2-y^2), d_(3z^2-r^2)
Cr       ...
Cr
Cr     In the case of (nsp=2) all of the spin up lm indices appear
Cr     consecutively followed by all of the spin down lm indices
Cr     (i.e. if nl=2 then lm=6 corresponds to spin down p_x).
Cr
Cr   The DOS channels are symmetrized if needed [ng > 0 and
Cr     mod(mull,10) = 1] for (mull=1,11,21).  No symmetrization is done
Cr     for (mull=30,31) and for (mull=1, L > 2) and therefore the full
Cr     BZ should be used (no special points) or symmetry-equivalent
Cr     channels should be averaged.
Cu Updates
Cu   10 Apr 02 Bug fix --- nfstg now set correctly
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nfileb,nbas,nl,nsp,nclass,ldim,mull,imfil,nchan,nlm,natm,
     .  ifi,ng,nev
      integer ipc(nbas),nrclas(nclass),idxdn(0:nl-1,1),indxsh(*),
     .  istab(nbas,ng)
      double precision g(3,3,ng),gd(5,5,ng),zll(ldim,ldim,2),
     .  sll(ldim,ldim,2),pdos(nlm,natm),bc(nlm,nlm,natm,natm),
     .  bcos(nlm,nlm,natm),pwk(nlm,natm),bwk(nlm,nlm,natm,natm),
     .  eband(nev),doswt(nchan,nev)
      logical lov
C Local variables
      integer ib,id,i,j,ic,l,m,lm,lsp,mull0,mull1,ichan,i2,id2,ic2,
     .  lm2,l2,m2,j2,nls,ion,ig,jc,jc2,ll,fopno,iomoms,nfstg
      double precision zj1,zj2,zj3,zj4,wgt,wgt2,wgt3,wgt4
      logical lmsum,dsym,lonsit

C --- Initialization ---
      mull0 = mod(mull,10)
      mull1 = mull / 10
      call rxx(mull0 .ne. 0 .and. mull0 .ne. 1,
     .  'MKWTTB: bad ones digit in mull switch')
      nls = nlm / nsp
      if (mull1 .eq. 1) nls = 0
      lsp = ldim / nsp
      lmsum = .not. (mull .eq. 1 .or. mull1 .eq. 3)
C     dsym = ng .gt. 0 .and. mull0 .eq. 1
      dsym = ng .gt. 1 .and. mull0 .eq. 1

C --- Sum over bands ---
      do  240  ib = 1, nev
        if (mull1 .eq. 0) then
          call dpzero(pdos,nlm*natm)
        else
          call dpzero(bc,  nlm*nlm*natm*natm)
          if (mull0 .eq. 0) call dpzero(bcos,nlm*nlm*natm)
        endif

C --- First loop over orbitals ---
        id = 0
        do  60  i = 1, nbas
          ic = i
          if (mull0 .eq. 0) ic = ipc(i)
          lm = 0
          do  50  l = 0, nl-1
            do  40  m = 1, 2*l+1
              lm = lm + 1
              if (lmsum) lm = l+1
              if (mull1 .eq. 1) lm = 1
              id = id + 1
              j = indxsh(id)
C --- Exclude "higher" orbitals ---
              if (j .gt. lsp) goto 40
              zj1 = zll(j,ib,1)
              zj2 = zll(j,ib,2)
              if (nsp .eq. 2) then
                zj3 = zll(j+lsp,ib,1)
                zj4 = zll(j+lsp,ib,2)
              endif

C --- Skip second loop if no overlap and partial DOS ---
              if (.not. lov .and. mull1 .eq. 0) then
                wgt = zj1*zj1 + zj2*zj2
                pdos(lm,ic) = pdos(lm,ic) + wgt
                if (nsp .eq. 2) then
                  wgt2 = zj3*zj3 + zj4*zj4
                  pdos(lm+nls,ic) = pdos(lm+nls,ic) + wgt2
                endif
              else

C --- Second loop over orbitals ---
              id2 = 0
              do  30  i2 = 1, nbas
                lonsit = i .eq. i2 .and. mull0 .eq. 0
                ic2 = i2
                if (mull0 .eq. 0) ic2 = ipc(i2)
                lm2 = 0
                do  20  l2 = 0, nl-1
                  do  10  m2 = 1, 2*l2+1
                    lm2 = lm2 + 1
                    if (lmsum) lm2 = l2+1
                    if (mull1 .eq. 1) lm2 = 1
                    id2 = id2 + 1
                    j2 = indxsh(id2)
C --- Exclude "higher" orbitals ---
                    if (j2 .gt. lsp) goto 10

C --- Accumulate Mulliken DOS weights ---
                    if (lov) then
                     wgt =
     .                (zll(j2,ib,1)*zj1 + zll(j2,ib,2)*zj2)*sll(j2,j,1)+
     .                (zll(j2,ib,2)*zj1 - zll(j2,ib,1)*zj2)*sll(j2,j,2)
                      if (nsp .eq. 2) then
                        wgt2 =(zll(j2+lsp,ib,1)*zj3
     .                       + zll(j2+lsp,ib,2)*zj4)*sll(j2+lsp,j+lsp,1)
     .                       +(zll(j2+lsp,ib,2)*zj3
     .                       - zll(j2+lsp,ib,1)*zj4)*sll(j2+lsp,j+lsp,2)
                        wgt3 = (zll(j2+lsp,ib,1)*zj1
     .                       +  zll(j2+lsp,ib,2)*zj2)*sll(j2+lsp,j,1)
     .                       + (zll(j2+lsp,ib,2)*zj1
     .                       -  zll(j2+lsp,ib,1)*zj2)*sll(j2+lsp,j,2)
                        wgt4 = (zll(j2,ib,1)*zj3
     .                       +  zll(j2,ib,2)*zj4)*sll(j2,j+lsp,1)
     .                       + (zll(j2,ib,2)*zj3
     .                       -  zll(j2,ib,1)*zj4)*sll(j2,j+lsp,2)
                      endif
                    else
                      wgt = zll(j2,ib,1)*zj1 + zll(j2,ib,2)*zj2
                      if (nsp .eq. 2) then
                        wgt2 = zll(j2+lsp,ib,1)*zj3
     .                       + zll(j2+lsp,ib,2)*zj4
                        wgt3 = zll(j2+lsp,ib,1)*zj1 
     .                       + zll(j2+lsp,ib,2)*zj2
                        wgt4 = zll(j2,ib,1)*zj3 + zll(j2,ib,2)*zj4
                      endif
                    endif
                    if (mull1 .eq. 0) then
C --- Partial DOS ---
                      pdos(lm,ic) = pdos(lm,ic) + wgt
                      if (nsp .eq. 2) then
                        pdos(lm,ic) = pdos(lm,ic) + wgt3
                        pdos(lm+nls,ic) = pdos(lm+nls,ic) + wgt2 + wgt4
                      endif
                    else
                      if (lonsit) then
C --- On-site contributions ---
                        bcos(lm,lm2,ic) = bcos(lm,lm2,ic) + wgt
                        if (nsp .eq. 2) then
                          bcos(lm+nls,lm2,ic) = bcos(lm+nls,lm2,ic)
     .                                        + wgt4
                          bcos(lm,lm2+nls,ic) = bcos(lm,lm2+nls,ic)
     .                                        + wgt3
                          bcos(lm+nls,lm2+nls,ic) = wgt2
     .                                      + bcos(lm+nls,lm2+nls,ic)
                        endif
                      else
C --- Inter-site contributions ---
                        bc(lm,lm2,ic,ic2) = bc(lm,lm2,ic,ic2) + wgt
                        if (nsp .eq. 2) then
                          bc(lm+nls,lm2,ic,ic2) = bc(lm+nls,lm2,ic,ic2)
     .                                          + wgt4
                          bc(lm,lm2+nls,ic,ic2) = bc(lm,lm2+nls,ic,ic2)
     .                                          + wgt3
                          bc(lm+nls,lm2+nls,ic,ic2) = wgt2
     .                                      + bc(lm+nls,lm2+nls,ic,ic2)
                        endif
                      endif
                    endif
   10             continue
   20           continue
   30         continue
              endif
   40       continue
   50     continue
   60   continue

C --- Symmetrize DOS if needed ---
        if (dsym .and. mull0 .eq. 1 .and. (mull1 .eq. 0
     .    .or. mull1 .eq. 1 .or. mull1 .eq. 2)) then
          if (mull1 .eq. 0) then
            call dpcopy(pdos,pwk,1,nlm*natm,1d0/ng)
            call dpzero(pdos,nlm*natm)
          else
            call dpcopy(bc,bwk,1,nlm*nlm*natm*natm,1d0/ng)
            call dpzero(bc,nlm*nlm*natm*natm)
          endif
C --- Loop over all symmetry operations and atoms ---
          do  130  ig = 1, ng
            do  120  ic = 1, natm
              jc = istab(ic,ig)
              if (mull1 .eq. 0) then
C --- Symmetrize partial DOS ---
                do  90  lm = 1, nlm
                  l = ll(lm)
                  m = lm - l*l
                  if (lm .gt. nls) then
                    l = ll(lm-nls)
                    m = lm - nls - l*l
                  endif
                  if (idxdn(l,ipc(ic)) .gt. 1) goto 90
                  if (l .eq. 0) then
C --- Symmetrize L = 0 (sites only) ---
                    pdos(lm,ic) = pdos(lm,ic) + pwk(lm,jc)
                  elseif (l .eq. 1) then
C --- Symmetrize L = 1 ---
                    lm2 = 1
                    if (lm .gt. nls) lm2 = 1 + nls
                    do  70  m2 = 1, 3
                      lm2 = lm2 + 1
                      pdos(lm,ic) = pdos(lm,ic)
     .                            + g(m,m2,ig)*g(m,m2,ig)*pwk(lm2,jc)
   70               continue
                  elseif (l .eq. 2) then
C --- Symmetrize L = 2 ---
                    lm2 = 4
                    if (lm .gt. nls) lm2 = 4 + nls
                    do  80  m2 = 1, 5
                      lm2 = lm2 + 1
                      pdos(lm,ic) = pdos(lm,ic)
     .                            + gd(m,m2,ig)*gd(m,m2,ig)*pwk(lm2,jc)
   80               continue
                  elseif (l .gt. 2) then
C --- No symmetrization for L > 2 ---
                    pdos(lm,ic) = pdos(lm,ic) + pwk(lm,ic)
                  endif
   90           continue
              else
C --- Symmetrize bond charge DOS ---
                do  110  ic2 = 1, natm
                  jc2 = istab(ic2,ig)
C --- Loop over all orbital pairs ---
                  do  100  lm = 1, nlm
                    do  100  lm2 = 1, nlm
C --- Symmetrize sites only ---
                      bc(lm2,lm,ic2,ic) = bc(lm2,lm,ic2,ic)
     .                                  + bwk(lm2,lm,jc2,jc)
  100             continue
  110           continue
              endif
  120       continue
  130     continue
        endif

C --- Accumulate Mulliken DOS weights ---
        ichan = 0
        if (mull1 .eq. 0) then
C --- Partial DOS ---
          if (imfil .eq. 0) then
C --- Generate all channels ---
            do  150  i = 1, natm
              do  140  lm = 1, nlm
                if (mull0 .eq. 0) then
                  l = lm - 1
                  if (lm .gt. nls) l = lm - nls - 1
                  ic = i
                else
                  l = ll(lm)
                  if (lm .gt. nls) l = ll(lm-nls)
                  ic = ipc(i)
                endif
                if (idxdn(l,ic) .gt. 1) goto 140
                ichan = ichan + 1
                doswt(ichan,ib) = pdos(lm,i)
  140         continue
  150       continue
          else
C --- Read channels from file MULL ---
            ifi = fopno('MULL')
            rewind ifi

C --- Start loop to read indices of DOS channels ---
  160       read(ifi,*,err=170,end=170) lm,i
            ichan = ichan + 1

C --- Check for bad atom index ---
            if (i .gt. natm .or. i .le. 0) then
              write(*,500) lm,i
  500         format(/' MKWTTB: lm=',i2,'   i=',i3)
              call rx('MKWTTB: bad atom index')
            endif

C --- Check for bad angular momentum index ---
            if (mull0 .eq. 0) then
              l = lm - 1
              if (lm .gt. nls) l = lm - nls - 1
              ic = i
            else
              l = ll(lm)
              if (lm .gt. nls) l = ll(lm-nls)
              ic = ipc(i)
            endif
            if (lm .gt. nlm .or. lm .le. 0 .or. idxdn(l,ic) .gt. 1) then
              write(*,500) lm,i
              call rx('MKWTTB: bad angular momentum index')
            endif

C --- Save appropriate DOS channel ---
            doswt(ichan,ib) = pdos(lm,i)
            goto 160
          endif
  170     continue

        elseif (mull1 .eq. 1) then
C --- Bond charge DOS by atom pairs ---
          if (imfil .eq. 0) then
C --- Generate all channels ---
            do  190  i = 1, natm
              do  180  i2 = i, natm
                ichan = ichan + 1
                if (i .ne. i2) then
                  doswt(ichan,ib) = bc(1,1,i,i2) + bc(1,1,i2,i)
                else
                  if (mull0 .eq. 0) then
                    doswt(ichan,ib) = bcos(1,1,i)
                    if (nrclas(i) .gt. 1) then
                      ichan = ichan + 1
                      doswt(ichan,ib) = bc(1,1,i,i)
                    endif
                  else
                    doswt(ichan,ib) = bc(1,1,i,i)
                  endif
                endif
  180         continue
  190       continue
          else
C --- Read channels from file MULL ---
            ifi = fopno('MULL')
            rewind ifi

C --- Start loop to read indices of DOS channels ---
  200       read(ifi,*,err=210,end=210) i,i2,ion
            ichan = ichan + 1

C --- Check for bad atom indices ---
            if (i .gt. natm .or. i .le. 0 .or. i2 .gt. natm
     .        .or. i2 .le. 0) then
              write(*,510) i,i2,ion
  510         format(/' MKWTTB: i=',i3,'   i2=',i3,'   ionsit=',i2)
              call rx('MKWTTB: bad atom index')
            endif

C --- Check for bad <ionsit> switch ---
            if (mull0 .eq. 0) then
              ic = i
              ic2 = i2
            else
              ic = ipc(i)
              ic2 = ipc(i2)
            endif
            if (((mull0 .eq. 1 .or. ic .ne. ic2) .and. ion .ne. 0)
     .        .or. (mull0 .eq. 0 .and. ic .eq. ic2
     .        .and. nrclas(ic) .eq. 1 .and. ion .eq. 0)) then
              write(*,510) i,i2,ion
              call rx('MKWTTB: bad <ionsit> switch')
            endif

C --- Save appropriate DOS channel ---
            if (i .ne. i2) then
              doswt(ichan,ib) = bc(1,1,i,i2) + bc(1,1,i2,i)
            else
              if (mull0 .eq. 0 .and. ion .ne. 0) then
                doswt(ichan,ib) = bcos(1,1,i)
              else
                doswt(ichan,ib) = bc(1,1,i,i)
              endif
            endif
            goto 200
          endif
  210     continue

        elseif (mull1 .eq. 2 .or. mull1 .eq. 3) then
C --- Bond charge DOS by orbital pairs ---
          ifi = fopno('MULL')
          rewind ifi

C --- Get indices for DOS channels from file MULL ---
  220     read(ifi,*,err=230,end=230) lm,ic,lm2,ic2,ion
          ichan = ichan + 1

C --- Check for bad atom indices ---
          if (ic .gt. natm .or. ic .le. 0 .or. ic2 .gt. natm
     .      .or. ic2 .le. 0) then
            write(*,520) lm,ic,lm2,ic2,ion
  520       format(/' MKWTTB: lm=',i2,'   ic=',i4,'   lm2=',i2,
     .        '   ic2=',i4,'   ion=',i2)
            call rx('MKWTTB: bad atom index')
          endif

C --- Check for bad angular momentum indices ---
          if (mull0 .eq. 0) then
            i = ic
            i2 = ic2
          else
            i = ipc(ic)
            i2 = ipc(ic2)
          endif
          if (mull1 .eq. 2) then
            l = lm - 1
            l2 = lm2 - 1
            if (lm .gt. nls) l = lm - nls - 1
            if (lm2 .gt. nls) l2 = lm2 - nls - 1
          elseif (mull1 .eq. 3) then
            l = ll(lm)
            l2 = ll(lm2)
            if (lm .gt. nls) l = ll(lm-nls)
            if (lm2 .gt. nls) l2 = ll(lm2-nls)
          endif
          if (lm .gt. nlm .or. lm .le. 0 .or. lm2 .gt. nlm
     .      .or. lm2 .le. 0 .or. idxdn(l,i) .gt. 1
     .      .or. idxdn(l2,i2) .gt. 1) then
            write(*,520) lm,ic,lm2,ic2,ion
            call rx('MKWTTB: bad angular momentum index')
          endif

C --- Check for bad on-site switch ---
          if (((mull0 .eq. 1 .or. i .ne. i2) .and. ion .ne. 0)
     .      .or. (mull0 .eq. 0 .and. i .eq. i2 .and. nrclas(i) .eq. 1
     .      .and. ion .eq. 0)) then
            write(*,520) lm,ic,lm2,ic2,ion
            call rx('MKWTTB: bad <ionsit> switch')
          endif

C --- Save appropriate DOS channel ---
          if (lm .ne. lm2 .or. ic .ne. ic2) then
            if (ic .ne. ic2 .or. mull0 .eq. 1 .or. ion .eq. 0) then
              doswt(ichan,ib) =
     .          bc(lm,lm2,ic,ic2) + bc(lm2,lm,ic2,ic)
            else
              doswt(ichan,ib) =
     .          bcos(lm,lm2,ic) + bcos(lm2,lm,ic)
            endif
          else
            if (mull0 .eq. 0 .and. ion .ne. 0) then
              doswt(ichan,ib) = bcos(lm,lm,ic)
            else
              doswt(ichan,ib) = bc(lm,lm,ic,ic)
            endif
          endif
          goto 220
  230     continue
        else
          call rx('MKWTTB: bad tens digit in mull switch')
        endif
        call rxx(ichan .ne. nchan,'MKWTTB: ichan ne nchan')
  240 continue

C --- Write dimensions, bands, and DOS weights to file BAND ---
      nfstg = 11
      j = iomoms(-nfileb,nl,1,1,2,ldim,nfstg,1,1,1,nev,nev,nchan,nchan,
     .  nev,eband,doswt,doswt,doswt,wgt,wgt)
      end
