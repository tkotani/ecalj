      subroutine dosspl(nfilem,nbmx,nsp,nspc,nchan,n,w,nkp,wgts,eband,
     .doswt,npts,emin,emax,sw,wk,zos)
C- Density of states by BZ sampling
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nfilem, file handle for decomposition of the norm (MOMS or BAND)
Ci   nbmx, dimension of eband;
Ci   nsp=1 spin degenerate, =2 non-deg;
Ci   nspc: 2 for spin up and down coupled, otherwise 1
Ci   nchan, no. of channels;
Ci   n, w : order and broadening; nkp no. of k-points;
Ci   wgts : sampling BZ weights (from BZMESH);
Ci   eband, work array for bands; doswt, work array for weights;
Ci   npts, no of points in energy range: [emin, emax];
Ci   sw:  see output zos; wk, work array of length npts.
Co Outputs:
Co   zos : DOS (or integrated DOS for sw=T) for each spin and nchan
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nkp,nchan,nsp,nspc,nbmx,npts,n,nfilem
      double precision w,eband(nbmx),emin,emax,wk(npts),wgts(nkp),
     .zos(npts,nsp,nchan),doswt(nchan,nbmx,nspc)
      logical sw
C local variables
      integer nspx,ikp,isp,ib,ichan,meshpt,i,nevmx,nfstg,iomoms,jsp,ksp
      double precision bin,wt,bwt,e,emesh,x,d,s,w1,w2,xx

C --- Read past header in moments file ---
      rewind nfilem
      read(nfilem)

      nspx = nsp / nspc
      call dpzero(zos,npts*nsp*nchan)
      bin = npts - 1
      bin = (emax - emin) / bin
      do  6  ikp = 1, nkp
        do  6  isp = 1, nspx
          nfstg = 11
          if (iomoms(nfilem,i,nsp,nspc,nkp,i,nfstg,nspc,1,1,nbmx,nbmx,
     .    nchan,nchan,nevmx,eband,doswt,doswt,doswt,xx,xx) .lt. 0)
     .    call rx('DOSSPL: failed to read moments')

C --- Sum over bands ---
          do  5  ib = 1, nevmx
            e = eband(ib)

C    ...  Integrated DOS for each mesh point in energy window ...
            do  4  meshpt = 1, npts
              x = (emin - e) / w
              call delstp(n,x,d,s,xx)
              w1 = 1d0 - s
              emesh = emin + (meshpt-1)*bin
              x = (emesh - e) / w
              call delstp(n,x,d,s,xx)
              w2 = 1d0 - s
              wt = abs(wgts(ikp))*(w2 - w1) / nsp
              if (dabs(wt) .lt. 1d-10) then
                if (emesh .gt. e) goto 5
                goto 4
              endif

              do  3  jsp = 1, nspc
C       ... ksp is isp for uncoupled spins, and jsp for coupled spins
                ksp = max(jsp,isp)

C       ... Loop over channels ...
                do  3  ichan = 1, nchan
                  bwt = wt*doswt(ichan,ib,jsp)
                  zos(meshpt,ksp,ichan) = zos(meshpt,ksp,ichan) + bwt
    3         continue
    4       continue
    5     continue
    6 continue

      if (sw) return

C --- DOS from finite difference of NOS ---
      bin = 2d0 * bin
      do  11  isp  = 1, nsp
        do  11  ichan = 1, nchan
          do  10  i = 2, npts - 1
            wk(i) = (zos(i+1,isp,ichan) - zos(i-1,isp,ichan)) / bin
   10     continue
          wk(1)    = wk(2)
          wk(npts) = wk(npts-1)
          call dcopy(npts,wk,1,zos(1,isp,ichan),1)
   11 continue
      end

