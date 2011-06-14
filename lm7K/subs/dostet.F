      subroutine dostet(nbmx,nsp,nspx,nevmx,nchan,n1,n2,n3,ntet,
     .idtet,eband,doswt,npts,emin,emax,lidos,wk,zos)
C- Density of states to third order by tetrahedron method
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbmx, first dim. of eband;
Ci   nsp=1 spin degenerate, =2 non-deg;
Ci   nspx: 1 for spin up and down coupled, otherwise nsp
Ci   nevmx, no. of bands;
Ci   nchan, no. of DOS channels; n1,n2,n3;
Ci   ntet, idtet, o/p from tetirr
Ci   eband, work array for bands; doswt, work array for weights;
Ci   npts, no of points in energy range: [emin, emax];
Ci   lidos :F zos = dos
Ci         :T zos = energy integral of dos
Co Outputs:
Co   zos : DOS (or integrated DOS for lidos=T) for each spin and nchan
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nchan,nsp,nspx,nbmx,npts,ntet,idtet(0:4,*),n1,n2,n3,nevmx
      double precision eband(nbmx,nspx,1),emin,emax,wk(npts),
     .zos(npts,nsp,nchan),doswt(nchan,nbmx,nsp,1)
      logical lidos
C ... Local parameters
      integer isp,ib,i,itet,ichan,iq1,iq2,iq3,iq4,nspc,jsp,ksp
      double precision bin,eigen(4),v,wt,ebot,dmin1,eg0(4)
C ... External calls
      external dcopy,dpzero,slinz

      if (npts .le. 1 .or. npts .le. 2 .and. .not. lidos) call rx1(
     .'dostet: npts(=%i) too small for DOS : require npts>2',npts)
      nspc = nsp / nspx
      call dpzero(zos,npts*nsp*nchan)
      bin = npts - 1
      bin = (emax - emin) / bin
      v = ( 3d0  -  nsp ) / ( n1 * n2 * n3 * 6d0 )

C --- Loop over tetrahedra ---
      do  5  itet = 1, ntet
        iq1 = idtet(1,itet)
        iq2 = idtet(2,itet)
        iq3 = idtet(3,itet)
        iq4 = idtet(4,itet)

C --- Loop over spins and sum over bands ---
        do  4  isp = 1, nspx
          do  3  ib = 1, nevmx
            eigen(1) = eband(ib,isp,iq1)
            eigen(2) = eband(ib,isp,iq2)
            eigen(3) = eband(ib,isp,iq3)
            eigen(4) = eband(ib,isp,iq4)
            call dcopy(4,eigen,1,eg0,1)
            ebot = dmin1(eigen(1),eigen(2),eigen(3),eigen(4))
            if (ebot .gt. emax) goto 3

            do  2  jsp = 1, nspc
C       ... ksp is isp for uncoupled spins, and jsp for coupled spins
              ksp = max(jsp,isp)

C       ... Accumulate no. states assuming constant wt from this tet
              do  2  ichan = 1, nchan
                wt = doswt(ichan,ib,ksp,iq1)
     .          + doswt(ichan,ib,ksp,iq2)
     .          + doswt(ichan,ib,ksp,iq3)
     .          + doswt(ichan,ib,ksp,iq4)
                wt = wt * idtet(0,itet) * v / 4d0
                call slinz(wt,eigen,emin,emax,zos(1,ksp,ichan),npts)
    2       continue
    3     continue
    4   continue
    5 continue

      if (lidos) return

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

