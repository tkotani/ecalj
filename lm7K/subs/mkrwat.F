      subroutine mkrwat(aavw,bas,nbas,iax,ntab,ips,nl,npr,plat,hcr,
     .  deltr,rwats)
C- Make Watson sphere radii
C ----------------------------------------------------------------------
Ci Inputs:
Ci   aavw  :scale for plat,bas that puts it in units of hcr and deltr
Ci         :(normally alat/avw)
Ci   bas   :basis vectors (units of alat)
Ci   nbas  :size of basis
Ci   deltr :radius of Watson sphere = rmx + deltr
Ci   iax   :structure-constant information relating to s (see STRRS)
Ci   ips   :the jth atom belongs to species ips(j)
Ci   jbas  :center of current cluster
Ci   nl    :number of 's
Ci   npr   :number of neighbors within rmaxs
Ci   plat  :primitive lattice vectors (units of alat)
Ci   hcr   :hard core screening radius, scaled by 1/avw
Co Outputs:
Co   :screening constants for Watson-sphere orbitals
Co   adotom:(kappa*avw)^2-derivative of 
Co   rwats :Watson sphere radius, scaled by 1/avw
Cr Remarks:
Cr   First, the distance from central atom to outermost atom is found.
Cr   Then the screening constants of the Watson orbitals are calculated.
Cr   N^al_omega(rwats)=N^0_omega(rwats)-*J^0_omega(rwats)=!0
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer niax,nbas,ntab(nbas),ips(*),nl,npr
      parameter (niax=10)
      integer iax(niax,npr)
      double precision bas(3,*),deltr,plat(3,3),rwats(nbas),aavw,
     .  hcr(nl,*)
C Local variables:
      integer ibas,is,ipr,jbas
      double precision drr2,dr(3),rmx,rout

C ... For each site, make Watson sphere radii
      do  10  jbas = 1, nbas

C ... Watson sphere radius = range of the cluster + deltr
      rmx = 0d0
      do  20  ipr = ntab(jbas)+1, ntab(jbas+1)
        ibas = iax(2,ipr)
        is   = ips(ibas)
        rout = aavw*dsqrt(drr2(plat,bas(1,ibas),bas(1,jbas),
     .    -iax(3,ipr),-iax(4,ipr),-iax(5,ipr),dr)) + hcr(1,is)
        rmx = max(rmx,rout)
   20 continue
      rwats(jbas) = rmx + deltr
   10 continue

      end

