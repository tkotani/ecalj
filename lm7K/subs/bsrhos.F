      subroutine bsrhos(nbas,nl,ipc,rhos,nrhos,qnu,pp,sop,eula,neul,
     .  bxc,bsite,nbf,lihdim,indxsh,mode,bsigr)
C- Double-counting term <B.sigma.rho>
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   rhos  :spin density matrix (for mode=0), by class
Ci         :rhos should be hermitian in spin space, but may not be owing
Ci         :to energy integration errors in the complex plane.
Ci         :bsrhos uses a symmetrized form to minimize errors.
Ci   qnu   :moments (for mode=1)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l- and m-independent,
Ci         :nl if Euler are l-dependent and m-independent
Ci         :nl**2 if Euler are l- and m-dependent
Ci   bsite :magnetic field by site
Ci   nbf   :1 if bsite is l- and m-independent,
Ci         :nl if bsite is l-dependent and m-independent
Ci         :nl**2 if bsite is l- and m-dependent
Ci         :99 : flag that there exists no external field
Ci   mode  :0, use spin density matrix to make moments along T
Ci         :1, use qnus to make moments along qnu
Co Outputs
Cr Remarks
Cr   Definition of rho in terms of M: (standard definition of sigma)
Cr      rho = M . sigma/2
Cr   Pauli matrices sigma: (NB: normal def is twice this)
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr
Cr   Given rho, M can be obtain from:
Cr     M_x =  2 Re(rho21) = Re (rho12+rho21)
Cr     M_y =  2 Im(rho21) = Im (rho21-rho12)
Cr            M_z =  (rho11)-(rho22)
Cr
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be quite hermitian; e.g. when rhos is generated
Cr   by a Green's function technique.
Cr
Cr   Double counting term is
Cr     Tr <(B.sigma)(rho)>
Cr
Cr   Input B.sigma is
Cr                1 (Bz   Bx-iBy)
Cr    B.sigma =   - (           )
Cr                2 (Bx+iBy  -Bz)
Cr
Cr   Then (CHECK wrong factors of 2 in both b and sigma)
Cr      Tr <(B.sigma)(rho)>
Cr      = 1/2 (Bz(rho11-rho22) + (Bx-iBy) rho21 + (Bx+iBy) rho12)
Cr      = Bz(rho11-rho22)/2 + (Bx-iBy)(Mx+iMy)/4 + (Bx+iBy)(Mx-iMy)/4
Cr      = Bz(rho11-rho22)/2 + Bx Mx / 2 + By My /2
Cr      = 1/2 (B . M)
Cr   This formula can be computed either with the moments qnu
Cr   or from the spin-density matrix.
Cu Updates
Cu   07 Apr 04 First created
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer nbas,nrhos,nl,nsp,ipc(nbas),mode,lihdim,indxsh(lihdim),
     .  neul,nbf
      parameter (nsp=2)
      double precision rhos(2,3,nrhos,2,2,1),eula(nbas,neul,3),bxc(3,*),
     .  qnu(3,0:nl-1,2,*),pp(6,0:nl-1,nsp,*),bsite(nbas,nbf,3),bsigr(2),
     .  sop(0:nl-1,nsp,nsp,9,*)
C Local variables
      end

