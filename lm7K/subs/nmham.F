      subroutine nmham(job,lidim,ldim,nmto,kmto,vmtz,g,Gd,g0n,hk,ok)
C- NMTO hamiltonian and overlap
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :not used now.
Ci   lidim :leading dimension of G,Gd
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   nmto :number of NMTO energies
Ci   kmto  :NMTO energies
Ci   vmtz  :muffin-tin zero, needed for the hamiltonian
Ci   G0n   :a work array, to hold G[0,..,N]
Ci Inputs/Outputs
Cio   G    :On input, kink matrix K for each of nmto energies
Cio        :On output, G = K^-1
Cio   Gd   :On input, energy derivative Kdot of K
Cio        :On output, Gdot = K^-1 Kdot K^-1
Co Outputs
Ci   hk    :NMTO hamilton matrix
Ci   ok    :NMTO overlap matrix
Cr Remarks
Cr   NMTO Hamiltonian, overlap matrix by Lagrange interpolation
Cr   See OK Andersen, Mt. St. Odile book, section 3.2.
Cr
Cr   Following the notation in the book,
Cr          G_i = the Green matrix for energy eps_i = K_i^-1
Cr       G[i,j] = the divided difference (G_i - G_j)/(eps_i - eps_j)
Cr                Note that G[i,i] = Gdot_i = G_i Kdot_i G_i
Cr      G[0..N] = Nth divided difference of the Green matrix
Cr    gam[0..N] = Nth divided difference of contracted Green's function
Cr
Cr   The Hamiltonian and overlap matrices both expressible in terms of
Cr   the following integrals, and which Andersen showed are related to
Cr   Hermite divided differences G[[0..M]..N] (92,161) :
Cr
Cr      <gam[0..N] | gam[0..M]> = - G[[0..M]..N]
Cr
Cr   f[[0..M]..N] is a polynomial of order M+N+1 which function value
Cr   and slope at the first M points and value at all points and can be
Cr   obtained by Hermite interpolation (see Appendix).
Cr
Cr
Cr   The Overlap matrix is (see Eq. 94)
Cr     <Chi^N | Chi^N > =  G[0..N]+^-1 <gam[0..N] gam[0..N]> G[0..N]^-1
Cr                      = -G[0..N]+^-1 G[[0..M]..N] G[0..N]^-1
Cr
Cr   and the Hamiltonian matrix is (see Eqn 95)
Cr     G[0..N]+^-1 <gam[0..N] | H | gam[0..N]> G[0..N]+^-1
Cr
Cr   which we can obtain from (Eqn 93)
Cr     <gam[0..N] | (H - eps_N) | gam[0..N]> = G[[0..N-1]N]
Cr
Cr   The computationally useful formula for G[0..N] is (Eqns 73,149)
Cr      G[0..N] = sum_i=0^N f^N_i G_i
Cr   where
Cr     f^N_i = prod_j=0^N,j<>N (eps_i-eps_j)^-1
Cr
Cr   A second computationally useful formula is (Eqns 92,161)
Cr     G[[0..M]..N] =  sum_i=0^N sum_j=0^M (-G[i,j]) f^N_i f^M_j
Cr
Cr   For the hamiltonian we need
Cr     <gam[0..N] | H | gam[0..N]> = G[[0..N-1]N] + eps_N G[[0..N]]
Cr
Cr   Using the relation f^(N-1)_j = f^(N)_j (eps_j-eps_N) we obtain
Cr    <gam[0..N] | H | gam[0..N]> = -(G[[0..N-1]N] + eps_N G[[0..N]])
Cr
Cr     =  - sum_i=0^N sum_j=0^N G[i,j] f^N_i f^N_j (eps_j-eps_N)
Cr        + eps_N *
Cr        - sum_i=0^N sum_j=0^N G[i,j] f^N_i f^N_j
Cr        (first sum can j=0..N because eps_j-eps_N=0 for j=N)
Cr
Cr     =    sum_i=0^N sum_j=0^N G[i,j] f^N_i f^N_j eps_j
Cr
Cr     =    G[0..N]+^-1 <Chi^N | H | Chi^N > G[0..N]-1
Cr
Cr   Since both H and Overlap are decorated by G[0..N]-1, i.e.
Cr     G[0..N]+-1 <..> G[0..N]-1
Cr
Cr   we defer this decoration for the eigenvectors.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,lidim,ldim,nmto
      double precision kmto(nmto),vmtz
      double complex G(lidim,lidim,nmto),Gd(lidim,lidim,nmto)
      double complex hk(ldim,ldim),ok(ldim,ldim),G0n(ldim,ldim)
C ... Local parameters
      integer nkapmx,ik,jk,l2,owk,ierr,iprint,i,j
      parameter (nkapmx=8)
      double precision fn(nkapmx),fac
      double complex alpha,beta
C ... Heap
      integer w(1)
      common /w/ w

      l2 = 2*ldim**2
      alpha = 1
      beta = 0
      call dpzero(hk,l2)
      call dpzero(ok,l2)

C --- f_n = prod_{n'=0, \neq n}^{N} (kmto_n-kmto_n')^-1 ---
      call nmefac(nmto,kmto,fn)

C --- Scale kink matrices to Aleph rep'sn ---
C      call defcc(owk,ldim)
C      call naleph(nmto,ldim,G,Gd,w(owk))
C      call rlse(owk)

C --- For each energy, overwrite K with G = K^-1, Kdot with Gdot ---
      call kink2g(lidim,ldim,nmto,G,Gd)

C ------ Hamiltonian and Overlap, Lagrange interp (Eqs 94, 95) --------

C --- Assemble G[[0..N]] (Eqn 92) for overlap and corresponding for h---
      do  ik = 1, nmto
        do  jk = 1, nmto
C         G[ik,jk] * (e_ik - e_jk), or Gd[ik,ik] if ik = jk
          if (ik .eq. jk) then
            fac = fn(ik)*fn(jk)
            call zmscop(0,ldim,ldim,lidim,ldim,0,0,0,0,Gd(1,1,ik),G0n)
          else
            do  j = 1, ldim
            do  i = 1, ldim
              G0n(i,j) = G(i,j,ik) - G(i,j,jk)
            enddo
            enddo
            fac = fn(ik)*fn(jk)/(kmto(ik) - kmto(jk))
          endif
C         Overlap, Eq. 92
          call daxpy(l2,fac,G0n,1,ok,1)
C         For Hamiltonian, see Remarks above
          call daxpy(l2,fac*(kmto(jk)+vmtz),G0n,1,hk,1)
        enddo
      enddo

      if (iprint() .ge. 120) then
        call zprm('<gam gam>',2,ok,ldim,ldim,ldim)
        call zprm('<gam e gam>',2,hk,ldim,ldim,ldim)
      endif

C --- Decorate ham,overlap by G[0,..,N]^-1 ---
      if (job .eq. 1) then

C   ... G[0,..,N]^-1 (Eqns 73 and 149)
        call dpzero(g0n,l2)
        do  ik = 1, nmto
C         call daxpy(l2,fn(ik),G(1,1,ik),1,g0n,1)
          do  j = 1, ldim
          do  i = 1, ldim
            G0n(i,j) = G0n(i,j) + fn(ik) * G(i,j,ik)
          enddo
          enddo
        enddo
        if (iprint() .ge. 120) then
          call zprm('G[0..N]',2,G0n,ldim,ldim,ldim)
        endif

C   ... G[0,..,N]^-1 from G[0,..,N]
        call defcc(owk,l2)
        call zqinv('h',G0n,ldim,0,ldim,w(owk),ldim,ierr)
        if (ierr .ne. 0) call rx('nmham:  G[0,..,N] is singular')
        call rlse(owk)

C   ... Scale Hamiltonian, Overlap by G[0,..,N]^-1
        call zgemm('C','N',ldim,ldim,ldim,alpha,G0n,ldim,
     .    ok,ldim,beta,w(owk),ldim)
        call zgemm('N','N',ldim,ldim,ldim,alpha,w(owk),ldim,
     .    G0n,ldim,beta,ok,ldim)

        call defcc(owk,l2)
        call zgemm('C','N',ldim,ldim,ldim,alpha,G0n,ldim,
     .    hk,ldim,beta,w(owk),ldim)
        call zgemm('N','N',ldim,ldim,ldim,alpha,w(owk),ldim,
     .    G0n,ldim,beta,hk,ldim)

        if (iprint() .ge. 110) then
          call zprm('Overlap',2,ok,ldim,ldim,ldim)
          call zprm('Hamiltonian',2,hk,ldim,ldim,ldim)
        endif
        call rlse(owk)
      endif

      end

      subroutine naleph(nmto,ldim,k,kd,diagkd)
C- Rescale kink matrices to aleph representation
C ----------------------------------------------------------------------
Ci Inputs
Ci   nmto
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Co Outputs
Co   diagkd:sqrt(diagonal part of kd)
Ci   k     :diagkdd by
Ci   kd
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nmto,ldim
      double complex k(ldim,ldim,nmto),kd(ldim,ldim,nmto),diagkd(ldim)
C ... Local parameters
      integer i,j,ik

      do  i = 1, ldim
        diagkd(i) = 1/sqrt(kd(i,i,1))
      enddo

      do  ik = 1, nmto
      do  i = 1, ldim
      do  j = 1, ldim
        kd(i,j,ik) = diagkd(i) * kd(i,j,ik) * diagkd(j)
        k (i,j,ik) = diagkd(i) * k (i,j,ik) * diagkd(j)
      enddo
      enddo
      enddo

      end

