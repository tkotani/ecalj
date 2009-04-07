      subroutine kink2g(lidim,ldim,nmto,k,kd)
C- Overwrites Kink matrix K with G=K^-1 and Kdot with Gdot
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   nmto  :number of NMTO energies
Ci Inputs/Outputs
Cio   k    :On input, kink matrix K for each of nmto energies
Cio        :On output, G = K^-1
Cio   kd   :On input, energy derivative Kdot of K
Cio        :On output, Gdot = K^-1 Kdot K^-1
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lidim,ldim,nmto
      double complex k(lidim,lidim,nmto),kd(lidim,lidim,nmto)
C ... Local parameters
      integer ierr,ik,owk,iprint
      double complex alpha,beta
C ... Heap
      integer w(1)
      common /w/ w

      alpha = 1
      beta = 0

C --- f_n = prod_{n'=0, \neq n}^{N} (kaph_n-kaph_n')^-1
C     call nmefac(nmto,kaph,fn,gnm)

C --- For each energy, overwrite K with G = K^-1, Kdot with Gdot ---
      call defcc(owk,ldim**2)
      do  ik = 1, nmto

        if (iprint() .ge. 120) then
          call zprm('K(ik)',2,k(1,1,ik),lidim,ldim,ldim)
          call zprm('Kd(ik)',2,kd(1,1,ik),lidim,ldim,ldim)
        endif

        call zqinv('h',k(1,1,ik),lidim,0,ldim,w(owk),ldim,ierr)
        call rxx(ierr.ne.0,'kink2g: Kink matrix singular')
C       Gdot = -G Kdot G
        call zgemm('N','N',ldim,ldim,ldim,alpha,kd(1,1,ik),lidim,
     .    k(1,1,ik),lidim,beta,w(owk),ldim)
        call dscal(2*ldim**2,-1d0,w(owk),1)
        call zgemm('N','N',ldim,ldim,ldim,alpha,k(1,1,ik),lidim,w(owk),
     .    ldim,beta,kd(1,1,ik),lidim)

        if (iprint() .ge. 120) then
          call zprm('G(ik)',2,k(1,1,ik),lidim,ldim,ldim)
          call zprm('Gd(ik)',2,kd(1,1,ik),lidim,ldim,ldim)
        endif
      enddo
      call rlse(owk)

      end

