      subroutine makpph(nl,nsp,nbas,lihdim,ipc,indxsh,pp,pph)
C- Create a vector of site-dependent potential parameters
C ----------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   lihdim:size of lower + downfolding block
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci          (makidx.f)
Ci   pp    :potential parameters (atomsr.f)
Co Outputs
Co   pph   :potential parameters in downfolding order (makpph.f)
Cr Remarks
Cr   Mapping of composite RL index to this vector is same as that of
Cr   eigenvectors permutated according to downfolding rules in indxsh
Cr   pph(1) : enu
Cr   pph(2) : calpha
Cr   pph(3) : sqrdel
Cr   pph(4) : palpha
Cr   pph(5) : oalp
Cu Updates
Cu   14 Sep 99 discard uneeded argument lmx; use lihdim; make both spins
Cb Bugs
Cb   makpph assumes that all m's corresponding to a given l
Cb   are folded down in the same way.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nsp,nbas,lihdim,ipc(nbas),indxsh(lihdim)
      double precision pp(6,nl,nsp,*),pph(5,lihdim,nsp)
C Local parameters
      integer ibas,ic,l,m,lmr,isp,k
      double precision oalp,gmina

      do  3  isp = 1, nsp
      lmr = 0
      do  3  ibas = 1, nbas
        ic = ipc(ibas)
        do  2  l = 0, nl-1
          k = l+1
          if (indxsh(lmr+1) .gt. lihdim) then
            lmr = lmr + 2*l+1
            goto 2
          endif
          gmina = pp(6,k,isp,ic) - pp(5,k,isp,ic)
          oalp = gmina/
     .    (pp(3,k,isp,ic)**2 + gmina*(pp(1,k,isp,ic) - pp(2,k,isp,ic)))
          do  1  m = -l, l
            lmr = lmr+1
            pph(1,indxsh(lmr),isp) = pp(1,k,isp,ic)
            pph(2,indxsh(lmr),isp) = pp(2,k,isp,ic)
            pph(3,indxsh(lmr),isp) = pp(3,k,isp,ic)
            pph(4,indxsh(lmr),isp) = pp(4,k,isp,ic)
            pph(5,indxsh(lmr),isp) = oalp
    1     continue
    2   continue
    3 continue
C     call yprm('pph',1,pph,0,5,5,lihdim*nsp)

      end

