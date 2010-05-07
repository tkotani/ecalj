      subroutine makipt(nl,nbas,nsp,ipc,lihdim,indxsh,pp,ndim,pti)
C- Make inverse potential function vector, 2nd generation LMTO
C-----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   lihdim:size of lower + downfolding block
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci   pp    :potential parameters in gamma representation (atomsr.f)
Ci   ndim  :leading dimension of pti
Co Outputs
Co   pti   :inverse potential functions
Cu Updates
Cu   22 Jul 07 patch: when C is exactly enu, pti=NULLI
Cu   14 Sep 99 discard uneeded argument lmx; use lihdim
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nbas,nsp,ipc(nbas),lihdim,ndim,indxsh(ndim)
      double precision pp(6,0:nl-1,nsp,1),pti(ndim,nsp)
C Local parameters
      integer lmr,ibas,l,m,ic,isp,iprint,i
      double precision alpha,c,delta,enu,gamma,xx         
      integer NULLI
      parameter (NULLI=-99999)

      if (iprint() .ge. 100) print *,
     .  'Makipt: Inverse potential functions, RL order'
      do  4  isp = 1, nsp
        lmr = 0
        do  3  ibas = 1, nbas
          ic = ipc(ibas)
          do  2  l = 0, nl-1
            do  1  m = -l, l
              lmr = lmr + 1
              if (indxsh(lmr) .le. lihdim) then
                enu   = pp(1,l,isp,ic)
                gamma = pp(5,l,isp,ic)
                alpha = pp(6,l,isp,ic)
                delta = pp(3,l,isp,ic)
                c     = pp(2,l,isp,ic)
                if (gamma .ne. alpha) then
                  xx = 1 + (c-enu)*(gamma-alpha)/delta**2
                  delta = delta*xx
                  c = enu + (c-enu)*xx
                endif
                if (enu .eq. c) then
                  pti(indxsh(lmr),isp) = NULLI
                else
                  pti(indxsh(lmr),isp) = pp(5,l,isp,ic) +
     .              delta**2 / ( enu - c )
                endif
              endif
    1       continue
    2     continue
    3   continue
        if (iprint() .ge. 100) then
          if (nsp .eq. 2) write(*,10) isp
          write (*,20) (pti(indxsh(i),isp),i=1,ndim)
        endif
    4 continue
   10 format(1x,'Spin ',i1)
   20 format(1x,9f8.4)
      end

