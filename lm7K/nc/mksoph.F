      subroutine mksoph(nl,nsp,nbas,lihdim,ipc,indxsh,sop,soph)
C- Create a vector of site-dependent s-o coupling parameters
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   lihdim:size of lower + downfolding block
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   indxsh:permutation indices ordering orbitals in downfolding order 
Ci   sop   :spin-orbit parameters by class (atomsr.f)
Co Outputs
Co   lix,mix: l's and m's for each RL channel, permuted according bto
Co   soph  :soph(1..6,js1,js2,n) =  spin orbit parameters 
Co         :sop(l,js1,js2,1..6,ic) ordered in downfolding order.
Co         :soph(7,1,1,*) l for this RL channel
Co         :soph(7,2,1,*) m for this RL channel
Cu Updates
Cu   08 Feb 03 altered soph(3..6) to accomodate B-field
Cu   14 Sep 99 discard uneeded argument lmx; use lihdim
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nsp,nbas,lihdim,ipc(nbas),indxsh(*)
      double precision sop(0:nl-1,nsp,nsp,9,*),soph(7,2,2,lihdim)
C Local parameters
      integer ibas,ic,l,m,lmr,is1,is2,js1,js2,k,i

      lmr = 0
      do  10  ibas = 1, nbas
        ic = ipc(ibas)
        do  8  l = 0, nl-1
          do  5  m = -l, l
            lmr = lmr+1
            if (indxsh(lmr) .gt. lihdim) goto 5
            i = indxsh(lmr)
            soph(7,1,1,i) = l
            soph(7,2,1,i) = m
C           if (l .eq. 0) goto 5
            do  2  is1 = 1, 2
            js1 = min(is1,nsp)
            do  2  is2 = 1, 2
              js2 = min(is2,nsp)
              if (l .eq. 0) then
                do  4  k = 4, 6
                  soph(k,is1,is2,i) = sop(l,js1,js2,k,ic)
    4           continue
              else
                do  1  k = 1, 6
                  soph(k,is1,is2,i) = sop(l,js1,js2,k,ic)
    1           continue
              endif
    2       continue
    5     continue
    8   continue
   10 continue
      end

