      subroutine addoff(nl,nbas,lmx,ipc,indxsh,ldim,pot0,sk,hk,indH)
C- Add self consistent off-site terms to H_k
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nl,nbas,lmx,ipc,indxsh,ldim,
Ci   indH: work length nbas
Ci   pot0: monopole potentials, from tbesel
Ci   pot0: monopole potentials, from tbesel
Ci   sk,hk: Overlap and H in k-representation (real part only)
Co Outputs:
Co   hk: off site terms added (eq 7.81 Mike's book)
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nl,nbas,indxsh(nl**2*nbas),lmx(1),ipc(nbas),ldim,
     .        indH(nbas)
      double precision pot0(nbas),hk(ldim,ldim),sk(ldim,ldim)
C Local Variables
      integer nlm,i,j,ib,jb,ixi,ixj,nlmi,nlmj,lmi,lmj
C Intrinsic function
      nlm(i) = (1+lmx(i))**2

C --- Form table of indices that mark offset of Rth block in H(k) ---
      indH(1) = 0
      do  ib = 2, nbas
        indH(ib) = indH(ib-1) + nlm(ipc(ib-1))
      enddo

      do  ib = 1, nbas
        do  jb = 1, nbas
          if (ib .ne. jb) then
            ixi = indH(ib)
            ixj = indH(jb)
            nlmi = nlm(ipc(ib))
            nlmj = nlm(ipc(jb))
            do lmi = 1, nlmi
              do  lmj = 1, nlmj
                i = indxsh(ixi+lmi)
                j = indxsh(ixj+lmj)
                if (i .le. ldim .and. j .le. ldim) then
                  hk(i,j) = hk(i,j) + (pot0(ib) + pot0(jb)) * sk(i,j)
     .                      * 0.5d0
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      end 
