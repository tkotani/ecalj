      subroutine addtos(nds,nda,iax,nenv,nkap,npr,salph,nprs,s)
C- Appends structure constant vector for single sphere to s
C ----------------------------------------------------------------
Ci Inputs
Ci   nds   :leading dimensions of s
Ci   nda   :leading dimension of salph
Ci   nenv  :number of envelope functions; second dimension of salph
Ci   iax   :structure-constant information relating to s
Ci   nkap  :number of families for a given orbital
Ci   npr   :number of neighbors in current cluster
Ci   sc    :structure constants connecting one site to all its neighbors
Ci   nprs  :total number of strux accumulated before this one
Co Outputs
Co   salph :is appended to s
Cb Bugs
Cb   No check is made to ensure that s is sufficiently large.
Cr Remarks
Cr   Converts s to format s(L',L,R'), the indices corresponding
Cr   to S_ R'L', RL.  Indices RL correspond to sphere and L channel
Cr   of the head of the screened orbital and R'L' to the sphere and
Cr   L channel of which S is the one-center expansion for the screened
Cr   function.  Expansion theorem is:
Cr     \sum_R''L'' K_R''L''(r) (1 + alpha_R''L'' S_R''L'',RL) =
Cr     \sum_L' (J_R'L'(r)  - alpha_R'L' K_R'L'(r))  S_R'L',RL
Cr   which is analagous to the expansion theorem Eqn (3) in mstrux.
Cr
Cu Updates
Cu   06 Aug 06 Extra argument nkap
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nds,nda,nenv,nkap,npr,nprs,niax
      parameter (niax=10)
      integer iax(niax,npr)
      double precision salph(nda,nkap,nenv,nkap),
     .  s(nds,nds,nkap,nkap,nprs+npr)
C Local parameters
      integer ipr,ilm,klm,ii,nlma,nlmb,ikap,jkap

      nlmb = iax(9,1)
      ii = 0
      do  ipr = 1, npr
        nlma = iax(9,ipr)
        do  ikap = 1, nkap
        do  jkap = 1, nkap
          do  klm = 1, nlmb
          do  ilm = 1, nlma
            s(ilm,klm,ikap,jkap,ipr+nprs) = salph(ii+ilm,ikap,klm,jkap)
          enddo
          enddo
        enddo
        enddo
        ii = ii+nlma
      enddo

      end

