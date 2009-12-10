#define F90 1
      subroutine prjrsp(nsp,nspc,icplx,nRLc,nkp,P0)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   icplx :1 if P0 is real, 2 if P0 is complex
Ci   nRLc  :dimensions P0
Ci   nkp   :number of k-points which there is a P0
Cio Inputs/Outputs
Ci   P0     :some rows and columns of P0 are projected out to conserve
Ci          charge neutrality; see Remarks
Cl Local variables
Cl   nspx   :number of separate spin channels
Cl           =1 except in collinear spin-polarized case
Cl   nRLx   :nRLc * nspc
Cr Remarks
Cr  Renormalization:  the response function P0 as calculated from Dyson's
Cr  equation does not preserve charge neutrality.  Thus for the last
Cr  energy point (last=.true.), gfp0ft projects out a portion of the
Cr  response function P0, so that any single potential shift preserves
Cr  charge neutrality, and a uniform potential shift conserves charge in
Cr  each of n channels.  Let T be the following projector matrix
Cr
Cr                (1 1 ... 1)
Cr      T =  n^-1 (1 1 ... 1)   with n = number of rows in T
Cr                (...      )
Cr                (1 1 ... 1)
Cr
Cr  All eigenvalues of T are zero except except one whose eigenvalue is
Cr  1, and whose eigenvector is z=(1,1,...,1)/sqrt(n).  Thus,
Cr  1-T = z+ eps z   operating on a matrix P0 subtracts out the
Cr  projection of P0 onto this eigenvector.  The dual projection
Cr      P0' = (1-T) P0 (1-T)
Cr  projects out the sum of charges for a potential shift in any channel
Cr  and also the induced charge in any channel arising from a uniform
Cr  shift in all channels.  It is readily seen that
Cr      (P0' - P0)_ij = 1/n sum_k (-P0_kj - P0_ik + 1/n sum_l P0_lk)
Cr  The response function P0_ij, defined as dq_i = sum_j P0_ij dV_j
Cr  is merely the imaginary part of P0_ij.
Cu Updates
Cu   24 May 02 Extended to include nkp q-points
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,nspc,nRLc,icplx,nkp
      double precision P0(icplx,nRLc,nspc,nRLc,nkp,nsp)
C ... Local parameters
      integer  k,i,j,isp,nspx,nRLc2,ikp
#if F90 | AUTO_ARRAY
      double precision sum,rsum(nRLc*2),csum(nRLc*2)
#else
      double precision sum,rsum(1),csum(1)
      stop 'prjrsp requires automatic arrays'
#endif

C     call prm('input P0',k,P0,nRLc*nsp,nRLc*nsp,nRLc*nsp)

      nspx = nsp
      if (nspc .eq. 2) nspx = 1
      nRLc2 = nRLc*nspc
      do  10  ikp = 1, nkp
      do  10  isp = 1, nspx
      do  10  i = 1, icplx

      sum = 0
      do  22  k = 1, nRLc2
        rsum(k) = 0
        csum(k) = 0
        do  24  j = 1, nRLc2
          rsum(k) = rsum(k) + P0(i,j,1,k,ikp,isp)
          csum(k) = csum(k) + P0(i,k,1,j,ikp,isp)
          sum     = sum + P0(i,j,1,k,ikp,isp)
   24   continue
        rsum(k) = rsum(k)/nRLc2
        csum(k) = csum(k)/nRLc2
   22 continue
      sum = sum/nRLc2**2
      do  26  k = 1, nRLc2
      do  26  j = 1, nRLc2
        P0(i,k,1,j,ikp,isp) = P0(i,k,1,j,ikp,isp) - rsum(k)
        P0(i,j,1,k,ikp,isp) = P0(i,j,1,k,ikp,isp) - csum(k)
        P0(i,k,1,j,ikp,isp) = P0(i,k,1,j,ikp,isp) + sum
   26 continue
   10 continue

C     call prm('final P0',k,P0,nRLc*nsp,nRLc*nsp,nRLc*nsp)
      end

