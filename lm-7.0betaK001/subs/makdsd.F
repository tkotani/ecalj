      subroutine makdsd(opts,ldr,ldc,nrow,ncol,ofpr,ofpl,pph,s,dsd)
C- Calculate sqrdel*s*sqrdel from structure constants and pot pars
C ----------------------------------------------------------------------
Ci Inputs
Ci   opts  :a compound of digits specifying scaling options
Ci         :1s digit
Ci         :0 add nothing to diagonal
Ci         :1 add pph(2)-pph(1) to diagonal
Ci         :10s digit
Ci         :0 scale both left and right of s by pph(3)
Ci         :1 scale right only
Ci         :2 scale left only
Ci         :3 scale neither
Ci   ldr   :leading (row) dimension of s and dsd
Ci   ldc   :second  (column) dimension of s and dsd
Ci   nrow  :number of rows to scale
Ci   ncol  :number of columns to scale
Ci   ofpr  :right offset to pph
Ci   ofpl  :left offset to pph
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci   s     :structure constants
Co Outputs
Co   dsd   :pph(3) * S * pph(3) ... + pph(2) - pph(1)
Co         :left, right factors, or diagonal terms may be suppressed,
Co         :depending on opts
Cr Remarks
Cr   s and dsd may point to the same address space
Cu Updates
Cu   19 Jan 01 revised and extended; argument list changed
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer opts,ldr,ldc,ofpr,ofpl,nrow,ncol
      double precision pph(5,*)
      double precision s(ldr,ldc,2),dsd(ldr,ldc,2)
C local parameters
      integer i,j,opt0,opt1

      opt0 = mod(opts,10)
      opt1 = mod(opts/10,10)

      if (opt1 .eq. 0) then
#if SGI_PARALLEL
C$DOACROSS local(i,j)
#endif
        do  i = 1, nrow
          do  j = 1, ncol
            dsd(i,j,1) = pph(3,i+ofpl)*s(i,j,1)*pph(3,j+ofpr)
            dsd(i,j,2) = pph(3,i+ofpl)*s(i,j,2)*pph(3,j+ofpr)
          enddo
        enddo
      elseif (opt1 .eq. 1) then
#if SGI_PARALLEL
C$DOACROSS local(i,j)
#endif
        do  i = 1, nrow
          do  j = 1, ncol
            dsd(i,j,1) = s(i,j,1)*pph(3,j+ofpr)
            dsd(i,j,2) = s(i,j,2)*pph(3,j+ofpr)
          enddo
        enddo
      elseif (opt1 .eq. 2) then
#if SGI_PARALLEL
C$DOACROSS local(i,j)
#endif
        do  i = 1, nrow
          do  j = 1, ncol
            dsd(i,j,1) = pph(3,i+ofpl)*s(i,j,1)
            dsd(i,j,2) = pph(3,i+ofpl)*s(i,j,2)
          enddo
        enddo
      elseif (opt1 .eq. 3) then
#if SGI_PARALLEL
C$DOACROSS local(i,j)
#endif
        do  i = 1, nrow
          do  j = 1, ncol
            dsd(i,j,1) = s(i,j,1)
            dsd(i,j,2) = s(i,j,2)
          enddo
        enddo
      endif

C --- Add diagonal ---
      if (opt0 .eq. 0) return
      do  i = 1, min(nrow,ncol)
        dsd(i,i,1) = dsd(i,i,1) + pph(2,i) - pph(1,i)
      enddo
      end

