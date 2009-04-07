      subroutine broydn (jac,df,dx,beta,gamma,bx,btx,namx,wc)
C- Calculate the jacobian matrix J^(m+1) according to eq. A13-A16 of
C  D. Vanderbilt, S. G. Louie in Phys. Rev B 30, 6118 (1984)
C ------------------------------------------------------------------
Ci Inputs:
Ci   df    :df^(m+1) = F^(m+1) - F^(m)  (unnormalized)
Ci   dx    :dx^(m+1) = x^(m+1) - x^(m)  (unnormalized)
Cio  beta  :defines jac
Cio  gamma :defines jac
Cw   bx    :work array
Cw   btx   :work array
Ci   namx  :upper limit to number of elements in the vector
Ci   wc    :mixing parameter (wc large: output vector strongly weighted)
Co Outputs:
Co   jac   :Jacobian matrix =  gamma * beta**(-1)
Cio  beta  :is updated
Cio  gamma :is updated
C ----------------------------------------------------------------
C     implicit none
C  Passed variables:
      integer namx
      double precision jac(namx,namx), beta(namx,namx),
     .                 gamma(namx,namx),df(namx),dx(namx),
     .                 bx(namx),btx(namx),wc

C  Local variables
      integer i,j,iprint
      double precision xn,sigma,wcx,tol
C External calls
      external dmpy,iprint

C --- Calculate the norm of x^(m+1)-x^(m) ---
      xn = 0
      do  10  i = 1, namx
   10 xn = xn + dx(i)**2

      tol = 1d-12
      if (xn .lt. tol) then
        if (iprint().ge.40) write(*,101)
  101   format(' BROYDN:  dx very small, matrix not updated')
      else
        wcx = wc / xn
        call dmpy(beta,namx,1,dx,namx,1,bx,namx,1,namx,1,namx)
        call dmpy(beta,1,namx,dx,namx,1,btx,namx,1,namx,1,namx)
        call dmpy(bx,1,namx,dx,namx,1,sigma,1,1,1,1,namx)
        sigma = sigma*wcx + 1
        sigma = wcx / sigma

C ---   Update gamma (eq A15), beta (implicit inverse) (eq. A14,A16) ---
C       NB: beta is beta inverse
        do  20   i = 1, namx
        do  20  j = 1, namx
          gamma(i,j) = gamma(i,j) - wc*df(i)*dx(j)/xn
          beta (i,j) = beta (i,j) - sigma*bx(i)*btx(j)
   20   continue
      endif

C --- J^(m+1) = gamma^(m+1) * beta^(m+1)**(-1) (eq. A13) ---
      call dmpy(gamma,namx,1,beta,namx,1,jac,namx,1,namx,namx,namx)

      end

