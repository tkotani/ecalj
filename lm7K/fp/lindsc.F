      subroutine lindsc(job,ng,gv,tpiba,elind,cn)
C- Make screened density, using a Lindhard response function
C ----------------------------------------------------------------------
Ci Inputs
Ci   job  1s digit specifies what cn contains on output;
Ci          0 cn is untouched
Ci          1 cn is overwritten by eps^-1 cn
Ci          2 cn is overwritten by (eps^-1-1) cn
Ci          3 cn is overwritten by eps^-1.  Input cn is irrelevant
Ci   ng    :number of G-vectors
Ci   gv    :list of dimensionless reciprocal lattice vectors (gvlist.f)
Ci   tpiba :2*pi/lattice constant
Ci   elind :Lindhard screening parameter
Co Outputs
Co   cn:     overwritten by
Co           cn untouched                (job = 0)
Co           eps^-1 (input cn)           (job = 1)
Co           (eps^-1 - 1)(input cn)      (job = 2)
Co           eps                         (job = 3)
Co           eps^-1                      (job = 4)
Cr Remarks
Cr   The Thomas-Fermi dielectric response is:
Cr     eps_TF(q) = 1 + 4 (k_F a_0) / pi (q a_0)^2
Cr               = 1 + 4 k_F / pi q^2 (Rydberg units)
Cr   The Lindhard function is
Cr     eps(q) = 1 + 4 k_F / pi q^2 * [...] , where
Cr              [...] = 1/2 + (1-x*x)/4x ln((1+x)/(1-x)) and
Cr              x = q / 2 k_F.
Cu Updates
Cu   28 Oct 01 routine revamped.  Original was ridiculously complicated.
C ----------------------------------------------------------------------
C     implicit none
      integer job,ng
      double precision gv(ng,3),tpiba,elind
      double complex cn(ng)
C Local variables
      integer i
      double precision g2,xx,yy,eps,pi

      pi = 4d0*datan(1d0)

Ckino isanrg is logical function,       call isanrg(job,0,4,'lindsc:','job', .true.)
      l_dummy_isanrg=isanrg(job,0,4,'lindsc:','job', .true.)
      if (job .eq. 0) return

      cn(1) = 0

C ... Early exit if elind is zero => eps = 1
      if (elind .eq. 0) then
        if (job .eq. 1) then
        elseif (job .eq. 2) then
          call dpzero(cn,ng*2)
        elseif (job .eq. 3) then
          do  i = 2, ng
            cn(i) = 1
          enddo
        endif
        return
      endif

      do  22  i = 2, ng
        g2 = tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
        xx = sqrt(g2/elind)/2
        yy = 0.5d0 + (1-xx**2)/(4*xx)*dlog(dabs((1+xx)/(1-xx)))
        eps = 1 + 4*dsqrt(elind)/(pi*g2)*yy
        if (job .eq. 1) then
          cn(i) = cn(i) / eps
        elseif (job .eq. 2) then
          cn(i) = cn(i) * (1/eps-1)
        elseif (job .eq. 3) then
          cn(i) = eps
        elseif (job .eq. 4) then
          cn(i) = 1 / eps
        endif
   22 continue

      end


