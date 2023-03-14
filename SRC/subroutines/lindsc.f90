subroutine lindsc(job,ng,gv,tpiba,elind,cn)! Make screened density, using a Lindhard response function
  !i Inputs
  !i   job  1s digit specifies what cn contains on output;
  !i          0 cn is untouched
  !i          1 cn is overwritten by eps^-1 cn
  !i          2 cn is overwritten by (eps^-1-1) cn
  !i          3 cn is overwritten by eps^-1.  Input cn is irrelevant
  !i   ng    :number of G-vectors
  !i   gv    :list of dimensionless reciprocal lattice vectors (gvlist.f)
  !i   tpiba :2*pi/lattice constant
  !i   elind :Lindhard screening parameter
  !o Outputs
  !o   cn:     overwritten by
  !o           cn untouched                (job = 0)
  !o           eps^-1 (input cn)           (job = 1)
  !o           (eps^-1 - 1)(input cn)      (job = 2)
  !o           eps                         (job = 3)
  !o           eps^-1                      (job = 4)
  !r Remarks
  !r   The Thomas-Fermi dielectric response is:
  !r     eps_TF(q) = 1 + 4 (k_F a_0) / pi (q a_0)^2
  !r               = 1 + 4 k_F / pi q^2 (Rydberg units)
  !r   The Lindhard function is
  !r     eps(q) = 1 + 4 k_F / pi q^2 * [...] , where
  !r              [...] = 1/2 + (1-x*x)/4x ln((1+x)/(1-x)) and
  !r              x = q / 2 k_F.
  !u Updates
  !u   28 Oct 01 routine revamped.  Original was ridiculously complicated.
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: job,ng
  double precision :: gv(ng,3),tpiba,elind
  double complex cn(ng)
  ! Local variables
  integer :: i
  double precision :: g2,xx,yy,eps,pi
  logical:: l_dummy_isanrg,isanrg

  pi = 4d0*datan(1d0)

  ! ino isanrg is logical function,       call isanrg(job,0,4,'lindsc:','job', .true.)
  l_dummy_isanrg=isanrg(job,0,4,'lindsc:','job', .true.)
  if (job == 0) return

  cn(1) = 0

  ! ... Early exit if elind is zero => eps = 1
  if (elind == 0) then
     if (job == 1) then
     elseif (job == 2) then
        call dpzero(cn,ng*2)
     elseif (job == 3) then
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
     if (job == 1) then
        cn(i) = cn(i) / eps
     elseif (job == 2) then
        cn(i) = cn(i) * (1/eps-1)
     elseif (job == 3) then
        cn(i) = eps
     elseif (job == 4) then
        cn(i) = 1 / eps
     endif
22 enddo

end subroutine lindsc


