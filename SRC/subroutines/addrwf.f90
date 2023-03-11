module m_addrwf !- Add constant * radial wave function to another radial wave function
  private
  public addrwf
contains
  subroutine addrwf(mode,z,l,v,ndg,n1,nr,rofi,rwgt,eadd,ev,fac,gadd, g,s)
    use m_vxtrap,only: rwftai
    use m_lmfinit,only: cc
    ! ----------------------------------------------------------------------
    !i Inputs
    !i  mode   :0 use both large and small components of radial w.f.
    !i         :1 use both large component of radial w.f. only
    !i   z     :nuclear charge
    !i         :(used only to compute overlap s, mode=0)
    !i   l     :l-quantum number
    !i   v     :spherical potential, without nuclear part
    !i         :(used only to compute overlap s, mode=0)
    !i   ndg   :leading dimension of g and gadd
    !i   n1    :if 0<n1<=nr, rwgt(n1) is scaled by 2
    !i         :(see vxtrap.f)
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   rwgt  :radial mesh weights for numerical integration
    !i   ev    :eigenvalue of input wave function g
    !i         :(used only to compute overlap s, mode=0)
    !i   eadd  :eigenvalue of wave function gadd
    !i         :(used only to compute overlap s, mode=0)
    !i   fac   :Add fac*gadd into g
    !i   gadd  :See fac
    !o Inputs/Outputs
    !o   g     :g is overwritten by g + fac*g
    !o   s     :overlap between gadd and new g
    !r Remarks
    !r   Input g and gadd are assumed to be solutions of the Schrodinger
    !r   equation with eigenvalues ev and eadd.  (For the scalar
    !r   relativistic case, the inner product depends slightly
    !r   on z,v, and eigenvalue)
    !u Updates
    !u   04 Sep 04 Adapted to extended local orbitals
    !u   12 Jul 04 ndg,n1 arguments (altered argument list)
    !u   14 Feb 02 New routine
    ! ----------------------------------------------------------------------
    implicit none
    integer :: l,ndg,n1,nr,mode
    real(8):: fac,rofi(nr),rwgt(nr), gadd(ndg,2),g(ndg,2),s
    real(8),optional:: z,v(nr),ev,eadd
    integer :: ir
    double precision :: vi,fllp1,gf11,gf22,gf12,r,tmc
    fllp1 = l*(l+1)
    s = 0
    if (n1 > 0 .AND. n1 < nr) rwgt(n1) = 2*rwgt(n1)
    if (mode == 0) then
       do  ir = 2, nr
          r = rofi(ir)
          if (fac /= 0) then
             g(ir,1) = g(ir,1) + fac*gadd(ir,1)
             g(ir,2) = g(ir,2) + fac*gadd(ir,2)
          endif
          !       Rest of the loop computes overlap between new g and gadd
          vi = v(ir) - 2d0*z/r
          tmc = cc - (vi-ev)/cc
          gf11 = 1d0 + fllp1/(tmc*r)**2
          tmc = cc - (vi-eadd)/cc
          gf22 = 1d0 + fllp1/(tmc*r)**2
          gf12 = (gf11 + gf22)/2
          s = s + rwgt(ir)*(gf12*g(ir,1)*gadd(ir,1) + g(ir,2)*gadd(ir,2))
       enddo
    else
       do  ir = 2, nr
          r = rofi(ir)
          if (fac /= 0) then
             g(ir,1) = g(ir,1) + fac*gadd(ir,1)
          endif
          s = s + rwgt(ir)*g(ir,1)*gadd(ir,1)
       enddo
    endif
    if (n1 > 0 .AND. n1 < nr) rwgt(n1) = rwgt(n1)/2
  end subroutine addrwf
end module m_addrwf

