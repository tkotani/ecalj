subroutine fradpk(kmax,rsma,lmxa,nr,rofi,fp,xp,vp,dp)
  !- Set up radial pkl functions, grad of same, vals and slopes.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   kmax  :polynomial cutoff in augmentation expansion
  !i   rsma  :augmentation smoothing radius
  !i   lmxa  :augmentation L-cutoff
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !o Outputs
  !o   fp    :(Pkl) * r
  !o   xp    :(Pkl') * r
  !o   vp    :Pkl at rofi(nr) = fp/rofi(nr)
  !o   dp    :Pkl' at rofi(nr)
  !r Remarks
  !r   Radial functions are multiplied by r, but vp and dp are not.
  !u Updates
  !u  16 May 00 Adapted from nfp frapkl.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: kmax,lmxa,nr
  double precision :: rofi(nr),fp(nr,0:lmxa,0:kmax),rsma, &
       xp(nr,0:lmxa,0:kmax),vp(0:lmxa,0:kmax),dp(0:lmxa,0:kmax)
  ! ... Local parameters
  integer :: k0,l0,l,k,i
  parameter ( k0=30, l0=7 )
  double precision :: pkl(0:k0,0:l0),rmt,fac,xx,r,rlp1,rr,rr2

  if (kmax+1 > k0) call rxi('fradpk: need k0 >= ',kmax+1)
  if (lmxa+1 > l0) call rxi('fradpk: need l0 >= ',lmxa+1)
  if (lmxa == -1) return

  ! ... Values and slopes of the pkl
  rmt = rofi(nr)
  call radpkl(rmt,rsma,kmax,lmxa+1,k0,pkl)
  do  l = 0, lmxa
     do  k = 0, kmax
        fac = (2*k+2*l+3d0)/(2*l+3d0)
        xx = pkl(k,l)
        dp(l,k) = 2*rmt**(l+1)/rsma**2*(xx-pkl(k,l+1)*fac*rsma) &
             + l*rmt**(l-1)*xx
        vp(l,k) = rmt**l*pkl(k,l)
     enddo
  enddo

  ! ... Set up the functions pkl(r) and derivative on radial mesh
  do  i = 1, nr
     r = rofi(i)
     call radpkl(r,rsma,kmax,lmxa+1,k0,pkl)
     do  l = 0, lmxa
        rlp1 = 0d0
        if ( r > 0d0 ) rlp1 = r**(l+1)
        do  k = 0, kmax
           fp(i,l,k) = pkl(k,l)*rlp1
           fac = (2*k+2*l+3d0)/(2*l+3d0)
           xx = pkl(k,l)
           rr = 0d0
           if ( r > 0d0 ) rr = r**l
           rr2 = rr * r * r
           !           xp(i,l,k) = 2*r**(l+2)/rsma**2*(xx-pkl(k,l+1)*fac*rsma)
           !    .         + l*r**l*xx
           xp(i,l,k) = 2*rr2/rsma**2*(xx-pkl(k,l+1)*fac*rsma) &
                + l*rr*xx
        enddo
     enddo
  enddo

end subroutine fradpk

