subroutine fradhd(nkaps,eh,rsmh,lh,lmxh,nr,rofi,fh,xh,vh,dh) !Radial envelope head functions and gradients
  use m_hansr,only: hansr
  !i Inputs
  !i   nkaps :1s digit: number of envelope function types per l q.n.
  !i   eh    :energies of smoothed Hankel (l- and kappa- dependent)
  !i   rsmh  :hankel smoothing radii (l- and kappa- dependent)
  !i         :rsmh=0 -> channel has no envelope function
  !i   lh    :basis l-cutoff for each of nkape groups
  !i   lmxh  :dimensions xh,dh,fh,vh
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !o Outputs
  !o   fh    :(envelope head functions) * r
  !o   xh    :(radial derivative of head functions) * r = d(fh)/dr - fh/r
  !o   vh    :function values at rofi(nr) = fh/rofi(nr)
  !o   dh    :function slopes at rofi(nr)
  !r Remarks
  !r   Radial functions are multiplied by r, but vh and dh are not.
  !r   Derivative calculate analytically, as
  !r     xh = h' = l/r h_l - h_l+1
  !u Updates
  !u   12 Aug 04 First implementation of extended local orbitals;
  !u             substitute call to hansmr with call to hansr.
  !u   10 Apr 02 Redimensioned eh,rsmh to accomodate larger lmax
  !u   16 May 00 Adapted from nfp fradhd.f
  implicit none
  integer :: nkaps,lh(nkaps),lmxh,nr,n0
  parameter (n0=10)
  double precision :: eh(n0,nkaps),rsmh(n0,nkaps),rofi(nr), &
       fh(nr,0:lmxh,nkaps),vh(0:lmxh,nkaps), &
       xh(nr,0:lmxh,nkaps),dh(0:lmxh,nkaps)
  integer :: ik,l,i,idx
  double precision :: rsm,e,rmt,r,vhh,dhh,xi(0:20),rl !,wk(2)
  fh=0d0; xh=0d0; vh=0d0; dh=0d0
  do  ik = 1, nkaps
     do  l = 0, lh(ik)
        rsm = rsmh(l+1,ik)
        e = eh(l+1,ik)
        if (rsm > 0) then
           rmt = rofi(nr)
           rl = rmt**l
           if (l == 0) rl = 1

           !           asm = 1d0/rsm
           !           call hansmr(rmt,e,asm,xi,l+1)
           !           vh(l,ik) = xi(l) * rmt**l
           !           dh(l,ik) = vh(l,ik)*l/rmt - xi(l+1)*rmt**(l+1)

           !  call hansr(rsm,0,l+1,1,l+1,e,rmt**2,1,1,idx,wk,10,xi) !we calculate xi at nr.
           call hansr(rsm,0,l+1,1,[l+1],[e],[rmt**2],1,1,[idx],10,xi) !we calculate xi at nr.
           vh(l,ik) = xi(l) * rl
           dh(l,ik) = vh(l,ik)*l/rmt - xi(l+1)*rl*rmt

           fh(1,l,ik) = 0
           xh(1,l,ik) = 0
           do  i = 2, nr

              r = rofi(i)
              rl = r**l
              if (l == 0) rl = 1

              !     call hansmr(r,e,asm,xi,l+1)
              !     call hansr(rsm,0,l+1,1,l+1,e,r**2,1,1,idx,wk,10,xi) !we calculate xi at r**2.
              call hansr(rsm,0,l+1,1,[l+1],[e],[r**2],1,1,[idx],10,xi) !we calculate xi at r**2.

              !             g0 = dexp(-asm*asm*r*r)
              vhh = xi(l) * rl
              dhh = vhh*l/r - xi(l+1)*rl*r
              fh(i,l,ik) = vhh*r
              xh(i,l,ik) = dhh*r

              !             call hansmd(1,r,e,rsm,l,hs,dhs,ddhs,xx,xx,xx)
              !             fh(i,l,ik) = hs(l)*r
              !             fh(i,l,ik) = dhs(l)*r

           enddo
        endif
     enddo
  enddo
end subroutine fradhd

