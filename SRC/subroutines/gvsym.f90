subroutine gvsym(ng,gv,ips0,bgv,c,csym)
  !- Symmetrize a function c, given in the form of a list
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   gv,ng   Lattice vectors, and number
  !i   ips0    pointer to first vector in star of this vector; see sgvsym
  !i   bgv     phase factor sum; see sgvsym
  !i   c       unsymmetrized function
  !o Outputs:
  !o   csym    symmetrized function
  !r Remarks:
  !u    7 Sep 98 Adapted from nfp gvsym.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: ng,ips0(ng)
  double precision :: gv(ng,3)
  double complex bgv(ng),c(ng),csym(ng)
  integer :: i,j,i0,kstar,ipr,iprint,nstar

  ! ... Sum up coefficients for first vector in each star
  do  10  i = 1, ng
     csym(i) = (0d0,0d0)
10 enddo
  do  12  i = 1, ng
     j = ips0(i)
     csym(j) = csym(j) + bgv(i)*c(i)
12 enddo

  ! ... Normalize
  do  20  i0 = 1, ng
     if (ips0(i0) == i0) then
        kstar = 0
        do  22  i = i0, ng
           if (ips0(i) == i0) kstar = kstar+1
22      enddo
        csym(i0) = csym(i0)/kstar
     endif
20 enddo

  ! ... Make all the coefficients
  do  30  i = 1, ng
     j = ips0(i)
     csym(i) = csym(j)*dconjg(bgv(i))
30 enddo

  ! ... Printout
  ipr = iprint()
  if (ipr <= 55) return
  print 255
  nstar = 0
  do  40  i0 = 1, ng
     if (ips0(i0) == i0) then
        nstar = nstar+1
        if (ipr >= 60) print *, ' '
        do  44  i = i0, ng
           if (ips0(i) == i0) then
              if (i == i0) then
                 print 251, nstar,i,gv(i,1),gv(i,2),gv(i,3), &
                      c(i),csym(i)
              else
                 if (ipr >= 60) &
                      print 252, i,gv(i,1),gv(i,2),gv(i,3), &
                      c(i),csym(i)
              endif
251           format(i4,i5,3f6.1,2f12.8,1x,2f12.8)
252           format(4x,i5,3f6.1,2f12.8,1x,2f12.8)
255           format(/' star  ig',8x,'recip',17x,'c_in',20x,'c_sym')
           endif
44      enddo
     endif
40 enddo

end subroutine gvsym

