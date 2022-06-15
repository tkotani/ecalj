subroutine symprj(nrclas,nlml,ngrp,nbas,istab,g,ag,plat,qlat, &
     pos,sym)
  !- Set up symmetry projectors for one class
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nrclas:number of atoms in this class
  !i   nlml  :L-cutoff to which symmetry projectors are calculated
  !i   ngrp  :size of space group
  !i   nbas  :size of basis
  !i   istab :not used
  !i   g     :point group operations
  !i   ag    :translation part of space group
  !i   plat  :primitive lattice vectors, in units of alat
  !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
  !i   pos   :basis vectors, in units of alat
  !o Outputs
  !o   sym   :symmetry projectors for each site within this class
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nlml,nrclas,ngrp,nbas,istab(nbas,ngrp)
  double precision :: sym(nlml,nlml,nrclas),plat(3,3),qlat(3,3), &
       g(3,3,ngrp),ag(3,ngrp),pos(3,nrclas),d(3)
  integer:: lmxl , ll , ig , ja , m , ia=99999 , iprint , ilm , l , jlm,  jlm1 , jlm2
  real(8),allocatable :: rmat_rv(:,:)
  double precision :: wgt
  real(8):: rfrac(3)
  lmxl = ll(nlml)
  allocate(rmat_rv(nlml,nlml))
  wgt = 1d0/ngrp
  sym = 0d0
  ! --- For each group operation, do ---
  do  10  ig = 1, ngrp
     ! ...   Find site mapped into first site under this operation
     do  12  ja = 1, nrclas
        do  14  m = 1, 3
           d(m) = g(m,1,ig)*pos(1,ja) + g(m,2,ig)*pos(2,ja) &
                + g(m,3,ig)*pos(3,ja) + ag(m,ig) - pos(m,1)
14      enddo
!        call shorbz(d,d,plat,qlat)
       ia = ja
       rfrac = matmul(d,qlat)
!       write(6,*)ig,ja,sum(abs(rfrac-nint(rfrac)))
       if(sum((rfrac-nint(rfrac))**2)< 1d-7) goto 80
!        if (d(1)*d(1)+d(2)*d(2)+d(3)*d(3) < 1d-7) goto 80
        !          if (d(1)*d(1)+d(2)*d(2)+d(3)*d(3) .lt. 1d-9) goto 80
12   enddo
     call rxi('symprj: no site mapped into first under op',ig)
80   continue
     ! ...   Make and add transformation matrix
     call ylmrtg ( nlml , g ( 1 , 1 , ig ) , rmat_rv )
     sym(:,:,ia)= sym(:,:,ia) + wgt*rmat_rv
     !        call dpadd ( sym ( 1 , 1 , ia ) , rmat_rv , 1 , nlml * nlml
     !     .  , wgt )
10 enddo
  if( iprint() >= 60 ) then
     do 20  ja = 1, nrclas
        write(6,727) ja
727     format(/' projection matrices for ja=',i3)
        ilm = 0
        do  22  l = 0, lmxl
           jlm1 = l*l+1
           jlm2 = (l+1)**2
           do  24  m = 1, 2*l+1
              ilm = ilm+1
              write(6,210) (sym(ilm,jlm,ja),jlm=jlm1,jlm2)
210           format(1x,9f12.6)
24         enddo
22      enddo
20   enddo
  endif
  deallocate(rmat_rv)
end subroutine symprj

