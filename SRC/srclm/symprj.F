      subroutine symprj(nrclas,nlml,ngrp,nbas,istab,g,ag,plat,qlat,
     .pos,sym)
C- Set up symmetry projectors for one class
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrclas:number of atoms in this class
Ci   nlml  :L-cutoff to which symmetry projectors are calculated
Ci   ngrp  :size of space group
Ci   nbas  :size of basis
Ci   istab :not used
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   plat  :primitive lattice vectors, in units of alat
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   pos   :basis vectors, in units of alat
Co Outputs
Co   sym   :symmetry projectors for each site within this class
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlml,nrclas,ngrp,nbas,istab(nbas,ngrp)
      double precision sym(nlml,nlml,nrclas),plat(3,3),qlat(3,3),
     .g(3,3,ngrp),ag(3,ngrp),pos(3,nrclas),d(3)
      integer:: lmxl , ll , ig , ja , m , ia , iprint , ilm , l , jlm,  jlm1 , jlm2
      real(8),allocatable :: rmat_rv(:,:)
      double precision wgt
      lmxl = ll(nlml)
      allocate(rmat_rv(nlml,nlml))
      wgt = 1d0/ngrp
      sym = 0d0
C --- For each group operation, do ---
      do  10  ig = 1, ngrp
C ...   Find site mapped into first site under this operation
        do  12  ja = 1, nrclas
          do  14  m = 1, 3
            d(m) = g(m,1,ig)*pos(1,ja) + g(m,2,ig)*pos(2,ja)
     .      + g(m,3,ig)*pos(3,ja) + ag(m,ig) - pos(m,1)
   14     continue
          call shorbz(d,d,plat,qlat)
          ia = ja
          if (d(1)*d(1)+d(2)*d(2)+d(3)*d(3) .lt. 1d-7) goto 80
c          if (d(1)*d(1)+d(2)*d(2)+d(3)*d(3) .lt. 1d-9) goto 80
   12   continue
        call rxi('symprj: no site mapped into first under op',ig)
   80   continue
C ...   Make and add transformation matrix
        call ylmrtg ( nlml , g ( 1 , 1 , ig ) , rmat_rv )
        sym(:,:,ia)= sym(:,:,ia) + wgt*rmat_rv 
c        call dpadd ( sym ( 1 , 1 , ia ) , rmat_rv , 1 , nlml * nlml 
c     .  , wgt )
   10 continue
      if( iprint >= 60 ) then
      do 20  ja = 1, nrclas
        write(6,727) ja
  727   format(/' projection matrices for ja=',i3)
        ilm = 0
        do  22  l = 0, lmxl
          jlm1 = l*l+1
          jlm2 = (l+1)**2
          do  24  m = 1, 2*l+1
            ilm = ilm+1
            write(6,210) (sym(ilm,jlm,ja),jlm=jlm1,jlm2)
  210       format(1x,9f12.6)
   24     continue
   22   continue
 20   continue
      endif
      deallocate(rmat_rv)
      end

