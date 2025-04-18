subroutine psymr0(lmxl,ic,nbas,ipc,pos0,pos,ipa,nrclas)  !- Make list of atoms in one class
  !i   lmxl  :lmxl(1..nbas) = l-cutoffs dimensioning offsets.
  !i         :Used in creating table offsets for table ipa, which see
  !i         :lmxl(1) < -1 => ipa is created differently and
  !i         :lmxl(2..) are not referenced; see description for ipa.
  !i   ic    :abs. val of ic = class for which to make list of members
  !i         :sign of ic<0 : flag to suppress generation of pos.
  !i   nbas  :size of basis
  !i   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
  !i   pos0  :basis vectors for each site (not used if ic<0)
  !o Outputs
  !o   ipa   :ipa(1..nrclas) is created in one of two modes:
  !o         :if lmxl(1)<-1 :
  !o         : ipa = list of sites belonging to class ic
  !o         :       In this mode lmxl is not used.
  !o         :if lmxl(1)>=-1 :
  !o         : ipa = table of offsets to an array for each site in class
  !o         :       ic.  It is intended for some array A containing
  !o         :       elements for each site 1..nbas, and that offset
  !o         :       to the A for site ib+1, offs(ib+1) is
  !o         :       (lmxl(ib)+1)**2 + offs(ib).
  !o   pos   :basis vectors for each member of class (made if ic>0)
  !o   nrclas:number of sites in class ic
  implicit none
  integer :: nbas,lmxl(nbas),ic,nrclas,ipc(nbas),ipa(nbas)
  double precision :: pos0(3,nbas),pos(3,nbas)
  logical :: loff
  integer :: ib,iclass,ioff
  logical :: lpos
  loff = lmxl(1) .ge. -1
  iclass = iabs(ic)
  lpos = ic .gt. 0
  nrclas = 0
  ioff = 0
  do  ib = 1, nbas
     if (ipc(ib) == iclass) then
        nrclas = nrclas+1
        if ( .NOT. loff) ioff = ib
        ipa(nrclas) = ioff
        if (lpos) pos(:,nrclas) = pos0(:,ib)
     endif
     ioff = ioff + (lmxl(ib)+1)**2
  enddo
end subroutine psymr0
subroutine pysmr1(wgt,nr,nlml,nsp,sym,rho1, rho2,nn)! Add wgt*sym*rho1 into rho2, transposed
  implicit none
  integer :: nr,nlml,nsp,nn
  double precision :: sym(nlml,nlml),rho1(nr,nlml,nsp), rho2(nr,nlml,nsp),wgt
  integer :: i,ilm,jlm,isp
  double precision :: ccc
  nn = 0
  do  isp = 1, nsp
     do  ilm = 1, nlml
        do  jlm = 1, nlml
           ccc = wgt*sym(jlm,ilm)
           if (dabs(ccc) > 1d-9) then
              nn = nn+1
              do  i = 1, nr
                 rho2(i,ilm,isp) = rho2(i,ilm,isp) + ccc*rho1(i,jlm,isp)
              enddo
           endif
        enddo
     enddo
  enddo
end subroutine pysmr1
subroutine symprj(nrclas,nlml,ngrp,nbas,istab,g,ag,plat,qlat,pos, sym) !- Set up symmetry projectors for one class
  use m_ll,only:ll
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
  implicit none
  integer :: nlml,nrclas,ngrp,nbas,istab(nbas,ngrp)
  double precision :: sym(nlml,nlml,nrclas),plat(3,3),qlat(3,3), g(3,3,ngrp),ag(3,ngrp),pos(3,nrclas),d(3)
  integer:: lmxl , ig , ja , m , ia=99999 , iprint , ilm , l , jlm,  jlm1 , jlm2
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
        d(:) = matmul(g(:,:,ig),pos(:,ja)) + ag(:,ig) - pos(:,1)
        ia = ja
        rfrac = matmul(d,qlat)
        if(sum((rfrac-nint(rfrac))**2)< 1d-7) goto 80
12   enddo
     call rxi('symprj: no site mapped into first under op',ig)
80   continue
     ! ...   Make and add transformation matrix
     call ylmrtg ( nlml , g ( 1 , 1 , ig ) , rmat_rv )
     sym(:,:,ia)= sym(:,:,ia) + wgt*rmat_rv
10 enddo
  if( iprint() >= 60 ) then
     do ja = 1, nrclas
        write(6,"(/' projection matrices for ja=',i3)") ja
        ilm = 0
        do  l = 0, lmxl
           jlm1 = l*l+1
           jlm2 = (l+1)**2
           do  m = 1, 2*l+1
              ilm = ilm+1
              write(6,"(1x,9f12.6)") (sym(ilm,jlm,ja),jlm=jlm1,jlm2)
           enddo
        enddo
     enddo
  endif
  deallocate(rmat_rv)
end subroutine symprj
subroutine pxsmr1(wgt,nr,nlml,nsp,sym,rho1, rho2,nn)! Add wgt*sym*rho1 into rho2
  implicit none
  integer :: nr,nlml,nsp,nn,ilm,jlm,i,isp
  real(8) :: sym(nlml,nlml),rho1(nr,nlml,nsp), rho2(nr,nlml,nsp),wgt
  nn = 0
  do  ilm = 1, nlml
     do  jlm = 1, nlml
        nn = nn+1
        rho2(:,ilm,:) = rho2(:,ilm,:) + wgt*sym(ilm,jlm)*rho1(:,jlm,:)
     enddo
  enddo
end subroutine pxsmr1
