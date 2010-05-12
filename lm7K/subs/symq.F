      subroutine psymr0(lmxl,ic,nbas,ipc,pos0,pos,ipa,nrclas)
C- Make list of atoms in one class
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmxl  :lmxl(1..nbas) = l-cutoffs dimensioning offsets.
Ci         :Used in creating table offsets for table ipa, which see
Ci         :lmxl(1) < -1 => ipa is created differently and
Ci         :lmxl(2..) are not referenced; see description for ipa.
Ci   ic    :abs. val of ic = class for which to make list of members
Ci         :sign of ic<0 : flag to suppress generation of pos.
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   pos0  :basis vectors for each site (not used if ic<0)
Co Outputs
Co   ipa   :ipa(1..nrclas) is created in one of two modes:
Co         :if lmxl(1)<-1 :
Co         : ipa = list of sites belonging to class ic
Co         :       In this mode lmxl is not used.
Co         :if lmxl(1)>=-1 :
Co         : ipa = table of offsets to an array for each site in class
Co         :       ic.  It is intended for some array A containing
Co         :       elements for each site 1..nbas, and that offset
Co         :       to the A for site ib+1, offs(ib+1) is
Co         :       (lmxl(ib)+1)**2 + offs(ib).
Co   pos   :basis vectors for each member of class (made if ic>0)
Co   nrclas:number of sites in class ic
Cu Updates
Cu   17 May 01 Add ic<0 as flag to suppress generation of pos
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,lmxl(nbas),ic,nrclas,ipc(nbas),ipa(nbas)
      double precision pos0(3,nbas),pos(3,nbas)
C ... Local parameters
      logical loff
      integer ib,iclass,ioff
      logical lpos

      loff = lmxl(1) .ge. -1
      iclass = iabs(ic)
      lpos = ic .gt. 0
      nrclas = 0
      ioff = 0
      do  ib = 1, nbas
        if (ipc(ib) .eq. iclass) then
          nrclas = nrclas+1
          if (.not. loff) ioff = ib
          ipa(nrclas) = ioff
          if (lpos) then
            pos(1,nrclas) = pos0(1,ib)
            pos(2,nrclas) = pos0(2,ib)
            pos(3,nrclas) = pos0(3,ib)
          endif
        endif
        ioff = ioff + (lmxl(ib)+1)**2
      enddo
      end

      subroutine psymq0(nrclas,nsp,ipa,wk,nvec,vec)
C- Symmetrize vector of objects with l=0 symmetry for one class of atoms
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrclas:number of atoms in the present class
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co   ipa   :list of site indices for each member of class
Ci   wk    :work array of length nvec*nsp
Ci   nvec  :size of vector for each member of class
Co Inputs/Outputs
Cio   vec  :On input,  unsymmetrized vector of quantities
Cio        :On output, elements of vec in class list ipa(1..nrclas)
Cio        :are symmetrized.
Cu Updates
Cu   17 May 01
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nrclas,nsp,nvec,ipa(nrclas)
      double precision wk(nvec,nsp),vec(nvec,nsp,1)
C ... Local parameters
      integer ia,ib,isp,i
      double precision fac

      call dpzero(wk,nsp*nvec)
      fac = 1d0/nrclas
      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
          do  i = 1, nvec
            wk(i,isp) = wk(i,isp) + vec(i,isp,ib)*fac
          enddo
        enddo
      enddo

      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
          do  i = 1, nvec
            vec(i,isp,ib) = wk(i,isp)
          enddo
        enddo
      enddo

      end

      subroutine symqmp(nrclas,nlml,nlmx,plat,posc,ngrp,g,ag,qwk,
     .  ipa,sym,qmp,nn)
C- Symmetrize multipole moments for a single class
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrclas:number of atoms in the ith class
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nlmx  :dimensions sym: sym is generated for ilm=1..nlmx
Ci   plat  :primitive lattice vectors, in units of alat
Ci   posc  :work array holding basis vectors for this class
Ci   ngrp  :number of group operations
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   qwk   :work array of dimension nlml
Ci   ipa   :ipa(1..nrclas) = table of offsets to qmp corresponding
Ci         :to each member of the class
Cio Inputs/Outputs
Cio  qmp   :On input,  unsymmetrized multipole moments
Cio        :On output, multipole moments are symmetrized
Co Outputs
Co   sym   :symmetry projectors for each member of the class
Co   nn    :number of elements symmetrized
Cu Updates
Cu   23 Aug 01 adapted from psymql
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nrclas,nlmx,ipa(nrclas),nlml,ngrp,nn
      double precision plat(3,3),posc(3,nrclas),g(3,3,ngrp),ag(3,ngrp)
      double precision sym(nlmx,nlmx,nrclas),qwk(nlml),qmp(*)
C ... Local parameters
      integer ia,ilm
      double precision wgt,xx,qlat(3,3)

      if (nlml .gt. nlmx) call rxi('symqmp: increase nlmx to',nlml)
      call dinv33(plat,1,qlat,xx)

C ... Make the symmetry projectors
      call symprj(nrclas,nlmx,ngrp,xx,xx,g,ag,plat,qlat,posc,sym)

C ... Accumulate symmetrized qmpol on first site
      call dpzero(qwk, nlml)
      do  ia = 1, nrclas
        call pxsmr1(1d0,1,nlml,1,sym(1,1,ia),qmp(1+ipa(ia)),qwk,nn)
      enddo

C ... Rotate and copy to all sites in class
      wgt = nrclas
      do  ia = 1, nrclas
        call dpzero(qmp(1+ipa(ia)), nlml)
        call pysmr1(wgt,1,nlml,1,sym(1,1,ia),qwk,qmp(1+ipa(ia)),nn)
      enddo

      nn = 0
      do  ilm = 1, nlml
        if (dabs(qmp(ilm+ipa(1))) .gt. 1d-6) nn = nn+1
      enddo

      end

      subroutine pxsmr1(wgt,nr,nlml,nsp,sym,rho1,rho2,nn)
C- Add wgt*sym*rho1 into rho2
C     implicit none
C ... Passed parameters
      integer nr,nlml,nsp,nn
      double precision sym(nlml,nlml),rho1(nr,nlml,nsp),
     .  rho2(nr,nlml,nsp),wgt
C ... Local parameters
      integer ilm,jlm,i,isp
      double precision ccc

      nn = 0
      do  isp = 1, nsp
        do  ilm = 1, nlml
          do  jlm = 1, nlml
            ccc = wgt*sym(ilm,jlm)
            if (dabs(ccc) .gt. 1d-9) then
              nn = nn+1
              do  i = 1, nr
                rho2(i,ilm,isp) = rho2(i,ilm,isp) + ccc*rho1(i,jlm,isp)
              enddo
            endif
          enddo
        enddo
      enddo
      end

      subroutine pysmr1(wgt,nr,nlml,nsp,sym,rho1,rho2,nn)
C- Add wgt*sym*rho1 into rho2, transposed
C     implicit none
C ... Passed parameters
      integer nr,nlml,nsp,nn
      double precision sym(nlml,nlml),rho1(nr,nlml,nsp),
     .  rho2(nr,nlml,nsp),wgt
C ... Local parameters
      integer i,ilm,jlm,isp
      double precision ccc

      nn = 0
      do  isp = 1, nsp
        do  ilm = 1, nlml
          do  jlm = 1, nlml
            ccc = wgt*sym(jlm,ilm)
            if (dabs(ccc) .gt. 1d-9) then
              nn = nn+1
              do  i = 1, nr
                rho2(i,ilm,isp) = rho2(i,ilm,isp) + ccc*rho1(i,jlm,isp)
              enddo
            endif
          enddo
        enddo
      enddo
      end

