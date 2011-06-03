      subroutine soprm(mode,lpzi,phi,phid,phiz,nr,nsp,lmxs,lmx,v,dv,enu,
     .ez,z,ri,wi,wk,sop,sopz)
C- Radial matrix elements between orbitals of different spin
C ---------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make spin-orbit parameters, i.e. matrix elements
Ci         :  <phi|so|phi> <phi|so|phidot> <phidot|so|phidot>
Ci         :  and for local orbitals that are present
Ci         :  <phiz|so|phiz> <phiz|so|phi> <phiz|so|phidot>
Ci         :2 make matrix elements for input magnetic field B
Ci         :  <phi|B|phi> <phi|B|phidot> <phidot|B|phidot>
Ci         :4 orthonormalize phi,phidot in separate spin channels
Ci         :5 1+4 above
Ci         :6 2+4 above
Ci   lpzi  :flags which channels have local orbitals
Ci   phi   :radial wave function * r
Ci   phid  :energy derivative of phi
Ci   phiz  :local orbital wave function
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmxs  :phi,phid,phiz,sop,sopz are dimensioned 0:lmxs
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   v     :electron part of potential: complete potential is v(r)-2*z/r.
Ci   dv    :(mode=1) radial derivative of potential, dV/dr with V=(V+ + V-)/2
Ci         :(mode=2) magnetic field
Ci         :(mode=4) not used
Ci   enu   :enu's for making charge density
Ci   ez    :enu's for local orbitals
Ci   z     :nuclear charge
Ci   ri    :radial mesh
Ci   wi    :weights for integration on mesh
Ci   wk    :work array of length nr*4
Co  Outputs
Co   sop   :sop(l,is1,is2,i=1..3) : matrix elements between orbitals
Co         :of spin is1 and is2 for quantum number l.
Co         :Three types of integrals are calculated for i=1..3:
Co         :<phi pert phi>  <phi pert phidot>  <phidot pert phidot>
Co   sopz  :sopz(l,is1,is2,i=1..3) : matrix elements between local orbitals
Co         :and orbitals of spin is1 and is2 for quantum number l.
Co         :Three types of integrals are calculated for i=1..3:
Co         :<phiz SO phiz>  <phiz SO phi>  <phiz SO phidot>.
Cr  Remarks
Cr   so = 2/(c^2) dV/dr*(1/r), V(r)=-2*z/r+v(r)
Cr   Note: so=1/(2*m^2*c^2)*(dV/dr*1/r), m=.5, c=274 (at. Rydberg units)
Cr   H_so = so*s^ dot l^, s^=0.5d0*sigma (Pauli matrix).
Cu Updates
Cu   17 Jan 07 Set ME for l=0 to zero (weren't necesssarily calc before)
Cu   11 Jul 05 Merged with sofp to make one routine
Cu   25 Apr 05 A. Chantis added local orbitals
Cu   05 Jan 04 leading dimensions of sop distinct from lmx
Cu   07 Feb 03 Added ability to compute matrix elements of external
Cu             field B.  New arg list and definition of mode.
C ---------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmx,mode,lpzi(0:lmx),lmxs,nr,nsp
      double precision z,
     .phi(nr,0:lmxs,nsp),phid(nr,0:lmxs,nsp),phiz(nr,0:lmxs,nsp),
     .ri(nr),sop(0:lmxs,nsp,nsp,3),v(nr,nsp),wk(nr,4),dv(nr),
     .wi(nr),sopz(0:lmxs,nsp,nsp,3),enu(0:8,nsp),ez(0:8,nsp)
C ... Local parameters
      integer l,ir,is,is1,is2,ipr,mode0,stdo,lgunit,lmin
      double precision c,pa,r,r1,r2,dot3,vavg,eavg,eavgz,dva,xx,
     .xxz,xxavg,wkz(nr,4)
C ... External calls
Ckino 2011.05 del external,       external daxpy,dscal,getpr
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c
      data pa /1d0/

C --- Setup ---
      call getpr(ipr)
      stdo = lgunit(1)
      mode0 = mod(mode,4)
C     c = 274.071979d0
      lmin = 1
      if (mode0 .eq. 2) lmin=0

C --- Orthonormalize phi, phidot, neglecting small component ---
      if (mode/4 .ne. 0) then
        if (ipr .gt. 50)
     .  print '(/'' soprm: overlaps  phi*phi     phi*phidot'')'
        do  1  is = 1, nsp
          do  1  l = 0, lmx
            r1 = dot3(nr,phi(1,l,is),phi(1,l,is),wi)
            call dscal(nr,1/dsqrt(r1),phi(1,l,is),1)
            r2 = dot3(nr,phi(1,l,is),phid(1,l,is),wi)
            call daxpy(nr,-r2,phi(1,l,is),1,phid(1,l,is),1)
            if (ipr .gt. 50) write(stdo,334) is,l,r1,r2
  334       format('  spin',i2,'  l=',i1,2f13.6)
    1   continue
      endif
      if (mode0 .eq. 0) return

C --- Matrix elements for each l ---
      do  is = 1, 4
        wk(1,is) = 0d0
        wkz(1,is) = 0d0
      enddo
      if (mode0 .eq. 2) then
        do  ir = 2, nr
          wk(ir,1) = wi(ir) * dv(ir)
        enddo
      endif

C .. Initialize matrix elements for s orbitals, in case not calculated
      l = 0
      do  is2 = 1, nsp
        do  is1 = 1, nsp
          do is = 1, 3
            sop(l,is1,is2,is) = 0d0
            if (lpzi(l) .ne. 0) then
              sopz(l,is1,is2,is) = 0d0
            endif
          enddo
        enddo
      enddo

      do  l = lmin, lmx
        eavg = (enu(l,1)+enu(l,nsp))/2
        if (lpzi(l) .ne. 0) then
          eavgz = (ez(l,1)+ez(l,nsp))/2
        endif
        do  is2 = 1, nsp
          do  is1 = 1, nsp
            if (mode0 .eq. 1) then
              do  ir = 2, nr
                r = ri(ir)
                vavg = (v(ir,1)+v(ir,nsp))/2 - 2*z/r
                dva  = dv(ir) + 2*z/r**2
                xx = 1/r/(1d0+pa*(eavg-vavg)/c**2)**2
                if (lpzi(l) .ne. 0) then
                  xxz = 1/r/(1d0+pa*(eavgz-vavg)/c**2)**2
                  xxavg = 0.5d0*(xx+xxz)
                  wkz(ir,1) = phiz(ir,l,is1)*dva
                  wkz(ir,2) = phiz(ir,l,is1)*xxz
                  wkz(ir,3) = phi(ir,l,is2)*xxavg
                  wkz(ir,4) = phid(ir,l,is2)*xxavg
                endif
                wk(ir,1) = phi(ir,l,is1)*dva
                wk(ir,3) = phid(ir,l,is1)*dva
                wk(ir,2) = phi(ir,l,is2)*xx
                wk(ir,4) = phid(ir,l,is2)*xx
              enddo
              sop(l,is1,is2,1) = dot3(nr,wk,wk(1,2),wi)*2d0/c**2
              sop(l,is1,is2,2) = dot3(nr,wk,wk(1,4),wi)*2d0/c**2
              sop(l,is1,is2,3) = dot3(nr,wk(1,3),wk(1,4),wi)*2d0/c**2
              if (lpzi(l) .ne. 0) then
                sopz(l,is1,is2,1) = dot3(nr,wkz,wkz(1,2),wi)*2d0/c**2
                sopz(l,is1,is2,2) = dot3(nr,wkz,wkz(1,3),wi)*2d0/c**2
                sopz(l,is1,is2,3) = dot3(nr,wkz,wkz(1,4),wi)*2d0/c**2
              endif
            elseif (mode0 .eq. 2) then
              sop(l,is1,is2,1) = dot3(nr,phi(1,l,is1),phi(1,l,is2),wk)
              sop(l,is1,is2,2) = dot3(nr,phi(1,l,is1),phid(1,l,is2),wk)
              sop(l,is1,is2,3) = dot3(nr,phid(1,l,is1),phid(1,l,is2),wk)
            else
              call rxi('soprm: bad mode:',mode)
            endif

          enddo
        enddo
      enddo

C --- Printout ---
      if (ipr .le. 50) return
      if (mode0 .eq. 1) write(stdo,332) 'spin-orbit coupling'
      if (mode0 .eq. 2) write(stdo,332) 'external field'
  332 format(' soprm:  matrix elements for perturbation from ',a/
     .13x,'l',4x,'<phi || phi>',2x,'<dot || phi>',2x,'<dot || dot>')
      if (nsp .eq. 1) then
        do  l = lmin, lmx
          write(stdo,333) '          ',
     .    l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
          if (lpzi(l) .ne. 0) then
            write(stdo,333) '          ',
     .      l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
          endif
        enddo
      else
        do  l = lmin, lmx
          write(stdo,333) 'up   up   ',
     .    l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
          write(stdo,333) 'down down ',
     .    l,sop(l,2,2,1),sop(l,2,2,2),sop(l,2,2,3)
          write(stdo,333) 'up   down ',
     .    l,sop(l,1,2,1),sop(l,1,2,2),sop(l,1,2,3)
          write(stdo,333) 'down up   ',
     .    l,sop(l,2,1,1),sop(l,2,1,2),sop(l,2,1,3)
          write(stdo,333)
          if (lpzi(l) .ne. 0) then
            write(stdo,335) 'up   up   ',
     .      l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
            write(stdo,335) 'down down ',
     .      l,sopz(l,2,2,1),sopz(l,2,2,2),sopz(l,2,2,3)
            write(stdo,335) 'up   down ',
     .      l,sopz(l,1,2,1),sopz(l,1,2,2),sopz(l,1,2,3)
            write(stdo,335) 'down up   ',
     .      l,sopz(l,2,1,1),sopz(l,2,1,2),sopz(l,2,1,3)
            write(stdo,335)
          endif
        enddo
      endif
C  333 format(1x,a,i3,3f20.15)
  333 format(1x,a,i3,1x, 3f14.8)
  335 format(1x,a,i3,'l',3f14.8)

      end


