      subroutine symtab(nbas,iclass,bas,g,ag,ng,rb,qb,iwk,istab)
C- Make symmetry transformation table for basis atoms; check classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas:  number of atoms in basis
Ci   iclass:  jth atom belongs to class iclass(j)
Ci   bas:     bas(i,j) are Cartesian coordinates of jth atom in basis
Ci   g,ag,ng: defines space group (see sgroup.f)
Ci            NB: sign of ng used as a flag to suppress checking of
Ci            atom classes.
Ci   rb:  real space lattice vectors
Ci   qb:  reciprocal space lattice vectors
Ci   iwk: work array.  Not used if input ng<0
Co Outputs
Co   istab:  site ib is transformed into istab(ib,ig) by grp op ig
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,ng
      integer iclass(1),iwk(1),istab(nbas,1)
      double precision bas(3,1),g(9,1),ag(3,1),rb(9),qb(9)
C Local variables
      integer ib,ic,ig,jb,jc

      if (ng .eq. 0) return

C --- Make atom transformation table ---
      do  20  ig = 1, iabs(ng)
        do  10  ib = 1, nbas
C         call gpfndx(g(1,ig),ag(1,ig),ib,jb,bas,nbas,rb,qb)
          call grpfnd(1d-5,g,ag,ig,bas,nbas,qb,ib,jb)
          if (jb .eq. 0) then
            write(*,100) ib,ig
            call rx('SYMTAB: bad group operations')
          endif
          if (ng .gt. 0) then
            ic = iclass(ib)
            jc = iclass(jb)
            if (ic .ne. jc) then
              write(*,110) ib,jb,ig
              call rx('SYMTAB: invalid atom classes')
            endif
          endif
          istab(ib,ig) = jb
   10   continue
   20 continue

  100 format(/' ***ERROR*** Map of atom ib=',i4,' not found for ig=',i3)
  110 format(/' ***ERROR*** Atom ib=',i4,
     .  ' not in same class as mapped atom jb=',i4,' for ig=',i3)

      if (ng .lt. 0) return

C --- Check atom classes ---
      do  50  ib = 1, nbas
        ic = iclass(ib)
        call iinit(iwk,nbas)
        do  30  ig = 1, ng
          iwk(istab(ib,ig)) = 1
   30   continue
        do  40  jb = 1, nbas
          if (iwk(jb) .eq. 1) goto 40
          jc = iclass(jb)
          if (ic .eq. jc) then
            write(*,120) ib,jb
            call rx('SYMTAB: bad atom classes')
          endif
   40   continue
   50 continue

  120 format(/' ***ERROR*** Atom ib=',i4,
     .  ' in same class as symmetry-inequivalent atom jb=',i4)

      end

