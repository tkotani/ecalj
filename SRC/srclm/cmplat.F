      integer function cmplat(nbas,plat,plat2,ng,g,ips,pos,ips2,pos2)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   plat2 :primitive lattice vectors of second lattice (not used now)
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   pos   :basis vectors
Ci   ips2  :species table of second lattice
Ci   pos2  :basis vectors of second lattice
Co Outputs
Co   cmplat:-1 species mismatch
Co         : 1 lattices have equivalent numbers of species, but
Co         :   no symmetry operation was found that aligns bases
Co         : 0 a symop was found that aligns the basis.
Cl Local variables
Cl   nks   :number of sites belonging to 'special' species
Cl         :used to find set of possible translations
Cl   nrs   :number of sites belonging to current species
Cr Remarks
Cr   At present, it is assumed both bases have the same lattice vectors.
Cu Updates
Cu   12 Sep 01 Originally created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,ips(nbas),ips2(nbas),ng
      double precision plat(3,3),plat2(3,3),pos(3,nbas),pos2(3,nbas),
     .g(9,ng)
C ... Local parameters
      logical latvec
      integer nspec,mxint,ipr,ks,nks,is,nrs,ns2,stdo,lgunit
      integer ipa(nbas),ipa2(nbas),iprm(nbas),ig,k,i,j
      double precision xx,posa(3,nbas),posa2(3,nbas),
     .qlat(3,3),trans(3),dbas(3),tol0,tol1,tol
      parameter (tol0=1d-5)

      call getpr(ipr)
      stdo = lgunit(1)
      call dinv33(plat,1,qlat,xx)
      tol = 0
      tol1 = 1d-6
      if (tol .eq. 0) tol1 = tol0

C ... Count the number of species; check for match
      cmplat = -1
      nspec = mxint(nbas,ips)
      if (nspec .ne. mxint(nbas,ips2)) then
        call info(20,0,0,
     .  ' ... species mismatch; lattices are not equivalent',0,0)
        return
      endif
C ... Check for eq. no. sites for each species;
C     Find species ks holding the fewest sites.
      ks = 0
      nks = 0
      do  is = 1, nspec
        call psymr0(-2,-is,nbas,ips,xx,xx,ipa,nrs)
        call psymr0(-2,-is,nbas,ips2,xx,xx,ipa,ns2)
        if (nrs .ne. ns2) then
          call awrit1(' ... basis differ in the number of species of'//
     .    ' type %i ... lattices inequivalent',' ',80,stdo,is)
          return
        endif
        if (nrs .lt. nks .or. nks .eq. 0) then
          ks = is
          nks = nrs
        endif
      enddo

C --- Loop over all possible (rotations+translations) ---
C     The set of allowed translations is obtained by matching one
C     of the posa2(spec ks), after rotation by ig, to the first
C     of the posa(spec ks).  A match is found if a
C     (rotation+translation) fits all species
      do  ig = 1, ng
C     write(6,'(3f12.6)') g(1:9,ig)
        do  k = 1, nks

C  ...  Get translation for this k
          call psymr0(-2,ks,nbas,ips,pos,posa,ipa,nks)
          call psymr0(-2,ks,nbas,ips2,pos2,posa2,ipa2,nks)
C       Rotate posa2(k) by g(ig)
          call dmpy(g(1,ig),3,1,posa2(1,k),3,1,trans,3,1,3,1,3)
C        call shorbz(trans,trans,plat,qlat)
C        if (k .eq. 2) print '(3f12.6)', trans
C       Subtract posa(1)
          call daxpy(3,-1d0,posa,1,trans,1)

C   ... See whether this translation is valid for all species
          do  is = 1, nspec
C         Next two lines create posa2(is), rotated by g(ig)
            call psymr0(-2,is,nbas,ips2,pos2,posa,ipa2,nrs)
            call dmpy(g(1,ig),3,1,posa,3,1,posa2,3,1,3,nrs,3)
C         Create posa(is)
            call psymr0(-2,is,nbas,ips,pos,posa,ipa,nrs)

C         Try to align by considering every possible pair.
C         Inefficient, but easy to program.
C         If there is a match, all iprm will be nonzero
            call iinit(iprm,nrs)
            do  i = 1, nrs
              do  j = 1, nrs
                dbas(1) = posa2(1,j) - posa(1,i) - trans(1)
                dbas(2) = posa2(2,j) - posa(2,i) - trans(2)
                dbas(3) = posa2(3,j) - posa(3,i) - trans(3)
                if (latvec(1,tol1,qlat,dbas)) then
                  if (iprm(i) .ne. 0) call fexit2(-1,111,' Exit -1 '//
     .            'CMPLAT: site %i apparently coincident with site %i',
     .            j,iprm(i))
                  iprm(i) = j
                  goto 20
                endif
              enddo
   20         continue
            enddo

C     ... No match unless all iprm are nonzero
            do  i = 1, nrs
              if (iprm(i) .eq. 0) then
                if (ipr .ge. 60) write(stdo,333) ig,trans,is
  333           format('  try grp op',i3,' trans =', 3f11.6,
     .          ' ... mismatched spec',i3)
                goto 10
              endif
            enddo

C       A match was found for this species ... try next species
          enddo
C       A match was found for all species ... return with match info
          call shorbz(trans,trans,plat,qlat)
          cmplat = 0
          call awrit2(' ... match for group op %i, translation%'//
     .    '3:1,6;6d',' ',80,stdo,ig,trans)
          return

C       No match found; try next translation
   10     continue
        enddo
      enddo

   99 continue
      cmplat = 1
      call awrit2(' ... no match found (%i group ops,'//
     .' %i translations)',' ',80,stdo,ng,nks)
      end

