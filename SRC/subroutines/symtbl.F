      subroutine symtbl(mode,tol,nbas,ipc,pos,g,ag,ng,qlat,istab)
C- Make symmetry transformation table for posis atoms; check classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1st digit
Ci         :0  site ib is transformed into istab(ib,ig) by grp op ig
Ci         :1  site istab(i,ig) is transformed into site i by grp op ig
Ci         :   (NB: old routine gpfndx used this convention)
Ci         :10s digit
Ci         :1  check atom classes
Ci   tol   :tol for which atoms are considered to be at the same site
Ci         :use 0 for symtbl to pick internal default
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   pos   :pos(i,j) are Cartesian coordinates of jth atom in basis
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   ng    :number of group operations
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Co Outputs
Co   istab :table of site permutations for each group op; see mode
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,ng,mode
      integer ipc(1),istab(nbas,1)
      double precision pos(3,1),g(9,1),ag(3,1),qlat(9),tol
C Local variables
      integer ib,ic,ig,jb,jc,mode1,mode10 !,oiwk
      double precision tol0,tol1
      character(200)::aaa
      parameter (tol0=1d-5)
C ... Heap
c      integer w(1)
c      common /w/ w
      integer,allocatable:: w_oiwk(:)

      if (ng .eq. 0) return
      mode1 = mod(mode,10)
      mode10 = mod(mode/10,10)
      tol1 = tol
      if (tol .eq. 0) tol1 = tol0

C     --- Make atom transformation table ---
      do  20  ig = 1, ng
        do  10  ib = 1, nbas
          call grpfnd(tol1,g,ag,ig,pos,nbas,qlat,ib,jb)
          if (jb .eq. 0) then
c     .    call fexit2(-1,111,' Exit -1 SYMTBL: no map for atom '//
c     .    'ib=%i, ig=%i',ib,ig)
            write(aaa,"('SYMTBL: no map for atom ib=',i0,' ig=',i0)") ib,ig
            call rx(aaa)
          endif
          if (mode10 .ne. 0) then
            ic = ipc(ib)
            jc = ipc(jb)
            if (ic .ne. jc) then
               write(aaa,"('SYMTBL: site ',i0,' not in same class as mapped site ',i0,', ig=',i0)")
     &              ib,jb,ig
               call rx(aaa)
c               call fexit3(-1,111,' Exit -1 SYMTBL: '//
c     .      'site %i not in same class as mapped site %i, ig=%i',
c     .      ib,jb,ig)
           endif    
          endif
          if (mode1 .eq. 0) then
            istab(ib,ig) = jb
          else
            istab(jb,ig) = ib
          endif
   10   continue
   20 continue

      if (mode10 .eq. 0) return

C --- Check atom classes ---
c      call defi(oiwk, nbas)
      allocate(w_oiwk(nbas))
      do  50  ib = 1, nbas
        ic = ipc(ib)
c        call iinit(w(oiwk),nbas)
        w_oiwk=0
        do  30  ig = 1, ng
c          w(oiwk+istab(ib,ig)-1) = 1
          w_oiwk(istab(ib,ig)) = 1
   30   continue
        do  40  jb = 1, nbas
c          if (w(oiwk+jb-1) .eq. 1) goto 40
          if (w_oiwk(jb) .eq. 1) goto 40
          jc = ipc(jb)
          if (ic .eq. jc) call fexit2(-1,111,' Exit -1 SYMTBL:  '//
     .    'site ib=%i in same class as inequivalent site jb=%i',ib,jb)
   40   continue
   50 continue
c      call rlse(oiwk)
      deallocate(w_oiwk)
      end

      subroutine istbpm(istab,nbas,ng,istab2)
C- Makes inverse of istab
C     implicit none
C Passed parameters
      integer nbas,ng
      integer istab(nbas,ng),istab2(nbas,ng)
C Local
      integer ib,ig,ibp

      do  ig = 1, ng
        do  ib = 1, nbas

          ibp = istab(ib,ig)
          istab2(ibp,ig) = ib
        enddo
      enddo

      end


