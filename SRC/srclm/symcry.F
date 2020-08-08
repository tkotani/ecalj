      subroutine symcry(tol,bas,bast,ipc,nbas,nclass,nrclas,
     .ng,plat,qlat,g,ag,istab)
C- Generates the symmetry ops of the crystal from those of the lattice
C ----------------------------------------------------------------------
Ci Inputs:
Ci   tol:   tol for which atoms are considered to be at the same site
Ci          use 0 for symtbl to pick internal default
Ci   bas   :basis vectors (scaled by alat)
Ci   ipc   :the jth atom belongs to class ipc(j)
Ci   nbas  :number of atoms in the basis
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   nrclas:number of atoms in the ith class
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   bast  :work array of same dimension as bas
Cio Inputs/Outputs:
Cio  ng    :number of allowed symmetry operations (see Remarks)
Cio         on input  number of symmetry operations of the lattice
Cio         on output number of symmetry operations of the crystal
Cio  g     :symmetry operation matrices
Co Outputs:
Co   ag    :symmetry operation vector
Co   istab :site ib is transformed into istab(ib,ig) by operation ig
Cr Remarks:
Cr   symcry finds the subset of the allowed ng point operations of the
Cr   lattice without a basis (see symlat.f) that are valid for
Cr   the crystal.
Cr
Cr   This routine is based on ATFTMT written by Worlton and Warren,
Cr   CPC 3, 88 (1972).
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer nbas,ng,ipc(nbas),nclass,istab(nbas,ng),nrclas(nclass)
      double precision plat(9,*),qlat(*),bas(3,*),bast(3,*),
     .g(9,*),ag(3,*),tol
C Local parameters:
      integer ibas,ic,iclbsj,icmin,ig,ipr,jbas,kbas,kc,
     .lgunit,m,mbas,nj,nm,ng0,stdo
C     integer mode(3)
      double precision dbas(3),tol0,tol1
      parameter (tol0=1d-5)
      logical latvec
      character sg*35

      call getpr(ipr)
      stdo = lgunit(1)
      tol1 = tol
      if (tol .eq. 0) tol1 = tol0

C --- Find the class with minimum number of atoms ---
      icmin = 1
      do  5  ic = 1, nclass
        if (nrclas(ic).lt.nrclas(icmin) .and. nrclas(ic).gt.0) icmin = ic
    5 continue
      ibas = iclbsj(icmin,ipc,nbas,1)

C --- For each group op, see whether it only shifts basis by some T ---
      ng0 = ng
      ng = 0
      do  30  ig = 1, ng0
C   ... Rotate the basis by g
        call dmpy(g(1,ig),3,1,bas,3,1,bast,3,1,3,nbas,3)
        do  20  nj = 1, nrclas(icmin)
          jbas = iclbsj(icmin,ipc,nbas,nj)
C     ... This is a candidate for translation ag
          do  22  m = 1, 3
            ag(m,ng+1) = bas(m,jbas)-bast(m,ibas)
   22     continue
          call shorbz(ag(1,ng+1),ag(1,ng+1),plat,qlat)
C          mode(1) = 2
C          mode(2) = 2
C          mode(3) = 2
C          call shorps(1,plat,mode,ag(1,ng+1),ag(1,ng+1))
C     ... See whether candidate works for all sites; also make istab
          do  10  kbas = 1, nbas
            kc = ipc(kbas)
            do  12  nm = 1, nrclas(kc)
              mbas = iclbsj(kc,ipc,nbas,nm)
              do  14  m = 1,3
                dbas(m) = bas(m,mbas)-bast(m,kbas)-ag(m,ng+1)
   14         continue
              if (latvec(1,tol1,qlat,dbas)) then
                istab(kbas,ng+1) = mbas
                goto 10
              endif
   12       continue
C       ... Candidate not valid
            if (ipr .ge. 90) then
              call asymop(g(1,ig),ag(1,ng+1),' ',sg)
              call awrit1(' symcry: excluded candidate ig=%,2i  '//sg
     .        //'%a',' ',80,stdo,ig)
            endif
            goto 20
   10     continue

C     --- Valid ag found; add g to list ---
          ng = ng+1
          if (ig .gt. ng) call dcopy(9,g(1,ig),1,g(1,ng),1)
          if (ipr .ge. 70) then
            call asymop(g(1,ng),ag(1,ng),' ',sg)
            call awrit1(' symcry: accepted candidate ig=%,2i  '//sg
     .      //'%a',' ',80,stdo,ig)
          endif
          goto 30
   20   continue
   30 continue

      if (ipr .ge. 30) call awrit2(' SYMCRY: crystal invariant under '//
     .'%i symmetry operations for tol=%;3g',' ',80,stdo,ng,tol1)
      if (ipr .ge. 60 .and. ng .gt. 1) then
        write(stdo,'('' ig  group op'')')
        do  60  ig = 1, ng
          call asymop(g(1,ig),ag(1,ig),' ',sg)
          write(stdo,'(i4,2x,a)') ig,sg
   60   continue
      endif

      end

