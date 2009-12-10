      subroutine pgfset(sspec,nbas,bas,plat,lorder,lrat,ips,vshft,
     .  pgfsl,pgfvl,pgord,pgfn,npl,npadl,npadr,pgplp)
C- Make arrays related to principal layers
C ----------------------------------------------------------------
Ci Inputs
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: eref name
Ci     Stored:    *
Ci     Passed to: *
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   plat  :principal lattice vectors
Ci   lorder:if .false., no final check is made on the ordering of sites
Ci         :within each PL.  Else, a check is made.  pgfset aborts if
Ci         :the sites are badly ordered.
Ci   lrat  :if T, rationalize pgfsl, numbering PL 0, 1, 2, ...
Ci   ips   :ips(ib) is species corresponding to site ib
Ci   vshft :(for printout only)
Ci   pfgsl :pgfsl(ib) is PL to which site ib corresponds
Ci   pgfvl :pgfsl(ib) is the PL potential group
Co Outputs
Co   npl   :number of principal layers
Co   npadl,npadr: number of padding sites L and R
Co   pgord: table of indices that sort basis (dimensioned 2*nbas)
Co   pgfn:  projection normal to transverse vectors, increasing with PL.
Co   pgplp: PL-dependent indices.  pgfset returns:
Co          1: cumulative number of basis atoms in this PL=0 ... this PL
Co          2: index to PL with inequivalent potential
Co          (the following are not made here)
Co          3: source (column) dimension of GF for this PL
Co          4: row (field) matrix dimension for this PL
Co          5: matrix dimension for this PL including i-waves
Co          6: offset to diagonal part of g
Co          PL are numbered 0...npl-1 here
Cio Inputs/Outputs
Cr Local variables
Cl   pgord  :permutation of sites that orders them by increasing PL.
Cl          :For a given PL, original order is preserved.
Cl   pgfn   :vector normal to plane Plat(1) x Plat(2)
Cr Remarks
Cr   The PL must be ordered so that their projection along plat(3)
Cr   increases with increasing PL, and the bulk left and right layers
Cr   have the smallest and largest projections.
Cr   pgfset does the following:
Cr     1. Finds the permutation pgord
Cr     2. Determines pgfn = vector normal to the plane P1xP2, 
Cr        and the rightmost projection onto pgfn in the L layer
Cr     3. Shift all PL that are ordered improperly relative to layer L
Cr     4. Optionally check that the set of PL are properly ordered
Cu Updates
Cu   09 Aug 05 Better check of sites ordering within PL
Cu   20 Sep 04 Better printout when PL are badly ordered
C ----------------------------------------------------------------
C     implicit none
      logical lorder,lrat
      integer nbas,npl,npadl,npadr,pgfsl(nbas),pgfvl(nbas),
     .  pgord(nbas,2),pgplp(6,0:1),ips(nbas)
      double precision vshft(-7:nbas),sspec(1)
C local variables
      logical flag
      integer i,j,ib,il,jl,jb,ipr,iv,nv,ib1,isw,isum,ibs,
     .  ils,ilL,stdo,nglob,ibold,ilold
      double precision bas(3,nbas),pgfn(3),plat(3,3),
     .  xx(4),xl,xpl,xr,xpr,x0,ddot,xold,xls,xlmax
      equivalence (xx(1),xl), (xx(2),xpl), (xx(3),xpr), (xx(4),xr)
      character*80 outs, fmt*20, spid*8

      call getpr(ipr)
      stdo = nglob('stdo')
      if (ipr .ge. 20) print '(1x)'

C --- Sort basis by increasing PL, preserving order for constant PL ---
      call ivshel(1,nbas,pgfsl,pgord,.true.)
      do  4  i = 1, nbas
    4 pgord(i,1) = pgord(i,1)+1
      j = 1
      jl = pgfsl(pgord(1,1))
      do  5  i = 1, nbas
        ib = pgord(i,1)
        il = pgfsl(ib)
        if (il .ne. jl) then
          call ivshel(1,i-j,pgord(j,1),pgord(1,2),.false.)
C         call awrit1('%32:1i',' ',80,stdo,pgord)
          j  = i
          jl = il
        endif
    5 continue
      call ivshel(1,nbas+1-j,pgord(j,1),pgord(1,2),.false.)

C --- Vector pgfn, and xl=Rightmost proj. in L layer ---
      call cross(plat,plat(1,2),pgfn)
      x0 = ddot(3,pgfn,1,plat(1,3),1)
      if (x0 .lt. 0) call dscal(3,-1d0,pgfn,1)
      x0 = ddot(3,pgfn,1,pgfn,1)
      call dscal(3,1/sqrt(x0),pgfn,1)
      xl = -9999d0
      ilL = pgfsl(pgord(1,1))
      do  i = 1, nbas
        ib = pgord(i,1)
        il = pgfsl(ib)
        if (il .eq. ilL) then
          xl = max(xl,ddot(3,pgfn,1,bas(1,ib),1))
        endif
      enddo

C --- Add plat(3) to any layer not to the right of the L layer ---
C     call prmx('starting pos',bas,3,3,nbas)
      flag = .false.
      do  jl = ilL+1, pgfsl(pgord(nbas,1))
        xr =  9999d0
        do  i = 1, nbas
          ib = pgord(i,1)
          il = pgfsl(ib)
          if (il .eq. jl) then
            if (ddot(3,pgfn,1,bas(1,ib),1) .lt. xr) then
              xr = ddot(3,pgfn,1,bas(1,ib),1)
              jb = ib
            endif
          endif
        enddo

        if (xr .lt. xl) then
C         First-time printout if needed
          if (.not. flag) then
            call info2(1,0,0,' PGFSET (warning): PL %i%-1j '//
     .        'extends beyond sites in other PL:'//
     .        '%N max P1xP2 for PL %i is %,6d ... add plat(3) for '//
     .        'these layers:%N'//'   PL  min P1xP2  ib   New P1xP2',
     .        ilL,xl)
            flag = .true.
          endif
          do  i = 1, nbas
            ib = pgord(i,1)
            il = pgfsl(ib)
            if (il .eq. jl) then
              call daxpy(3,1d0,plat(1,3),1,bas(1,ib),1)
            endif
          enddo
          xold = xr
          xr = ddot(3,pgfn,1,bas(1,jb),1)
          if (ipr .ne. 0) write(stdo,333) jl,xold,jb,xr
  333     format(i4,f12.6,i4,f12.6)
        endif
      enddo
C     call prmx('pos after shift bulk',bas,3,3,nbas)

C --- Check for increasing projection with increasing PL  ---
      flag = .false.
      ibs = pgord(1,1)
      ils = pgfsl(ibs)
      xold  = -1d10
      ibold = 0
      ilold = 0
      xlmax = -1d10
      jl = pgfsl(pgord(1,1))
      npl = 1
      pgplp(1,npl) = 0
      pgplp(2,npl) = pgfvl(pgord(1,1))
      do  20  i = 1, nbas
        ib = pgord(i,1)
        il = pgfsl(ib)
        xl = ddot(3,pgfn,1,bas(1,ib),1)
C       Stay within a layer ... update xlmax
        if (il .eq. jl) then
          pgplp(1,npl) = pgplp(1,npl)+1
          xls = xlmax
          xlmax = max(xlmax,xl)
C         for error tracking only ... save site index that changes xlmax
          if (xlmax .ne. xls) ibs = ib
          if (xlmax .ne. xls) ils = il
C       New layer: update xold and jl
        elseif (il .gt. jl) then
          xold = xlmax
          ibold = ibs
          ilold = ils
          xlmax = xl
          npl = npl+1
          pgplp(1,npl) = 1
          pgplp(2,npl) = pgfvl(ib)
          jl = il
        else
          call rx('bug in pgfset')
        endif
        if (lorder .and. xl .lt. xold) then
          call info8(10,1,0,' PGFSET:  program will abort ...%N'
     .    //' site %i (in PL %i) has projection along P1 x P2 = %;4d%N'
     .    //' site %i (in PL %i) has projection along P1 x P2 = %;4d',
     .      ibold,ilold,xold,ib,il,xl,0,0)
          flag = .true.
        endif
        if (lrat) pgfsl(ib) = npl
   20 continue
      npadl = pgplp(1,1)
      npadr = pgplp(1,npl)

C --- Printout ---
      if (ipr .ge. 30) then
        call awrit5(' PGFSET: %i principal layers, %i+%i+%i sites.'//
     .    '  Normal:%3:1;5d',' ',80,stdo,npl,nbas,npadl,npadr,pgfn)
        j = 18
        fmt = '(1x,a,20I4)'
        do  38  i = 0, npl-1, j
          if (i .gt. 0) print *, ' '
          print fmt,'  PL |',(ib, ib = i, min(j+i-1,npl-1))
          print fmt,'size |',(pgplp(1,ib+1), ib = i, min(j+i-1,npl-1))
          print fmt,' PLV |',(pgplp(2,ib+1), ib = i, min(j+i-1,npl-1))
   38   continue
        do  34  i = 2, nbas
   34   if (pgord(i,1) .le. pgord(i-1,1)) goto 36
        goto 37
   36   if (lrat) print '('' Basis reordered by increasing PL:'')'
   37   continue
C       j = 18
C       fmt = '(1x,a,20I4)'
        if (nbas .lt. 100) fmt = '(1x,a,26I3)'
        if (ipr .ge. 40) then
            print 334
  334       format(' New ib Old   Spec      PL  PLV',15x,'Pos',
     .        16x,'Pos.h',4x,'Vshift')
          il = -1
          do  33  jb = 1, nbas
            ib = pgord(jb,1)
            outs = ' '
            call spacks(0,'spec name',sspec,spid,ips(ib),ips(ib))
            write(outs,335) jb, ib, spid, pgfsl(ib)-1,
     .        pgfvl(ib), (bas(ib1,ib), ib1=1,3),
     .        ddot(3,pgfn,1,bas(1,ib),1)
            call awrit3('%?#n#%23p        ##%a%?#n#%;10,6D##',outs,
     .        len(outs),-stdo,isw(pgfsl(ib).eq.il),
     .        vshft(ib),vshft(ib))
  335       format(i4,i6,4x,a8,2i4,1x,3f10.5,f10.5,f10.6)
            il = pgfsl(ib)
   33     continue
        else
          do  30  jb = 1, nbas, j
            print fmt, '  PL',
     .      (pgfsl(pgord(ib,1)), ib = 1+jb-1, min(j+jb-1,nbas))
            print fmt, 'Site',
     .      (pgord(ib,1), ib = 1+jb-1, min(j+jb-1,nbas))
   30     continue
        endif
      endif

      if (flag) call fexit(-1,111,' PGFSET: badly ordered sites',0)

C --- Convert pgplp(1) into accumulated number.  isum avoids compiler bug
      do  40  jl = npl, 1, -1
   40 pgplp(1,jl) =  isum(jl,pgplp(1,1),6)
      pgplp(1,npl+1) = pgplp(1,npl) + npadr
C     pgplp(1,-1) = 0
C      pgplp(1,npl+2) = pgplp(1,npl+1) + npadr
C      pgplp(1,npl+2) = pgplp(1,npl+1) + npadr
C      print '(a,20i3)','size |',(pgplp(1,jl), jl = 0, npl+1)
C      stop

C --- Map pgfvl(isite) into pgplp(2,PL) ---
      if (lrat) then
        il = 0
        iv = pgfvl(pgord(1,1))
C ...   Check that all sites with a given pgfsl have constant pgfvl
        do  51  jl = 1, nbas
          ib = pgord(jl,1)
C         print *, pgfsl(ib),il, pgfvl(ib),iv
          if (pgfsl(ib) .eq. il .and. pgfvl(ib) .ne. iv)
     .      call fexit3(-1,111,' Exit -1 PGFSET: improper pgfvl(%i):'//
     .      ' expected %i but found %i',ib,iv,pgfvl(ib))
          il = pgfsl(ib)
          iv = pgfvl(ib)
   51   continue
        il = -1
        nv = 0
        do  52  jl = 1, nbas
          ib = pgord(jl,1)
          if (pgfsl(ib) .ne. il) then
            nv = nv+1
C           pgplp(2,nv) = pgfvl(ib)
          endif
          il = pgfsl(ib)
   52   continue
      endif

C     call prmx('end of pgfset, pos',bas,3,3,nbas)

      end

