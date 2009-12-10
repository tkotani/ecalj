      subroutine mksym(mode,slabl,ssymgr,sctrl,slat,ssite,sarray)
C- Setup for symmetry group
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : 1s digit
Ci           0  make space group only
Ci           1  add inversion to point group (see aginv below)
Ci           2  Same as 1, but make additionally ag,istab for extra
Ci              operations, using -g for rotation part; see Remarks
Ci           10s digit
Ci           0  do nothing about classes
Ci           1  split species into classes
Ci           2  Also assign class labels to each class
Ci           4  Assume class and species are the same.
Ci   slabl : species labels
Ci   ssymgr: string containing symmetry group generators.
Ci           if ssymgr contains 'find', mksym will add basis atoms as
Ci           needed to guarantee generators are valid, and generate
Ci           internally any additonal group operations needed to
Ci           complete the space group.
Ci
Cio Inputs/outputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nspec nsite lpgf
Ci     Stored:    nclass
Ci     Passed to: *
Cio  slat  :struct containing information about the lattice
Cio    Elts read: plat ldist dist
Cio    Stored:    npgrp nsgrp oistab oag osymgr
Cio          :nsgrp= number of space group operations.
Cio          :npgrp= number of point group operations.  It may differ
Cio                  from nsgrp if inversion is added.
Cio          :istab= table of site permutations for each group op
Cio          :ag   = translation part the group ops
Cio          :symgr= point group (rotation) part of each group op
Cio  ssite :struct containing site-specific information
Cio    Elts read: *
Ci     Stored:    class clabel
Cio  sarray:structure containing offsets to various arrays
Cio    Elts read: opos oips npadl npadr
Cio    Stored:    onrcp oipcp oics oipc oclabl nclasp onrc
Co           oipc,oipcp,onrc,oics and oclabl, and
Co           oipc:  pointer to class table, pad equiv to double pad
Co           oipcp: pointer to class table, pad nequiv to double pad
Co Outputs
Co  ssymgr:  may be enlarged if it contains 'find' initially
Co  slat     the following are set: npgrp nsgrp oistab oag osymgr
Cr Remarks
Cr   In certain cases the inversion operation may be added to the space
Cr   group, for purposes of k integration.  This is permissible when the
Cr   hamiltonian has the form h(-k) = h*(k).  In that case, the
Cr   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
Cr
Cr   Also, the Green's functions are related G(-k) = Gtranspose(k).
Cr   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
Cr   then G(-k) = Gtranspose(k), and the same (g,ag) information is
Cr   needed for either rotation.
Cu Updates
Cu   04 Jan 06 Enabled keyword 'point' in ssymgr
Cu             Generates extra group ops when generators create more 
Cu             than the maximum expected
Cu   07 Oct 05 Put altered positions into site structure
Cu   03 Nov 01 Adapted to altered gensym
Cu   26 Jan 01 Added new '2' option to 1s digit of mode
Cu   31 May 00 revised meaning of mode
C ----------------------------------------------------------------------
C     implicit none
      character*8 slabl(1),ssymgr*(*)
      integer mode,nsgrp,npgrp
      double precision sctrl(1),slat(1),ssite(1),sarray(1)
C ... Heap
      integer w(1)
      common /w/ w
C Local variables
      logical T,F,cmdopt,a2bin,ltmp
      integer idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nbas,
     .  nbas0,nbasp,nsite,nclass,nclasp,nclspp,ngen,ngnmx,nspec,usegen,
     .  npadl,npadr,ldist,nggen,ngmx,stdo,nglob,incli
      integer oclabl,oics,oips,oistab,onrc,onrcp,
     .  onrspc,opos,oipc,oipcp,oiwk,oag,osymgr,opos2,oips2,aginv
      parameter (T=.true., F=.false., ngnmx=10)
      character*120 gens,strn*72
      double precision gen(9,ngnmx),plat(3,3),dist(3,3),qlat(3,3),xx,
     .  fptol

      call upack('ctrl nbas nspec nsite',sctrl,nbas,nspec,nsite,0,0)
      call upack('lat plat ldist dist',slat,plat,ldist,dist,0,0)
      ngmx = 48
      stdo = nglob('stdo')
      if (iprint() .ge. 30) write(stdo,'(1x)')

C ... Re-entry when ngmx was increased
    5 continue

C ... For memory allocation of class data, assume initially nclass=nbas
      call defrr(oag,   3*ngmx)
      call defrr(osymgr,9*ngmx)
      call defi(oipc,  nsite)
      call defi(oics,  nsite)
      call defi(onrspc,nbas)
      call upack('array opos oips',sarray,opos,oips,0,0,0)
C     call defi(oips,nbas)
C     call defrr(opos,3*nbas)
C     call spackv(10,'site spec',ssite,1,nbas,w(oips))
C     call spackv(10,'site pos',ssite,1,nbas,w(opos))

*     call prmx('basis',w(opos),3,3,nbas)
*     print *, (w(oips-1+i), i=1,nbas)

C ... Extract keywords from input string; copy rest to gens
C     keyword 'point' is appended to gens
      call words(ssymgr,ngen)
      j1 = 1
      idest = 1
      usegen = 2
      gens = ' '
      ltmp = .false.
      do  ig = 1, ngen
        call word(ssymgr,ig,j1,j2)
        if (ssymgr(j1:j2) .eq. 'find') then
          usegen = 0
        else if (ssymgr(j1:j2) .eq. 'point') then
          ltmp = .true.
        else
          call strncp(gens,ssymgr,idest,j1,j2-j1+2)
          idest = idest+j2-j1+2
        endif
      enddo
      if (ltmp) then
        call strncp(gens,'point',idest,1,5)
      endif

C --- Generate space group ---
      nbas0 = nbas
      if (cmdopt('--fixpos',8,0,strn)) then
        j1 = 8+1
        if (strn(9:13) .eq. ':tol=') then
          j1 = 13
        endif
        if (strn(9:9) .ne. ':' .or.
     .    .not. a2bin(strn,fptol,4,0,' ',j1,len(strn))) fptol = 1d-5
      else
        fptol = 0
      endif

C ... When generating the group the basis may become enlarged ... 
C     copy larger files relevant arrays to larger workspace
      call defi (oistab,(ngmx+1)*nbas)
C     call defi (oistb2,(ngmx+1)*nbas)
      call defi (oips2, ngmx*nbas)
      call defdr(opos2, 3*ngmx*nbas)
      call icopy(nbas,w(oips),1,w(oips2),1)
      call dcopy(3*nbas,w(opos),1,w(opos2),1)
      call gensym(slabl,gens,usegen,T,F,fptol,F,nbas,nspec,ngmx,plat,
     .  plat,ldist,dist,w(opos2),w(oips2),w(onrspc),nsgrp,w(osymgr),
     .  w(oag),ngen,gen,ssymgr,nggen,isym,w(oistab))
      if (nbas .gt. nbas0) call rxs('gensym: the basis was enlarged.',
     .  ' Check group operations.')
      if (nggen .gt. nsgrp) then
        call info2(1,0,0,' MKSYM (warning): generators create more than'
     .    //' ngmx=%i group ops ...',ngmx,0)
        ngmx = ngmx*16
        call rlse(oag)
        goto 5
      endif
      if (fptol .ne. 0) then
        call spackv(11,'site pos',ssite,1,nbas,w(opos2))
        call spackv(10,'site pos',ssite,1,nbas,w(opos))
      endif
      call rlse(oips2)

C --- Add inversion to point group ---
      incli = -1
      npgrp = nsgrp
      if (mod(mode,10) .ne. 0) then
        ngen = ngen+1
        call dpzero(gen(1,ngen),9)
        gen(1,ngen) = -1
        gen(5,ngen) = -1
        gen(9,ngen) = -1
        call pshpr(iprint()-40)
        call grpgen(gen(1,ngen),1,w(osymgr),npgrp,ngmx)
        call poppr
        incli = npgrp-nsgrp        
      endif
C     Not necessary, because it will be remade anyway
C     call redfi(oistab,npgrp*nbas)

C --- Printout of symmetry operations ---
      call info5(20,0,0,' MKSYM:  found %i space group operations'//
     .  '%?#n#; adding inversion generated %i ops#%j#'//
     .  '%?#(n==0)# ... includes inversion##',nsgrp,npgrp-nsgrp,npgrp,
     .  incli,0)
      if (incli .eq. -1) then
      call info0(30,0,0,'         no attempt to add inversion symmetry')
      endif
      if (mod(mode/10,10) .eq. 0) goto 100

C --- Split species into classes ---
      call defi(onrc, nspec)
      call icopy(nsite,w(oips),1,w(oipc),1)

      call splcls(mod(mode/10,10).eq.4,w(opos),nbas,nsgrp,w(oistab),
     .  nspec,slabl,nclass,w(oipc),w(oics),w(onrc))

      call upack('array npadl npadr',sarray,npadl,npadr,0,0,0)
      nclasp = nclass
      nclspp = nclass
      nbasp  = nbas + npadl + npadr
C ... (pgf) Extra classes for padding layers
      lpgf = igets('ctrl lpgf',sctrl)
      if (lpgf .gt. 0 .and. lpgf .lt. 10) then
C   ... Remake ipc,ics for padding and double-padding layers
        call pgfpp(nbas,npadl,npadr,1,1,.true.,.false.,w(oipc),w(oipc))
        call redfi(oics,nsite)
        call defi(oiwk,nsite)
        call pvsym1(nbas,npadl+npadr,nclasp,w(oipc),w(oics))
        nclspp = nclasp
        call pvsym1(nbasp,npadl+npadr,nclspp,w(oipc),w(oics))
        if (iprint() .ge. 30) then
          call awrit2('         %i classes for embedding region + %i '//
     .      'for padding layers',' ',80,stdo,nclass,nclasp-nclass)
        endif
        call redfi(oics,   nclspp)
C   ... Separate ipc, ipcd (dble padded ic eqiv =sngl pad vs not equiv)
        call defi (oipcp,  nsite)
        call icopy(nsite,w(oipc),1,w(oipcp),1)
        call icopy(npadl+npadr,w(oipc+nbas),1,w(oipc+nbasp),1)
      else
C   ... Reallocate arrays as permanent arrays
        call redfi(oics,   nclspp)
        oipcp = oipc
      endif

C ... Remake istab
      if (mod(mode,10) .eq. 2) then
        call defi(oistab,  npgrp*nbas)
      else
        call defi(oistab,  nsgrp*nbas)
      endif
      call dinv33(plat,1,qlat,xx)
C     call prmx('pos',w(opos),3,3,nbas)
C      call symtab(nbas,w,w(opos),w(osymgr),w(oag),-nsgrp,plat,qlat,w,
C     .  w(oistab))
C      call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
C      call symtbl(0,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
C     .  w(oistab))
C      call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
      call symtbl(1,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
     .  w(oistab))
C     Debugging check
C     call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
C     call symtbl(0,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
C     .   w(oistb2))
C     call shoist(w(oistb2),nbas,w(oag),w(osymgr),nsgrp)
C     call istbpm(w(oistab),nbas,nsgrp,w(oistb2))
C     call shoist(w(oistb2),nbas,w(oag),w(osymgr),nsgrp)
C     stop

C   ... For each i*grpop fix translation ag to be that of grpop
        if (mod(mode,10) .eq. 2) then
          if (aginv(w(osymgr),nsgrp,npgrp,nbas,w(oag),w(oistab)) .ne. 0)
     .      call rx('bug in mksym')
        endif

C ... poke ipc into ssite
      call spackv(11,'site class',ssite,1,nbasp,w(oipc))

C --- Create nrc and class labels ---
      if (mod(mode/10,10) .eq. 2 .or. mod(mode/10,10) .eq. 4) then
        call defrr(oclabl, nclspp)
        call defi (onrcp,  nclspp)
        call pvsym2(3,2*nbasp-nbas,nclspp,w(oics),w(oipcp),nspec,slabl,
     .    ssite,w(oclabl),w(onrcp))
        if (nclasp .eq. nclass) then
          onrc = onrcp
        else
          call defi (onrc,   nclasp)
          call icopy(nclasp,w(onrcp),1,w(onrc),1)
          call iaxpy(nclasp-nclass,1,w(onrcp+nclasp),1,w(onrc+nclass),1)
C          call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
C          call awrit2('%n:1i',' ',100,6,nclspp,w(onrcp))
        endif
      else
        oclabl = 0
        onrcp = 0
      endif

C       call awrit2('%n:1i',' ',100,6,nsite,w(oipc))
C       call awrit2('%n:1i',' ',100,6,nsite,w(oipcp))
C       pause

C --- Cleanup: poke class and symmetry info into structures ---
      call pack1('ctrl nclass',sctrl,nclass)
      call pack5('array onrcp oipcp',sarray,onrcp,oipcp,0,0,0)
      call pack5('array oics oipc oclabl nclasp onrc',sarray,oics,oipc,
     .  oclabl,nclasp,onrc)
c     call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
  100 continue
      call pack5('lat npgrp nsgrp oistab oag osymgr',slat,npgrp,nsgrp,
     .  oistab,oag,osymgr)

      end
      subroutine pvsym1(nbas,npad,nclass,ipc,ics)
C- Fix classes for padding layers
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas   number of basis atoms without padding layers
Ci   nclass number of classes without padding layers
Ci Outputs
Co   nclass updated
Co   ipc(nbas..nbas+npad) updated
Co   ics updated
C ----------------------------------------------------------------------
C     implicit none
      integer nbas,nclass,ipc(nbas),ics(nclass)
      integer i,npad,nc,is,ic,j
      if (npad .eq. 0) return

C ... ipc->-ipc to flag which classes have been renumbered
      do  5  i = nbas+1, nbas+npad
    5 ipc(i) = -ipc(i)

C ... For each site in padding layer, assign class
      nc = nclass
C ... Look for next occurrence of class not assigned
      do  10  i = nbas+1, nbas+npad

        if (ipc(i) .gt. 0) goto 10
C   ... Class ipc(i) need to be assigned:  do for each site in class
C       ic is (-)original class, is is species of original class
        ic = ipc(i)
        is = ics(-ic)

        nc = nc+1
        do  12  j = i, nbas+npad
   12   if (ipc(j) .eq. ic) ipc(j) = nc

C   ... Species index is same as that of original class
        ics(nc) = is

   10 continue

      nclass = nc
C      call awrit2('%n:1i',' ',100,6,nbas+npad,ipc)
C      call awrit2('%n:1i',' ',100,6,nclass,ics)
C      pause

      end
      subroutine pvsym2(mode,nbas,nclass,ics,ipc,nspec,slabl,ssite,
     .  dclabl,nrc)
C- Create class labels from species labels (double precision format)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nosplt: T copy class and species
Ci     mode: 0 do nothing
Ci           1 create class labels clabl
Ci           2 create number of sites in each class nrc
Ci      ipc: for padding sites ib ipc(ib) = class
Ci   nclass: number of classes
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    clabel
Co Outputs:
Co   dclabl: class labels in double precision format
Co      nrc: number of sites in each class
Cu Updates
Cu   18 Dec 01 Packs class label into ssite->clabel
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nbas,nclass,nspec,ics(1),ipc(nbas),nrc(1)
      double precision dclabl(nclass),ssite(1)
      character*8 slabl(nspec)
      integer ic,iclbsj,idx,is,ib
      character*8 clabl

C --- Make class labels from species labels ---
      if (mod(mode,2) .eq. 1) then
      do  10  is = 1, nspec
        do  12  idx = 1, nbas
          ic = iclbsj(is,ics,-nclass,idx)
          if (ic .lt. 0) goto 13
          call clabel(slabl,is,idx,clabl)
          call s8tor8(clabl,dclabl(ic))
   12   continue
   13   continue
   10 continue
      endif

      do  20  ib = 1, nbas
        ic = ipc(ib)
        call r8tos8(dclabl(ic),clabl)
        call spacks(1,'site clabel',ssite,clabl,ib,ib)
   20 continue

C      do  ib = 1, nbas
C        call spacks(0,'site clabel',ssite,clabl,ib,ib)
C        print *, ib, clabl
C      enddo
C      pause

C --- Create nrc ---
      if (mod(mode/2,2) .eq. 1) then
      call iinit(nrc,nclass)
      do  30  ib = 1, nbas
        ic = ipc(ib)
        nrc(ic) = nrc(ic)+1
   30 continue
      endif

C ... debugging : check class labels
C      call awrit2('%n:1i',' ',100,6,nclass,nrc)
C      do  50  ic = 1, nclass
C        call r8tos8(dclabl(ic),clabl)
C        print 333, ic, clabl
C  333   format(i4,1x,a)
C   50 continue
C      pause

      end

      integer function aginv(g,npg,nsg,nbas,ag,istab)
C- Finds translation part for extra inversion-generated group ops
C ----------------------------------------------------------------------
Ci Inputs
Ci   g     :point group operations
Ci   npg   :number of point group operations
Ci   nsg   :number of space group operations
ci   nbas  :dimensions istab.  If zero, istab is not altered.
Co Outputs
Co   ag    :for each group op ig in range nsg+1..npg,
Co         :the corresponding group op jg is found such that
Co         :g(jg) = -g(ig), and ag(jg) is copied to ag(ig); see Remarks
Co   istab :istab(jg) is copied into istab(ig) just in the same way as
Co         :ag(jg) is copied into ag(ig).
Co         :The sign of istab(1) is set negative to indicate that this
Co         :operation is derived from a true space group operation.
Cr Remarks
Cr   In certain cases the inversion operation may be added to the space
Cr   group, for purposes of k integration.  This is permissible when the
Cr   hamiltonian has the form h(-k) = h*(k).  In that case, the
Cr   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
Cr
Cr   Also, the Green's functions are related G(-k) = Gtranspose(k).
Cr   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
Cr   then G(-k) = Gtranspose(k), and the same (g,ag) information is
Cr   needed for either rotation.  Let g~ = -g:
Cr     G0(g~^-1 -k) into G(k),
Cr
Cr
Cr   Routine aborts with aginv=ig if for group ig it cannot find
Cr   corresponding jg such that g(jg) = -g(ig)
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npg,nsg,nbas,istab(nbas,1)
      double precision g(3,3,npg),ag(3,npg)
C ... Local parameters
      integer ig,jg
      double precision ginv(3,3)
      logical grpeql

      do  ig = npg+1, nsg
        call dpcopy(g(1,1,ig),ginv,1,9,-1d0)
        do  jg = 1, nsg
          if (grpeql(g(1,1,jg),ginv)) then
            call dcopy(3,ag(1,jg),1,ag(1,ig),1)
            if (nbas .gt. 0) then
              call icopy(nbas,istab(1,jg),1,istab(1,ig),1)
              istab(1,ig) = -istab(1,jg)
            endif
C           call awrit2(' aginv mapping %,2i ->%,2i',' ',80,6,ig,jg)
            goto 10
          endif
        enddo
        aginv = ig
        return
   10   continue
      enddo
      aginv = 0
      end

