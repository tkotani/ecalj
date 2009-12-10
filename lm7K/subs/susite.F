      subroutine susite(sctrl,sham,spot,sarray,slat,sspec,ossite)
C- Sets up permanent arrays related to basis
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nsite nbas nl lpgf
Ci     Stored:    npl sdxsi nsite nbasp
Ci     Passed to: lgors
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: *
Ci     Stored:    neula oeula
Ci     Passed to: *
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: *
Ci     Stored:    ovshf oaamom
Ci     Passed to: *
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: *
Ci     Stored:    ncl ocllst oclp oclssl opgfsl opgfvl opgplp opos oips
Ci                npadl npadr
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat platl platr
Ci     Stored:    opos
Ci     Passed to: clsset
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: eref name
Ci     Stored:    *
Ci     Passed to: *
Ci   w(ossite): array containing site data
Ci     Elts read: *
Ci     Stored:    clabel pos spec pl plv eula
Ci     Passed to: dvprm dvperm pgfpp pvsub2
Co Outputs
Cr Remarks
Cr  *Arrays which susite sets:
Cr      w(opos),w(oips),w(opgplp),w(oeula)
Cr      w(opgfsl),w(opgfvl),w(opgord) for PGF
Cr
Cr  *pgplp holds information about crystal subblocks, and its meaning
Cr   depends on the context.  pgplp is dimensioned pgplp(6,-1:*)
Cr   pgplp(1,-1) contains information about the context-dependence,
Cr   to permit certain routines that depend on the context to be
Cr   applied across different ones.
Cr
Cr    *For the usual crystal case, there are no subblocks.
Cr     In this case, pgplp(1,-1) = -1.  susite creates and sets pgplp.
Cr
Cr    *For the planar Green's function, subblocks are principal layers.
Cr     In this case, it two additional subblocks are added to treat
Cr     the boundaries: one to the "left" of the first PL, and one to the
Cr     "right" of the last.  Here, pgplp(1,-1) = 0.  Some rows
Cr     of pgplp are initially set by pgfset; others are set in pgfasa.
Cr     Information is kept about npl+2 layers.
Cr
Cr    *For the order-N Green's function, subblocks are the crystal
Cr     subblocks.  Here, pgplp(1,-1) = -2.
Cr
Cu Updates
Cu   24 May 08 eula and magf read from master only, MPI mode
Cu   07 Feb 03 looks for, and reads magnetic field from file
Cu   23 Jan 02 susite packs seref into sham.  Altered argument list.
Cu   18 Dec 01 susite packs species labels into ssite
Cu   13 Oct 98 convention for pgplp extended to non-layer case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ossite
      double precision sctrl(1),sham(1),spot(1),sarray(1),slat(1),
     .  sspec(1)
C ... Local parameters
      integer procid,master,mpipid
      logical lgors,pass1,lorder,mlog
      integer lpgf,igets,nsite,nbas,nl,nbasp,npadl,npadr,ncl,npl,i,ib,
     .  ix(3),ssize,neul,nbf,opgplp,opgord,opos,oips,ovshft,opgfsl,
     .  opgfvl,owk,oeula,omagf,oaamom,oclssl,ocllst,oclord,oclp,nclp,
     .  igetss,iprint
      parameter (nclp=9)
      double precision plat(9),platl(9),platr(9),pgfn(3),xsi(5)
      double precision dval,dglob,seref,dgetss
      character*8 spid
C ... Heap
      integer w(1)
      common /w/ w

      call upack('ctrl nsite nbas nl',sctrl,nsite,nbas,nl,0,0)

C     call rxx(nsite.ne.nbas,'susite not set up for nsite.ne.nbas')
      ssize = nint(dval(w(ossite),1))
      pass1 = .false.
      nbasp = nbas
      npl   = nbas
      master = 0
      procid = mpipid(1)
      mlog = .false.

C ... Pack species labels; make sumeref
      seref = 0
      do  ib = 1, nsite
        i = igetss('site spec',ib,w(ossite))
C       Skip over sites with no species (e.g. multipoles)
        if (i .gt. 0) then
          call spacks(0,'spec name',sspec,spid,i,i)
          call spacks(1,'site clabel',w(ossite),spid,ib,ib)
          seref = seref + dgetss('spec eref',i,sspec)
        endif
      enddo
      call pack1('ham seref',sham,seref)

C ... Start again if ssite was reallocated
   10 pass1 = .not. pass1
      call upack('lat plat platl platr',slat,plat,platl,platr,0,0)
      call defrr(opos,3*nsite)
      call spackv(10,'site pos',w(ossite),1,nsite,w(opos))
      call defi(oips,nsite)
      call spackv(10,'site spec',w(ossite),1,nsite,w(oips))

C --- Order-N: group all sites in a crystal into clusters ---
      if (lgors('ctrl lordn,-1',sctrl)) then
        call defi(oclssl, nsite)
        call defi(ocllst, nsite+1)
        call defi(oclp,  -(nsite+1)*nclp)
        call defi(oclord, nsite)
        call clsset(11,slat,nsite,w(opos),ncl,w(oclssl),w(ocllst),
     .    w(oclord),w(oclp))

C   ... Reorder the site structure according to clord
        call defdr(owk, ssize*nsite)
        call dvprm(ssize,nsite,w(ossite),w(owk),w(oclord),.true.)
        call rlse(owk)
        call dvset(w(ossite),1,1,dble(ssize))
C   ... Unpack arrays, now permuted; shorten pos
        ix(1) = 1
        ix(2) = 1
        ix(3) = 1
        call spackv(10,'site pos',w(ossite),1,nsite,w(opos))
        call shorps(nbas,plat,ix,w(opos),w(opos))
        call spackv(11,'site pos',w(ossite),1,nsite,w(opos))
        call spackv(10,'site spec',w(ossite),1,nsite,w(oips))
C       call shstru('site',w(ossite),1,nsite)
        call redfi(oclp,   nclp*(ncl+1))
        call pack5('array ncl ocllst oclp oclssl',sarray,ncl,ocllst,
     .    oclp,oclssl,0)
      endif

C --- Check and order principal layers, and sites by PL ---
      lpgf = igets('ctrl lpgf',sctrl)
      if (lpgf .ne. 0) then
        call defi(opgfsl,-nsite)
        call spackv(10,'site pl',w(ossite),1,nsite,w(opgfsl))
        call defi(opgfvl,-nsite)
        call spackv(10,'site plv',w(ossite),1,nsite,w(opgfvl))
        call defdr(ovshft, -8-2*nbasp)
        call defi(opgplp,  6*(npl+2))
C        print *, 'debuging susite'
C        call ivset(w(opgplp),1,6*(npl+2),-99)
C   ... no accumulated sites in -1 PL
        w(opgplp) = 0
        call defi(opgord,  2*nbas)
        if (pass1) call pshpr(min(iprint(),1))
C       if (pass1) call pshpr(0)
C       call prmx('pos before pgfset',w(opos),3,3,nbas)
        lorder = .not. pass1
        call pgfset(sspec,nbas,w(opos),plat,lorder,.true.,w(oips),
     .    w(ovshft),w(opgfsl),w(opgfvl),w(opgord),pgfn,npl,npadl,npadr,
     .    w(opgplp))
        nbasp = nbas + (npadl + npadr)
        ix(1) = dglob('nbasp',dble(nbasp),1)
        ix(2) = dglob('npl',dble(npl),1)
        nsite = nbas + 2*(npadl + npadr)
        if (pass1) then
          call poppr
          call redfrr(ossite,nsite*ssize)
C         Zero out enlarged last part of it
          call defrr(owk,-nsite*ssize)
          call dcopy(nbas*ssize,w(ossite),1,w(owk),1)
          call dcopy(nsite*ssize,w(owk),1,w(ossite),1)
          call rlse(owk)
          goto 10
        endif
        call defdr(owk, ssize*nbas)
C   ... Reset pos,pfgsl, since they may be shifted by pgfset
        call spackv(11,'site pos',w(ossite),1,nsite,w(opos))
        call spackv(11,'site pl',w(ossite),1,nsite,w(opgfsl))
C   ... Reorder the site structure according to pgord
        do  12  i = 1, nbas
   12   w(opgord+i-1) = w(opgord+i-1)-1
C       call awrit2('%n:1i',' ',100,6,nbas,w(opgord))
C       call prmx('ssite before perm',w(ossite),ssize,ssize,nbas)
        call dvperm(ssize,nbas,w(ossite),w(owk),w(opgord),.true.)
        call dvset(w(ossite),1,1,dble(ssize))
C       call prmx('ssite after perm',w(ossite),ssize,ssize,nbas)
C   ... Doubly pad ssite
        call pgfpp(nbas,npadl,npadr,ssize,1,.false.,.false.,
     .    w(ossite),w(ossite))
C   ... Unpack arrays, now permuted and padded
        call spackv(10,'site pl',w(ossite),1,nsite,w(opgfsl))
        call spackv(10,'site plv',w(ossite),1,nsite,w(opgfvl))
        call spackv(10,'site pos',w(ossite),1,nsite,w(opos))
        call spackv(10,'site spec',w(ossite),1,nsite,w(oips))

C   ... Shift doubly padded bas; repack
        if (nbasp .gt. nbas) then
        call pgbasp(nbas,npadl,npadr,w(opos),plat,platl,platr,w(owk))
        call dscal(9,2d0,platl,1)
        call pgbasp(nbasp,npadl,npadr,w(owk),plat,platl,platr,w(opos))
        call dscal(9,.5d0,platl,1)
        call spackv(11,'site pos',w(ossite),1,nsite,w(opos))
C       call prmx('pos, double pad',w(opos),3,3,nsite)
C       call prmx('ssite, double pad',w(ossite),ssize,ssize,nsite)
C   ... Shift pgfsl in padding layers, and repack
        call pvsub1(w(opgfsl),nbas,npadl,npadr)
        call spackv(11,'site pl',w(ossite),1,nsite,w(opgfsl))
C        call awrit2('%n:1i',' ',100,6,nsite,w(opgfsl))
C        call awrit2('%n:1i',' ',100,6,nbas,w(opgfvl))
C   ... Pick up padded ips,pgfsl,pgfvl
        call spackv(10,'site spec',w(ossite),1,nsite,w(oips))
        call spackv(10,'site plv',w(ossite),1,nsite,w(opgfvl))
        call redfi(opgplp,   6*(npl+2))
        call pack5('array opgfsl opgfvl opgplp',sarray,opgfsl,opgfvl,
     .    opgplp,0,0)
        call pack1('ctrl npl',sctrl,npl)
        call pack1('pot ovshf',spot,ovshft)
C       call shstru('site',w(ossite),1,nsite)
        endif
        call rlse(opgord)
        call pack5('ctrl nsite nbasp',sctrl,nsite,nbasp,0,0,0)
      else
C   ... No partitioning of hamiltonian: make global pgplp
        npadl = 0
        npadr = 0
        call defi(opgplp, -12)
        w(opgplp+0) = -1
        w(opgplp+6) =  nbas
        call defdr(ovshft, -8-2*nbas)
        call pack5('array opgfsl opgfvl opgplp',sarray,1,1,opgplp,0,0)
        call pack1('ctrl npl',sctrl,1)
        call pack1('pot ovshf',spot,ovshft)
      endif

C --- Euler angles and external magnetic fields ---
      if (lgors('ctrl lncol,11',sctrl)) then
        call defdr(oeula, -nbasp*nl**2*3+3)
        call dpzero(xsi,5)
        if (procid .eq. master) then
          call pvsub2(0,w(ossite),nbas,nbasp,nl,w(oeula),neul,xsi)
        endif
        call mpibc1(neul,1,2,mlog,'susite','neul')
        call mpibc1(xsi,5,4,mlog,'susite','xsi')
        call redfrr(oeula, nbasp*neul*3+3)
        call mpibc1(w(oeula),nbas*neul*3,4,mlog,'susite','eula')
        call pack1('ctrl sdxsi',sctrl,xsi)
        call pack2('ham neula oeula',sham,neul,oeula)

        nbf = 1
        omagf = 1
        if (lgors('ctrl lncol,8',sctrl)) then
          call defdr(omagf, -nbasp*nl**2*3)
          if (procid .eq. master) then
            call pvsub2(1,w(ossite),nbas,nbasp,nl,w(omagf),nbf,xsi)
          endif
          call mpibc1(nbf,1,2,mlog,'susite','nbf')
          call mpibc1(w(omagf),nbas*neul*3,4,mlog,'susite','magf')
          call redfrr(omagf, nbasp*nbf*3)
C         call prmx('b field',w(omagf),nbasp*nbf,nbasp*nbf,3)
        endif
        call pack2('ham nbf omagf',sham,nbf,omagf)

        call defdr(oaamom, -nbasp)
        call pack1('pot oaamom',spot,oaamom)
      else
        call defdr(oeula,1)
        call pack2('ham neula oeula',sham,0,oeula)
      endif

      call pack5('lat opos',slat,opos,0,0,0,0)
      call pack5('array opos oips npadl npadr',sarray,opos,oips,npadl,
     .  npadr,0)
C     call pack5('ctrl nsite nbasp',sctrl,nsite,nbasp,0,0,0)
C     call shstru('site',w(ossite),1,nsite)

      end
      subroutine pvsub1(pgfsl,nbas,npadl,npadr)
C- Sets PL indices for padding layers
C     implicit none
      integer nbas,npadl,npadr,pgfsl(nbas)
      integer nbasp,i

      nbasp = nbas + (npadl + npadr)
      do  10  i = 1, npadl
      pgfsl(nbas+i) = pgfsl(nbas+i) - 1
   10 pgfsl(nbasp+i) = pgfsl(nbasp+i) - 2
      do  20  i = 1, npadr
      pgfsl(nbas+npadl+i) = pgfsl(nbas+npadl+i) + 1
   20 pgfsl(nbasp+npadl+i) = pgfsl(nbasp+npadl+i) + 2

      end
      subroutine pvsub2(mode,ssite,nbas,nbasp,nl,eula,neul,xsi)
C- Unpacks Euler angles or magnetic fields; reads from disk
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 for euler angles, 1 for bfield
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read:
Ci     Stored:    eula
Ci     Passed to: spackv
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   nl    :(global maximum l) + 1
Co Outputs
Co   eula  :(mode 0) Euler angles for noncollinear spins
Co         :(mode 1) Magnetic field read
Co   neul  :1 if Euler angles (Bfield) are l-independent, nl otherwise
Co   xsi   :global deamon parameters for spin dynamics (not read mode 1)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   29 Jan 03 Added mode to also read b-field
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,nbasp,nl,neul
      double precision ssite(1),eula(nbasp,*)
C ... Local parameters
      integer oeula,ifi,fopn,fxst
      double precision xsi(5)
      character name*6
C ... Heap
      integer w(1)
      common /w/ w

      call isanrg(mode,0,1,'pvsub2','mode',.true.)
      if (mode .eq. 0) name = 'eula  '
      if (mode .eq. 1) name = 'bfield'

      neul = 1
      call defr(oeula,3*nbasp)
      call spackv(10,'site '//name,ssite,1,nbasp,w(oeula))
C     call prmx(name,w(oeula),3,3,nbasp)
      call dmcpy(w(oeula),1,3,eula,nbasp,1,nbasp,3)
C     call prmx(name//' transposed',eula,nbasp,nbasp,3)

C ... Read angles from disk, if available
      if (fxst(name) .eq. 1) then
        ifi = fopn(name)
        rewind ifi
        if (mode .eq. 0) call ioeula(nbasp,nl,eula,neul,xsi,ifi)
        if (mode .eq. 1) call iomagf(nbasp,nl,eula,neul,ifi)
        call fclose(ifi)
      endif
C     call prmx(name//' yet again',eula,nbasp*neul,nbasp*neul,3)

C     At present, no noncollinearity in padding layers
      if (nbasp .gt. nbas)
     .  call dmcpy(0d0,0,0,eula(1+nbas,1),nbasp,1,nbasp-nbas,neul*3)

C     Repack if neula is 1
      if (neul .eq. 1) then
        call dmcpy(eula,nbasp,1,w(oeula),1,3,nbasp,3)
        call spackv(11,'site '//name,ssite,1,nbasp,w(oeula))
      endif

      call rlse(oeula)

      end

