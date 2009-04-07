      subroutine mmag(sctrl,sarray,sspec,slat,smove,sstrn)
C- Micromagnetics simulations
C ----------------------------------------------------------------------
Ci   sctrl :struct containing parameters governing program flow
Ci   sarray:structure containing offsets to various arrays
Ci   sspec :struct containing species-specific information
Ci   slat  :struct containing information about the lattice
Ci   sstrn :struct containing global strings
Co Outputs
Cr Remarks
Cr   Code was originally written for an earlier version of LMTO.
Cr   It has been adapted to lm-6.9 up to the point of generating
Cr   the hamiltonian.  The rest is commented out, as the
Cr   structure information needs be recast into smove.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mxcsiz
      double precision sctrl(1),sspec(1),slat(1),sarray(1),smove(1)
      character*(*) sstrn
C Local
      character htype(2)*5
      double precision alat,plat(3,3),avw
      integer i,i1,i2,igets,j,lpbc,lpgf,nbas,nbasp,nspec,nttab,ltype(2),
     .  ntype(2),oepair,oesite,oiax,oipair,oips,oisite,ontab,opgfsl,
     .  opos,orham,ortab
C ... Heap
      integer w(1)
      common /w/ w
      data htype /'h1:','h2:'/

      call upack('ctrl nbas nbasp nspec',sctrl,nbas,nbasp,nspec,0,0)
      call upack('lat alat plat avw',slat,alat,plat,avw,0,0)
      call upack('array opgfsl oips opos',sarray,opgfsl,oips,opos,0,0)
      call defrr(orham,nspec)
      call spackv(10,'spec rham',sspec,1,nspec,w(orham))
      lpgf = igets('ctrl lpgf',sctrl)
      lpbc = 0
      if (lpgf .gt. 0) lpbc = 11

C --- Make iax and rtab ---
      if (lpbc .eq. 0) then
        i = 3
        j = -1
      elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
        i = 2
        j = 1
      else
        call rx('MMHAM: not implemented for lpbc>1')
      endif
      mxcsiz = 0
      call pairs(nbas,nbasp,alat,plat,w(orham),w(opos),
     .  w(oips),i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
      call defrr(ortab,3*nttab)
      call mkrtab(000,alat,plat,w(opos),w(oiax),nttab,w(opos),w(ortab))

C --- Make the Hamiltonian ---
      call defrr(oesite, nbasp)
      call defi(oisite,  nbasp)
      call defrr(oepair, nttab)
      call defi(oipair,  nttab)
      call upacks('strn mmham',i1,i2)
      ntype(1) = 1
      ntype(2) = 1
      ltype(1) = -1
      ltype(2) = -1
      call suemph(sspec,sstrn(i1:i2),w(ontab),w(oiax),w(ortab),nbasp,
     .  nttab,w(oips),ntype,ltype,htype,w(oesite),w(oisite),w(oepair),
     .  w(oipair))

C      call defrr(owk,-(3*nbas+3))
C      call u_bsi(s_sdyn,intgmd,nmodt,modt,ct,xx,xx,xx,mx,-1)
C      neq = nbas*2 + nmodt
C      if (intgmd .eq. 1) call defrr(obswk,(neq*(4+mx)+5))
C      nstep = -1
C      call defrr(ofrc, -(nbas*3+3))
CC  ... For historical compatibility use -1, but should be 1
C       call dvset(w(oaamom),1,nbas,1d0)
C
CC --- Micromagnetic dynamics ---
C      ifi = fopn('save')
C      call poseof(ifi)
C   62 call mm_dyn(nbas,xsi,w(oeula),
C     .  qspirl,w(owk),w(oebil),w(oibil),nstep,irmmd,w(oiax),nttab,
C     .  w(obswk),w(ofrc),amag,w(oaamom),etot,s_sdyn,nvarms,ifi)
C
CC --- Calculate the pairwise forces and continue ---
C      if (nstep .gt. 0) then
C        etot = 0
C        call dpzero(w(ofrc),nbas*3+3)
C        call rx('mmag: need update mmpair')
C        call mmpair(nbas,nttab,w(oiax),w(oebil),w(oibil),w(owk),
C     .   etot,amag,w(ofrc))
C        goto 62
C      endif

      end

