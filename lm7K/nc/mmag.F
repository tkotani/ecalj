      subroutine mmag(sctrl,sham,sarray,sspec,slat,smove,sstrn)
C- Micromagnetics simulations
C ----------------------------------------------------------------------
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nbasp nspec sdxsi nvario lpgf
Ci     Stored:
Ci     Passed to:
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: qss oeula
Ci     Stored:
Ci     Passed to:
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: opgfsl oips opos
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:
Ci     Stored:    rham
Ci     Passed to: spackv suemph
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat avw
Ci     Stored:
Ci     Passed to:
Ci   smove :struct for dynamics information; see routine umove
Ci     Elts read: nmodt modt ct prmint:1.15.5
Ci     Duplicate: prmint:
Ci     Stored:
Ci     Passed to: mmdyn
Ci   sstrn :struct for global strings
Ci     Elts read: mmham
Ci     Stored:
Ci     Passed to:
Cl   nmodt :number of thermostat modes
Cl   modt  :modt(1..nmodt) : thermostat modes
Cl   ct    :ct(1..nmodt)   : thermostat mode coefficients
Cl   neq   :total number of integration variables
Cl   intgmd:integration mode
Cl         :1 Bulirsch-Stoer integration
Cl   amag  :average magnetization
Co Outputs
Cr Remarks
Cr   mmag integrates the Landau-Lifshitz equations of motion
Cr   with global daemons.  It operates in one of two modes:
Cr
Cr 1.  An empirical (heisenberg) hamiltonian H is supplied.
Cr
Cr     It is supplied EITHER from the contents of a file, OR in the
Cr     form of a set of empirical rules.  In the former case,
Cr     structure element 'strn mmham' contains a string
Cr     FILE=filename; the contents of H are supplied by file
Cr     filename.  The syntax of that file are desribed in (and data
Cr     are read by) subroutine iorsj.  Alternatively, structure
Cr     element 'strn mmham' contains a set of rules that define the
Cr     Heisenberg coefficients.  In this latter case, H is assembled
Cr     by suemph.f, which also documents the syntax for the rules.
Cr
Cr     In this mode mmag calls mmdyn, which computes the forces
Cr     internally and integrates the equations of motion over the
Cr     specified interval (see struc smove).
Cr
Cr 2.  The magnetic forces are supplied by the calling program.
Cr
Cr     In this mode mmag calls mmdyn for a single time, supplying
Cr     the Eular angles and forces for that time.  mmdyn returns
Cr     with a new time for which it requires updated Eular angles
Cr     and forces.
Cr
Cu Updates
Cu   02 Oct 03 Implement --wcvmj to write CVM file
Cu   03 Jun 03 Writes to simpler sv file, rather than save file
Cu   24 Nov 02 Micromagnetics hamiltonian may be read from a file
Cu   17 Aug 01 updated from ASA package version 5.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mxcsiz
      double precision sctrl(1),sham(1),sspec(1),slat(1),sarray(1),
     .  smove(1)
      character sstrn*(*)
C Local
      character outs*256,jfilnm*72
      logical cmdopt,a2bin,parstr,rdstrn
      character htype(2)*5,dc*1
      double precision alat,plat(3,3),avw,xsi(3),qspirl(4),ct(3),
     .  amag(3),etot,tnow,xx,ttot,sclrsj,tolrsj
      integer oaamom,oepair,oesite,oeula,ofrc,ohpair,oiax,ointwk,oipair,
     .  oips,oisite,ontab,opgfsl,opos,orham,ortab,owk
      integer fopna,i,i1,i2,ifi,jfi,igets,j,j1,j2,k,lpbc,lpgf,ltype(2),
     .  modt(3),nbas,nbasp,neq,nmodt,nspec,nstep,nttab,ntype(2),nvarms,
     .  iorsj,iprint,parg,nfilm,nfilj,istyle
C     Symmetry operations
      integer nsgrp,og,oag,oistab,oistb2,ohuns
C     Specific to integration mode
      integer intgmd,mx
C ... Heap
      integer w(1)
      common /w/ w
      data htype /'h1:','h2:'/

      call upack('ctrl nbas nbasp nspec sdxsi nvario',sctrl,nbas,nbasp,
     .  nspec,xsi,nvarms)
      call upack('lat alat plat avw',slat,alat,plat,avw,0,0)
      call upack('array opgfsl oips opos',sarray,opgfsl,oips,opos,0,0)
      call defrr(orham,nspec)
      call spackv(10,'spec rham',sspec,1,nspec,w(orham))
      call upack('ham qss oeula',sham,qspirl,oeula,0,0,0)
      lpgf = igets('ctrl lpgf',sctrl)
      lpbc = 0
      if (lpgf .gt. 0) lpbc = 11

C --- Make or read the Hamiltonian ---
      ntype(1) = 1
      ntype(2) = 1
      ltype(1) = 0
      ltype(2) = -1
      call upacks('strn mmham',i1,i2)

      call defrr(oaamom,nbasp)
      call dvset(w(oaamom),1,nbasp,1d0)
      call defrr(oesite, nbasp)
      call defi(oisite,  nbasp)

C --- Read R.S. hamiltonian from file ---
      i = 0
      if (parstr(sstrn(i1:i2),'FILE=',i2-i1+1-5,5,'=',i,j)) then
        call word(sstrn(i1+j:),1,j1,j2)
        jfilnm = sstrn(i1+j+j1-1:i1+j+j2-1)
        ifi = fopna(jfilnm,-1,0)
        call info0(20,1,0,' ... Read r.s. hamiltonian from file '//
     .    jfilnm(1:j2-j1+1))

C   ... Setup call for file read
        nfilj = iorsj(0,ifi,nbasp,alat,plat,w(oaamom),w,w,w,w,xx)
        call rxx(nfilj.lt.0,'MMAG: missing or incompatible file'//
     .    ' lattice vectors')

C   ... Make a neighbor table that encompasses all file entries
        if (lpbc .eq. 0) then
          i = 3
          j = -1
        elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
          i = 2
          j = 1
        else
          call rx('MMHAM: not implemented for lpbc>1')
        endif
        call info2(35,0,0,' ... Make neighbor table: rcut = %d * alat',
     .    xx,0)
        call pshpr(iprint()-10)
        mxcsiz = 0
C       call tcn('pairs')
        call pairs(nbas,nbasp,alat,plat,alat*xx/2*(1+1d-6),w(opos),
     .    -1,i,j,w,nttab,ontab,oiax,mxcsiz)
C       call tcx('pairs')
        call defrr(ortab,3*nttab)
        call mkrtab(0,alat,plat,w(opos),w(oiax),nttab,w(opos),w(ortab))
        call poppr

        call defrr(oepair, -nttab)
        call defi(oipair,  -nttab)

C   ... Read Heisenberg coefficients from file
C       NB: note file's rtab does not scale by alat
        call dscal(3*nttab,1/alat,w(ortab),1)
        call defrr(ohpair,-nttab)
        nfilm = iorsj(3,ifi,nbasp,alat,plat,w(oaamom),w(oiax),w(ortab),
     .    w(ontab),w(ohpair),xx)
        if (nfilm .le. 0) call rx('mmag: no hamiltonian read')

        if (nfilj.ne.nfilm) call info2(20,0,0,'%8f(warning) only'//
     .    ' matched %i out of %i file entries',nfilm,nfilj)

C   ... Purge iax table of unused entries
        if (nfilm.ne.nttab) then
          call info0(20,0,0,' ... Reduce iax table and'//
     .      ' re-read pair hamiltonian')

C         Purge iax table of unused entries
          call symiax(1,plat,nbas,w,w,w,0,w(ontab),w(oiax),nttab,i)
          if (nttab .ne. nfilm) then
            call fexit2(-1,111,' Exit -1 : mmag: bug in symiax:'//
     .        ' nttab=%i, nfilm=%i',nttab,nfilm)
          endif
C         Remake rtab
          call mkrtab(000,1d0,plat,w(opos),w(oiax),nttab,w(opos),
     .      w(ortab))

C         Re-read pair hamiltonian with new iax table
          nfilm = iorsj(2,ifi,nbasp,alat,plat,w(oaamom),w(oiax),
     .      w(ortab),w(ontab),w(ohpair),xx)
          if (nttab .ne. nfilm) call rx('bug in iorsj')

C         Close input file rsj
          call fclose(ifi)

          call upack('lat nsgrp oistab osymgr oag',slat,nsgrp,oistab,og,
     .      oag,0)
          call defi(oistb2,nsgrp*nbas)
          call istbpm(w(oistab),nbas,nsgrp,w(oistb2))
          call defrr(ohuns,-nttab)
          call dcopy(nttab,w(ohpair),1,w(ohuns),1)
          i = 0 + 10*2*0 + 1000
C          call ppair1(40,1,nbas,nbas,w(oips),alat,plat,w(opos),xx,
C     .      nttab,w(ontab),w(oiax),mxcsiz)
          call rsmsym(i,plat,1,w,0,nbas,w(opos),1,1,
     .      1,nttab,w(ontab),w(oiax),w(og),w(oistb2),nsgrp,1,
     .      w(ohuns),w(ohpair))

        endif

        call dscal(3*nttab,alat,w(ortab),1)
C       Distribute into epair and esite
        call pmmag1(nbas,w(ontab),w(ortab),
     .    w(ohpair),w(oesite),w(oepair))
        call rlse(ohpair)

C --- Use empirical rules to make hamiltonian ---
      else

C   ... Make iax and rtab
        if (lpbc .eq. 0) then
          i = 3
          j = -1
        elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
          i = 2
          j = 1
        else
          call rx('MMHAM: not implemented for lpbc>1')
        endif

        call pshpr(iprint()-10)
        mxcsiz = 0
        call pairs(nbas,nbasp,alat,plat,w(orham),w(opos),
     .    w(oips),i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
        call defrr(ortab,3*nttab)
        call mkrtab(0,alat,plat,w(opos),w(oiax),nttab,w(opos),w(ortab))
        call poppr

        call defrr(oepair, -nttab)
        call defi(oipair,  -nttab)

        call suemph(sspec,sstrn(i1:i2),w(ontab),w(oiax),w(ortab),nbasp,
     .  nttab,w(oips),ntype,ltype,htype,w(oesite),w(oisite),w(oepair),
     .  w(oipair))
      endif

C --- Write R.S. Jij to file.  Switch syntax: --wrsj[:fn=name][:scl=#] ---
      if (cmdopt('--wrsj',6,0,outs)) then
        i = 6
        sclrsj = 1
        tolrsj = 1d-6
        jfilnm = 'rsj2'
        dc = outs(i+1:i+1)
        if (dc .ne. ' ') then
C     ... Return here to resume parsing for arguments
          j2 = i
   50     continue
          j2 = j2+1
          if (outs(j2:j2) .eq. dc) goto 50
          j1 = min(len(outs),j2)
          call nwordg(outs,0,dc//' ',1,j1,j2)
          if (j2 .ge. j1) then
            if (.false.) then
            elseif (outs(j1:j1+3) .eq. 'scl=')  then
              j = 0
              i = parg('scl=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,sclrsj)
              if (i .le. 0) goto 59
            elseif (outs(j1:j1+3) .eq. 'tol=')  then
              j = 0
              i = parg('tol=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,tolrsj)
              if (i .le. 0) goto 59
            elseif (outs(j1:j1+2) .eq. 'fn=')  then
              jfilnm = outs(j1+3:j2)
            else
              goto 59
            endif
            goto 50
   59       call rxs('mmag: failed to parse --wrsj options:',outs)
          endif
        endif

        call info0(20,0,0,' ... writing r.s. hamiltonian to file '//
     .    jfilnm)
        ifi = fopna(jfilnm,-1,0)
        call dscal(3*nttab,1/alat,w(ortab),1)
        i = iorsj(1,-ifi,nbasp,alat,plat,w(oaamom),w(oiax),
     .      w(ortab),w(ontab),w(oepair),xx)
        call dscal(3*nttab,alat,w(ortab),1)
        call fclose(ifi)
        call rx0('mmag after writing file '//jfilnm)

C --- Write Jij to file suitable to CVM input ---
      elseif (cmdopt('--wcvmj',7,0,outs)) then
        i = 7
        istyle = 6
        sclrsj = 1
        tolrsj = 1d-6
        jfilnm = 'cvmj'
        dc = outs(i+1:i+1)
        if (dc .ne. ' ') then
C     ... Return here to resume parsing for arguments
          j2 = i
  150     continue
          j2 = j2+1
          if (outs(j2:j2) .eq. dc) goto 150
          j1 = min(len(outs),j2)
          call nwordg(outs,0,dc//' ',1,j1,j2)
          if (j2 .ge. j1) then
            if (.false.) then
            elseif (outs(j1:j1+4) .eq. 'mode2')  then
              istyle = 7
            elseif (outs(j1:j1+3) .eq. 'scl=')  then
              j = 0
              i = parg('scl=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,sclrsj)
              if (i .le. 0) goto 159
            elseif (outs(j1:j1+3) .eq. 'tol=')  then
              j = 0
              i = parg('tol=',4,outs(j1:),j,len(outs(j1:)),
     .          dc//' ',1,1,k,tolrsj)
              if (i .le. 0) goto 159
            elseif (outs(j1:j1+2) .eq. 'fn=')  then
              jfilnm = outs(j1+3:j2)
            else
              goto 159
            endif
            goto 150
  159       call rxs('mmag: failed to parse --wcvmj options: ',outs)
          endif
        endif

        call info0(20,0,0,' ... writing r.s. hamiltonian to file '//
     .    jfilnm)
        ifi = fopna(jfilnm,-1,0)
        call dscal(3*nttab,1/alat,w(ortab),1)
        call dscal(nttab,sclrsj,w(oepair),1)
        i = iorsj(istyle,-ifi,nbasp,alat,plat,w(oaamom),w(oiax),
     .      w(ortab),w(ontab),w(oepair),xx)
        call dscal(3*nttab,alat,w(ortab),1)
        call fclose(ifi)
        call rx0('mmag after writing file '//jfilnm)

      endif

C ... Printout of initial conditions
      call pmmag2(nbas,w(oaamom),w(ontab),w(oiax),w(oepair))

      call defrr(owk,-(3*nbas+3))
      call upack('move nmodt modt ct',smove,nmodt,modt,ct,0,0)

      neq = nbas*2 + nmodt
      intgmd = igets('move prmint:1.1',smove)
      if (intgmd .eq. 1) then
        mx = igets('move prmint:5.5',smove)
        call defrr(ointwk,(neq*(4+mx)+5))
      else
        call rxi('mmag: unknown integration method',intgmd)
      endif
      nstep = -1
      call defrr(ofrc, -(nbas*3+3))

C --- Switches for restarting ---
      ifi = fopna('SV',-1,0)
      call poseof(ifi)
      tnow = 0
      j = 8
      if (cmdopt('-tstart=',j,0,outs))
     .  call rxx(.not.a2bin(outs,tnow,4,0,' ',j,len(outs)),
     .  'mmag: failed to parse -tstart=')
      j = 4

C     Switch --cont tells mmag to do the following:
C     1. Set tnow to last entry in save file and increment ttot
C     2. Read euler angles from eula-sv file
      if (cmdopt('--cont',j,0,outs)) then
        call rx('mmag not ready for new save')
        backspace ifi
        if (rdstrn(ifi,outs,len(outs),.false.)) then
        i = 0
        if (parstr(outs,'time=',len(outs)-5,5,'=',i,j)) then
        if (parg('time=',4,outs,i,len(outs),' ',1,1,k,tnow).eq.1) then
          call upack1('move tstot',smove,ttot)
          ttot = ttot + tnow
          call pack1('move tstot',smove,ttot)
          call info2(20,1,0,
     .      ' ... read Euler angles from file eula-sv;'//
     .      ' set start time to %d',tnow,0)
          goto 30
        endif
        endif
        endif
        call info0(20,0,0,
     .    ' ... read Euler angles from file eula-sv;'//
     .    ' unable to read start time from save file')
   30   continue

        jfi = fopna('eula-sv',-1,0)
        rewind jfi
        call ioeula(nbas,3,w(oeula),j,xsi,jfi)
        call fclose(jfi)
        call rxx(j.ne.1,'file eula-sv has neul>1')

      endif

C --- Micromagnetics dynamics ---
C  62 continue
      call mmdyn(smove,nbas,xsi,w(oeula),qspirl,w(owk),w(oepair),
     .  w(oipair),nstep,w(oiax),nttab,w(ointwk),w(ofrc),amag,
     .  w(oaamom),etot,tnow,nvarms,ifi)

C --- Calculate the pairwise forces and continue ---
C      if (nstep .gt. 0) then
C        etot = 0
C        call dpzero(w(ofrc),nbas*3+3)
C        call rx('mmag: need update mmpair')
C        call mmpair(nbas,nttab,w(oiax),w(oebil),w(oibil),w(oaamom),w(owk),
C     .   etot,amag,w(ofrc))
C        goto 62
C      endif

      end
      subroutine pmmag1(nbas,ntab,rtab,hpair,esite,ebil)
C- Distribute hpair into one-body and bilinear terms
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   rtab  :table of connecting vectors
Ci   hpair :pair hamiltonian.  On-site terms are those for which
Ci         :connecting vectors are zero
Co Outputs
Ci   esite :on-site terms
Ci   ebil  :bilinear terms
Cl Local variables
Cl         :
Cr Remarks
Cr   Assumes first entry in iax table is on-site term, if it exists
Cu Updates
Cu   12 Apr 03 Allow there to be no on-site term
Cu   24 Nov 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,ntab(nbas+1)
      double precision hpair(*),esite(nbas),ebil(*),rtab(3,*)
C ... Local parameters
      integer ib,it

      do  ib = 1, nbas
        it = ntab(ib)+1
        ebil(it) = hpair(it)
        if (rtab(1,it)**2+rtab(2,it)**2+rtab(3,it)**2.eq.0) ebil(it)=0
        esite(ib) = hpair(ntab(ib)+1)
        do  it = ntab(ib)+2, ntab(ib+1)
          ebil(it) = hpair(it)
        enddo
      enddo

      end

      subroutine pmmag2(nbas,amom,ntab,iax,hpair)
C- Printout of energy for collinear spins, also decp'n by site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   amom  :vector of moments
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   hpair :pair hamiltonian.  On-site terms, if present, are
Ci         :assumed to have zero amplitude.
Ci         :connecting vectors are zero
Cr Remarks
Cr   Assumes first entry in iax table is on-site term
Cu Updates
Cu   24 Nov 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,ntab(nbas+1),niax
      parameter (niax=10)
      integer iax(niax,*)
      double precision hpair(*),amom(nbas)
C ... Local parameters
      integer ib,jb,it,nttab,iprint
      double precision dsum,ei,kboltz,facJ0,e0,amag0
C     kboltz is in units eV K^-1
      parameter (kboltz=8.617d-5)

      nttab = ntab(nbas+1)

C      x = dsum(nttab,hpair,1) / 2
C      call info2(20,1,0,' Collinear energy 1/2*sum_ij J_ij = %,8;8g',x,
C     .  0)

      e0 = 0
      amag0 = dsum(nbas,amom,1)
      do  it = 1, nttab
        ib = iax(1,it)
        jb = iax(2,it)
        e0 = e0 + dsign(1d0,amom(ib)*amom(jb))*hpair(it)/2
C        print 333, it, e0
C  333   format(i5,f17.7)
      enddo
      call info2(20,1,0,
     .  ' Collinear energy 1/2*sum_ij sgn(s_i s_j) J_ij = %,8;8g'//
     .  '  sum s_i = %;6d',e0,amag0)

      if (iprint() .lt. 50) return

      facJ0 = 2d0/3d0 * 13.6d0 / kboltz
      call info0(0,0,0,
     .  '   site      J_00(mRy)     J_0(mRy)    2/3 J_0 (K)')

      do  ib = 1, nbas
C        e0 = hpair(ntab(ib)+1)
        e0 = 0
        ei = 0
        do  it = ntab(ib)+1, ntab(ib+1)
          ei = ei + hpair(it)
        enddo
        e0 = -0.5d0*e0
        ei = -0.5d0*ei
        call info5(0,0,0,' %,6i %;12,3D %;12,3D %;12,1D ',
     .    ib,e0*1000,ei*1000,facJ0*ei,0)
      enddo

      end

