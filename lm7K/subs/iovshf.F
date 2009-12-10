      subroutine iovshf(nbas,mode,commnt,ef0,ef,vconst,vshft,ifi)
C- I/O of GF or PGF potential array containing potential shifts
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   mode  :0 do nothing
Ci         :1 I/O Ef,vconst(1)
Ci         :2 I/O Ef,vconst(1..3)
Ci         :4 I/O vshft(1..nbas)
Ci         :5 I/O combination 1+4
Ci         :6 I/O combination 2+4
Ci   commnt:(file read) if commnt is nonblank, it is used
Ci         :as a format string for printout.
Ci         :(file write) if commnt is nonblank, it is
Ci         :included in the printout of global quantities
Ci   ef0   :(use in file read) reference Fermi level;
Ci         :used to adjust vconst when file value differs.
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  vconst:array of global potential shifts
Cio        :vconst(1) = potential shift
Cio        :vconst(2) = potential shift of L end region (PGF)
Cio        :vconst(3) = potential shift of R end region (PGF)
Cio  vshft :array of site potential shifts
Cio  ef    :Fermi level
Cr Remarks
Cr   Potential shifts are stored in two ways: global shifts
Cr   applying to all sites within a region, and site-specific
Cr   shifts.  Regions are defined in one of two contexts:
Cr   1.  Crystal case (one region)
Cr       shift is stored in vconst(1)
Cr   2.  Layer code (many PL), three regions.
Cr       Region C : global shift of active stored in vconst(1)
Cr       Region L : global shift stored in vconst(2)
Cr       Region R : global shift stored in vconst(3)
Cr   This subroutine does not mix the global and site-specific
Cr   regions and thus doesn't need to distinguish between the
Cr   two contexts.
Cr
Cr   iovshf reads in potential shifts (and class groupings) for every
Cr   site (class), and parameters relating to the shifts.
Cr   There is a header line which holds global quantities.
Cr      ef=# vconst=#   (crystal code)
Cr   or
Cr      ef=# vconst=# vconst(L)=# vconst(R)=# (layer code)
Cr   Following the header site- dependent shifts can be incorporated
Cr   in subseqent lines.  First a header line is required; then
Cr   lines of this form 
Cr      site shifts           <- header line
Cr       #1  #2               <- #1=site index ib, #2=shift
Cu Updates
Cu   21 Mar 03 completely redesigned; new argument list.
Cu   19 Feb 02 vshft(-5,-3,-1) no longer used for DOS
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) commnt
      integer nbas,mode,ifi
      double precision vconst(3),vshft(nbas),ef0,ef
      integer ib
C ... Local parameters
      integer i,j1,j2,ix(7),a2vec
      double precision xx(2),dasum
      integer k,swt(7),rdtok
      logical rdstrn
      character strn*128
C ... for rdfiln
C      integer mxchr,mxlev,lstsiz,ctlen
C      parameter (mxchr=20,mxlev=4,lstsiz=200,ctlen=120)
C      character*120 s,strn2,vnam(mxlev)*16,ctbl(mxchr,2)*(ctlen)
C      logical loop0(0:mxlev)
C      integer nr,nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
C     .  nlist(0:mxlev)

      if (mode .eq. 0) return

C --- File read ---
      if (ifi .gt. 0) then
        rewind ifi
        if (mode .ge. 4) call dpzero(vshft,nbas)

C        nr = 0
C        s = ' '
C   10   call rdfiln(ifi,'#{}%',mxlev,loop0,nlin,list,lstsiz,
C     .    ilist,nlist,vnam,ctbl,mxchr,strn2,s,len(s),nr)
C        if (nr .le. 0) goto 20
C        if (s .eq. ' ') goto 10

        if (.not. rdstrn(ifi,strn,len(strn),.false.)) return
        call stswt('init io,read stop',len(strn),swt)
        call stswt('token,opt,fixlen',0,swt)
        ef = ef0
        k = rdtok('ef=',strn,' ',' ',' ,',4,swt,1,1,ef,' ')
        if (k .eq. 0) backspace ifi
        if (k .ne. 0) then
        if (k .lt. 0) call rxs('iovsf: bad ef= in',strn)
        k = rdtok('vconst=',strn,' ',' ',' ,',4,swt,1,1,vconst,' ')
        if (k .lt. 0) call rxs('iovsf: bad vconst= in',strn)
        k =rdtok('vconst(L)=',strn,' ',' ',' ,',4,swt,1,1,vconst(2),' ')
        if (k .lt. 0) call rxs('iovsf: bad vconst(L)= in',strn)
        k =rdtok('vconst(R)=',strn,' ',' ',' ,',4,swt,1,1,vconst(3),' ')
        if (k .lt. 0) call rxs('iovsf: bad vconst(R)= in',strn)

        strn = ' IOVSHF: ' // commnt
        call strip(strn,j1,j2)
        if (j1+7 .le. j2) then
          call info(20,1,0,strn(1:j2),ef,vconst)
        endif
        endif
        if (abs(ef-ef0) .ne. 0) then
          call info2(1,0,0,' IOVSHF (warning): fermi level mismatch:'//
     .      ' input = %,7;7d file value=%,7;7d',ef0,ef)
        endif

        if (mode .ge. 4) then
          if (.not. rdstrn(ifi,strn,len(strn),.false.)) return
          call word(strn,1,j1,j2)
          if (j2 .lt. j1) return
          if (strn(j1:j2) .ne. 'site') return
          call word(strn,2,j1,j2)
          if (j2 .lt. j1) return
          if (strn(j1:j2) .ne. 'shifts') return
   10     continue
          if (.not. rdstrn(ifi,strn,len(strn),.false.)) return
          i = 0
          if (a2vec(strn,len(strn),i,4,' ',1,2,2,ix,xx) .ne. 2)
     .      call rx('IOVSHF : failed to parse line "'//strn//'%a"')
          ib = nint(xx(1))
          vshft(ib) = xx(2)
          if (ib .gt. nbas) call rxs('VSHFT : ib>nbas, line :',strn)
          goto 10
        endif

C --- File write to disk ---
      else
        strn = commnt

C       Write global shift data
        if (mod(mode,4) .eq. 1) then
          call awrit2(strn(1:10)//'%a ef=%,7;7d  vconst=%,7;7d',' ',80,
     .      -ifi,ef,vconst)
        elseif (mod(mode,4) .eq. 2) then
          call awrit4(strn(1:10)//'%a ef=%,7;7d  vconst=%,7;7d '//
     .      'vconst(L)=%,7;7d vconst(R)=%,7;7d',
     .      strn,120,-ifi,ef,vconst,vconst(2),vconst(3))
        endif

C       Write site shift data
        if (mode .ge. 4 .and. dasum(nbas,vshft,1) .ne. 0d0) then
          write(-ifi,'('' site shifts'')')
          do  ib = 1, nbas
            if (vshft(ib) .ne. 0) write(-ifi,334) ib, vshft(ib)
  334       format(i6,2f12.6)
          enddo
        endif
      endif

      end
C      subroutine vshfav(mode,nbl,pgplp,vshfi,vshfo)
CC- Collect average average vshfos, adjust vshfo(ib), or disperse avg
CC ----------------------------------------------------------------
CCi mode  :See Remarks about distinctions between crystal subblocks
CCi       :1s digit
CCi         1  Disperse average shift = vshfo(-6,-4,-2) into individual
CCi            sites (which of vshfo(-6,-4,-2) to use depends on which
CCi            region site ib belongs to; see Remarks) and subsequently
CCi            set vshfo(-6,-4,-2) to zero.
CCi         2  This mode is used to update vshfo for potential shifts
CCi            generated by a band program.
CCi            It computes change <vshfo(ib)-vshfi(ib)> averaging over
CCi            ib within a region (see Remarks); add < > to
CCi            vshfo(-6,-4,-2) (depending on region) and subtract it
CCi            from individual sites.
CCi         0  This mode is similar to mode 2 (see description above).
CCi            The difference with mode 2 is that instead of averaging
CCi            <vshfo(ib)-vshfi(ib)> over a region,
CCi            <vshfo(ib)> is averaged.
CCi       10s digit
CCi         1  initially copy o to i
CCi         2  initially copy i to o
CCi nbl   :number of crystal subblocks (susite.f)
CCi pgplp :index and dimensioning information for crystal subblocks.
CCi        The meaning of pgplp depends on the context; see Remarks
CCi        and susite.f
CCio vshfi:input potential shift
CCio vshfo:output potential shift
CCr Remarks
CCr   This code is used in two distinct contexts:
CCr   1.  Crystal case (one PL)
CCr       There is one PL; nbl=1 and a constant average is stored in
CCr       vshft(-2).
CCr   2.  Layer code (many PL, three regions.
CCr       Region L corresponds to PL=-1; average is stored in vshft(-6)
CCr       Region R corresponds to PL=nbl; average is stored in vshft(-4)
CCr       Region C corresponds to -1<PL<nbl;  average stored in vshft(-2)
CCr   The contexts are distinguished by whether pgplp(1,-1) .ge. 0
CC ----------------------------------------------------------------
C      implicit none
CC Passed variables
C      integer mode,nbl,pgplp(6,-1:nbl)
C      double precision vshfi(-7:18),vshfo(-7:18)
CC Local variables
C      logical lpad
C      integer ibl,ib,ib1,ib2,nb,k,mode0,mode1
C      double precision vshfa
C
C      stop 'oops!'
C      mode0 = mod(mode,10)
C      mode1 = mod(mode/10,10)
C
C      if (mode1 .ne. 0) then
CC   ... Call to gtibpl picks up nbasp
C        call gtibpl(nbl,nbl,pgplp,ib1,ib2)
C        if (mode1 .eq. 1) call dcopy(8+ib2,vshfo,1,vshfi,1)
C        if (mode1 .eq. 2) call dcopy(8+ib2,vshfi,1,vshfo,1)
C      endif
C
C      lpad = pgplp(1,-1) .ge. 0
C      vshfa = 0
C      nb = 0
C      do  10  ibl = -1, nbl
C        if (.not. lpad .and. (ibl .eq. -1 .or. ibl .eq. nbl)) goto 10
C        call gtibpl(ibl,nbl,pgplp,ib1,ib2)
C        k = 0
C        if (ibl .eq. -1) k = -6
C        if (ibl .eq. nbl) k = -4
C        if (ibl .eq. nbl-1) k = -2
C        do  12  ib = ib1, ib2
C          nb = nb+1
C          vshfa = vshfa + vshfo(ib)
C          if (mode0 .eq. 2) then
C            vshfa = vshfa - vshfi(ib)
C          endif
C   12   continue
C
C        if (k .lt. 0) then
CC     ... Collect average, combined with vshfo(-6,-4,-2)
C          if (mode0 .eq. 0 .or. mode0 .eq. 2) then
C            vshfa = vshfa / nb
C            vshfo(k) = vshfo(k) + vshfa
C            do  20  ib = ib2-nb+1, ib2
C   20       vshfo(ib) = vshfo(ib) - vshfa
C            vshfa = 0
C            nb = 0
CC     ... Collect average vshfo(-6,-4,-2) into individual sites
C          else
C            vshfa = vshfo(k)
C            vshfo(k) = 0
C            do  30  ib = ib2-nb+1, ib2
C   30       vshfo(ib) = vshfo(ib) + vshfa
C            nb = 0
C          endif
C        endif
C   10 continue
C
C      end

