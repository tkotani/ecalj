
      subroutine slist(lstyle,strlst,slabl,z,nspec,nlist,list)
C- Generates a list of species from a string specification
C ----------------------------------------------------------------
Ci Inputs
Ci   lstyle: style of specification; see Remarks
Ci   strlst:  string specifying list of species
Ci   slabl :list of species labels
Ci   nspec :number of species.
Ci   z     :table nuclear charge by species
Co Outputs
Co   nlist :number of species in list
Co   list  :a list of species satisfying specifications
Cr Remarks
Cr *Syntax of strlst: depends on one of three styles (lstyle)
Cr
Cr *lstyle=1 : a list of integers; see mkilst.f for complete syntax.
Cr             Example: '1,4:6,11' generates a list of five numbers,
Cr             1,4,5,6,11.
Cr
Cr *lstyle=2 : the list is specified according to an expression.
Cr             The expression can involve the species index is and
Cr             atomic number z.  Any species satisfying expression is
Cr             included in the list.  Example:  'is<6&z==14'
Cr
Cr *lstyle=3 : strlst is a list of species by name; entries are
Cr           : separated by commas.  Each entry must be in slabl.
Cu Updates
Cu   13 Sep 01 adapted from clist
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lstyle,nlist,nspec,list(1)
      character*(*) strlst
      character*8 slabl(nspec)
      double precision z(nspec)
C ... Local parameters
      integer iv0,ival,is,i,j,ls,j1,j2
      logical a2bin,sw
      character alabl*8
      integer,allocatable:: ilst(:)
C ... Heap
c      integer w(1)
c      common /w/ w
C ... External calls
ckino del defi,rlse
      external clrsyv,hunti,ishell,lodsyv,mkils0,mkilst,numsyv,
     .nwordg,rxi,rxs,tokmat

      ls = len(strlst)
      nlist = 0
C     call iinit(list,nspec)
c      goto (10,20,30) lstyle
      select case(lstyle*10)
C -- lstyle=1 ---
c   10 continue
      case(10)
      call mkils0(strlst,nlist,i)
c        call defi(oilst, nlist)
      allocate(ilst(nlist))
      call mkilst(strlst,nlist,ilst) !w(oilst))
      if (nlist .eq. 0) return
      call ishell(nlist,ilst) !w(oilst))
      list(1) = ilst(1) !w(oilst)
      j = 1
      do  i = 2, nlist
c          if (w(oilst+i-1) .gt. list(j)
c     .    .and. w(oilst+i-1) .le. nspec) then
c            j = j+1
c            list(j) = w(oilst+i-1)
        if (ilst(i) .gt. list(j)
     .    .and. ilst(i) .le. nspec) then
          j = j+1
          list(j) = ilst(i)
        endif
      enddo
      nlist = j
c        call rlse(oilst)
      deallocate(ilst)
      return

C --- lstyle=2 ---
c   20 continue
      case(20)
      call numsyv(iv0)
      nlist = 0
      do  42  is = 1, nspec
        call lodsyv('is',1,dble(is),ival)
        call lodsyv('z',1,z(is),ival)
        i = 0
        if (a2bin(strlst,sw,0,0,strlst(ls:ls),i,ls)) then
          if (sw) then
            nlist = nlist+1
            list(nlist) = is
          endif
C   ... Abort if a2bin can't parse expression
        else
          call rxs('slist: failed to parse',strlst)
        endif
   42 continue
      call clrsyv(iv0)
      return

C --- lstyle=3 ---
c   30 continue
      case(30)
      nlist = 0
      j2 = -1
C ... Return here to parses next species
   31 continue
      j1 = j2+2
      call nwordg(strlst,1,', ',1,j1,j2)
      if (j1 .gt. j2) return
      alabl = strlst(j1:j2)
      call tokmat(alabl,slabl,nspec,8,' ',is,i,.false.)
      if (is .eq. -1 .or. is .ge. nspec)
     .  call rxs('slist: failed to parse species list ',strlst)
      is = is+1
      i = 0
      call hunti(list,nlist,is,0,i)
      if (i .ge. nlist .or. list(i+1) .ne. is) then
        nlist = nlist+1
        list(nlist) = is
        call ishell(nlist,list)
      endif
C     print *, j1,j2,alabl,is,'   ',list(1:nlist)
      goto 31

      case default
      call rxi('slist: bad style',lstyle)
      return

      end select
      end


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs
c$$$  subroutine sumlst(lopt,nchmx,nbas,ng,ssite,sspec,sopts,mode,
c$$$     .nsites,lsites,lmxch,nchan,lchan,nll)
c$$$      use m_struc_def  !Cgetarg
c$$$C- Set up list of sites and mode for Mulliken or partial DOS projection
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs:
c$$$Ci   lopt  :1s digit: choice of l channels
c$$$Ci         :  0 for partial DOS (lcutoff to lmxa; see Remarks)
c$$$Ci         :  1 for Mulliken analysis (lcutoff to lmxb; see Remarks)
c$$$Ci         :10s digit what to return
c$$$Ci         :  0 return all outputs
c$$$Ci         :  1 return nsites and lmxch only
c$$$Ci              (NB: if sopts contains site list, lsites also returned)
c$$$Ci              (NB: if sopts contains nl=, value of nll is altered)
c$$$Ci   nchmx :(lopt=1 only) dimensions array lchan
c$$$Ci   nbas  :size of basis
c$$$Ci   ng    :# of space group ops --- checks the m-resolved case
c$$$Ci   ssite :struct for site-specific information; see routine usite
c$$$Ci     Elts read: spec
c$$$Ci     Stored:    *
c$$$Ci     Passed to: *
c$$$Ci   sspec :struct for species-specific information; see routine uspec
c$$$Ci     Elts read: lmxb lmxa
c$$$Ci     Stored:    *
c$$$Ci     Passed to: *
c$$$Ci   sopts  from command line (see Remarks) with
c$$$Ci           the --mull or --pdos stripped
c$$$Co Outputs:
c$$$Co   mode  :(see remarks)
c$$$Co   nsites:(number of sites)
c$$$Co   lsites:(list of sites)
c$$$Co   nchan :(number of DOS channels)
c$$$Co   lmxch :(maximum basis l for all of the atoms in the list)
c$$$Co   lchan :(array of channel numbers for each site,l,m) see Remarks
c$$$Co         :(lopt=1 only; not used for lopt=0)
c$$$Co   nll   :if nonzero, use lmax = nll-1 for all sites in making
c$$$Co         :dos channels.  NB: input value of nll passed by caller
c$$$Co         :is its default value.
c$$$Cr Remarks
c$$$Cr   sumlst:options passed in sopts
c$$$Cr          [:mode=#][nl=#][:sites=list]]
c$$$Cr   modes are 0: all sites atom-resolved (default --mull)
c$$$Cr             1: all sites l-resolved
c$$$Cr             2: all sites lm-resolved (default --pdos)
c$$$Cr             Adding 3 to mode flags sumlst that a site-list is
c$$$Cr             to follow.  This is not needed, however, since supplying
c$$$Cr             a site-list automatically will cause 3 to be added.
c$$$Cr
c$$$Cr   The --mull (lopt=1) and --pdos (lopt=0) modes are similar;
c$$$Cr   however lopt=1 case is complicated by the possibility of multiple
c$$$Cr   basis functions per l channel, and also different lmxb in different
c$$$Cr   kappa channels.  To handle that case, sumlst creates a list marking
c$$$Cr   which DOS channel each site, l and m belongs.  sumlst returns
c$$$Cr   this array in the form that mullmf can reference as
c$$$Cr   lchan(ib) lchan(l,ib) or lchan(lm,ib) depending on the mode
c$$$Cr   and lchan will be dimensioned using the maximum l, lmxch:
c$$$Cr   lchan(ib) lchan(lmxch+1,ib) or lchan((lmxch+1)**2,ib)
c$$$Cr   depending on the mode. The ordering of ib is either as the
c$$$Cr   internal list of sites (mode .lt. 3) or as the list returned
c$$$Cr   from the cmd-line string by mkilst. (see mchan.f)
c$$$Cu Updates
c$$$Cu   30 Aug 04 Changed first argument to integer lopt
c$$$Cu   18 Jan 02 (MvS) redesigned switches for compatibility with site
c$$$Cu                   lists in other context
c$$$Cu   27 Mar 01 (MvS) extended to handle both --mull and --pdos
c$$$Cu   20 Mar 01 Written by ATP
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C Passed Parameters
c$$$      character*(*) sopts
c$$$      integer lopt,nchmx,nbas,ng,nsites,mode,lsites(*),nchan,lmxch,
c$$$     .lchan(*),nll
c$$$      type(s_site)::ssite(*)
c$$$      type(s_spec)::sspec(*)
c$$$
c$$$C Local Variables
c$$$      integer mxmode
c$$$      parameter (mxmode = 6)
c$$$      character dc*1,modstr(0:mxmode-1)*128
c$$$      integer lmax,isite,lmdim,ltop,ilm
c$$$      integer i,j,ipr,iprmin,ib,is,iprint,lgunit,stdo,igetss,isw,j2,j1,
c$$$     .m,iv(10),parg,lstyle,lopt0,lopt1
c$$$      double precision z
c$$$      data modstr /'all sites atom-projected',
c$$$     .'all sites l-projected','all sites lm-projected',
c$$$     .'site list atom-projected',
c$$$     .'site list l-projected','site list lm-projected' /
c$$$      logical:: isanrg, l_dummy_isanrg
c$$$      character(256)::ttt
c$$$      lopt0 = mod(lopt,10)
c$$$      lopt1 = mod(lopt/10,10)
c$$$      ipr = iprint()
c$$$      iprmin = 30
c$$$      mode = -1
c$$$      dc = sopts(1:1)
c$$$      stdo = lgunit(1)
c$$$      mode = 2
c$$$      nsites = 0
c$$$      if (lopt0 .eq. 1) mode = 0
c$$$      if (dc .ne. ' ') then
c$$$C   ... Return here to resume parsing for arguments
c$$$        j2 = 0
c$$$   10   continue
c$$$        j2 = j2+1
c$$$        if (sopts(j2:j2) .eq. dc) goto 10
c$$$        j1 = min(len(sopts),j2)
c$$$        call nwordg(sopts,0,dc//' ',1,j1,j2)
c$$$        if (j2 .ge. j1) then
c$$$          if (.false.) then
c$$$C         DOS mode; see Remarks
c$$$          elseif (sopts(j1:j1+4) .eq. 'mode=')  then
c$$$            m = 0
c$$$            i = parg('mode=',2,sopts(j1:),m,len(sopts(j1:)),
c$$$     .      dc//' ',1,1,iv,mode)
c$$$            if (i .le. 0) goto 999
c$$$Ckino isanrg is logical function,             call isanrg(mode,0,5,' sumlst:','mode',.true.)
c$$$            l_dummy_isanrg=isanrg(mode,0,5,' sumlst:','mode',.true.)
c$$$            if (ipr .ge. 10) then
c$$$              if ((mode .eq. 2 .or. mode .eq. 5) .and. ng .gt. 1)
c$$$     .        write(stdo,*)' **WARNING** sumlst: for lm-decomposition suppress symops'
c$$$            endif
c$$$C         Number of energy points.  Not used now
c$$$          elseif (sopts(j1:j1+2) .eq. 'nl=')  then
c$$$            m = 0
c$$$            i = parg('nl=',2,sopts(j1:),m,len(sopts(j1:)),
c$$$     .      dc//' ',1,1,iv,nll)
c$$$            if (i .le. 0) goto 999
c$$$C         Site list
c$$$          elseif (sopts(j1:j1+5) .eq. 'sites=')  then
c$$$            if (mode .lt. 3) mode = mode+3
c$$$C           We need to make slabl,z before using other lstyle
c$$$            lstyle = 1
c$$$            call slist(lstyle,sopts(j1+6:j2),' ',z,nbas,nsites,lsites)
c$$$          else
c$$$            goto 999
c$$$          endif
c$$$          goto 10
c$$$        endif
c$$$      endif
c$$$      if (ipr .ge. iprmin+1) write(stdo,1)
c$$$    1 format (' sumlst: Site   l_max')
c$$$      if (mode .lt. 3) nsites = nbas
c$$$
c$$$      if (nsites .eq. 0) call rxi
c$$$     .('sumlst: --sites=list is required input for mode =',mode)
c$$$
c$$$C --- Generate lsites,lmxch ---
c$$$      lmxch = -1
c$$$      do  isite = 1, nsites
c$$$        if (mode .lt. 3) then
c$$$          ib = isite
c$$$          if (lopt1 .eq. 0) then
c$$$            lsites(isite) = ib
c$$$          endif
c$$$        else
c$$$          ib = lsites(isite)
c$$$        endif
c$$$Ckino isanrg is logical function,         call isanrg(ib,1,nbas,' sumlst:','site',.true.)
c$$$        l_dummy_isanrg=isanrg(ib,1,nbas,' sumlst:','site',.true.)
c$$$
c$$$        is=ssite(ib)%spec
c$$$
c$$$        if ( lopt0 .eq. 1 ) lmax = int(sspec(is)%lmxb)
c$$$
c$$$        if ( lopt0 .eq. 0 ) lmax = int(sspec(is)%lmxa)
c$$$
c$$$        if (nll .gt. 0) lmax = min(nll-1,lmax)
c$$$        lmxch = max(lmxch,lmax)
c$$$      enddo
c$$$      if (lopt1 .eq. 1) return
c$$$
c$$$C --- Generate nchan,lchan ---
c$$$      i = 0
c$$$      do  isite = 1, nsites
c$$$        ib = lsites(isite)
c$$$
c$$$        is=ssite(ib)%spec
c$$$
c$$$        if ( lopt0 .eq. 1 ) lmax = int(sspec(is)%lmxb)
c$$$
c$$$        if ( lopt0 .eq. 0 ) lmax = int(sspec(is)%lmxa)
c$$$
c$$$        if (nll .gt. 0) lmax = min(nll-1,lmax)
c$$$        if (ipr .gt. iprmin) write(stdo,2) ib, lmax
c$$$    2   format (9x,i4,5x,i1)
c$$$        if (mode .eq. 0 .or. mode .eq. 3) then
c$$$          ltop = 1
c$$$          lmdim = 1
c$$$        elseif (mode .eq. 1 .or. mode .eq. 4) then
c$$$          ltop = lmax+1
c$$$          lmdim = lmxch+1
c$$$        elseif (mode .eq. 2 .or. mode .eq. 5) then
c$$$          ltop = (lmax+1)**2
c$$$          lmdim = (lmxch+1)**2
c$$$        endif
c$$$        if (lopt0 .eq. 1) then
c$$$          do  ilm = 1, ltop
c$$$            i = i + 1
c$$$            if (i.gt.nchmx) call rxi(' sumlst: bndfp needs nchmx > ',i)
c$$$            call mchan(lmdim,0d0,0d0,0,nsites,0,isite,ilm,-1,i,lchan)
c$$$          enddo
c$$$        else
c$$$          i = i+ltop
c$$$        endif
c$$$      enddo
c$$$      nchan = i
c$$$C --- Printout --
c$$$      if (ipr .ge. iprmin) then
c$$$         call strip(modstr(mode),i,j)
c$$$         if(lopt0==0) ttt=' sumlst:  Mulliken mode '
c$$$         if(lopt0==1) ttt=' sumlst:  Partial DOS mode '
c$$$         write(stdo,"(a,i5,a,i5,' sites ',i5,' channels')")ttt,mode,modstr(mode)(i:j),nsites,nchan
c$$$      endif
c$$$      return
c$$$  999 continue
c$$$      call rxs('sumlst: failed to parse options in ',sopts(j1:))
c$$$      end subroutine sumlst


c$$$      logical function pars1v(strn,size,nam,term,cast,val)
c$$$C- Parses a string for a specified variable
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   recrd(0:*): string recrd is parsed from i to size-1
c$$$Ci   iopt:   0: var declaration of existing variable is supressed
c$$$Ci           1: var declaration of existing variable supersedes
c$$$Ci   mxdecl: parsing ends if number of declarations exceeds mxdecl
c$$$Co Outputs
c$$$Co   i:      last character parsed
c$$$Cr Remarks
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$C Passed parameters
c$$$      integer size,cast
c$$$      character*(*) strn,nam,term*1
c$$$      double precision val
c$$$C Local parameters
c$$$      integer i,j
c$$$      logical parstr,a2bin
c$$$
c$$$      pars1v = .false.
c$$$      i = 0
c$$$      if (.not. parstr(strn,nam,size,len(nam),term,i,j)) return
c$$$      if (a2bin(strn,val,cast,0,' ',j,size)) pars1v = .true.
c$$$      end
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

      
c$$$#define USEOPTIONAL
c$$$
c$$$      module m_ioden2
c$$$      contains
c$$$      subroutine ioden2(mode,sopts,nsp,plat,alat,nbas,pos,z,k1,k2,k3,
c$$$     .ngabc,smrho,wk,nwk,modrhg)
c$$$
c$$$C- Kernel called by ioden
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   mode  :0 just return these arguments parsed from sopts:
c$$$Ci         :  modrhg
c$$$Ci         :1 copy and save density in appropriate plane
c$$$Ci   sopts :string containing options; see ioden above.
c$$$Ci   plat  :primitive lattice vectors, in units of alat (for printout)
c$$$Ci   alat  :lattice constant
c$$$Ci   nbas  :number of atoms in unit cell
c$$$Ci   pos   :atomic positions in unit cell
c$$$Ci   z     : atomic number of atoms in cell
c$$$Ci   nsp   :number of spin channels
c$$$Ci   k1..3 :dimensions smrho
c$$$Ci   ngabc :number of divisions in smrho along plat(1..3)
c$$$Ci   smrho :smooth density on uniform mesh
c$$$Ci   wk    :2D work array
c$$$Ci   nwk   :dimensions wk
c$$$Co Outputs
c$$$Co   modrhg:controls what part of core density is added
c$$$Co   The projection of smrho onto a plane is written to disk
c$$$Co   or if option is g3d the full 3d grid is written to disk
c$$$Cl Local variables
c$$$Cl   lic   :lic(1..3,1..2) = step length in each vector axes 1..2
c$$$Cl         :lic(4,1..2)    = number of points along axes 1..2
c$$$Cl         :lic(5..7,1)    = starting element for each lattice vector
c$$$Cr Remarks
c$$$Cr   Choosing p1, p2, p3 sets the origin for that plat and also
c$$$Cr   creates a default set of vectors, choosing the other lattice
c$$$Cr   vectors in the system.  Thus:
c$$$Cr   p1=# => lic(1..7,1) = (0 1 0 nb 1 1 1) lic(1..7,2) = (0 0 1 nc **)
c$$$Cr   p2=# => lic(1..7,1) = (0 0 1 nc 1 1 1) lic(1..7,2) = (1 0 0 na **)
c$$$Cr   p3=# => lic(1..7,1) = (1 0 0 na 1 1 1) lic(1..7,2) = (0 1 0 nb **)
c$$$Cu Updates
c$$$Cu   09 Feb 01 Revised how a plane is defined.
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character*(*) sopts
c$$$      integer mode,ngabc(3),k1,k2,k3,nwk,nsp,modrhg,nbas
c$$$#ifdef USEOPTIONAL
c$$$      double complex ,optional:: smrho(k1,k2,k3,nsp)
c$$$      double precision plat(3,3),alat,pos(3,nbas)
c$$$      real(8),optional:: wk(nwk,nwk,nsp)
c$$$#else
c$$$      double complex smrho(k1,k2,k3,nsp)
c$$$      double precision plat(3,3),wk(nwk,nwk,nsp),alat,pos(3,nbas)
c$$$#endif
c$$$      double precision z(nbas)
c$$$C ... Local parameters
c$$$      logical lok
c$$$      integer lgunit,stdo,j2,j1,fopn,ifi,i
c$$$      integer lic(7,2),j,k,iv(3),jv(3),kv(4),a2vec,isp
c$$$      double precision vecs(3,2),vv,ddot
c$$$      character*120 dc*1, fn, prfmt*40, sout
c$$$      logical lg3d
c$$$      integer i1,i2,i3
c$$$
c$$$      logical:: l_dummy_isanrg,isanrg
c$$$
c$$$ccccccccccccccccccccccccccccc
c$$$c      print *,' zzxxx=',z
c$$$c      stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzz'
c$$$ccccccccccccccccccccccccccccc
c$$$
c$$$      stdo = lgunit(1)
c$$$      dc = sopts(1:1)
c$$$      fn = 'smrho'
c$$$C     prfmt = '(8f10.6)'
c$$$      prfmt = '(8e14.6)'
c$$$      lg3d = .false.
c$$$C     Default: origin at (1,1,1), (na+1,nb+1) points along 1st,2nd axes
c$$$      call iinit(lic,14)
c$$$      call ivset(lic,5,7,1)
c$$$      lic(1,1) = 1
c$$$      lic(4,1) = ngabc(1)+1
c$$$      lic(2,2) = 1
c$$$      lic(4,2) = ngabc(2)+1
c$$$
c$$$      if (dc .ne. ' ') then
c$$$C ... Return here to resume parsing for arguments
c$$$        j2 = 0
c$$$   10   continue
c$$$        j2 = j2+1
c$$$        if (sopts(j2:j2) .eq. dc) goto 10
c$$$        j1 = min(len(sopts),j2)
c$$$        call nwordg(sopts,0,dc//' ',1,j1,j2)
c$$$        if (j2 .ge. j1) then
c$$$          if (.false.) then
c$$$          elseif (sopts(j1:j1+2) .eq. 'fn=')  then
c$$$            if (j1+3 .gt. j2) goto 99
c$$$            fn = sopts(j1+3:j2)
c$$$C   ... option g3d
c$$$          elseif (sopts(j1:j1+2).eq. 'g3d') then
c$$$            lg3d = .true.
c$$$C   ... option o=#,#,#
c$$$          elseif (sopts(j1:j1+1) .eq. 'o=')  then
c$$$            if (j1+2 .gt. j2) goto 99
c$$$            i = j1+1
c$$$            if (a2vec(sopts,j2,i,2,', '//dc,3,2,3,kv,iv) .ne. 3)
c$$$     .      goto 99
c$$$            lic(5,1) = mod(iv(1)+ngabc(1),ngabc(1))+1
c$$$            lic(6,1) = mod(iv(2)+ngabc(2),ngabc(2))+1
c$$$            lic(7,1) = mod(iv(3)+ngabc(3),ngabc(3))+1
c$$$C   ... option l[12]=#,#,#,#
c$$$          elseif (sopts(j1:j1+2) .eq. 'l1=' .or.
c$$$     .    sopts(j1:j1+2) .eq. 'l2=')  then
c$$$            if (j1+3 .gt. j2) goto 99
c$$$            i = 0
c$$$            call chrps2(sopts(j1+1:),'12',2,0,i,iv)
c$$$C         this check should never be necessary
c$$$Ckino isanrg is logical function, C         call isanrg(iv,1,2,' ','iv in ioden ... bug ..',.true.)
c$$$C         l_dummy_isanrg=isanrg(iv,1,2,' ','iv in ioden ... bug ..',.true.)
c$$$            i = j1+2
c$$$            j = iv(1)
c$$$            k = a2vec(sopts,j2,i,2,', '//dc,3,2,4,kv,lic(1,j))
c$$$            if (k .ne. 3 .and. k .ne. 4) goto 99
c$$$          elseif (sopts(j1:j1+4) .eq. 'core=') then
c$$$            if (j1+5 .gt. j2) goto 99
c$$$            i = j1+4
c$$$            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,modrhg) .ne. 1)
c$$$     .      goto 99
c$$$          else
c$$$            call rxs('ioden: unrecognised option ... ',sopts(j1:j2))
c$$$          endif
c$$$          goto 10
c$$$        endif
c$$$      endif
c$$$      if (mode .eq. 0) return
c$$$
c$$$C --- Printout ---
c$$$      if (lg3d) then
c$$$        ifi = fopn(fn)
c$$$        rewind ifi
c$$$        call awrit3('      Writing smooth density to file'//fn//
c$$$     .  ' %a : full 3d grid (%i,%i,%i).',' ',80,stdo,k1,k2,k3)
c$$$        do isp = 1, nsp
c$$$C           call awrit4('%% 3D block data %i %i %i, spin=%i',' ',80,ifi,
c$$$C     .     k1,k2,k3,isp)
c$$$          write(ifi,'("CRYSTAL")')
c$$$          write(ifi,'("PRIMVEC")')
c$$$          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
c$$$     .    ,i2=1,3)
c$$$          write(ifi,'("PRIMCOORD")')
c$$$          write(ifi,'(2i5)') nbas,1
c$$$          do i = 1, nbas
c$$$            write(ifi,'(i4,2x,3f10.5)') int(z(i)),
c$$$     .      (pos(i2,i)*alat*0.529177208,i2=1,3)
c$$$          enddo
c$$$          write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
c$$$          write(ifi,'("charge_density_spin_",i1)') isp
c$$$          write(ifi,'("BEGIN_DATAGRID_3D_isp_",i1)') isp
c$$$          write(ifi,'(3i4)') k1,k2,k3
c$$$          write(ifi,'(3f10.5)') 0.,0.,0.
c$$$          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
c$$$     .    ,i2=1,3)
c$$$          write(ifi,'(8e14.6)')  
c$$$     .    (((dble(smrho(i1,i2,i3,isp)),i1=1,k1),i2=1,k2),i3=1,k3)
c$$$          write(ifi,'("END_DATAGRID_3D_isp_",i1)') isp
c$$$          write(ifi,'("END_BLOCK_DATAGRID_3D")')
c$$$        enddo
c$$$        call fclose(ifi)
c$$$        return
c$$$      else
c$$$        call awrit3('         Writing smooth density to file '//fn//
c$$$     .  '%a : origin at (%i,%i,%i).',' ',80,stdo,lic(5,1)-1,
c$$$     .  lic(6,1)-1,lic(7,1)-1)
c$$$        call dpzero(vecs,3)
c$$$        do  i = 1, 3
c$$$          call daxpy(3,dble(lic(4+i,1)-1)/ngabc(i),plat(1,i),1,vecs,1)
c$$$        enddo
c$$$        call ivset(kv,1,3,2)
c$$$        call shorps(1,plat,kv,vecs,vecs)
c$$$        write(stdo,368) (vecs(i,1), i=1,3)
c$$$  368   format(9x,'In cartesian coordinates, origin is',3f10.6)
c$$$
c$$$        call dpzero(vecs,6)
c$$$        do  k = 1, 2
c$$$          do  i = 1, 3
c$$$            vv = dble((lic(4,k)-1)*lic(i,k))/dble(ngabc(i))
c$$$            call dpadd(vecs(1,k),plat(1,i),1,3,vv)
c$$$          enddo
c$$$          vv = dsqrt(ddot(3,vecs(1,k),1,vecs(1,k),1))
c$$$          sout = ' '
c$$$          call awrit8('%9fv%i: (%i pts) = (%i,%i,%i)p -> '//
c$$$     .    '(%,6;6d,%,6;6d,%,6;6d)',sout,80,0,k,lic(4,k),
c$$$     .    lic(1,k),lic(2,k),lic(3,k),vecs(1,k),vecs(2,k),vecs(3,k))
c$$$          call awrit1('%a l=%,6;6d',sout,-80,-stdo,vv)
c$$$        enddo
c$$$
c$$$C ... Sanity checks
c$$$        lok = .true.
c$$$        lok = lok .and. lic(4,1) .gt. 0
c$$$        lok = lok .and. lic(4,2) .gt. 0
c$$$        if (.not. lok) call fexit2(-1,1,' Exit -1 ioden: number of '//
c$$$     .  'points along axes (%i,%i) are not > 0',lic(4,1),lic(4,2))
c$$$        if (nwk.lt.lic(4,1) .or. nwk.lt.lic(4,2)) call rx('increase nwk')
c$$$
c$$$C --- Copy points to wk from plane ---
c$$$        call icopy(3,lic(5,1),1,kv,1)
c$$$        do  k = 1, lic(4,2)
c$$$          call icopy(3,kv,1,iv,1)
c$$$C         PRINT *, K, KV(1), KV(2), KV(3)
c$$$          do  i = 1, lic(4,1)
c$$$C         PRINT *, ' ',I, IV(1), IV(2), IV(3)
c$$$            do  j = 1, 3
c$$$              jv(j) = mod(iv(j)-1,ngabc(j))+1
c$$$   44         if (jv(j) .le. 0) then
c$$$                jv(j) = jv(j) + ngabc(j)
c$$$                goto 44
c$$$              endif
c$$$            enddo
c$$$            if (jv(1) .le. 0 .or. jv(2) .le. 0 .or. jv(3) .lt. 0) then
c$$$              call rx('bug in ioden')
c$$$            endif
c$$$            do  isp = 1, nsp
c$$$              wk(i,k,isp) = dble(smrho(jv(1),jv(2),jv(3),isp))
c$$$            enddo
c$$$C          if (i .eq. 1 .or. i .eq. 16) then
c$$$C          PRINT *, I,K, JV(1), JV(2), JV(3),WK(I,K)
c$$$C          endif
c$$$            do  j = 1, 3
c$$$              iv(j) = iv(j) + lic(j,1)
c$$$            enddo
c$$$          enddo
c$$$          do  j = 1, 3
c$$$            kv(j) = kv(j) + lic(j,2)
c$$$          enddo
c$$$        enddo
c$$$
c$$$C --- Write the file ---
c$$$        ifi = fopn(fn)
c$$$        rewind ifi
c$$$        do  isp = 1, nsp
c$$$          call ywrm(0,' spin 1',1,ifi,prfmt,wk(1,1,isp),0,nwk,lic(4,1),
c$$$     .    lic(4,2))
c$$$        enddo
c$$$        call fclose(ifi)
c$$$        return
c$$$      endif
c$$$
c$$$   99 continue
c$$$      call rxs('ioden: failed to parse option ... ',sopts(j1:j2))
c$$$
c$$$      end subroutine ioden2
c$$$      end module m_ioden2
c$$$
c$$$      subroutine ioden ( sopts , slat , ssite , sspec , sv_p_orhat 
c$$$     . , smrho )
c$$$
c$$$
c$$$      use m_struc_def  !Cgetarg
c$$$      use m_globalvariables
c$$$      use m_ioden2
c$$$C- File I/O charge density on a uniform mesh in a plane or full 3d mesh
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   sopts :character string specifying plane and other options.
c$$$Ci         :Specifications and options are separated by a delimiter,
c$$$Ci         :which is the first character in sopts.
c$$$Ci         :
c$$$Ci         :if option is g3d then the density is written to
c$$$Ci         :disk on the full 3d grid.   Otherwise:
c$$$Ci         :The density is written to disk for a uniform of mesh of
c$$$Ci         :points in one plane.  This information is specified by three
c$$$Ci         :groups of numbers: the origin, a first direction vector with
c$$$Ci         :its number of points, and a second direction vector with its
c$$$Ci         :number of points.
c$$$Ci         :
c$$$Ci         :At present, these points must be taken from the points on
c$$$Ci         :the smooth mesh, smrho.  In this case, all three groups of
c$$$Ci         :information are sets of integers.  For example,
c$$$Ci         :specify the origin by three numbers:
c$$$Ci         :    o=#1,#2,#3
c$$$Ci         :The point (#1,#2,#3) corresponds to the Cartesian coordinate
c$$$Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
c$$$Ci         :where (n1,n2,n3) are the number of divisions in the
c$$$Ci         :mesh along the three lattice vectors (p1,p2,p3).
c$$$Ci         :o=0,0,0 corresponds to the origin.
c$$$Ci         :
c$$$Ci         :Specify the direction vectors by
c$$$Ci         :    l1=#1,#2,#3[,#4]
c$$$Ci         :    l2=#1,#2,#3[,#4]
c$$$Ci         :
c$$$Ci         :l1 and l2 specify the first and second direction vectors,
c$$$Ci         :respectively.  #1,#2,#3 select the
c$$$Ci         :increments in mesh points along each of the three lattice
c$$$Ci         :vectors that define the direction vector.  Thus in Cartesian
c$$$Ci         :coordinates a direction vector is
c$$$Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
c$$$Ci         :where as before (n1,n2,n3) are the number of divisions in
c$$$Ci         :the mesh along the three lattice vectors (p1,p2,p3).
c$$$Ci         :The last number (#4) specifies how many points to take
c$$$Ci         :in that direction.
c$$$Ci         :
c$$$Ci         :Other options:
c$$$Ci         :  g3d       output of 3d grid plus headers in xsf format ready for xcrysden
c$$$Ci         :  fn=name   specifies the file name for file I/O
c$$$Ci         :
c$$$Ci         :   core=#   specifies how local rho is to be included
c$$$Ci         :            #=0 include core densities - nuclear charge
c$$$Ci         :            #=1 include core densities
c$$$Ci         :            #=2 (default) exclude core densities
c$$$Ci         :            #=-1 no local densities to be included
c$$$Ci         :            #=-2 true local density, no smoothed part
c$$$Ci         :            #=-3 istl-local sm densities, no true part
c$$$Ci         :
c$$$Ci         :Example: using '~' as the delimiter, and suppose
c$$$Ci         :n1=n2=48 and n3=120, the specification
c$$$Ci         :  ~fn=myrho~o=0,0,60~l1=1,1,0,49~l2=0,0,1,121
c$$$Ci         :writes 'myrho.ext' a mesh (49,121) points.
c$$$Ci         :The origin (first point) lies at (p3/2)
c$$$Ci         :The first vector points along (p1+p2), and has that length;
c$$$Ci         :the second vector points along p3, and has that length.
c$$$Ci   slat  :struct for lattice information; see routine ulat
c$$$Ci     Elts read: plat nabc ng okv vol
c$$$Ci     Stored:    *
c$$$Ci     Passed to: rhgcmp
c$$$Ci   ssite :struct for site-specific information; see routine usite
c$$$Ci     Elts read: *
c$$$Ci     Stored:    *
c$$$Ci     Passed to: rhgcmp
c$$$Ci   sspec :struct for species-specific information; see routine uspec
c$$$Ci     Elts read: *
c$$$Ci     Stored:    *
c$$$Ci     Passed to: rhgcmp
c$$$Ci   orhat :offsets to arrays containing local atom densities
c$$$Ci   smrho :smooth density on uniform mesh
c$$$Co Outputs
c$$$Co   The sum local gaussian densities and smrho is written to disk
c$$$Cl Local variables
c$$$Cl   modrhg:controls what part of core density is added
c$$$Cl         :-1 Add no local densities
c$$$Cl         : 0 include core density - sm nuclear charge
c$$$Cl         : 1 include core density, but no nuclear charge
c$$$Cl         : 2 exclude core density
c$$$Cr Remarks
c$$$Cr    sopts specifies which plane(s) are written to disk
c$$$Cr Bugs
c$$$Cr   Routines create smoothe approximation to density, not true density
c$$$Cu Updates
c$$$Cu   July 7 05 (Walter Lambrecht) New option g3d
c$$$Cu   25 Aug 04 New modes -2, -3
c$$$Cu   24 May 03 Corrections to errors in messages
c$$$Cu   23 Oct 01 Local densities are more accurately represented
c$$$Cu             in G_kL expansion:  k=0..kmax.  Added core= option.
c$$$Cu   25 Apr 01 Simplified the plane specification
c$$$Cu   02 Mar 01 Spin polarized
c$$$Cu   09 Feb 01 Added local gaussian densities to file I/O
c$$$Cu             and extended how a plane may be specified
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character*(*) sopts
c$$$      type(s_lat)::slat
c$$$      type(s_site)::ssite(*)
c$$$      type(s_spec)::sspec(*)
c$$$
c$$$      type(s_rv1) :: sv_p_orhat
c$$$
c$$$      double complex smrho(*)
c$$$C ... Local parameters
c$$$      integer ngabc(3),n1,n2,n3,k1,k2,k3,kmax
c$$$      integer kkk,lgunit,ng,nglob,nsp,nwk,stdo,modrhg,nbas,i
c$$$ckino Dec.9.2011:          real(8),pointer :: rv_p_opos(:) =>NULL()
c$$$
c$$$ckino Dec.9.2011:         integer,pointer :: iv_p_okv(:) =>NULL()
c$$$
c$$$      complex(8) ,allocatable :: cn_zv(:)
c$$$      complex(8) ,allocatable :: psrho_zv(:)
c$$$ckino not complex but real      complex(8) ,allocatable :: wk_zv(:)
c$$$      real(8) ,allocatable :: wk_rv(:)
c$$$
c$$$      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
c$$$      double precision vol,xx0,xx,dval,plat(3,3),alat
c$$$      integer ib,is,igetss
c$$$C ... Heap
c$$$      integer ::iwdummy
c$$$
c$$$      logical:: l_dummy_isanrg,isanrg
c$$$
c$$$C ... External calls
c$$$ckino del upack,upack2
c$$$      external awrit2,dcopy,fclose,fftz3,fftz30,gvgetf,
c$$$     .gvputf,icopy,ivset,mkils0,mkilst,nwordg,poppr,
c$$$     .pshpr,rhgcmp,rhomom,rx
c$$$
c$$$      real(8),allocatable:: zz(:)
c$$$C ... Unpack and setup
c$$$Changenglob      nsp = nglob('nsp')
c$$$      nsp = globalvariables%nsp
c$$$Changenglob      nbas = nglob('nbas')
c$$$      nbas = globalvariables%nbas
c$$$      stdo = lgunit(1)
c$$$
c$$$      i_copy_size=size(slat%plat)
c$$$      call dcopy(i_copy_size,slat%plat,1,plat,1)
c$$$      i_copy_size=size(slat%nabc)
c$$$      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
c$$$      ng=slat%ng
c$$$ckino Dec.9.2011:         iv_p_okv => slat%iv_p_okv
c$$$
c$$$      vol=slat%vol
c$$$
c$$$ccccccccccccccccccccccccccccc
c$$$c      allocate(zz(nbas))
c$$$c      zz=sspec(1:nbas)%z
c$$$c      print *,' zz1=',zz
c$$$c      print *,' zz2=',sspec(1:nbas)%z
c$$$c      stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzz'
c$$$ccccccccccccccccccccccccccccc
c$$$
c$$$
c$$$      alat=slat%alat
c$$$ckino Dec.9.2011:         rv_p_opos => slat%rv_p_opos
c$$$
c$$$
c$$$c      call defrr(oz,nbas)
c$$$c      do ib = 1,nbas
c$$$c        is = int(ssite(ib)%spec)
c$$$cc        w(oz-2+ib*2)=sspec(is)%z !but Kino's email on 5June2010.
c$$$c        call dcopy(1,sspec(is)%z,1,w(oz-2+ib*2),1)
c$$$c      enddo
c$$$      call fftz30(n1,n2,n3,k1,k2,k3)
c$$$      kkk = k1*k2*k3
c$$$      kmax = 3
c$$$      modrhg = 2
c$$$c      call ioden2(0,sopts,nsp,plat,alat,nbas,w(opos),w(oz),k1,k2,k3,
c$$$c     .ngabc,w,w,nwk,modrhg)
c$$$Ckino  mode=0 checks format
c$$$#ifdef USEOPTIONAL
c$$$      call ioden2 ( mode=0 , sopts=sopts ,nsp= nsp ,plat= plat ,
c$$$ckino Dec.9.2011:       .   alat=alat , nbas=nbas , pos=rv_p_opos
c$$$ckino Jan.04.2012:       .   alat=alat , nbas=nbas , pos=slat%rv_p_opos 
c$$$     .   alat=alat , nbas=nbas , pos=slat%rv_a_opos 
c$$$     .  , z=sspec ( ssite ( 1:nbas ) %spec ) %z , k1=k1 , k2=k2 , k3=k3
c$$$     .   , ngabc=ngabc 
c$$$     .  , nwk=nwk , modrhg=modrhg )
c$$$#else
c$$$ckino Dec.9.2011:          call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , rv_p_opos
c$$$ckino Jan.04.2012:          call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , slat%rv_p_opos
c$$$      call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , slat%rv_a_opos 
c$$$     .  , sspec ( ssite ( 1:nbas ) %spec ) %z , k1 , k2 , k3 , ngabc 
c$$$     .  , iwdummy , iwdummy , nwk , modrhg )
c$$$#endif
c$$$
c$$$
c$$$Ckino isanrg is logical function,       call isanrg(modrhg,-3,2,'ioden','core option',.true.)
c$$$      l_dummy_isanrg=isanrg(modrhg,-3,2,'ioden','core option',.true.)
c$$$
c$$$C ... Overwrite smrho+, smrho- with smrho, smrho+ - smrho-
c$$$Ckino  smrho is complex(8), dsumdf() is for real(8)
c$$$Ckino  so kkk*2 is necessary to use dsumdf().
c$$$Ckino
c$$$Ckino  return values: smrho(:kkk)=spin(1)+spin(2)
c$$$Ckino  and            smrho(kkk+1:kkk*2)=spin(1)-spin(2)
c$$$      if (nsp .eq. 2) then
c$$$        call dsumdf(kkk*2,1d0,smrho,0,1,smrho(1+kkk),0,1)
c$$$      endif
c$$$
c$$$C ... Put n0(G) into psrho and cn
c$$$      allocate(psrho_zv(kkk*nsp))
c$$$
c$$$      allocate(cn_zv(ng*nsp))
c$$$
c$$$Ckino smrho is complex(8) , so kkk*2 is necessary.
c$$$      call dcopy ( kkk * 2 * nsp , smrho , 1 , psrho_zv , 1 )
c$$$
c$$$      call fftz3 ( psrho_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp 
c$$$     ., 0 , - 1 )
c$$$
c$$$ckino Dec.9.2011:         call gvgetf ( ng , nsp , iv_p_okv , k1 , k2 , k3 , psrho_zv ,
c$$$ckino Dec.28.2011:         call gvgetf ( ng , nsp , slat%iv_p_okv , k1 , k2 , k3 , psrho_zv ,
c$$$      call gvgetf ( ng , nsp , slat%iv_a_okv , k1 , k2 , k3 , psrho_zv , 
c$$$     . cn_zv )
c$$$Ckino now cn_zv(G,nsp)
c$$$
c$$$
c$$$Ckino xx0 = real and G=0 component
c$$$      xx0 = dval ( cn_zv , 1 )
c$$$
c$$$
c$$$C ... Add sum of local gaussian densities to mesh density
c$$$      print *,' modrhg=',modrhg
c$$$Ckino core=#   specifies how local rho is to be included
c$$$Ckino         :            #=0 include core densities - nuclear charge
c$$$Ckino         :            #=1 include core densities
c$$$Ckino         :            #=2 (default) exclude core densities
c$$$Ckino         :            #=-1 no local densities to be included
c$$$Ckino         :            #=-2 true local density, no smoothed part
c$$$Ckino         :            #=-3 istl-local sm densities, no true part
c$$$      if (modrhg .ge. 0 .or. modrhg .le. -2) then
c$$$        if (modrhg .eq.  0) i = 131
c$$$        if (modrhg .eq.  1) i =  31
c$$$        if (modrhg .eq.  2) i =   1
c$$$        if (modrhg .eq. -2) i =   2
c$$$        if (modrhg .eq. -3) i =   3
c$$$        if ( i .eq. 2 ) call dscal ( ng * nsp * 2 , 0d0 , cn_zv , 
c$$$     .  1 )
c$$$Ckino                      -> cn_zv=0 if i==2
c$$$
c$$$        if ( i .eq. 3 ) call dscal ( ng * nsp * 2 , - 1d0 , cn_zv 
c$$$     .  , 1 )
c$$$Ckino                      -> cn_zv=-cn_zv if i==3
c$$$
c$$$        call rhgcmp ( i , 1 , nbas , ssite , sspec , slat , sv_p_orhat 
c$$$     .   , kmax , ng , cn_zv )
c$$$
c$$$
c$$$        if ( i .eq. 3 ) call dscal ( ng * nsp * 2 , - 1d0 , cn_zv 
c$$$     .  , 1 )
c$$$
c$$$      endif
c$$$      xx = dval ( cn_zv , 1 )
c$$$
c$$$
c$$$C ... FFT (n0 + gaussians) (G) to real-space mesh
c$$$      if (modrhg .ge. 0) then
c$$$ckino Dec.9.2011:           call gvputf ( ng , nsp , iv_p_okv , k1 , k2 , k3 , cn_zv , psrho_zv
c$$$ckino Dec.28.2011:           call gvputf ( ng , nsp , slat%iv_p_okv , k1 , k2 , k3 , cn_zv , psrho_zv
c$$$        call gvputf ( ng , nsp , slat%iv_a_okv , k1 , k2 , k3 , cn_zv , psrho_zv 
c$$$     .   )
c$$$
c$$$
c$$$        call awrit2('%N ioden : local densities + envelope density,'
c$$$     .  //' Qloc=%,6;6d  Q=%,6;6d',' ',80,stdo,(xx-xx0)*vol,xx*vol)
c$$$        if (modrhg .lt. 2) then
c$$$          call awrit1('%9fLocal densities include core'//
c$$$     .    '%?#n==0#+nuclear## contributions.',' ',80,stdo,modrhg)
c$$$        endif
c$$$      elseif (modrhg .eq. -2) then
c$$$        call info2(0,0,0,'%N ioden : local'//
c$$$     .  ' densities (true-smooth terms),  Qloc=%,6;6d',xx*vol,0)
c$$$      elseif (modrhg .eq. -3) then
c$$$        call info2(0,0,0,'%N ioden : smooth (envelope - local)'//
c$$$     .  ' density,  Qs=%,6;6d  Qs-Qs(loc)=%,6;6d',
c$$$     .  xx0*vol,xx*vol)
c$$$      else
c$$$        call info2(0,0,0,'%N ioden : smooth density only'//
c$$$     .  ' (no local densities added) Qsm=%,6;6d',xx0*vol,0)
c$$$      endif
c$$$
c$$$      call fftz3 ( psrho_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp 
c$$$     ., 0 , 1 )
c$$$Ckino now psrho_zv is real space mesh.
c$$$
c$$$
c$$$      if (allocated(cn_zv)) deallocate(cn_zv)
c$$$
c$$$
c$$$C ... File I/O
c$$$      nwk = 12*max(k1,k2,k3)
c$$$
c$$$ckino not complex but real,      allocate(wk_zv(nwk**2*nsp))
c$$$      allocate(wk_rv(nwk**2*nsp))
c$$$
c$$$Ckino mode=1 write smrho to the file
c$$$ckino Dec.9.2011:         call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , rv_p_opos
c$$$ckino Jan.04.2012:         call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , slat%rv_p_opos
c$$$      call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , slat%rv_a_opos 
c$$$     . , sspec ( ssite ( 1:nbas ) %spec ) %z , k1 , k2 , k3 , ngabc
c$$$ckino not complex but real     . , psrho_zv , wk_zv , nwk , modrhg ) 
c$$$     . , psrho_zv , wk_rv , nwk , modrhg )
c$$$
c$$$
c$$$ckino not complex but real      if (allocated(wk_zv)) deallocate(wk_zv)
c$$$      if (allocated(wk_rv)) deallocate(wk_rv)
c$$$      if (allocated(psrho_zv)) deallocate(psrho_zv)
c$$$
c$$$c      call rlse(oz)
c$$$
c$$$
c$$$C ... Restore smrho+, smrho-
c$$$      if (nsp .eq. 2) then
c$$$        call dsumdf(kkk*2,.5d0,smrho,0,1,smrho(1+kkk),0,1)
c$$$      endif
c$$$      end subroutine ioden
c$$$
c$$$
c$$$
c$$$
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
c$$$  #define MULLIKEN 1
c$$$      subroutine suqlst(sopts,iop,nband,efermi,nsp,evl,nfbn,ifblst,nq,
c$$$     .qp, onespdummy) !onesp removed 
c$$$      use m_globalvariables
c$$$C- Set up a list of q-points in various modes for energy bands
c$$$C  See also entry suqlsw(nband,jsp,nsp,evl)
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   sopts :character string describing special options that can
c$$$Ci         :selects the mode for which qp are read.
c$$$Ci         :*default mode: the qp file consists of a list of
c$$$Ci         :               lines and the number of qp in each line.
c$$$Ci         :*list mode   : the qp file consists of a list of qp
c$$$Ci         :*contour mode: input for specifying qp on a uniform 2D
c$$$Ci                       : mesh, for contour plots.
c$$$Ci         :See Remarks for the syntax of each mode.
c$$$Ci
c$$$Ci         :Options are separated by delimiters; the first character is
c$$$Ci         :the delimiter.  The following list the strings declaring
c$$$Ci         :options assuming the delimiter is '/'. (Space as the first
c$$$Ci         :implies that there are no options)
c$$$Ci         :/fn=fnam     read qp from file 'fnam'
c$$$Ci         :/con         input file for contour plot mode (mode 3)
c$$$Ci         :/qp          input file specifies a list of qp (mode 2)
c$$$Ci         :             (default assumes a list of symmetry lines.)
c$$$Ci         :/long        write bands with extra digits precision
c$$$Ci         :             (has no effect for symmetry line mode)
c$$$Ci         :/spin1       generate bands only for first spin
c$$$Ci         :/lst=list    write only those bands specified in a list.
c$$$Ci         :             For syntax of list, see slatsm/mkilst.f
c$$$Ci         :/evn=#       keep track of smallest, largest eval for
c$$$Ci                       #th band, and printout at close.
c$$$Ci         :/ef=#        change efermi to #.
c$$$Ci         :Example: --band~long~qp~lst=2:5
c$$$Ci   iop   :options passed by the caller
c$$$Ci         :1s digit is suppress actions for parallel mode
c$$$Ci         : mode 1:
c$$$Ci         : 1: suppress writing line header info file
c$$$Ci         : 2: Return total number of qp to be generated in all lines
c$$$Ci         :    and also suppress writing line header info file
c$$$Ci         :    Thus iop=2 => all modes return total no qp to be generated
c$$$Ci         : mode 3:
c$$$Ci         : nonzero: do not allocate evsav until iq>nq
c$$$Ci
c$$$Ci   nband :(suqlst) maximum number of energy bands to write
c$$$Ci         :(suqlsw) actual number of energy bands to write
c$$$Ci
c$$$Ci   efermi:Fermi energy (written to bnds file)
c$$$Ci
c$$$Ci   nsp   :2 for spin-polarized case, otherwise 1
c$$$Ci         :NB used only in file write to indicate how many times
c$$$Ci         :nsp is ALTERED to nsp=1 if spin1 option is set
c$$$Ci
c$$$Ci   evl   :eigenvalues for current qp (used only in suqlsw)
c$$$Cio Inputs/Outputs (see Remarks)
c$$$Cio  nq    :Input nq=0 :if zero, flags suqlst to set up bands mode:
c$$$Cio                    :some initializations; then returns ready
c$$$Cio                    :for first block of qp
c$$$Cio        :           :if nonzero, setup for for next qp
c$$$Cio        :Output nq  :if zero, suqlst has no more qp to calculate.
c$$$Cio                    :if nonzero, nq is the number of k-points in the
c$$$Cio                    :current block.  Caller then calls suqlst nq
c$$$Cio                    :times (one for each k-point), which will return
c$$$Cio                    :qp for each k-point.  After exactly nq calls,
c$$$Cio                    :suqlst will start another block, if it exists.
c$$$Cio                    :See Remarks for schematic of calling sequence.
c$$$Co Outputs
c$$$Co   nfbn  :(color weights) number of elements in iblst
c$$$Co   iblst :(color weights) list of orbital indices for color weights:
c$$$Co         :color weight is the sum of contributions from eigenvector
c$$$Co         :components in iblst
c$$$Co   qp    :not set, if suqlst is called starting a new block.
c$$$Co         :(see Remarks).  Else, k-point at which to generate bands
c$$$Co    onesp(removed by takao) :
c$$$Co    if spin1 flag is encountered, onesp is set to 1 and nsp is set to 1
c$$$Cl Local variables
c$$$Cl   iq    :current qp in this block
c$$$Cl   mode  :1 symmetry-line mode
c$$$Cl         :2 list-of-qp mode
c$$$Cl         :3 contour mode
c$$$Cl   q1    :starting qp for symmetry mode; only meaningful in that mode
c$$$Cl   q2    :ending   qp for symmetry mode; only meaningful in that mode
c$$$Cl   nevn  :(optional) band index; routine monitors largest, smallest value
c$$$Cl         :for that index.
c$$$Cl   evmnn :smallest value found for qp of specified band index
c$$$Cl   evmxn :largest  value found for qp of specified band index
c$$$Cl   ifiq  :file logical unit for input qp file
c$$$Cl   ifib  :file logical unit for output bands file
c$$$Cl   ql    :local copy of current qp
c$$$Cl   oqp
c$$$Cr Remarks
c$$$Cr   suqlst is designed to be called to generate qp in groups or
c$$$Cr   blocks.  The calling sequence is:
c$$$Cr
c$$$Cr     nq = 0  <- flags that first call to suqlst, to set up mode
c$$$Cr     do  iblock = 1, forever
c$$$Cr       This call generates nq, the number of points in this block
c$$$Cr       call suqlst(nband,efermi,nsp,evl,nq,qp,onesp)
c$$$Cr       if (nq .eq. 0) stop
c$$$Cr       do  iq = 1, nq
c$$$Cr         This call generates qp for current block
c$$$Cr         call suqlst(ndimh,ef0,nsp,w,nkp,qp,onesp) <- returns qp
c$$$Cr         do  isp = 1, nsp
c$$$Cr         call suqlsw(ndimh,qp,evl(1,isp)) <- saves evl for this qp
c$$$Cr                                             (call is optional)
c$$$Cr         enddo
c$$$Cr       enddo
c$$$Cr     enddo
c$$$Cr
c$$$Cr   The following modes are implemented:
c$$$Cr     mode=1 reads qp from syml file, and generates qp along each
c$$$Cr            specified symmetry line.  Structure of qp file:
c$$$Cr            file has one line for each symmetry line as follows:
c$$$Cr               nq      q1x   q1y   q1z      q2x   q2y   q2z
c$$$Cr               ...
c$$$Cr            the line entries have meanings:
c$$$Cr            --# qp-   ---starting qp---    --- ending qp ---
c$$$Cr            Any line with nq=0 implies no more lines.
c$$$Cr     mode=2 reads qp from specified file and generates qp for each
c$$$Cr            specified qp.  File consists sets of qpx,qpy,qpz for each
c$$$Cr            qp sought.  Typically a file would consists of lines like
c$$$Cr               q1x   q1y   q1z
c$$$Cr               q2x   q2y   q2z
c$$$Cr               ...
c$$$Cr     mode=3 generates qp for a uniform mesh in a plane (contour plot)
c$$$Cr            The file supplies information describing a rectangle in
c$$$Cr            the Brillouin zone.  It consists of a single line,
c$$$Cr            which contains the following:
c$$$Cr             v1    range  n     v2    range  n   height  list-of-bands
c$$$Cr
c$$$Cr            v1 and v2 are two vectors specifying the plane of the
c$$$Cr            contour.  range and n (one each for v1 and v2) are the
c$$$Cr            starting and final amplitudes of those vectors, and the
c$$$Cr            the number of points within the vector.  list-of-bands
c$$$Cr            is a list of integers which specify which bands are to
c$$$Cr            be written to the output file.  'height' is the 'z' axis.
c$$$Cr            For example,
c$$$Cr             v1    range  n     v2    range  n   height  list-of-bands
c$$$Cr            1 0 0  -1 1   51   0 1 0  -1 1   51   0.00    4,5
c$$$Cr            creates a file of 51x51 points, with the four corners
c$$$Cr            (-1,-1,0),  (1,-1,0),  (-1,1,0),  (1,1,0)
c$$$Cr            going through the gamma-point. Two bands (4,5) are stored.
c$$$Cu Updates
c$$$Cu   08 Jul 08 Extend to case where number of bands can be q dependent
c$$$Cu             modes 1,2: suqlsw writes out number of bands with qp
c$$$Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1 
c$$$c                                                          -->forbidden 2014takao
c$$$Cu   02 Jul 06 Color mode extended to two colors
c$$$Cu   02 Jul 06 New color mode (one color weight only)
c$$$Cu   14 Feb 05 contour mode saves both spins in spin-polarized case
c$$$Cu   20 Oct 03 suqlst works properly in contour mode
c$$$Cu   28 Aug 01 added ef switch
c$$$Cu   23 Jan 01 first written
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character*(*) sopts
c$$$      integer nq,nband,nsp,nfbn(2),ifblst(nband,*),onespdummy,iop!onesp
c$$$      double precision efermi,qp(3),evl(nband)
c$$$C ... Local variables
c$$$      character strn*120, strn2*120, dc*1, fn*120, prfmt*40
c$$$      logical rdstrn,fmlong
c$$$      integer i,iprint,j,j1,j2,k,nglob,rdm,
c$$$     .parg,a2vec,iv(7),stdo,jsp,op1
c$$$      double precision xx,xx1,xx2,sclp,plndst,xxv(7)
c$$$      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
c$$$      integer oqp,nqall,imode
c$$$      double precision x1,x2,y1,y2,evmxn,evmnn,
c$$$     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)
c$$$      real(8),allocatable :: evsav(:,:,:,:)
c$$$C ... MPI
c$$$      integer procid,master,mpipid
c$$$
c$$$      logical:: l_dummy_isanrg,isanrg
c$$$
c$$$ctakao
c$$$      real(8)::qlx(3)
c$$$C ... Heap
c$$$c      integer w(1)
c$$$c      common /w/ w
c$$$      common /suqlsd/
c$$$     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
c$$$     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c$$$c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy
c$$$      save evsav
c$$$      real(8),allocatable:: w_oqp(:,:)
c$$$
c$$$      procid = mpipid(1)
c$$$      master = 0
c$$$Changenglob      stdo = nglob('stdo')
c$$$      stdo = globalvariables%stdo
c$$$      op1 = mod(iop,10)
c$$$      nqall = 0
c$$$
c$$$C --- First call ... setup and parse options ---
c$$$      if (nq .eq. 0) then
c$$$C   ... Defaults
c$$$        fmlong = .false.
c$$$        mode = 1
c$$$        fn = 'syml'
c$$$        nblst = 0
c$$$        nevn = 0
c$$$        evmxn = -99d9
c$$$        evmnn = 99d9
c$$$        nfbn(1) = 0
c$$$        nfbn(2) = 0
c$$$
c$$$        dc = sopts(1:1)
c$$$        if (dc .ne. ' ') then
c$$$C   ... Return here to resume parsing for arguments
c$$$          j2 = 0
c$$$   10     continue
c$$$          j2 = j2+1
c$$$          if (sopts(j2:j2) .eq. dc) goto 10
c$$$          j1 = min(len(sopts),j2)
c$$$          call nwordg(sopts,0,dc//' ',1,j1,j2)
c$$$          if (j2 .ge. j1) then
c$$$            if (.false.) then
c$$$            elseif (sopts(j1:j1+2) .eq. 'fn=')  then
c$$$              if (j1+3 .gt. j2) call rx('suqlst: bad file name')
c$$$              fn = sopts(j1+3:j2)
c$$$            elseif (sopts(j1:j2) .eq. 'qp')  then
c$$$              mode = 2
c$$$            elseif (sopts(j1:j2) .eq. 'spin1')  then
c$$$               call rx('suqlst:spin removed! takao')
c$$$c              onesp = 1
c$$$c              nsp = 1
c$$$            elseif (sopts(j1:j2) .eq. 'con')  then
c$$$              mode = 3
c$$$            elseif (sopts(j1:j2) .eq. 'long')  then
c$$$              fmlong = .true.
c$$$            elseif (sopts(j1:j1+3) .eq. 'col=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list, col=..')
c$$$              call mkils0(sopts(j1+4:j2),nfbn,ifblst)
c$$$              call mkilst(sopts(j1+4:j2),nfbn,ifblst)
c$$$              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
c$$$            elseif (sopts(j1:j1+4) .eq. 'col2=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list, col2=..')
c$$$              call mkils0(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
c$$$              call mkilst(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
c$$$              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
c$$$            elseif (sopts(j1:j1+3) .eq. 'lst=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list, lst=..')
c$$$              call mkils0(sopts(j1+4:j2),nblst,iblst)
c$$$              if (nblst .gt. 100) call rx('increase size of iblst')
c$$$              call mkilst(sopts(j1+4:j2),nblst,iblst)
c$$$            elseif (sopts(j1:j1+2) .eq. 'ef=')  then
c$$$              j = 0
c$$$              i = parg('ef=',4,sopts(j1:),j,len(sopts(j1:)),
c$$$     .        dc//' ',1,1,i,efermi)
c$$$              if (i .le. 0) call rxs(
c$$$     .        'suqlst: failed to parse string for fermi level:  ',
c$$$     .        sopts(j1:))
c$$$            elseif (sopts(j1:j1+3) .eq. 'evn=')  then
c$$$              if (j1+4 .gt. j2) call rx('suqlst: bad list')
c$$$              i = j1+3
c$$$              xxv(1) = 0
c$$$              j = a2vec(sopts,j2,i,2,dc//' ',2,3,1,iv,nevn)
c$$$              if (j .ne. 1 .or. nevn .gt. nband)
c$$$     .        call rx('suqlst: bad value for evn')
c$$$            else
c$$$              call rxs('suqlst: failed to parse argument, ',sopts(j1:j2))
c$$$            endif
c$$$            goto 10
c$$$          endif
c$$$        endif
c$$$
c$$$        if (nfbn(1) .eq. 0) then
c$$$          call info2(20,0,0,' suqlst:  generate bands, mode %i',mode,0)
c$$$        elseif (nfbn(2) .eq. 0) then
c$$$          call ilst2a(ifblst,nfbn,strn)
c$$$          call info2(20,0,0,
c$$$     .    ' suqlst:  generate bands with color weights, mode %i.'//
c$$$     .    ' %N%10f%i components: '//strn//'%a',mode,nfbn)
c$$$        else
c$$$          call ilst2a(ifblst,nfbn,strn)
c$$$          call info2(20,0,0,' suqlst:  generate bands '//
c$$$     .    'with two color weights, mode %i.',mode,0)
c$$$          call ilst2a(ifblst,nfbn,strn)
c$$$          call info2(20,0,0,
c$$$     .    '%10f%i components, first color:  '//strn//'%a',nfbn,0)
c$$$          call ilst2a(ifblst(1,2),nfbn(2),strn)
c$$$          call info2(20,0,0,
c$$$     .    '%10f%i components, second color: '//strn//'%a',nfbn(2),0)
c$$$        endif
c$$$        if (nfbn(1) .eq. 0 .and. nfbn(2) .ne. 0) call info0(20,0,0,
c$$$     .  ' suqlst:  (warning) col2 specified, but col not')
c$$$
c$$$C   ... Open qp file
c$$$        if (procid .eq. master) then
c$$$          ifiq = fopno(fn)
c$$$          rewind ifiq
c$$$
c$$$C   ... open bands file
c$$$          ifib = fopnn('BNDS')
c$$$          rewind ifib
c$$$
c$$$C   ... Write header
c$$$          if (mode .eq. 1) then
c$$$            i = nblst
c$$$            if (nblst .eq. 0) i = nband
c$$$            if (nfbn(1) .eq. 0) then
c$$$C           Use separate format statment to circumvent gfortran bug
c$$$              write(ifib,335) i,efermi,0
c$$$  335         format(i5,f10.5,i6)
c$$$            elseif (nfbn(2) .eq. 0) then
c$$$              call ilst2a(ifblst,nfbn,strn)
c$$$              call strip(strn,j1,j2)
c$$$              write(ifib,336) i,efermi,1,strn(j1:j2)
c$$$  336         format(i5,f10.5,i6:'  col= ',a:'  col2= ',a)
c$$$            else
c$$$              call ilst2a(ifblst,nfbn,strn)
c$$$              call strip(strn,j1,j2)
c$$$              call ilst2a(ifblst(1,2),nfbn(2),strn2)
c$$$              call strip(strn2,j,k)
c$$$              write(ifib,336) i,efermi,2,strn(j1:j2),strn2(j:k)
c$$$            endif
c$$$          endif
c$$$        endif
c$$$
c$$$C   ... Other initializations
c$$$        iq = 0
c$$$        nq = -1
c$$$      endif
c$$$
c$$$C --- Setup for a new block of k-points, depending on mode ---
c$$$C     This branch occurs on completion of the last qp of the current block
c$$$C     which is marked by iq>nq
c$$$C     At the completion of this block: nq must be computed and:
c$$$C     (mode=1) q1,q2 set up.
c$$$C              Note: if 1s digit of iop is set in this mode,
c$$$C              this branch returns sum of all qp in all lines.
c$$$C              No setup for q1,q2; no
c$$$C     (mode=2) oqp allocated and loaded
c$$$C     (mode=3) nblst = number of bands to save
c$$$      if (iq .gt. nq) then
c$$$        iq = 1
c$$$C       Bands along specified symmetry lines
c$$$        if (mode .eq. 1) then
c$$$  725     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 999
c$$$          if (strn(1:1) .eq. '#') goto 725
c$$$C         Exit if first entry in line is zero
c$$$          i = 0
c$$$          xxv(1) = 0
c$$$          i = a2vec(strn,len(strn),i,4,', ',2,3,1,iv,xxv)
c$$$          if (i .eq. 1 .and. xxv(1) .eq. 0) goto 999
c$$$          i = 0
c$$$          i = a2vec(strn,len(strn),i,4,', ',2,3,7,iv,xxv)
c$$$          if (i .ne. 7 .and. iprint().ge.10) then
c$$$            write(stdo,
c$$$     .      '(/'' suqlst (warning) skipping line:''/''  '',a)') strn
c$$$            goto 725
c$$$          endif
c$$$          nq = xxv(1)
c$$$C         1 qp is nonsensical for a line
c$$$          if (nq .eq. 1) nq = 0
c$$$C         No qp: exit
c$$$          if (nq .le. 0) goto 999
c$$$C         setup q1,q2
c$$$          call dcopy(3,xxv(2),1,q1,1)
c$$$          call dcopy(3,xxv(5),1,q2,1)
c$$$          if (iprint().ge.10) write(stdo,785) nq,q1,q2
c$$$  785     format(/' suqlst:  nq=',i3,'   q1=',3f7.4,'   q2=',3f7.4)
c$$$C         Write line header information to disk
c$$$          if (op1 .eq. 0) then
c$$$            write(ifib,337) nq*nsp
c$$$  337       format(2i5)
c$$$C         Accumulate all qp and cycle until all lines are read
c$$$          elseif (op1 .eq. 2) then
c$$$            nqall = nqall + nq
c$$$            goto 725
c$$$          endif
c$$$C       Bands for a list of specified qp
c$$$        elseif (mode .eq. 2) then
c$$$C         Only one block for this mode.
c$$$C         Flag that prior block already completed: nq>0
c$$$          if (nq .gt. 0) goto 999
c$$$          nq = 0
c$$$          j = 3
c$$$          if (rdm(ifiq,10000,0,' ',xx,j,nq) .ne. 1) call
c$$$     .    rx('suqlst:  bad qp file')
c$$$
c$$$c          call defrr(oqp,3*nq)
c$$$          allocate(w_oqp(3,nq))
c$$$
c$$$          rewind ifiq
c$$$          j = rdm(ifiq,10000,3*nq,' ',w_oqp,j,nq)
c$$$          call awrit1('%N suqlst: read %i qp from file '//fn//
c$$$     .    '%a',' ',80,stdo,nq)
c$$$          if (j .ne. 1) call rx('suqlst: failed to read qp')
c$$$          i = nblst
c$$$          if (nblst .eq. 0) i = nband
c$$$          if (nfbn(1) .ne. 0) i = i*2
c$$$          call awrit2('%% rows %i cols %i',' ',80,ifib,nq,i+4)
c$$$C       Bands on a uniform mesh in a specified plane (contour plot)
c$$$        elseif (mode .eq. 3 .and. nq .eq. -1) then
c$$$  825     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 998
c$$$          if (strn(1:1) .eq. '#') goto 825
c$$$          call words(strn,i)
c$$$          if (i .ne. 14) call rxi('suqlst con mode: expected 14 '//
c$$$     .    'arguments from input file but read',i)
c$$$          iblst(1) = -1
c$$$          nblst = 1
c$$$          backspace ifiq
c$$$          read(ifiq,*) vectx,x1,x2,nqx,vecty,y1,y2,nqy,plndst
c$$$          call word(strn,14,j1,j2)
c$$$          call mkilss(11,strn(j1:j2),nblst,iblst)
c$$$          if (nblst .le. 0) call rx('suqlst: no bands in list')
c$$$          call dscal(3,1/dsqrt(sclp(vectx,vectx)),vectx,1)
c$$$          call dscal(3,1/dsqrt(sclp(vecty,vecty)),vecty,1)
c$$$C         Subtract from vecty projection onto vectx
c$$$C         call daxpy(3,-sclp(vectx,vecty),vectx,1,vecty,1)
c$$$          call cross(vectx,vecty,vect3)
c$$$          call dscal(3,plndst/dsqrt(sclp(vect3,vect3)),vect3,1)
c$$$          nq = nqx*nqy
c$$$          if (iprint() .ge. 10) then
c$$$            write(stdo,717) vectx,x1,x2,nqx,plndst,vecty,y1,y2,nqy,nq
c$$$  717       format(' vx=',3f9.6,'  x1,x2=',2f9.6,'  nx=',i3,'   h=',f9.6/
c$$$     .      ' vy=',3f9.6,'  y1,y2=',2f9.6,'  ny=',i3,'  np=',i5)
c$$$            write(strn,'('' save %i bands: %'',i2,'':1i'')') nblst
c$$$            if (strn(17:18) .eq. '% ') strn(17:18) = ' %'
c$$$            call awrit2(strn,strn,80,stdo,nblst,iblst)
c$$$            if (dabs(sclp(vectx,vecty)) .gt. 1d-6 .and. iprint().gt.20)
c$$$     .      print *,'suqlst (warning): input plane vectors not orthogonal'
c$$$          endif
c$$$
c$$$          if (iop .eq. 0) then
c$$$            allocate(evsav(nqx,nqy,nblst,nsp))
c$$$          endif
c$$$
c$$$C       Contour plot, cleanup.
c$$$C       Note: cleanup handled by suqlsw when last qp is called.
c$$$        elseif (mode .eq. 3) then
c$$$C          if (op1 .ne. 0) then
c$$$C            allocate(evsav(nqx,nqy,nblst,nsp))
c$$$C            return
c$$$C          endif
c$$$          call rx('suqlst: caller should never reach this branch')
c$$$        endif
c$$$
c$$$C --- Generate qp for this iq, depending on mode ---
c$$$      else
c$$$        if (mode .eq. 1) then
c$$$          xx = dble(iq-1)/dble(nq-1)
c$$$          qp(1) = xx*q2(1) + (1-xx)*q1(1)
c$$$          qp(2) = xx*q2(2) + (1-xx)*q1(2)
c$$$          qp(3) = xx*q2(3) + (1-xx)*q1(3)
c$$$        elseif (mode .eq. 2) then
c$$$          call dpscop(w_oqp,qp,3,iq*3-2,1,1d0)
c$$$        elseif (mode .eq. 3) then
c$$$C         Inner Loop:  excursions in y; outer loop: excursions in x
c$$$          j = mod(iq-1,nqy)
c$$$          i = (iq-1-j)/nqy
c$$$          if (nqx .le. 1) then
c$$$            xx1 = x1
c$$$          else
c$$$            xx1 =i*(x2-x1)/(nqx-1)+x1
c$$$          endif
c$$$          if (nqy .le. 1) then
c$$$            xx2 = y1
c$$$          else
c$$$            xx2 =j*(y2-y1)/(nqy-1)+y1
c$$$          endif
c$$$          do  93  k = 1, 3
c$$$            qp(k) = xx1*vectx(k) + xx2*vecty(k) + vect3(k)
c$$$   93     continue
c$$$          if (j .eq. 0 .and. iprint().ge.20) write(stdo,718) i+1,nqx,qp
c$$$  718     format(' line',i3,' of',i3,'   q(1)=',3f10.6)
c$$$        else
c$$$          call rx('suqlst: bad mode')
c$$$        endif
c$$$        iq = iq+1
c$$$C       Hold onto local copy of qp
c$$$        call dcopy(3,qp,1,ql,1)
c$$$      endif
c$$$      if(allocated(w_oqp)) deallocate(w_oqp)
c$$$      return
c$$$
c$$$C --- No more qp blocks: cleanup ---
c$$$  999 continue
c$$$      nq = 0
c$$$      if (nevn .ne. 0) call awrit3(' eval no. %i:  minimum eval'//
c$$$     .' = %;8F  maximum eval = %;8F',' ',80,stdo,nevn,evmnn,evmxn)
c$$$
c$$$      if (mode .eq. 1) then
c$$$        if (op1 .eq. 0) then
c$$$          write(ifib,337) 0
c$$$        elseif (op1 .eq. 2) then
c$$$          nq = nqall
c$$$        endif
c$$$      elseif (mode .eq. 3) then
c$$$        call rx('not ready')
c$$$C       call xxxbnd(w(oev),nblst,nqx,nqy,ifib)
c$$$C       return
c$$$      endif
c$$$      if(allocated(w_oqp)) deallocate(w_oqp)
c$$$      return
c$$$
c$$$C --- Error exit ---
c$$$  998 call rxs('suqlst: failed to read file contents, file ',fn)
c$$$
c$$$      entry suqlsm(imode)
c$$$C- Return qlist mode
c$$$      imode = mode
c$$$      return
c$$$
c$$$      entry suqlsw(nband,jsp,nsp,evl)
c$$$C- Write or store the energy bands to file for this qp
c$$$
c$$$      if (mode .eq. 1) then
c$$$        if (nblst .eq. 0) then
c$$$          prfmt = '(3f10.5,i6/(10f8.4))'
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.5/(10f8.4))'
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 2) then
c$$$        if (fmlong .and. nblst .eq. 0) then
c$$$          prfmt = '(3f15.10,i6/(5f15.10))'
c$$$        elseif (nblst .eq. 0) then
c$$$          prfmt = '(3f10.6,i6/(8f10.6))'
c$$$        elseif (fmlong .and. nblst .gt. 0) then
c$$$          prfmt = '(3f15.10/(5f15.10))'
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.6/(8f10.6))'
c$$$C         if (nblst .le. 5) prfmt = '(8f10.6)'
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 3) then
c$$$        j = mod(iq-2,nqy)
c$$$        i = (iq-2-j)/nqy
c$$$Changenglob        stdo = nglob('stdo')
c$$$        stdo = globalvariables%stdo
c$$$        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c$$$  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
c$$$
c$$$        do  k = 1, nblst
c$$$          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
c$$$        enddo
c$$$
c$$$C   ... If last qp generated, dump to file and exit
c$$$        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
c$$$          if (fmlong) then
c$$$            prfmt = '(5f15.10/(5f15.10))'
c$$$          else
c$$$            prfmt = '(8f10.6/(8f10.6))'
c$$$          endif
c$$$          rewind ifib
c$$$          do  j1 = 1, nsp
c$$$            do  k  = 1, nblst
c$$$              if (nsp .eq. 1)
c$$$     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
c$$$              if (nsp .eq. 2)
c$$$     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
c$$$     .        ifib,nqx,nqy,j1)
c$$$
c$$$              do  i = 1, nqx
c$$$                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          call rx0('finished generating bands on q-mesh')
c$$$        endif
c$$$      endif
c$$$
c$$$C ... Keep running tab on smallest,largest eval
c$$$      if (nevn .ne. 0) then
c$$$        evmxn = max(evmxn,evl(nevn))
c$$$        evmnn = min(evmnn,evl(nevn))
c$$$      endif
c$$$
c$$$
c$$$      entry suqlsw2(nband,jsp,nsp,evl,qlx) !takao added ql
c$$$C- Write or store the energy bands to file for this qp
c$$$      ql=qlx
c$$$      if (mode .eq. 1) then
c$$$        if (nblst .eq. 0) then
c$$$          prfmt = '(3f10.5,i6/(10f8.4))'
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.5/(10f8.4))'
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$      elseif (mode .eq. 2) then
c$$$        if (fmlong .and. nblst .eq. 0) then
c$$$          prfmt = '(3f15.10,i6/(5f15.10))'
c$$$        elseif (nblst .eq. 0) then
c$$$          prfmt = '(3f10.6,i6/(8f10.6))'
c$$$        elseif (fmlong .and. nblst .gt. 0) then
c$$$          prfmt = '(3f15.10/(5f15.10))'
c$$$        elseif (nblst .gt. 0) then
c$$$          prfmt = '(3f10.6/(8f10.6))'
c$$$C         if (nblst .le. 5) prfmt = '(8f10.6)'
c$$$        else
c$$$          call rx('bug in suqlsw')
c$$$        endif
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 3) then
c$$$        j = mod(iq-2,nqy)
c$$$        i = (iq-2-j)/nqy
c$$$Changenglob        stdo = nglob('stdo')
c$$$        stdo = globalvariables%stdo
c$$$        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c$$$c  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
c$$$
c$$$        do  k = 1, nblst
c$$$          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
c$$$        enddo
c$$$
c$$$C   ... If last qp generated, dump to file and exit
c$$$        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
c$$$          if (fmlong) then
c$$$            prfmt = '(5f15.10/(5f15.10))'
c$$$          else
c$$$            prfmt = '(8f10.6/(8f10.6))'
c$$$          endif
c$$$          rewind ifib
c$$$          do  j1 = 1, nsp
c$$$            do  k  = 1, nblst
c$$$              if (nsp .eq. 1)
c$$$     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
c$$$              if (nsp .eq. 2)
c$$$     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
c$$$     .        ifib,nqx,nqy,j1)
c$$$
c$$$              do  i = 1, nqx
c$$$                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          call rx0('finished generating bands on q-mesh')
c$$$        endif
c$$$      endif
c$$$
c$$$C ... Keep running tab on smallest,largest eval
c$$$      if (nevn .ne. 0) then
c$$$        evmxn = max(evmxn,evl(nevn))
c$$$        evmnn = min(evmnn,evl(nevn))
c$$$      endif
c$$$      end
c$$$
c$$$      subroutine suqlse(nband,jsp,nsp,ndimhx,ifbn,nfbn,ifblst,ndlst,
c$$$     .evec,wk)
c$$$      use m_globalvariables
c$$$C- Write to file the projection of eigenvector subblock for this qp
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   nband :number of energy bands to write
c$$$Ci   jsp   :current spin index (not used now)
c$$$Ci   nsp   :number of spins (not used now)
c$$$Ci   ndimhx:dimensions evec
c$$$Ci   ifbn  :index to color list (nfbn,ifblst)
c$$$Ci   nfbn  :number of elements for color weights projection
c$$$Ci   ifblst:list of elements for  color weights projection
c$$$Ci   ldlst :leading dimension of ifblst
c$$$Ci   evec  :eigenvectors
c$$$Ci   wk    :work array of same dimension as evl
c$$$Co Outputs
c$$$Cl Local variables
c$$$Cl         :
c$$$Cr Remarks
c$$$Cr   Inefficient, but it works
c$$$Cu Updates
c$$$Cu   08 Jul 08 New argument ndlst so ifblst can be dimensioned
c$$$Cu             independently from nband
c$$$Cu   05 Jun 06 First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer jsp,nsp,nband,ifbn,nfbn(2),ndlst,ifblst(ndlst,*)
c$$$      integer ndimhx
c$$$      double precision wk(1)
c$$$      double complex evec(ndimhx,ndimhx)
c$$$C ... Local parameters
c$$$      double complex zdotc
c$$$      character prfmt*40
c$$$      logical fmlong
c$$$      integer i,iprint,j,j1,k,nglob,stdo
c$$$      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
c$$$c      integer oqp
c$$$      double precision x1,x2,y1,y2,evmxn,evmnn,
c$$$     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)
c$$$
c$$$      common /suqlsd/
c$$$     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
c$$$     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c$$$c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy
c$$$
c$$$      integer ipiv(ndimhx)
c$$$      real(8),allocatable :: evsav(:,:,:,:)
c$$$      complex(8),allocatable:: evecc(:,:),work(:,:)
c$$$
c$$$      allocate(evecc(ndimhx,ndimhx),work(ndimhx,ndimhx))
c$$$
c$$$#if MULLIKEN
c$$$      call zcopy(ndimhx**2,evec,1,evecc,1)
c$$$      call zgetrf(ndimhx,ndimhx,evecc,ndimhx,ipiv,j)
c$$$      if (j .ne. 0) call rx('mullmf: failed to generate overlap')
c$$$      call zgetri(ndimhx,evecc,ndimhx,ipiv,work,ndimhx**2,j)
c$$$#endif
c$$$      do  i = 1, ndimhx
c$$$        wk(i) = 0
c$$$        do  j = 1, nfbn(ifbn)
c$$$          k = ifblst(j,ifbn)
c$$$          if (k .le. 0 .or. k .gt. ndimhx) then
c$$$            call fexit2(-1,111,' Exit -1 : suqlst: component %i'//
c$$$     .      ' outside range (1:%i)',k,ndimhx)
c$$$          endif
c$$$#if MULLIKEN
c$$$          wk(i) = wk(i) + evecc(i,k)*evec(k,i)
c$$$#else
c$$$          wk(i) = wk(i) + dconjg(evec(k,i))*evec(k,i)
c$$$#endif
c$$$        enddo
c$$$
c$$$#ifndef MULLIKEN
c$$$        wk(i) = wk(i) / zdotc(ndimhx,evec(1,i),1,evec(1,i),1)
c$$$        if (wk(i) .lt. 0 .or. wk(i) .gt. 1+1d-14) then
c$$$          call rx('bug in suqlse')
c$$$        endif
c$$$#endif
c$$$
c$$$      enddo
c$$$
c$$$      if (mode .eq. 1) then
c$$$        prfmt = '(3f10.5/(10f8.4))'
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) ql, (wk(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) ql, (wk(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 2) then
c$$$        if (fmlong) then
c$$$          prfmt = '(5f15.10)'
c$$$        else
c$$$          prfmt = '(8f10.6)'
c$$$C         if (nblst .gt. 0 .and. nblst .le. 5) prfmt = '(8f10.6)'
c$$$        endif
c$$$        if (nblst .eq. 0) then
c$$$          write(ifib,prfmt) (wk(i),i=1,nband)
c$$$        else
c$$$          write(ifib,prfmt) (wk(iblst(i)),i=1,nblst)
c$$$        endif
c$$$
c$$$      elseif (mode .eq. 3) then
c$$$        call rx('need copy weights into different place, mode=3')
c$$$        j = mod(iq-2,nqy)
c$$$        i = (iq-2-j)/nqy
c$$$Changenglob        stdo = nglob('stdo')
c$$$        stdo = globalvariables%stdo
c$$$        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c$$$  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)
c$$$
c$$$        do  k = 1, nblst
c$$$          evsav(i+1,j+1,k,jsp) = wk(iblst(k))
c$$$        enddo
c$$$
c$$$C   ... If last qp generated, dump to file and exit
c$$$        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
c$$$          if (fmlong) then
c$$$            prfmt = '(5f15.10/(5f15.10))'
c$$$          else
c$$$            prfmt = '(8f10.6/(8f10.6))'
c$$$          endif
c$$$          rewind ifib
c$$$          do  j1 = 1, nsp
c$$$            do  k  = 1, nblst
c$$$              if (nsp .eq. 1)
c$$$     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
c$$$              if (nsp .eq. 2)
c$$$     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
c$$$     .        ifib,nqx,nqy,j1)
c$$$
c$$$              do  i = 1, nqx
c$$$                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
c$$$              enddo
c$$$
c$$$            enddo
c$$$          enddo
c$$$          call rx0('finished generating bands on q-mesh')
c$$$        endif
c$$$      endif
c$$$
c$$$      deallocate(evecc,work)
c$$$
c$$$      end
c$$$
c$$$      subroutine suqlsr(mode,ifi,nsp,nbf,lde,nb1,nb2,ls1,ls2,nq,qp,eb)
c$$$C- Count qp, optionally read them and energy bands from file
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   mode  :Any combination of the following is allowed:
c$$$Ci         :1 count number of qp in file, return as nq
c$$$Ci         :  if 1s bit mode is 0, nq is an input,
c$$$Ci         :  the file nq cannot exceed input nq
c$$$Ci         :2 return qp read in file into array qp
c$$$Ci         :4 return bands read in file into array eb
c$$$Ci   ifi   :read from logical unit ifi
c$$$Ci   nsp   :number of spins in band file
c$$$Ci   nbf   :number of bands in band file
c$$$Ci   lde   :leading dimension of eb
c$$$Ci   nb1,nb2: read bands nb1..nb2 into eb
c$$$Ci   ls1   :(nsp=2): read first spin only
c$$$Ci   ls2   :(nsp=2): read second spin only
c$$$Cio Inputs/Outputs
c$$$Cio  nq    :numbefr of k-points (returned if 1s bit of mode is set)
c$$$Co Outputs
c$$$Co   qp    :k-points, returned if 4s bit of mode is set
c$$$Co   eb    :energy bands, returned if 4s bit of mode is set
c$$$Cl Local variables
c$$$Cl         :
c$$$Cr Remarks
c$$$Cr
c$$$Cu Updates
c$$$Cu   06 Jun 07  First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer mode,ifi,nsp,nq,nbf,lde,nb1,nb2
c$$$      logical ls1,ls2
c$$$      double precision qp(3,nq),eb(lde,nsp,nq)
c$$$C ... Local parameters
c$$$      integer nq1,nql,mode1,mode2,mode4,iq1,i,nline,isp
c$$$      double precision ql(3),evl(3000)
c$$$      logical :: l_dummy_isanrg,isanrg
c$$$
c$$$      if (mode .eq. 0) return
c$$$      mode1 = mod(mode,2)
c$$$      mode2 = mod(mode/2,2)
c$$$      mode4 = mod(mode/4,2)
c$$$      if (mode4 .ne. 0) then
c$$$Ckino isanrg is logical function,         call isanrg(nb2,1,min(lde,nbf),'suqlsr:','top band index',
c$$$        l_dummy_isanrg=isanrg(nb2,1,min(lde,nbf),'suqlsr:','top band index',
c$$$     .  .true.)
c$$$Ckino isanrg is logical function,         call isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
c$$$        l_dummy_isanrg=isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
c$$$      endif
c$$$
c$$$C --- For each panel, do ---
c$$$      nql = 0
c$$$      nline = 0
c$$$   91 continue
c$$$      read(ifi,*) nq1
c$$$      if (nq1 .le. 0) goto 90
c$$$      isp = 0
c$$$      do  iq1 = 1, nq1
c$$$        isp = mod(isp,2)+1
c$$$        nql = nql+1
c$$$        if (mode1 .eq. 0 .and. nql .gt. nq) call rxi(
c$$$     .  'suqlsr: file more q-points than allocated: nqmx=',nq)
c$$$        read(ifi,*,END=999,ERR=999) ql(1),ql(2),ql(3)
c$$$        if (mode2 .ne. 0) call dcopy(3,ql,1,qp(1,nql),1)
c$$$        read(ifi,*,END=999,ERR=999) (evl(i),i=1,nbf)
c$$$        if (mode4 .ne. 0) then
c$$$C         Copy only if appropriate spin
c$$$          if (isp .eq. 1 .and. ls2) then
c$$$          elseif (isp .eq. 2 .and. ls1) then
c$$$          else
c$$$            call dcopy(nb2-nb1+1,evl(nb1),1,eb(1,1,nql),1)
c$$$          endif
c$$$        endif
c$$$      enddo
c$$$      nline = nline+1
c$$$      goto 91
c$$$C     End of loop over lines
c$$$   90 continue
c$$$
c$$$
c$$$      if (mode1 .eq. 0) then
c$$$        call info2(30,1,0,' suqlsr: found %i qp in %i lines from file',
c$$$     .  nql,nline)
c$$$      else
c$$$        nq = nql/nsp
c$$$        call info2(30,1,0,' suqlsr: read %i qp in %i lines from file',
c$$$     .  nql,nline)
c$$$      endif
c$$$
c$$$      if (mod(nql,2) .ne. 0 .and. nsp .eq. 2) then
c$$$        call info0(20,0,0, ' suqlsr (warning): '//
c$$$     .  'odd number of bands encountered but nsp=2')
c$$$      endif
c$$$
c$$$      return
c$$$  999 continue
c$$$      call rxi('suqlsr: failed to read bands file, nq=',nql)
c$$$
c$$$      end
c$$$
c$$$      subroutine xxxbnd(evl,nbnd,nqx,nqy,ifib)
c$$$C     implicit none
c$$$      integer nbnd,nqx,nqy,ifib,ix,iy,ib,ifi,fopna
c$$$      double precision evl(nbnd,nqy,nqx)
c$$$      character*10 strn
c$$$
c$$$      do  10  ib = 1, nbnd
c$$$        ifi = ifib
c$$$        if (nbnd .ne. 1) then
c$$$          strn = ' '
c$$$          call awrit1('bnd%i',strn,10,0,ib)
c$$$          ifi = fopna(strn,30,0)
c$$$        endif
c$$$        write(ifi,337) nqx,nqy
c$$$  337   format('% rows ',i5,' cols ',i5)
c$$$        do  20  ix = 1, nqx
c$$$          write(ifi,333) (evl(ib,iy,ix), iy=1,nqy)
c$$$  333     format(6f12.6)
c$$$   20   continue
c$$$        if (nbnd .ne. 1) call fclose(ifi)
c$$$   10 continue
c$$$      end
c$$$      double precision function sclp(v1,v2)
c$$$      double precision v1(3),v2(3)
c$$$      sclp  = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
c$$$      end
c$$$
      
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs
c$$$  integer function rdtok(token,instr,descr,sep0,sep,cast,swt,
c$$$     .n1,count,resn,resc)
c$$$C- Seek a token within a string and read associated data into array
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   token: string marking data to follow token
c$$$Ci   instr: string in which to seek token
c$$$Ci   descr: descriptor string composed of one or more parts.
c$$$Ci        : It is used to format printouts generated by rdtok,
c$$$Ci        : in 'print' modes
c$$$Ci        :*The first part, if present, is a formatting string
c$$$Ci        : for writing the contents of a token, awrite conventions
c$$$Ci        : It is delimited by a double-colon, e.g.
c$$$Ci        :    descr = '%;6d,,string-for-remaining-parts'
c$$$Ci        : The format statement is used when the token's contents
c$$$Ci        : are written out (swt(5)>2, below).  A default format
c$$$Ci        : is used if one is not supplied
c$$$Ci        :*The second part, if present, is a character consisting
c$$$Ci        : of a single digit '0'..'9', and is the number of spaces
c$$$Ci        : to indent before writing the token, in all modes swt(5).
c$$$Ci        : It is delimited by a double-semicolon, e.g.
c$$$Ci        :    descr = '%;6d,,4;;string-for-remaining-part'
c$$$Ci        : If not present, a default indentation of 2 is used.
c$$$Ci        :*The third part, if present, describe the token's purpose.
c$$$Ci        : It is used only for input mode 0 (swt(5)=0).
c$$$Ci        : It consists of any part of the string after the first or
c$$$Ci        : second parts.
c$$$Ci    sep0: the character immediately to the left of the token
c$$$Ci          must be one of the characters in sep0
c$$$Ci    sep:  The token's arguments are delineated by characters in
c$$$Ci       :  this list.
c$$$Ci    cast: data following token is one of the following.
c$$$Ci          cast    class of token
c$$$Ci           0      logical
c$$$Ci           1      char
c$$$Ci           2      int
c$$$Ci           3      real
c$$$Ci           4      double
c$$$Ci        : result is stored into resn for cast=0,2,3,4
c$$$Ci        :           stored into resc for cast=1
c$$$Ci     n1 : Data associated with token are put into:
c$$$Ci          resn(n1),resn(n1+1),... for cast=0,2,3,4    and
c$$$Ci          resc(n1),resc(n1+1),... for cast=1
c$$$Ci   count: The number of elements (numbers or strings) to try and parse
c$$$Ci        : It may be permissible to read fewer than `count' elements;
c$$$Ci        : see swt(3)
c$$$Ci Inputs/Outputs
c$$$Cio  swt  : a vector of switches governing how gtcat and rdtok behave.
c$$$Cio       : Its internal structure is described in stswt, below.
c$$$Cio       : User is advised to set these parameters via routine stswt.
c$$$Cio       : and query the state of using quswt; then there is no need
c$$$Cio       : know the internal structure of swt.
c$$$Cio       : As described in detail in stst, quantities that can be set
c$$$Cio       : and probed are the following:
c$$$Cio       : structure of swt.
c$$$Cio       : swt()  function (see stswt)
c$$$Ci        :  1     offset to instr where parsing begins
c$$$Cio       :  2     parse up to this character
c$$$Ci        :  3     governs behavior when token or some data is missing,
c$$$Ci        :        whether tokens are optional or required, if
c$$$Ci        :        there may be multiple occurences of this token, and
c$$$Ci        :        whether number of elements may be less than count.
c$$$Ci        :  5     read and printout mode
c$$$Ci        :        0=> do not attempt to read token but show
c$$$Ci        :            what token would have been sought
c$$$Ci        :        1=> attempt to find token and read contents
c$$$Ci        :        2=> print token name
c$$$Ci        :        3=> combination of 1 and 2
c$$$Ci        :  6     file logical unit for printout
c$$$Co        :  7     index to last character parsed
c$$$Co Outputs
c$$$Co   rtdok: -999   No token was sought
c$$$Co        :   0    No match to token is found
c$$$Co        :   n    token matched and converted sans error n elements
c$$$Co        :        (for count=0, returns 1 if token matched)
c$$$Co        :  -n    if error on conversion of argument n
c$$$Cl Local variables
c$$$Cl   lopt :  0 if token is optional
c$$$Cl        :  1 if token is required
c$$$Cl        :  2 alternate token may be subst.
c$$$Cl        :  3 if token is ignored
c$$$Cb Bugs
c$$$Cb   when n1>0, data not placed in proper location for integer, logical
c$$$Cb   arrays
c$$$Cu Updates
c$$$Cu   06 Nov 01 Initially created
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer swt(7),cast,n1,count
c$$$      character*(1) instr(*)
c$$$      character*(*) token,descr,resc(1),sep0,sep
c$$$      double precision resn(1)
c$$$C ... Local parameters
c$$$      logical logval,parstr,l1,l2,lhelp,lshow,lread
c$$$      integer a2vec,awrite,i,ifipr,ifmt,io,ip,is,isave,it(99),ival,k,ks,
c$$$     .linstr,rdtk2,recl,retval,ltok,quswt,lfixl,lopt,lmult,lio
c$$$      parameter (recl=500)
c$$$      character   aa*(recl),prs*(recl)
c$$$      character*8 nmcast(0:4),fmt*20
c$$$      save ip,prs
c$$$      data nmcast /'logical','char','integer','real','double'/
c$$$      data ip /0/
c$$$
c$$$C ... Setup
c$$$      lopt  = quswt('token,opt',swt)
c$$$      if (lopt .eq. 3) then
c$$$        rdtok = -999
c$$$        return
c$$$      endif
c$$$      ifipr = quswt('lunit',swt)
c$$$      ltok = len(token)
c$$$      lhelp = quswt('io,help',swt) .eq. 1
c$$$      lshow = quswt('io,w',swt) .eq. 1
c$$$      lread = quswt('io,r',swt) .eq. 1
c$$$      lio   = quswt('io',swt)
c$$$      lfixl = quswt('token,fixlen',swt)
c$$$      lmult = quswt('token,mult',swt)
c$$$
c$$$C --- Show what token would have been sought ---
c$$$      if (lhelp) then
c$$$
c$$$C       Find indentation, offset to descriptor string
c$$$C       Default indentation
c$$$        ks = 2
c$$$C       Default offset to start of descriptor
c$$$        i = 0
c$$$        if (descr .ne. ' ') then
c$$$C         Strip off first parts, if either exists
c$$$          l2 = parstr(descr,';;',len(descr)-1,2,'a',i,k)
c$$$          if (l2) then
c$$$            ks = ichar(descr(i:i)) - ichar('0')
c$$$          else
c$$$            i = 0
c$$$            l1 = parstr(descr,',,',len(descr)-1,2,'a',i,k)
c$$$            if (.not. l1) i = -2
c$$$          endif
c$$$        endif
c$$$        fmt = ' '
c$$$        is= ltok
c$$$        if (count .eq. 0) then
c$$$          aa = fmt(1:ks) // ' *token  ' // token(1:is) //
c$$$     .    '%a (no arguments; %?#(n==0)#optional#required#)'
c$$$          call awrit1(aa,' ',72,ifipr,lopt)
c$$$        elseif (count .eq. 1) then
c$$$          aa = fmt(1:ks) // ' *token  ' // token(1:is) //
c$$$     .    '%a takes 1 argument of cast '//nmcast(cast)//
c$$$     .    '%a %?#(n==0)#(optional)##%-1j'//
c$$$     .    '%a %?#(n==1)#(required)##%-1j'//
c$$$     .    '%a %?#(n==2)#%N'//fmt(1:ks)//
c$$$     .    '   (choose either this token or the next one)##'
c$$$          call awrit1(aa,' ',128,ifipr,lopt)
c$$$        else
c$$$          aa = fmt(1:ks) // ' *token  ' // token(1:is) //
c$$$     .    '%a takes %?#(n==0)#up to ##'//
c$$$     .    '%i arguments of cast '//nmcast(cast)//
c$$$     .    '%a %?#(n==0)#(optional)##%-1j'//
c$$$     .    '%a %?#(n==1)#(required)##%-1j'//
c$$$     .    '%a %?#(n==2)#%N'//fmt(1:ks)//
c$$$     .    '   (choose either this token or the next one)##'
c$$$          call awrit3(aa,' ',128,ifipr,lfixl,count,lopt)
c$$$        endif
c$$$        if (descr .ne. ' ') then
c$$$          fmt = ' '
c$$$          aa = '   ' // fmt(1:ks) // descr(i+3:)
c$$$          call awrit0(aa,' ',-len(aa),ifipr)
c$$$        endif
c$$$        rdtok = -999
c$$$        ip = 0
c$$$        return
c$$$
c$$$C --- No token to read ---
c$$$      elseif (.not. lread) then
c$$$        retval = -999
c$$$        goto 999
c$$$      endif
c$$$
c$$$C --- Try and find token ---
c$$$      linstr = swt(2)
c$$$      if (token .eq. ' ') then
c$$$        is = swt(1)
c$$$      else
c$$$        is = rdtk2(instr,linstr,token,sep0,lopt.eq.1,swt(1))
c$$$      endif
c$$$C ... No token found ... exit
c$$$      if (is .le. 0) then
c$$$        swt(7) = linstr
c$$$        rdtok = 0
c$$$        return
c$$$      endif
c$$$
c$$$C --- Parse for token arguments ---
c$$$      swt(7) = is
c$$$
c$$$C ... Case count=0
c$$$      if (count .eq. 0) then
c$$$        retval = 1
c$$$        goto 999
c$$$      endif
c$$$
c$$$C ... Case further delimit range of search
c$$$      if (lmult .ne. 0) then
c$$$        ks = rdtk2(instr,linstr,token,sep0,.false.,is)
c$$$        if (ks .gt. 0) linstr = ks-ltok-1
c$$$        swt(2) = linstr
c$$$      endif
c$$$
c$$$C ... Parse for vector of character arguments
c$$$      if (cast .eq. 1) then
c$$$        isave = is - ltok
c$$$        do  i = 1, count
c$$$
c$$$C         Skip past leading blanks
c$$$          is = is-1
c$$$          call skipbl(instr,linstr,is)
c$$$          is = is+1
c$$$
c$$$C         End-of-category encountered before loop over count exhausted
c$$$          if (is .ge. linstr) then
c$$$            if (lfixl .ne. 0) then
c$$$              if (count .eq. 1) then
c$$$                aa = 'Exit -1 rdtok: no string after token '//
c$$$     .          token(1:ltok)
c$$$                call fexit(-1,1,aa,0)
c$$$              else
c$$$                aa = 'Exit -1 rdtok: failed to parse %i strings, token'
c$$$                call strncp(aa,instr,53,isave,min(recl-53,linstr-isave)-1)
c$$$              endif
c$$$              call fexit(-1,1,aa,count)
c$$$            endif
c$$$            retval = -i
c$$$            swt(7) = linstr+1
c$$$            goto 999
c$$$          endif
c$$$
c$$$          call cpstr(instr,linstr,1001,sep,is,io,resc(i+n1-1))
c$$$          is = is+1
c$$$        enddo
c$$$        swt(7) = is-1
c$$$        retval = count
c$$$        goto 999
c$$$      endif
c$$$
c$$$C ... Parse for vector of numerical arguments
c$$$      isave = is - ltok
c$$$      swt(7) = is
c$$$      is = is-1
c$$$      k = a2vec(instr,linstr,is,cast,sep,len(sep),-99,-count,it,
c$$$     .resn(n1))
c$$$      retval = k
c$$$
c$$$C ... a2vec encountered an illegitimate expression
c$$$      if (k .lt. 0) then
c$$$        if (lfixl .ne. 0) then
c$$$          aa = ' Exit -1 rdtok: error parsing expr no %i, token'
c$$$          call strncp(aa,instr,52,isave,min(recl-52,linstr-isave)-1)
c$$$          call fexit(-1,1,aa,-k)
c$$$        endif
c$$$C ... a2vec read fewer than count expressions
c$$$      elseif (k .lt. count) then
c$$$        if (lfixl .ne. 0) then
c$$$          aa = ' Exit -1 rdtok: read only %i of %i values, token'
c$$$          call strncp(aa,instr,52,isave,min(recl-52,linstr-isave)-1)
c$$$          call fexit2(-1,1,aa,k,count)
c$$$        endif
c$$$        retval = -k
c$$$      endif
c$$$
c$$$C --- Cleanup: check for printout ---
c$$$  999 continue
c$$$      rdtok = retval
c$$$      if (.not. lshow) return
c$$$
c$$$C     Find format, if supplied, and indentation
c$$$C     Default indentation
c$$$      ks = 2
c$$$      fmt = ' '
c$$$      if (descr .ne. ' ') then
c$$$        i = 0
c$$$        l2 = parstr(descr,';;',len(descr)-1,2,'a',i,k)
c$$$        if (l2) then
c$$$          ks = ichar(descr(i:i)) - ichar('0')
c$$$        endif
c$$$        ifmt = 0
c$$$        l1 = parstr(descr,',,',len(descr)-1,2,'a',ifmt,k)
c$$$        if (l1) fmt = descr(1:ifmt)
c$$$      endif
c$$$
c$$$C ... Write token name to string
c$$$      if (ip .eq. 0) then
c$$$        prs = ' '
c$$$      else
c$$$        ks = min(ks,2)
c$$$      endif
c$$$      ip = ks + ip + awrite(token,prs(ks+ip+1:),128,0,0,0,0,0,0,0,0,0)
c$$$      if (count .eq. 0) then
c$$$C   ... format applied after token write
c$$$C        if (fmt .ne. ' ')
c$$$C     .  ip = awrite(fmt,prs,128,0,0,0,0,0,0,0,0,0)
c$$$        goto 998
c$$$      endif
c$$$
c$$$C ... Write token contents to string, vec of length |retval|
c$$$C     ip = ip+1
c$$$
c$$$C     Number of elements to write
c$$$      k = retval
c$$$      if (retval .lt. 0) k = -retval-1
c$$$      if (lio .eq. 2) k = count
c$$$C     Write the elements
c$$$      do  i = 1, k
c$$$        if (cast .eq. 1) then
c$$$          ip = ip + awrite(resc(n1+i-1),prs(ip+1:),128,
c$$$     .    0,0,0,0,0,0,0,0,0)
c$$$          call skpblb(prs,ip,ip)
c$$$          ip = ip+1
c$$$          if (fmt .ne. ' ')
c$$$     .    ip = awrite(fmt(1:ifmt),prs,128,0,0,0,0,0,0,0,0,0)
c$$$        elseif (cast .eq. 0) then
c$$$          if (fmt .eq. ' ') then
c$$$            fmt = '%l'
c$$$            ifmt = 2
c$$$          endif
c$$$          ip = ip + awrite(fmt(1:ifmt),prs(ip+1:),128,0,
c$$$     .    logval(resn,n1+i-1),0,0,0,0,0,0,0)
c$$$        elseif (cast .eq. 2) then
c$$$          if (fmt .eq. ' ') then
c$$$            fmt = '%i'
c$$$            ifmt = 2
c$$$          endif
c$$$          ip = ip + awrite(fmt(1:ifmt),prs(ip+1:),128,0,
c$$$     .    ival(resn,n1+i-1),0,0,0,0,0,0,0)
c$$$        elseif (cast .eq. 4) then
c$$$          if (fmt .eq. ' ') then
c$$$            fmt = '%;10F'
c$$$            ifmt = 5
c$$$          endif
c$$$          ip = ip + awrite(fmt(1:ifmt),prs(ip+1:),128,0,resn(n1+i-1),
c$$$     .    0,0,0,0,0,0,0)
c$$$        endif
c$$$C       ip = ip+2
c$$$        if (i .lt. k) then
c$$$          if (prs(ip:ip) .ne. sep(1:1)) then
c$$$            ip = ip+1
c$$$            write(prs(ip:ip),'(a)') sep(1:1)
c$$$          endif
c$$$        elseif (retval .lt. 0 .and. lio .ne. 2) then
c$$$          ip = ip+1
c$$$          write(prs(ip:),'(a,''***'')') sep(1:1)
c$$$          ip = ip+3
c$$$        endif
c$$$      enddo
c$$$
c$$$  998 continue
c$$$c     call skpblb(prs,ip,ip)
c$$$C     ip = ip+1
c$$$
c$$$C ... Printout, or just preserve in string prs if ifipr is 0
c$$$      if (ifipr .gt. 0) then
c$$$        write(ifipr,'(a)') prs(1:ip)
c$$$        ip = 0
c$$$      else
c$$$      endif
c$$$
c$$$      end
c$$$
c$$$      integer function rdtk2(instr,linstr,token,sep0,lreqd,i1)
c$$$C- Kernel called by rdtok to find a token
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   instr :string in which to seek token
c$$$Ci   linstr:length of instr
c$$$Ci   token :token
c$$$Ci    sep0 :the character immediately to the left of the token
c$$$Ci         :must be one of the characters in sep0
c$$$Ci   lreqd :T, abort with error if token is not found
c$$$Ci    i1   :start search at i1
c$$$Co Outputs
c$$$Co   rdtk2 :if token found, rdtk2 points to first character past token
c$$$Co         :if no token found, rdtok = -linstr
c$$$Cu Updates
c$$$Cu   02 Nov 01 First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      logical lreqd
c$$$      integer linstr,i1
c$$$      character*(*) token,sep0
c$$$      character*(1) instr(linstr)
c$$$C ... Local parameters
c$$$      integer ltok,k,iterm,is,lgunit,stderr
c$$$      character ct*1
c$$$
c$$$      ltok = len(token)
c$$$      ct = token
c$$$      stderr = lgunit(2)
c$$$
c$$$C --- Loop over length of string ---
c$$$C     Match by:
c$$$C     1.  First character matches at instr(k)
c$$$C     2.  character to left matches on of sep0
c$$$C     3.  token matches
c$$$      do  k = i1, linstr
c$$$
c$$$        if (instr(k) .eq. ct) then
c$$$          iterm = 1
c$$$          if (k .gt. 1) then
c$$$            is = 0
c$$$            call chrps2(instr(k-1),sep0,len(sep0),0,is,iterm)
c$$$          endif
c$$$          if (iterm .gt. 0) then
c$$$            do  is = 1, ltok
c$$$              if (instr(k+is-1) .ne. token(is:is)) goto 31
c$$$            enddo
c$$$C           A match was found
c$$$            rdtk2 = ltok+k
c$$$            return
c$$$          endif
c$$$        endif
c$$$   31   continue
c$$$      enddo
c$$$C ... We have fallen through the loop, so no token found
c$$$      rdtk2 = -linstr
c$$$      if (lreqd) then
c$$$        write(stderr,'(1x)')
c$$$        call rxs2('rdtok: failed to find token `',token,'''')
c$$$      endif
c$$$      end
c$$$

cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
c$$$  subroutine a2rotm(strn,linv,ipr,rotm)
c$$$C- Generate Euler angles from a sequence of rotations, ascii input
c$$$Ci strn: Sequence of rotation matrices.  Syntax: rot1[,rot2...],
c$$$Ci   where each rotj looks like (x,y,z)angle.  x:, y: or z: may be
c$$$Ci   substituted for (x,y,z).  Example:  strn that makes rotation matrix
c$$$Ci   corresponding to the Euler angles alpha=pi/4 beta=pi/3 gamma=pi/2:
c$$$Ci   (0,0,1)pi/4,(0,1,0)pi/3,(0,0,1)pi/2   or just  z:pi/4,y:pi/3,z:pi/2
c$$$Ci linv: T, return transpose of rotm, corresponding to rotation of
c$$$Ci   positions rather than the coordinates (default).
c$$$Ci ipr: verbosity
c$$$Co rotm: rotation matrix.  Follow with a call to rm2eua to find Euler
c$$$Co   angles.
c$$$C     implicit none
c$$$      logical linv
c$$$      character*(*) strn, rchr*7, ss*80
c$$$      double precision rotm(3,3)
c$$$      double precision alpha,beta,gamma,rotj(3,3),r2(3,3),angle,
c$$$     .r(3),pi,ddot,fuzz,ct,theta,phi
c$$$
c$$$      integer ls,itrm,ich,a2vec,ix(3),i,j,irot,ipr,i1mach
c$$$      parameter (fuzz = 1d-12)
c$$$      data rchr /'(XxYyZz'/
c$$$
c$$$
c$$$      pi = 4*datan(1d0)
c$$$      call dpzero(rotm,9)
c$$$      rotm(1,1) = 1
c$$$      rotm(2,2) = 1
c$$$      rotm(3,3) = 1
c$$$      ss = strn
c$$$      ls = len(ss)
c$$$      irot = 0
c$$$      ich = 0
c$$$C --- Entry point to accumulate new rotation ---
c$$$    5 continue
c$$$      irot = irot+1
c$$$      call chrps2(ss,rchr,len(rchr),ich,ich,itrm)
c$$$      if (itrm .eq. 0) goto 999
c$$$      if (itrm .gt. 1) then
c$$$        call dpzero(r,3)
c$$$        if (ss(ich+2:ich+2) .ne. ':') goto 999
        itrm = itrm/2
        if (itrm .gt. 3) call rx('bug in a2rotm')
   11   r(itrm) = 1
        ich = ich+2
      else
        ich = ich+1
        if (a2vec(ss,ls,ich,4,',)',2,2,3,ix,r) .ne. 3) goto 999
      endif
      if (a2vec(ss,ls,ich,4,', ',2,2,1,ix,angle) .ne. 1) goto 999
      call dpzero(rotj,9)
C ... Rotation about r
      call dscal(3,1/dsqrt(ddot(3,r,1,r,1)),r,1)
      ct = r(3)
      phi = 0
      if ( ct .gt.  1 .and.  ct .lt. 1-fuzz) ct = 1
      if (-ct .gt.  1 .and. -ct .lt. 1-fuzz) ct = -1
      if (dabs(ct) .lt. 1) phi = datan2(r(2),r(1))
      theta = dacos(ct)
      if (theta .lt. 0) theta = 2*pi - theta
      call rotma(phi,theta,angle,rotj)
C ... rotm <- rotj*rotm
      call dpcopy(rotm,r2,1,9,1d0)
      call dmpy(rotj,3,1,r2,3,1,rotm,3,1,3,3,3)
      if (ipr .gt. 40) then
        call rm2eua(rotm,alpha,beta,gamma)
        call awrit4(' a2rotm:  rotation %i, alpha = %1;6d  '//
     .  'beta = %1;6d  gamma = %1;6d',' ',80,i1mach(2),
     .  irot,alpha,beta,gamma)
        print 335, ((rotm(i,j),j=1,3),i=1,3)
  335   format((3f15.9))

C ... code to make beta, phi om from Euler angles
C        st = dsqrt(1-ct**2)
C        print *, 'beta=',2d0*dasin(st*dsin(angle))
C        print *, 'theta,phi=',theta,phi
C        print *, 'phi=',phi,(alpha-gamma-pi)/2
C        if (dabs((phi-(alpha-gamma-pi)/2)) .gt. 1d-10) then
C          print *, 'warning: phi off', need fix this
C          phi = phi - pi
C        endif
C        if (dabs((phi-(alpha-gamma-pi)/2)) .gt. 1d-10)  stop 'phi off'
C        gp = (alpha+gamma)/2
C        print *, 'gp=',gp
C        print *, 'angle=', angle,
C     .    dasin(dsqrt((datan(gp)**2+dsin(beta/2)**2)/(datan(gp)**2+1)))
C        tg = dcos(theta)*dtan(angle)
C        tg = dtan(gp)
C        print *, dcos(angle)**2*tg**2+dsin(beta/2)**2-dsin(angle)**2,
C     .    dtan(gp),tg
C        print *, tg**2+dsin(beta/2)**2-dsin(angle)**2*(1+tg**2)
C        print *, (tg**2+dsin(beta/2)**2)/(1+tg**2)-dsin(angle)**2
C        print *, dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2))-dsin(angle)
C        print *, 'angle=', angle,
C     .    dasin(dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2)))
CC ...   sin and cos angle, each undetermined up to a sign
C        sw = dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2))
C        cw = dsqrt(1-sw**2)
C        print *, 'sin angle=', sin(angle),sw
C        if (tg .lt. 0) cw = -cw
CC ...   theta, assuming ??? cw positive
C        th = datan2(dsin(beta/2),cw*tg)
C        print *, 'theta=', theta, datan2(dsin(beta/2),cw*tg)
C        print *, 'theta=', theta, datan2(dsin(beta/2),-cw*tg)
C        st = dsin(th)
C        ct = dcos(th)
C        om = datan2(tg,ct)
C
Cc ...   check
C        print *, 'ck tg', tg,tan(om)*ct
C        print *, 'ck sb', sin(beta/2),sin(th)*sin(om)
C
CC ok for ROT=(.1,.2,.3).9*pi
CC but backwards for ROT=(.1,.2,-.3)-.9*pi ROT=(.1,-.2,-.3)-.9*pi
C        stop
C
C        call euler(r(1),r(2),r(3),angle,alpha,beta,gamma)

      endif

C --- Pick up another rotation matrix, or quit if done ---
      if (ix(1) .eq. 1) goto 5
      if (linv) then
        call dpcopy(rotm,r2,1,9,1d0)
        do    i = 1, 3
          do   j = 1, 3
             rotm(i,j) = r2(j,i)
          enddo
          enddo
        if (ipr .ge. 50) then
          call rm2eua(rotm,alpha,beta,gamma)
          call awrit3(' a2rotm:  inverse, alpha = %1;6d  '//
     .    'beta = %1;6d  gamma = %1;6d',' ',80,i1mach(2),
     .    alpha,beta,gamma)
          print 335, ((rotm(i,j),j=1,3),i=1,3)
        endif
      endif
      return

  999 call rx('a2rotm: rotation must be one of x:#,y:#,z:#,(vec)#')
      end
      subroutine eua2rm(alpha,beta,gamma,r)
C- Generate the rotation matrix corresponding to Euler angles
C ----------------------------------------------------------------------
C  Definition: R = Rz(gamma)*Ry(beta)*Rz(alpha)
Cr This definition has the property that a vector v rotated, i.e.
Cr    R v = v'
Cr is rotated to v' = zhat 
Cr when v is defined by the polar coordinates (alpha,beta)
Cr    vx = cos(alpha) sin(beta)
Cr    vy = sin(alpha) sin(beta)
Cr    vx = cos(beta)
Cr Note also that NB: r^-1 = r+
C ----------------------------------------------------------------------
C     implicit none
      double precision r(3,3),alpha,beta,gamma
C     integer i,j
      double precision ca,cb,cg,sa,sb,sg

      ca = dcos(alpha)
      sa = dsin(alpha)
      cb = dcos(beta)
      sb = dsin(beta)
      cg = dcos(gamma)
      sg = dsin(gamma)

C     Euler angles: R3(=R_z(gamma)) R2(=R_y(beta)) R1(=R_z(alpha))
      r(1,1) =  cg*cb*ca - sg*sa
      r(2,1) = -sg*cb*ca - cg*sa
      r(3,1) =  sb*ca
      r(1,2) =  cg*cb*sa + sg*ca
      r(2,2) = -sg*cb*sa + cg*ca
      r(3,2) =  sb*sa
      r(1,3) = -cg*sb
      r(2,3) =  sg*sb
      r(3,3) =  cb

C      print 335, ((r(i,j),j=1,3),i=1,3)
C  335 format((3f15.9))
      end
      subroutine rm2eua(r,alpha,beta,gamma)
C- Generate Euler angles from a rotation matrix
C     implicit none
      double precision r(3,3),alpha,beta,gamma,fuzz,pi
      integer i1mach,i,j
      double precision ca,cb,cg,sa,sb,sg,amg,apg,xx,yy,d1mach
      parameter (fuzz = 1d-8)

      pi = 4*datan(1d0)
      cb = r(3,3)
      if (dabs(cb) .gt. 1 .and.  dabs(cb) .lt. 1+fuzz) cb = sign(1d0,cb)
C --- Handle case beta is near zero or pi ---
      if (dabs(dabs(cb)-1) .lt. dsqrt(fuzz)) then
        xx = r(1,3)**2 + r(2,3)**2
        yy = r(3,1)**2 + r(3,2)**2
        beta = datan2(dsqrt((xx+yy)/2),r(3,3))
      else
        beta = dacos(cb)
      endif
      if (dsqrt(r(3,2)**2 + r(3,1)**2) .lt. fuzz .or.
     .dsqrt(r(2,3)**2 + r(1,3)**2) .lt. fuzz) then
        if (cb .gt. 0) then
          apg = datan2(r(1,2),r(2,2))
          if (r(3,1)**2+r(3,2)**2 .ne. 0d0) then
            alpha = datan2(r(3,2),r(3,1))
            gamma = apg - alpha
          elseif (r(1,3)**2+r(2,3)**2 .ne. 0d0) then
            gamma = datan2(r(2,3),-r(1,3))
            alpha = apg - gamma
          else
            gamma = 0d0
            alpha = apg
          endif
        else
          amg = datan2(-r(1,2),r(2,2))
          if (r(3,1)**2+r(3,2)**2 .ne. 0d0) then
            alpha = datan2(r(3,2),r(3,1))
            if (dsqrt(r(3,1)**2+r(3,2)**2) .lt. d1mach(3)*2) alpha = 0
            gamma = alpha - amg
          elseif (r(1,3)**2+r(2,3)**2 .ne. 0d0) then
            gamma = datan2(r(2,3),-r(1,3))
            alpha = amg + gamma
          else
            gamma = 0
            alpha = amg
          endif
        endif
      else
        alpha = datan2(r(3,2),r(3,1))
        gamma = datan2(r(2,3),-r(1,3))
      endif
      if (alpha .gt.  pi) alpha = alpha - 2*pi
      if (alpha .lt. -pi) alpha = alpha + 2*pi
      if (gamma .gt.  pi) gamma = gamma - 2*pi
      if (gamma .lt. -pi) gamma = gamma + 2*pi
C      call awrit3(' rm2eua:  alpha = %1;6d  beta = %1;6d'//
C     .  '  gamma = %1;6d',' ',80,i1mach(2),alpha,beta,gamma)
C ... Check validity of entire rotation matrix
      ca = dcos(alpha)
      sa = dsin(alpha)
      cb = dcos(beta)
      sb = dsin(beta)
      cg = dcos(gamma)
      sg = dsin(gamma)
C      print *,   dabs(r(1,1)-(ca*cb*cg-sa*sg)) .gt. fuzz*10 ,
C     .    dabs(r(1,2)-(sa*cb*cg+ca*sg)) .gt. fuzz*10 ,
C     .    dabs(r(1,3)-(-sb*cg)) .gt. fuzz*10 ,
C     .    dabs(r(2,1)-(-ca*cb*sg-sa*cg)) .gt. fuzz*10 ,
C     .    dabs(r(2,2)-(-sa*cb*sg+ca*cg)) .gt. fuzz*10 ,
C     .    dabs(r(2,3)-(sb*sg)) .gt. fuzz*10 ,
C     .    dabs(r(3,1)-(ca*sb)) .gt. fuzz*10 ,
C     .    dabs(r(3,2)-(sa*sb)) .gt. fuzz*10 ,
C     .    dabs(r(3,3)-cb) .gt. fuzz*10
      if (dabs(r(1,1)-(ca*cb*cg-sa*sg)) .gt. fuzz*10 .or.
     .dabs(r(1,2)-(sa*cb*cg+ca*sg)) .gt. fuzz*10 .or.
     .dabs(r(1,3)-(-sb*cg)) .gt. fuzz*10 .or.
     .dabs(r(2,1)-(-ca*cb*sg-sa*cg)) .gt. fuzz*10 .or.
     .dabs(r(2,2)-(-sa*cb*sg+ca*cg)) .gt. fuzz*10 .or.
     .dabs(r(2,3)-(sb*sg)) .gt. fuzz*10 .or.
     .dabs(r(3,1)-(ca*sb)) .gt. fuzz*10 .or.
     .dabs(r(3,2)-(sa*sb)) .gt. fuzz*10 .or.
     .dabs(r(3,3)-cb) .gt. fuzz*10) then
        call awrit3(' rm2eua found alpha=%1;9d  beta=%1;9d'//
     .  '  gamma=%1;9d from rot:',' ',80,i1mach(2),alpha,beta,gamma)
        print 335, ((r(i,j),j=1,3),i=1,3)
  335   format((3f20.15))
        call rx('rm2eua: matrix not consistent with Euler angles: numerical error-->fuzz too large?')
      endif
      end
      subroutine rotma(phi,theta,gamma,rotm)
C- Generate rotation matrix for rotation about specified polar angle
C  Rotates coordinates by gamma about (phi,theta)
C     implicit none
      double precision rotm(3,3),r2(3,3),r1(3,3),phi,theta,gamma,xx
C     integer i,j

      call dpzero(r1,9)
      r1(1,1) = dcos(phi)
      r1(2,2) = r1(1,1)
      r1(2,1) = -dsin(phi)
      r1(1,2) = -r1(2,1)
      r1(3,3) = 1
      call dpzero(r2,9)
      r2(2,2) = 1
      r2(1,1) = dcos(theta)
      r2(3,3) = r2(1,1)
      r2(3,1) = dsin(theta)
      r2(1,3) = -r2(3,1)
c ... rotm <- Rtheta*Rphi
      call dmpy(r2,3,1,r1,3,1,rotm,3,1,3,3,3)
      call dpzero(r2,9)
      r2(1,1) = dcos(gamma)
      r2(2,2) = r2(1,1)
      r2(2,1) = -dsin(gamma)
      r2(1,2) = -r2(2,1)
      r2(3,3) = 1
c ... r1 <- Rgamma*(Rtheta*Rphi)
      call dmpy(r2,3,1,rotm,3,1,r1,3,1,3,3,3)
C ... r2 <- (Rtheta*Rphi)^-1, rotm <- (Rtheta*Rphi)^-1*Rg*(Rtheta*Rphi)
      call dinv33(rotm,0,r2,xx)
      call dmpy(r2,3,1,r1,3,1,rotm,3,1,3,3,3)
C      print 336, ((rotm(i,j),j=1,3),i=1,3)
C  336 format(' rot matrix (Rtheta*Rphi)^-1*Rg*(Rtheta*Rphi):'/(3f12.6))
      end
C      subroutine euler(x,y,z,om,al,be,ga)
C note this convention rotates backwards from usual one.
C      implicit none
C      double precision x,y,z,om,al,be,ga
C      double precision pi,fi,st,ct,ay,tg,om2,r,atg
C      data pi /3.141592653589793d0/
C      om2 = om
C      r = dsqrt(x*x + y*y + z*z)
C      x = x/r
C      y = y/r
C      z = z/r
C      ay = dabs(y)
C      ct = z
C      st = dsqrt(x*x + y*y)
C      if (st .eq. 0d0) then
C        fi = 0d0
C      else
C        fi = dasin(ay/st)
C      endif
C      if (x .gt. 0d0 .and. y .lt. 0d0) fi = 2d0*pi - fi
C      if (x .lt. 0d0 .and. y .gt. 0d0) fi =      pi - fi
C      if (x .lt. 0d0 .and. y .lt. 0d0) fi =      pi + fi
C      be = 2d0*dasin(st*dsin(om2))
C      if (dabs(om/pi-0.5d0) .lt. 0.001d0) then
C        if (ct .ne. 0d0) then
C          atg = pi*dsign(1d0,(0.5d0-om/pi)*ct)/2d0
C        else
C          atg = fi - pi / 2d0
C        end if
C      else
C        tg = ct*dsin(om2)/dcos(om2)
C        atg = datan(tg)
C      endif
C      al = atg + fi - pi/2d0
C      ga = atg - fi + pi/2d0
C      if (ga .lt. -pi) ga = ga+2*pi
C
C      print *, 'fi,be,gp,=',fi,be,atg
C      tg = dtan(atg)
C
C      print *, 'om=',om,datan(tg/ct)
C      print *, dcos(om)**2*tg**2+sin(be/2)**2-dsin(om)**2
C
C      print 1,al,be,ga
C    1 format(' al=',f10.6,' be=',f10.6,' gam=',f10.6)
C      stop
C      end



cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs      
  module m_pvioeu
      contains
      subroutine pvioeu(mode,ssite,eula,nbas,neul)
      use m_struc_def  
C- Writes the Euler angles to a file or to stdout
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 write in format suitable for rdm
Ci         :1 write in printout format
Ci         :10s digit
Ci         :0 printout is for Euler angles
Ci         :1 printout is for Bfield
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    clabel
Ci   eula  :Euler angles for noncollinear spins
Ci   nbas  :size of basis
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Co Outputs
Co   Writes to stdout the Euler angles
Cr Remarks
Cr
Cu Updates
Cu   14 Feb 03 output can be for euler angles or for B-field
Cu   22 May 02 Writes label information if passed through in ssite
Cu             Altered argument list
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,neul
      real(8):: eula(nbas,neul,3)
      type(s_site),optional::ssite(*)
C ... Local parameters
      integer i,j,k,lgunit,stdo,mode0,mode1,lbl
      double precision rotm(3,3)
      character*72 outs, clabl*8
      character*255:: formatc(344:346),lbll
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      stdo = lgunit(1)
      clabl = ' '
      if (mode0 .ne. 0) then
        lbl=344
        if (mode1 .ne. 0) lbl=345
        if (mode1 .eq. 0 .and. mode0 .ne. 0) lbl=346
        formatc(344)="('   ib     alpha        beta       gamma':6x,a)"
        formatc(345)="('   ib      bx          by          bz  ':6x,a)"
        formatc(346)="('   ib   alpha      beta     gamma':5x,a,"//
     .  "11x,'z axis rotated to:')"
        lbll=formatc(lbl)

        if ( neul .gt. 1)  write(outs,trim(lbll))
c        if (ssite(1)%size .eq. 0 .or. neul .gt. 1)  write(outs,trim(lbll))
c        if (ssite(1)%size .ne. 0 .and. neul .eq. 1) write(outs,trim(lbll)) 'class'
      endif
c  344 format('   ib     alpha        beta       gamma':6x,a)
c  345 format('   ib      bx          by          bz  ':6x,a)
c  346 format('   ib   alpha      beta     gamma':5x,a,
c     .  11x,'z axis rotated to:')
      if (neul .gt. 4) then
        call awrit0('%a     lm',outs,len(outs),0)
      elseif (neul .gt. 1) then
        call awrit0('%a      l',outs,len(outs),0)
      endif
      call awrit0(outs,' ',-len(outs),stdo)
      do  21  i = 1, nbas
c        if (ssite(1)%size .ne. 0) then
c          do i_spacks=i,i
c            call spacks_copy('u',ssite(i_spacks)%clabel,i,i,clabl,i_spacks)
c          enddo
c        endif
        if (neul .gt. 1 .and. mode0 .ne. 0) then
          if (clabl .eq. ' ') call awrit1('# site %,4i',' ',80,stdo,i)
          if (clabl .ne. ' ') call awrit1('# site %,4i  class '//clabl,
     .    ' ',80,stdo,i)
          write(stdo,'(i5,3f12.6,i4)')
     .    (i,(eula(i,j,k),k=1,3),j, j=1,neul)
        elseif (neul.eq.1 .and. mode0.ne.0 .and. clabl.eq.' ') then
          write(stdo,'(i5,3f12.6)')
     .    (i,(eula(i,j,k),k=1,3), j=1,neul)
        elseif (neul .eq. 1 .and. mode0 .ne. 0) then
          call eua2rm(eula(i,1,1),eula(i,1,2),eula(i,1,3),rotm)
          write(stdo,'(i5,3f10.6,3x,a,1x,3f10.6)')
     .    i,(eula(i,1,k),k=1,3), clabl, (rotm(3,k),k=1,3)
        else
          write(stdo,'(3f16.12)') ((eula(i,j,k),k=1,3),j=1,neul)
        endif
   21 continue

      end subroutine pvioeu
      end module m_pvioeu

      subroutine ioeula(nbas,nl,eula,neul,xsi,ifi)
C- I/O of Euler angles
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: <0 for write, >0 for read
Ci   nbas: number of basis atoms
Co Outputs
Co   Euler angles are read in or written out
Co   neul  is read in or written out (either 1, nl or nl**2)
Co   xsi is read in if xsi=# is available on the first line
Co   xsi is written out if xsi is nonzero
Cr Remarks
Cr   Aborts on read when nbas does not match file
Cu Updates
Cu   24 May 08 Only read/write master (MPI)
Cu             Note: it is the caller's responsibility to broadcast
Cu             This is because neul and the dimensions of eula
Cu             are only known after reading.  Typical Broadcast:
Cu             call mpibc1(neul,1,2,mlog,'susite','neul')
Cu             call mpibc1(xsi,5,4,mlog,'susite','xsi')
Cu             call mpibc1(eula,nbas*neul*3,4,mlog,'susite','eula')
Cu   19 Nov 97 use rdfiln for reading first line
C ----------------------------------------------------------------
      use m_pvioeu
C     implicit none
      integer nbas,neul,nl,ifi
      double precision eula(nbas,neul,3),tmp(125),xsi(3)
      integer i,j,k,ipr,lgunit,nc,rdm,nxsi,a2vec,ix(10)
      logical parstr,a2bin
      character*80 ss,ss2*4
      integer master,mpipid
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),
     .vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .nlist(0:mxlev)
      equivalence (a,ss)
      data rdarg /'#{}% c'/


C --- I/O through master node only (MPI) ---
      master = 0
      if (mpipid(1) .eq. master) then

C --- Write ---
        if (ifi .lt. 0) then
C       rewind (-ifi)
          if (neul .ne. 1 .and. neul .ne. nl .and. neul .ne. nl*nl)
     .    call rxi('ioeula: bad dim for eula: neula=',neul)
          call awrit3('%% rows %i cols %i neula %i',ss,80,0,nbas*neul,3,
     .    neul)
          nxsi = 0
          do  11  i = 1, 3
            if (xsi(i) .eq. 0) goto 12
            nxsi = nxsi+1
   11     continue
   12     continue
          if (nxsi .gt. 0)
     .    call awrit3('%a nxsi=%i %n:1;6,6d',ss,80,0,nxsi,nxsi,xsi)
          call awrit0('%a',ss,80,ifi)
          do  10  i = 1, nbas
            if (neul .gt. 1) call awrit1('# ib %,4i',' ',80,-ifi,i)
            write(-ifi,'(3f16.12)') ((eula(i,j,k),k=1,3),j=1,neul)
   10     continue
C --- Read ---
        else
          call getpr(ipr)
          rewind ifi
          nr = 0
          call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,ss,recrd,recl,nr)
C   ... No error if file is empty, but don't bother reading
          if (nr .eq. 0) then
            if (ipr .gt. 30) print '('' IOEULA:   empty file'')'
            return
          endif
          if (ss(1:1) .eq. '%') then
C   ... Read xsi if it is there
            i = 0
            if (parstr(ss,'nxsi=',len(ss)-5,5,' ',i,j)) then
              j = j-1
              i = j
              if (a2bin(ss,nxsi,2,0,' ',j,len(ss)-5)) then
                j = a2vec(ss,len(ss)-5,i,4,' #',2,2,nxsi+1,ix,tmp)
                if (j .ne. nxsi+1) call rx('ioeula:  failed to read xsi')
                call dcopy(nxsi,tmp(2),1,xsi,1)
              endif
            endif
            neul = 1
            i = 0
            if (parstr(ss,'neula ',len(ss)-6,6,' ',i,j)) then
              if (.not. a2bin(ss,neul,2,0,' ',j,len(ss)))
     .        call rx('ioeula: failed to parse for neula')
            endif
          else
            neul = 1
            rewind ifi
          endif

C   ... Expect nbas*neul*3 rows in file ... skip reading if not
          if (neul .ne. 1 .and. neul .ne. nl .and. neul .ne. nl*nl)
     .    call rxi('ioeula: illegal dimension: neula =',neul)
          nc = 3
          ss2 = 'in'
          if (neul .eq. nl) ss2 = ' '
          if (neul .eq. nl*nl) ss2 = '%bm-'
          if (ipr .gt. 30) call awrit1(' IOEULA:  reading l-'//ss2//
     .    '%adependent Euler angles (file %i)',' ',80,lgunit(1),ifi)
C        if (nr .ne. neul*nbas .or. nc .ne. 3) then
C          if (ipr .ge. 20) call awrit3(' ioeula: input skipped --'//
C     .      ' expected %i rows (3 cols) but found %i (%i cols)',
C     .      ' ',80,lgunit(1),neul*nbas,nr,nc)

C   ... Read the data, 1 atom at a time ...
          call dpzero(eula,nbas*neul*3)
          do  20  i = 1, nbas
            call rxx(rdm(ifi,0,neul*nc,' ',tmp,neul,nc).ne.1,
     .      'ioeula:  bad file')
            call dmcpy(tmp,neul,1,eula(i,1,1),nbas*neul,nbas,neul,nc)
   20     continue
          if (ipr .ge. 50) call pvioeu(mode=1,eula=eula,nbas=nbas,neul=neul)
        endif
      endif
      end subroutine ioeula
