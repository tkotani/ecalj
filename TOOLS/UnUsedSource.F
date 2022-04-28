      subroutine suclst(nsitmx,nbas,nsp,ssite,sspec,clsopt,isite,iclsl,iclsn,nsites)
      use m_struc_def
      use m_ext,only:sname
C- Set up list of sites for core level spectra
C ----------------------------------------------------------------------
Ci Inputs: nsitmx (dimension of isite,iclsl,iclsn: set in bndfp)
Ci         nbas,ssite,sspec; clsopt (command line options for parsing)
Co Outputs:
Co         isite,iclsl,iclsn (list of site,core-l and core-n)
Co         nsites (total number of site triples returned)
Cr Remarks
Cr parses the string
Cr  --cls[.ib,l,n[.ib,l,n][.lst=list,l,n[.lst=list,l,n]][.fn=file]]
Cr       where '.' is any separator
Cr                 except a ':' or a ',' if lst= is used
Cr                 except a ',' if ib,l,n is used
Cr       ib,l,n is a triple of site,core-l,core-n
Cr       list is a list of sites having l and n in common
Cr       file is a file of nsites records of triples
Cr example
Cr       --cls/2,1,2/4,2,3/lst=3,5:7,0,1
Cr  produces
Cr         site  l  n
Cr           2   1  2
Cr           4   2  3
Cr           3   0  1
Cr           5   0  1
Cr           6   0  1
Cr           7   0  1
Cu Updates
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nsitmx
      integer isite(nsitmx),iclsl(nsitmx),iclsn(nsitmx),nsites,nbas,nsp
      character*(*) clsopt
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C Local Variables
      integer i,j,j1,j2,iprmin,nblst,ib,is,l,n,ifi,p,q
      integer iprint,lgunit,igetss
      character clabl*8,fn*16,dc*1,outstr*64,cc2*2
      logical a2bin
      character(10)::i2char
      iprmin = 10

      dc = clsopt(1:1)
      if (dc .eq. ' ') call rx('suclst: found --cls but no site')
      nsites = 0
      j2 = 0
      fn =  ' '
C   ... Return here to resume parsing for arguments
    1 continue
      j2 = j2 + 1
      if (clsopt(j2:j2) .eq. dc) goto 1
      j1 = min(len(clsopt),j2)
      call nwordg(clsopt,0,dc//' ',1,j1,j2)
      if (j2 .ge. j1) then
        if (clsopt(j1:j1+2) .eq. 'fn=')  then
          if (j1+3 .gt. j2) call rx('suclst: bad file name')
          fn = clsopt(j1+3:j2)
C   ... read in data from file
          open(newunit=ifi,file=trim(fn)//'.'//trim(sname)) 
          i = 0
    2     continue
          read (ifi,*,end=3,err=3) isite(i+1),iclsl(i+1),iclsn(i+1)
          i = i + 1
          if (i .gt. nsitmx) call rxi('bndfp needs nsitmx',nsites)
          goto 2
    3     continue
          close(ifi)
          nsites = i
          goto 4
        elseif (clsopt(j1:j1+3) .eq. 'lst=') then
          if (j1+4 .gt. j2) call rx('suclst: bad list')
C   ... expect ,l,n at end of word; move pointer back to end of list
          j = j2 - 4
          call mkils0(clsopt(j1+4:j),nblst,i)
          call mkilst(clsopt(j1+4:j),nblst,isite(1+nsites))
          j = j + 1
          if (clsopt(j:j) .ne. ',') goto 5
          i = 0
          j = j + 1
          cc2 = clsopt(j:j)//' '
          if (.not. a2bin(cc2,l,2,0,' ',i,-1)) goto 5
          j = j + 1
          if (clsopt(j:j) .ne. ',') goto 5
          i = 0
          j = j + 1
          cc2 = clsopt(j:j)//' '
          if (.not. a2bin(cc2,n,2,0,' ',i,-1)) goto 5
          if (j .ne. j2) call rx('bug in suclst')
          nsites = nsites + nblst
          if (nsites.gt.nsitmx) call rxi('bndfp needs nsitmx',nsites)
          call icopy(nblst,l,0,iclsl(nsites+1-nblst),1)
          call icopy(nblst,n,0,iclsn(nsites+1-nblst),1)
        else
          j = j1
          i = 0
          cc2 = clsopt(j:j)//' '
          if (.not. a2bin(cc2,ib,2,0,' ',i,-1)) goto 6
          j = j + 1
          if (clsopt(j:j) .ne. ',') goto 6
          j = j + 1
          i = 0
          cc2 = clsopt(j:j)//' '
          if (.not. a2bin(cc2,l,2,0,' ',i,-1)) goto 6
          j = j + 1
          if (clsopt(j:j) .ne. ',') goto 6
          j = j + 1
          i = 0
          cc2 = clsopt(j:j)//' '
          if (.not. a2bin(cc2,n,2,0,' ',i,-1)) goto 6
          if (j .ne. j2) call rx('bug in suclst')
          nsites = nsites + 1
          if (nsites.gt.nsitmx) call rxi('bndfp needs nsitmx',nsites)
          isite(nsites) = ib
          iclsl(nsites) = l
          iclsn(nsites) = n
        endif
        goto 1
      endif
    4 continue
      if (nsites .eq. 0) call rx('suclst: found --cls but no site')
C --- print table of sites, l and n
      if (nsites .gt. 0 .and. iprint() .ge. iprmin) then
        if (nsp .eq. 1) then
          write(lgunit(1),"(a,/,a)") ' suclst: set up table of site, l and n for CLS',
     .    '   site   label     l  n   DOS channels'
        else
          write(lgunit(1),"(a,/,a)") ' suclst: set up table of site, l and n for CLS',
     .    '   site   label     l  n   DOS channels up(down)'
        endif
        j = 1
        do  i = 1, nsites
          ib = isite(i)
          if (ib .gt. nbas) call rx('suclst: ib > nbas')
          is = int(ssite(ib)%spec)
          clabl=sspec(is)%name
          if (nsp .eq. 1) then
c            call awrit6('%:3,3i     '//clabl//'%,2i%:1,2i%3f%i,%i,%i',
c     .      ' ',128,lgunit(1),isite(i),iclsl(i),iclsn(i),j,j+1,j+2)
             write(lgunit(1),"(a,2i3,2x,3i5)")
     &            i2char(isite(i))//'     '//trim(clabl),iclsl(i),iclsn(i),j,j+1,j+2
          else
c            call awrit3('%:3,3i     '//clabl//'%,2i%:1,2i',
c     .      outstr,128,0,isite(i),iclsl(i),iclsn(i))
             write(lgunit(1),"(a,2i3)")
     &            i2char(isite(i))//'     '//trim(clabl),iclsl(i),iclsn(i)
            call strip(outstr,p,q)
c            call awrit6(outstr(1:q)//'%3f%i(%i),%i(%i),%i(%i)',' ',128,
c     .      lgunit(1),j,j+1,j+2,j+3,j+4,j+5)
            write(lgunit(1),"(a,6i4)") outstr(1:q),j,j+1,j+2,j+3,j+4,j+5
          endif
          j = j + 3*nsp
        enddo
      endif
      return
    5 continue
C --- bad list ERROR
      call fexit0(-1,'suclst: bad list, expecting --cls'//dc//'lst=list,l,n')
    6 continue
C --- bad triple ERROR
      call fexit0(-1,'suclst: bad triple, expecting --cls'//dc//'site,l,n')
      end subroutine suclst


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
  subroutine sumlst(lopt,nchmx,nbas,ng,ssite,sspec,sopts,mode,
     .nsites,lsites,lmxch,nchan,lchan,nll)
      use m_struc_def  !Cgetarg
C- Set up list of sites and mode for Mulliken or partial DOS projection
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lopt  :1s digit: choice of l channels
Ci         :  0 for partial DOS (lcutoff to lmxa; see Remarks)
Ci         :  1 for Mulliken analysis (lcutoff to lmxb; see Remarks)
Ci         :10s digit what to return
Ci         :  0 return all outputs
Ci         :  1 return nsites and lmxch only
Ci              (NB: if sopts contains site list, lsites also returned)
Ci              (NB: if sopts contains nl=, value of nll is altered)
Ci   nchmx :(lopt=1 only) dimensions array lchan
Ci   nbas  :size of basis
Ci   ng    :# of space group ops --- checks the m-resolved case
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxb lmxa
Ci     Stored:    *
Ci     Passed to: *
Ci   sopts  from command line (see Remarks) with
Ci           the --mull or --pdos stripped
Co Outputs:
Co   mode  :(see remarks)
Co   nsites:(number of sites)
Co   lsites:(list of sites)
Co   nchan :(number of DOS channels)
Co   lmxch :(maximum basis l for all of the atoms in the list)
Co   lchan :(array of channel numbers for each site,l,m) see Remarks
Co         :(lopt=1 only; not used for lopt=0)
Co   nll   :if nonzero, use lmax = nll-1 for all sites in making
Co         :dos channels.  NB: input value of nll passed by caller
Co         :is its default value.
Cr Remarks
Cr   sumlst:options passed in sopts
Cr          [:mode=#][nl=#][:sites=list]]
Cr   modes are 0: all sites atom-resolved (default --mull)
Cr             1: all sites l-resolved
Cr             2: all sites lm-resolved (default --pdos)
Cr             Adding 3 to mode flags sumlst that a site-list is
Cr             to follow.  This is not needed, however, since supplying
Cr             a site-list automatically will cause 3 to be added.
Cr
Cr   The --mull (lopt=1) and --pdos (lopt=0) modes are similar;
Cr   however lopt=1 case is complicated by the possibility of multiple
Cr   basis functions per l channel, and also different lmxb in different
Cr   kappa channels.  To handle that case, sumlst creates a list marking
Cr   which DOS channel each site, l and m belongs.  sumlst returns
Cr   this array in the form that mullmf can reference as
Cr   lchan(ib) lchan(l,ib) or lchan(lm,ib) depending on the mode
Cr   and lchan will be dimensioned using the maximum l, lmxch:
Cr   lchan(ib) lchan(lmxch+1,ib) or lchan((lmxch+1)**2,ib)
Cr   depending on the mode. The ordering of ib is either as the
Cr   internal list of sites (mode .lt. 3) or as the list returned
Cr   from the cmd-line string by mkilst. (see mchan.f)
Cu Updates
Cu   30 Aug 04 Changed first argument to integer lopt
Cu   18 Jan 02 (MvS) redesigned switches for compatibility with site
Cu                   lists in other context
Cu   27 Mar 01 (MvS) extended to handle both --mull and --pdos
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      character*(*) sopts
      integer lopt,nchmx,nbas,ng,nsites,mode,lsites(*),nchan,lmxch,
     .lchan(*),nll
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C Local Variables
      integer mxmode
      parameter (mxmode = 6)
      character dc*1,modstr(0:mxmode-1)*128
      integer lmax,isite,lmdim,ltop,ilm
      integer i,j,ipr,iprmin,ib,is,iprint,lgunit,stdo,igetss,isw,j2,j1,
     .m,iv(10),parg,lstyle,lopt0,lopt1
      double precision z
      data modstr /'all sites atom-projected',
     .'all sites l-projected','all sites lm-projected',
     .'site list atom-projected',
     .'site list l-projected','site list lm-projected' /
      logical:: isanrg, l_dummy_isanrg
      character(256)::ttt
      lopt0 = mod(lopt,10)
      lopt1 = mod(lopt/10,10)
      ipr = iprint()
      iprmin = 30
      mode = -1
      dc = sopts(1:1)
      stdo = lgunit(1)
      mode = 2
      nsites = 0
      if (lopt0 .eq. 1) mode = 0
      if (dc .ne. ' ') then
C   ... Return here to resume parsing for arguments
        j2 = 0
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) .eq. dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 .ge. j1) then
          if (.false.) then
C         DOS mode; see Remarks
          elseif (sopts(j1:j1+4) .eq. 'mode=')  then
            m = 0
            i = parg('mode=',2,sopts(j1:),m,len(sopts(j1:)),
     .      dc//' ',1,1,iv,mode)
            if (i .le. 0) goto 999
Ckino isanrg is logical function,             call isanrg(mode,0,5,' sumlst:','mode',.true.)
            l_dummy_isanrg=isanrg(mode,0,5,' sumlst:','mode',.true.)
            if (ipr .ge. 10) then
              if ((mode .eq. 2 .or. mode .eq. 5) .and. ng .gt. 1)
     .        write(stdo,*)' **WARNING** sumlst: for lm-decomposition suppress symops'
            endif
C         Number of energy points.  Not used now
          elseif (sopts(j1:j1+2) .eq. 'nl=')  then
            m = 0
            i = parg('nl=',2,sopts(j1:),m,len(sopts(j1:)),
     .      dc//' ',1,1,iv,nll)
            if (i .le. 0) goto 999
C         Site list
          elseif (sopts(j1:j1+5) .eq. 'sites=')  then
            if (mode .lt. 3) mode = mode+3
C           We need to make slabl,z before using other lstyle
            lstyle = 1
            call slist(lstyle,sopts(j1+6:j2),' ',z,nbas,nsites,lsites)
          else
            goto 999
          endif
          goto 10
        endif
      endif
      if (ipr .ge. iprmin+1) write(stdo,1)
    1 format (' sumlst: Site   l_max')
      if (mode .lt. 3) nsites = nbas

      if (nsites .eq. 0) call rxi
     .('sumlst: --sites=list is required input for mode =',mode)

C --- Generate lsites,lmxch ---
      lmxch = -1
      do  isite = 1, nsites
        if (mode .lt. 3) then
          ib = isite
          if (lopt1 .eq. 0) then
            lsites(isite) = ib
          endif
        else
          ib = lsites(isite)
        endif
Ckino isanrg is logical function,         call isanrg(ib,1,nbas,' sumlst:','site',.true.)
        l_dummy_isanrg=isanrg(ib,1,nbas,' sumlst:','site',.true.)

        is=ssite(ib)%spec

        if ( lopt0 .eq. 1 ) lmax = int(sspec(is)%lmxb)

        if ( lopt0 .eq. 0 ) lmax = int(sspec(is)%lmxa)

        if (nll .gt. 0) lmax = min(nll-1,lmax)
        lmxch = max(lmxch,lmax)
      enddo
      if (lopt1 .eq. 1) return

C --- Generate nchan,lchan ---
      i = 0
      do  isite = 1, nsites
        ib = lsites(isite)

        is=ssite(ib)%spec

        if ( lopt0 .eq. 1 ) lmax = int(sspec(is)%lmxb)

        if ( lopt0 .eq. 0 ) lmax = int(sspec(is)%lmxa)

        if (nll .gt. 0) lmax = min(nll-1,lmax)
        if (ipr .gt. iprmin) write(stdo,2) ib, lmax
    2   format (9x,i4,5x,i1)
        if (mode .eq. 0 .or. mode .eq. 3) then
          ltop = 1
          lmdim = 1
        elseif (mode .eq. 1 .or. mode .eq. 4) then
          ltop = lmax+1
          lmdim = lmxch+1
        elseif (mode .eq. 2 .or. mode .eq. 5) then
          ltop = (lmax+1)**2
          lmdim = (lmxch+1)**2
        endif
        if (lopt0 .eq. 1) then
          do  ilm = 1, ltop
            i = i + 1
            if (i.gt.nchmx) call rxi(' sumlst: bndfp needs nchmx > ',i)
            call mchan(lmdim,0d0,0d0,0,nsites,0,isite,ilm,-1,i,lchan)
          enddo
        else
          i = i+ltop
        endif
      enddo
      nchan = i
C --- Printout --
      if (ipr .ge. iprmin) then
         call strip(modstr(mode),i,j)
         if(lopt0==0) ttt=' sumlst:  Mulliken mode '
         if(lopt0==1) ttt=' sumlst:  Partial DOS mode '
         write(stdo,"(a,i5,a,i5,' sites ',i5,' channels')")ttt,mode,modstr(mode)(i:j),nsites,nchan
      endif
      return
  999 continue
      call rxs('sumlst: failed to parse options in ',sopts(j1:))
      end subroutine sumlst


      logical function pars1v(strn,size,nam,term,cast,val)
C- Parses a string for a specified variable
C ----------------------------------------------------------------
Ci Inputs
Ci   recrd(0:*): string recrd is parsed from i to size-1
Ci   iopt:   0: var declaration of existing variable is supressed
Ci           1: var declaration of existing variable supersedes
Ci   mxdecl: parsing ends if number of declarations exceeds mxdecl
Co Outputs
Co   i:      last character parsed
Cr Remarks
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer size,cast
      character*(*) strn,nam,term*1
      double precision val
C Local parameters
      integer i,j
      logical parstr,a2bin

      pars1v = .false.
      i = 0
      if (.not. parstr(strn,nam,size,len(nam),term,i,j)) return
      if (a2bin(strn,val,cast,0,' ',j,size)) pars1v = .true.
      end
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

      
#define USEOPTIONAL

      module m_ioden2
      contains
      subroutine ioden2(mode,sopts,nsp,plat,alat,nbas,pos,z,k1,k2,k3,
     .ngabc,smrho,wk,nwk,modrhg)

C- Kernel called by ioden
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 just return these arguments parsed from sopts:
Ci         :  modrhg
Ci         :1 copy and save density in appropriate plane
Ci   sopts :string containing options; see ioden above.
Ci   plat  :primitive lattice vectors, in units of alat (for printout)
Ci   alat  :lattice constant
Ci   nbas  :number of atoms in unit cell
Ci   pos   :atomic positions in unit cell
Ci   z     : atomic number of atoms in cell
Ci   nsp   :number of spin channels
Ci   k1..3 :dimensions smrho
Ci   ngabc :number of divisions in smrho along plat(1..3)
Ci   smrho :smooth density on uniform mesh
Ci   wk    :2D work array
Ci   nwk   :dimensions wk
Co Outputs
Co   modrhg:controls what part of core density is added
Co   The projection of smrho onto a plane is written to disk
Co   or if option is g3d the full 3d grid is written to disk
Cl Local variables
Cl   lic   :lic(1..3,1..2) = step length in each vector axes 1..2
Cl         :lic(4,1..2)    = number of points along axes 1..2
Cl         :lic(5..7,1)    = starting element for each lattice vector
Cr Remarks
Cr   Choosing p1, p2, p3 sets the origin for that plat and also
Cr   creates a default set of vectors, choosing the other lattice
Cr   vectors in the system.  Thus:
Cr   p1=# => lic(1..7,1) = (0 1 0 nb 1 1 1) lic(1..7,2) = (0 0 1 nc **)
Cr   p2=# => lic(1..7,1) = (0 0 1 nc 1 1 1) lic(1..7,2) = (1 0 0 na **)
Cr   p3=# => lic(1..7,1) = (1 0 0 na 1 1 1) lic(1..7,2) = (0 1 0 nb **)
Cu Updates
Cu   09 Feb 01 Revised how a plane is defined.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) sopts
      integer mode,ngabc(3),k1,k2,k3,nwk,nsp,modrhg,nbas
#ifdef USEOPTIONAL
      double complex ,optional:: smrho(k1,k2,k3,nsp)
      double precision plat(3,3),alat,pos(3,nbas)
      real(8),optional:: wk(nwk,nwk,nsp)
#else
      double complex smrho(k1,k2,k3,nsp)
      double precision plat(3,3),wk(nwk,nwk,nsp),alat,pos(3,nbas)
#endif
      double precision z(nbas)
C ... Local parameters
      logical lok
      integer lgunit,stdo,j2,j1,fopn,ifi,i
      integer lic(7,2),j,k,iv(3),jv(3),kv(4),a2vec,isp
      double precision vecs(3,2),vv,ddot
      character*120 dc*1, fn, prfmt*40, sout
      logical lg3d
      integer i1,i2,i3

      logical:: l_dummy_isanrg,isanrg

ccccccccccccccccccccccccccccc
c      print *,' zzxxx=',z
c      stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzz'
ccccccccccccccccccccccccccccc

      stdo = lgunit(1)
      dc = sopts(1:1)
      fn = 'smrho'
C     prfmt = '(8f10.6)'
      prfmt = '(8e14.6)'
      lg3d = .false.
C     Default: origin at (1,1,1), (na+1,nb+1) points along 1st,2nd axes
      call iinit(lic,14)
      call ivset(lic,5,7,1)
      lic(1,1) = 1
      lic(4,1) = ngabc(1)+1
      lic(2,2) = 1
      lic(4,2) = ngabc(2)+1

      if (dc .ne. ' ') then
C ... Return here to resume parsing for arguments
        j2 = 0
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) .eq. dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 .ge. j1) then
          if (.false.) then
          elseif (sopts(j1:j1+2) .eq. 'fn=')  then
            if (j1+3 .gt. j2) goto 99
            fn = sopts(j1+3:j2)
C   ... option g3d
          elseif (sopts(j1:j1+2).eq. 'g3d') then
            lg3d = .true.
C   ... option o=#,#,#
          elseif (sopts(j1:j1+1) .eq. 'o=')  then
            if (j1+2 .gt. j2) goto 99
            i = j1+1
            if (a2vec(sopts,j2,i,2,', '//dc,3,2,3,kv,iv) .ne. 3)
     .      goto 99
            lic(5,1) = mod(iv(1)+ngabc(1),ngabc(1))+1
            lic(6,1) = mod(iv(2)+ngabc(2),ngabc(2))+1
            lic(7,1) = mod(iv(3)+ngabc(3),ngabc(3))+1
C   ... option l[12]=#,#,#,#
          elseif (sopts(j1:j1+2) .eq. 'l1=' .or.
     .    sopts(j1:j1+2) .eq. 'l2=')  then
            if (j1+3 .gt. j2) goto 99
            i = 0
            call chrps2(sopts(j1+1:),'12',2,0,i,iv)
C         this check should never be necessary
Ckino isanrg is logical function, C         call isanrg(iv,1,2,' ','iv in ioden ... bug ..',.true.)
C         l_dummy_isanrg=isanrg(iv,1,2,' ','iv in ioden ... bug ..',.true.)
            i = j1+2
            j = iv(1)
            k = a2vec(sopts,j2,i,2,', '//dc,3,2,4,kv,lic(1,j))
            if (k .ne. 3 .and. k .ne. 4) goto 99
          elseif (sopts(j1:j1+4) .eq. 'core=') then
            if (j1+5 .gt. j2) goto 99
            i = j1+4
            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,modrhg) .ne. 1)
     .      goto 99
          else
            call rxs('ioden: unrecognised option ... ',sopts(j1:j2))
          endif
          goto 10
        endif
      endif
      if (mode .eq. 0) return

C --- Printout ---
      if (lg3d) then
        ifi = fopn(fn)
        rewind ifi
        call awrit3('      Writing smooth density to file'//fn//
     .  ' %a : full 3d grid (%i,%i,%i).',' ',80,stdo,k1,k2,k3)
        do isp = 1, nsp
C           call awrit4('%% 3D block data %i %i %i, spin=%i',' ',80,ifi,
C     .     k1,k2,k3,isp)
          write(ifi,'("CRYSTAL")')
          write(ifi,'("PRIMVEC")')
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
     .    ,i2=1,3)
          write(ifi,'("PRIMCOORD")')
          write(ifi,'(2i5)') nbas,1
          do i = 1, nbas
            write(ifi,'(i4,2x,3f10.5)') int(z(i)),
     .      (pos(i2,i)*alat*0.529177208,i2=1,3)
          enddo
          write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
          write(ifi,'("charge_density_spin_",i1)') isp
          write(ifi,'("BEGIN_DATAGRID_3D_isp_",i1)') isp
          write(ifi,'(3i4)') k1,k2,k3
          write(ifi,'(3f10.5)') 0.,0.,0.
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
     .    ,i2=1,3)
          write(ifi,'(8e14.6)')  
     .    (((dble(smrho(i1,i2,i3,isp)),i1=1,k1),i2=1,k2),i3=1,k3)
          write(ifi,'("END_DATAGRID_3D_isp_",i1)') isp
          write(ifi,'("END_BLOCK_DATAGRID_3D")')
        enddo
        call fclose(ifi)
        return
      else
        call awrit3('         Writing smooth density to file '//fn//
     .  '%a : origin at (%i,%i,%i).',' ',80,stdo,lic(5,1)-1,
     .  lic(6,1)-1,lic(7,1)-1)
        call dpzero(vecs,3)
        do  i = 1, 3
          call daxpy(3,dble(lic(4+i,1)-1)/ngabc(i),plat(1,i),1,vecs,1)
        enddo
        call ivset(kv,1,3,2)
        call shorps(1,plat,kv,vecs,vecs)
        write(stdo,368) (vecs(i,1), i=1,3)
  368   format(9x,'In cartesian coordinates, origin is',3f10.6)

        call dpzero(vecs,6)
        do  k = 1, 2
          do  i = 1, 3
            vv = dble((lic(4,k)-1)*lic(i,k))/dble(ngabc(i))
            call dpadd(vecs(1,k),plat(1,i),1,3,vv)
          enddo
          vv = dsqrt(ddot(3,vecs(1,k),1,vecs(1,k),1))
          sout = ' '
          call awrit8('%9fv%i: (%i pts) = (%i,%i,%i)p -> '//
     .    '(%,6;6d,%,6;6d,%,6;6d)',sout,80,0,k,lic(4,k),
     .    lic(1,k),lic(2,k),lic(3,k),vecs(1,k),vecs(2,k),vecs(3,k))
          call awrit1('%a l=%,6;6d',sout,-80,-stdo,vv)
        enddo

C ... Sanity checks
        lok = .true.
        lok = lok .and. lic(4,1) .gt. 0
        lok = lok .and. lic(4,2) .gt. 0
        if (.not. lok) call fexit2(-1,1,' Exit -1 ioden: number of '//
     .  'points along axes (%i,%i) are not > 0',lic(4,1),lic(4,2))
        if (nwk.lt.lic(4,1) .or. nwk.lt.lic(4,2)) call rx('increase nwk')

C --- Copy points to wk from plane ---
        call icopy(3,lic(5,1),1,kv,1)
        do  k = 1, lic(4,2)
          call icopy(3,kv,1,iv,1)
C         PRINT *, K, KV(1), KV(2), KV(3)
          do  i = 1, lic(4,1)
C         PRINT *, ' ',I, IV(1), IV(2), IV(3)
            do  j = 1, 3
              jv(j) = mod(iv(j)-1,ngabc(j))+1
   44         if (jv(j) .le. 0) then
                jv(j) = jv(j) + ngabc(j)
                goto 44
              endif
            enddo
            if (jv(1) .le. 0 .or. jv(2) .le. 0 .or. jv(3) .lt. 0) then
              call rx('bug in ioden')
            endif
            do  isp = 1, nsp
              wk(i,k,isp) = dble(smrho(jv(1),jv(2),jv(3),isp))
            enddo
C          if (i .eq. 1 .or. i .eq. 16) then
C          PRINT *, I,K, JV(1), JV(2), JV(3),WK(I,K)
C          endif
            do  j = 1, 3
              iv(j) = iv(j) + lic(j,1)
            enddo
          enddo
          do  j = 1, 3
            kv(j) = kv(j) + lic(j,2)
          enddo
        enddo

C --- Write the file ---
        ifi = fopn(fn)
        rewind ifi
        do  isp = 1, nsp
          call ywrm(0,' spin 1',1,ifi,prfmt,wk(1,1,isp),0,nwk,lic(4,1),
     .    lic(4,2))
        enddo
        call fclose(ifi)
        return
      endif

   99 continue
      call rxs('ioden: failed to parse option ... ',sopts(j1:j2))

      end subroutine ioden2
      end module m_ioden2

      subroutine ioden ( sopts , slat , ssite , sspec , sv_p_orhat 
     . , smrho )


      use m_struc_def  !Cgetarg
      use m_globalvariables
      use m_ioden2
C- File I/O charge density on a uniform mesh in a plane or full 3d mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :character string specifying plane and other options.
Ci         :Specifications and options are separated by a delimiter,
Ci         :which is the first character in sopts.
Ci         :
Ci         :if option is g3d then the density is written to
Ci         :disk on the full 3d grid.   Otherwise:
Ci         :The density is written to disk for a uniform of mesh of
Ci         :points in one plane.  This information is specified by three
Ci         :groups of numbers: the origin, a first direction vector with
Ci         :its number of points, and a second direction vector with its
Ci         :number of points.
Ci         :
Ci         :At present, these points must be taken from the points on
Ci         :the smooth mesh, smrho.  In this case, all three groups of
Ci         :information are sets of integers.  For example,
Ci         :specify the origin by three numbers:
Ci         :    o=#1,#2,#3
Ci         :The point (#1,#2,#3) corresponds to the Cartesian coordinate
Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
Ci         :where (n1,n2,n3) are the number of divisions in the
Ci         :mesh along the three lattice vectors (p1,p2,p3).
Ci         :o=0,0,0 corresponds to the origin.
Ci         :
Ci         :Specify the direction vectors by
Ci         :    l1=#1,#2,#3[,#4]
Ci         :    l2=#1,#2,#3[,#4]
Ci         :
Ci         :l1 and l2 specify the first and second direction vectors,
Ci         :respectively.  #1,#2,#3 select the
Ci         :increments in mesh points along each of the three lattice
Ci         :vectors that define the direction vector.  Thus in Cartesian
Ci         :coordinates a direction vector is
Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
Ci         :where as before (n1,n2,n3) are the number of divisions in
Ci         :the mesh along the three lattice vectors (p1,p2,p3).
Ci         :The last number (#4) specifies how many points to take
Ci         :in that direction.
Ci         :
Ci         :Other options:
Ci         :  g3d       output of 3d grid plus headers in xsf format ready for xcrysden
Ci         :  fn=name   specifies the file name for file I/O
Ci         :
Ci         :   core=#   specifies how local rho is to be included
Ci         :            #=0 include core densities - nuclear charge
Ci         :            #=1 include core densities
Ci         :            #=2 (default) exclude core densities
Ci         :            #=-1 no local densities to be included
Ci         :            #=-2 true local density, no smoothed part
Ci         :            #=-3 istl-local sm densities, no true part
Ci         :
Ci         :Example: using '~' as the delimiter, and suppose
Ci         :n1=n2=48 and n3=120, the specification
Ci         :  ~fn=myrho~o=0,0,60~l1=1,1,0,49~l2=0,0,1,121
Ci         :writes 'myrho.ext' a mesh (49,121) points.
Ci         :The origin (first point) lies at (p3/2)
Ci         :The first vector points along (p1+p2), and has that length;
Ci         :the second vector points along p3, and has that length.
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat nabc ng okv vol
Ci     Stored:    *
Ci     Passed to: rhgcmp
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: rhgcmp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: rhgcmp
Ci   orhat :offsets to arrays containing local atom densities
Ci   smrho :smooth density on uniform mesh
Co Outputs
Co   The sum local gaussian densities and smrho is written to disk
Cl Local variables
Cl   modrhg:controls what part of core density is added
Cl         :-1 Add no local densities
Cl         : 0 include core density - sm nuclear charge
Cl         : 1 include core density, but no nuclear charge
Cl         : 2 exclude core density
Cr Remarks
Cr    sopts specifies which plane(s) are written to disk
Cr Bugs
Cr   Routines create smoothe approximation to density, not true density
Cu Updates
Cu   July 7 05 (Walter Lambrecht) New option g3d
Cu   25 Aug 04 New modes -2, -3
Cu   24 May 03 Corrections to errors in messages
Cu   23 Oct 01 Local densities are more accurately represented
Cu             in G_kL expansion:  k=0..kmax.  Added core= option.
Cu   25 Apr 01 Simplified the plane specification
Cu   02 Mar 01 Spin polarized
Cu   09 Feb 01 Added local gaussian densities to file I/O
Cu             and extended how a plane may be specified
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) sopts
      type(s_lat)::slat
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      type(s_rv1) :: sv_p_orhat

      double complex smrho(*)
C ... Local parameters
      integer ngabc(3),n1,n2,n3,k1,k2,k3,kmax
      integer kkk,lgunit,ng,nglob,nsp,nwk,stdo,modrhg,nbas,i
ckino Dec.9.2011:          real(8),pointer :: rv_p_opos(:) =>NULL()

ckino Dec.9.2011:         integer,pointer :: iv_p_okv(:) =>NULL()

      complex(8) ,allocatable :: cn_zv(:)
      complex(8) ,allocatable :: psrho_zv(:)
ckino not complex but real      complex(8) ,allocatable :: wk_zv(:)
      real(8) ,allocatable :: wk_rv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision vol,xx0,xx,dval,plat(3,3),alat
      integer ib,is,igetss
C ... Heap
      integer ::iwdummy

      logical:: l_dummy_isanrg,isanrg

C ... External calls
ckino del upack,upack2
      external awrit2,dcopy,fclose,fftz3,fftz30,gvgetf,
     .gvputf,icopy,ivset,mkils0,mkilst,nwordg,poppr,
     .pshpr,rhgcmp,rhomom,rx

      real(8),allocatable:: zz(:)
C ... Unpack and setup
Changenglob      nsp = nglob('nsp')
      nsp = globalvariables%nsp
Changenglob      nbas = nglob('nbas')
      nbas = globalvariables%nbas
      stdo = lgunit(1)

      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      ng=slat%ng
ckino Dec.9.2011:         iv_p_okv => slat%iv_p_okv

      vol=slat%vol

ccccccccccccccccccccccccccccc
c      allocate(zz(nbas))
c      zz=sspec(1:nbas)%z
c      print *,' zz1=',zz
c      print *,' zz2=',sspec(1:nbas)%z
c      stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzz'
ccccccccccccccccccccccccccccc


      alat=slat%alat
ckino Dec.9.2011:         rv_p_opos => slat%rv_p_opos


c      call defrr(oz,nbas)
c      do ib = 1,nbas
c        is = int(ssite(ib)%spec)
cc        w(oz-2+ib*2)=sspec(is)%z !but Kino's email on 5June2010.
c        call dcopy(1,sspec(is)%z,1,w(oz-2+ib*2),1)
c      enddo
      call fftz30(n1,n2,n3,k1,k2,k3)
      kkk = k1*k2*k3
      kmax = 3
      modrhg = 2
c      call ioden2(0,sopts,nsp,plat,alat,nbas,w(opos),w(oz),k1,k2,k3,
c     .ngabc,w,w,nwk,modrhg)
Ckino  mode=0 checks format
#ifdef USEOPTIONAL
      call ioden2 ( mode=0 , sopts=sopts ,nsp= nsp ,plat= plat ,
ckino Dec.9.2011:       .   alat=alat , nbas=nbas , pos=rv_p_opos
ckino Jan.04.2012:       .   alat=alat , nbas=nbas , pos=slat%rv_p_opos 
     .   alat=alat , nbas=nbas , pos=slat%rv_a_opos 
     .  , z=sspec ( ssite ( 1:nbas ) %spec ) %z , k1=k1 , k2=k2 , k3=k3
     .   , ngabc=ngabc 
     .  , nwk=nwk , modrhg=modrhg )
#else
ckino Dec.9.2011:          call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , rv_p_opos
ckino Jan.04.2012:          call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , slat%rv_p_opos
      call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , slat%rv_a_opos 
     .  , sspec ( ssite ( 1:nbas ) %spec ) %z , k1 , k2 , k3 , ngabc 
     .  , iwdummy , iwdummy , nwk , modrhg )
#endif


Ckino isanrg is logical function,       call isanrg(modrhg,-3,2,'ioden','core option',.true.)
      l_dummy_isanrg=isanrg(modrhg,-3,2,'ioden','core option',.true.)

C ... Overwrite smrho+, smrho- with smrho, smrho+ - smrho-
Ckino  smrho is complex(8), dsumdf() is for real(8)
Ckino  so kkk*2 is necessary to use dsumdf().
Ckino
Ckino  return values: smrho(:kkk)=spin(1)+spin(2)
Ckino  and            smrho(kkk+1:kkk*2)=spin(1)-spin(2)
      if (nsp .eq. 2) then
        call dsumdf(kkk*2,1d0,smrho,0,1,smrho(1+kkk),0,1)
      endif

C ... Put n0(G) into psrho and cn
      allocate(psrho_zv(kkk*nsp))

      allocate(cn_zv(ng*nsp))

Ckino smrho is complex(8) , so kkk*2 is necessary.
      call dcopy ( kkk * 2 * nsp , smrho , 1 , psrho_zv , 1 )

      call fftz3 ( psrho_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp 
     ., 0 , - 1 )

ckino Dec.9.2011:         call gvgetf ( ng , nsp , iv_p_okv , k1 , k2 , k3 , psrho_zv ,
ckino Dec.28.2011:         call gvgetf ( ng , nsp , slat%iv_p_okv , k1 , k2 , k3 , psrho_zv ,
      call gvgetf ( ng , nsp , slat%iv_a_okv , k1 , k2 , k3 , psrho_zv , 
     . cn_zv )
Ckino now cn_zv(G,nsp)


Ckino xx0 = real and G=0 component
      xx0 = dval ( cn_zv , 1 )


C ... Add sum of local gaussian densities to mesh density
      print *,' modrhg=',modrhg
Ckino core=#   specifies how local rho is to be included
Ckino         :            #=0 include core densities - nuclear charge
Ckino         :            #=1 include core densities
Ckino         :            #=2 (default) exclude core densities
Ckino         :            #=-1 no local densities to be included
Ckino         :            #=-2 true local density, no smoothed part
Ckino         :            #=-3 istl-local sm densities, no true part
      if (modrhg .ge. 0 .or. modrhg .le. -2) then
        if (modrhg .eq.  0) i = 131
        if (modrhg .eq.  1) i =  31
        if (modrhg .eq.  2) i =   1
        if (modrhg .eq. -2) i =   2
        if (modrhg .eq. -3) i =   3
        if ( i .eq. 2 ) call dscal ( ng * nsp * 2 , 0d0 , cn_zv , 
     .  1 )
Ckino                      -> cn_zv=0 if i==2

        if ( i .eq. 3 ) call dscal ( ng * nsp * 2 , - 1d0 , cn_zv 
     .  , 1 )
Ckino                      -> cn_zv=-cn_zv if i==3

        call rhgcmp ( i , 1 , nbas , ssite , sspec , slat , sv_p_orhat 
     .   , kmax , ng , cn_zv )


        if ( i .eq. 3 ) call dscal ( ng * nsp * 2 , - 1d0 , cn_zv 
     .  , 1 )

      endif
      xx = dval ( cn_zv , 1 )


C ... FFT (n0 + gaussians) (G) to real-space mesh
      if (modrhg .ge. 0) then
ckino Dec.9.2011:           call gvputf ( ng , nsp , iv_p_okv , k1 , k2 , k3 , cn_zv , psrho_zv
ckino Dec.28.2011:           call gvputf ( ng , nsp , slat%iv_p_okv , k1 , k2 , k3 , cn_zv , psrho_zv
        call gvputf ( ng , nsp , slat%iv_a_okv , k1 , k2 , k3 , cn_zv , psrho_zv 
     .   )


        call awrit2('%N ioden : local densities + envelope density,'
     .  //' Qloc=%,6;6d  Q=%,6;6d',' ',80,stdo,(xx-xx0)*vol,xx*vol)
        if (modrhg .lt. 2) then
          call awrit1('%9fLocal densities include core'//
     .    '%?#n==0#+nuclear## contributions.',' ',80,stdo,modrhg)
        endif
      elseif (modrhg .eq. -2) then
        call info2(0,0,0,'%N ioden : local'//
     .  ' densities (true-smooth terms),  Qloc=%,6;6d',xx*vol,0)
      elseif (modrhg .eq. -3) then
        call info2(0,0,0,'%N ioden : smooth (envelope - local)'//
     .  ' density,  Qs=%,6;6d  Qs-Qs(loc)=%,6;6d',
     .  xx0*vol,xx*vol)
      else
        call info2(0,0,0,'%N ioden : smooth density only'//
     .  ' (no local densities added) Qsm=%,6;6d',xx0*vol,0)
      endif

      call fftz3 ( psrho_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp 
     ., 0 , 1 )
Ckino now psrho_zv is real space mesh.


      if (allocated(cn_zv)) deallocate(cn_zv)


C ... File I/O
      nwk = 12*max(k1,k2,k3)

ckino not complex but real,      allocate(wk_zv(nwk**2*nsp))
      allocate(wk_rv(nwk**2*nsp))

Ckino mode=1 write smrho to the file
ckino Dec.9.2011:         call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , rv_p_opos
ckino Jan.04.2012:         call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , slat%rv_p_opos
      call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , slat%rv_a_opos 
     . , sspec ( ssite ( 1:nbas ) %spec ) %z , k1 , k2 , k3 , ngabc
ckino not complex but real     . , psrho_zv , wk_zv , nwk , modrhg ) 
     . , psrho_zv , wk_rv , nwk , modrhg )


ckino not complex but real      if (allocated(wk_zv)) deallocate(wk_zv)
      if (allocated(wk_rv)) deallocate(wk_rv)
      if (allocated(psrho_zv)) deallocate(psrho_zv)

c      call rlse(oz)


C ... Restore smrho+, smrho-
      if (nsp .eq. 2) then
        call dsumdf(kkk*2,.5d0,smrho,0,1,smrho(1+kkk),0,1)
      endif
      end subroutine ioden




cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  #define MULLIKEN 1
      subroutine suqlst(sopts,iop,nband,efermi,nsp,evl,nfbn,ifblst,nq,
     .qp, onespdummy) !onesp removed 
      use m_globalvariables
C- Set up a list of q-points in various modes for energy bands
C  See also entry suqlsw(nband,jsp,nsp,evl)
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :character string describing special options that can
Ci         :selects the mode for which qp are read.
Ci         :*default mode: the qp file consists of a list of
Ci         :               lines and the number of qp in each line.
Ci         :*list mode   : the qp file consists of a list of qp
Ci         :*contour mode: input for specifying qp on a uniform 2D
Ci                       : mesh, for contour plots.
Ci         :See Remarks for the syntax of each mode.
Ci
Ci         :Options are separated by delimiters; the first character is
Ci         :the delimiter.  The following list the strings declaring
Ci         :options assuming the delimiter is '/'. (Space as the first
Ci         :implies that there are no options)
Ci         :/fn=fnam     read qp from file 'fnam'
Ci         :/con         input file for contour plot mode (mode 3)
Ci         :/qp          input file specifies a list of qp (mode 2)
Ci         :             (default assumes a list of symmetry lines.)
Ci         :/long        write bands with extra digits precision
Ci         :             (has no effect for symmetry line mode)
Ci         :/spin1       generate bands only for first spin
Ci         :/lst=list    write only those bands specified in a list.
Ci         :             For syntax of list, see slatsm/mkilst.f
Ci         :/evn=#       keep track of smallest, largest eval for
Ci                       #th band, and printout at close.
Ci         :/ef=#        change efermi to #.
Ci         :Example: --band~long~qp~lst=2:5
Ci   iop   :options passed by the caller
Ci         :1s digit is suppress actions for parallel mode
Ci         : mode 1:
Ci         : 1: suppress writing line header info file
Ci         : 2: Return total number of qp to be generated in all lines
Ci         :    and also suppress writing line header info file
Ci         :    Thus iop=2 => all modes return total no qp to be generated
Ci         : mode 3:
Ci         : nonzero: do not allocate evsav until iq>nq
Ci
Ci   nband :(suqlst) maximum number of energy bands to write
Ci         :(suqlsw) actual number of energy bands to write
Ci
Ci   efermi:Fermi energy (written to bnds file)
Ci
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci         :NB used only in file write to indicate how many times
Ci         :nsp is ALTERED to nsp=1 if spin1 option is set
Ci
Ci   evl   :eigenvalues for current qp (used only in suqlsw)
Cio Inputs/Outputs (see Remarks)
Cio  nq    :Input nq=0 :if zero, flags suqlst to set up bands mode:
Cio                    :some initializations; then returns ready
Cio                    :for first block of qp
Cio        :           :if nonzero, setup for for next qp
Cio        :Output nq  :if zero, suqlst has no more qp to calculate.
Cio                    :if nonzero, nq is the number of k-points in the
Cio                    :current block.  Caller then calls suqlst nq
Cio                    :times (one for each k-point), which will return
Cio                    :qp for each k-point.  After exactly nq calls,
Cio                    :suqlst will start another block, if it exists.
Cio                    :See Remarks for schematic of calling sequence.
Co Outputs
Co   nfbn  :(color weights) number of elements in iblst
Co   iblst :(color weights) list of orbital indices for color weights:
Co         :color weight is the sum of contributions from eigenvector
Co         :components in iblst
Co   qp    :not set, if suqlst is called starting a new block.
Co         :(see Remarks).  Else, k-point at which to generate bands
Co    onesp(removed by takao) :
Co    if spin1 flag is encountered, onesp is set to 1 and nsp is set to 1
Cl Local variables
Cl   iq    :current qp in this block
Cl   mode  :1 symmetry-line mode
Cl         :2 list-of-qp mode
Cl         :3 contour mode
Cl   q1    :starting qp for symmetry mode; only meaningful in that mode
Cl   q2    :ending   qp for symmetry mode; only meaningful in that mode
Cl   nevn  :(optional) band index; routine monitors largest, smallest value
Cl         :for that index.
Cl   evmnn :smallest value found for qp of specified band index
Cl   evmxn :largest  value found for qp of specified band index
Cl   ifiq  :file logical unit for input qp file
Cl   ifib  :file logical unit for output bands file
Cl   ql    :local copy of current qp
Cl   oqp
Cr Remarks
Cr   suqlst is designed to be called to generate qp in groups or
Cr   blocks.  The calling sequence is:
Cr
Cr     nq = 0  <- flags that first call to suqlst, to set up mode
Cr     do  iblock = 1, forever
Cr       This call generates nq, the number of points in this block
Cr       call suqlst(nband,efermi,nsp,evl,nq,qp,onesp)
Cr       if (nq .eq. 0) stop
Cr       do  iq = 1, nq
Cr         This call generates qp for current block
Cr         call suqlst(ndimh,ef0,nsp,w,nkp,qp,onesp) <- returns qp
Cr         do  isp = 1, nsp
Cr         call suqlsw(ndimh,qp,evl(1,isp)) <- saves evl for this qp
Cr                                             (call is optional)
Cr         enddo
Cr       enddo
Cr     enddo
Cr
Cr   The following modes are implemented:
Cr     mode=1 reads qp from syml file, and generates qp along each
Cr            specified symmetry line.  Structure of qp file:
Cr            file has one line for each symmetry line as follows:
Cr               nq      q1x   q1y   q1z      q2x   q2y   q2z
Cr               ...
Cr            the line entries have meanings:
Cr            --# qp-   ---starting qp---    --- ending qp ---
Cr            Any line with nq=0 implies no more lines.
Cr     mode=2 reads qp from specified file and generates qp for each
Cr            specified qp.  File consists sets of qpx,qpy,qpz for each
Cr            qp sought.  Typically a file would consists of lines like
Cr               q1x   q1y   q1z
Cr               q2x   q2y   q2z
Cr               ...
Cr     mode=3 generates qp for a uniform mesh in a plane (contour plot)
Cr            The file supplies information describing a rectangle in
Cr            the Brillouin zone.  It consists of a single line,
Cr            which contains the following:
Cr             v1    range  n     v2    range  n   height  list-of-bands
Cr
Cr            v1 and v2 are two vectors specifying the plane of the
Cr            contour.  range and n (one each for v1 and v2) are the
Cr            starting and final amplitudes of those vectors, and the
Cr            the number of points within the vector.  list-of-bands
Cr            is a list of integers which specify which bands are to
Cr            be written to the output file.  'height' is the 'z' axis.
Cr            For example,
Cr             v1    range  n     v2    range  n   height  list-of-bands
Cr            1 0 0  -1 1   51   0 1 0  -1 1   51   0.00    4,5
Cr            creates a file of 51x51 points, with the four corners
Cr            (-1,-1,0),  (1,-1,0),  (-1,1,0),  (1,1,0)
Cr            going through the gamma-point. Two bands (4,5) are stored.
Cu Updates
Cu   08 Jul 08 Extend to case where number of bands can be q dependent
Cu             modes 1,2: suqlsw writes out number of bands with qp
Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1 
c                                                          -->forbidden 2014takao
Cu   02 Jul 06 Color mode extended to two colors
Cu   02 Jul 06 New color mode (one color weight only)
Cu   14 Feb 05 contour mode saves both spins in spin-polarized case
Cu   20 Oct 03 suqlst works properly in contour mode
Cu   28 Aug 01 added ef switch
Cu   23 Jan 01 first written
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) sopts
      integer nq,nband,nsp,nfbn(2),ifblst(nband,*),onespdummy,iop!onesp
      double precision efermi,qp(3),evl(nband)
C ... Local variables
      character strn*120, strn2*120, dc*1, fn*120, prfmt*40
      logical rdstrn,fmlong
      integer i,iprint,j,j1,j2,k,nglob,rdm,
     .parg,a2vec,iv(7),stdo,jsp,op1
      double precision xx,xx1,xx2,sclp,plndst,xxv(7)
      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
      integer oqp,nqall,imode
      double precision x1,x2,y1,y2,evmxn,evmnn,
     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)
      real(8),allocatable :: evsav(:,:,:,:)
C ... MPI
      integer procid,master,mpipid

      logical:: l_dummy_isanrg,isanrg

ctakao
      real(8)::qlx(3)
C ... Heap
c      integer w(1)
c      common /w/ w
      common /suqlsd/
     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy
      save evsav
      real(8),allocatable:: w_oqp(:,:)

      procid = mpipid(1)
      master = 0
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      op1 = mod(iop,10)
      nqall = 0

C --- First call ... setup and parse options ---
      if (nq .eq. 0) then
C   ... Defaults
        fmlong = .false.
        mode = 1
        fn = 'syml'
        nblst = 0
        nevn = 0
        evmxn = -99d9
        evmnn = 99d9
        nfbn(1) = 0
        nfbn(2) = 0

        dc = sopts(1:1)
        if (dc .ne. ' ') then
C   ... Return here to resume parsing for arguments
          j2 = 0
   10     continue
          j2 = j2+1
          if (sopts(j2:j2) .eq. dc) goto 10
          j1 = min(len(sopts),j2)
          call nwordg(sopts,0,dc//' ',1,j1,j2)
          if (j2 .ge. j1) then
            if (.false.) then
            elseif (sopts(j1:j1+2) .eq. 'fn=')  then
              if (j1+3 .gt. j2) call rx('suqlst: bad file name')
              fn = sopts(j1+3:j2)
            elseif (sopts(j1:j2) .eq. 'qp')  then
              mode = 2
            elseif (sopts(j1:j2) .eq. 'spin1')  then
               call rx('suqlst:spin removed! takao')
c              onesp = 1
c              nsp = 1
            elseif (sopts(j1:j2) .eq. 'con')  then
              mode = 3
            elseif (sopts(j1:j2) .eq. 'long')  then
              fmlong = .true.
            elseif (sopts(j1:j1+3) .eq. 'col=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list, col=..')
              call mkils0(sopts(j1+4:j2),nfbn,ifblst)
              call mkilst(sopts(j1+4:j2),nfbn,ifblst)
              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
            elseif (sopts(j1:j1+4) .eq. 'col2=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list, col2=..')
              call mkils0(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
              call mkilst(sopts(j1+5:j2),nfbn(2),ifblst(1,2))
              if (nfbn(1) .lt. 0) call rx('suqlst: bad list, col=..')
            elseif (sopts(j1:j1+3) .eq. 'lst=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list, lst=..')
              call mkils0(sopts(j1+4:j2),nblst,iblst)
              if (nblst .gt. 100) call rx('increase size of iblst')
              call mkilst(sopts(j1+4:j2),nblst,iblst)
            elseif (sopts(j1:j1+2) .eq. 'ef=')  then
              j = 0
              i = parg('ef=',4,sopts(j1:),j,len(sopts(j1:)),
     .        dc//' ',1,1,i,efermi)
              if (i .le. 0) call rxs(
     .        'suqlst: failed to parse string for fermi level:  ',
     .        sopts(j1:))
            elseif (sopts(j1:j1+3) .eq. 'evn=')  then
              if (j1+4 .gt. j2) call rx('suqlst: bad list')
              i = j1+3
              xxv(1) = 0
              j = a2vec(sopts,j2,i,2,dc//' ',2,3,1,iv,nevn)
              if (j .ne. 1 .or. nevn .gt. nband)
     .        call rx('suqlst: bad value for evn')
            else
              call rxs('suqlst: failed to parse argument, ',sopts(j1:j2))
            endif
            goto 10
          endif
        endif

        if (nfbn(1) .eq. 0) then
          call info2(20,0,0,' suqlst:  generate bands, mode %i',mode,0)
        elseif (nfbn(2) .eq. 0) then
          call ilst2a(ifblst,nfbn,strn)
          call info2(20,0,0,
     .    ' suqlst:  generate bands with color weights, mode %i.'//
     .    ' %N%10f%i components: '//strn//'%a',mode,nfbn)
        else
          call ilst2a(ifblst,nfbn,strn)
          call info2(20,0,0,' suqlst:  generate bands '//
     .    'with two color weights, mode %i.',mode,0)
          call ilst2a(ifblst,nfbn,strn)
          call info2(20,0,0,
     .    '%10f%i components, first color:  '//strn//'%a',nfbn,0)
          call ilst2a(ifblst(1,2),nfbn(2),strn)
          call info2(20,0,0,
     .    '%10f%i components, second color: '//strn//'%a',nfbn(2),0)
        endif
        if (nfbn(1) .eq. 0 .and. nfbn(2) .ne. 0) call info0(20,0,0,
     .  ' suqlst:  (warning) col2 specified, but col not')

C   ... Open qp file
        if (procid .eq. master) then
          ifiq = fopno(fn)
          rewind ifiq

C   ... open bands file
          ifib = fopnn('BNDS')
          rewind ifib

C   ... Write header
          if (mode .eq. 1) then
            i = nblst
            if (nblst .eq. 0) i = nband
            if (nfbn(1) .eq. 0) then
C           Use separate format statment to circumvent gfortran bug
              write(ifib,335) i,efermi,0
  335         format(i5,f10.5,i6)
            elseif (nfbn(2) .eq. 0) then
              call ilst2a(ifblst,nfbn,strn)
              call strip(strn,j1,j2)
              write(ifib,336) i,efermi,1,strn(j1:j2)
  336         format(i5,f10.5,i6:'  col= ',a:'  col2= ',a)
            else
              call ilst2a(ifblst,nfbn,strn)
              call strip(strn,j1,j2)
              call ilst2a(ifblst(1,2),nfbn(2),strn2)
              call strip(strn2,j,k)
              write(ifib,336) i,efermi,2,strn(j1:j2),strn2(j:k)
            endif
          endif
        endif

C   ... Other initializations
        iq = 0
        nq = -1
      endif

C --- Setup for a new block of k-points, depending on mode ---
C     This branch occurs on completion of the last qp of the current block
C     which is marked by iq>nq
C     At the completion of this block: nq must be computed and:
C     (mode=1) q1,q2 set up.
C              Note: if 1s digit of iop is set in this mode,
C              this branch returns sum of all qp in all lines.
C              No setup for q1,q2; no
C     (mode=2) oqp allocated and loaded
C     (mode=3) nblst = number of bands to save
      if (iq .gt. nq) then
        iq = 1
C       Bands along specified symmetry lines
        if (mode .eq. 1) then
  725     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 999
          if (strn(1:1) .eq. '#') goto 725
C         Exit if first entry in line is zero
          i = 0
          xxv(1) = 0
          i = a2vec(strn,len(strn),i,4,', ',2,3,1,iv,xxv)
          if (i .eq. 1 .and. xxv(1) .eq. 0) goto 999
          i = 0
          i = a2vec(strn,len(strn),i,4,', ',2,3,7,iv,xxv)
          if (i .ne. 7 .and. iprint().ge.10) then
            write(stdo,
     .      '(/'' suqlst (warning) skipping line:''/''  '',a)') strn
            goto 725
          endif
          nq = xxv(1)
C         1 qp is nonsensical for a line
          if (nq .eq. 1) nq = 0
C         No qp: exit
          if (nq .le. 0) goto 999
C         setup q1,q2
          call dcopy(3,xxv(2),1,q1,1)
          call dcopy(3,xxv(5),1,q2,1)
          if (iprint().ge.10) write(stdo,785) nq,q1,q2
  785     format(/' suqlst:  nq=',i3,'   q1=',3f7.4,'   q2=',3f7.4)
C         Write line header information to disk
          if (op1 .eq. 0) then
            write(ifib,337) nq*nsp
  337       format(2i5)
C         Accumulate all qp and cycle until all lines are read
          elseif (op1 .eq. 2) then
            nqall = nqall + nq
            goto 725
          endif
C       Bands for a list of specified qp
        elseif (mode .eq. 2) then
C         Only one block for this mode.
C         Flag that prior block already completed: nq>0
          if (nq .gt. 0) goto 999
          nq = 0
          j = 3
          if (rdm(ifiq,10000,0,' ',xx,j,nq) .ne. 1) call
     .    rx('suqlst:  bad qp file')

c          call defrr(oqp,3*nq)
          allocate(w_oqp(3,nq))

          rewind ifiq
          j = rdm(ifiq,10000,3*nq,' ',w_oqp,j,nq)
          call awrit1('%N suqlst: read %i qp from file '//fn//
     .    '%a',' ',80,stdo,nq)
          if (j .ne. 1) call rx('suqlst: failed to read qp')
          i = nblst
          if (nblst .eq. 0) i = nband
          if (nfbn(1) .ne. 0) i = i*2
          call awrit2('%% rows %i cols %i',' ',80,ifib,nq,i+4)
C       Bands on a uniform mesh in a specified plane (contour plot)
        elseif (mode .eq. 3 .and. nq .eq. -1) then
  825     if (.not. rdstrn(ifiq,strn,len(strn),.false.)) goto 998
          if (strn(1:1) .eq. '#') goto 825
          call words(strn,i)
          if (i .ne. 14) call rxi('suqlst con mode: expected 14 '//
     .    'arguments from input file but read',i)
          iblst(1) = -1
          nblst = 1
          backspace ifiq
          read(ifiq,*) vectx,x1,x2,nqx,vecty,y1,y2,nqy,plndst
          call word(strn,14,j1,j2)
          call mkilss(11,strn(j1:j2),nblst,iblst)
          if (nblst .le. 0) call rx('suqlst: no bands in list')
          call dscal(3,1/dsqrt(sclp(vectx,vectx)),vectx,1)
          call dscal(3,1/dsqrt(sclp(vecty,vecty)),vecty,1)
C         Subtract from vecty projection onto vectx
C         call daxpy(3,-sclp(vectx,vecty),vectx,1,vecty,1)
          call cross(vectx,vecty,vect3)
          call dscal(3,plndst/dsqrt(sclp(vect3,vect3)),vect3,1)
          nq = nqx*nqy
          if (iprint() .ge. 10) then
            write(stdo,717) vectx,x1,x2,nqx,plndst,vecty,y1,y2,nqy,nq
  717       format(' vx=',3f9.6,'  x1,x2=',2f9.6,'  nx=',i3,'   h=',f9.6/
     .      ' vy=',3f9.6,'  y1,y2=',2f9.6,'  ny=',i3,'  np=',i5)
            write(strn,'('' save %i bands: %'',i2,'':1i'')') nblst
            if (strn(17:18) .eq. '% ') strn(17:18) = ' %'
            call awrit2(strn,strn,80,stdo,nblst,iblst)
            if (dabs(sclp(vectx,vecty)) .gt. 1d-6 .and. iprint().gt.20)
     .      print *,'suqlst (warning): input plane vectors not orthogonal'
          endif

          if (iop .eq. 0) then
            allocate(evsav(nqx,nqy,nblst,nsp))
          endif

C       Contour plot, cleanup.
C       Note: cleanup handled by suqlsw when last qp is called.
        elseif (mode .eq. 3) then
C          if (op1 .ne. 0) then
C            allocate(evsav(nqx,nqy,nblst,nsp))
C            return
C          endif
          call rx('suqlst: caller should never reach this branch')
        endif

C --- Generate qp for this iq, depending on mode ---
      else
        if (mode .eq. 1) then
          xx = dble(iq-1)/dble(nq-1)
          qp(1) = xx*q2(1) + (1-xx)*q1(1)
          qp(2) = xx*q2(2) + (1-xx)*q1(2)
          qp(3) = xx*q2(3) + (1-xx)*q1(3)
        elseif (mode .eq. 2) then
          call dpscop(w_oqp,qp,3,iq*3-2,1,1d0)
        elseif (mode .eq. 3) then
C         Inner Loop:  excursions in y; outer loop: excursions in x
          j = mod(iq-1,nqy)
          i = (iq-1-j)/nqy
          if (nqx .le. 1) then
            xx1 = x1
          else
            xx1 =i*(x2-x1)/(nqx-1)+x1
          endif
          if (nqy .le. 1) then
            xx2 = y1
          else
            xx2 =j*(y2-y1)/(nqy-1)+y1
          endif
          do  93  k = 1, 3
            qp(k) = xx1*vectx(k) + xx2*vecty(k) + vect3(k)
   93     continue
          if (j .eq. 0 .and. iprint().ge.20) write(stdo,718) i+1,nqx,qp
  718     format(' line',i3,' of',i3,'   q(1)=',3f10.6)
        else
          call rx('suqlst: bad mode')
        endif
        iq = iq+1
C       Hold onto local copy of qp
        call dcopy(3,qp,1,ql,1)
      endif
      if(allocated(w_oqp)) deallocate(w_oqp)
      return

C --- No more qp blocks: cleanup ---
  999 continue
      nq = 0
      if (nevn .ne. 0) call awrit3(' eval no. %i:  minimum eval'//
     .' = %;8F  maximum eval = %;8F',' ',80,stdo,nevn,evmnn,evmxn)

      if (mode .eq. 1) then
        if (op1 .eq. 0) then
          write(ifib,337) 0
        elseif (op1 .eq. 2) then
          nq = nqall
        endif
      elseif (mode .eq. 3) then
        call rx('not ready')
C       call xxxbnd(w(oev),nblst,nqx,nqy,ifib)
C       return
      endif
      if(allocated(w_oqp)) deallocate(w_oqp)
      return

C --- Error exit ---
  998 call rxs('suqlst: failed to read file contents, file ',fn)

      entry suqlsm(imode)
C- Return qlist mode
      imode = mode
      return

      entry suqlsw(nband,jsp,nsp,evl)
C- Write or store the energy bands to file for this qp

      if (mode .eq. 1) then
        if (nblst .eq. 0) then
          prfmt = '(3f10.5,i6/(10f8.4))'
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.5/(10f8.4))'
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        else
          call rx('bug in suqlsw')
        endif

      elseif (mode .eq. 2) then
        if (fmlong .and. nblst .eq. 0) then
          prfmt = '(3f15.10,i6/(5f15.10))'
        elseif (nblst .eq. 0) then
          prfmt = '(3f10.6,i6/(8f10.6))'
        elseif (fmlong .and. nblst .gt. 0) then
          prfmt = '(3f15.10/(5f15.10))'
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.6/(8f10.6))'
C         if (nblst .le. 5) prfmt = '(8f10.6)'
        else
          call rx('bug in suqlsw')
        endif
        if (nblst .eq. 0) then
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        else
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 3) then
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
Changenglob        stdo = nglob('stdo')
        stdo = globalvariables%stdo
        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)

        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
            do  k  = 1, nblst
              if (nsp .eq. 1)
     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
              if (nsp .eq. 2)
     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
     .        ifib,nqx,nqy,j1)

              do  i = 1, nqx
                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
              enddo
            enddo
          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif

C ... Keep running tab on smallest,largest eval
      if (nevn .ne. 0) then
        evmxn = max(evmxn,evl(nevn))
        evmnn = min(evmnn,evl(nevn))
      endif


      entry suqlsw2(nband,jsp,nsp,evl,qlx) !takao added ql
C- Write or store the energy bands to file for this qp
      ql=qlx
      if (mode .eq. 1) then
        if (nblst .eq. 0) then
          prfmt = '(3f10.5,i6/(10f8.4))'
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.5/(10f8.4))'
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        else
          call rx('bug in suqlsw')
        endif
      elseif (mode .eq. 2) then
        if (fmlong .and. nblst .eq. 0) then
          prfmt = '(3f15.10,i6/(5f15.10))'
        elseif (nblst .eq. 0) then
          prfmt = '(3f10.6,i6/(8f10.6))'
        elseif (fmlong .and. nblst .gt. 0) then
          prfmt = '(3f15.10/(5f15.10))'
        elseif (nblst .gt. 0) then
          prfmt = '(3f10.6/(8f10.6))'
C         if (nblst .le. 5) prfmt = '(8f10.6)'
        else
          call rx('bug in suqlsw')
        endif
        if (nblst .eq. 0) then
          write(ifib,prfmt) ql, nband, (evl(i),i=1,nband)
        else
          write(ifib,prfmt) ql, (evl(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 3) then
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
Changenglob        stdo = nglob('stdo')
        stdo = globalvariables%stdo
        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
c  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)

        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = evl(iblst(k))
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
            do  k  = 1, nblst
              if (nsp .eq. 1)
     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
              if (nsp .eq. 2)
     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
     .        ifib,nqx,nqy,j1)

              do  i = 1, nqx
                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
              enddo
            enddo
          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif

C ... Keep running tab on smallest,largest eval
      if (nevn .ne. 0) then
        evmxn = max(evmxn,evl(nevn))
        evmnn = min(evmnn,evl(nevn))
      endif
      end

      subroutine suqlse(nband,jsp,nsp,ndimhx,ifbn,nfbn,ifblst,ndlst,
     .evec,wk)
      use m_globalvariables
C- Write to file the projection of eigenvector subblock for this qp
C ----------------------------------------------------------------------
Ci Inputs
Ci   nband :number of energy bands to write
Ci   jsp   :current spin index (not used now)
Ci   nsp   :number of spins (not used now)
Ci   ndimhx:dimensions evec
Ci   ifbn  :index to color list (nfbn,ifblst)
Ci   nfbn  :number of elements for color weights projection
Ci   ifblst:list of elements for  color weights projection
Ci   ldlst :leading dimension of ifblst
Ci   evec  :eigenvectors
Ci   wk    :work array of same dimension as evl
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr   Inefficient, but it works
Cu Updates
Cu   08 Jul 08 New argument ndlst so ifblst can be dimensioned
Cu             independently from nband
Cu   05 Jun 06 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer jsp,nsp,nband,ifbn,nfbn(2),ndlst,ifblst(ndlst,*)
      integer ndimhx
      double precision wk(1)
      double complex evec(ndimhx,ndimhx)
C ... Local parameters
      double complex zdotc
      character prfmt*40
      logical fmlong
      integer i,iprint,j,j1,k,nglob,stdo
      integer iblst(100),ifib,ifiq,iq,mode,nblst,nevn,nqx,nqy
c      integer oqp
      double precision x1,x2,y1,y2,evmxn,evmnn,
     .q1(3),q2(3),ql(3),vectx(3),vecty(3),vect3(3)

      common /suqlsd/
     .q1,q2,ql,evmxn,evmnn,ifiq,ifib,mode,iq,oqp,nevn,nblst,
     .iblst,fmlong,nqx,nqy,x1,x2,y1,y2,vectx,vecty,vect3
c      common /suqlsd2/ql,ifib,mode,iq,nblst,iblst,fmlong,nqx,nqy

      integer ipiv(ndimhx)
      real(8),allocatable :: evsav(:,:,:,:)
      complex(8),allocatable:: evecc(:,:),work(:,:)

      allocate(evecc(ndimhx,ndimhx),work(ndimhx,ndimhx))

#if MULLIKEN
      call zcopy(ndimhx**2,evec,1,evecc,1)
      call zgetrf(ndimhx,ndimhx,evecc,ndimhx,ipiv,j)
      if (j .ne. 0) call rx('mullmf: failed to generate overlap')
      call zgetri(ndimhx,evecc,ndimhx,ipiv,work,ndimhx**2,j)
#endif
      do  i = 1, ndimhx
        wk(i) = 0
        do  j = 1, nfbn(ifbn)
          k = ifblst(j,ifbn)
          if (k .le. 0 .or. k .gt. ndimhx) then
            call fexit2(-1,111,' Exit -1 : suqlst: component %i'//
     .      ' outside range (1:%i)',k,ndimhx)
          endif
#if MULLIKEN
          wk(i) = wk(i) + evecc(i,k)*evec(k,i)
#else
          wk(i) = wk(i) + dconjg(evec(k,i))*evec(k,i)
#endif
        enddo

#ifndef MULLIKEN
        wk(i) = wk(i) / zdotc(ndimhx,evec(1,i),1,evec(1,i),1)
        if (wk(i) .lt. 0 .or. wk(i) .gt. 1+1d-14) then
          call rx('bug in suqlse')
        endif
#endif

      enddo

      if (mode .eq. 1) then
        prfmt = '(3f10.5/(10f8.4))'
        if (nblst .eq. 0) then
          write(ifib,prfmt) ql, (wk(i),i=1,nband)
        else
          write(ifib,prfmt) ql, (wk(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 2) then
        if (fmlong) then
          prfmt = '(5f15.10)'
        else
          prfmt = '(8f10.6)'
C         if (nblst .gt. 0 .and. nblst .le. 5) prfmt = '(8f10.6)'
        endif
        if (nblst .eq. 0) then
          write(ifib,prfmt) (wk(i),i=1,nband)
        else
          write(ifib,prfmt) (wk(iblst(i)),i=1,nblst)
        endif

      elseif (mode .eq. 3) then
        call rx('need copy weights into different place, mode=3')
        j = mod(iq-2,nqy)
        i = (iq-2-j)/nqy
Changenglob        stdo = nglob('stdo')
        stdo = globalvariables%stdo
        if (iprint() .ge. 60) write(stdo,345) i+1,j+1,ql
  345   format(' saving point, iq=',i5,' jq=',i5,'  qp=',3f12.5)

        do  k = 1, nblst
          evsav(i+1,j+1,k,jsp) = wk(iblst(k))
        enddo

C   ... If last qp generated, dump to file and exit
        if (i+1 .eq. nqx .and. j+1 .eq. nqy .and. jsp .eq. nsp) then
          if (fmlong) then
            prfmt = '(5f15.10/(5f15.10))'
          else
            prfmt = '(8f10.6/(8f10.6))'
          endif
          rewind ifib
          do  j1 = 1, nsp
            do  k  = 1, nblst
              if (nsp .eq. 1)
     .        call awrit2('%% rows %i cols %i',' ',80,ifib,nqx,nqy)
              if (nsp .eq. 2)
     .        call awrit3('%% rows %i cols %i spin %i',' ',80,
     .        ifib,nqx,nqy,j1)

              do  i = 1, nqx
                write(ifib,prfmt) (evsav(i,j,k,j1), j=1,nqy)
              enddo

            enddo
          enddo
          call rx0('finished generating bands on q-mesh')
        endif
      endif

      deallocate(evecc,work)

      end

      subroutine suqlsr(mode,ifi,nsp,nbf,lde,nb1,nb2,ls1,ls2,nq,qp,eb)
C- Count qp, optionally read them and energy bands from file
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :Any combination of the following is allowed:
Ci         :1 count number of qp in file, return as nq
Ci         :  if 1s bit mode is 0, nq is an input,
Ci         :  the file nq cannot exceed input nq
Ci         :2 return qp read in file into array qp
Ci         :4 return bands read in file into array eb
Ci   ifi   :read from logical unit ifi
Ci   nsp   :number of spins in band file
Ci   nbf   :number of bands in band file
Ci   lde   :leading dimension of eb
Ci   nb1,nb2: read bands nb1..nb2 into eb
Ci   ls1   :(nsp=2): read first spin only
Ci   ls2   :(nsp=2): read second spin only
Cio Inputs/Outputs
Cio  nq    :numbefr of k-points (returned if 1s bit of mode is set)
Co Outputs
Co   qp    :k-points, returned if 4s bit of mode is set
Co   eb    :energy bands, returned if 4s bit of mode is set
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   06 Jun 07  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ifi,nsp,nq,nbf,lde,nb1,nb2
      logical ls1,ls2
      double precision qp(3,nq),eb(lde,nsp,nq)
C ... Local parameters
      integer nq1,nql,mode1,mode2,mode4,iq1,i,nline,isp
      double precision ql(3),evl(3000)
      logical :: l_dummy_isanrg,isanrg

      if (mode .eq. 0) return
      mode1 = mod(mode,2)
      mode2 = mod(mode/2,2)
      mode4 = mod(mode/4,2)
      if (mode4 .ne. 0) then
Ckino isanrg is logical function,         call isanrg(nb2,1,min(lde,nbf),'suqlsr:','top band index',
        l_dummy_isanrg=isanrg(nb2,1,min(lde,nbf),'suqlsr:','top band index',
     .  .true.)
Ckino isanrg is logical function,         call isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
        l_dummy_isanrg=isanrg(nb1,1,nb2,'suqlsr:','bottom band index',.true.)
      endif

C --- For each panel, do ---
      nql = 0
      nline = 0
   91 continue
      read(ifi,*) nq1
      if (nq1 .le. 0) goto 90
      isp = 0
      do  iq1 = 1, nq1
        isp = mod(isp,2)+1
        nql = nql+1
        if (mode1 .eq. 0 .and. nql .gt. nq) call rxi(
     .  'suqlsr: file more q-points than allocated: nqmx=',nq)
        read(ifi,*,END=999,ERR=999) ql(1),ql(2),ql(3)
        if (mode2 .ne. 0) call dcopy(3,ql,1,qp(1,nql),1)
        read(ifi,*,END=999,ERR=999) (evl(i),i=1,nbf)
        if (mode4 .ne. 0) then
C         Copy only if appropriate spin
          if (isp .eq. 1 .and. ls2) then
          elseif (isp .eq. 2 .and. ls1) then
          else
            call dcopy(nb2-nb1+1,evl(nb1),1,eb(1,1,nql),1)
          endif
        endif
      enddo
      nline = nline+1
      goto 91
C     End of loop over lines
   90 continue


      if (mode1 .eq. 0) then
        call info2(30,1,0,' suqlsr: found %i qp in %i lines from file',
     .  nql,nline)
      else
        nq = nql/nsp
        call info2(30,1,0,' suqlsr: read %i qp in %i lines from file',
     .  nql,nline)
      endif

      if (mod(nql,2) .ne. 0 .and. nsp .eq. 2) then
        call info0(20,0,0, ' suqlsr (warning): '//
     .  'odd number of bands encountered but nsp=2')
      endif

      return
  999 continue
      call rxi('suqlsr: failed to read bands file, nq=',nql)

      end

      subroutine xxxbnd(evl,nbnd,nqx,nqy,ifib)
C     implicit none
      integer nbnd,nqx,nqy,ifib,ix,iy,ib,ifi,fopna
      double precision evl(nbnd,nqy,nqx)
      character*10 strn

      do  10  ib = 1, nbnd
        ifi = ifib
        if (nbnd .ne. 1) then
          strn = ' '
          call awrit1('bnd%i',strn,10,0,ib)
          ifi = fopna(strn,30,0)
        endif
        write(ifi,337) nqx,nqy
  337   format('% rows ',i5,' cols ',i5)
        do  20  ix = 1, nqx
          write(ifi,333) (evl(ib,iy,ix), iy=1,nqy)
  333     format(6f12.6)
   20   continue
        if (nbnd .ne. 1) call fclose(ifi)
   10 continue
      end
      double precision function sclp(v1,v2)
      double precision v1(3),v2(3)
      sclp  = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
      end

      
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs
  integer function rdtok(token,instr,descr,sep0,sep,cast,swt,
     .n1,count,resn,resc)
C- Seek a token within a string and read associated data into array
C ----------------------------------------------------------------
Ci Inputs
Ci   token: string marking data to follow token
Ci   instr: string in which to seek token
Ci   descr: descriptor string composed of one or more parts.
Ci        : It is used to format printouts generated by rdtok,
Ci        : in 'print' modes
Ci        :*The first part, if present, is a formatting string
Ci        : for writing the contents of a token, awrite conventions
Ci        : It is delimited by a double-colon, e.g.
Ci        :    descr = '%;6d,,string-for-remaining-parts'
Ci        : The format statement is used when the token's contents
Ci        : are written out (swt(5)>2, below).  A default format
Ci        : is used if one is not supplied
Ci        :*The second part, if present, is a character consisting
Ci        : of a single digit '0'..'9', and is the number of spaces
Ci        : to indent before writing the token, in all modes swt(5).
Ci        : It is delimited by a double-semicolon, e.g.
Ci        :    descr = '%;6d,,4;;string-for-remaining-part'
Ci        : If not present, a default indentation of 2 is used.
Ci        :*The third part, if present, describe the token's purpose.
Ci        : It is used only for input mode 0 (swt(5)=0).
Ci        : It consists of any part of the string after the first or
Ci        : second parts.
Ci    sep0: the character immediately to the left of the token
Ci          must be one of the characters in sep0
Ci    sep:  The token's arguments are delineated by characters in
Ci       :  this list.
Ci    cast: data following token is one of the following.
Ci          cast    class of token
Ci           0      logical
Ci           1      char
Ci           2      int
Ci           3      real
Ci           4      double
Ci        : result is stored into resn for cast=0,2,3,4
Ci        :           stored into resc for cast=1
Ci     n1 : Data associated with token are put into:
Ci          resn(n1),resn(n1+1),... for cast=0,2,3,4    and
Ci          resc(n1),resc(n1+1),... for cast=1
Ci   count: The number of elements (numbers or strings) to try and parse
Ci        : It may be permissible to read fewer than `count' elements;
Ci        : see swt(3)
Ci Inputs/Outputs
Cio  swt  : a vector of switches governing how gtcat and rdtok behave.
Cio       : Its internal structure is described in stswt, below.
Cio       : User is advised to set these parameters via routine stswt.
Cio       : and query the state of using quswt; then there is no need
Cio       : know the internal structure of swt.
Cio       : As described in detail in stst, quantities that can be set
Cio       : and probed are the following:
Cio       : structure of swt.
Cio       : swt()  function (see stswt)
Ci        :  1     offset to instr where parsing begins
Cio       :  2     parse up to this character
Ci        :  3     governs behavior when token or some data is missing,
Ci        :        whether tokens are optional or required, if
Ci        :        there may be multiple occurences of this token, and
Ci        :        whether number of elements may be less than count.
Ci        :  5     read and printout mode
Ci        :        0=> do not attempt to read token but show
Ci        :            what token would have been sought
Ci        :        1=> attempt to find token and read contents
Ci        :        2=> print token name
Ci        :        3=> combination of 1 and 2
Ci        :  6     file logical unit for printout
Co        :  7     index to last character parsed
Co Outputs
Co   rtdok: -999   No token was sought
Co        :   0    No match to token is found
Co        :   n    token matched and converted sans error n elements
Co        :        (for count=0, returns 1 if token matched)
Co        :  -n    if error on conversion of argument n
Cl Local variables
Cl   lopt :  0 if token is optional
Cl        :  1 if token is required
Cl        :  2 alternate token may be subst.
Cl        :  3 if token is ignored
Cb Bugs
Cb   when n1>0, data not placed in proper location for integer, logical
Cb   arrays
Cu Updates
Cu   06 Nov 01 Initially created
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer swt(7),cast,n1,count
      character*(1) instr(*)
      character*(*) token,descr,resc(1),sep0,sep
      double precision resn(1)
C ... Local parameters
      logical logval,parstr,l1,l2,lhelp,lshow,lread
      integer a2vec,awrite,i,ifipr,ifmt,io,ip,is,isave,it(99),ival,k,ks,
     .linstr,rdtk2,recl,retval,ltok,quswt,lfixl,lopt,lmult,lio
      parameter (recl=500)
      character   aa*(recl),prs*(recl)
      character*8 nmcast(0:4),fmt*20
      save ip,prs
      data nmcast /'logical','char','integer','real','double'/
      data ip /0/

C ... Setup
      lopt  = quswt('token,opt',swt)
      if (lopt .eq. 3) then
        rdtok = -999
        return
      endif
      ifipr = quswt('lunit',swt)
      ltok = len(token)
      lhelp = quswt('io,help',swt) .eq. 1
      lshow = quswt('io,w',swt) .eq. 1
      lread = quswt('io,r',swt) .eq. 1
      lio   = quswt('io',swt)
      lfixl = quswt('token,fixlen',swt)
      lmult = quswt('token,mult',swt)

C --- Show what token would have been sought ---
      if (lhelp) then

C       Find indentation, offset to descriptor string
C       Default indentation
        ks = 2
C       Default offset to start of descriptor
        i = 0
        if (descr .ne. ' ') then
C         Strip off first parts, if either exists
          l2 = parstr(descr,';;',len(descr)-1,2,'a',i,k)
          if (l2) then
            ks = ichar(descr(i:i)) - ichar('0')
          else
            i = 0
            l1 = parstr(descr,',,',len(descr)-1,2,'a',i,k)
            if (.not. l1) i = -2
          endif
        endif
        fmt = ' '
        is= ltok
        if (count .eq. 0) then
          aa = fmt(1:ks) // ' *token  ' // token(1:is) //
     .    '%a (no arguments; %?#(n==0)#optional#required#)'
          call awrit1(aa,' ',72,ifipr,lopt)
        elseif (count .eq. 1) then
          aa = fmt(1:ks) // ' *token  ' // token(1:is) //
     .    '%a takes 1 argument of cast '//nmcast(cast)//
     .    '%a %?#(n==0)#(optional)##%-1j'//
     .    '%a %?#(n==1)#(required)##%-1j'//
     .    '%a %?#(n==2)#%N'//fmt(1:ks)//
     .    '   (choose either this token or the next one)##'
          call awrit1(aa,' ',128,ifipr,lopt)
        else
          aa = fmt(1:ks) // ' *token  ' // token(1:is) //
     .    '%a takes %?#(n==0)#up to ##'//
     .    '%i arguments of cast '//nmcast(cast)//
     .    '%a %?#(n==0)#(optional)##%-1j'//
     .    '%a %?#(n==1)#(required)##%-1j'//
     .    '%a %?#(n==2)#%N'//fmt(1:ks)//
     .    '   (choose either this token or the next one)##'
          call awrit3(aa,' ',128,ifipr,lfixl,count,lopt)
        endif
        if (descr .ne. ' ') then
          fmt = ' '
          aa = '   ' // fmt(1:ks) // descr(i+3:)
          call awrit0(aa,' ',-len(aa),ifipr)
        endif
        rdtok = -999
        ip = 0
        return

C --- No token to read ---
      elseif (.not. lread) then
        retval = -999
        goto 999
      endif

C --- Try and find token ---
      linstr = swt(2)
      if (token .eq. ' ') then
        is = swt(1)
      else
        is = rdtk2(instr,linstr,token,sep0,lopt.eq.1,swt(1))
      endif
C ... No token found ... exit
      if (is .le. 0) then
        swt(7) = linstr
        rdtok = 0
        return
      endif

C --- Parse for token arguments ---
      swt(7) = is

C ... Case count=0
      if (count .eq. 0) then
        retval = 1
        goto 999
      endif

C ... Case further delimit range of search
      if (lmult .ne. 0) then
        ks = rdtk2(instr,linstr,token,sep0,.false.,is)
        if (ks .gt. 0) linstr = ks-ltok-1
        swt(2) = linstr
      endif

C ... Parse for vector of character arguments
      if (cast .eq. 1) then
        isave = is - ltok
        do  i = 1, count

C         Skip past leading blanks
          is = is-1
          call skipbl(instr,linstr,is)
          is = is+1

C         End-of-category encountered before loop over count exhausted
          if (is .ge. linstr) then
            if (lfixl .ne. 0) then
              if (count .eq. 1) then
                aa = 'Exit -1 rdtok: no string after token '//
     .          token(1:ltok)
                call fexit(-1,1,aa,0)
              else
                aa = 'Exit -1 rdtok: failed to parse %i strings, token'
                call strncp(aa,instr,53,isave,min(recl-53,linstr-isave)-1)
              endif
              call fexit(-1,1,aa,count)
            endif
            retval = -i
            swt(7) = linstr+1
            goto 999
          endif

          call cpstr(instr,linstr,1001,sep,is,io,resc(i+n1-1))
          is = is+1
        enddo
        swt(7) = is-1
        retval = count
        goto 999
      endif

C ... Parse for vector of numerical arguments
      isave = is - ltok
      swt(7) = is
      is = is-1
      k = a2vec(instr,linstr,is,cast,sep,len(sep),-99,-count,it,
     .resn(n1))
      retval = k

C ... a2vec encountered an illegitimate expression
      if (k .lt. 0) then
        if (lfixl .ne. 0) then
          aa = ' Exit -1 rdtok: error parsing expr no %i, token'
          call strncp(aa,instr,52,isave,min(recl-52,linstr-isave)-1)
          call fexit(-1,1,aa,-k)
        endif
C ... a2vec read fewer than count expressions
      elseif (k .lt. count) then
        if (lfixl .ne. 0) then
          aa = ' Exit -1 rdtok: read only %i of %i values, token'
          call strncp(aa,instr,52,isave,min(recl-52,linstr-isave)-1)
          call fexit2(-1,1,aa,k,count)
        endif
        retval = -k
      endif

C --- Cleanup: check for printout ---
  999 continue
      rdtok = retval
      if (.not. lshow) return

C     Find format, if supplied, and indentation
C     Default indentation
      ks = 2
      fmt = ' '
      if (descr .ne. ' ') then
        i = 0
        l2 = parstr(descr,';;',len(descr)-1,2,'a',i,k)
        if (l2) then
          ks = ichar(descr(i:i)) - ichar('0')
        endif
        ifmt = 0
        l1 = parstr(descr,',,',len(descr)-1,2,'a',ifmt,k)
        if (l1) fmt = descr(1:ifmt)
      endif

C ... Write token name to string
      if (ip .eq. 0) then
        prs = ' '
      else
        ks = min(ks,2)
      endif
      ip = ks + ip + awrite(token,prs(ks+ip+1:),128,0,0,0,0,0,0,0,0,0)
      if (count .eq. 0) then
C   ... format applied after token write
C        if (fmt .ne. ' ')
C     .  ip = awrite(fmt,prs,128,0,0,0,0,0,0,0,0,0)
        goto 998
      endif

C ... Write token contents to string, vec of length |retval|
C     ip = ip+1

C     Number of elements to write
      k = retval
      if (retval .lt. 0) k = -retval-1
      if (lio .eq. 2) k = count
C     Write the elements
      do  i = 1, k
        if (cast .eq. 1) then
          ip = ip + awrite(resc(n1+i-1),prs(ip+1:),128,
     .    0,0,0,0,0,0,0,0,0)
          call skpblb(prs,ip,ip)
          ip = ip+1
          if (fmt .ne. ' ')
     .    ip = awrite(fmt(1:ifmt),prs,128,0,0,0,0,0,0,0,0,0)
        elseif (cast .eq. 0) then
          if (fmt .eq. ' ') then
            fmt = '%l'
            ifmt = 2
          endif
          ip = ip + awrite(fmt(1:ifmt),prs(ip+1:),128,0,
     .    logval(resn,n1+i-1),0,0,0,0,0,0,0)
        elseif (cast .eq. 2) then
          if (fmt .eq. ' ') then
            fmt = '%i'
            ifmt = 2
          endif
          ip = ip + awrite(fmt(1:ifmt),prs(ip+1:),128,0,
     .    ival(resn,n1+i-1),0,0,0,0,0,0,0)
        elseif (cast .eq. 4) then
          if (fmt .eq. ' ') then
            fmt = '%;10F'
            ifmt = 5
          endif
          ip = ip + awrite(fmt(1:ifmt),prs(ip+1:),128,0,resn(n1+i-1),
     .    0,0,0,0,0,0,0)
        endif
C       ip = ip+2
        if (i .lt. k) then
          if (prs(ip:ip) .ne. sep(1:1)) then
            ip = ip+1
            write(prs(ip:ip),'(a)') sep(1:1)
          endif
        elseif (retval .lt. 0 .and. lio .ne. 2) then
          ip = ip+1
          write(prs(ip:),'(a,''***'')') sep(1:1)
          ip = ip+3
        endif
      enddo

  998 continue
c     call skpblb(prs,ip,ip)
C     ip = ip+1

C ... Printout, or just preserve in string prs if ifipr is 0
      if (ifipr .gt. 0) then
        write(ifipr,'(a)') prs(1:ip)
        ip = 0
      else
      endif

      end

      integer function rdtk2(instr,linstr,token,sep0,lreqd,i1)
C- Kernel called by rdtok to find a token
C ----------------------------------------------------------------------
Ci Inputs
Ci   instr :string in which to seek token
Ci   linstr:length of instr
Ci   token :token
Ci    sep0 :the character immediately to the left of the token
Ci         :must be one of the characters in sep0
Ci   lreqd :T, abort with error if token is not found
Ci    i1   :start search at i1
Co Outputs
Co   rdtk2 :if token found, rdtk2 points to first character past token
Co         :if no token found, rdtok = -linstr
Cu Updates
Cu   02 Nov 01 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lreqd
      integer linstr,i1
      character*(*) token,sep0
      character*(1) instr(linstr)
C ... Local parameters
      integer ltok,k,iterm,is,lgunit,stderr
      character ct*1

      ltok = len(token)
      ct = token
      stderr = lgunit(2)

C --- Loop over length of string ---
C     Match by:
C     1.  First character matches at instr(k)
C     2.  character to left matches on of sep0
C     3.  token matches
      do  k = i1, linstr

        if (instr(k) .eq. ct) then
          iterm = 1
          if (k .gt. 1) then
            is = 0
            call chrps2(instr(k-1),sep0,len(sep0),0,is,iterm)
          endif
          if (iterm .gt. 0) then
            do  is = 1, ltok
              if (instr(k+is-1) .ne. token(is:is)) goto 31
            enddo
C           A match was found
            rdtk2 = ltok+k
            return
          endif
        endif
   31   continue
      enddo
C ... We have fallen through the loop, so no token found
      rdtk2 = -linstr
      if (lreqd) then
        write(stderr,'(1x)')
        call rxs2('rdtok: failed to find token `',token,'''')
      endif
      end


cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
  subroutine a2rotm(strn,linv,ipr,rotm)
C- Generate Euler angles from a sequence of rotations, ascii input
Ci strn: Sequence of rotation matrices.  Syntax: rot1[,rot2...],
Ci   where each rotj looks like (x,y,z)angle.  x:, y: or z: may be
Ci   substituted for (x,y,z).  Example:  strn that makes rotation matrix
Ci   corresponding to the Euler angles alpha=pi/4 beta=pi/3 gamma=pi/2:
Ci   (0,0,1)pi/4,(0,1,0)pi/3,(0,0,1)pi/2   or just  z:pi/4,y:pi/3,z:pi/2
Ci linv: T, return transpose of rotm, corresponding to rotation of
Ci   positions rather than the coordinates (default).
Ci ipr: verbosity
Co rotm: rotation matrix.  Follow with a call to rm2eua to find Euler
Co   angles.
C     implicit none
      logical linv
      character*(*) strn, rchr*7, ss*80
      double precision rotm(3,3)
      double precision alpha,beta,gamma,rotj(3,3),r2(3,3),angle,
     .r(3),pi,ddot,fuzz,ct,theta,phi

      integer ls,itrm,ich,a2vec,ix(3),i,j,irot,ipr,i1mach
      parameter (fuzz = 1d-12)
      data rchr /'(XxYyZz'/


      pi = 4*datan(1d0)
      call dpzero(rotm,9)
      rotm(1,1) = 1
      rotm(2,2) = 1
      rotm(3,3) = 1
      ss = strn
      ls = len(ss)
      irot = 0
      ich = 0
C --- Entry point to accumulate new rotation ---
    5 continue
      irot = irot+1
      call chrps2(ss,rchr,len(rchr),ich,ich,itrm)
      if (itrm .eq. 0) goto 999
      if (itrm .gt. 1) then
        call dpzero(r,3)
        if (ss(ich+2:ich+2) .ne. ':') goto 999
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
