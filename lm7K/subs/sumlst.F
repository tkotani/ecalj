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
            call isanrg(mode,0,5,' sumlst:','mode',.true.)
            if (ipr .ge. 10) then
              if ((mode .eq. 2 .or. mode .eq. 5) .and. ng .gt. 1)
     .        call awrit0(' **WARNING** sumlst: '//
     .        'for lm-decomposition suppress symops',' ',128,stdo)
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
        call isanrg(ib,1,nbas,' sumlst:','site',.true.)

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
        call awrit4(' sumlst:  %?#n#Mulliken#Partial DOS# mode %i ('//
     .  modstr(mode)(i:j)//')  %i sites %i channels',' ',128,stdo,
     .  lopt0,mode,nsites,nchan)
      endif
      return

  999 continue
      call rxs('sumlst: failed to parse options in ',sopts(j1:))

      end subroutine sumlst


