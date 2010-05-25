      subroutine suemph(sspec,hamrul,ntab,iax,rtab,nbas,nttab,ips,
     .ntype,ltype,htype,esite,isite,epair,ipair)

      use m_struc_def  !Cgetarg

C- Assemble the coefficients to an empirical hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   sspec :struct containing species-specific information
Ci     Elts read: z
Ci     Stored:    name
Ci   hamrul: a string containing set of rules defining the
Ci         : the hamiltonian.  lhrule is the length of string hamrul.
Ci   nbas  : number of sites
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   ntype : ntype(1) number types of on-site interactions
Ci         : ntype(2) number types of pairwise interactions
Ci   ltype : a vector of integers specifying which interactions sought.
Ci         : The vector consists of ntype(1) integers describing the
Ci         : on-site interactions, followed by ntype(2) integers
Ci         : describing the pairwise interactions.
Ci         : The vector corresponds to the hamiltonian types
Ci         : ltype = 0  => interaction is not used
Ci         : ltype = 1  => interaction is not required
Ci         : ltype =-1  => interaction is required
Ci   htype : a vector of strings naming each of the ntype interactions.
Ci         : The vector consists of ntype(1) strings naming the
Ci         : on-site interactions, followed by ntype(2) strings naming
Ci         : the pairwise interactions.
Co Outputs
Co   esite : onsite interaction
Co   isite : contains which rule was used for onsite interaction
Co   epair : pairwise interaction
Co   ipair : contains which rule was used for pair interaction
Cr Remarks
Cr  *hamrul is a string containing a sequence of rules defining the
Cr   hamiltonian.  A particular rule has a syntax
Cr        htype value [expr1] [expr2] [...]
Cr   h-type is a string specifying what kind of interaction, eg on-site,
Cr   pairwise, biquadratic.  These strings are passed in array htype.
Cr  *For now, the interactions must be either of the one-center or
Cr   the two-center type, and only one kind of one-center and one
Cr   kind of two-center hamiltonian is allowed.
Cr  *The expressions expr1, etc. are integer expressions that are evaluated
Cr   for each site(pair), and if any of the expr evaluate to nonzero,
Cr   that particular rule is used for that site(pair).  expr1 etc are
Cr   evaluated in a site(pair)-specific way by setting variables in the
Cr   symbolic variables specific to each site(pair) (see below).
Cr   Thus, the first rule suemph encounters in which one of expr1, expr2,
Cr   etc. evaluates to nonzero, is the rule used to define the
Cr   interaction for that pair.  The interaction is taken to be zero for
Cr   a particular pair if no rules meet this criterion.  Also, pairwise
Cr   exchange interactions connecting the i-j pair implicitly defines
Cr   the both the j-i interaction at the same time.
Cr  *For each site, the following are loaded into the variables table:
Cr     is1:    species index for the site
Cr     zs1:    atomic number for the site
Cr     ib1:    site index
Cr  *For each pair, the following are loaded into the variables table:
Cr     is1,is2,zs1,zs2: species index and atomic number
Cr     ib1,ib2:         site index
Cr     d:               pair separation length
Cr  *Example: suppose there is a single one-center hamiltonian and
Cr            a single two-center hamiltonian, and hamrul is:
Cr    'h1: -3 zs1==8 h2: .1 zs1==8&zs2==8 h2: .2 is1==2&ib2>=6'
Cr   All sites with atomic number z=8 have a one-center interaction
Cr   coefficient of -3; there are no other 1-center interactions.
Cr   All pairs both sites having atomic number z=8 are assigned a
Cr   coefficient .1.  All remaining pairs with the first site belong
Cr   to species 2 and the second site having index 6 are assigned a
Cr   coefficient 0.2.  All other interactions are zero.
Cl Local variables
Cl   Program works internally by looping over rules, and for each
Cl   rule looping over any pairs which have no interaction defined.
Cl   It breaks hamrul into substrings, calling 'word' and 'nword'.
Cl  *jtypj is an index to what kind of interaction the current rule
Cr   refers to (see strings htype).
Cl  *iw1 holds the index to the string in hamrul currently being parsed.
Cl  *iwj holds indices to starting and ending substrings for each
Cl   expression of the current rule.
Cl  *iw2 counts how many expr are associated with the current rule.
Cl  *i1,i2 hold positions of current word being parsed in hamrul.
Cl  *nrule is a counter for the number of rules parsed so far.
Cl  *ib1,ib2,is1,is2 are the site and class indices of the current pair.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,ntab(1),niax,nttab,ips(1),
     .isite(nbas),ipair(nttab),ltype(2),ntype(2)
      parameter (niax=10)
      integer iax(niax,1)
      real(8):: rtab(3,1) , esite(nbas) , epair(nttab)
      type(s_spec)::sspec(*)

      character hamrul*(*),htype(2)*5
C Local parameters
      integer mxexpr
      character*80 outs,strn*120,alabl*8,alabl2*8
      double precision d,eint
      parameter (mxexpr=10)
      integer i,i1,i1mach,i2,ib1,ib2,is1,is2,ipr,iprint,ishell,iv0,ival,
     .iw,iw1,iw2,iwj(2,mxexpr),j,j1,j2,jtypj,lgunit,nrule,nshell,nw,
     .n1,n2,n12,ih1,ih2
      double precision dshell,dgetss,zs1,zs2
      logical a2bin,lsw,lquit

C --- Setup ---
      ipr = iprint()
      call words(hamrul, nw)
      call word(hamrul,nw,i1,i2)
      if (nw .le. 1) call rx('suemph: missing hamiltonian')
      nw = nw-1
      call word(hamrul,1,i1,i2)
      iw1 = 0
      nrule = 0
      call iinit(ipair,nttab)
      call dpzero(epair,nttab)
      n1 = ntype(1)
      n2 = ntype(2)
      n12 = n1+n2
      if (n1 .gt. 1 .or. n2 .gt. 1) call rx('suemph: not implemented')

C --- Get the sequence of expressions for next hamiltonian rule ---
    2 continue
      iw1 = iw1+1
C   ... There need to be at least two more words
      if (iw1+1 .gt. nw) goto 60
      nrule = nrule+1
      i1 = i2+1
      call word(hamrul,iw1,i1,i2)
C        print *, hamrul(i1:i2)
      call tokmat(hamrul(i1:i2),htype,2,5,' ',jtypj,i,.false.)
      jtypj = jtypj+1
C   ... Error if we don't recognize the interaction
      if (jtypj .le. 0) then
        outs = 'unrecognized interaction "'//hamrul(i1:i2)//'"'
        call rx('suemph: '//outs)
      endif
      if (jtypj .gt. 2) then
        outs = 'interaction "'//hamrul(i1:i2)//'" not implemented'
        call rx('suemph: '//outs)
      endif
C   ... Find the word specifying magnitude of the interaction
      iw1 = iw1+1
      i1 = i2+1
      call nword(hamrul,1,i1,i2)
      ih1 = i1
      ih2 = i2
C       This generates interaction, but it may depend on z,d, etc,
C       so we postpone evaluating until we have loaded variables
C        j = i1-1
C        if (.not. a2bin(hamrul,eint,4,0,' ',j,i2)) then
C          outs = 'suemph: failed to parse "'//hamrul(i1:i2)//'"'
C          call rx(outs)
C        endif

C  ---  Find all subsequent words associated with this interaction ---
      iw2 = 0
      do  4  iw = 1, nw-iw1+1
        i1 = i2+1
        call nword(hamrul,1,i1,i2)
        if (iw .gt. mxexpr) call rx('suemph: too many expr')
        iwj(1,iw) = i1
        iwj(2,iw) = i2
C     ... Quit if start of new interaction
        call tokmat(hamrul(i1:i2),htype,2,5,' ',j,i,.false.)
        if (j .ge. 0) goto 6
        iw2 = iw2+1
    4 continue
    6 continue

C   --- On-site: assign this hamiltonian to sites matching rule ---
      call numsyv(iv0)
      if (jtypj .le. n1 .and. ltype(jtypj) .ne. 0) then
        do  10  ib1 = 1, nbas
C     ... The interaction is already defined for this site
          if (isite(ib1) .ne. 0) goto 10
          is1 = ips(ib1)
          zs1 = (sspec(is1)%z)

          call lodsyv('ib1',1,dble(ib1),ival)
          call lodsyv('is1',1,dble(is1),ival)
          call lodsyv('zs1',1,zs1,ival)
C         call shosyv(0,0,0,6)
C     ... New interaction if any expression is nonzero
          do  12  iw = 1, iw2
            j = 0
            lsw = a2bin(hamrul(iwj(1,iw):iwj(2,iw)),
     .      ival,2,0,' ',j,iwj(2,iw)-iwj(1,iw)+1)
            if (lsw .and. ival .ne. 0) then
              isite(ib1) = nrule
              j = ih1-1
              if (.not. a2bin(hamrul,eint,4,0,' ',j,ih2)) then
                outs = 'suemph: failed to parse "'//hamrul(ih1:ih2)//'"'
                call rx(outs)
              endif
              esite(ib1) = eint
            endif
   12     continue
   10   continue

C   --- Pairs: assign this hamiltonian to pairs matching rule ---
      elseif (jtypj .gt. n1 .and. jtypj .le. n12
     ..and. ltype(jtypj) .ne. 0) then
        do  20  i = 1, nttab
C     ... The interaction is already defined for this pair
          if (ipair(i) .ne. 0) goto 20
C     ... Load the variables table with pair-specific information
          ib1 = iax(1,i)
          ib2 = iax(2,i)
          d   = dsqrt(rtab(1,i)**2 + rtab(2,i)**2 + rtab(3,i)**2)
          if (d .lt. 1d-6) goto 20
          is1 = ips(ib1)
          is2 = ips(ib2)
          zs1 = (sspec(is1)%z)

          zs2 = (sspec(is2)%z)

          call lodsyv('ib1',1,dble(ib1),ival)
          call lodsyv('ib2',1,dble(ib2),ival)
          call lodsyv('is1',1,dble(is1),ival)
          call lodsyv('is2',1,dble(is2),ival)
          call lodsyv('zs1',1,zs1,ival)
          call lodsyv('zs2',1,zs2,ival)
          call lodsyv('d',1,d,ival)
C         call shosyv(0,0,0,6)
C     ... New interaction if any expression is nonzero
          do  22  iw = 1, iw2
            j = 0
            lsw = a2bin(hamrul(iwj(1,iw):iwj(2,iw)),
     .      ival,2,0,' ',j,iwj(2,iw)-iwj(1,iw)+1)
            if (lsw .and. ival .ne. 0) then
              j = ih1-1
              if (.not. a2bin(hamrul,eint,4,0,' ',j,ih2)) then
                outs = 'suemph: failed to parse "'//hamrul(ih1:ih2)//'"'
                call rx(outs)
              endif
              ipair(i) = nrule
              epair(i) = eint
              j = iax(6,i)
              epair(j) = eint
              ipair(j) = nrule
            endif
   22     continue
   20   continue
      endif
      call clrsyv(iv0)
      iw1 = iw1+iw2

C --- Look for next rule ---
      goto 2

C --- Cleanup and printout ---
   60 continue

      lquit = .false.
      if (ltype(1) .ne. 0) then
        j1 = 0
        do  61  i = 1, nbas
          if (isite(i) .eq. 0) j1=j1+1
   61   continue
        lquit = lquit .or. j1 .gt. 0 .and. ltype(1) .lt. 0
      else
        j1 = 0
      endif
      if (ltype(n1+1) .ne. 0) then
        j2 = -nbas
        do  62  i = 1, nttab
          if (ipair(i) .eq. 0) j2=j2+1
   62   continue
        lquit = lquit .or. j2 .gt. 0 .and. ltype(n1+1) .lt. 0
      else
        j2 = 0
      endif
      call awrit5('%x suemph:  %i rule%-1j%?#n>1#s##,'//
     .'  %i sites%?#n# (missing h1 for %-1j%i)##,'//
     .'  %i pairs%?#n# (missing h2 for %-1j%i)##',
     .outs,80,-lgunit(2),nrule,nbas,j1,nttab-nbas,j2)
      if (ipr .ge. 20) then
        print *
        call awrit0('%a',outs,-80,-lgunit(1))
      endif

      if (ltype(1) .ne. 0) then
        if (ipr .ge. 40) then
          print 331
  331     format(/' One-center hamiltonian:'/'  site           h1')
          do  114  ib1 = 1, nbas
            do i_spacks=ips(ib1),ips(ib1)
              call spacks_copy('u',sspec(i_spacks)%name,ips(ib1),ips(ib1),alabl,i_spacks)
            enddo

            call awrit3('%x%,4i('//alabl//'%a)%10p%;12,6D(%i)',
     .      outs,80,0,ib1,esite(ib1),isite(ib1))
            call locase(outs)
            call awrit0('%a',outs,-80,-i1mach(2))
  114     continue
        endif
      endif

      if (ltype(n1+1) .ne. 0) then
        if (ipr .ge. 40) then
          print 332
  332     format(/' Two-center hamiltonian:')
          if (ipr .ge. 45) print 333
  333     format('  pair   ib        jb          d          h2')
          do  124  ib1 = 1, nbas
            do i_spacks=ips(ib1),ips(ib1)
              call spacks_copy('u',sspec(i_spacks)%name,ips(ib1),ips(ib1),alabl,i_spacks)
            enddo

            if (ipr .ge. 40) then
              j = 0
              dshell = 0
              nshell = 0
              ishell = 1
              strn = ' dist(n) ='
              do  126  i = ntab(ib1)+2, ntab(ib1+1)
                d   = dsqrt(rtab(1,i)**2 + rtab(2,i)**2 + rtab(3,i)**2)
                if (dabs(d-dshell) .gt. 1d-6) then
                  nshell = nshell+1
                  if (nshell .gt. 1) call awrit2('%a %,3;3d(%i)',
     .            strn,len(strn),0,dshell,i-ishell)
                  ishell = i
                  dshell = d
                endif
                if (ipair(i) .eq. 0) j=j+1
  126         continue
              if (nshell .gt. 1) call awrit2('%a %,3;3d(%i)',
     .        strn,len(strn),0,dshell,i-ishell)
              call awrit5('%x site %i('//alabl//'%a)  %i pairs  %i shells'//
     .        '%?#n# (J1 missing in %i pairs)#%j#',outs,80,0,ib1,
     .        (ntab(ib1+1)-ntab(ib1)-1),nshell,j,j)
              call awrit0('%a',outs,-80,-i1mach(2))
              if (ipr .gt. 40) call awrit0('%a',strn,-80,-i1mach(2))
            endif
            if (ipr .ge. 45) then
              do  128  i = ntab(ib1)+1, ntab(ib1+1)
                ib2 = iax(2,i)
                do i_spacks=ips(ib2),ips(ib2)
                  call spacks_copy('u',sspec(i_spacks)%name,ips(ib2),ips(ib2),alabl2,i_spacks)
                enddo

                d   = dsqrt(rtab(1,i)**2 + rtab(2,i)**2 + rtab(3,i)**2)
                outs = ' '
                call awrit6('%,4i%,4i('//alabl//'%a)%14p%,5i('//alabl2//
     .          '%a)%25p%;12,6D%;12,6D(%i)',outs,80,0,i,ib1,ib2,d,epair(i),
     .          ipair(i))
                call locase(outs)
                call awrit0('%a',outs,-80,-i1mach(2))
  128         continue
            endif
  124     continue
        endif
      endif

      call rxx(lquit,'suemph: incomplete hamiltonian')

      end subroutine suemph


