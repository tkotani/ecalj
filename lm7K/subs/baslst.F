Cdelw1 var ['olist', '1', 'olist', 'nspec', 'integer', 'iv_w_', 'baslst']
Cdelw1 undel [['wref'], 'olist', 'baslst']
      subroutine baslst(lstyld,iopt,slst,ilast,ips,nbas,slabl,z,nopt,
     .  optlst,optarg,nlist,list)
C- Generates a list of sites from a string specification
C ----------------------------------------------------------------
Ci Inputs
Ci   lstyld:default list style
Ci   iopt  :1s digit
Ci         : 1  print out list
Ci         :10s digit
Ci         : 0, do not sort list
Ci         : 1, sort list, paring duplicates
Ci   slst  :string bearing information for list
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   nbas  :size of basis
Ci   z     :Table of nuclear charges by species
Ci     ... baslst is designed to convert switches (substring of slst)
Ci         into numbers; not implemented
Ci   nopt  :number of special modifiers
Ci   optlst:string of nopt tokens for modifiers
Co Outputs
Co   nlist :number of elements in list
Co   list  :list of sites
Co   ilast :index to last char in string parsed by baslst
Co   optarg:numbers possibly generated from optional switches
Cr Remarks
Cu Updates
Cu   07 Apr 03 New arguments lstyld and ilast
Cu   13 Sep 01 Supersedes old optget.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lstyld,iopt,nbas,ips(1),nlist,list(*),nopt,ilast
      character slst*(*),optlst(nopt)*(*)
      character*8 slabl(*)
      double precision z(*),optarg(1)
C ... Local parameters
Cdelw1 do not change  [['wref'], 'olist', 'baslst']
      integer i,j,j1,j2,ls,m,lstyle,k1,k2,
     .  iv,parg,olist,nlists,mxint,nspec,ib,is,iclbsj
      character dc*1
      integer lgunit
C ... Heap
      integer w(1)
      common /w/ w
C ... External calls
      external awrit2,defi,ishell,mkilst,nwordg,rxs,slist,word

      nlist = 0
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [olist]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (slst .eq. ' ') return
      ls = len(slst)
      j1 = 1
      dc = slst(j1:j1)
      j1 = j1+1

C ... Return here to resume parsing for arguments
C  40 continue
      call nwordg(slst,0,dc//' ',1,j1,j2)

C ... Parse special arguments (not implemented)
      if (slst(j2+1:j2+1) .ne. ' ')  then
        do  8  i = 1, nopt
          call word(optlst(i),1,k1,k2)
          m = j1-1
          j = parg(optlst(i)(k1:k2),4,slst,m,ls,dc,1,1,iv,optarg(i))
    8   continue
      endif

C --- Parse style for specifiying sites list ---
      lstyle = lstyld
      m = j1-1
      i = parg('style=',2,slst,m,ls,dc,1,1,iv,lstyle)
      if (i .gt. 0) then
        j1 = m+2
        call nwordg(slst,0,dc//' ',1,j1,j2)
      endif

C --- List of sites ---
      if (lstyle .gt. 1) then
        nspec = mxint(nbas,ips)
Cdelw1 do not change olist because of ['wref']
        call defi(olist, nspec)
        if (lstyle .eq. 2) j2=j2+1
Cdelw1 do not change olist because of ['wref']
        call slist(lstyle,slst(j1:j2),slabl,z,nspec,nlists,w(olist))
        if (lstyle .ne. 2) j2=j2+1
        nlist = 0
        do  i = 1, nlists
Cdelw1 1 is digit, dropped
Cdelw1 do not change olist because of ['wref']
          is = w(olist+i-1)
          do  j = 1, nbas
            ib = iclbsj(is,ips,-nbas,j)
            if (ib .ge. 0) then
              nlist = nlist+1
              list(nlist) = ib
            endif
          enddo
        enddo
        if (nlist .eq. 0 .or. mod(iopt/10,10) .eq. 0) goto 99
        call ishell(nlist,list)
C        j = 1
C        do  i = 2, nlist
C          if (list(i) .gt. list(j) .and. list(i) .le. nbas) then
C            j = j+1
C            list(j) = list(i)
C          endif
C        enddo
C        nlist = j
C        j2=j2+1

      elseif (slst(j1:j1+1) .eq. 'z ' .or.
     .        slst(j1:j1+1) .eq. 'Z ') then
        nlist = 0
        do  10  ib = 1, nbas
          is = ips(ib)
          if (z(is) .eq. 0) then
            nlist = nlist+1
            list(nlist) = ib
          endif
   10   continue

        j2=j2+1

      else
        j = j1
        call nwordg(slst,0,dc//' ',1,j,j2)
        call mkilst(slst(j1:j2),nlist,list)
        if (nlist .lt. 0)
     .    call rxs('baslst: failed to parse site list ',slst(j1:))
        if (nlist .eq. 0 .or. mod(iopt/10,10) .eq. 0) goto 99
        call ishell(nlist,list)
        j = 1
        do  16  i = 2, nlist
          if (list(i) .gt. list(j) .and. list(i) .le. nbas) then
            j = j+1
            list(j) = list(i)
          endif
   16   continue
        nlist = j
        j2=j2+1
      endif

      ilast = j2
Cki#error, have return with len(w_varlist)>0 at line 139
Cdelw1 w_varlist remains: olist 
Cdelw1 w_varlistundel: [olist]
Cdelw1 w_varlist (undel), remains: [olist]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate olist because of [wref]
   99 if (mod(iopt,10) .eq. 0) return
      call awrit2(' baslst: %n:1i',slst,255,lgunit(1),nlist,list)
Cdelw1 w_varlist remains: olist 
Cdelw1 w_varlistundel: [olist]
Cdelw1 w_varlist (undel), remains: [olist]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate olist because of [wref]
      end

