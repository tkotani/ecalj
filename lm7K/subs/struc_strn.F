      subroutine packs(snamein,sstrn,packingw)
C pack paking word(packingw) into sstrn. Use upacks to unpack.
Ci snamein: sname defined in sstrn_sname.
Ci sstrn: long string. sstrings.
Ci packingw 
      use m_sstrnsname
      character*(*):: snamein,sstrn,packingw
      integer:: i,lll
      if(initialized==0) then
        call initialize()
        initialized=1
      endif
      do i=1,nsname
        lll=len_trim(packingw) 
        if(lll>nx) then
          call rxs('packs: too large packing word.'//
     &              'Need to modify code',packingw)
        endif
        if(trim(sname(i))==trim(snamein)) then
          iend(i) = init(i) +lll -1
          sstrn(init(i):iend(i))= packingw
          return
        endif
      enddo
      call rxs('packs: no sname is found in sstrn_sname. '//
     &          'bug or need to modify codes',snamein)
      end
      
      subroutine upacks(snamein,init_,iend_)
c return init and iend for sstring(init:iend) for given snamein
      use m_sstrnsname
      integer::init_,iend_
      character*(*):: snamein
      if(initialized==0) then
        call initialize()
        initialized=1
      endif
      do i=1,nsname
        if(trim(sname(i))==trim(snamein)) then
          init_=init(i)
          iend_=iend(i)
          return
        endif
      enddo
      call rxs('upacks: no sname is found for snamein. '//
     &          'bug or need to modify codes',snamein)
      end

C$$$      subroutine lstra(sname,indxs,offe,noffe)
C$$$C- Returns a list of indices corresponding to names in struc
C$$$C ----------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   sname   struc type: must be one of those listed in styps, below,
C$$$Ci           followed by a ist of elements sought in struc;
C$$$Ci           see Remarks, below
C$$$Co Outputs
C$$$Co   indxs   index to which entry in styps sname corresponds
C$$$Co   offe    list of indices in struc corresponding to names
C$$$Co   noffe   number of entries in offe
C$$$Ci   ... entry lstrs
C$$$Ci   offe(1) index in struc for which name is sought
C$$$Co   sname   element name associated with index offe(1)
C$$$Co   noffe   number of elements in structure
C$$$Ci   ... entry lstrp
C$$$Ci   offi     offset in work array to this structure
C$$$Cr Remarks
C$$$Cr   lstra makes indices associated with names to entries in
C$$$Cr   structures it knows about.  For lstra to work, the indices
C$$$Cr   elements in a given structure are sorted alphabetically by
C$$$Cr   name.  Internal tables within each structure keep track
C$$$Cr   of the offset within the structure.
C$$$Cr   NB: the structures named here should be consistent with
C$$$Cr       those in spack and spackv, below.
C$$$Cr
C$$$Cr  *To create a new structure (e.g. 'ugw')
C$$$Cr     1.  Add name to 'styps'.  Must add new nlist,list,lists.
C$$$Cr         list must be alphabetically ordered,
C$$$Cr         so you must shift all indices list[k,...] to
C$$$Cr         list[k+1,...] in the data statements fo list[k]
C$$$Cr         and similarly parameter statements in nlist[k]
C$$$Cr         Update nlists; add entry in data statement for 'offs'
C$$$Cr     2.  Create a new subroutine, e.g. 'ugw'  The number of
C$$$Cr         elements defined there must match nlist[k]
C$$$Cr     3.  Create a new data statement for list[k]
C$$$Cr     4.  Add a new entry in spack.  Again, shift the
C$$$Cr         computed 'goto' statements for k+1...
C$$$Cr         If the structure has a species index, first argument
C$$$Cr         in 'call pack(...') should be 10+lpack, otherwise
C$$$Cr         just lpack.
C$$$Cr     5.  In routine spackv, make sure that the indices
C$$$Cr         corresponding to those structures with species
C$$$Cr         (eg uspec, usite) correspond with those in spack,
C$$$Cr         and if your new struc has species, add lines there.
C$$$Cu Updates
C$$$Cu   29 Jul 08 New ctrl->pfloat
C$$$Cu   05 Jul 08 New ham->(lmxax,ndham,npwmin,npwpad,pwemax,pwemin,pwmode,oveps)
C$$$Cu   27 Jun 08 New str->rmaxg,ivl spec->ehvl
C$$$Cu   07 Dec 07 Add sham->qpoff
C$$$Cu   17 Nov 07 LDA+U parameters sham->nlibu,lmaxu,udiag; remove sctrl->nlibu
C$$$Cu   08 Jun 07 New stb and sspec parameters for JEK's additions to tbe
C$$$Cu   01 May 07 New pot->opdel
C$$$Cu   27 Mar 07 New ctrl->nlibu
C$$$Cu   27 Jun 06 New sham->pmax
C$$$Cu   21 Mar 06 New ctrl->defm
C$$$Cu   01 Feb 06 New spec->rcfa
C$$$Cu   21 Dec 05 Increase spec->coreq to length 2
C$$$Cu    1 Aug 05 Parameters for ES finder
C$$$Cu   27 Apr 05 (WRL) parameters for LDA+U
C$$$Cu   12 Sep 04 Increment size of bz->semsh
C$$$Cu   25 Jun 04 Add uspec->eh3,rs3,vmtz
C$$$Cu   19 May 04 Add sstr->loka
C$$$Cu   13 May 04 Add spec->(colxbs,lxi,radxbs,rcut,rint)
C$$$Cu             Add sham->alfsi,dabc
C$$$Cu             Add ssite->mpole,dpole
C$$$Cu             add str->nbisi,nalf,ncupl,ndust,adec,wx
C$$$Cu    6 Apr 04 Add pot->nrhos and pot->obxc
C$$$Cu    5 Mar 04 Enlarge sham->eterms to size 20
C$$$Cu    2 Mar 04 Added rpad to LAT structure
C$$$Cu    6 Jan 04 Added rsstol to HAM structure
C$$$Cu   24 Jun 03 Added ecuts to GW structure
C$$$Cu   24 Jun 03 (Chantis) relativistic arrays in spot
C$$$Cu   21 May 03 Added sham->sigp
C$$$Cu   09 Mar 03 Enlarged sham->eterms
C$$$Cu   28 Feb 03 Potential function arrays for fully rel case (upot)
C$$$Cu   27 Jan 03 parms for magnetic field in (ham, site, spot) strux
C$$$Cu   22 Jul 02 More parameters to sham for r.s. sigma
C$$$Cu   22 Jul 02 Added parameters to sham and sgw for r.s. sigma
C$$$Cu   10 Apr 02 Redimensionsed sspec to accomodate larger lmax
C$$$Cu   12 Oct 00 new structure sgw.
C$$$Cu    2 Feb 00 igets returns nint(strux entries with double cast)
C$$$Cu             if specified range, lgors uses first entry
C$$$C ----------------------------------------------------------------
C$$$C     implicit none
C$$$      character*(*) sname
C$$$      integer indxs,offe(1),noffe
C$$$C Local variables
C$$$      integer nmax,ntyps,
C$$$     .  nlistf
C$$$C     ntyps = number of structures; nmax = max length of structure
C$$$      parameter (ntyps=1,nmax=9)
C$$$      parameter (nlistf=9)
C$$$      integer nlists(ntyps),offs(ntyps),i,j,offi
C$$$      character*6 lists(nmax,ntyps),styps(ntyps), listf(nlistf)
C$$$      equivalence
C$$$     .  (lists(1,1),listf)
C$$$      save styps,nlists,offs
C$$$C     styps = structure names
C$$$      data styps /'strn'/
C$$$      data offs  /0/
C$$$      data listf /'amix','gemb','gfopt','jobid','map','mix','mmham',
C$$$     .  'sxopt','symg'/
C$$$      data nlists /nlistf/

C$$$C      print 333, listd
C$$$C  333 format(a)
C$$$C      stop

C$$$      call word(sname,1,i,j)
C$$$C     Find index indxs to structure with ascii repsn snam(i:j)
C$$$      call strxls(styps,ntyps,sname(i:j),indxs,noffe)
C$$$      if (noffe.lt.0) call rxs('lstra:  unknown structure: ',sname)

C$$$      if (j+1 .lt. len(sname)) then
C$$$        call strxls(lists(1,indxs),nlists(indxs),sname(j+2:),offe,noffe)
C$$$      else
C$$$        noffe = 0
C$$$      endif

C$$$      end



C$$$      subroutine ustrn(sstrn,offe,noffe,off,cast,nelt)
C$$$C- Returns offsets, casts and number for string structure
C$$$C ----------------------------------------------------------------
C$$$Ci Inputs (see Remarks)
C$$$Ci   sstrn  not used; present only for compatibility with other
C$$$Ci                    unpacking routines.
C$$$Ci   offe   list of indices offe(i) for i=1..offe (see Remarks)
C$$$Ci          offe(1) = 0 => do nothing
C$$$Ci          offe(i)>0 unpacks offset and length of string in off,nelt
C$$$Ci                    Here offe(i) are the list of i=1..noffe indices
C$$$Ci          offe(i)<0 assigns offset and length of string from off,nelt
C$$$Ci                    Here -offe(i) are the list of i=1..noffe indices
C$$$Ci   noffe  number of entries in offe
C$$$Ci   nelt   (for each offe, if offe<0) string length
C$$$Co  Outputs
C$$$Co   For each i, i=1..noffe, return the following. If offe(i) > 0:
C$$$Co     off    for each parameter in offe, returns offset to sstrn
C$$$Co     cast   for each parameter in offe, returns cast=1 of entry
C$$$Co     nelt   for each parameter in offe, returns length of string
C$$$Co   If offe(i) < 0:
C$$$Co     off    for each parameter in offe, returns offset to sstrn
C$$$Co     cast   for each parameter in offe, returns cast=1 of entry
C$$$Cr Remarks.
C$$$Cr   sstrn holds a set of starting and ending points in a large char
C$$$Cr   array (not passed here) corresponding to strings labeled below.
C$$$Cr   lstra associates string labels with index (eg 'mix' with '6')
C$$$Cr
C$$$Cr   The string structure is different from the others, because
C$$$Cr   of the screwy way fortran handles string arguments.
C$$$Cr
C$$$Cr   ustrn takes a list of noffe indices (passed in offe)
C$$$Cr   and, if offe>0
C$$$Cr       returns a list of starting positions (in off(1..noffe)) and
C$$$Cr       string lengths (in nelt(1..offe)) in the global character array
C$$$Cr       (not passed here) that consitute strings associated w/ indices
C$$$Cr   or, if offe<0
C$$$Cr       assigns offsets and lengths (in off(1..noffe) and nelt(1...noffe)
C$$$Cr       based for future storage.
C$$$Cr
C$$$Cr   Caller must actually copy the string into or out of the structure.
C$$$Cr
C$$$Cr   Example: copy string instr(i..strsiz) to 'mmham' slot
C$$$Cr   Here sstrng is name of char array holding all strings.
C$$$Cr       ... Get offs(1) = index to mmham
C$$$Cr       call lstra('strn mmham',tmp,offs,i)
C$$$Cr       ... Allocate string of length 1 to get offset into offi (i1
C$$$Cr       call ustrn(tmp,-offs,1,offi,tmp,1)
C$$$Cr       ... copy instr(i..strsiz) to sstrng, compressing spaces.
C$$$Cr       call cpstr(instr,strsiz,101,'}',i,j,sstrng(offi(1):))
C$$$Cr       ... j characters were copied; reallocate string of proper length
C$$$Cr       call ustrn(tmp,-offs,1,offi,i1,1)
C$$$Cr
C$$$Cr   The index to the elements are arranged alphabetically by name:
C$$$Cr   The following maps the internal location to the index:
C$$$Cr   off offe  name    purpose
C$$$Cr         1   amix    amix mode
C$$$Cr         2   gemb    file names for embedding GF
C$$$Cr         3   gfopt   GF options
C$$$Cr         4   jobid   job id
C$$$Cr         5   map     mapping
C$$$Cr         6   mix     mixing mode
C$$$Cr         7   mmham   micromagnetics hamiltonian
C$$$Cr         8   sxopt   SX options
C$$$Cr         9   symg    symmetry group
C$$$C ----------------------------------------------------------------
C$$$C     implicit none
C$$$      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C$$$C Local variables
C$$$      integer i,nlists
C$$$      parameter (nlists=9)
C$$$      double precision sstrn(1)
C$$$      integer ilists(nlists),ilen(nlists),il,offmx,k
C$$$      save ilists,ilen,offmx
C$$$      data offmx /0/ ilists /nlists*0/

C$$$      if (offe(1) .eq. 0) then
C$$$C       size = noffe
C$$$        return
C$$$      endif

C$$$      do  10  i = 1, noffe
C$$$        il = offe(i)
C$$$        cast(i) = 1
C$$$        if (il .lt. 0) then
C$$$C     ... See whether to re-allocate string.
C$$$          if (ilists(-il) .ne. 0) then
C$$$            do  12  k = 1, nlists
C$$$C       ... Branch if not the top one; make a new allocation
C$$$   12       if (ilists(-il) .lt. ilists(k)) goto 14
C$$$C       ... It is the top one; just remake ilen and offmx
C$$$            goto 16
C$$$          endif
C$$$C     ... Entry point for new string allocation
C$$$   14     ilists(-il) = offmx+1
C$$$C     ... Entry point for last string re-allocation
C$$$   16     ilen(-il)   = nelt(i)
C$$$          offmx  =  ilists(-il)+nelt(i)
C$$$          off(i) =  ilists(-il)
C$$$        else
C$$$          off(i) =  ilists(il)
C$$$          nelt(i) = ilen(il)
C$$$        endif
C$$$   10 continue

C$$$      end


C$$$      subroutine upacks(sname,is1,is2)
C$$$C- Unpacks one entry in string structure
C$$$C     implicit none
C$$$      character*(*) sname
C$$$      integer is1,is2,i,offs,j,xx,leni,casti

C$$$      call lstra(sname,i,offs,j)
C$$$      if (j .eq. -1)
C$$$     .  call rxs('upacks: unrecognized element "',sname)
C$$$      call ustrn(xx,offs,1,is1,casti,leni)
C$$$      is2 = is1 + leni-1

C$$$      end
