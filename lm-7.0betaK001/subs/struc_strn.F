
      subroutine lstra(sname,indxs,offe,noffe)
C- Returns a list of indices corresponding to names in struc
C ----------------------------------------------------------------
Ci Inputs
Ci   sname   struc type: must be one of those listed in styps, below,
Ci           followed by a ist of elements sought in struc;
Ci           see Remarks, below
Co Outputs
Co   indxs   index to which entry in styps sname corresponds
Co   offe    list of indices in struc corresponding to names
Co   noffe   number of entries in offe
Ci   ... entry lstrs
Ci   offe(1) index in struc for which name is sought
Co   sname   element name associated with index offe(1)
Co   noffe   number of elements in structure
Ci   ... entry lstrp
Ci   offi     offset in work array to this structure
Cr Remarks
Cr   lstra makes indices associated with names to entries in
Cr   structures it knows about.  For lstra to work, the indices
Cr   elements in a given structure are sorted alphabetically by
Cr   name.  Internal tables within each structure keep track
Cr   of the offset within the structure.
Cr   NB: the structures named here should be consistent with
Cr       those in spack and spackv, below.
Cr
Cr  *To create a new structure (e.g. 'ugw')
Cr     1.  Add name to 'styps'.  Must add new nlist,list,lists.
Cr         list must be alphabetically ordered,
Cr         so you must shift all indices list[k,...] to
Cr         list[k+1,...] in the data statements fo list[k]
Cr         and similarly parameter statements in nlist[k]
Cr         Update nlists; add entry in data statement for 'offs'
Cr     2.  Create a new subroutine, e.g. 'ugw'  The number of
Cr         elements defined there must match nlist[k]
Cr     3.  Create a new data statement for list[k]
Cr     4.  Add a new entry in spack.  Again, shift the
Cr         computed 'goto' statements for k+1...
Cr         If the structure has a species index, first argument
Cr         in 'call pack(...') should be 10+lpack, otherwise
Cr         just lpack.
Cr     5.  In routine spackv, make sure that the indices
Cr         corresponding to those structures with species
Cr         (eg uspec, usite) correspond with those in spack,
Cr         and if your new struc has species, add lines there.
Cu Updates
Cu   29 Jul 08 New ctrl->pfloat
Cu   05 Jul 08 New ham->(lmxax,ndham,npwmin,npwpad,pwemax,pwemin,pwmode,oveps)
Cu   27 Jun 08 New str->rmaxg,ivl spec->ehvl
Cu   07 Dec 07 Add sham->qpoff
Cu   17 Nov 07 LDA+U parameters sham->nlibu,lmaxu,udiag; remove sctrl->nlibu
Cu   08 Jun 07 New stb and sspec parameters for JEK's additions to tbe
Cu   01 May 07 New pot->opdel
Cu   27 Mar 07 New ctrl->nlibu
Cu   27 Jun 06 New sham->pmax
Cu   21 Mar 06 New ctrl->defm
Cu   01 Feb 06 New spec->rcfa
Cu   21 Dec 05 Increase spec->coreq to length 2
Cu    1 Aug 05 Parameters for ES finder
Cu   27 Apr 05 (WRL) parameters for LDA+U
Cu   12 Sep 04 Increment size of bz->semsh
Cu   25 Jun 04 Add uspec->eh3,rs3,vmtz
Cu   19 May 04 Add sstr->loka
Cu   13 May 04 Add spec->(colxbs,lxi,radxbs,rcut,rint)
Cu             Add sham->alfsi,dabc
Cu             Add ssite->mpole,dpole
Cu             add str->nbisi,nalf,ncupl,ndust,adec,wx
Cu    6 Apr 04 Add pot->nrhos and pot->obxc
Cu    5 Mar 04 Enlarge sham->eterms to size 20
Cu    2 Mar 04 Added rpad to LAT structure
Cu    6 Jan 04 Added rsstol to HAM structure
Cu   24 Jun 03 Added ecuts to GW structure
Cu   24 Jun 03 (Chantis) relativistic arrays in spot
Cu   21 May 03 Added sham->sigp
Cu   09 Mar 03 Enlarged sham->eterms
Cu   28 Feb 03 Potential function arrays for fully rel case (upot)
Cu   27 Jan 03 parms for magnetic field in (ham, site, spot) strux
Cu   22 Jul 02 More parameters to sham for r.s. sigma
Cu   22 Jul 02 Added parameters to sham and sgw for r.s. sigma
Cu   10 Apr 02 Redimensionsed sspec to accomodate larger lmax
Cu   12 Oct 00 new structure sgw.
Cu    2 Feb 00 igets returns nint(strux entries with double cast)
Cu             if specified range, lgors uses first entry
C ----------------------------------------------------------------
C     implicit none
      character*(*) sname
      integer indxs,offe(1),noffe
C Local variables
      integer nmax,ntyps,
     .  nlistf
C     ntyps = number of structures; nmax = max length of structure
      parameter (ntyps=1,nmax=9)
      parameter (nlistf=9)
      integer nlists(ntyps),offs(ntyps),i,j,offi
      character*6 lists(nmax,ntyps),styps(ntyps), listf(nlistf)
      equivalence
     .  (lists(1,1),listf)
      save styps,nlists,offs
C     styps = structure names
      data styps /'strn'/
      data offs  /0/
      data listf /'amix','gemb','gfopt','jobid','map','mix','mmham',
     .  'sxopt','symg'/
      data nlists /nlistf/

C      print 333, listd
C  333 format(a)
C      stop

      call word(sname,1,i,j)
C     Find index indxs to structure with ascii repsn snam(i:j)
      call strxls(styps,ntyps,sname(i:j),indxs,noffe)
      if (noffe.lt.0) call rxs('lstra:  unknown structure: ',sname)

      if (j+1 .lt. len(sname)) then
        call strxls(lists(1,indxs),nlists(indxs),sname(j+2:),offe,noffe)
      else
        noffe = 0
      endif

      end



      subroutine ustrn(sstrn,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for string structure
C ----------------------------------------------------------------
Ci Inputs (see Remarks)
Ci   sstrn  not used; present only for compatibility with other
Ci                    unpacking routines.
Ci   offe   list of indices offe(i) for i=1..offe (see Remarks)
Ci          offe(1) = 0 => do nothing
Ci          offe(i)>0 unpacks offset and length of string in off,nelt
Ci                    Here offe(i) are the list of i=1..noffe indices
Ci          offe(i)<0 assigns offset and length of string from off,nelt
Ci                    Here -offe(i) are the list of i=1..noffe indices
Ci   noffe  number of entries in offe
Ci   nelt   (for each offe, if offe<0) string length
Co  Outputs
Co   For each i, i=1..noffe, return the following. If offe(i) > 0:
Co     off    for each parameter in offe, returns offset to sstrn
Co     cast   for each parameter in offe, returns cast=1 of entry
Co     nelt   for each parameter in offe, returns length of string
Co   If offe(i) < 0:
Co     off    for each parameter in offe, returns offset to sstrn
Co     cast   for each parameter in offe, returns cast=1 of entry
Cr Remarks.
Cr   sstrn holds a set of starting and ending points in a large char
Cr   array (not passed here) corresponding to strings labeled below.
Cr   lstra associates string labels with index (eg 'mix' with '6')
Cr
Cr   The string structure is different from the others, because
Cr   of the screwy way fortran handles string arguments.
Cr
Cr   ustrn takes a list of noffe indices (passed in offe)
Cr   and, if offe>0
Cr       returns a list of starting positions (in off(1..noffe)) and
Cr       string lengths (in nelt(1..offe)) in the global character array
Cr       (not passed here) that consitute strings associated w/ indices
Cr   or, if offe<0
Cr       assigns offsets and lengths (in off(1..noffe) and nelt(1...noffe)
Cr       based for future storage.
Cr
Cr   Caller must actually copy the string into or out of the structure.
Cr
Cr   Example: copy string instr(i..strsiz) to 'mmham' slot
Cr   Here sstrng is name of char array holding all strings.
Cr       ... Get offs(1) = index to mmham
Cr       call lstra('strn mmham',tmp,offs,i)
Cr       ... Allocate string of length 1 to get offset into offi (i1
Cr       call ustrn(tmp,-offs,1,offi,tmp,1)
Cr       ... copy instr(i..strsiz) to sstrng, compressing spaces.
Cr       call cpstr(instr,strsiz,101,'}',i,j,sstrng(offi(1):))
Cr       ... j characters were copied; reallocate string of proper length
Cr       call ustrn(tmp,-offs,1,offi,i1,1)
Cr
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr         1   amix    amix mode
Cr         2   gemb    file names for embedding GF
Cr         3   gfopt   GF options
Cr         4   jobid   job id
Cr         5   map     mapping
Cr         6   mix     mixing mode
Cr         7   mmham   micromagnetics hamiltonian
Cr         8   sxopt   SX options
Cr         9   symg    symmetry group
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,nlists
      parameter (nlists=9)
      double precision sstrn(1)
      integer ilists(nlists),ilen(nlists),il,offmx,k
      save ilists,ilen,offmx
      data offmx /0/ ilists /nlists*0/

      if (offe(1) .eq. 0) then
C       size = noffe
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = 1
        if (il .lt. 0) then
C     ... See whether to re-allocate string.
          if (ilists(-il) .ne. 0) then
            do  12  k = 1, nlists
C       ... Branch if not the top one; make a new allocation
   12       if (ilists(-il) .lt. ilists(k)) goto 14
C       ... It is the top one; just remake ilen and offmx
            goto 16
          endif
C     ... Entry point for new string allocation
   14     ilists(-il) = offmx+1
C     ... Entry point for last string re-allocation
   16     ilen(-il)   = nelt(i)
          offmx  =  ilists(-il)+nelt(i)
          off(i) =  ilists(-il)
        else
          off(i) =  ilists(il)
          nelt(i) = ilen(il)
        endif
   10 continue

      end


      subroutine upacks(sname,is1,is2)
C- Unpacks one entry in string structure
C     implicit none
      character*(*) sname
      integer is1,is2,i,offs,j,xx,leni,casti

      call lstra(sname,i,offs,j)
      if (j .eq. -1)
     .  call rxs('upacks: unrecognized element "',sname)
      call ustrn(xx,offs,1,is1,casti,leni)
      is2 = is1 + leni-1

      end
