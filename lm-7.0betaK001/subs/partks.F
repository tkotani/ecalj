      integer function prtkss(categ,token,term,sname,sstrn,mxlen,chr,
     .                  istart,iopt)
C- Read into a string token into sstrn
Ci mxlen: maximum length to allow for string
C     implicit none
      integer istart,iopt,mxlen
      character*(*) sname,categ(0:1)*1,term,chr,token,sstrn

C Local variables
      double precision xx
C iolib common block
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      integer offs,j,i,offi,casti,leni,partok
      call lstra(sname,i,offs,j)
      if (j .eq. -1) call rxs('prtkss: unrecognized element "',sname)

C ... Pick up, or allocate initial space string
      if (optio .ne. 2) call ustrn(xx,-offs,1,offi,casti,mxlen)
      if (optio .eq. 2) call ustrn(xx,offs,1,offi,casti,leni)

      prtkss = partok(categ,token,term,mxlen,sstrn(offi:),1,1,istart,
     .  iopt)

      if (optio .eq. 1) then
        call skpblb(sstrn(offi:),iabs(mxlen),leni)
        call ustrn(xx,-offs,1,offi,casti,leni+1)
      endif
      end
      integer function partks(categ,token,term,sname,struc,
     .  is,chr,count,cast,istart,iopt)
C- Read into structure a vector of numbers parsed from a token
C ----------------------------------------------------------------
Ci Inputs
Ci   sname  'struc-name elt[?$|?val][;alias][,mask][:range[.number]]
Ci           See salias for syntax of modifiers.
Ci           ?$ => whether element should be read is conditional on
Ci                 whether prior element is read
Ci   chr
Ci   struc  structure array
Ci  Rules for termination of a token
Ci
Ci  1. It is enclosed in quotes '...'
Ci
Ci  2. The terminator is the last element in term
Ci   is:    0 if struc has no species index, else species index
Ci   categ: string to parse for token.
Ci          Parsing is between istart and internal variable subsiz.
Ci   token: token, eg 'nabc='
Ci   term:  term(1:1) = last char of token.
Ci          If len(term)=1, expressions are terminated by spaces.
Ci          If len(term)>1, partks will use any of the term(2:*)
Ci          the separator in parsing the list of elements; see Remarks
Ci   count: maximum number of elements to read 
Ci          count<0 => partks will read exactly -count elements
Ci                     or abort w/ error)
Ci          count=0 => see Outputs, result
Ci   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
Ci   istart:offset to first character of category where search begins
Ci   In the case of a character field, result should be the max size.
Ci   optio: ones digit
Ci          0: show token that would have been sought
Ci          1: attempt to find token and read contents
Ci          2: print out token and contents
Ci          tens digit
Ci          0: require ' ' preceding token
Ci          1: allow any of term to precede token
Ci          2: allow any of term to terminate token, in addition
Ci             to being separators of element list
Ci          3: combination of 1 and 2
Ci   iopt  (tells partks whether token is optional, required, for ignored)
Ci          0  token is optional
Ci          1  token is required
Ci          2, do nothing but set noerr to F, return partks=-1
Co Outputs
Co   res:array into which elements are read.
Co          if count=0 partks does parse for elements following.
Ci          if also cast=0 partks returns noerr in res
Co   noerr: see remarks
Co    iend: offset to first character past match
Co  partks: number of elements actually read (-1 if no attempt to read
Co          an element -- if optio is 0 or 2 or iopt is 2.)
Cr Remarks
Cr   on return, noerr is always false unless a token was actually
Cr   sought (optio=1) and found (note that this differs from getcat).
Cr   Passing common-block variable iend in place of istart has the
Cr   effect of updating istart to point beyond the end of the match.
Cr
Cr   ichoos is automatically incremented on each call if nchoos .ne. 0
Cr   ichoos and nchoos are automatically reset to 0 if a match is found
Cr   or if ichoos=nchoos
Cr
Cr   Caution when using multiple-terminators ie len(term)>0:
Cr   partks attempts to parse the expression for each term(2:len(term)),
Cr   and will use any one that succeeds.  To avoid ambiguity, DO NOT use
Cr   terminators that may be part of an allowed arithmetic expression.
Cu Updates
Cu  2 Feb 00  handle setting first element d.p. entry to logical
C ----------------------------------------------------------------
C     implicit none
      integer cast,is,count,istart,iopt
      character*(*) sname,categ(0:*)*1,term,chr,token
      double precision struc(1)

C Local variables
      integer procid,master,mpipid
      logical lres(100),lgors,last,pvpar2
      integer i,j,k,partok,off,i1mach,j1,j2,nw,i1,i2,strsiz,mask(10),
     .  a2vec,is1,is2,bitlow,bitor,bitand,k1,k2
      integer offi(10),casti(10),nelti(10),ix(100),ires(100),range(2,10)
      character*1 strn*1024,strn2*80,strn3*80
      character*8 ename*80,tsname*80,alias*80,switch*80
      double precision xx
C iolib common block
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      master = 0
      procid = mpipid(1)

      if (iopt .eq. 2) then
        noerr = .false.
        partks = -1
        return
      endif

C ... Cull true struc name, and names for purposes of input
      call salias(sname,tsname,alias,switch,range,mask)
C      print *, 'sname ',sname
C      print *, 'tsname ',tsname
C      print *, 'mask=',mask
C      print *, 'alias ',alias
C      stop

C ... Pick up offsets in structure
      call spacki(2,tsname,struc,is,offi,casti,nelti,i,i)

C ... Print out information sought
      if (mod(optio,10) .eq. 0 .and. cast .lt. 5) then
        off = offi(1)
        i=partok(categ,token,term,struc(off),chr,count,cast,istart,iopt)

C --- Get contents of token, put into appropriate entry ---
      else
        if (cast .eq. 0) then
          lres(1) = lgors(sname,struc)
          partks =
     .      partok(categ,token,term,lres,chr,count,cast,istart,iopt)
          if (mod(optio,10) .eq. 1) then
C           No range specified; use lsets to set
            if (range(1,1) .eq. -1) then
              call lsets(tsname,struc,lres,mask)
C           Range of length 1 specified; do it here
            elseif (range(1,1) .eq. range(2,1)) then
              mask(1) = bitlow(mask(1))
              i = offi(1) + range(1,1) - 1
              if (lres(1)) then
                struc(i) = bitor(int(struc(i)),mask)
              else
                struc(i) = int(struc(i)) - bitand(int(struc(i)),mask)
              endif
            else
              call rx('partks not ready to set multiple log bits')
            endif
         endif
        elseif (cast .eq. 1) then
          call spacks(0,tsname,struc,strn,is,is)
C     ... Assume 8 characters fit into strn
          partks =
     .      partok(categ,token,term,8,strn,count,cast,istart,iopt)
          call spacks(1,tsname,struc,strn,is,is)
C   ... Poke doubles directly into struc
        elseif (cast .lt. 5) then
          off = offi(1) + nint(struc(1))*max(is-1,0)
          call word(tsname,1,is1,is2)
          call word(switch,1,j1,j2)
          if (pvpar2(tsname(is1:is2+1),switch(j1:j2),
     .      range,.true.,is,struc,mask,offi,casti,nelti,off)) then
            if (mask(1) .ne. 0) then
              xx = struc(off)
              if (mod(optio,10) .eq. 1) xx = xx - bitand(int(xx),mask)
              if (mod(optio,10) .eq. 2) struc(off) =bitand(int(xx),mask)
            endif
            partks =
     .      partok(categ,token,term,struc(off),chr,count,4,istart,iopt)
            if (mask(1) .ne. 0) then
              if (mod(optio,10) .eq. 1) then
                struc(off) = bitand(int(struc(off)),mask) + xx
              else
                struc(off) = xx
              endif
            endif
          endif
        elseif (cast .eq. 5) then
          last = .false.
          call word(tsname,1,is1,is2)
          if (mod(optio,10) .eq. 0) then
            if (procid .eq. master) then
            call word(token,1,j1,j2)
            strn = '   token  ' // token(j1:j2) //
     .             ' of cast struc, elements: ' // alias
            call awrit1('%a%?;n==0; (optional);;',strn,len(strn),
     .        -i1mach(2),iopt)
            strn = '          '//chr//'%a'
            call awrit0(strn,' ',-len(strn),i1mach(2))
            endif
          else
            call words(alias,nw)
            strn = '   ' // token
C       ... Display
            partks = -1
            if (optio .eq. 2) then
            if (procid .eq. master) then
              do  20  i = 1, nw
C           ... Local modification of offset,range,cond. continuation
                call word(switch,i,j1,j2)
                if (.not. pvpar2(tsname(is1:is2+1),switch(j1:j2),
     .            range(1,i),last,is,struc,mask(i),offi(i),casti(i),
     .            nelti(i),off)) then
                  call awrit0('%a (*)',strn,len(strn),0)
                  goto 22
                endif
                call word(alias,i,j1,j2)
                if (i .eq. 1 .and. mask(i) .ne. 0 .and.
     .            casti(i) .eq. 2 .and. nelti(i) .eq. 1) then
                  lres(1) = lgors(sname,struc)
                  last = lres(1)
                  strn2 = '%a '//alias(j1:j2)//'=%l'
                  call awrit1(strn2,strn,len(strn),0,lres)
                else
                  strn2 = '%a '//alias(j1:j2)//'=%g'
                  do  24  j = 1, nelti(i)
                    call awrit1(strn2,strn,len(strn),0,struc(off+j-1))
                    strn2 = '%a %g'
   24             continue
                  last = int(struc(off+nelti(i)-1)) .ne. 0
                endif
*               strn2 = '%a %g'
   20         continue
   22         continue
              call awrip(j)
              print *, strn(1:j)
            endif
C       ... Read
            else
              i = partok(categ,token,term,lres,strn,0,0,istart,iopt)
              if (.not. lres(1)) return
              partks = 1
              call skipbl(categ,subsiz,iend)
              if (iend .ge. subsiz) return
              call pvpar1(categ,iend,token,term(2:),subsiz,strn,strsiz)

C         ... Now parse the input from 'strn'
              i1 = 0
              do  30  i = 1, nw
C           ... Get number of elements again, in case the number changed
                call word(tsname,i+1,j1,j2)
                ename = tsname(is1:is2+1)//tsname(j1:j2)
                call spacki(2,ename,struc,is,offi(i),casti(i),nelti(i),
     .            i,i)
C           ... Local modification of offset,range,cond. continuation
                call word(switch,i,k1,k2)
                if (.not.pvpar2(tsname(is1:is2+1),switch(k1:k2),
     .            range(1,i),last,is,struc,mask(i),offi(i),casti(i),
     .            nelti(i),off)) goto 32
                k = len(term)-1
                last = .false.
C           ... Read a logical, mask it with already existing integer
                if (casti(i) .eq. 2 .and. mask(i) .ne. 0) then
                  i2 = a2vec(strn,strsiz,i1,0,term(2:),
     .              k,-k,nelti(i),ix,lres)
                  call discop(struc(off),ires,nelti(i),1,1,0)
                  mask(i) = bitlow(mask(i))
                  if (i2 .gt. 0) last = lres(i2)
                  do  34  k = 1, nelti(i)
                    if (k .le. i2 .and. lres(k)) then
                      ires(k) = bitor(ires(k),mask(i))
                    else
                      ires(k) = ires(k) - bitand(ires(k),mask(i))
                    endif
   34             continue
                  call idscop(ires,struc(off),nelti(i),1,1)
                else
                  i2 = a2vec(strn,strsiz,i1,4,term(2:),
     .              k,-k,nelti(i),ix,struc(off))
                  if (i2 .gt. 0) last = int(struc(off+i2-1)) .ne. 0
                endif
                if (i2 .lt. nelti(i) .and. count .lt. 0) then
                  print '('' partks: parse failed at '',a,'' ... '',a)'
     .              ,strn(1:i1),strn(i1+1:min(i1+20,strsiz))
                  strn3 = ' '
C                 print 'failed to parse ... or 'sought # elements ...'
                  call awrit3(': %?;n<0;'//
     .              'failed to parse arg %j%i;'//
     .              'sought %i elements but found %i;',strn3,
     .              len(strn3),0,i2,nelti(i),iabs(i2))
                  call word(alias,i,j1,j2)
                  call rxs4('partks token ',token,' element ',
     .              alias(j1:j2),strn3)
                endif
   30         continue
   32         continue
            endif
          endif
        endif
      endif

      end

      subroutine pvpar1(categ,ic,token,term,maxsiz,strn,strsiz)
Ci ic starting position
C     implicit none
      integer ic,maxsiz,strsiz
      character categ(0:1)*1,strn*1024,ct*1,term*(*),token*(*)
      integer iend,i,j,k,is
      logical last
      iend = ic
      strn = ' '

C ... Case string in quotation marks.  Take insides of marks
      if (categ(ic) .eq. '"' .or. categ(ic) .eq. '''') then
        ct = categ(ic)
        ic = ic+1
        call strcop(strn,categ(ic),maxsiz-ic,ct,strsiz)
        if (strn(strsiz:strsiz) .ne. ct .and. strsiz .eq. maxsiz-ic)
     .    call rxs4('partks token ',token,' missing quotation "',
     .    strn(1:strsiz),'%a " ...')
        strn(strsiz:strsiz) = ' '
        iend = iend+strsiz+1

C ... Otherwise, continue until terminator is last entry in term
      else
        i = iend
        is = 1
   23   continue
        call chrps2(categ,term,len(term),maxsiz,ic,k)
        iend = ic
        strsiz = ic-i+1
        last = k .eq. len(term)
        ic = ic+1

C   ... Compress whitespace into one character
        j = i
        call skipbl(categ,maxsiz,i)
        if (i .gt. j) i = i-1

        call strncp(strn,categ,is,i+1,ic-i)
        is = is+ic-i
        i = ic
C       print *, strn(1:is)
        call skipbl(categ,maxsiz,ic)
        if (ic .lt. maxsiz .and. .not. last) goto 23
      endif
      strsiz = is
      ic = iend
C      print *, 'final strn=',strn(1:strsiz)
C      print *, categ(iend-1),categ(iend),categ(iend+1)
C      stop
      end

      logical function pvpar2(tsname,switch,range,last,is,struc,mask,
     .  offi,cast,nelti,off)
C- Make local modification of offsets, ranges
C ----------------------------------------------------------------------
Ci Inputs
Ci   tsname:structure name
Ci   switch:for conditional read/write of elements:
Ci         :'none' => there are no conditions
Ci         :'$'    => conditional on whether
Ci         :          prior element was read
Ci         :else:  entry in a structure
Ci   range :first,last entries in current element
Ci         :(applicable for vector members of a structure)
Ci         :range(2)<=0 => last entry = last vector element of member
Ci   last  :T if prior element was read; for conditional read
Ci   is    :species index, if structure has species
Ci   struc :structure
Ci   mask  :masks marking individual bits within an element
Ci   offi  :offset to first entry in structure
Co Outputs
Co   cast  :cast of elements
Co   nelti :number of elements to read
Co   off   :offset to first element in structure
Co  pvpar2 :F flag to skip reading or writing these elements; occurs
Co         :  when there is some unsatisfied condition set by switch
Co         :T flag to read or write these elements
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
      character*(*) tsname,switch
      logical last
      integer is,nelti,offi,range(2),off,cast,mask
      double precision struc(1)
      integer j,igetss
      character*40 ename

      pvpar2 = .true.

C ... Local modification of cast
      if (mask .gt. 0) cast = 2

C ... Local modification of offsets, ranges
      off = offi + nint(struc(1))*max(is-1,0)
      if (range(1) .gt. 0) then
        off = off + range(1)-1
        nelti = nelti - (range(1)-1)
        if (range(2) .ge. range(1)) nelti = range(2)-range(1)+1
      endif
C ... Check for conditional output
      if (switch .eq. '$') then
        if (.not. last) pvpar2 = .false.
      elseif (switch .ne. 'none') then
        ename = tsname//switch
        j = igetss(ename,is,struc)
        pvpar2 = j .ne. 0
      endif
      end

