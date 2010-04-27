      subroutine salias(sname,trunam,alias,switch,range,mask)

C- Unpack modifiers to structure returning true names alias, mask, range
C ----------------------------------------------------------------------
Ci Inputs
Ci   sname :structure name followed by one or more elements, e.g.
Ci         :'bz nkabc efmax'
Ci         :Each element can be accompanied by modifiers that
Ci         :specify some portion of the structure,  e.g.:
Ci         :'ctrl lncol;force,16'
Ci         ;See remarks.
Co Outputs
Co   trunam:same as sname, but without the modifiers.  e.g.
Co         :sname 'ctrl lncol;force,16' => trunam => 'ctrl lncol'
Co   alias :alias for name (see Remarks)
Co   switch:for conditional I/O of an element.  salias just returns
Co         :string following '?'
Co   range :a starting offset (and optionally upper bound) for a name
Co         :--- useful for names corresponding to vectors only
Co   mask  :masks marking individual bits within an element
Co         :--- useful for names corresponding to scalar integers only
Cr Remarks
Cr   salias unpacks modifiers to elements in a structure, to
Cr   specify some subportion of the element, or some alias
Cr   associated with the name.
Cr
Cr   The syntax of sname is 'struc-name elt1 elt2 ...'
Cr   where each of elt1 elt2 ... has the syntax with optional modifiers
Cr     element-name[?string][;alias][,mask][:range[.upper-bound]]
Cr   e.g.
Cr     'str iinv;ncut:1.1 iinv;nit:2.2 iinv;tol:3.3'
Cr  *'element-name' is the entry name of the structure
Cr  *';alias'  is a string naming a alias for this part of element-name
Cr   (used in input/output)
Cr  *',mask' is a mask (for single integer entries only) that singles
Cr   out individual bits within the integer
Cr  *':range[.upper-bound]' points to a starting offset for names
Cr   that hold vector entries, and optionally a vector upper bound
C ----------------------------------------------------------------------
C     implicit none
      integer mask(1),range(2,1)
      character*(*) sname,trunam,alias,switch
      integer j1,j2,is1,is2,ia,is,it,nw,i,j,ix(2),a2vec
      logical a2bin
      character ct*1,strn*40

      call words(sname,nw)
      call word(sname,1,j1,j2)
      trunam = sname(j1:j2)
      it = j2+1
      ia = 0
      is = 0

      do  12  i = 1, nw-1

        call ivset(range(1,i),1,2,-1)
        j1 = j2+1
        mask(i) = 0
        call nword(sname,1,j1,j2)
        is1 = j1
        call nwordg(sname,0,';:?, ',1,is1,is2)
        trunam(it+1:) = sname(is1:is2)
        it = it + is2-is1+2
        alias(ia+1:) = sname(is1:is2)
        switch(is+1:) = 'none'
*       print *, 'alias is',alias
        if (is2 .ne. j2) then
   10     call nwordg(sname,0,';:?, ',2,is1,is2)
          ct = sname(is1-1:is1-1)
          if (ct .eq. ';') then
            alias(ia+1:) = sname(is1:is2)
          elseif (ct .eq. ':') then
            strn = sname(is1:is2)
            j = 0
            j = a2vec(strn,len(strn),j,2,'.;:?, ',6,2,2,ix,range(1,i))
            if (j.lt.0) call rxs('salias: failed to parse',sname(1:is2))
          elseif (ct .eq. '?') then
            switch(is+1:) = sname(is1:is2)
          elseif (ct .eq. ',') then
            j = is1-1
            if (.not. a2bin(sname,mask(i),2,0,' ',j,is2-1)) call
     .        rxs('salias: failed to parse mask in ',sname)
          endif
          if (is2 .lt. j2) goto 10
        endif

        call skp2bl(alias,len(alias),ia)
        ia = ia+1
        call skp2bl(switch,len(switch),is)
        is = is+1

   12 continue


C      print *, 'sname ',sname
C      print *, 'trunam ',trunam
C      print *, 'lres=',lres,' mask=',(mask(i),i=1,nw-1)
C      print *, 'alias ',alias
C      print *, 'switch ',switch
C      stop

Cgetarg       end
       end subroutine salias 


