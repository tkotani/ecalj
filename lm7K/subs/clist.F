Cdelw1 var ['oilst', '0', 'oilst', 'nlist', 'integer', 'iv_w_', 'clist']
Cdelw1 undel [['wref'], 'oilst', 'clist']
      subroutine clist(lstyle,slist,dclabl,z,nclass,nlist,list)
C- Generates a list of classes from a string specification
C ----------------------------------------------------------------
Ci Inputs
Ci   slist:  string specifying list of classes
Ci   lstyle: style of slist specification; see Remarks
Ci   nclass  number of classes.
Co Outputs
Co   nlist,list a list of classes satisfying specifications
Cr Remarks
Cr *Syntax of slist: depends on one of three styles (lstyle)
Cr
Cr *lstyle=1 : a list of integers; see mkilst.f for complete syntax.
Cr             Example: '1,4:6,11' generates a list of five numbers,
Cr             1,4,5,6,11.
Cr
Cr *lstyle=2 : the list is specified according to an expression.
Cr             The expression can involve the class index ic and 
Cr             atomic number z.  Any class satisfying expression is
Cr             included in the list.  Example:  'ic<6&z==14'
Cr
Cr *lstyle=3 : is specifically for unix systems.  slist is a filename
Cr             with the usual unix wildcards, eg a[1-6].  'clist'
Cr             makes a system 'ls' call for that string; any class
Cr             which 'ls' finds is included in the list.
C ----------------------------------------------------------------
C     implicit none
      integer lstyle,nlist,nclass,list(1)
      character*(*) slist
      double precision z(nclass),dclabl(nclass)
C Local variables
Cdelw1 do not change  [['wref'], 'oilst', 'clist']
      integer iv0,ic,ival,is,i,j,ls,ls1,oilst
      logical a2bin,sw
      character strn*120,filnam*72,cnam*72,clabl*8
C Heap:
      integer w(1)
      common /w/ w

      ls = len(slist)
      nlist = 0
      goto (10,20,30) lstyle
      call rxi('clist: bad style',lstyle)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oilst]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C -- lstyle=1 ---
   10 continue
      call mkils0(slist,nlist,i)
Cdelw1 do not change oilst because of ['wref']
      call defi(oilst, nlist)
Cdelw1 do not change oilst because of ['wref']
      call mkilst(slist,nlist,w(oilst))
Cki#error, have return with len(w_varlist)>0 at line 50
Cdelw1 w_varlist remains: oilst 
Cdelw1 w_varlistundel: [oilst]
Cdelw1 w_varlist (undel), remains: [oilst]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate oilst because of [wref]
      if (nlist .eq. 0) return
Cdelw1 do not change oilst because of ['wref']
      call ishell(nlist,w(oilst))
Cdelw1 do not change oilst because of ['wref']
      list(1) = w(oilst)
      j = 1
      do  12  i = 2, nlist
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 do not change oilst because of ['wref']
        if (w(oilst+i-1) .gt. list(j)
     .      .and. w(oilst+i-1) .le. nclass) then
          j = j+1
Cdelw1 1 is digit, dropped
Cdelw1 do not change oilst because of ['wref']
          list(j) = w(oilst+i-1)
        endif
   12 continue
      nlist = j
Cdelw1 rlse name= oilst old_list= oilst 
Cdelw1 rlse name= oilst new_list= (None)
Cdelw1 not deallocate oilst because of [wref]
Cdelw1       call rlse(oilst)
       call rlse(oilst)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oilst]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C --- lstyle=2 ---
   20 continue
      call numsyv(iv0)
      nlist = 0
      do  42  ic = 1, nclass
        call lodsyv('ic',1,dble(ic),ival)
        call lodsyv('z',1,z(ic),ival)
        is = 0
        if (a2bin(slist,sw,0,0,slist(ls:ls),is,ls)) then
          if (sw) then
            nlist = nlist+1
            list(nlist) = ic
          endif
C   ... Abort if a2bin can't parse expression
        else
          call rxs('clist: failed to parse',slist)
        endif
   42 continue
      call clrsyv(iv0)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oilst]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C --- lstyle=3 ---
   30 continue
      nlist = 0
      call skpblb(slist,ls,ls1)
      call ffnam(slist(1:ls1+1),filnam)
      do  44  ic = 1, nclass
        call pshpr(0)
        call r8tos8(dclabl(ic),clabl)
        call ffnam(clabl,cnam)
        call poppr
        call awrit0('%xls ' // filnam //'%a|grep -s '
     .    // cnam // '%a>/dev/null',strn,len(strn),0)
        call locase(strn)
        call fsystm(strn,j)
        if (j .eq. 0) then
          nlist = nlist+1
          list(nlist) = ic
        endif
   44 continue

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oilst]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

