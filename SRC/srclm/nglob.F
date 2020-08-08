      integer function nglob(name)
C- Returns value of integer global variable
C ----------------------------------------------------------------------
Ci Inputs
Ci   name  :string specifying variable name
Co Outputs
Co  nglob  :nearest integer of variable's value
Co         :nglob aborts if name is missing from the table.
Cr Remarks
Cr  Variables may be set with dglob, below.
C ----------------------------------------------------------------------
C     implicit none
      character *(*) name
C ... Local parameters
      integer nigvmx,ngv,iset,it,i
      parameter (nigvmx=30)
      double precision gv(nigvmx),dglob,val
      character*10 sgv(nigvmx)
      save gv,sgv,ngv
      data ngv/0/

      do  10  i = 1, ngv
        it = i
        if (name .eq. sgv(it)) goto 11
   10 continue
      call rxs2('nglob: no variable "',name,'" in table')
   11 continue
      nglob = nint(gv(it))
      return

      entry dglob(name,val,iset)
C- Set global variable and value, or retrieve double precision value
C ----------------------------------------------------------------------
Ci Inputs
Ci   name  :string specifying variable name
Ci     val :variable's value, used when iset=1
Ci    iset :0 return value of variable with name 'name'
Ci          dglob aborts if name is missing.
Ci         :1 set variable's value to val.  If variable is missing from
Ci          the table, add it to the table, and assign its value
Co Outputs
Co  dglob  :(iset=0) the index to the table entry of named variable
Co          (iset=1) value of variable 'name'
Cr Remarks
C ----------------------------------------------------------------------

      if (iset .eq. 1) then
        do  20  i = 1, ngv
          it = i
          if (name .eq. sgv(it)) goto 21
   20   continue
        it = -1
   21   continue
        if (it .eq. -1) then
          ngv = ngv+1
          if (ngv .gt. nigvmx) call rx('dglob increase nigvmx')
          sgv(ngv)= name
          it = ngv
        endif
        gv(it) = val
        dglob = it

      else
        do  30  i = 1, ngv
          it = i
          if (name .eq. sgv(it)) goto 31
   30   continue
        call rxs2('dglob: no variable "',name,'" in table')
   31   continue
        dglob = gv(it)
      endif

      end


C      subroutine fmain
C
C      integer nglob,i
C      double precision dglob,xx
C
C      xx = dglob('nspin',2.2d0,1)
C      print *, 'should be 1',xx
C      xx = dglob('nsp',3.3d0,1)
C      print *, 'should be 2',xx
C      xx = dglob('nbas',4.4d0,1)
C      print *, 'should be 3',xx
C      xx = dglob('nsp',5.45d0,1)
C      print *, 'should be 2',xx
C      xx = dglob('nsp',0d0,0)
C      print *, 'should be 5.45',xx
C      xx = nglob('nsp')
C      print *, 'should be 5',xx
C      xx = nglob('nspin')
C      print *, 'should be 2',xx
C      xx = nglob('nbas')
C      print *, 'should be 4',xx
C      xx = dglob('nbas',0d0,0)
C      print *, 'should be 4.4',xx
C      end

