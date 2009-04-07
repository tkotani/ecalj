      integer function iclbas(class,iclass)
C- Returns an index to iclbas atom in basis given class
C ----------------------------------------------------------------
Ci Inputs
Ci   class,iclass
Co Outputs
Co   iclbas
Cb Bugs
Cr    No check is made on the length of class
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer class,iclass(1)
C Local parameters
      integer ibas

      do  10  ibas = 1, 1000000
        iclbas = ibas
        if (iclass(ibas) .eq. class) return
   10 continue
      call rx('ICLBAS: bad input')
      end
      integer function iclbsj(ic,ipc,nbas,nrbas)
C- Returns an index to nrbas atom in basis given the class
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ic    :class index
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nbas  : abs  = number of atoms in the basis
Ci         : sign <0 to return with -n if there are fewer than nrbas
Ci         :  members of class ic, where n=number members of class ic
Ci   nrbas :the nrbas-th basis atom of class ic is sought
Co Outputs:
Co   iclbsj:the nrbas-th atom belonging to class ic
Cu Updates
Cu   10 May 01 Returns -n rather than -1 when nrbas-th atom not found
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ic,nbas,ipc(nbas),nrbas
C Local parameters
      integer ibas,n,nbasa

      iclbsj = 1
      n = 0
      nbasa = iabs(nbas)
      do  10  ibas = 1, nbasa
        if (ipc(ibas) .eq. ic) n = n+1
        if (n .eq. nrbas) then
          iclbsj = ibas
          return
        endif
   10 continue

      if (nbas .lt. 0) then
        iclbsj = -n
        return
      endif

      call fexit3(-1,111,' Exit -1 ICLBSJ: sought atom no. %i'//
     .  ' in class %i but only %i atoms exist',nrbas,ic,n)

      end

