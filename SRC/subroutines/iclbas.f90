integer function iclbsj(ic,ipc,nbas,nrbas)  !- Returns an index to nrbas atom in basis given the class
  !i   ic    :class index
  !i   ipc   :class index: site ib belongs to class ipc(ib) 
  !i   nbas  : abs  = number of atoms in the basis
  !i         : sign <0 to return with -n if there are fewer than nrbas
  !i         :  members of class ic, where n=number members of class ic
  !i   nrbas :the nrbas-th basis atom of class ic is sought
  !o Outputs:
  !o   iclbsj:the nrbas-th atom belonging to class ic
  !u Updates
  !u   10 May 01 Returns -n rather than -1 when nrbas-th atom not found
  ! ----------------------------------------------------------------------
  implicit none
  integer :: ic,nbas,ipc(*),nrbas
  integer :: ibas,n,nbasa
  iclbsj = 1
  n = 0
  nbasa = iabs(nbas)
  do  10  ibas = 1, nbasa
     if (ipc(ibas) == ic) n = n+1
     if (n == nrbas) then
        iclbsj = ibas
        return
     endif
10 enddo
  if (nbas < 0) then
     iclbsj = -n
     return
  endif
  call rx3('ICLBSJ: sought atom no.#1 in class #2 but only #3 atoms exist. #1#2#3=',nrbas,ic,n)
end function iclbsj
