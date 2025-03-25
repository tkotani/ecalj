module m_nvfortran !this is because nvfortran24.1 have no implementation to findloc,count
  public findloc,count
  interface findloc 
     module procedure findlocl,findloci
  endinterface findloc
  interface count 
     module procedure countl
  endinterface count
  private
contains
  pure integer function findlocl(lll,value,dim,back) !this is for
    implicit none
    intent(in):: lll,value,dim,back
    integer:: nnn,ns,ne,inc,i
    logical:: lll(:),value
    logical,optional:: back
    integer,optional:: dim
    nnn=size(lll)
    if(present(back)) then
       if(back) then
          ns=nnn
          ne=1
          inc=-1
       endif
    else
       ns=1
       ne=nnn
       inc=1
    endif     !   nn1 = findloc([(rsmh1(i,j)>0d0,i=1,nnx)],value=.true.,dim=1,back=.true.)! 1st MTO sets
    findlocl=0
    do i=ns,ne,inc
       if(intl(value)-intl(lll(i))==0) then
          findlocl=i
          return
       endif
    enddo
  contains
    pure integer function intl(lll)
      logical,intent(in):: lll
      if(lll) intl=1
      if(.not.lll)intl=0
    endfunction intl
  end function findlocl
  pure integer function findloci(lll,value,dim,back) 
    implicit none
    intent(in):: lll,value,dim,back
    integer:: nnn,ns,ne,inc,i, lll(:),value
    logical,optional:: back
    integer,optional:: dim
    nnn=size(lll)
    if(present(back)) then
       if(back) then
          ns=nnn
          ne=1
          inc=-1
       endif
    else
       ns=1
       ne=nnn
       inc=1
    endif
    findloci=0
    do i=ns,ne,inc
       if(lll(i)==value) then
          findloci=i
          return
       endif
    enddo
  end function findloci
  integer function countl(lll)
    integer:: i,m
    logical:: lll(:)
    m=0
    do i=1,size(lll)
       if(lll(i)) m=m+1
    enddo
    countl=m
  end function countl
endmodule m_nvfortran
