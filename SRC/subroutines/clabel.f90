subroutine clabel(slabl,is,idx,clabl)
  !- Make class label from species label
  !i Inputs
  !i   slabl,is: species label is slabl(is)
  !i   idx:      class to make is the idx_th class for species is
  !i             From the ics table, the class corrsponding to idx is
  !i             iclbsj(is,ics,-nclass,idx)
  !o   clabl:    class label
  !     implicit none
  integer :: is,idx
  character(8) :: slabl(is),clabl
  clabl = slabl(is)
  if (idx > 1) write(clabl,"(i5)") idx
end subroutine clabel

