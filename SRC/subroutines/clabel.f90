      subroutine clabel(slabl,is,idx,clabl)
C- Make class label from species label
Ci Inputs
Ci   slabl,is: species label is slabl(is)
Ci   idx:      class to make is the idx_th class for species is
Ci             From the ics table, the class corrsponding to idx is
Ci             iclbsj(is,ics,-nclass,idx)
Co   clabl:    class label
C     implicit none
      integer is,idx
      character*8 slabl(is),clabl
      clabl = slabl(is)
      if (idx .gt. 1) write(clabl,"(i5)") idx 
      end

