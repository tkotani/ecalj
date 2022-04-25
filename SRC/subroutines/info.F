      subroutine info(jpr,l1,l2,string,a1,a2)
C- Printout when ipr>=jpr, with some blank lines
C ----------------------------------------------------------------------
Ci Inputs
Ci   jpr   :print with verbosity is at least abs(jpr)
Ci         :jpr<0 -> print only if master node (MPI)
Ci   l1    :number of blank lines to print out before string
Ci   l2    :number of blank lines to print out after string
Ci         :l2<0 => write without carriage return
Ci   string:string to print
Co Outputs
Co   string is printed on stdo 
Cr Remarks
Cu Updates
Cu   24 Nov 02 New info2, info5, info8
Cu   26 Jun 02 Add arguments a1,a2 so info can write data
Cu   02 Nov 01 Use awrite instead of write statement
Cu             Special handling of l2=-1
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer jpr,l1,l2,a1,a2,a3,a4,a5,a6,a7,a8
      character*(*) string
C ... Local parameters
      integer i,iprint,stdo,lgunit,recl,awrite,mpipid
      parameter (recl=1024)
      character*(recl) lstr
      entry info0(jpr,l1,l2,string)
      entry info2(jpr,l1,l2,string,a1,a2)
      entry info5(jpr,l1,l2,string,a1,a2,a3,a4,a5)
      entry info8(jpr,l1,l2,string,a1,a2,a3,a4,a5,a6,a7,a8)
      if (iprint() .lt. iabs(jpr)) return
      if (jpr .lt. 0) then
        i = mpipid(1)
        if (i .ne. 0) return
      endif
      stdo = lgunit(1)
      do  10  i = 1, l1
         write(stdo,100)
 10   continue
      if (l2 .ge. 0) then
        call awrit8(string,' ',recl,stdo,a1,a2,a3,a4,a5,a6,a7,a8)
C       call awrit2(string,' ',recl,stdo,a1,a2)
      else
        i = awrite(string,lstr,recl,0,a1,a2,a3,a4,a5,a6,a7,a8)
        call cwrite(lstr,0,i-1,0)
      endif
      do  20  i = 1, l2
         write(stdo,100)
 20   continue
  100 format(1x)
      return
      end

