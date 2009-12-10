      subroutine getujh(sspec,nl,nbas,ipc,dclabl,idu,uh,jh)
C- Unpack Hubbard U's and J's from spec structure
C ----------------------------------------------------------------------
Ci Inputs:
Ci   sspec,nl,nbas,ipc
Co Outputs:
Co   idu,uh,jh
Cr Remarks
Cr   self explanatory
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nl,nbas,ipc(1),idu(4,nbas)
      double precision sspec(1),dclabl(1),uh(4,nbas),jh(4,nbas)
C Local Variables
      integer ib,is,ipr,iprint,i1mach
      character*8 clabl

      ipr = iprint()
      if (ipr .gt. 30) then
        call awrit0(' GETUJH: unpacking U and J ..',' ',120,i1mach(2))
      endif
      do  ib = 1, nbas
        is = ipc(ib)
        call upack('spec idu uh jh',sspec,is,
     .             idu(1,ib),uh(1,ib),jh(1,ib),0)
        if (ipr .gt. 30) then
          call r8tos8(dclabl(is),clabl)
          call awrit7('Atom %#5i '//clabl//' IDU=%n:2i U=%n:2d J=%n:2d',
     .      ' ',256,i1mach(2),ib,nl,idu(1,ib),nl,uh(1,ib),nl,jh(1,ib))
        endif
      enddo
      end
