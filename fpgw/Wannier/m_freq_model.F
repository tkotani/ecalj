!>Frequency mesh generator
!! - OUTPUT
!!   - fhris :histgram bins to accumlate im part
!!   - freq_r: omega along real axis
!!   - freq_i: omega along imag axis
!!   - wiw: integration weight along im axis
!!   - npm: npm=1 means only positive omega;npm=2 means positive and negative omega.
!! - NOTE: change of frequency mesh defined here may destroy consistency or not. Need check
      module m_freq_model
      real(8),allocatable:: frhis(:)
      integer:: nwhis,npm,nw_i,nw 

      contains
!> Get data set for m_freq_model. All arguments are input.
!! - This read "Input.dat" (dw,omg_c) and TimeReversal() 
!! - All arguments are input
      subroutine getfreq_model(omg2max)
      use m_keyvalue,only:getkeyvalue

      implicit none
      logical:: tetra
      real(8),intent(in):: omg2max

      real(8),allocatable:: freqx(:),wx(:),expa(:)
      logical:: timereversal,onceww
      integer:: nw2,iw,ihis
      real(8)::omg_c,dw,omg2,wemax
      real(8), allocatable :: freqr2(:)  
      real(8)::  pi       = 4d0*datan(1d0)

      logical,save:: done=.false.


      tetra = .True.

      if(done) call rx('gerfreq is already done') !sanity check
      done =.true.
      nw=-99999 !for sanity check

c      nw = nw_input
!! We first accumulate Imaginary parts. 
!! Then it is K-K transformed to obtain real part.
      call getkeyvalue("Input.dat","dw",dw )
      call getkeyvalue("Input.dat","omg_c",omg_c )
      write(6,"('dw, omg_c= ',2f13.5)") dw, omg_c
!! histogram bin divisions
      nw2=int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
      allocate(freqr2(nw2))     !+1 b/c (iw-1)
      do iw=1,nw2
         freqr2(iw)=dw*(iw-1)+dw**2/2./omg_c*(iw-1)**2
      enddo                     !linear + quadratic term
      if (nw2 < 2 ) call rx( "m_freq: nw2 < 2")
      if (dw*(nw-2) > freqr2(nw2-1)) call rx("m_freq: dw*(nw-2) > freqr2(nw2-1)")
      nwhis = nw2-1
      allocate(frhis(1:nwhis+1))
      frhis = freqr2(1:nwhis+1)
      write(6,*)' we set frhis nwhis=',nwhis
c      nw  = idint (wemax/2d0/dw) + 3

!! Timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
!      if(.not.timereversal()) then
!         write(6,"('TimeReversal off mode')")
!         npm=2
!         nw_i=-nw
!         if(.not.tetra)      call rx( ' tetra=T for timereversal=off')
!      endif
!      write(6,*)'Timereversal=',Timereversal()

      end subroutine getfreq_model
      end module m_freq_model
