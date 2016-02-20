!>Frequency mesh generator
!! - OUTPUT
!!   - fhris :histgram bins to accumlate im part
!!   - freq_r: omega along real axis
!!   - freq_i: omega along imag axis
!!   - wiw: integration weight along im axis
!!   - npm: npm=1 means only positive omega;npm=2 means positive and negative omega.
!! - NOTE: change of frequency mesh defined here may destroy consistency or not. Need check
      module m_freq
      real(8),allocatable:: frhis(:),freq_r(:),freq_i(:),wiw(:)
      integer:: nwhis,npm,nw_i,nw 

      contains
!> Get data set for m_freq. All arguments are input.
!! - This read GWinput (dw,omg_c) and TimeReversal() 
!! - All arguments are input
      subroutine getfreq(epsmode,realomega,imagomega,tetra,omg2max,wemax,niw,ua,mpi__root)
      use keyvalue,only: getkeyvalue
      implicit none
      integer,intent(in):: niw !,nw_input
      logical,intent(in):: realomega,imagomega,tetra,mpi__root,epsmode
      real(8),intent(in):: omg2max,ua

      real(8),allocatable:: freqx(:),wx(:),expa(:)
      logical:: timereversal,onceww
      integer:: nw2,iw,ihis
      real(8)::omg_c,dw,omg2,wemax
      real(8), allocatable :: freqr2(:)  
      real(8)::  pi       = 4d0*datan(1d0)

      logical,save:: done=.false.
      if(done) call rx('gerfreq is already done') !sanity check
      done =.true.
      nw=-99999 !for sanity check

c      nw = nw_input
!! We first accumulate Imaginary parts. 
!! Then it is K-K transformed to obtain real part.
      call getkeyvalue("GWinput","dw",dw )
      call getkeyvalue("GWinput","omg_c",omg_c )
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

!! frhis_m
      nw=nw2-1       ! nw+1 is how many points of real omega we use
                     ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
                     ! maximum nw=nw2-1 because nwhis=nw2-1
      do iw=3,nw2-1
                     !nw is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
                     !Here frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
                     !nw was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
                     ! for correlation Coulomb Wc(omg),
                     ! and one more point omg=dw*(nw-1) needed for extrapolation.
                     ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
                     ! and one more point omg=frhis_m(nw) needed for extropolation
                     ! used in subroutine alagr3z in  sxcf.f.
         omg2 = (freqr2(iw-2)+freqr2(iw-1))/2d0
         if (omg2 > wemax/2d0 ) then !>dw*(nw_input-3)) then !omg is in unit of Hartree
            nw=iw
            exit     ! 'nw_input' is only used to get maximum frequency for
                     !  dressed coulomb line
         endif
      enddo

      if(epsmode) then
        nw  = nwhis-1
c        niw = 0
      endif

      allocate(freq_r(0:nw))
      freq_r(0)=0d0
      do iw=1,nw
        freq_r(iw)=(frhis(iw)+frhis(iw+1))/2d0
      enddo

!! Plot frhis ----------
      if(onceww(1)) then
        write(6,*)' --- Frequency bins to accumulate Im part  (a.u.) are ---- '
        do ihis= 1, min(10,nwhis)
          write(6,"(' ihis Init  End=', i5,2f13.6)") ihis,frhis(ihis),frhis(ihis+1)
        enddo
        write(6,*) 'ihis  ...'
        do ihis= max(min(10,nwhis),nwhis-10), nwhis
          write(6,"(' ihis Init  End=', i5,2f13.6)") ihis,frhis(ihis),frhis(ihis+1)
        enddo
      endif
!! Timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
      if(.not.timereversal()) then
         write(6,"('TimeReversal off mode')")
         npm=2
         nw_i=-nw
         if(.not.tetra)      call rx( ' tetra=T for timereversal=off')
      endif
      write(6,*)'Timereversal=',Timereversal()
      if(realomega .and. mpi__root) then
         open(UNIT=3111,file='freq_r') !write number of frequency
                                !points nwp and frequensies in 'freq_r' file
         write(3111,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
         do iw= nw_i,-1
            write(3111,"(d23.15,2x,i6)") -freq_r(-iw),iw !This file is reffere by hsfp0 and so.
         enddo
         do iw= 0,nw
            write(3111,"(d23.15,2x,i6)") freq_r(iw),iw !This file is reffere by hsfp0 and so.
         enddo
         close(3111)
      endif
!! set freq_i
      if (imagomega) then
         write(6,*)' freqimg: niw =',niw
         allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )
         call freq01 (niw,ua, 
     o        freqx,freq_i,wx,expa)
         ! Generate gaussian frequencies x between (0,1) and w=(1-x)/x      
         allocate(wiw(niw))
         do iw=1,niw
           wiw(iw)=wx(iw)/(2d0*pi*freqx(iw)*freqx(iw))
         enddo
         deallocate(freqx,wx,expa)
      endif
      end subroutine getfreq
      end module m_freq
