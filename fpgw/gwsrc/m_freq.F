!>Frequency mesh generator
!! - OUTPUT
!!   - fhris :histgram bins to accumlate im part
!!   - freq_r: omega along real axis
!!   - freq_i: omega along imag axis
!!   - wiw: integration weight along im axis
!!   - npm: npm=1 means only positive omega;npm=2 means positive and negative omega.
!! - NOTE: change of frequency mesh defined here may destroy consistency or not. Need check
      module m_freq
      real(8),allocatable,protected:: frhis(:),freq_r(:),freq_i(:),wiw(:)
      integer,protected:: nwhis,npm,nw_i,nw 

c      real(8),allocatable,protected:: frhis0(:)
c      integer,protected:: nwhis0
      contains
!> Get data set for m_freq. All arguments are input.
!! - This read GWinput (dw,omg_c) and TimeReversal() 
!! - All arguments are input
      subroutine getfreq(epsmode,realomega,imagomega,tetra,omg2max,wemax,niw,ua,mpi__root)
      use m_keyvalue,only:getkeyvalue

      implicit none
      integer,intent(in):: niw !,nw_input
      logical,intent(in):: realomega,imagomega,tetra,mpi__root,epsmode
      real(8),intent(in):: omg2max,ua

      real(8),allocatable:: freqx(:),wx(:),expa(:)
      logical:: timereversal,onceww
      integer:: nw2,iw,ihis
      real(8)::omg_c,dw,omg2,wemax
      real(8), allocatable :: freqr2(:)  ,frhis_tmp(:)
      real(8)::  pi = 4d0*datan(1d0), aa,bb,ratio,oratio,daa
      integer::nee,noo,ifif,ifile_handle

      logical,save:: done=.false.
      if(done) call rx('gerfreq is already done') !sanity check
      done =.true.
      nw=-99999 !for sanity check

c      nw = nw_input
!! Histogram bin divisions
!! We first accumulate Imaginary parts. 
!! Then it is K-K transformed to obtain real part.

c      call getkeyvalue("GWinput","dw",dw )
c      call getkeyvalue("GWinput","omg_c",omg_c )
c      write(6,"('dw, omg_c= ',2f13.5)") dw, omg_c
      call getkeyvalue("GWinput","HistBin_ratio",oratio, default=1.03d0) 
      call getkeyvalue("GWinput","HistBin_dw",dw, default=1d-5) !a.u.
      aa = oratio-1d0
      bb = dw/aa
      iw = 0d0
      do 
        iw=iw+1
        if( bb*( exp(aa*(iw-1)) - 1d0 ) >omg2max+1d-6) exit
      enddo  
      nwhis = iw+2 !+2 for margin. Necessary?
      allocate(frhis(1:nwhis+1))
      do iw = 1,nwhis+1
        frhis(iw) = bb*( exp(aa*(iw-1)) - 1d0 )
      enddo
      write(6,"('dw, omg_ratio, nwhis= ',d9.2,f13.5,i6)") dw, aa,nwhis

!! Determine nw. Is this correct?
      do iw=3,nwhis
         omg2 = (frhis(iw-2)+frhis(iw-1))/2d0
         if (omg2 > wemax/2d0 ) then !>dw*(nw_input-3)) then !omg is in unit of Hartree
            nw=iw
            exit 
         endif
      enddo
!! document need to be fixed...
c      nw=nw2-1      ! nw+1 is how many points of real omega we use
                     ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
                     ! maximum nw=nw2-1 because nwhis=nw2-1
!! document need to be fixed...
                     !nw is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
                     !Here frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
                     !nw was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
                     ! for correlation Coulomb Wc(omg),
                     ! and one more point omg=dw*(nw-1) needed for extrapolation.
                     ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
                     ! and one more point omg=frhis_m(nw) needed for extropolation
                     ! used in subroutine alagr3z in  sxcf.f.

!! Determine freq_r
      if(epsmode) then
        nw  = nwhis-1
      endif
      allocate(freq_r(0:nw))
      freq_r(0)=0d0
      do iw=1,nw
        freq_r(iw)=(frhis(iw)+frhis(iw+1))/2d0
      enddo

!! Timereversal=F is implimented only for tetra=T and sergeyv=T
!! nw_i and npm
      npm=1
      nw_i=0
      if(.not.timereversal()) then
         write(6,"('TimeReversal off mode')")
         npm=2
         nw_i=-nw
         if(.not.tetra)   call rx( ' tetra=T for timereversal=off')
      endif
      write(6,*)'Timereversal=',Timereversal()

!! Write freq_r
      if(realomega .and. mpi__root) then
        ifif=ifile_handle()
         open(UNIT=ifif,file='freq_r') !write number of frequency points nwp and frequensies in 'freq_r' file
         write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
         do iw= nw_i,-1
            write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw 
         enddo
         do iw= 0,nw
            write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw 
         enddo
         close(ifif)
      endif

!! Determine freq_i  : gaussian frequencies x between (0,1) and w=(1-x)/x      
      if (imagomega) then
         write(6,*)' freqimg: niw =',niw
         allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )
         call freq01 (niw,ua, 
     o        freqx,freq_i,wx,expa) 
         allocate(wiw(niw))
         do iw=1,niw
           wiw(iw)=wx(iw)/(2d0*pi*freqx(iw)*freqx(iw))
         enddo
         deallocate(freqx,wx,expa)
      endif

!! Plot frhis 
      if(onceww(1)) then
        write(6,*)' we set frhis nwhis noo-->nee=',nwhis,noo,nee
        write(6,*)' --- Frequency bins to accumulate Im part  (a.u.) are ---- '
        do ihis= 1, nwhis !min(10,nwhis)
          write(6,"(' ihis Init  End=', i5,2f18.11)") ihis,frhis(ihis),frhis(ihis+1)
        enddo
      endif
      end subroutine getfreq
      end module m_freq
