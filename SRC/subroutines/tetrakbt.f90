logical function usetetrakbt()
  usetetrakbt=.false.
end function usetetrakbt

!! ---------------------
subroutine tetrakbt(ebfin,eafin,v,voltet, frhis, nwhis, kbt, &
     wtthis ) !this is accumlating variable
  !! kbt vertion
  !! Histgram weights for each bin.
  !! The i-th bin of the Histgrams is [frhis(i) frhis(i+1)].
  !! wtthis(ihis) += \int_{frhis(ihis)^{frhis(ihis+1)} d \omega  \int d^3k (fb-fa) * \delta(\omega +v(k) )
  !!  Total sum of the histgram weight is equal to voltet.
  !!  wtthis(ihis): Note this is accumulating variable.
  implicit none
  real(8),intent(in):: ebfin(4),eafin(4),v(4), voltet,kbt
  real(8):: ebft(4),eaft(4),ebf(3),eaf(3)
  integer(4)::  inihis, iedhis, nwhis,ichk=2313
  integer(4) :: ieaord(1:4),i,isig,idif(3),idf,ix,n,itmp,ihis
  real(8)   ::  WW(4),  norm,s1,s2,pvn, &
       integb3p, integb3m , integb2p, integb2m,ww2p,ww3p
  real(8)::frhis(nwhis+1),wtthis(nwhis,4),intega,integb,stot,xxx,wx,fafac,fbfac
  logical ::chkwrt=.false., matrix_linear
  real(8):: www(1:4)=.25d0,ec,wcg(4),wttt,el,eh !kkvkin(4,4),
  !      if(chkwrt) then
  !        print *, ' tetrakbt: ==========================='
  !        write(6,"(' i=',i3,' x ea eb =',3f10.5)")
  !     &    (i, v(i), eafin(i), ebfin(i),i=4,1,-1)
  !      endif

  !! sorted so that ww(1)<ww(2)<ww(3)<ww(4)
  call sortea( -v,ieaord,4,isig)
  WW(1:4) = -v( ieaord(1:4) )
  ebft = ebfin( ieaord(1:4) )
  eaft = eafin( ieaord(1:4) )

  if(chkwrt) then
     write(6,"(a,f10.5,i5)")' tetrakbt: ===========================',voltet,nwhis
     write(6,"(' i=',i3,' ww ebft eaft =',3f10.5)") (i, ww(i), ebft(i), eaft(i),i=4,1,-1)
  endif
  if(( .NOT. (WW(1)<=WW(2))) .OR. ( .NOT. (WW(2)<=WW(3))) .OR. ( .NOT. (WW(3)<=WW(4))) ) then
     write(6,"(/,' --- tetrakbt: wrong order WW=',4d14.6)") WW
     call rx( 'tetrakbt: wrong order of WW')
  endif

  inihis= -999
  iedhis= -999
  ix=1
  do ihis = 1,nwhis
     if(ix==1 .AND. WW(ix)<frhis(ihis+1)) then
        inihis = ihis
        ix=4
     endif
     if(ix==4 .AND. WW(ix)<frhis(ihis+1)) then
        iedhis = ihis
        exit
     endif
  enddo
  if(iedhis==-999 .OR. inihis==-999) then
     call rx( ' tetrakbt: can not find inihis iedhis')
  endif
  !      if(chkwrt)
  !      print *,' inihis iedhis=', inihis,iedhis
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      print *,' ###########integtetn test #########'
  !      ww(1:4)=(/1d0,1.5d0,9.5d0,10d0/)
  !      do ix=1,101
  !        wx = (ix-1)/10d0
  !        call integtetn(WW, Wx, xxx)
  !        write(61,"(i3,2d23.16)")ix,wx,xxx
  !      enddo
  !      stop 'test end integtetn---'
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  call integtetn(WW, WW(4), norm)
  !      if(chkwrt) write(6,"(' norm voltet=',2d24.16)") norm,voltet
  pvn = voltet/norm
  !      write(6,"(' norm voltet pvn=',3d14.6)") norm,voltet,pvn
  intega = 0d0
  do ihis = inihis, iedhis
     if( frhis(ihis+1)>ww(4) ) then
        integb = norm
     else
        call integtetn(WW, frhis(ihis+1), integb)
     endif

     el = frhis(ihis)
     eh = frhis(ihis+1)
     ! Center energy for histgram bin= [frhis(ihis), frhis(ihis+1)]
     call eaf_triangle(el,eh,ww, eaft, fafac)  ! ww should satisfy ww(1)<=ww(2)<=ww(3)<=ww(4)
     call eaf_triangle(el,eh,ww, ebft, fbfac)
     !        write(6,"('fbfac-fafac=',2f9.4,2x,f9.4)") fbfac, fafac, fbfac-fafac
     wtthis(ihis,1) = wtthis(ihis,1) &
          + (fafac-fbfac) * pvn*(integb - intega) !(fafac-fbfac) is introduced.

     if(chkwrt) then
        write(6,"(' ihis [Init End] wtt ',i5,2f9.4,d13.4 &
             ,' fbfa=', 3f9.4)") &
             ihis, frhis(ihis), frhis(ihis+1),pvn*(integb-intega), fbfac,fafac, fafac-fbfac
     endif

     intega = integb
  enddo
  !      if(chkwrt) then
  !        dini = 0d0
  !        do ihis=inihis,iedhis
  !          if(ihis>1) dini=wtthis(ihis-1)
  !          write(6,"(' ihis [Init  End]', i4,2f13.6,' wtt=',f13.6)")
  !     &    ihis,dini,frhis(ihis),wtthis(ihis)
  !          dini = frhis(ihis)
  !        enddo
  !      endif
  !      if(chkwrt) write(ichk,*) ' end of intttvc6'
end subroutine tetrakbt

!-----------------------------------------------------
subroutine eaf_triangle(el,eh,ww,eaft, fafac)
  !! We assume ww are soreted as ww(1)<ww(2)<ww(3)<ww(4)
  implicit none
  real(8),intent(in) :: el,eh, ww(4),eaft(4)
  real(8),intent(out):: fafac
  real(8):: eafww,kbt,eaf(3),factri,ecenter,fafac2,fafac1,s2,s3,s2e,s3e,eps=1d-10
  kbt=0d0 !dummy now
  ecenter=(el+eh)/2d0

  !! For cases when tetrahderon is neighter touching eh nor el.
  if( ww(1) <= ecenter .AND. ecenter <=ww(4) ) then
     continue
  elseif( ecenter <  ww(1) .AND.  eh< ww(4)) then
     ecenter = eh
  elseif( ecenter >  ww(4) .AND.  el> ww(1) ) then
     ecenter = el
  else
     ecenter= (ww(2)+ww(3)) /2d0 ! When tetrahderon is inside between [el, ecenter] or [ecenter,eh]
  endif

  if(     ecenter <  ww(1)) then
     write(6,*) 'el eh',el,eh
     write(6,*) 'eaf_triangle',ecenter,ww(1)
     call rx('eaf_triangle: ecenter < ww(1)')
  elseif( ecenter <= ww(2)) then
     eaf(1) = eafww(ecenter,1,2,ww,eaft) !line 12
     eaf(2) = eafww(ecenter,1,3,ww,eaft) !line 13
     eaf(3) = eafww(ecenter,1,4,ww,eaft) !line 14
     fafac= factri(eaf,kbt)
  elseif( ecenter <= ww(3)) then
     s2= (ww(2)-ww(1))/(ww(3)-ww(1)+eps) !s2 is the iso-energy surface area at e=ww(2). 13-14-2
     s3= (ww(4)-ww(3))/(ww(4)-ww(2)+eps) !s3 at e=ww(3) as well. 3-24-14
     eaf(1) = eafww(ecenter,1,3,ww,eaft) !Triangle 13-14-23
     eaf(2) = eafww(ecenter,1,4,ww,eaft)
     eaf(3) = eafww(ecenter,2,3,ww,eaft)
     s2e= s2*(ecenter-ww(3))**2/(ww(3)-ww(2))**2 ! area of 13-14-23
     fafac1 = factri(eaf,kbt)* s2e
     eaf(1) = eafww(ecenter,2,4,ww,eaft) ! Triangle 24-14-23
     s3e= s3*(ecenter-ww(2))**2/(ww(3)-ww(2))**2 ! area of 24-14-23
     fafac2 = factri(eaf,kbt)* s3e
     fafac=  (fafac1+fafac2)/(s2e+s3e)
  elseif( ecenter <= ww(4)) then
     eaf(1) = eafww(ecenter,1,4,ww,eaft)
     eaf(2) = eafww(ecenter,2,4,ww,eaft)
     eaf(3) = eafww(ecenter,3,4,ww,eaft)
     fafac= factri(eaf,kbt)
  else
     call rx('eaf_triangle: ecenter else')
  endif
  !$$$      write(6,"('xxx center=',4d19.4)") ecenter
  !$$$      write(6,"('xxx www000=',4d19.4)") ww(1:4)
  !$$$      write(6,"('xxx eaft00=',4d19.4)") eaft(1:4)
  !$$$      write(6,"('xxx eaf000=',3d19.4)") eaf(1:3)
end subroutine eaf_triangle

real(8) function eafww(ecenter,i,j, ww, eaft)
  implicit none
  integer,intent(in):: i,j
  real(8),intent(in):: ww(4),eaft(4)
  real(8)::wgt,eps=1d-10,ecenter
  wgt = (ecenter-ww(i))/(ww(j)-ww(i)+eps)
  eafww = eaft(i)* (1d0-wgt) + eaft(j)* wgt
END PROGRAM

real(8) function factri(eaf, kbt)
  !! Triangle weight for the domain of eaf<0. ---> we will make it for finite temperature version
  !! Temperature factor kbt will be implemented by H.Okumura soon!
  !! Return f(eaf)-f(efb) factor
  !!  input: kbt: Temperature in Ry.
  !!         eaf:  ea - efermi in Ry.
  implicit none
  real(8):: eps=1d-10!to avoid ZeroDIV
  real(8),intent(in):: eaf(3)
  real(8),intent(in):: kbt
  integer:: isig,ieaord(3)
  real(8) :: e1,e2,e3
  call sortea( eaf,ieaord, 3 ,isig)
  e1= eaf(ieaord(1)) ! Note that e1,e2,e3 are just nicknames for variables.
  e2= eaf(ieaord(2)) ! e1 < e2 <e3
  e3= eaf(ieaord(3))
  if(     0d0 <= e1 ) then
     factri = 0d0
  elseif( 0d0 < e2  ) then
     factri = e1**2/(e2-e1+eps)/(e3-e1+eps)
  elseif( 0d0 < e3  ) then
     factri = 1d0 - e3**2/(e3-e2+eps)/(e3-e1+eps)
  else
     factri = 1d0
  endif
END PROGRAM

