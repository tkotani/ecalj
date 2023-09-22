!>finite-temperature tetrahedron method
module m_tetrakbt 
  use m_keyvalue,only: getkeyvalue
  implicit none
  public:: tetrakbt_init, tetrakbt, kbt,integtetn
  private
  real(8):: tt,kbt, eps=1d-26 !!! tt: temperature [K].  eps is to avoid ZeroDIV, severe condition 
  integer:: i, ndiv !for numerical integration; h=(b-a)/ndiv
  logical:: init = .true.
  real(8):: kb=8.6171d-5 ![eV/K]
  !! developer
  real(8):: err_sum=0d0, err_sum2=0d0
  integer:: nn
contains
  !-----------------------------------------------------
  subroutine tetrakbt_init()
    real(8):: temperature, rydberg
    call getkeyvalue("GWinput","t_tetrakbt",temperature,default=3d+2)
    tt = temperature+1d-12 !avoid 0
    kbt=kb*tt/rydberg()
    if (init) then
       ndiv = ndiv_tt(tt)
       ! ndiv = 10000 ! tmp
       write(6,"(' tetrakbt_init: T[K], kbt[Ry], kbt[eV]',I5,2E13.5)") int(tt),kbt,kbt*rydberg()
       init = .false.
       write(6,*) "ndiv",ndiv
    endif
  end subroutine tetrakbt_init
  !-----------------------------------------------------
  subroutine tetrakbt( ebfin,eafin,v,voltet,frhis,nwhis, wtthis) !this is accumlating variable
!    use m_tetwt5,only: integtetn
    intent(in)::       ebfin,eafin,v,voltet,frhis,nwhis
    intent(out)::                                         wtthis
    !! kbt version
    !! Histgram weights for each bin.
    !! The i-th bin of the Histgrams is [frhis(i) frhis(i+1)].
    !! wtthis(ihis) += \int_{frhis(ihis)^{frhis(ihis+1)} d \omega  \int d^3k (fb-fa) * \delta(\omega +v(k) )
    !!  Total sum of the histgram weight is equal to voltet.
    !!  wtthis(ihis): Note this is accumulating variable.
    real(8):: ebfin(4),eafin(4),v(4), voltet
    real(8):: ebft(4),eaft(4),ebf(3),eaf(3)
    integer::  inihis, iedhis, nwhis,ichk=2313
    integer :: ieaord(1:4),i,isig,idif(3),idf,ix,n,itmp,ihis
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
    ! sorted so that ww(1)<ww(2)<ww(3)<ww(4)

    !! this block is equivalent to sortea(ea,ieaord,n,isig) which gives ieaord (ordering index of ea)
    ! ifdef EXPAND_SORTEA
    n=4
    do i = 1,n
       ieaord(i) = i
    enddo
    do ix= 2,n
       do i=ix,2,-1
          if( -v(ieaord(i-1)) > -v(ieaord(i) ) ) then
             itmp = ieaord(i-1)
             ieaord(i-1) = ieaord(i)
             ieaord(i) = itmp
             !            isig= -isig
             cycle
          endif
          exit
       enddo
    enddo
    ! else
    !      call sortea( -v,ieaord,4,isig)
    ! endif
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
    ! cccccccccccccccccccccccccccccccccc
!!! eaft + ebft = WW ??? --> No
    !$$$      do i = 1,4
    !$$$         write(6,"('--- oku: eaft,ebat,eaft+ebft,WW',i5,4E13.5)") i,eaft(i),ebft(i),eaft(i)+ebft(i),WW(i)
    !$$$      enddo
    ! cccccccccccccccccccccccccccccccccc
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
    intent(in) ::          el,eh, ww,eaft
    intent(out)::                          fafac
    !! We assume ww are soreted as ww(1)<ww(2)<ww(3)<ww(4)
    real(8):: el,eh, ww(4),eaft(4)
    real(8):: fafac
    real(8):: eaf(3),ecenter,fafac2,fafac1,s2,s3,s2e,s3e
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
       fafac= factri(eaf)
    elseif( ecenter <= ww(3)) then
       s2= (ww(2)-ww(1))/(ww(3)-ww(1)+eps) !s2 is the iso-energy surface area at e=ww(2). 13-14-2
       s3= (ww(4)-ww(3))/(ww(4)-ww(2)+eps) !s3 at e=ww(3) as well. 3-24-14
       eaf(1) = eafww(ecenter,1,3,ww,eaft) !Triangle 13-14-23
       eaf(2) = eafww(ecenter,1,4,ww,eaft)
       eaf(3) = eafww(ecenter,2,3,ww,eaft)
       s2e= s2*(ecenter-ww(3))**2/(ww(3)-ww(2))**2 ! area of 13-14-23
       fafac1 = factri(eaf)* s2e
       eaf(1) = eafww(ecenter,2,4,ww,eaft) ! Triangle 24-14-23
       s3e= s3*(ecenter-ww(2))**2/(ww(3)-ww(2))**2 ! area of 24-14-23
       fafac2 = factri(eaf)* s3e
       fafac=  (fafac1+fafac2)/(s2e+s3e)
    elseif( ecenter <= ww(4)) then
       eaf(1) = eafww(ecenter,1,4,ww,eaft)
       eaf(2) = eafww(ecenter,2,4,ww,eaft)
       eaf(3) = eafww(ecenter,3,4,ww,eaft)
       fafac= factri(eaf)
    else
       call rx('eaf_triangle: ecenter else')
    endif
    !$$$      write(6,"('xxx center=',4d19.4)") ecenter
    !$$$      write(6,"('xxx www000=',4d19.4)") ww(1:4)
    !$$$      write(6,"('xxx eaft00=',4d19.4)") eaft(1:4)
    !$$$      write(6,"('xxx eaf000=',3d19.4)") eaf(1:3)
  end subroutine eaf_triangle

  !-------------------------------------------------------------------------
  real(8) function eafww(ecenter,i,j, ww, eaft)
    integer,intent(in):: i,j
    real(8),intent(in):: ww(4),eaft(4)
    real(8)::wgt,ecenter
    wgt = (ecenter-ww(i))/(ww(j)-ww(i)+eps)
    eafww = eaft(i)* (1d0-wgt) + eaft(j)* wgt
  END function eafww

  !-------------------------------------------------------------------------
  real(8) function factri(eaf)
    !! Triangle weight for the domain of eaf<0. ---> we will make it for finite temperature version
    !! Temperature factor kbt will be implemented by H.Okumura soon!
    !! Return f(eaf)-f(efb) factor
    !!  input: kbt: Temperature in Ry.
    !!         eaf:  ea - efermi in Ry.
    real(8),intent(in):: eaf(3)
    ! eal(8),intent(in):: kbt
    real(8):: voltp, s1tp, s2tp ! return s1tp+s2tp
    real(8):: s1int1, s1int2, stp
    integer:: isig,ieaord(3)
    real(8) :: e1,e2,e3, factri_ref
    logical(8):: fdfunc = .false.
    logical(8):: debug = .false.

    call sortea( eaf,ieaord, 3 ,isig)
    e1= eaf(ieaord(1)) ! Note that e1,e2,e3 are just nicknames for variables.
    e2= eaf(ieaord(2)) ! e1 < e2 <e3
    e3= eaf(ieaord(3))

    !! write(6,"('--oku e1, e2-e1, e3-e2, e3-e1: ',4E13.5)") e1,(e2-e1)/e1,(e3-e2)/e1,(e3-e1)/e1
    !! save calculation time
    if     ( 10*kbt <= e1 ) then
       factri = 0d0
    elseif (e3 <= -10*kbt ) then
       factri = 1d0
    else
       !! norm factor for triangle prism (volume of triangle prism)
       voltp=(e2-e1)*(e3-e1)/2d0
       if ( .TRUE. ) then
          !! all integratino
          !!   int_e1^e2 (z-e1)f(z)dz + (e2-e1)/(e3-e2)*int_e2^e3 (e3-z)f(z)dz
          !! = int_e1^e2 (z-e1)f(z)dz - (e2-e1)/(e3-e2)*int_e2^e3 (z-e3)f(z)dz
          !! s1tp=numintall(e1,e2,e1)
          !! s2tp=numintall(e2,e3,e3)
          s1tp=int_simpson(e1,e2,e1)
          s2tp=int_simpson(e2,e3,e3)
          stp = s1tp - s2tp*(e2-e1)/(e3-e2)
          stp = stp/(voltp+eps)
          factri=stp
          !! divided I0 + I1
       else
          !! analitically calculated: int_a^b dx/(exp(x)+1)
          s2tp = (e2-e1)*(e3-e1) + kbt*( e1*( log(exp(e2/kbt)+1) - log(exp(e1/kbt)+1) ) &
               - e3*(e2-e1)*( log(exp(e3/kbt)+1) - log(exp(e2/kbt)+1) )/(e3-e2+eps) )
          s1int1 = integral_0t(e1,e2)
          s1int2 = integral_0t(e2,e3)
          s2tp = s1int1 - s1int2*(e2-e1)/(e3-e2+eps)
          !! numerical integration
          s1int1 = integral_1t(e1,e2)
          s1int2 = integral_1t(e2,e3)
          s1tp = s1int1 - s1int2*(e2-e1)/(e3-e2+eps)
          stp  = (s1tp+s2tp)/(voltp+eps)
          factri=stp
       endif
       !! write(6,"('voltp,s1,s2,s1+s2/voltp',4E13.5)") voltp, s1tp, s2tp, stp
       ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    endif
    !! !! check comparision
    if (debug) then
       factri_ref=factri0(e1,e2,e3)
       factri = factri_ref    !! check Feb13, 2020
       if ( log10(abs(factri-factri_ref)) > -7d0 ) then
          write(6,"('factri,factri0,diff',3E13.5)") factri,factri_ref,factri-factri_ref
       endif
       err_sum  = err_sum  + (factri-factri_ref)
       err_sum2 = err_sum2 + (factri-factri_ref)**2
       nn = nn + 1
       if (nn == 10000 .AND. tt < 1d-6 ) then !! only T=0 K
          write(6,"('nn, error(average)',I8,2E13.5)") nn,err_sum/nn,sqrt(err_sum2/nn-(err_sum/nn)**2)
       endif
       if (nn < 10000) then
          write(6,"('factri,factri0,diff',3E13.5)") factri,factri_ref,factri-factri_ref
       endif
    endif

  end function factri

  !--------------------------------------------------------------
  !! available only 0 K
  real(8) function factri0(e1,e2,e3)
    implicit none
    real(8),intent(in):: e1, e2, e3
    if (    0d0 <= e1) then
       factri0 = 0d0
    elseif( 0d0 < e2  ) then
       factri0 = e1**2/(e2-e1+eps)/(e3-e1+eps)
    elseif( 0d0 < e3  ) then
       factri0 = 1d0 - e3**2/(e3-e2+eps)/(e3-e1+eps)
    elseif ( e3 <= 0d0) then
       factri0 = 1d0
    endif
  end function factri0

  !--------------------------------------------------------------
  !! numerical integration: int x/(exp(x)+1) dx
  real(8) function integral_1t(e_start,e_end)
    !! trapezoidal rule
    implicit none
    real(8),intent(in):: e_start, e_end
    real(8):: hdiv, sumval, x
    real(8) :: intsum !! the value of integral

    hdiv=(e_end-e_start)/float(ndiv)
    sumval=0d0
    integral_1t=0d0
    do i = 1,ndiv
       x=e_start + i*hdiv
       sumval = sumval + 2d0*funcgx(x,1)
    enddo
    integral_1t = (funcgx(e_start,1)+sumval+funcgx(e_end,1))*hdiv/2d0
  end function integral_1t

  !--------------------------------------------------------------
  !! numerical integration: int 1/(exp(x)+1) dx
  real(8) function integral_0t(e_start,e_end)
    !! trapezoidal rule
    implicit none
    real(8),intent(in):: e_start, e_end
    real(8):: hdiv, sumval, x
    real(8) :: intsum !! the value of integral

    hdiv=(e_end-e_start)/float(ndiv)
    sumval=0d0
    integral_0t=0d0
    do i = 1,ndiv
       x=e_start + i*hdiv
       sumval = sumval + 2d0*funcgx(x,0)
    enddo
    integral_0t = (funcgx(e_start,0)+sumval+funcgx(e_end,0))*hdiv/2d0
  end function integral_0t

  !--------------------------------------------------------------
  real(8) function funcgx(x,n)
    implicit none
    integer:: n
    real(8):: x!, kbt
    funcgx = x**n/(exp(x/kbt)+1)
  end function funcgx

  !--------------------------------------------------------------
  real(8) function numintall(e_start,e_end,eeref)
    !! trapezoidal rule
    implicit none
    real(8),intent(in):: e_start, e_end, eeref
    real(8):: hdiv, sumval, x
    real(8) :: intsum !! the value of integral

    hdiv=(e_end-e_start)/float(ndiv)
    sumval=0d0
    numintall=0d0
    do i = 1,ndiv
       x=e_start + i*hdiv
       sumval = sumval + 2d0*funcgx2(x,eeref)
    enddo
    numintall = (funcgx2(e_start,eeref)+sumval+funcgx2(e_end,eeref))*hdiv/2d0
  end function numintall

  !--------------------------------------------------------------
  !! (x-e)f(x)
  real(8) function funcgx2(x,ee)
    implicit none
    real(8):: x, ee
    funcgx2 = (x-ee)/(exp(x/kbt)+1)
  end function funcgx2

  !--------------------------------------------------------------
  integer function ndiv_tt(ttt)
    implicit none
    real(8):: ttt
    integer:: n0 = 100000
    if (ttt < 1d0) then
       ndiv_tt = n0/2d0
    else
       !! ex: ndiv_tt = 10^4 for T=10  K
       !! ex: ndiv_tt = 10^3 for T=100 K
       ndiv_tt = int(n0/ttt)
       !! ndiv_tt must be an even number (for Simpson's rule)
       if (mod(ndiv_tt, 2) == 1) then
          ndiv_tt = ndiv_tt + 1
       endif
    endif
  end function ndiv_tt

  !--------------------------------------------------------------
  real(8) function int_simpson(e_start,e_end,eeref)
    !! Simpson's rule for integration
    implicit none
    real(8),intent(in):: e_start, e_end, eeref
    real(8):: hdiv, sumval, xeven, xodd
    real(8) :: intsum !! the value of integral
    integer:: n2div
    hdiv=(e_end-e_start)/float(ndiv)
    sumval=0d0
    int_simpson=0d0

    ! ndiv must be even for applying Simpson's rule
    n2div=int(ndiv/2d0)
    do i = 1,n2div
       xeven=e_start + 2*i*hdiv
       xodd =e_start + (2*i-1)*hdiv
       if (i==n2div) then !! only odd
          sumval = sumval + 4d0*funcgx2(xodd,eeref)
       else !! even + odd
          sumval = sumval + 4d0*funcgx2(xodd,eeref) + 2d0*funcgx2(xeven,eeref)
       endif
    enddo
    int_simpson = (funcgx2(e_start,eeref)+sumval+funcgx2(e_end,eeref))*hdiv/3d0
  end function int_simpson
  subroutine integtetn(e, ee, integb)  !> Calculate primitive integral of integb = 1/pi Imag[\int^ee dE' 1/(E' -e(k))] = \int^ee dE' S[E']
    !! \remark
    !!  S[E] : is area of the cross-section between the omega-constant plane and the tetrahedron. [here we assumee e1<e2<e3<e4].
    !!  Normalization is not considered! Rath&Freeman Integration of imaginary part on Eq.(17)
    implicit none
    real(8)::  e(1:4), ee, integb,a,b !,e1,e2,e3,e4 V1,V2,V3,V4 , ,D1,D2,D3,D4
    associate( e1=>e(1)-3d-8, e2=>e(2)-2d-8, e3=>e(3)-1d-8, e4=>e(4))
      associate(V1=>ee-e1, V2=>ee-e2, V3=>ee-e3,V4=>ee-e4)
        if(ee<e1) then
           integb=0d0
        elseif( e1<=ee .AND. ee<e2 ) then
           integb = V1**3/((e4-e1)*(e3-e1)*(e2-e1))
        elseif( e2<=ee .AND. ee<e3 ) then
           a  =  V1/ ((e4-e1)*(e3-e1)*(e2-e1))**(1d0/3d0)
           b  =  V2/(-(e4-e2)*(e3-e2)*(e1-e2))**(1d0/3d0)
           integb = (a-b) * (a**2+a*b+b**2)
        elseif( e3<=ee .AND. ee<e4 ) then
           integb = 1d0 - V4**3/((e1-e4)*(e2-e4)*(e3-e4))
        elseif( ee==e4 ) then
           integb = 1d0
        else
           call rx( ' integtetn: ee>e4')
        endif
      endassociate
    endassociate
  end subroutine integtetn
end module m_tetrakbt


