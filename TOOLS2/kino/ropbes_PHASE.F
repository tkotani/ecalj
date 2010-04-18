C ----------------------------------------------------------------
C   spherical Bessel function at x=r(i)*sqrt(e) divided by x**l 
c
Ci Inputs
Ci   r    list of points
Ci   e    energy
Ci   y,h  work vectors of length n each 
Co Outputs
Co   xi   J(r,l)/r**l, according to standard definition
C
c origianl: Feb. 15, 2010, Hiori Kino
      subroutine ropbes(r,e,lmax,y,h,xi,n)
      implicit none
      integer lmax,n
      double precision e,r(n),xi(n,0:lmax),h(n),y(n)

      double precision:: e2,eps=1.0d-10
      integer:: i,l
      e2=sqrt(e)
      do i=1,n
         h(i)= r(i)*e2
         if (h(i)< eps) h(i)=eps 
      enddo
      do l=0,lmax
        call dsjnv(l,n,h,y)
        do i=1,n
        xi(i,l)=y(i)/h(i)**l 
        enddo
      enddo   
 
      end subroutine ropbes

C
C dsjnv and dsjnvn is a part of PHASE v7.00.
C
C --------------------------------------------------------------------------------
C 
C COPYRIGHT of the program codes
C Copyright (C) 1993-2007 Hideki Katagiri, Koichi Kato, Tsuyoshi Miyazaki,
C Yoshitada Morikawa, Hideaki Sawada, Toshihiro Uchiyama, Tsuyoshi Uda,
C Takahiro Yamasaki, Noriaki Hamada, Akira Yanase, Takenori Yamamoto, 
C Hideaki Tsukioka, Masakuni Okamoto, Hideo Mizouchi, Kiyoshi Betsuyaku 
C and Kazuki Mae.
C 
C It is understood by the authors that the Institute of Industrial Science (IIS),
C the University of Tokyo, distributes this program as "RSS21 Free Software" with
C users' agreement with the terms and conditions written in the file, LICENSE.pdf
C or LICENSE_J.pdf (in Japanese).
C 
C 
C HISTORY
C The original version of this set of the computer programs "PHASE"
C was developed by the members of the Theory Group of Joint Research
C Center for Atom Technology (JRCAT), based in Tsukuba, in the period
C 1993-2001. The names of the contributors to the original version are Hideki 
C Katagiri, K. Kato, T. Miyazaki, Y. Morikawa, H. Sawada, T. Uchiyama,
C T. Uda and T. Yamasaki.
C Since 2002, this set has been tuned and new functions have been
C added to it as a part of the national project "Frontier Simulation
C Software for Industrial Science (FSIS)",  which is supported by
C the IT program of the Ministry of Education, Culture, Sports,
C Science and Technology (MEXT) of Japan. The program was developed further
C mainly by T. Yamasaki. T. Uda, T. Yamamoto, H. Tsukioka, M. Okamoto,
C H. Mizouchi, K. Betsuyaku and K. Mae contributed to the improvement of the code.
C The tetrahedron interpolation codes developed by N. Hamada, A. Yanase 
C and Kiyoyuki Terakura was included. 
C The symmetrization code developed by A. Yanase and N. Hamada was also included. 
C The manual and tutorial were written by Makoto Itoh with the cooperation by 
C Mineo Saito, H. Tsukioka, T. Yamamoto and T. Yamasaki.
C The sample calculations were prepared by T. Yamamoto, H. Tsukioka 
C and Hiroyoshi Momida. 
C Since 2006, this program set has been developed as a part of the national
C project "Revolutionary Simulation Software (RSS21)", which is supported by
C the next-generation IT program of MEXT of Japan.
C The activity of "Multiscale Simulation System for Function Analysis of
C Nanomaterials", RSS21, is supervised by Takahisa Ohno.
C 
C 
C CONTACT ADDRESS
C "Revolutionary Simulation Software" project,
C The Institute of Industrial Science (IIS), The University of Tokyo
C 4-6-1 Komaba, Meguro-ku, Tokyo 153-8505, Japan
C FAX +81-(0)3-5452-6662
C E-mail software@rss21.iis.u-tokyo.ac.jp
C URL http://www.rss21.iis.u-tokyo.ac.jp
C 
C * When distributing RSS21 Software duplications, the user must attach
C the full text in this file.
C 
C --------------------------------------------------------------------------------

#define NONE
        subroutine dsjnv(n,k,x,y)
!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 7.00)
!        
!  SUBROUINE: dsjnv
!     
!  AUTHOR(S): M. Okamoto and T. Yamamoto   December/24/2004
!     
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!     
!=======================================================================
        implicit none
        integer,intent(in)  :: n,k
        real(8),intent(in)  :: x(k)
        real(8),intent(out) :: y(k)
        integer :: ik

        real(8) :: x2,x4
        real(8),parameter :: EPS = 1.0d0
        real(8),parameter :: 
     & cc10 =  1.0000000000000000000d0, 
     & cc20 = -1.6666666666666666667d-1, 
     & cc30 =  8.3333333333333333333d-3, 
     & cc40 = -1.9841269841269841270d-4, 
     & cc50 =  2.7557319223985890653d-6, 
     & cc60 = -2.5052108385441718775d-8, 
     & cc70 =  1.6059043836821614599d-10, 
     & cc80 = -7.6471637318198164759d-13
        real(8),parameter :: 
     & cc11 =  3.3333333333333333333d-1, 
     & cc21 = -3.3333333333333333333d-2, 
     & cc31 =  1.1904761904761904762d-3, 
     & cc41 = -2.2045855379188712522d-5, 
     & cc51 =  2.5052108385441718775d-7, 
     & cc61 = -1.9270852604185937519d-9, 
     & cc71 =  1.0706029224547743066d-11, 
     & cc81 = -4.4983316069528332211d-14
        real(8),parameter :: 
     & cc12 =  6.6666666666666666667d-2, 
     & cc22 = -4.7619047619047619048d-3, 
     & cc32 =  1.3227513227513227513d-4, 
     & cc42 = -2.0041686708353375020d-6, 
     & cc52 =  1.9270852604185937519d-8, 
     & cc62 = -1.2847235069457291680d-10, 
     & cc72 =  6.2976642497339665096d-13, 
     & cc82 = -2.3675429510278069585d-15
        real(8),parameter :: 
     & cc13 =  9.5238095238095238095d-3, 
     & cc23 = -5.2910052910052910053d-4, 
     & cc33 =  1.2025012025012025012d-5, 
     & cc43 = -1.5416682083348750015d-7, 
     & cc53 =  1.2847235069457291680d-9, 
     & cc63 = -7.5571970996807598115d-12, 
     & cc73 =  3.3145601314389297419d-14, 
     & cc83 = -1.1274014052513366469d-16
        real(8),parameter :: 
     & cc14 =  1.0582010582010582011d-3, 
     & cc24 = -4.8100048100048100048d-5, 
     & cc34 =  9.2500092500092500093d-7, 
     & cc44 = -1.0277788055565833344d-8, 
     & cc54 =  7.5571970996807598115d-11, 
     & cc64 = -3.9774721577267156903d-13, 
     & cc74 =  1.5783619673518713057d-15, 
     & cc84 = -4.9017452402232028126d-18
        real(8),parameter :: 
     & cc15 =  9.6200096200096200096d-5, 
     & cc25 = -3.7000037000037000037d-6, 
     & cc35 =  6.1666728333395000062d-8, 
     & cc45 = -6.0457576797446078492d-10, 
     & cc55 =  3.9774721577267156903d-12, 
     & cc65 = -1.8940343608222455668d-14, 
     & cc75 =  6.8624433363124839376d-17, 
     & cc85 = -1.9606980960892811250d-19
        real(8),parameter :: 
     & cc16 =  7.4000074000074000074d-6, 
     & cc26 = -2.4666691333358000025d-7, 
     & cc36 =  3.6274546078467647095d-9, 
     & cc46 = -3.1819777261813725522d-11, 
     & cc56 =  1.8940343608222455668d-13, 
     & cc66 = -8.2349320035749807252d-16, 
     & cc76 =  2.7449773345249935751d-18, 
#ifdef NONE
     & cc86 = -7.2618448003306708335d-21
        real(8),parameter :: 
     & cc17 =  4.9333382666716000049d-7, 
     & cc27 = -1.4509818431387058838d-8, 
     & cc37 =  1.9091866357088235313d-10, 
     & cc47 = -1.5152274886577964534d-12, 
     & cc57 =  8.2349320035749807252d-15, 
     & cc67 = -3.2939728014299922901d-17, 
     & cc77 =  1.0166582720462939167d-19, 
     & cc87 = -2.5040844139071278736d-22
        real(8),parameter :: 
     & cc18 =  2.9019636862774117676d-8, 
     & cc28 = -7.6367465428352941253d-10, 
     & cc38 =  9.0913649319467787206d-12, 
     & cc48 = -6.5879456028599845801d-14, 
     & cc58 =  3.2939728014299922901d-16, 
     & cc68 = -1.2199899264555527000d-18, 
     & cc78 =  3.5057181794699790231d-21, 
     & cc88 = -8.0776916577649286246d-24
#else
     & cc86 = -7.2618448003306708335d-21
#endif

        select case (abs(n))
        case (0)
           do ik=1,k
              if (abs(x(ik)) < EPS) then
                 x2=x(ik)*x(ik)
                 y(ik) = (cc10+x2*(cc20+x2*(cc30+x2*(cc40 
     .                +x2*(cc50+x2*(cc60+x2*(cc70+x2*cc80 
     .                 )))))))
              else
                 y(ik) = sin(x(ik))/x(ik)
              end if
           end do
        case (1)
           do ik=1,k
              x2=x(ik)*x(ik)
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc11+x2*(cc21+x2*(cc31+x2*(cc41 
     .                +x2*(cc51+x2*(cc61+x2*(cc71+x2*cc81 
     .                 )))))))*x(ik)
              else
                 y(ik) = (sin(x(ik))-x(ik)*cos(x(ik)))/x2
              end if
           end do
        case (2)
           do ik=1,k
              x2=x(ik)*x(ik)
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc12+x2*(cc22+x2*(cc32+x2*(cc42 
     .                +x2*(cc52+x2*(cc62+x2*(cc72+x2*cc82 
     .                 )))))))*x2
              else
                 y(ik) = ((3.d0-x2)*sin(x(ik))-3.d0*x(ik)*cos(x(ik)))/(x(ik)*x2)
              end if
           end do
        case (3)
           do ik=1,k
              x2=x(ik)*x(ik)
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc13+x2*(cc23+x2*(cc33+x2*(cc43 
     .                +x2*(cc53+x2*(cc63+x2*(cc73+x2*cc83 
     .                 )))))))*x2*x(ik)
              else
                 y(ik) = ((15.d0-6.d0*x2)*sin(x(ik)) 
     .                  -(15.d0-x2)*x(ik)*cos(x(ik)))/(x2*x2)
              end if
           end do
        case (4)
           do ik=1,k
              x2=x(ik)*x(ik)
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc14+x2*(cc24+x2*(cc34+x2*(cc44 
     .                +x2*(cc54+x2*(cc64+x2*(cc74+x2*cc84 
     .                 )))))))*x2*x2
              else
                 y(ik) = ((105.d0+(-45.d0+x2)*x2)*sin(x(ik)) 
     .                  -(105.d0-10.d0*x2)*x(ik)*cos(x(ik)))/(x(ik)*x2*x2)
              end if
           end do
        case (5)
           do ik=1,k
              x2=x(ik)*x(ik)
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc15+x2*(cc25+x2*(cc35+x2*(cc45 
     .                +x2*(cc55+x2*(cc65+x2*(cc75+x2*cc85 
     .                 )))))))*x2*x2*x(ik)
              else
                 y(ik) = ((945.d0+(-420.d0+15.d0*x2)*x2)*sin(x(ik)) 
     .                  -(945.d0+(-105.d0+x2)*x2)*x(ik)*cos(x(ik)))/(x2*x2*x2)
              end if
           end do
        case (6)
           do ik=1,k
              x2=x(ik)*x(ik)
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc16+x2*(cc26+x2*(cc36+x2*(cc46 
     .                +x2*(cc56+x2*(cc66+x2*(cc76+x2*cc86 
     .                 )))))))*x2*x2*x2
              else
                 y(ik) = ((10395.d0+(-4725.d0+(210.d0-x2)*x2)*x2)*sin(x(ik)) 
     .                  -(10395.d0+(-1260.d0+21.d0*x2)*x2)*x(ik)*cos(x(ik)) 
     .                )/(x(ik)*x2*x2*x2)
              end if
           end do
#ifdef NONE
! spherical bessels with l=7,8 are not performed in double precision.
        case (7)
           do ik=1,k
              x2=x(ik)*x(ik)
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc17+x2*(cc27+x2*(cc37+x2*(cc47 
     .                +x2*(cc57+x2*(cc67+x2*(cc77+x2*cc87 
     .                 )))))))*x2*x2*x2*x(ik)
              else
                 x4=x2*x2
                 y(ik) = ((135135.d0+(-62370.d0+(3150.d0-28.d0*x2)*x2)*x2)*sin(x(ik)) 
     .                  -(135135.d0+(-17325.d0+(378.d0-x2)*x2)*x2)*x(ik)*cos(x(ik)) 
     .                 )/(x4*x4)
              end if
           end do
        case (8)
           do ik=1,k
              x2=x(ik)*x(ik)
              x4=x2*x2
              if (abs(x(ik)) < EPS) then
                 y(ik) = (cc18+x2*(cc28+x2*(cc38+x2*(cc48 
     .               +x2*(cc58+x2*(cc68+x2*(cc78+x2*cc88 
     .                 )))))))*x4*x4
              else
                 y(ik) = ( 
     .         (2027025.d0+(-945945.d0+(51975.d0+(-630.d0+x2)*x2)*x2)*x2)*sin(x(ik)) 
     .         - (2027025.d0+(-270270.d0+(6930.d0-36.d0*x2)*x2)*x2)*x(ik)*cos(x(ik)) 
     .                 )/(x(ik)*x4*x4)
              end if
           end do
#endif
        case default
           y(1:k) = 0.d0
           write(*,*) 'n=',n,', is not supported in dsjnv'
           stop
        end select

        end subroutine dsjnv

        subroutine dsjnvn(n,k,idp,x,y)
!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 7.00)
!        
!  SUBROUINE: dsjnvn
!     
!  AUTHOR(S): M. Okamoto and T. Yamamoto   December/24/2004
!     
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
!
!  FURTHER MODIFICATION: T. Yamasaki, March/01/2006
!   idp is introduced
!     
!=======================================================================
        implicit none
        integer,intent(in)  :: n,k,idp
        real(8),intent(in)  :: x(k)
        real(8),intent(out) :: y(k)
        integer :: ik, i

        real(8) :: x2,x4
        real(8),parameter :: EPS = 1.0d0
        real(8),parameter :: 
     & cc10 =  1.0000000000000000000d0, 
     & cc20 = -1.6666666666666666667d-1, 
     & cc30 =  8.3333333333333333333d-3, 
     & cc40 = -1.9841269841269841270d-4, 
     & cc50 =  2.7557319223985890653d-6, 
     & cc60 = -2.5052108385441718775d-8, 
     & cc70 =  1.6059043836821614599d-10, 
     & cc80 = -7.6471637318198164759d-13
        real(8),parameter :: 
     & cc11 =  3.3333333333333333333d-1, 
     & cc21 = -3.3333333333333333333d-2, 
     & cc31 =  1.1904761904761904762d-3, 
     & cc41 = -2.2045855379188712522d-5, 
     & cc51 =  2.5052108385441718775d-7, 
     & cc61 = -1.9270852604185937519d-9, 
     & cc71 =  1.0706029224547743066d-11, 
     & cc81 = -4.4983316069528332211d-14
        real(8),parameter :: 
     & cc12 =  6.6666666666666666667d-2, 
     & cc22 = -4.7619047619047619048d-3, 
     & cc32 =  1.3227513227513227513d-4, 
     & cc42 = -2.0041686708353375020d-6, 
     & cc52 =  1.9270852604185937519d-8, 
     & cc62 = -1.2847235069457291680d-10, 
     & cc72 =  6.2976642497339665096d-13, 
     & cc82 = -2.3675429510278069585d-15
        real(8),parameter :: 
     & cc13 =  9.5238095238095238095d-3, 
     & cc23 = -5.2910052910052910053d-4, 
     & cc33 =  1.2025012025012025012d-5, 
     & cc43 = -1.5416682083348750015d-7, 
     & cc53 =  1.2847235069457291680d-9, 
     & cc63 = -7.5571970996807598115d-12, 
     & cc73 =  3.3145601314389297419d-14, 
     & cc83 = -1.1274014052513366469d-16
        real(8),parameter :: 
     & cc14 =  1.0582010582010582011d-3, 
     & cc24 = -4.8100048100048100048d-5, 
     & cc34 =  9.2500092500092500093d-7, 
     & cc44 = -1.0277788055565833344d-8, 
     & cc54 =  7.5571970996807598115d-11, 
     & cc64 = -3.9774721577267156903d-13, 
     & cc74 =  1.5783619673518713057d-15, 
     & cc84 = -4.9017452402232028126d-18
        real(8),parameter :: 
     & cc15 =  9.6200096200096200096d-5, 
     & cc25 = -3.7000037000037000037d-6, 
     & cc35 =  6.1666728333395000062d-8, 
     & cc45 = -6.0457576797446078492d-10, 
     & cc55 =  3.9774721577267156903d-12, 
     & cc65 = -1.8940343608222455668d-14, 
     & cc75 =  6.8624433363124839376d-17, 
     & cc85 = -1.9606980960892811250d-19
        real(8),parameter :: 
     & cc16 =  7.4000074000074000074d-6, 
     & cc26 = -2.4666691333358000025d-7, 
     & cc36 =  3.6274546078467647095d-9, 
     & cc46 = -3.1819777261813725522d-11, 
     & cc56 =  1.8940343608222455668d-13, 
     & cc66 = -8.2349320035749807252d-16, 
     & cc76 =  2.7449773345249935751d-18, 
#ifdef NONE
     & cc86 = -7.2618448003306708335d-21
        real(8),parameter :: 
     & cc17 =  4.9333382666716000049d-7, 
     & cc27 = -1.4509818431387058838d-8, 
     & cc37 =  1.9091866357088235313d-10, 
     & cc47 = -1.5152274886577964534d-12, 
     & cc57 =  8.2349320035749807252d-15, 
     & cc67 = -3.2939728014299922901d-17, 
     & cc77 =  1.0166582720462939167d-19, 
     & cc87 = -2.5040844139071278736d-22
        real(8),parameter :: 
     & cc18 =  2.9019636862774117676d-8, 
     & cc28 = -7.6367465428352941253d-10, 
     & cc38 =  9.0913649319467787206d-12, 
     & cc48 = -6.5879456028599845801d-14, 
     & cc58 =  3.2939728014299922901d-16, 
     & cc68 = -1.2199899264555527000d-18, 
     & cc78 =  3.5057181794699790231d-21, 
     & cc88 = -8.0776916577649286246d-24
#else
     & cc86 = -7.2618448003306708335d-21
#endif

        select case (abs(n))
        case (0)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc10+x2*(cc20+x2*(cc30+x2*(cc40 
     & +x2*(cc50+x2*(cc60+x2*(cc70+x2*cc80 
     & )))))))
           end do
           do ik = idp, k
              y(ik) = sin(x(ik))/x(ik)
           end do
        case (1)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc11+x2*(cc21+x2*(cc31+x2*(cc41 
     & +x2*(cc51+x2*(cc61+x2*(cc71+x2*cc81 
     & )))))))*x(ik)
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              y(ik) = (sin(x(ik))-x(ik)*cos(x(ik)))/x2
           end do
        case (2)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc12+x2*(cc22+x2*(cc32+x2*(cc42 
     & +x2*(cc52+x2*(cc62+x2*(cc72+x2*cc82 
     & )))))))*x2
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              y(ik) = ((3.d0-x2)*sin(x(ik))-3.d0*x(ik)*cos(x(ik)))/(x(ik)*x2)
           end do
        case (3)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc13+x2*(cc23+x2*(cc33+x2*(cc43 
     & +x2*(cc53+x2*(cc63+x2*(cc73+x2*cc83 
     & )))))))*x2*x(ik)
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              y(ik) = ((15.d0-6.d0*x2)*sin(x(ik)) 
     .              -(15.d0-x2)*x(ik)*cos(x(ik)))/(x2*x2)
           end do
        case (4)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc14+x2*(cc24+x2*(cc34+x2*(cc44 
     & +x2*(cc54+x2*(cc64+x2*(cc74+x2*cc84 
     & )))))))*x2*x2
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              y(ik) = ((105.d0+(-45.d0+x2)*x2)*sin(x(ik)) 
     & -(105.d0-10.d0*x2)*x(ik)*cos(x(ik)))/(x(ik)*x2*x2)
           end do
        case (5)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc15+x2*(cc25+x2*(cc35+x2*(cc45 
     & +x2*(cc55+x2*(cc65+x2*(cc75+x2*cc85 
     & )))))))*x2*x2*x(ik)
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              y(ik) = ((945.d0+(-420.d0+15.d0*x2)*x2)*sin(x(ik)) 
     & -(945.d0+(-105.d0+x2)*x2)*x(ik)*cos(x(ik)))/(x2*x2*x2)
           end do
        case (6)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc16+x2*(cc26+x2*(cc36+x2*(cc46 
     & +x2*(cc56+x2*(cc66+x2*(cc76+x2*cc86 
     & )))))))*x2*x2*x2
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              y(ik) = ((10395.d0+(-4725.d0+(210.d0-x2)*x2)*x2)*sin(x(ik)) 
     & -(10395.d0+(-1260.d0+21.d0*x2)*x2)*x(ik)*cos(x(ik)) 
     & )/(x(ik)*x2*x2*x2)
           end do
#ifdef NONE
! spherical bessels with l=7,8 are not performed in double precision.
        case (7)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              y(ik) = (cc17+x2*(cc27+x2*(cc37+x2*(cc47 
     & +x2*(cc57+x2*(cc67+x2*(cc77+x2*cc87 
     & )))))))*x2*x2*x2*x(ik)
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              x4=x2*x2
              y(ik) = ((135135.d0+(-62370.d0+(3150.d0-28.d0*x2)*x2)*x2)*sin(x(ik)) 
     & -(135135.d0+(-17325.d0+(378.d0-x2)*x2)*x2)*x(ik)*cos(x(ik)) 
     & )/(x4*x4)
           end do
        case (8)
           do ik=1,idp-1
              x2=x(ik)*x(ik)
              x4=x2*x2
              y(ik) = (cc18+x2*(cc28+x2*(cc38+x2*(cc48 
     & +x2*(cc58+x2*(cc68+x2*(cc78+x2*cc88 
     & )))))))*x4*x4
           end do
           do ik=idp,k
              x2=x(ik)*x(ik)
              x4=x2*x2
              y(ik) = ( 
     & (2027025.d0+(-945945.d0+(51975.d0+(-630.d0+x2)*x2)*x2)*x2)*sin(x(ik)) 
     & - (2027025.d0+(-270270.d0+(6930.d0-36.d0*x2)*x2)*x2)*x(ik)*cos(x(ik)) 
     & )/(x(ik)*x4*x4)
           end do
#endif
        case default
           y(1:k) = 0.d0
           write(*,*) 'n=',n,', is not supported in dsjnvn'
           stop
        end select

      end subroutine dsjnvn

