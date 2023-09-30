!>utils for LDAU
module m_ldau_util
  use m_lgunit,only:stdo
  use m_MPItk,only: master_mpi
  use m_lmfinit,only: lmxa_i=>lmxa
  use m_ftox
  public ldau, sudmtu,chkdmu,mixmag,symdmu
  private
contains
  subroutine ldau(vrsion,l,iblu,UH,JH,dmatu,nsp,lmaxu,vorb,Eorb) !Makes Vorb and Eorb from dmatu for given site and l
    implicit none
    intent(out)                                       vorb,Eorb
    intent(in)    vrsion,l,iblu,UH,JH,dmatu,nsp,lmaxu
    !i   vrsion:LDA+U version 1 AMF; 2 FLL; 3 mixed Petukhov version
    !i         :see Remarks
    !i         :add 100's digit: make Eorb only; do not update vorb
    !i   l     :l block for which LDA+U is defined
    !i   iblu  :index to current LDA+U block
    !i   UH    :Hubbard U
    !i   JH    :Hubbard J
    !i   dmatu :density matrix for LDA+U
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmaxu :dimensioning parameter for U matrix
    !o Outputs
    !o   vorb  :orbital-dependent potential matrices
    !o   Eorb  :orbital energy
    !l Local variables
    !r Remarks
    !r   See Liechtenstein PRB 52, R5467 (1995) for FLL limit
    !r   See Petukhov      PRB 67, 153106 (2003) for AMF
    !r                                           Eq. 5 for mixed
    !u Updates
    !u   09 Nov 05 (wrl) Convert dmat to complex form
    !u   29 Oct 05 Switch to evaluate Etot without updating vorb
    !u   27 Apr 05 Lambrecht first created
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: l,iblu,nsp,lmaxu,vrsion,iprint
    real(8) :: UH,JH,Eorb
    complex(8) vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
    complex(8) dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
    integer :: m1,m2,ii,isp,lvrs,lnov
    real(8) :: aaa,num(2),aven(2),nnum,bot
    real(8) :: E1,E2
    complex(8) Vorb1(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
    complex(8) Vorb2(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
    complex(8) den1(-3:3,-3:3),den2(-3:3,-3:3)
    lvrs = mod(vrsion,100)
    lnov = mod(vrsion/100,10)
    !if(iprint()>0) write(stdo,ftox)' ldau: version(See ldau.F)= ',lvrs,'iblu=', iblu
    ! see Petukhov et al. PRB 67, 153106 (2003) construct aaa=alpha mixing
    ! of two ldau versions eq. 5
    if (lvrs == 3) then
       do  isp = 1, 2
          aven(isp) = 0.0d0
          do  m1 = -l,l
             aven(isp) = aven(isp) + dmatu(m1,m1,isp,iblu)
          enddo
          aven(isp) = aven(isp)/(2*l+1)
          do  m1 = -l, l
             do  m2 = -l, l
                den1(m1,m2) = dmatu(m1,m2,isp,iblu)
                den2(m1,m2) = 0.d0
             enddo
             den1(m1,m1) = dmatu(m1,m1,isp,iblu) - aven(isp)
          enddo
          do  m1 = -l,l
             do  m2 = -l,l
                do  ii = -l,l
                   den2(m1,m2) = den2(m1,m2) + den1(m1,ii)* &
                        den1(m2,ii)
                enddo
             enddo
          enddo
          num(isp) = 0.0d0
          do  m1 = -l,l
             num(isp) = num(isp) + den2(m1,m1)
          enddo
       enddo
       nnum = 0d0
       bot = 0d0
       do  isp = 1, 2
          nnum = nnum + num(isp)
          bot = bot + aven(isp)*(1d0 - aven(isp))
       enddo
       if (bot == 0d0) stop 'LDAU: divide by zero bot'
       aaa = nnum/((2*l+1)*bot)

       !       call two types of vorb  and average them weighted according to aaa
       call vldau(UH,JH,1,dmatu,l,nsp,lmaxu,iblu,E1,Vorb1)
       call vldau(UH,JH,2,dmatu,l,nsp,lmaxu,iblu,E2,Vorb2)
       if (lnov == 0) then
          do  isp = 1, 2
             do  m1 = -l, l
                do  m2 = -l, l
                   vorb(m1,m2,isp,iblu) = (1-aaa)*Vorb1(m1,m2,isp,iblu) + &
                        aaa*Vorb2(m1,m2,isp,iblu)
                enddo
             enddo
          enddo
       endif
       Eorb = (1-aaa)*E1 + aaa*E2
    else
       call vldau(UH,JH,vrsion,dmatu,l,nsp,lmaxu,iblu,Eorb,vorb)
    endif
  end subroutine ldau
  subroutine vldau(UH,JH,vrsion,dmatu,l,nsp,lmaxu,iblu,Eorb,vorb)
    intent(out)                                        Eorb,vorb
    !- Set up LDA+U potential from U and J for one l
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   UH    :Hubbard U
    !i   JH    :Hubbard J
    !i   vrsion:LDA+U version
    !i         :1 AMF; 2 FLL; see Remarks
    !i         :4 majority and minority shifted by U
    !i         :4 majority spin shifted by U, minority by J
    !i         :add 100's digit: make Eorb only; do not update vorb
    !i   dmatu :density matrix for LDA+U, spherical harmonics
    !i   l     :l block for which LDA+U is defined
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmaxu :dimensioning parameter for U matrix
    !i   iblu  :index to current LDA+U block
    !o Outputs
    !o   vorb  :orbital dependent-potential matrices
    !o   Eorb  :orbital energy
    !l Local variables
    !r Remarks
    !r   See Petukhov      PRB 67, 153106 (2003) for AMF, FLL in spherical approx
    !r   See Liechtenstein PRB 52, R5467 (1995) for FLL limit
    !u Updates
    !u   06 May 07 Bug fix: return dmatu unchanged in AFM limit
    !u   09 Nov 05 (wrl) Convert dmat to complex form
    !u   29 Oct 05 Switch to evaluate Etot without updating vorb
    !u   27 Apr 05 Lambrecht first created
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: l,nsp,lmaxu,iblu,vrsion
    real(8) :: UH,JH,Eorb
    complex(8) vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
    complex(8) dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
    ! ... Local parameters
    integer :: m,m1,m2,m3,isp,lvrs,lnov
    real(8) :: Eldau,Edc,Ueff,dmat4
    real(8) :: trace(2)=1d99,n0(2),ttrace
    complex(8) Vnew(-3:3,-3:3)
    real(8) :: Vee(-3:3,-3:3,-3:3,-3:3)
    complex(8) vtemp1,vtemp2
    integer :: iot(2)
    lvrs = mod(vrsion,100)
    lnov = mod(vrsion/100,10)
    !     iot(i) is spin complement to spin i
    iot(1) = 2
    iot(2) = 1
    !     Calculate n_sig = Tr(rho_sig)/(2l+1)
    do  isp = 1, nsp
       trace(isp) = 0d0
       do  m = -l, l
          trace(isp) = trace(isp) + dmatu(m,m,isp,iblu)
       enddo
    enddo
    ttrace = trace(1) + trace(2)
    do  isp = 1, nsp
       n0(isp) = trace(isp)/(2*l+1)
    enddo
    !     AMF
    ! see Petukhov  PRB 67, 153106 (2003) but generalize ala Liechtenstein
    ! for nonspherical case
    ! AMF means construct with delta dmatu instead of dmatu in what follows
    if (lvrs == 1) then
       do  isp = 1, nsp
          do m = -l, l
             dmatu(m,m,isp,iblu) = dmatu(m,m,isp,iblu) - n0(isp)
          enddo
       enddo
    endif
    !     End AMF
    !     Eq. 6, PRB 67, 153106 (2003)
    !     Note:  JH=0 => Vee(m,m2,m1,m3) = U delta(m,m1) delta(m2,m3)
    call veecomp(Vee,l,UH,JH)
    !      print *, 'vee'
    !      do m=-l,l
    !        do m1=-l,l
    !          print *,'m,m1',m,m1
    !          do  m2=-l,l
    !            print '(7f10.4)', (Vee(m,m2,m1,m3),m3=-l,l)
    !          enddo
    !        enddo
    !      enddo
    !      stop
    !     See Liechtenstein PRB 52, R5467 (1995) for FLL limit
    Eldau = 0d0
    do  isp = 1, 2
       do  m = -l, l
          do  m1 = -l, l
             Vnew(m,m1) = (0d0,0d0)
             do  m2 = -l, l
                do  m3 = -l, l

                   !             Case lvrs = 4:
                   !             Potential shift and dmat input, not U and dmat
                   !             Replace dmat with diagonal dmat4 => V diagonal, l independent
                   if (lvrs == 4) then
                      if (m /= m1 .OR. m2 /= m3 .OR. m /= m2) goto 10

                      dmat4 = n0(isp)/2 + n0(3-isp)/2
                      !               Spherical average from Petukhov
                      !               Vnew(m,m) = - Ueff (dmat_eff - 0.5d0)
                      Vnew(m,m1) = Vnew(m,m1) + UH
                      Ueff = -UH/(dmat4-0.5d0)
                      !               Petukhov Eq. 3.  Factor of 1/2 comes later
                      Eldau = Eldau - Ueff*(dmat4**2 - dmat4)

                      !               Mimic mode 2
                      !               U2 = 0
                      !               if (m .eq. m3) U2 = Ueff
                      !               print *, Ueff/2*(n0(isp)**2 - n0(isp))
                      !               print *,
                      !     .           Ueff*dmat4*dmat4 + (Ueff - U2)*dmat4*dmat4
                      !               Eldau = Eldau +
                      !     .           Ueff*dmat4*dmat4 + (Ueff - U2)*dmat4*dmat4
                   elseif (lvrs == 5) then
                      if (m == m1) Vnew(m,m1)=UH*(isp-1.5d0)*2
                      Eldau = 0
                   else
                      !             First line in Eq. 5, PRB 52, R5467
                      !             NB: J=0 => vtemp1 = U delta(m,m1) delta(m2,m3)
                      vtemp1 = Vee(m,m2,m1,m3)*dmatu(m2,m3,iot(isp),iblu)
                      !             Second and third lines in Eq. 5, PRB 52, R5467
                      vtemp2 = (Vee(m,m2,m1,m3) - Vee(m,m2,m3,m1))* &
                           dmatu(m2,m3,isp,iblu)
                      Vnew(m,m1) = Vnew(m,m1) + vtemp1 + vtemp2
                      Eldau = Eldau + Vee(m,m2,m1,m3)*dmatu(m,m1,isp,iblu)* &
                           dmatu(m2,m3,iot(isp),iblu) + (Vee(m,m2,m1,m3) - &
                           Vee(m,m2,m3,m1))*dmatu(m,m1,isp,iblu)* &
                           dmatu(m2,m3,isp,iblu)
                   endif
10                 continue
                enddo
             enddo
          enddo
          if (lvrs == 2) then !   FLL: see last line, Eq. 5, PRB 52, R5467
             Vnew(m,m) = Vnew(m,m) - UH*(ttrace-0.5d0) + JH*(trace(isp)-0.5d0)
          endif
       enddo
       if (lnov == 0) then
          do  m=-l,l
             do  m1 = -l,l
                Vorb(m,m1,isp,iblu) = Vnew(m,m1)
             enddo
          enddo
       endif
    enddo
    Eldau = Eldau/2d0
    Edc = 0d0
    if (lvrs == 2) then
       Edc = 0.5d0*UH*ttrace*(ttrace-1d0) - 0.5d0*JH* &
            (trace(1)*(trace(1)-1d0) + trace(2)*(trace(2)-1d0))
    endif
    Eorb = Eldau - Edc
    !if (lvrs == 4) then
    !   write(stdo,ftox)' vldau:  Eldau =',ftof(Eldau),'Ueff=',Ueff,'Eorb =',Eorb
    !else
    !   write(stdo,ftox)' vldau:  Eldau =',ftof(Eldau),'Edc =',Edc,'Eorb=',Eorb
    !endif
    if (lvrs == 1) then !!   Restore dmatu in AMF case
       do  isp = 1, nsp
          do m = -l, l
             dmatu(m,m,isp,iblu) = dmatu(m,m,isp,iblu) + n0(isp)
          enddo
       enddo
    endif
  end subroutine vldau

  subroutine Veecomp(Vee,l,U_H,J_H)
    !- Calculate Vee from U, J, Slater integrals, and Clebsch-Gordon coefficients.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   l     :l for which U_H and J_H are defined
    !i   U_H   :Screened direct coulomb integral (Hubbard U)
    !i   J_H   :Exchange integral (Hubbard J)
    !o Outputs
    !o   Vee   :Vee(m,m2,m3,m4) = <m,m2|Vee|m1,m3>
    !l Local variables
    !r Remarks
    !u Updates
    !u   02 Jun 05 Lambrecht bug fixes
    !u   27 Apr 05 Larson first created
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: l
    real(8) :: U_H,J_H
    real(8) :: Vee(-3:3,-3:3,-3:3,-3:3)
    ! ... Local parameters
    integer :: m,m1,m2,m3,k,jj,sg
    real(8) :: F(0:2*3),r42,r62,t1,t2,t3,ak
    !     Slater integrals from U_H, J_H (this will be updated later)

    call dpzero(vee,7**4)
    call dpzero(F,7)

    F(0) = U_H
    !     d-states -> r42 = F(4)/F(2) = 0.625
    if (l == 2) then
       r42 = 0.625d0
       F(2) = 14d0*J_H/(1d0+r42)
       F(4) = F(2)*r42
    endif
    !     f-states
    if (l == 3) then
       r42 = 451d0/675d0
       r62 = 1001d0/2025d0
       F(2) = 3d0*J_H/(2d0/15d0+1d0/11d0*r42+50d0/429d0*r62)
       F(4) = F(2)*r42
       F(6) = F(2)*r62
    endif
    !     Calculate <m,m2|Vee|m1,m3>
    do  m = -l, l
       do  m2 = -l, l
          do  m1 = -l, l
             do  m3 = -l, l
                do  k = 0, 2*l,2
                   !                     print *,'k=',k,'F=',F(k)
                   ak = 0d0
                   call t3j_all(l,k,l,0,0,0,t1)
                   !                     print *,'l=',l,'k=',k,'F=',F(k),'t1=',t1
                   do  jj = -k, k
                      call t3j_all(l,k,l,-m,jj,m1,t2)
                      call t3j_all(l,k,l,-m2,-jj,m3,t3)
                      sg = (-1)**(m+jj+m2)
                      ak = ak + sg*t2*t3
                   enddo
                   ak = (2*l+1)**2*t1**2*ak
                   !                     print '(5i3,3f8.2)',
                   !     .               m,m2,m1,m3,k,Vee(m,m2,m1,m3),ak,F(k)
                   Vee(m,m2,m1,m3) = Vee(m,m2,m1,m3) + ak*F(k)
                   !                     print '(5i3,3f9.4)',
                   !     .               m,m2,m1,m3,k,Vee(m,m2,m1,m3),ak,F(k)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine Veecomp

  subroutine t3j_all(j1,j2,j3,m1,m2,m3,t3j)
    ! ----------------------------------------------------------------
    !- Calculate 3-j symbols.
    !     implicit none
    integer :: j1,j2,j3,m1,m2,m3
    integer :: J,L1,L2,K0,K1,K2,KMAX,KMIN,trr,k
    real(8) :: sgJ,sgk,C,ab,t3j
    real(8) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,b
    real(8) :: a1a,a2a,a3a,a4a,a5a,a6a,a7a,a8a,a9a
    t3j = 0d0
    if ((m1+m2+m3) /= 0) goto 142
    call tri_rule(j1,j2,j3,trr)
    if (trr < 0) goto 142
    call fctor1( j1 + j2 - j3,a1)
    call fctor1( j1 - j2 + j3,a2)
    call fctor1(-j1 + j2 + j3,a3)
    call fctor1( j1 + m1,a4)
    call fctor1( j1 - m1,a5)
    call fctor1( j2 + m2,a6)
    call fctor1( j2 - m2,a7)
    call fctor1( j3 + m3,a8)
    call fctor1( j3 - m3,a9)
    call fctor1( j1 + j2 + j3 + 1,b)
    ab = b
    a1a = a1
    a2a = a2
    a3a = a3
    a4a = a4
    a5a = a5
    a6a = a6
    a7a = a7
    a8a = a8
    a9a = a9
    J   = j1 - j2 - m3
    sgJ = (-1)**(J)
    C   = sgJ*DSQRT(a1a*a2a*a3a*a4a*a5a*a6a*a7a*a8a*a9a/ab)

    K0   = j1 + j2 - j3
    K1   = j1 - m1
    K2   = j2 + m2
    L1   = j2 - j3 - m1
    L2   = j1 - j3 + m2
    KMAX = min0(K1,K2,K0)
    KMIN = max0(0,L1,L2)

    t3j = 0d0
    do  k = KMIN, KMAX
       call fctor1(k,a1)
       call fctor1(j1 +j2 - j3 - k,a2)
       call fctor1(j1 -m1 - k,a3)
       call fctor1(j2 +m2 - k,a4)
       call fctor1(j3 -j2 + m1 + k,a5)
       call fctor1(j3 -j1 - m2 + k,a6)
       sgk = (-1d0)**k
       a1a = a1
       a2a = a2
       a3a = a3
       a4a = a4
       a5a = a5
       a6a = a6
       t3j = t3j + sgk/(a1a*a2a*a3a*a4a*a5a*a6a)
    enddo
    t3j = C*t3j
142 continue
  end subroutine t3j_all
  subroutine fctor1(nfac,fctor2)
    ! ----------------------------------------------------------------
    !     implicit none
    real(8) :: f,f2(0:100),fctor2
    integer :: nfac,i
    f2(0) = 1d0
    f = 0d0
    do  i = 1,nfac
       f = f+1d0
       f2(i) = f2(i-1)*f
    enddo
    fctor2 = f2(nfac)
  end subroutine fctor1

  subroutine tri_rule(j1,j2,j3,tri_ru)
    ! ----------------------------------------------------------------
    !     implicit none
    integer :: j1,j2,j3,tri_ru
    !     check the triangular rule
    tri_ru = 1
    if ((j1+j2-j3) < 0)tri_ru = -1
    if ((j1-j2+j3) < 0)tri_ru = -1
    if ((-j1+j2+j3) < 0)tri_ru = -1
  end subroutine tri_rule

  ! ssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine mixmag(sss)
    use m_amix,only: amix
    !  subroutine pqmixa(nda,nmix,mmix,mxsav,beta,rms2,a,tj)
    !- Mixing routine for sigma. Modified from pqmixa in subs/pqmix.f
    !- Anderson mixing of a vector
    !i  mmix: number of iterates available to mix
    ! o nmix: nmix > 0: number of iter to try and mix
    !i        nmix < 0: use mmix instead of nmix.
    !o  nmix: (abs)  number of iter actually mixed.
    !o        (sign) <0, intended that caller update nmix for next call.
    !  MvS Feb 04 use sigin as input sigma if available (lsigin=T)
    !             Add mixnit as parameter
    implicit none
    !      logical lsigin
    integer,parameter:: nda=1
    integer:: nmix,mmix
    integer,parameter:: mxsav=10
    real(8) :: rms2,tj(mxsav),beta
    integer :: im,imix,jmix,onorm,okpvt,oa !,amix
    integer :: iprintxx,ifi,nitr,ndaf
    real(8)::sss(nda),sigin(nda)
    real(8):: tjmax
    real(8),allocatable::norm(:),a(:,:,:)
    integer,allocatable:: kpvt(:)
    integer::ret
    character(20) :: fff
    logical :: fexist
    real(8):: acc
    integer:: ido,ifile_handle
    iprintxx = 30
    beta=.3d0
    allocate ( a(nda,0:mxsav+1,2) )
    fff="mixmag.aftest"
    INQUIRE (FILE =fff, EXIST = fexist)
    if(fexist)      write(6,*)'... reading file mixsigma'
    if( .NOT. fexist) write(6,*)'... No file mixsigma'
    open(newunit=ifi,file=fff,form='unformatted')
    if(fexist) then
       read(ifi,err=903,end=903) nitr,ndaf
       if (ndaf /= nda) goto 903
       read(ifi,err=903,end=903) a
       goto 902
    endif
    goto 901
903 continue
    print 368
368 format(5x,'(warning) file mismatch ... mixing file not read')
901 continue
    nitr = 0
902 continue
    a(:,0,1) = sss   !output
    if( .NOT. fexist) a(:,0,2) = sss !input
    !     if input sigma available, use it instead of file a(:,0,2)
    !      if (lsigin) then
    !        write(6,*)'... using input sigma read from sigm file'
    !        a(:,0,2) = sigin  !input
    !      endif
    write(6,*)'sum sss=',sum(abs(sss))
    imix=9
    mmix = min(max(nitr-1,0),imix)
    if (mmix > mxsav) mmix = mxsav
    !     this information already printed out by amix
    !     write(6,*)'mixing parameters for amix are fixed in mixsigma'
    !     write(6,*)'   beta       =', beta
    !     write(6,*)'   tjmax      =', tjmax
    !     write(6,*)'   mmix mxsav =', mmix,mxsav
    !     call getkeyvalue("GWinput","mixtj",acc,default=0d0,status=ret)
    acc=0d0
    if(acc/=0d0) then
       write(6,*)' readin mixtj from GWinput: mixtj=',acc
       tjmax=abs(acc)+1d-3
       if(mmix==1) then
          tj(1)=1d0
       else
          tj(1)= acc
          tj(2)= 1-acc
          mmix=2
       endif
       ido=2
    else
       tjmax=5d0
       ido=0
    endif
    !      allocate(norm(mxsav**2),kpvt(mxsav))
    imix = amix(nda,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax, &
         a,tj,rms2) !norm,kpvt,
    !      deallocate(norm, kpvt)
    ! ... Restore PQ array, updating new x
    !      call dpscop(a,w(oa),nda,1+nda*(mxsav+2),1+nda*(mxsav+2),1d0)
    !      call dcopy(nda*(mxsav+2)*2,w(oa),1,a,1)
    ! ...
    sss = a(:,0,2)
    rewind(ifi)
    write(ifi) nitr+1,nda
    write(ifi) a
    close(ifi)
  end subroutine mixmag
  subroutine chkdmu(eks, dmatu,dmatuo,vorb,eorb)
    use m_lmfinit,only: stdl,nbas,nsp,nlibu,lmaxu,ispec,lldau,tolu=>mix_tolu,umix=>mix_umix,stdo,idu,uh,jh,ham_lsig,addinv
    use m_mksym,only: g=>symops,istab=>oistab, ng =>ngrp
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    !use m_ldauu,only: ldau
    implicit none
    intent(in)::       eks,       dmatuo
    intent(out)::                        vorb,eorb
    ! LDA+U potential vorb and total energy eorb.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlibu: number of U blocks
    !i   lmaxu :dimensioning parameter for U matrix
    !i   dmatu : dmatu produced in current iteration
    !i         : dmatu is passed in real harmonics
    !i   dmatuo: dmatu produced in prior iteration
    !i         : dmatuo is passed in real harmonics
    ! xxx   idvsh=0 dmatu, dmatuo, vorb input/output in real harmonics
    !i   umix  :linear mixing parameter for density matrix
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat in dmats(*,lldau(ib))
    !i   ng    :number of group operations
    !i   g     :point group operations
    !i   istab :site istab(i,ig) is transformed into site i by grp op ig
    !o  vorb  :orbital dependent potential matrices
    !o  eorb  : U contribution to LDA+U total energy
    ! ----------------------------------------------------------------------
    integer:: ierr
    include "mpif.h"
    integer:: idvsh=0
    real(8):: eks, eorbxxx,eorb
    complex(8) dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    complex(8) dmatuo(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    complex(8) vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    integer :: l,lmxa,ib,is,iblu,igetss,idmat,ivsiz,ifile_handle !,idu(4)
    integer :: iprint,ipl,havesh
    real(8) :: ddmat,eorbi,eterms(20),ddot,xx !,uh(4),jh(4)
    logical:: fexist,mmtargetx,eee
    real(8),allocatable:: uhall(:,:)
    real(8):: mmsite(nbas),uhxx,mmhist(10000),uhhist(10000),mmtarget,uhdiff,ddo
    real(8),save::uhxnew,uhx,alpha,alphax
    integer:: nn,ibas,ifx,key,i,nit
    integer,save::ncount=0
    complex(8):: dmatuav(-lmaxu:lmaxu,-lmaxu:lmaxu)
    real(8):: sss(1),sigin
    real(8):: fac,ssss
    if (nlibu == 0) return
    ssss=0d0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib) !ssite(ib)%spec
          ssss = ssss + sum(abs(uh(:,is)))+sum(abs(jh(:,is)))
       endif
    enddo
    havesh =0
    idvsh  =0  ! We assume real harmonics for i/o
    ipl = 1
    ivsiz = nsp*nlibu*(lmaxu*2+1)**2
    ! --- Symmetrize output dmatu (req. real harmonics); compare diff ---
    if(iprint()>=60)call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,'Unsymmetrized out dmats',dmatu)
    call symdmu(nlibu,dmatu, nbas,nsp, lmaxu, ng, g, istab, lldau, xx)
    if(addinv) dmatu= dreal(dmatu) !this is symmetrization: needed because we skip phi_-k for some
    ! k points because of extra symmetry of inversion for k points. 2023-jan
    if(master_mpi)write(stdo,ftox)
    if(master_mpi)write(stdo,ftox)'chkdmu: LDA+U. RMSdiff of dmat from symmetrization =',ftod(xx,2)
    ! --- Compute U contribution to total energy; make vorb ---
    call rotycs(1,dmatu,nbas,nsp,lmaxu,lldau) !mode=1, dmatu is from rh to sh
    havesh = 1
    eorb = 0
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib) 
          lmxa=lmxa_i(is)
          do  l = 0, min(lmxa,3)
             if (idu(l+1,is) /= 0) then
                iblu = iblu+1
                eorbi = 999
                call ldau(100+idu(l+1,is),l,iblu,uh(l+1,is),jh(l+1,is),dmatu,nsp,lmaxu,vorb,eorbi)
                eorb = eorb + eorbi
             endif
          enddo
       endif
    enddo
    ! --- LDA total energy terms ---
    if(master_mpi)write(stdo,ftox)'eks =',ftof(eks),'e[U]=',ftof(eorb),'Etot(LDA+U)=',ftof(eks+eorb)
    if(master_mpi)write(stdl,ftox)'ldau EHK',ftof(eks),'U',ftof(eorb),'ELDA+U',ftof(eks+eorb)
    ! --- Restore dmatu, vorb to real harmonics
    call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau) !-1, from sh to rh idvsh=0
    havesh = 0
    if(master_mpi)then
       ddmat = sum(abs(dmatu-dmatuo)**2)**.5
       write(stdo,ftox)'LDA+U update density matrix ... RMS diff in densmat',ftod(ddmat)
    endif
    ddo = sum(abs(dmatuo)**2) !dmatuo=0d0 if no dmatu.* occnum.*
    if(ddo>1d-10.and.ssss>1d-10) dmatu = umix*dmatu+(1d0-umix)*dmatuo ! new*umix + old*(1-umix)
    call rotycs(1,dmatu,nbas,nsp,lmaxu,lldau) !from rh to sh
    havesh = 1
    vorb=1d99
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib)
          lmxa=lmxa_i(is)
          do  l = 0, min(lmxa,3)
             if (idu(l+1,is) /= 0) then
                iblu = iblu+1
                call pshpr(iprint()-20)
                uhx = uh(l+1,is)
                call ldau(idu(l+1,is),l,iblu,uhx,jh(l+1,is),dmatu,nsp, lmaxu,vorb,eorbxxx) 
                call poppr
             endif
          enddo
       endif
    enddo
    !  ... Symmetrize vorb to check (symdmu requires real harmonics)
    call rotycs(-1,vorb,nbas,nsp,lmaxu,lldau) !-1 means that vorb is converted from sh to rh
    call symdmu(nlibu,vorb,nbas , nsp , lmaxu , ng , g , istab , lldau , xx )
    call rotycs(1,vorb,nbas,nsp,lmaxu,lldau) !1 means that vorb is converted from rh to sh
    !!=>  At this point, dmatu and vorb are in spherical harmonics
    if(Iprint()>20) write(stdo,ftox)'RMS change in vorb from symmetrization =',ftod(xx)
    if(xx>.0001d0 .AND. iprint()>30) write(stdo,'(a)')'(warning) RMS change unexpectely large'
    !    if(iprint()>0) write(6,ftox)'=== representation in spherical harmonics dmatu ==='
    !    call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau, 'Mixed dmats',dmatu)
    !    if(iprint()>0) write(6,ftox)'=== representation in spherical harmonics vorb ==='
    !    call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau, 'New vorb',vorb)
    if (master_mpi) then
       idmat = ifile_handle()
       open(idmat,file='dmats.'//trim(sname)) !havesh mush be 1
       if(iprint()>0) write(stdo,*)' ... Writing density matrix dmats.*... '
       !write(idmat,ftox)'# === dmats in sph harmonics (read by lmfp-m_ldau_init-sudmtu)==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,' dmats !spherical harmonics',dmatu)
       write(idmat,ftox)'# --- We do not read following data. Just for check ---'
       write(idmat,ftox)'# vorb in sph harmonics (not readin) ==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau, 'New vorb',vorb)
    endif
    call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau) !dmatu is converted from sh to rh
    call rotycs(-1,vorb,nbas,nsp,lmaxu,lldau)  !vorb is converted from sh to rh
    havesh=0                  !I recovered this 2022May8
    if (master_mpi) then ! write in real harmonics
       write(idmat,ftox)'# dmatu in real harmonics (not readin) ==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,'Mixed dmats',dmatu)
       write(idmat,ftox)'# vorb  in real harmonics (not readin) ==='
       call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,'New vorb',vorb)
       close(idmat)
    endif
  end subroutine chkdmu
  subroutine sudmtu(dmatu,vorb) !not touch module variables
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    use m_lmfinit,only: nbas,nsp,nlibu,lmaxu,lldau,ispec,stdo,slabl,idu,uh,jh
    use m_mksym,only: g=>symops,istab=>oistab, ng =>ngrp
    !- Initialize site density matrix and vorb  for LDA+U
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlibu : nlibu total number of U blocks
    !i   lmaxu :dimensioning parameter for U matrix
    !i   idvsh :0 dmatu and vorb returned in real harmonics
    !i         :1 dmatu and vorb returned in spherical harmonics
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat in dmats(*,lldau(ib))
    !i   ng    :number of group operations
    !i   g     :point group operations
    !i   istab :site istab(i,ig) is transformed into site i by grp op ig
    !o Outputs
    !o   dmatu :density matrix for LDA+U orbitals
    !o         :in real spherical harmonics basis
    !o   vorb  :orbital dependent potential matrices
    !o         :in real spherical harmonics basis
    !l Local variables
    !l   eorb  : U contribution to LDA+U total energy
    !r Remarks
    !r   Reads in diagonal occupation numbers from file occnum.ext or dmatu
    !r   given in order of m=-l,l, isp=1,2, and constructs initial vorb
    !    use m_ldauu,only: ldau
    implicit none
    integer:: idvsh=0
    complex(8):: dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),&
         Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    logical :: mmtargetx,eee
    integer :: i,isp,ib,l,lmxa,m,m2,foccn,havesh,ivsiz,ipr,ifx
    real(8):: nocc(-3:3,2),iv(7),eorb,xx
    integer:: igetss,is,idmat,fxst,iblu,nlm,nlmu,nn,m1 !idu(4),
    complex(8):: tmp(7,7),img=(0d0,1d0)
    real(8):: tempr(7,7),tempi(7,7)
    character str*80,spid*8,aaa*24, xn*8
    complex(8) :: dmwk_zv(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    real(8):: uhx,uhxx
    ! ... MPI
    include "mpif.h"
    integer :: procid,master,mpipid,ierr
    logical :: mlog,occe,dexist,readtemp,cmdopt0
    real(8)::sss
    character(128):: bbb
    call rxx(nsp.ne.2,'LDA+U must be spin-polarized!')
    procid = mpipid(1)
    master = 0
    call getpr(ipr)
    !! When LDAU is dummy (usually just in order to print our dmats file).

    !NOT set dmatu=0
    ! sss=0d0
    ! do  ib = 1, nbas
    !    if (lldau(ib) /= 0) then
    !       is = ispec(ib) 
    !       sss = sss + sum(abs(uh(:,is)))+sum(abs(jh(:,is)))
    !    endif
    ! enddo
    ! if(sss<1d-6) then
    !    dmatu = 0d0
    !    havesh = 1
    !    goto 1185
    ! endif

    ! Read in dmatu if file  dmats.ext  exists ---
    if(procid /= master) goto 1185
    inquire(file='dmats.'//trim(sname),exist=dexist)
    inquire(file='occnum.'//trim(sname),exist=occe) !if no dmats, try occnum.
    if(dexist) then
       open(newunit=idmat,file='dmats.'//trim(sname))
825    continue
       !       read(idmat,*)str    !bug at 2022-12-11.
       read(idmat,"(a)")str !bug recovered at 2023-01-28
       if(str(1:1) == '#') goto 825
       if(index(str,' sharm ')/=0) then
          havesh = 1 !spherical(complex) harmonics
       else
          havesh = 0 !real harmonics
       endif
       if(havesh ==1) bbb='spherical harmonics' !complex harmonics
       if(havesh ==0) bbb='real harmonics'      !real harmonics 
       if(master_mpi) write(stdo,*)' sudmtu: reading density matrix from file dmats in '//trim(bbb)
       rewind idmat
       iblu = 0
       do  ib = 1, nbas
          if (lldau(ib) /= 0) then
             is = ispec(ib) 
             lmxa=lmxa_i(is)
             do l = 0, min(lmxa,3)
                if (idu(l+1,is) /= 0) then
                   iblu = iblu+1
                   nlm = 2*l+1
                   nlmu = 2*lmaxu+1
                   do  isp = 1, 2
                      readtemp=.false.
                      read(idmat,*,end=9888,err=9888)
                      do m1=1,nlm
                         read(idmat,*,end=9888,err=9888) (tempr(m1,m2),m2=1,nlm)
                      enddo
                      read(idmat,*,end=9888,err=9888)
                      do m1=1,nlm
                         read(idmat,*,end=9888,err=9888) (tempi(m1,m2),m2=1,nlm)
                      enddo
                      readtemp=.true.
9888                  continue
                      if( .NOT. readtemp) call rxi('sudmtu failed to read dmats for site',ib)
                      dmatu(-l:l,-l:l,isp,iblu)=tempr(1:nlm,1:nlm)+img*tempi(1:nlm,1:nlm)
                   enddo
                endif
             enddo
          endif
       enddo
    elseif(occe) then         !if no occunum.* initial dmatu=0
       write(stdo,*) 'sudmtu: initial (diagonal) density-matrix from occ numbers'
       open(newunit=foccn,file='occnum.'//trim(sname))
       havesh = 1
12     continue
       !      read(foccn,*) str  !bug at 2022-12-11.
       read(foccn,"(a)") str  !bug recovered at 2023-01-28
       if (str(1:1) == '#') goto 12
       if (str(1:1) == '%') then
          if(index(str,' real ')/=0) havesh=0
       else
          backspace(foccn)
       endif
       iblu = 0
       dmatu=0d0
       do  ib = 1, nbas
          if (lldau(ib) /= 0) then
             is = ispec(ib)
             lmxa=lmxa_i(is)
             do l = 0,min(lmxa,3)
                if (idu(l+1,is) /= 0) then
                   iblu = iblu+1
                   do  isp = 1, 2
11                    continue
                      if (str(1:1) == '#') goto 11
                      read(foccn,*) nocc(-l:l,isp)
                      write(stdo,ftox)' occnum: site',ib,'l',l,'isp',isp,' ',ftof(nocc(-l:l,isp))
                   enddo
                   do isp = 1, 2
                      do m = -l, l
                         dmatu(m,m,isp,iblu) = nocc(m,isp)
                      enddo
                   enddo
                endif
             enddo
          endif
       enddo
       close(foccn)
    else
       dmatu=0d0
    endif
1185 continue
    ! ... Initial printout
    call praldm(0,51,51,havesh,nbas,nsp,lmaxu,lldau,' dmats read from disk',dmatu)
    ivsiz = nsp*nlibu*(lmaxu*2+1)**2
    call mpibc1(dmatu,2*ivsiz,4,mlog,'sudmtu','dmatu')
    call mpibc1(havesh,1,2,.false.,' ',' ')
    ! ... Density matrix in real or spherical harmonics (fixed by idvsh)
    if (havesh /= idvsh) then
       call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = idvsh
    endif
    ! ... Symmetrize dmatu (symdmu requires real harmonics)
    dmwk_zv=dmatu
    if (havesh == 1) then
       call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = 0
    endif
    call symdmu(nlibu,dmatu , nbas , nsp , lmaxu , ng , g , istab , lldau , xx )
    if (havesh /= idvsh) then
       call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu, lldau)
       call rotycs(2*idvsh-1,dmwk_zv,nbas,nsp,lmaxu, lldau )
       havesh = idvsh
    endif
    if (ng /= 0) then
       if(master_mpi)write(stdo,ftox)'sudmtu:  RMS change in dmats from symmetrization',ftof(xx)
       if (xx > .01d0) write(stdo,*)'(warning) RMS change unexpectely large'
       call daxpy ( ivsiz * 2 , - 1d0 , dmatu , 1 , dmwk_zv , 1 )
       if(ipr>=60) write(stdo,*)' change in dmat wrought by symmetrization'
       call praldm ( 0 , 60 , 60 , 0 , nbas , nsp , lmaxu , lldau ,' ' , dmwk_zv )
    endif
    !    Print dmats in specified harmonics
    dmwk_zv=dmatu
    if (havesh /= idvsh) then
       call rotycs ( 2 * idvsh - 1 , dmwk_zv , nbas , nsp , lmaxu, lldau )
    endif
    if(cmdopt0('--showdmat')) then
       if(master_mpi)write(stdo,*)
       call praldm(0,30,30,idvsh,nbas,nsp,lmaxu,lldau,' Symmetrized dmats' , dmwk_zv )
    endif
    !     Print dmats in complementary harmonics
    i = 1-idvsh
    call rotycs(2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu, lldau )
    if(cmdopt0('--showdmat')) then
       if(master_mpi)write(stdo,*)
       call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau, ' Symmetrized dmats' , dmwk_zv )
    endif
    ! ... Make Vorb (ldau requires spherical harmonics)
    if (havesh /= 1) then
       call rotycs(1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = 1
    endif
    if(master_mpi) write(stdo,*)
    iblu = 0
    do  20  ib = 1, nbas
       if (lldau(ib) == 0) goto 20
       is = ispec(ib) 
       lmxa=lmxa_i(is)
       spid=slabl(is) 
       i = min(lmxa,3)
       if(master_mpi) write(stdo,ftox)'Species '//spid//'mode',idu(1:i+1,is),'U',ftof(uh(1:i+1,is),2),'J',ftof(jh(1:i+1,is),2)
       do  22  l = 0, i
          if (idu(l+1,is) /= 0) then
             iblu = iblu+1
             uhx=uh(l+1,is)
             call ldau(idu(l+1,is),l,iblu,uhx,jh(l+1,is),dmatu,nsp,lmaxu,vorb,eorb)
          endif
22     enddo
20  enddo
    call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,' Unsymmetrized vorb',vorb)
!!!! ==>     At this point, dmatu and vorb are in spherical harmonics (complex)

    ! ... Symmetrize vorb to check (symdmu requires real harmonics)
    call rotycs(-1,vorb,nbas,nsp,lmaxu,lldau)
    call symdmu (nlibu, vorb, nbas , nsp , lmaxu , ng , g , istab , lldau , xx )
    !     EITHER: vorb =>  spherical harmonics OR dmatu => real harmonics
    if (idvsh == 1) then
       call rotycs(1,vorb,nbas,nsp,lmaxu,lldau)
       havesh = 1
    else
       call rotycs(-1,dmatu,nbas,nsp,lmaxu,lldau)
       havesh = 0
    endif
    if (master_mpi.and.ng /= 0) then
       write(stdo,ftox)' sudmtu:  RMS change in vorb from symmetrization = ',ftof(xx)
       if (xx > .01d0) write(stdo,*)'          (warning) RMS change unexpectely large'
       write(stdo,*)
    endif
    if(cmdopt0('--showdmat')) call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,' Symmetrized vorb',vorb)
    i = 1-idvsh
    dmwk_zv=vorb
    call rotycs( 2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , lldau )
    if(cmdopt0('--showdmat')) then
       write(stdo,*) ! vorb in complementary harmonics
       call praldm(0,30,30, i , nbas , nsp , lmaxu , lldau , ' Vorb' , dmwk_zv )
    endif
    eorb = 0d0
    return
99  continue 
    call rx('bad occnum file, site l= '//trim(xn(ib))//trim(xn(l)))
  end subroutine sudmtu
  subroutine rotycs(mode,a,nbas,nsp,lmaxu,lldau) !- Rotate matrix a from real to spherical harmonics for LDA+U objects densmat and vorb
    use m_lmfinit,only:idu,ispec
    !i mode =1 from real a to spherical a
    !i      -1 from spherical a to real a
    !i a: matrix to be transformed a(m,m,isp,iblu)  could be vorb or dmat
    !i nbas : number of sites
    !i nsp  : number of spins
    !i lmaxu: lmax for U
    !i lldau  :lldau(ib)=0 => no U on this site otherwise
    !i        :U on site ib with dmat in dmats(*,lldau(ib))
    !o a rotated in place
    !r Remarks
    !r order of cubic harmonics ls (l-1)s,ms...1s 0 1c mc... (l-1)c lc
    !r order of spherical harmonics -l:l
    !r Yl-m=(Ylmc-iYlms)/sqrt(2)  Ylm=(-1)**m*conjg(Yl-m)
    !u Updates
    !u   18 Jan 06 A. Chantis changed rotation matrices in accordance with
    !u             the definition of real harmonics used in the rest of
    !u             the code (Hund's rules satisfied as indicated by orb. moment)
    !u   09 Nov 05 (wrl) Convert dmat to complex form
    !u   30 Apr 05 Lambrecht first created
    !----------------------------------------------------------------
    implicit none
    integer :: nbas,lldau(nbas),mode,lmaxu,nsp
    complex(8),target:: a(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nbas)
    integer:: ib,m,l,is,igetss,i,j,k,ll,isp,iblu
    complex(8):: b(2,2),c(2,2),add,bb(2,2)
    complex(8),parameter:: s2= 1/dsqrt(2d0), img=(0d0,1d0)
    complex(8),pointer:: rot(:,:)
    complex(8),save,target::rott(2,2,5,-1:1)
    logical,save:: init=.true.
    logical:: cmdopt0
    if(mode/=1 .AND. mode/=-1) call rx('ROTYCS: mode must be 1 or -1')
    if(init)then
       do m=1,5
          rott(1,1:2,m,-1) = [s2,          s2*(-1d0)**m] !mode -1: spherical to cubic basis
          rott(2,1:2,m,-1) = [img*s2, -img*s2*(-1d0)**m]
          rott(:,:,  m, 1) = transpose(dconjg(rott(:,:,m,-1))) !mode 1 : cubic to spherical
       enddo
       init=.false.
    endif
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib)==0) cycle
       is  = ispec(ib) 
       do  l = 0, min(lmxa_i(is),3)
          if (idu(l+1,is) ==0) cycle
          iblu = iblu+1
          do  isp = 1, 2
             do  m = 1, l
                bb(1,:) = [a( m,m,isp,iblu), a( m,-m,isp,iblu)]
                bb(2,:) = [a(-m,m,isp,iblu), a(-m,-m,isp,iblu)]
                rot => rott(:,:,m,mode)
                c = matmul(matmul(rot,bb),transpose(dconjg(rot))) !c=rot*b*rot^+
                a(m,m,isp,iblu)   = c(1,1)
                a(m,-m,isp,iblu)  = c(1,2)
                a(-m,m,isp,iblu)  = c(2,1)
                a(-m,-m,isp,iblu) = c(2,2)
             enddo
          enddo
       enddo
    enddo
  end subroutine rotycs
  subroutine praldm(ifi,ipr1,ipr2,sharm,nbas,nsp,lmaxu,lldau,strn,dmatu) !- Writes out a site density-matrix-like object for all sites
    use m_struc_def
    use m_lmfinit,only:idu,ispec
    !i Inputs
    !i   ifi   :if zero, write to stdo, in screen style format
    !i         :else, write to file ifi in high-precision format
    !i   ipr1  :if verbosity ipr>ipr1, print header
    !i   ipr2  :if verbosity ipr>ipr2, print contents of dmats
    !i   sharm :0 if in real harmonics, 1 if in spherical harmonics
    !i   nbas  :size of basis
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmaxu :dimensioning parameter for U matrix
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i          U on site ib with dmat in dmats(*,lldau(ib))
    !i   strn  :string put into header
    !i   dmatu :density matrix for LDA+U
    !o Outputs
    !o   dmatu is written to file ifi
    implicit none
    integer :: nbas,nsp,lldau(nbas),ifi,lmaxu,ipr1,ipr2,sharm,i_copy_size
    complex(8)::   dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
    character strn*(*)
    integer :: iblu,ib,is,igetss,lmxa,l !,idu(4)
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib) 
          lmxa=lmxa_i(is)
          do  l = 0, min(lmxa,3)
             if (idu(l+1,is) /= 0) then
                iblu = iblu+1
                call prdmts(ifi,ipr1,ipr2,sharm,strn,ib,l,lmaxu,iblu,dmatu, nsp,1)
             endif
          enddo
       endif
    enddo
  end subroutine praldm
  subroutine prdmts(ifi,ipr1,ipr2,sharm,strn,ib,l,lmaxu,iblu,dmats, nsp,nspc)!- Writes out a site density-matrix-like object for a single l
    use m_ftox
    use m_lmfinit,only: stdo
    !i Inputs
    !i   ifi   :if zero, write to stdo, in screen style format
    !i         :else, write to file ifi in high-precision format
    !i   ipr1  :if verbosity ipr>ipr1, print header
    !i   ipr2  :if verbosity ipr>ipr2, print contents of dmats
    !i   sharm :0 if in real harmonics, 1 if in spherical harmonics
    !i   strn  :string put into header
    !i   ib    :site index (ib=0 suppresses printout)
    !i   l     :dmats defined for l block
    !i   lmaxu :dimensioning parameter for dmats
    !i   iblu  :index to current block
    !i   dmats :site density matrix
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
    !o Outputs
    !o  header and dmats are printed to stdo
    implicit none
    integer :: ifi,ipr1,ipr2,l,ib,lmaxu,nsp,nspc,iblu,sharm
    complex(8) :: dmats(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
    character strn*(*)
    integer :: isp,ipr,m1,m2,mpipid,nlm
    character strnl*120,strn1*30,lll*60
    if (nspc == 2) call rx('prdmts not ready for nspc=2')
    if (mpipid(1) /= 0) return
    call getpr(ipr)
    nlm = 2*l+1
    if (ifi /= 0) then
       if(sharm==0) lll=' real    rharm'
       if(sharm/=0) lll=' complex sharm'
       do  isp = 1, nsp
          if (ipr>=ipr1)write(ifi,ftox)'% ',trim(lll),trim(strn),'l=',l,'site',ib,'spin',isp
          if (ipr< ipr2)return
          do  m1 = -l, l
             write(ifi,'(7(f12.7,2x))') (dreal(dmats(m1,m2,isp,iblu)),m2=-l,l)
          enddo
          write(ifi,*)
          do  m1 = -l, l
             write(ifi,'(7(f12.7,2x))') (dimag(dmats(m1,m2,isp,iblu)),m2=-l,l)
          enddo
       enddo
    else
       if(sharm == 0) then
          strnl = strn // ' real harmonics'
       else
          strnl = strn // ' spherical harmonics'
       endif
       do  isp = 1, nsp
          if (ipr < ipr2) return
          write(stdo,ftox) trim(strnl),'l=',l,'ib',ib,'isp',isp
          do  m1 = -l, l
             write(stdo,'(7(f9.5,2x))')(dreal(dmats(m1,m2,isp,iblu)),m2=-l,l)
          enddo
          write(stdo,'(1x)')
          do  m1 = -l, l
             write(stdo,'(7(f9.5,2x))')(dimag(dmats(m1,m2,isp,iblu)),m2=-l,l)
          enddo
       enddo
    endif
  end subroutine prdmts
  subroutine symdmu(nlibu,dmatu,nbas,nsp,lmaxu,ng,g, istab,lldau,rms)!- Symmetrize LDA+U density matrix dmatu
    use m_struc_def
    use m_ftox
    use m_lmfinit,only:idu,ispec
    use m_lgunit,only:stdo
    !i Inputs
    !i   dmatu :density matrix for LDA+U
    !i   dmatw :work array of same dimension as dmatu.
    !i         :Returns original dmatu on output.
    !i   nbas  :size of basis
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmaxu :dimensioning parameter for U matrix
    !i   ng    :number of group operations.  Program does nothing if ng=0
    !i   g     :point group operations
    !i   istab :table of site permutations for each group op (mksym.f,symtbl.f)
    !i   istab :site istab(i,ig) is transformed into site i by grp op ig
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i          U on site ib with dmat beginning at dmats(*,lldau(ib))
    !l Local variables
    !l   ofjbl :number of dens mat arrays preceding current one for
    !l         :current site
    ! o Inputs/Outputs
    ! o dmatu  :density matrix or vorb symmetrized on output.
    ! o        :the output dmatu is the sum of the original dmatw
    ! o        :and the symmetrized dmatu.  If output dmatu is
    ! o        :to be a symmetrized version of the input,
    ! o        :dmatw MUST be initially zero on input
    !o Outputs
    !o   rms   :rms change in dmatu from symmetrization
    !r Notes
    !r   Routine uses real harmonics
    !u Updates
    !u   30 Jan 06 dmats now symmetrized across different sites.
    !u             Unsuccessful attmempt to include spinor rotations
    !u   09 Nov 05 (wrl) Convert dmat to complex form
    !u   30 Apr 05 Lambrecht first created
    !--------------------------------------------------------------
    implicit none
    integer :: nbas,lldau(nbas),ng,nsp,lmaxu,istab(nbas,ng),i_copy_size
    integer :: is,igetss,lmxa,m1,m2,ilm1,ilm2,ib,l,isp,m3,m4,ig,iblu,nlibu,jb,jblu,ofjbl,lwarn
    real(8):: rmat(16,16),r(-3:3,-3:3),ddot,g(9,*),rms,xx
    complex(8):: sdmat(-3:3,-3:3,2,2),&  ! ... for spinor rotations
         dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),&
         dmatw(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    logical :: cmdopt0
    character (40) :: str
    dmatw=0d0
    rms = 0
    if (cmdopt0('--nosymdm') .OR. ng == 0) return
    ! --- Setup for spinor part ---
    lwarn = 0
    do  ig = 1, ng
       xx = ddet33(g(1,ig))
       !       Extract to rdmat pure rotation part of g; make Euler angles
       if (xx > 0) then
          call dpcopy(g(1,ig),rmat,1,9,1d0)
       else
          call dpcopy(g(1,ig),rmat,1,9,-1d0)
       endif
       if (dabs(xx*g(9,ig)-1) > 1d-6) lwarn = lwarn+1
    enddo
    if (lwarn > 0) write(stdo,ftox)'symdmu (warning): ',lwarn,'symops rotate z axis'
    ! --- For each site density-matrix, do ---
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib) /= 0) then
          is = ispec(ib) !ssite(ib)%spec
          lmxa=lmxa_i(is)
          ofjbl = -1
          do  l = 0, min(lmxa,3)
             if (idu(l+1,is) /= 0) then
                iblu = iblu+1
                ofjbl = ofjbl+1
                do  ig = 1, ng
                   jb = istab(ib,ig)
                   jblu = lldau(jb) + ofjbl
                   !               Rotation matrices for spherical harmonics up to f orbitals
                   call ylmrtg(16,g(1,ig),rmat)
                   !               Pick out the one we need
                   ilm1 = l**2
                   do  m1 = -l, l
                      ilm1 = ilm1+1
                      ilm2 = l**2
                      do  m2 = -l, l
                         ilm2 = ilm2+1
                         r(m1,m2) = rmat(ilm1,ilm2)
                      enddo
                   enddo
                   !           ... Spatial rotation: dmatu(iblu) -> sdmat
                   do  isp = 1, nsp
                      do  m1 = -l, l
                         do  m2 = -l, l
                            sdmat(m1,m2,isp,isp) = 0
                            sdmat(m1,m2,isp,3-isp) = 0
                            do  m3 = -l, l
                               do  m4 = -l, l
                                  sdmat(m1,m2,isp,isp) = sdmat(m1,m2,isp,isp) &
                                       + r(m1,m3)*dmatu(m3,m4,isp,iblu)*r(m2,m4)
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                   do   isp = 1, nsp
                      do  m1 = -l, l
                         do  m2 = -l, l
                            dmatw(m1,m2,isp,jblu) = dmatw(m1,m2,isp,jblu) + &
                                 sdmat(m1,m2,isp,isp)/ng
                         enddo
                      enddo
                   enddo
                enddo
             endif
          enddo
       endif
    enddo
    !     Exchange original for symmetrized dmatu
    nlibu = iblu
    is = nsp*nlibu*(lmaxu*2+1)**2
    call dswap(2*is,dmatw,1,dmatu,1)
    !     RMS change in dmatu
    call daxpy(2*is,-1d0,dmatu,1,dmatw,1)
    rms = dsqrt(ddot(2*is,dmatw,1,dmatw,1)/(2*is))
    call daxpy(2*is,1d0,dmatu,1,dmatw,1)
  end subroutine symdmu
  double precision function ddet33(matrix)
    !- Calculates the determinant of a 3X3 matrix
    double precision :: matrix(9),ddot,m1cm2(3)
    call cross(matrix(4),matrix(7),m1cm2)
    ddet33 = ddot(3,matrix(1),1,m1cm2,1)
  end function ddet33
end module m_ldau_util
