module m_ldauu
  public ldau
  private
  contains
subroutine ldau(vrsion,l,iblu,UH,JH,dmatu,nsp,lmaxu,vorb,Eorb)
  use m_lgunit,only:stdo
  use m_ftox
  intent(out)                                         vorb,Eorb
  intent(in)      vrsion,l,iblu,UH,JH,dmatu,nsp,lmaxu
  !- Makes Vorb and Eorb from dmatu for given site and l
  ! ----------------------------------------------------------------------
  !i Inputs
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
  double precision :: UH,JH,Eorb
  double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
  double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
  integer :: m1,m2,ii,isp,lvrs,lnov
  double precision :: aaa,num(2),aven(2),nnum,bot
  double precision :: E1,E2
  double complex Vorb1(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
  double complex Vorb2(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
  double complex den1(-3:3,-3:3),den2(-3:3,-3:3)
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
  use m_ftox
  use m_lgunit,only:stdo
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
  double precision :: UH,JH,Eorb
  double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
  double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
  ! ... Local parameters
  integer :: m,m1,m2,m3,isp,lvrs,lnov
  double precision :: Eldau,Edc,Ueff,dmat4
  double precision :: trace(2)=1d99,n0(2),ttrace
  double complex Vnew(-3:3,-3:3)
  double precision :: Vee(-3:3,-3:3,-3:3,-3:3)
  double complex vtemp1,vtemp2
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
10               continue
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
  double precision :: U_H,J_H
  double precision :: Vee(-3:3,-3:3,-3:3,-3:3)
  ! ... Local parameters
  integer :: m,m1,m2,m3,k,jj,sg
  double precision :: F(0:2*3),r42,r62,t1,t2,t3,ak
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
  double precision :: sgJ,sgk,C,ab,t3j
  double precision :: a1,a2,a3,a4,a5,a6,a7,a8,a9,b
  double precision :: a1a,a2a,a3a,a4a,a5a,a6a,a7a,a8a,a9a
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
  double precision :: f,f2(0:100),fctor2
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

end module m_ldauu
