subroutine vldau(UH,JH,vrsion,dmatu,l,nsp,lmaxu,iblu,Eorb,vorb)
  intent(out)                                          vorb,Eorb
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
  double precision :: trace(2),n0(2),ttrace
  double complex Vnew(-3:3,-3:3)
  double precision :: Vee(-3:3,-3:3,-3:3,-3:3)
  double complex vtemp1,vtemp2
  integer :: iot(2)
  !      double complex tracev

  !      integer stdo,nglob
  !      stdo = nglob('stdo')

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

        !         FLL: see last line, Eq. 5, PRB 52, R5467
        if (lvrs == 2) then
           Vnew(m,m) = Vnew(m,m) - &
                UH*(ttrace-0.5d0) + JH*(trace(isp)-0.5d0)
           !            print *, vnew(m,m)
        endif
        !         End FLL
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
  if (lvrs == 4) then
     call info5(20,0,0,' vldau:  Eldau = %,6;6d  '// &
          'Ueff = %,6;6d  Eorb = %,6;6d  ',Eldau,Ueff,Eorb,0,0)
  else
     call info5(20,0,0,' vldau:  Eldau = %,6;6d  '// &
          'Edc = %,6;6d  Eorb = %,6;6d  ',Eldau,Edc,Eorb,0,0)
  endif

  !      if (l .eq. 2)  then
  !        print 333, (dble(vnew(m,m)), m=-l,l)
  !  333   format(5f12.5)
  !      stop
  !      endif


  !     Restore dmatu in AMF case
  if (lvrs == 1) then
     do  isp = 1, nsp
        do m = -l, l
           dmatu(m,m,isp,iblu) = dmatu(m,m,isp,iblu) + n0(isp)
        enddo
     enddo
  endif
  !     End AMF


end subroutine vldau


