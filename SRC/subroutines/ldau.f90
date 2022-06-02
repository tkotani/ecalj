subroutine ldau(vrsion,l,iblu,UH,JH,dmatu,nsp,lmaxu,vorb,Eorb)
  use m_lmfinit,only:stdo
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
  if(iprint()>0) write(stdo,ftox)' ldau: version(See ldau.F)= ',lvrs,'iblu=', iblu
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

