      subroutine ldau(vrsion,l,iblu,UH,JH,dmatu,nsp,lmaxu,vorb,Eorb)
C- Makes Vorb and Eorb from dmatu for given site and l
C ----------------------------------------------------------------------
Ci Inputs
Ci   vrsion:LDA+U version 1 AMF; 2 FLL; 3 mixed Petukhov version
Ci         :see Remarks
Ci         :add 100's digit: make Eorb only; do not update vorb
Ci   l     :l block for which LDA+U is defined
Ci   iblu  :index to current LDA+U block
Ci   UH    :Hubbard U
Ci   JH    :Hubbard J
Ci   dmatu :density matrix for LDA+U
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmaxu :dimensioning parameter for U matrix
Co Outputs
Co   vorb  :orbital-dependent potential matrices
Co   Eorb  :orbital energy
Cl Local variables
Cr Remarks
Cr   See Liechtenstein PRB 52, R5467 (1995) for FLL limit
Cr   See Petukhov      PRB 67, 153106 (2003) for AMF
Cr                                           Eq. 5 for mixed
Cu Updates
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   29 Oct 05 Switch to evaluate Etot without updating vorb
Cu   27 Apr 05 Lambrecht first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,iblu,nsp,lmaxu,vrsion
      double precision UH,JH,Eorb
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
C ... Local parameters
      integer m1,m2,ii,isp,lvrs,lnov
      double precision aaa,num(2),aven(2),nnum,bot
      double precision E1,E2
      double complex Vorb1(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
      double complex Vorb2(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,iblu)
      double complex den1(-3:3,-3:3),den2(-3:3,-3:3)
C      integer stdo,nglob
C      stdo = nglob('stdo')

      lvrs = mod(vrsion,100)
      lnov = mod(vrsion/100,10)
      call info2(30,0,0,' ldau:   version = %i  iblu = %i',lvrs,iblu)

C     call prmx('dmatu, spin 1',dmatu,2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

c see Petukhov et al. PRB 67, 153106 (2003) construct aaa=alpha mixing
C of two ldau versions eq. 5
      if (lvrs .eq. 3) then
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
                den2(m1,m2) = den2(m1,m2) + den1(m1,ii)*
     .          den1(m2,ii)
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
        if (bot .eq. 0d0) stop 'LDAU: divide by zero bot'
        aaa = nnum/((2*l+1)*bot)

C       call two types of vorb  and average them weighted according to aaa
        call vldau(UH,JH,1,dmatu,l,nsp,lmaxu,iblu,E1,Vorb1)
        call vldau(UH,JH,2,dmatu,l,nsp,lmaxu,iblu,E2,Vorb2)
        if (lnov .eq. 0) then
          do  isp = 1, 2
            do  m1 = -l, l
              do  m2 = -l, l
                vorb(m1,m2,isp,iblu) = (1-aaa)*Vorb1(m1,m2,isp,iblu) +
     .          aaa*Vorb2(m1,m2,isp,iblu)
              enddo
            enddo
          enddo
        endif
        Eorb = (1-aaa)*E1 + aaa*E2
      else
        call vldau(UH,JH,vrsion,dmatu,l,nsp,lmaxu,iblu,Eorb,vorb)
      endif
      end

