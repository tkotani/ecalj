      subroutine vldau(UH,JH,vrsion,dmatu,l,nsp,lmaxu,iblu,Eorb,vorb)
C- Set up LDA+U potential from U and J for one l
C ----------------------------------------------------------------------
Ci Inputs
Ci   UH    :Hubbard U
Ci   JH    :Hubbard J
Ci   vrsion:LDA+U version 
Ci         :1 AMF; 2 FLL; see Remarks
Ci         :4 majority and minority shifted by U
Ci         :4 majority spin shifted by U, minority by J
Ci         :add 100's digit: make Eorb only; do not update vorb
Ci   dmatu :density matrix for LDA+U, spherical harmonics
Ci   l     :l block for which LDA+U is defined
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmaxu :dimensioning parameter for U matrix
Ci   iblu  :index to current LDA+U block
Co Outputs
Co   vorb  :orbital dependent-potential matrices
Co   Eorb  :orbital energy
Cl Local variables
Cr Remarks
Cr   See Petukhov      PRB 67, 153106 (2003) for AMF, FLL in spherical approx
Cr   See Liechtenstein PRB 52, R5467 (1995) for FLL limit
Cu Updates
Cu   06 May 07 Bug fix: return dmatu unchanged in AFM limit
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   29 Oct 05 Switch to evaluate Etot without updating vorb
Cu   27 Apr 05 Lambrecht first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,nsp,lmaxu,iblu,vrsion
      double precision UH,JH,Eorb
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
C ... Local parameters
      integer m,m1,m2,m3,isp,lvrs,lnov
      double precision Eldau,Edc,Ueff,dmat4
      double precision trace(2),n0(2),ttrace
      double complex Vnew(-3:3,-3:3)
      double precision Vee(-3:3,-3:3,-3:3,-3:3)
      double complex vtemp1,vtemp2
      integer iot(2)
c      double complex tracev

C      integer stdo,nglob
C      stdo = nglob('stdo')

      lvrs = mod(vrsion,100)
      lnov = mod(vrsion/100,10)

c     iot(i) is spin complement to spin i
      iot(1) = 2
      iot(2) = 1

C     Calculate n_sig = Tr(rho_sig)/(2l+1)
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

C     AMF
C see Petukhov  PRB 67, 153106 (2003) but generalize ala Liechtenstein
C for nonspherical case
C AMF means construct with delta dmatu instead of dmatu in what follows
      if (lvrs .eq. 1) then
        do  isp = 1, nsp
          do m = -l, l
            dmatu(m,m,isp,iblu) = dmatu(m,m,isp,iblu) - n0(isp)
          enddo
        enddo
      endif
C     End AMF
C     Eq. 6, PRB 67, 153106 (2003)
C     Note:  JH=0 => Vee(m,m2,m1,m3) = U delta(m,m1) delta(m2,m3)
      call veecomp(Vee,l,UH,JH)
C      print *, 'vee'      
C      do m=-l,l
C        do m1=-l,l
C          print *,'m,m1',m,m1
C          do  m2=-l,l
C            print '(7f10.4)', (Vee(m,m2,m1,m3),m3=-l,l)
C          enddo
C        enddo
C      enddo
C      stop
C     See Liechtenstein PRB 52, R5467 (1995) for FLL limit
      Eldau = 0d0
      do  isp = 1, 2
        do  m = -l, l
          do  m1 = -l, l
            Vnew(m,m1) = (0d0,0d0)
            do  m2 = -l, l
              do  m3 = -l, l

C             Case lvrs = 4:
C             Potential shift and dmat input, not U and dmat
C             Replace dmat with diagonal dmat4 => V diagonal, l independent
              if (lvrs .eq. 4) then
                if (m.ne.m1 .or. m2.ne.m3 .or. m.ne.m2) goto 10

                dmat4 = n0(isp)/2 + n0(3-isp)/2
C               Spherical average from Petukhov
C               Vnew(m,m) = - Ueff (dmat_eff - 0.5d0)
                Vnew(m,m1) = Vnew(m,m1) + UH 
                Ueff = -UH/(dmat4-0.5d0)
C               Petukhov Eq. 3.  Factor of 1/2 comes later
                Eldau = Eldau - Ueff*(dmat4**2 - dmat4)

C               Mimic mode 2
C               U2 = 0
C               if (m .eq. m3) U2 = Ueff
C               print *, Ueff/2*(n0(isp)**2 - n0(isp))
C               print *, 
C     .           Ueff*dmat4*dmat4 + (Ueff - U2)*dmat4*dmat4
C               Eldau = Eldau +
C     .           Ueff*dmat4*dmat4 + (Ueff - U2)*dmat4*dmat4
              elseif (lvrs .eq. 5) then 
                if (m .eq. m1) Vnew(m,m1)=UH*(isp-1.5d0)*2 
                Eldau = 0
              else
C             First line in Eq. 5, PRB 52, R5467
C             NB: J=0 => vtemp1 = U delta(m,m1) delta(m2,m3)
              vtemp1 = Vee(m,m2,m1,m3)*dmatu(m2,m3,iot(isp),iblu)
C             Second and third lines in Eq. 5, PRB 52, R5467
              vtemp2 = (Vee(m,m2,m1,m3) - Vee(m,m2,m3,m1))*
     .                  dmatu(m2,m3,isp,iblu)
              Vnew(m,m1) = Vnew(m,m1) + vtemp1 + vtemp2

              Eldau = Eldau + Vee(m,m2,m1,m3)*dmatu(m,m1,isp,iblu)*
     .          dmatu(m2,m3,iot(isp),iblu) + (Vee(m,m2,m1,m3) -
     .          Vee(m,m2,m3,m1))*dmatu(m,m1,isp,iblu)*
     .          dmatu(m2,m3,isp,iblu)
              endif
   10         continue
              enddo
            enddo
          enddo

C         FLL: see last line, Eq. 5, PRB 52, R5467
          if (lvrs .eq. 2) then
            Vnew(m,m) = Vnew(m,m) -
     .                  UH*(ttrace-0.5d0) + JH*(trace(isp)-0.5d0)
c            print *, vnew(m,m)
          endif
C         End FLL
        enddo

        if (lnov .eq. 0) then
          do  m=-l,l
            do  m1 = -l,l
              Vorb(m,m1,isp,iblu) = Vnew(m,m1)
            enddo
          enddo
        endif

      enddo

      Eldau = Eldau/2d0
      Edc = 0d0
      if (lvrs .eq. 2) then
        Edc = 0.5d0*UH*ttrace*(ttrace-1d0) - 0.5d0*JH*
     .    (trace(1)*(trace(1)-1d0) + trace(2)*(trace(2)-1d0))
      endif
      Eorb = Eldau - Edc
      if (lvrs .eq. 4) then
        call info5(20,0,0,' vldau:  Eldau = %,6;6d  '//
     .    'Ueff = %,6;6d  Eorb = %,6;6d  ',Eldau,Ueff,Eorb,0,0)
      else
        call info5(20,0,0,' vldau:  Eldau = %,6;6d  '//
     .    'Edc = %,6;6d  Eorb = %,6;6d  ',Eldau,Edc,Eorb,0,0)
      endif

C      if (l .eq. 2)  then
C        print 333, (dble(vnew(m,m)), m=-l,l)
C  333   format(5f12.5)
C      stop
C      endif


C     Restore dmatu in AMF case
      if (lvrs .eq. 1) then
        do  isp = 1, nsp
          do m = -l, l
            dmatu(m,m,isp,iblu) = dmatu(m,m,isp,iblu) + n0(isp)
          enddo
        enddo
      endif
C     End AMF


      end


