      subroutine vlm2us(lmaxu,rmt,idu,lmxa,iblu,vorb,ppnl,vumm)
C- Rotate vorb from (phi,phidot) to (u,s) and store in vumm
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmaxu :dimensioning parameter for U matrix
Ci   lmxa  :augmentation l-cutoff
Ci   vorb  :orbital-dependent potential matrices
Ci   ppnl  :potential parameters
Cio Inputs/Outputs
Cio  iblu  :index to current LDA+U block
Cio        :on input, index to last LDA+U block that was accessed
Cio        :iblu will be incremented to from blocks at this site
Co Outputs
Co   vumm  :vorb for this site in (us) representation
Co         :vumm(m1,m2,1) = <u| vorb(m1,m2) |u>
Co         :vumm(m1,m2,2) = <u| vorb(m1,m2) |s>
Co         :vumm(m1,m2,3) = <s| vorb(m1,m2) |u>
Co         :vumm(m1,m2,4) = <s| vorb(m1,m2) |s>
Co         :vumm(m1,m2,5) = <u| vorb(m1,m2) |z>
Co         :vumm(m1,m2,6) = <s| vorb(m1,m2) |z>
Co         :vumm(m1,m2,7) = <z| vorb(m1,m2) |z>
Co         :vumm(m1,m2,8) = <z| vorb(m1,m2) |u>
Co         :vumm(m1,m2,9) = <z| vorb(m1,m2) |s>
Cu Updates
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   08 Jun 05 (MvS) extended to local orbitals
Cu   30 Apr 05 Lambrecht first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmaxu,lmxa,iblu,idu(4)
      double precision rmt
      integer n0,nppn,nab
      parameter (n0=10,nppn=12,nab=9)
      double precision ppnl(nppn,n0,2)
      double complex Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,*),
     .               vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)
C ... Local parameters
      integer m1,m2,l,i
      double precision phi,dlphi,phip,dlphip,dphi,dphip
      double precision r12,r21,r11,r22,det
      double precision phz,dphz
      double complex vzz,vuz,vsz,vzu,vzs

C     call prmx('vorb',vorb(1,1,2,1),2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

C ... Rotate Vorb from phi,phidot basis to u,s basis
      do  l = 0, min(lmxa,3)
        if (idu(l+1) .ne. 0) then
          iblu = iblu+1
          do  i = 1, 2
            dlphi  = ppnl(3,l+1,i)/rmt
            dlphip = ppnl(4,l+1,i)/rmt
            phi    = ppnl(5,l+1,i)
            phip   = ppnl(6,l+1,i)
            dphi   = phi*dlphi/rmt
            dphip  = dlphip/rmt*phip
            det = phi*dphip - dphi*phip
            r11 = dphip/det
C           r12 = -dphi/det
            r21 = -phip/det
C           r22 = phi/det
            do  m1 = -l, l
              do  m2 = -l, l
                vumm(m1,m2,1,i,l) = Vorb(m1,m2,i,iblu)*r11*r11
                vumm(m1,m2,2,i,l) = Vorb(m1,m2,i,iblu)*r11*r21
                vumm(m1,m2,3,i,l) = Vorb(m1,m2,i,iblu)*r21*r11
                vumm(m1,m2,4,i,l) = Vorb(m1,m2,i,iblu)*r21*r21
              enddo
            enddo

            phz  = ppnl(11,l+1,i)
            dphz = ppnl(12,l+1,i)

            if (phz .ne. 0) then
              do  m1 = -l, l
                do  m2 = -l, l
                  vzz = phz**2*vumm(m1,m2,1,i,l) +
     .              phz*dphz*(vumm(m1,m2,2,i,l)+vumm(m1,m2,3,i,l)) +
     .              dphz**2*vumm(m1,m2,4,i,l)
                  vuz = - phz*vumm(m1,m2,1,i,l) - dphz*vumm(m1,m2,2,i,l)
                  vsz = - phz*vumm(m1,m2,3,i,l) - dphz*vumm(m1,m2,4,i,l)
                  vzu = - phz*vumm(m1,m2,1,i,l) - dphz*vumm(m1,m2,3,i,l)
                  vzs = - phz*vumm(m1,m2,2,i,l) - dphz*vumm(m1,m2,4,i,l)

                  vumm(m1,m2,5,i,l) = vuz
                  vumm(m1,m2,6,i,l) = vsz
                  vumm(m1,m2,7,i,l) = vzz
                  vumm(m1,m2,8,i,l) = vzu
                  vumm(m1,m2,9,i,l) = vzs

                enddo
              enddo
            endif

C            if (l .eq. 3) then
C              print *, 'Vorb before rotation isp=',i,'l=',l
C              print ('(7f8.4)'),((Vorb(m1,m2,i,iblu),m2=-l,l),m1=-l,l)
C              print *, 'rotated vumm'
C              print ('(7(2f8.4,1x))'),
C     .          (((vumm(m1,m2,k,i,l),m2=-l,l),m1=-l,l),k=1,4)
C              stop
C            endif

          enddo
        endif
      enddo
      end

