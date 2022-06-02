subroutine vlm2us(lmaxu,rmt,idu,lmxa,iblu,vorb,ppnl,vumm)
  !- Rotate vorb from (phi,phidot) to (u,s) and store in vumm
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lmaxu :dimensioning parameter for U matrix
  !i   lmxa  :augmentation l-cutoff
  !i   vorb  :orbital-dependent potential matrices
  !i   ppnl  :potential parameters
  ! o Inputs/Outputs
  ! o  iblu  :index to current LDA+U block
  ! o        :on input, index to last LDA+U block that was accessed
  ! o        :iblu will be incremented to from blocks at this site
  !o Outputs
  !o   vumm  :vorb for this site in (us) representation
  !o         :vumm(m1,m2,1) = <u| vorb(m1,m2) |u>
  !o         :vumm(m1,m2,2) = <u| vorb(m1,m2) |s>
  !o         :vumm(m1,m2,3) = <s| vorb(m1,m2) |u>
  !o         :vumm(m1,m2,4) = <s| vorb(m1,m2) |s>
  !o         :vumm(m1,m2,5) = <u| vorb(m1,m2) |z>
  !o         :vumm(m1,m2,6) = <s| vorb(m1,m2) |z>
  !o         :vumm(m1,m2,7) = <z| vorb(m1,m2) |z>
  !o         :vumm(m1,m2,8) = <z| vorb(m1,m2) |u>
  !o         :vumm(m1,m2,9) = <z| vorb(m1,m2) |s>
  !u Updates
  !u   09 Nov 05 (wrl) Convert dmat to complex form
  !u   08 Jun 05 (MvS) extended to local orbitals
  !u   30 Apr 05 Lambrecht first created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lmaxu,lmxa,iblu,idu(4)
  double precision :: rmt
  integer :: n0,nppn,nab
  parameter (n0=10,nppn=12,nab=9)
  double precision :: ppnl(nppn,n0,2)
  double complex Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,*), &
       vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)
  ! ... Local parameters
  integer :: m1,m2,l,i
  double precision :: phi,dlphi,phip,dlphip,dphi,dphip
  double precision :: r12,r21,r11,r22,det
  double precision :: phz,dphz
  double complex vzz,vuz,vsz,vzu,vzs

  !     call prmx('vorb',vorb(1,1,2,1),2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

  ! ... Rotate Vorb from phi,phidot basis to u,s basis
  do  l = 0, min(lmxa,3)
     if (idu(l+1) /= 0) then
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
           !           r12 = -dphi/det
           r21 = -phip/det
           !           r22 = phi/det
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

           if (phz /= 0) then
              do  m1 = -l, l
                 do  m2 = -l, l
                    vzz = phz**2*vumm(m1,m2,1,i,l) + &
                         phz*dphz*(vumm(m1,m2,2,i,l)+vumm(m1,m2,3,i,l)) + &
                         dphz**2*vumm(m1,m2,4,i,l)
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

           !            if (l .eq. 3) then
           !              print *, 'Vorb before rotation isp=',i,'l=',l
           !              print ('(7f8.4)'),((Vorb(m1,m2,i,iblu),m2=-l,l),m1=-l,l)
           !              print *, 'rotated vumm'
           !              print ('(7(2f8.4,1x))'),
           !     .          (((vumm(m1,m2,k,i,l),m2=-l,l),m1=-l,l),k=1,4)
           !              stop
           !            endif

        enddo
     endif
  enddo
end subroutine vlm2us

