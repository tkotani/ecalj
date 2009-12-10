      subroutine hmfr2c(nbas,nl,indxsh,qsp,eula,neul,ldim,
     .                  lihdim,ipc,nsp,pprel,sk,hk,ok)
C-    Fully Relativistic 2c Hamiltonian in gamma-rep.
C--------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ldg   :number of lower+intermediate+higher orbitals
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   sk    :structure constants, s^beta
Ci   nsp   : 2 for spin polarized case
Cu Updates
Cu   18 Jun 04 (A Chantis) first created
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nl,lihdim,indxsh(*),nsp,neul,ldim,ipc(nbas),lmr,l,
     .        ic,ms1,ms2,imu,m1,m2,i1,j1,i,j,l2,ibas
      double precision sk(ldim,2,ldim,2*2),ok(ldim,2,ldim,2*2),
     . hk(ldim,2,ldim,2*2),qsp(4),eula(nbas,neul,3),
     . pprel(4,nl,2*nl,2,2,*),crel(2,2),delrel(2,2),
     . mu,u,clebsh(2,2),ctmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . wtmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),clms(lihdim,2,2),
     . wlms(lihdim,2,2),wlmst(lihdim,2,2)
C    . gtmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),glms(lihdim,2,2),gamrel(2,2)
C      double precision sfr(ldim,2,ldim,2*2)
      call tcn('hmfr2c')

       l2 = ldim**2*4

C        do  22  i1 = 1, 2
C        do  22  j1 = 1, 4
C        do  24  j = 1, ldim
C        do  24  i = 1, ldim
C          print*, i1,j1,sk(i,i1,j,j1)
C 24   continue
C 22   continue

C       call yprm('starting Srel',2,sk,(ldim)**2,ldim,ldim,ldim)

       call dcopy(ldim**2*2,sk,1,ok,1)

       call rotspn(0,nbas,nl,indxsh,eula,neul,qsp(4),0,0,
     .   0,ldim,ldim,ldim,ldim,ok,sk)

C        do  23  i1 = 1, 2
C        do  23  j1 = 1, 4
C        do  25  j = 1, ldim
C        do  25  i = 1, ldim
C          print*, i1,j1,sk(i,i1,j,j1)
C 25    continue
C 23   continue
C      stop

c ...   Make double complex format
c       call ztoyy(tmps,2*ldim,2*ldim,2*ldim,2*ldim,0,1)
c       call zprm('ending S',2,tmps,2*ldim,2*ldim,2*ldim)

c       do i = 1, 2
c       do j = 1, 2
c       do i1 = 1, ldim
c          do j1= 1, ldim
c             print*, tmps(i1,i,j1,j)
c          enddo
c       enddo
c       enddo
c       enddo


c     call yprm('starting Srel',12,sk,(2*ldim)**2,2*ldim,2*ldim,2*ldim)


C ... Nylify these matrices, the loop over lmr never
C     determines some elements (such as xlms(1,2,1),xlms(5,2,1), e.t.c)
C    ifc compiler does not like it

       do ms1 = 1, 2
          do ms2 = 1, 2
             do i1 = 1, ldim
              clms(i1,ms1,ms2) = 0.d0
C             glms(i1,ms1,ms2) = 0.d0
              wlms(i1,ms1,ms2) = 0.d0
             enddo
          enddo
        enddo

       lmr = 0
       do  ibas = 1, nbas
        ic = ipc(ibas)
        do  l = 0, nl-1
          if (indxsh(lmr+1) .gt. lihdim) then
            lmr = lmr + 2*l+1
            goto 2
          endif

        do  imu = 1, 2*(l+1)

C   ... Extract pprel-gamma

            do ms1 = 1, 2
            do ms2 = 1, 2
              crel(ms1,ms2)   = pprel(1,l+1,imu,ms1,ms2,ic)
C             gamrel(ms1,ms2) = pprel(2,l+1,imu,ms1,ms2,ic)
              delrel(ms1,ms2) = pprel(3,l+1,imu,ms1,ms2,ic)
            enddo
            enddo

C       ... Setup for rotation of kappa-mu to lms rep'sn
            mu= dble(imu-l) - 1.5d0
            u = mu/(dble(l)+0.5d0)
            clebsh(1,1) = dsqrt(0.5d0*(1d0+u))
            clebsh(1,2) = -dsqrt(0.5d0*(1d0-u))
            clebsh(2,2) = clebsh(1,1)
            clebsh(2,1) = -clebsh(1,2)

C       --- Rotate to lms rep'sn from kappa-mu rep'sn  ---
            do ms1 = 1, 2
            do ms2 = 1, 2
              m1 = int(mu - (ms1-1.5d0))
              m2 = int(mu - (ms2-1.5d0))
              if (iabs(m1) .le. l .and. iabs(m2) .le. l) then

C       --- Test rotate C,gamma,delt to lms rep'sn from kappa-mu rep'sn  ---
                ctmp(m1,m2,ms1,ms2) =
     .            clebsh(1,ms1)*crel(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*crel(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*crel(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*crel(2,2)*clebsh(2,ms2)

C                gtmp(m1,m2,ms1,ms2) =
C     .            clebsh(1,ms1)*gamrel(1,1)*clebsh(1,ms2) +
C     .            clebsh(1,ms1)*gamrel(1,2)*clebsh(2,ms2) +
C     .            clebsh(2,ms1)*gamrel(2,1)*clebsh(1,ms2) +
C     .            clebsh(2,ms1)*gamrel(2,2)*clebsh(2,ms2)

                wtmp(m1,m2,ms1,ms2) =
     .            clebsh(1,ms1)*delrel(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*delrel(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*delrel(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*delrel(2,2)*clebsh(2,ms2)

              else
                ctmp(m1,m2,ms1,ms2) = 0d0
C               gtmp(m1,m2,ms1,ms2) = 0d0
                wtmp(m1,m2,ms1,ms2) = 0d0
              endif

            enddo
            enddo

C       ... Poke into P (ms1 is -1/2 and ms2 is +1/2)

            do  ms1 = 1, 2
            do  ms2 = 1, 2
              m1 = int(mu - (ms1-1.5d0))
              m2 = int(mu - (ms2-1.5d0))
              if (iabs(m1) .le. l .and. iabs(m2) .le. l) then
              lmr = nl*nl*(ibas-1) + (l+1)**2 - (2*l+1) + (l+m1) + 1
C............fix lmr: order (1,0,-1) instead of (-1,0,1)
              lmr = lmr - 2*m1

              clms(indxsh(lmr),ms1,ms2) = ctmp(m1,m2,ms1,ms2)
C             glms(indxsh(lmr),ms1,ms2) = gtmp(m1,m2,ms1,ms2)
              wlms(indxsh(lmr),ms1,ms2) = wtmp(m1,m2,ms1,ms2)

              endif

             enddo
            enddo

C     ... Loop over imu
          enddo

C   ... Loop over l
        enddo
    2   continue
C ... Loop over ibas
      enddo

C        call prmx('c-rel',clms,lihdim,lihdim,4)
C        call prmx('g-rel',glms,lihdim,lihdim,4)
C        call prmx('w-rel',wlms,lihdim,lihdim,4)

C ... Construct Hamiltonian
C      print *, 'START BUILDING HAMILTONIAN'
      do i1 = 1, ldim
      do i = 1, 2
         do j = 1, 2
            wlmst(i1,i,j) = 0.d0
         enddo
      enddo
      enddo

c      do i = 1, 2
c      do j = 1, 2
c         do i1 = 1, ldim
c          do j1 = 1, ldim
c             sfr(i1,i,j1,j) = 0.0d0
c          enddo
c         enddo
c      enddo
c      enddo

      do j = 1, 2
         do i1 = 1, ldim
         wlmst(i1,j,j) = wlms(i1,j,j)
         enddo
      enddo

      do i1 = 2, ldim-1
         wlmst(i1+1,2,1) = wlms(i1,1,2)
         wlmst(i1,1,2) = wlms(i1+1,2,1)
      enddo

c      do i = 1, 2
c         do i1 = 1, ldim
c            sfr(i1,i,i1,i) = wlms(i1,i,i)
c         enddo
c      enddo
c      do i1 = 1, ldim-1
c         sfr(i1,1,i1+1,2) = wlms(i1,1,2)
c      enddo
c      do i1 = 2, ldim
c         sfr(i1,2,i1-1,1) = wlms(i1,2,1)
c      enddo

c      call yprm('S',2,sfr,2*ldim,2*ldim,2*ldim)

C ... 1-1 block
      do i = 1, ldim
       do j = 1, ldim

      if(i .eq. ldim .and. j .lt. ldim) then

      hk(i,1,j,1)=wlmst(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .        wlmst(i,1,1)*sk(i,1,j+1,2)*wlms(j+1,2,1)

      hk(i+l2,1,j,1)=wlmst(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .        wlmst(i,1,1)*sk(i+l2,1,j+1,2)*wlms(j+1,2,1)


      else

      if(j .eq. ldim .and. i .lt. ldim) then
      hk(i,1,j,1)=wlmst(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .            wlmst(i,1,2)*sk(i+1,2,j,1)*wlms(j,1,1)

      hk(i+l2,1,j,1)=wlmst(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .            wlmst(i,1,2)*sk(i+l2+1,2,j,1)*wlms(j,1,1)

      else

      if(i .eq. ldim .and. j .eq. ldim) then
      hk(i,1,j,1)=wlmst(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)

      hk(i+l2,1,j,1)=wlmst(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)

      else
      hk(i,1,j,1)=wlmst(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .         wlmst(i,1,2)*sk(i+1,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,1,1)*sk(i,1,j+1,2)*wlms(j+1,2,1)+
     .         wlmst(i,1,2)*sk(i+1,2,j+1,2)*wlms(j+1,2,1)

      hk(i+l2,1,j,1)=wlmst(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .         wlmst(i,1,2)*sk(i+l2+1,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,1,1)*sk(i+l2,1,j+1,2)*wlms(j+1,2,1)+
     .         wlmst(i,1,2)*sk(i+l2+1,2,j+1,2)*wlms(j+1,2,1)

        endif
       endif
      endif

C ... 2-2 block

      if(i .eq. 1 .and. j .gt. 1) then
      hk(i,2,j,2)=wlmst(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .           wlmst(i,2,2)*sk(i,2,j-1,1)*wlms(j-1,1,2)

      hk(i+l2,2,j,2)=wlmst(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .           wlmst(i,2,2)*sk(i+l2,2,j-1,1)*wlms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .gt. 1) then
      hk(i,2,j,2)=wlmst(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .           wlmst(i,2,1)*sk(i-1,1,j,2)*wlms(j,2,2)

      hk(i+l2,2,j,2)=wlmst(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .           wlmst(i,2,1)*sk(i+l2-1,1,j,2)*wlms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. 1) then
      hk(i,2,j,2)=wlmst(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)

      hk(i+l2,2,j,2)=wlmst(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)

      else

      hk(i,2,j,2)=wlmst(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .         wlmst(i,2,1)*sk(i-1,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,2,2)*sk(i,2,j-1,1)*wlms(j-1,1,2)+
     .         wlmst(i,2,1)*sk(i-1,1,j-1,1)*wlms(j-1,1,2)

      hk(i+l2,2,j,2)=wlmst(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .         wlmst(i,2,1)*sk(i+l2-1,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,2,2)*sk(i+l2,2,j-1,1)*wlms(j-1,1,2)+
     .         wlmst(i,2,1)*sk(i+l2-1,1,j-1,1)*wlms(j-1,1,2)

        endif
       endif
      endif

C ... 1-2 block
      if(i .eq. ldim .and. j .gt. 1) then
      hk(i,1,j,2)=wlmst(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,1)*sk(i,1,j-1,1)*wlms(j-1,1,2)

      hk(i+l2,1,j,2)=wlmst(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,1)*sk(i+l2,1,j-1,1)*wlms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .lt. ldim) then
      hk(i,1,j,2)=wlmst(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,2)*sk(i+1,2,j,2)*wlms(j,2,2)

      hk(i+l2,1,j,2)=wlmst(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,2)*sk(i+l2+1,2,j,2)*wlms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. ldim) then
      hk(i,1,j,2)=wlmst(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)

      hk(i+l2,1,j,2)=wlmst(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)

      else

      hk(i,1,j,2)=wlmst(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,1)*sk(i,1,j-1,1)*wlms(j-1,1,2)+
     .         wlmst(i,1,2)*sk(i+1,2,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,2)*sk(i+1,2,j-1,1)*wlms(j-1,1,2)

      hk(i+l2,1,j,2)=wlmst(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,1)*sk(i+l2,1,j-1,1)*wlms(j-1,1,2)+
     .         wlmst(i,1,2)*sk(i+l2+1,2,j,2)*wlms(j,2,2)+
     .         wlmst(i,1,2)*sk(i+l2+1,2,j-1,1)*wlms(j-1,1,2)

        endif
       endif
      endif

C ... 2-1 block
      if(i .eq. 1 .and. j .lt. ldim) then
      hk(i,2,j,1)=wlmst(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,2)*sk(i,2,j+1,2)*wlms(j+1,2,1)

      hk(i+l2,2,j,1)=wlmst(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,2)*sk(i+l2,2,j+1,2)*wlms(j+1,2,1)

      else

      if(j .eq. ldim .and. i .gt. 1) then
      hk(i,2,j,1)=wlmst(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,1)*sk(i-1,1,j,1)*wlms(j,1,1)

      hk(i+l2,2,j,1)=wlmst(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,1)*sk(i+l2-1,1,j,1)*wlms(j,1,1)

      else

        if(j .eq. ldim .and. i .eq. 1) then
      hk(i,2,j,1)=wlmst(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)

      hk(i+l2,2,j,1)=wlmst(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)

        else

      hk(i,2,j,1)=wlmst(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,1)*sk(i-1,1,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,2)*sk(i,2,j+1,2)*wlms(j+1,2,1)+
     .         wlmst(i,2,1)*sk(i-1,1,j+1,2)*wlms(j+1,2,1)

      hk(i+l2,2,j,1)=wlmst(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,1)*sk(i+l2-1,1,j,1)*wlms(j,1,1)+
     .         wlmst(i,2,2)*sk(i+l2,2,j+1,2)*wlms(j+1,2,1)+
     .         wlmst(i,2,1)*sk(i+l2-1,1,j+1,2)*wlms(j+1,2,1)

        endif
       endif
      endif

       enddo
      enddo
C         print*, 'MADE W^T*S*W'

c       call yprm('W^T*S*W',12,hk,(2*ldim)**2,2*ldim,2*ldim,2*ldim)

C ... Add C

      do i = 1, ldim
       do i1 = 1, 2
        hk(i,i1,i,i1) = hk(i,i1,i,i1) + clms(i,i1,i1)
       enddo
      enddo

C      print*,'ADDED DIAGONAL C'
      do i = 2, ldim
       hk(i-1,1,i,2) = hk(i-1,1,i,2) + clms(i-1,1,2)
       hk(i,2,i-1,1) = hk(i,2,i-1,1) + clms(i,2,1)
      enddo

C      print*, 'MADE THE HAMILTONIAN'
C ... Make Overlap matrix for gamma-rep.

       do i = 1, 2
       do j = 1, 2
       do i1 = 1, ldim
          do j1= 1, ldim
             ok(i1,i,j1,j) = 0.d0
             ok(i1+l2,i,j1,j) = 0.d0
             if (i1 .eq. j1 .and. i .eq. j) then
             ok(i1,i,j1,j) = 1.d0
             endif
C             print*, dcmplx(ok(i1,i,j1,j),ok(i1+l2,i,j1,j))
C             print*, dcmplx(hk(i1,i,j1,j),hk(i1+l2,i,j1,j))
          enddo
       enddo
       enddo
       enddo

c       do i = 1, 2
c       do j = 1, 4
c       do i1 = 1, ldim
c          do j1= 1, ldim
c             print*, i1,i,j1,j,hk(i1,i,j1,j)
c          enddo
c       enddo
c       enddo
c       enddo
c       stop

C      call yprm('hk relativ.',12,hk,(2*ldim)**2,2*ldim,2*ldim,2*ldim)
C      call yprm('ok relativ.',12,ok,(2*ldim)**2,2*ldim,2*ldim,2*ldim)

       call tcx('hmfr2c')

       end











