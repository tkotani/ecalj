      subroutine hmfr3c(nbas,nl,indxsh,qsp,eula,neul,ldim,
     .                  lihdim,ipc,nsp,pprel,pph,ksop,sk,hk,ok)
C- Fully Relativistic 3c Hamiltonian in gamma-rep.
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
     . pph(6,nl,nsp,*),
     . wtmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . clms(lihdim,2,2),wlms(lihdim,2,2),
     . wlmst(lihdim,2,2),
     . dlt(2,2),dlrelt(2,2),cme(2,2),penu(2,2),
     . hk1(ldim,2,ldim,2*2),hk2(ldim,2,ldim,2*2),hk3(ldim,2,ldim,2*2),
     . hk4(ldim,2,ldim,2*2),mk1(ldim,2,ldim,2),mk2(ldim,2,ldim,2),
     . mk3(ldim,2,ldim,2),mk4(ldim,2,ldim,2),
     . mx0(2,2),mx1(2,2),mx2(2,2),mx3(2,2),tmp(2,2),tmp1(2,2),
     . mx0tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . mx1tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . mx2tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . mx3tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . mx0lms(lihdim,2,2),mx1lms(lihdim,2,2),mx2lms(lihdim,2,2),
     . mx3lms(lihdim,2,2),ksop(0:nl-1,nsp,nsp,9,*),pg(2,2),
     . enum(2,2),enu
C     double precision glms(lihdim,2,2)
C     double precision gamrel(2,2),gtmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2)
C     double precision sfr(ldim,2,ldim,2*2)
C Overlap third order related matrices
      integer k
      double precision ox0(2,2),ox1(2,2),ox2(2,2),ox3(2,2),
     . ox0tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . ox1tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . ox2tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . ox3tmp(-(nl-1):(nl-1),-(nl-1):(nl-1),2,2),
     . ox0lms(lihdim,2,2),ox1lms(lihdim,2,2),ox2lms(lihdim,2,2),
     . ox3lms(lihdim,2,2),ok1(ldim,2,ldim,2*2),ok2(ldim,2,ldim,2*2),
     . ok3(ldim,2,ldim,2*2),ok4(ldim,2,ldim,2*2),
     . omk1(ldim,2,ldim,2),omk2(ldim,2,ldim,2),
     . omk3(ldim,2,ldim,2),omk4(ldim,2,ldim,2)

      call tcn('hmfr3c')

       l2 = ldim**2*4

C       do i = 1, 2
C       do j = 1, 2
C       do i1 = 1, ldim
C          do j1= 1, ldim
C             ok(i1,i,j1,j) = 0.d0
C             ok(i1+l2,i,j1,j) = 0.d0
C             if (i1 .eq. j1 .and. i .eq. j) then
C             ok(i1,i,j1,j) = 1.d0
C             endif
C             print*, dcmplx(ok(i1,i,j1,j),ok(i1+l2,i,j1,j))
C             print*, dcmplx(hk(i1,i,j1,j),hk(i1+l2,i,j1,j))
C          enddo
C       enddo
C       enddo
C       enddo

C  ... Kronecker's Symbol
       do ms1 = 1, 2
          do ms2 = 1, 2
             dlt(ms1,ms2) =0d0
             if (ms1 .eq. ms2) dlt(ms1,ms2) = 1d0
          enddo
       enddo

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
              mx0lms(i1,ms1,ms2) = 0.d0
              mx1lms(i1,ms1,ms2) = 0.d0
              mx2lms(i1,ms1,ms2) = 0.d0
              mx3lms(i1,ms1,ms2) = 0.d0
              ox0lms(i1,ms1,ms2) = 0.d0
              ox1lms(i1,ms1,ms2) = 0.d0
              ox2lms(i1,ms1,ms2) = 0.d0
              ox3lms(i1,ms1,ms2) = 0.d0
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

           enu = 0.5d0*(pph(1,l+1,1,ic)+pph(1,l+1,2,ic))

            do ms1 = 1, 2
            do ms2 = 1, 2
              crel(ms1,ms2)   = pprel(1,l+1,imu,ms1,ms2,ic)
C             gamrel(ms1,ms2) = pprel(2,l+1,imu,ms1,ms2,ic)
              delrel(ms1,ms2) = pprel(3,l+1,imu,ms1,ms2,ic)
              dlrelt(ms1,ms2) = pprel(3,l+1,imu,ms2,ms1,ic)
C   ...Construct p*Enu matrix, (C-Enu) matrix
              pg(ms1,ms2) = pprel(4,l+1,imu,ms1,ms2,ic)
              enum(ms1,ms2) = 0d0
C              cme(ms1,ms2) = crel(ms1,ms2) -
C     .        0.5d0*(pph(1,l+1,1,ic)+pph(1,l+1,2,ic))*dlt(ms1,ms2)
            enddo
            enddo


C            enum(1,1) = -0.5d0*dble((l+1))*ksop(l,1,1,1,ic) + enu
C            enum(2,2) = 0.5d0*dble(l)*ksop(l,1,1,1,ic) + enu

            if (imu .eq. 1) then
             enum(1,1) = -0.5d0*dble((l+1))*ksop(l,1,1,1,ic)
     .                    + pph(1,l+1,1,ic)
             enum(2,2) = 0.5d0*dble(l)*ksop(l,1,1,1,ic) +
     .                    pph(1,l+1,1,ic)
            else
             if (imu .eq. 2*(l+1)) then
             enum(1,1) = -0.5d0*dble((l+1))*ksop(l,2,2,1,ic)
     .                    + pph(1,l+1,2,ic)
             enum(2,2) = 0.5d0*dble(l)*ksop(l,2,2,1,ic) +
     .                    pph(1,l+1,2,ic)
             else
            enum(1,1) = -0.5d0*dble((l+1))*(ksop(l,1,1,1,ic) +
     .      ksop(l,2,2,1,ic))*0.5d0 + enu
            enum(2,2) = 0.5d0*dble(l)*(ksop(l,1,1,1,ic) +
     .      ksop(l,2,2,1,ic))*0.5d0 + enu
             endif
             endif

            do ms1 = 1, 2
               do ms2 = 1, 2
                  cme(ms1,ms2) = crel(ms1,ms2) - enum(ms1,ms2)
               enddo
            enddo

C   ... Construct mx0 =(C-Enu)*Enu*p*(C-Enu), tmp=(C-Enu)*Enu*p

            call mul22(enum,pg,penu)

            call mul22(cme,penu,tmp)
            call mul22(tmp,cme,mx0)

C   ... Construct mx1 = tmp*W^(T), mx2 = W^(T)*Enu*p*(C-Enu)

            call mul22(tmp,dlrelt,mx1)
            call mul22(dlrelt,penu,tmp1)
            call mul22(tmp1,cme,mx2)

C   ... Contrcut mx3  = W^(T)*Enu*p

            call mul22(dlrelt,penu,mx3)

C   ... Construct  the matrices needed for overlap matrix

C   ... Construct ox0 =(C-Enu)*p*(C-Enu), tmp=(C-Enu)*p

            call mul22(cme,pg,tmp)
            call mul22(tmp,cme,ox0)

C   ... Construct ox1 = tmp*W^(T), ox2 = W^(T)*p*(C-Enu)

            call mul22(tmp,dlrelt,ox1)
            call mul22(dlrelt,pg,tmp1)
            call mul22(tmp1,cme,ox2)

C   ... Contrcut ox3  = W^(T)*p

            call mul22(dlrelt,pg,ox3)
c----------------------------------------------------

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

C       --- Rotate C,gamma,delt to lms rep'sn from kappa-mu rep'sn  ---
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
C      --- Rotate mx(0 to 4)
                mx0tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*mx0(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*mx0(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*mx0(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*mx0(2,2)*clebsh(2,ms2)

C .... Test
C                mx0tmp(m1,m2,ms1,ms2) =
C     .           clebsh(1,ms1)*pg(1,1)*clebsh(1,ms2) +
C     .            clebsh(1,ms1)*pg(1,2)*clebsh(2,ms2) +
C     .            clebsh(2,ms1)*pg(2,1)*clebsh(1,ms2) +
C     .            clebsh(2,ms1)*pg(2,2)*clebsh(2,ms2)
C ....

                mx1tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*mx1(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*mx1(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*mx1(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*mx1(2,2)*clebsh(2,ms2)

                mx2tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*mx2(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*mx2(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*mx2(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*mx2(2,2)*clebsh(2,ms2)

                mx3tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*mx3(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*mx3(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*mx3(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*mx3(2,2)*clebsh(2,ms2)

C   ... Rotate the matrices from which the ovrelap is costructed

                ox0tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*ox0(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*ox0(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*ox0(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*ox0(2,2)*clebsh(2,ms2)

                ox1tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*ox1(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*ox1(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*ox1(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*ox1(2,2)*clebsh(2,ms2)

                ox2tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*ox2(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*ox2(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*ox2(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*ox2(2,2)*clebsh(2,ms2)

                ox3tmp(m1,m2,ms1,ms2) =
     .           clebsh(1,ms1)*ox3(1,1)*clebsh(1,ms2) +
     .            clebsh(1,ms1)*ox3(1,2)*clebsh(2,ms2) +
     .            clebsh(2,ms1)*ox3(2,1)*clebsh(1,ms2) +
     .            clebsh(2,ms1)*ox3(2,2)*clebsh(2,ms2)

              else
                ctmp(m1,m2,ms1,ms2) = 0d0
C               gtmp(m1,m2,ms1,ms2) = 0d0
                wtmp(m1,m2,ms1,ms2) = 0d0
                mx0tmp(m1,m2,ms1,ms2) = 0d0
                mx1tmp(m1,m2,ms1,ms2) = 0d0
                mx2tmp(m1,m2,ms1,ms2) = 0d0
                mx3tmp(m1,m2,ms1,ms2) = 0d0
                ox0tmp(m1,m2,ms1,ms2) = 0d0
                ox1tmp(m1,m2,ms1,ms2) = 0d0
                ox2tmp(m1,m2,ms1,ms2) = 0d0
                ox3tmp(m1,m2,ms1,ms2) = 0d0
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
              mx0lms(indxsh(lmr),ms1,ms2) = mx0tmp(m1,m2,ms1,ms2)
              mx1lms(indxsh(lmr),ms1,ms2) = mx1tmp(m1,m2,ms1,ms2)
              mx2lms(indxsh(lmr),ms1,ms2) = mx2tmp(m1,m2,ms1,ms2)
              mx3lms(indxsh(lmr),ms1,ms2) = mx3tmp(m1,m2,ms1,ms2)
              ox0lms(indxsh(lmr),ms1,ms2) = ox0tmp(m1,m2,ms1,ms2)
              ox1lms(indxsh(lmr),ms1,ms2) = ox1tmp(m1,m2,ms1,ms2)
              ox2lms(indxsh(lmr),ms1,ms2) = ox2tmp(m1,m2,ms1,ms2)
              ox3lms(indxsh(lmr),ms1,ms2) = ox3tmp(m1,m2,ms1,ms2)
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
C         call prmx('pg-rel',mx0lms,lihdim,lihdim,4)

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

C  ... Construct the thrird order terms

C  ... mx1*S*W

C ... 1-1 block
      do i = 1, ldim
       do j = 1, ldim

      if(i .eq. ldim .and. j .lt. ldim) then

      hk1(i,1,j,1)=mx1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .        mx1lms(i,1,1)*sk(i,1,j+1,2)*wlms(j+1,2,1)

      hk1(i+l2,1,j,1)=mx1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .        mx1lms(i,1,1)*sk(i+l2,1,j+1,2)*wlms(j+1,2,1)


      else

      if(j .eq. ldim .and. i .lt. ldim) then
      hk1(i,1,j,1)=mx1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .            mx1lms(i,1,2)*sk(i+1,2,j,1)*wlms(j,1,1)

      hk1(i+l2,1,j,1)=mx1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .            mx1lms(i,1,2)*sk(i+l2+1,2,j,1)*wlms(j,1,1)

      else

      if(i .eq. ldim .and. j .eq. ldim) then
      hk1(i,1,j,1)=mx1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)

      hk1(i+l2,1,j,1)=mx1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)

      else
      hk1(i,1,j,1)=mx1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .         mx1lms(i,1,2)*sk(i+1,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,1,1)*sk(i,1,j+1,2)*wlms(j+1,2,1)+
     .         mx1lms(i,1,2)*sk(i+1,2,j+1,2)*wlms(j+1,2,1)

      hk1(i+l2,1,j,1)=mx1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .         mx1lms(i,1,2)*sk(i+l2+1,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,1,1)*sk(i+l2,1,j+1,2)*wlms(j+1,2,1)+
     .         mx1lms(i,1,2)*sk(i+l2+1,2,j+1,2)*wlms(j+1,2,1)

        endif
       endif
      endif

C ... 2-2 block

      if(i .eq. 1 .and. j .gt. 1) then
      hk1(i,2,j,2)=mx1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .           mx1lms(i,2,2)*sk(i,2,j-1,1)*wlms(j-1,1,2)

      hk1(i+l2,2,j,2)=mx1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .           mx1lms(i,2,2)*sk(i+l2,2,j-1,1)*wlms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .gt. 1) then
      hk1(i,2,j,2)=mx1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .           mx1lms(i,2,1)*sk(i-1,1,j,2)*wlms(j,2,2)

      hk1(i+l2,2,j,2)=mx1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .           mx1lms(i,2,1)*sk(i+l2-1,1,j,2)*wlms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. 1) then
      hk1(i,2,j,2)=mx1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)

      hk1(i+l2,2,j,2)=mx1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)

      else

      hk1(i,2,j,2)=mx1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .         mx1lms(i,2,1)*sk(i-1,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,2,2)*sk(i,2,j-1,1)*wlms(j-1,1,2)+
     .         mx1lms(i,2,1)*sk(i-1,1,j-1,1)*wlms(j-1,1,2)

      hk1(i+l2,2,j,2)=mx1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .         mx1lms(i,2,1)*sk(i+l2-1,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,2,2)*sk(i+l2,2,j-1,1)*wlms(j-1,1,2)+
     .         mx1lms(i,2,1)*sk(i+l2-1,1,j-1,1)*wlms(j-1,1,2)

        endif
       endif
      endif

C ... 1-2 block
      if(i .eq. ldim .and. j .gt. 1) then
      hk1(i,1,j,2)=mx1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,1)*sk(i,1,j-1,1)*wlms(j-1,1,2)

      hk1(i+l2,1,j,2)=mx1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,1)*sk(i+l2,1,j-1,1)*wlms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .lt. ldim) then
      hk1(i,1,j,2)=mx1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,2)*sk(i+1,2,j,2)*wlms(j,2,2)

      hk1(i+l2,1,j,2)=mx1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,2)*sk(i+l2+1,2,j,2)*wlms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. ldim) then
      hk1(i,1,j,2)=mx1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)

      hk1(i+l2,1,j,2)=mx1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)

      else

      hk1(i,1,j,2)=mx1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,1)*sk(i,1,j-1,1)*wlms(j-1,1,2)+
     .         mx1lms(i,1,2)*sk(i+1,2,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,2)*sk(i+1,2,j-1,1)*wlms(j-1,1,2)

      hk1(i+l2,1,j,2)=mx1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,1)*sk(i+l2,1,j-1,1)*wlms(j-1,1,2)+
     .         mx1lms(i,1,2)*sk(i+l2+1,2,j,2)*wlms(j,2,2)+
     .         mx1lms(i,1,2)*sk(i+l2+1,2,j-1,1)*wlms(j-1,1,2)

        endif
       endif
      endif

C ... 2-1 block
      if(i .eq. 1 .and. j .lt. ldim) then
      hk1(i,2,j,1)=mx1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,2)*sk(i,2,j+1,2)*wlms(j+1,2,1)

      hk1(i+l2,2,j,1)=mx1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,2)*sk(i+l2,2,j+1,2)*wlms(j+1,2,1)

      else

      if(j .eq. ldim .and. i .gt. 1) then
      hk1(i,2,j,1)=mx1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,1)*sk(i-1,1,j,1)*wlms(j,1,1)

      hk1(i+l2,2,j,1)=mx1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,1)*sk(i+l2-1,1,j,1)*wlms(j,1,1)

      else

        if(j .eq. ldim .and. i .eq. 1) then
      hk1(i,2,j,1)=mx1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)

      hk1(i+l2,2,j,1)=mx1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)

        else

      hk1(i,2,j,1)=mx1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,1)*sk(i-1,1,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,2)*sk(i,2,j+1,2)*wlms(j+1,2,1)+
     .         mx1lms(i,2,1)*sk(i-1,1,j+1,2)*wlms(j+1,2,1)

      hk1(i+l2,2,j,1)=mx1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,1)*sk(i+l2-1,1,j,1)*wlms(j,1,1)+
     .         mx1lms(i,2,2)*sk(i+l2,2,j+1,2)*wlms(j+1,2,1)+
     .         mx1lms(i,2,1)*sk(i+l2-1,1,j+1,2)*wlms(j+1,2,1)

        endif
       endif
      endif

       enddo
      enddo

C ... W*S*mx2

C ... 1-1 block
      do i = 1, ldim
       do j = 1, ldim

      if(i .eq. ldim .and. j .lt. ldim) then

      hk2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx2lms(j,1,1)+
     .        wlms(i,1,1)*sk(i,1,j+1,2)*mx2lms(j+1,2,1)

      hk2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx2lms(j,1,1)+
     .        wlms(i,1,1)*sk(i+l2,1,j+1,2)*mx2lms(j+1,2,1)


      else

      if(j .eq. ldim .and. i .lt. ldim) then
      hk2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx2lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+1,2,j,1)*mx2lms(j,1,1)

      hk2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx2lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+l2+1,2,j,1)*mx2lms(j,1,1)

      else

      if(i .eq. ldim .and. j .eq. ldim) then
      hk2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx2lms(j,1,1)

      hk2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx2lms(j,1,1)

      else
      hk2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx2lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+1,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,1,1)*sk(i,1,j+1,2)*mx2lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+1,2,j+1,2)*mx2lms(j+1,2,1)

      hk2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx2lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,1,1)*sk(i+l2,1,j+1,2)*mx2lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j+1,2)*mx2lms(j+1,2,1)

        endif
       endif
      endif

C ... 2-2 block

      if(i .eq. 1 .and. j .gt. 1) then
      hk2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx2lms(j,2,2)+
     .           wlms(i,2,2)*sk(i,2,j-1,1)*mx2lms(j-1,1,2)

      hk2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx2lms(j,2,2)+
     .           wlms(i,2,2)*sk(i+l2,2,j-1,1)*mx2lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .gt. 1) then
      hk2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx2lms(j,2,2)+
     .           wlms(i,2,1)*sk(i-1,1,j,2)*mx2lms(j,2,2)

      hk2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx2lms(j,2,2)+
     .           wlms(i,2,1)*sk(i+l2-1,1,j,2)*mx2lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. 1) then
      hk2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx2lms(j,2,2)

      hk2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx2lms(j,2,2)

      else

      hk2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx2lms(j,2,2)+
     .         wlms(i,2,1)*sk(i-1,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,2,2)*sk(i,2,j-1,1)*mx2lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i-1,1,j-1,1)*mx2lms(j-1,1,2)

      hk2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx2lms(j,2,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,2,2)*sk(i+l2,2,j-1,1)*mx2lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j-1,1)*mx2lms(j-1,1,2)

        endif
       endif
      endif

C ... 1-2 block
      if(i .eq. ldim .and. j .gt. 1) then
      hk2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*mx2lms(j-1,1,2)

      hk2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*mx2lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .lt. ldim) then
      hk2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*mx2lms(j,2,2)

      hk2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*mx2lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. ldim) then
      hk2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx2lms(j,2,2)

      hk2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx2lms(j,2,2)

      else

      hk2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*mx2lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j-1,1)*mx2lms(j-1,1,2)

      hk2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*mx2lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*mx2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j-1,1)*mx2lms(j-1,1,2)

        endif
       endif
      endif

C ... 2-1 block
      if(i .eq. 1 .and. j .lt. ldim) then
      hk2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*mx2lms(j+1,2,1)

      hk2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*mx2lms(j+1,2,1)

      else

      if(j .eq. ldim .and. i .gt. 1) then
      hk2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*mx2lms(j,1,1)

      hk2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*mx2lms(j,1,1)

      else

        if(j .eq. ldim .and. i .eq. 1) then
      hk2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx2lms(j,1,1)

      hk2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx2lms(j,1,1)

        else

      hk2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*mx2lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i-1,1,j+1,2)*mx2lms(j+1,2,1)

      hk2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*mx2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*mx2lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j+1,2)*mx2lms(j+1,2,1)

        endif
       endif
      endif

       enddo
      enddo

C ... W*S*mx3

C ... 1-1 block
      do i = 1, ldim
       do j = 1, ldim

      if(i .eq. ldim .and. j .lt. ldim) then

      hk3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx3lms(j,1,1)+
     .        wlms(i,1,1)*sk(i,1,j+1,2)*mx3lms(j+1,2,1)

      hk3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx3lms(j,1,1)+
     .        wlms(i,1,1)*sk(i+l2,1,j+1,2)*mx3lms(j+1,2,1)


      else

      if(j .eq. ldim .and. i .lt. ldim) then
      hk3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx3lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+1,2,j,1)*mx3lms(j,1,1)

      hk3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx3lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+l2+1,2,j,1)*mx3lms(j,1,1)

      else

      if(i .eq. ldim .and. j .eq. ldim) then
      hk3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx3lms(j,1,1)

      hk3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx3lms(j,1,1)

      else
      hk3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*mx3lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+1,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,1,1)*sk(i,1,j+1,2)*mx3lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+1,2,j+1,2)*mx3lms(j+1,2,1)

      hk3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*mx3lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,1,1)*sk(i+l2,1,j+1,2)*mx3lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j+1,2)*mx3lms(j+1,2,1)

        endif
       endif
      endif

C ... 2-2 block

      if(i .eq. 1 .and. j .gt. 1) then
      hk3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx3lms(j,2,2)+
     .           wlms(i,2,2)*sk(i,2,j-1,1)*mx3lms(j-1,1,2)

      hk3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx3lms(j,2,2)+
     .           wlms(i,2,2)*sk(i+l2,2,j-1,1)*mx3lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .gt. 1) then
      hk3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx3lms(j,2,2)+
     .           wlms(i,2,1)*sk(i-1,1,j,2)*mx3lms(j,2,2)

      hk3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx3lms(j,2,2)+
     .           wlms(i,2,1)*sk(i+l2-1,1,j,2)*mx3lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. 1) then
      hk3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx3lms(j,2,2)

      hk3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx3lms(j,2,2)

      else

      hk3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*mx3lms(j,2,2)+
     .         wlms(i,2,1)*sk(i-1,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,2,2)*sk(i,2,j-1,1)*mx3lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i-1,1,j-1,1)*mx3lms(j-1,1,2)

      hk3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*mx3lms(j,2,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,2,2)*sk(i+l2,2,j-1,1)*mx3lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j-1,1)*mx3lms(j-1,1,2)

        endif
       endif
      endif

C ... 1-2 block
      if(i .eq. ldim .and. j .gt. 1) then
      hk3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*mx3lms(j-1,1,2)

      hk3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*mx3lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .lt. ldim) then
      hk3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*mx3lms(j,2,2)

      hk3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*mx3lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. ldim) then
      hk3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx3lms(j,2,2)

      hk3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx3lms(j,2,2)

      else

      hk3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*mx3lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j-1,1)*mx3lms(j-1,1,2)

      hk3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*mx3lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*mx3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j-1,1)*mx3lms(j-1,1,2)

        endif
       endif
      endif

C ... 2-1 block
      if(i .eq. 1 .and. j .lt. ldim) then
      hk3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*mx3lms(j+1,2,1)

      hk3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*mx3lms(j+1,2,1)

      else

      if(j .eq. ldim .and. i .gt. 1) then
      hk3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*mx3lms(j,1,1)

      hk3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*mx3lms(j,1,1)

      else

        if(j .eq. ldim .and. i .eq. 1) then
      hk3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx3lms(j,1,1)

      hk3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx3lms(j,1,1)

        else

      hk3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*mx3lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i-1,1,j+1,2)*mx3lms(j+1,2,1)

      hk3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*mx3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*mx3lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j+1,2)*mx3lms(j+1,2,1)

        endif
       endif
      endif

       enddo
      enddo

C  ... Make the Overlap matrix terms

C  ... ox1*S*W

C ... 1-1 block
      do i = 1, ldim
       do j = 1, ldim

      if(i .eq. ldim .and. j .lt. ldim) then

      ok1(i,1,j,1)=ox1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .        ox1lms(i,1,1)*sk(i,1,j+1,2)*wlms(j+1,2,1)

      ok1(i+l2,1,j,1)=ox1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .        ox1lms(i,1,1)*sk(i+l2,1,j+1,2)*wlms(j+1,2,1)


      else

      if(j .eq. ldim .and. i .lt. ldim) then
      ok1(i,1,j,1)=ox1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .            ox1lms(i,1,2)*sk(i+1,2,j,1)*wlms(j,1,1)

      ok1(i+l2,1,j,1)=ox1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .            ox1lms(i,1,2)*sk(i+l2+1,2,j,1)*wlms(j,1,1)

      else

      if(i .eq. ldim .and. j .eq. ldim) then
      ok1(i,1,j,1)=ox1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)

      ok1(i+l2,1,j,1)=ox1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)

      else
      ok1(i,1,j,1)=ox1lms(i,1,1)*sk(i,1,j,1)*wlms(j,1,1)+
     .         ox1lms(i,1,2)*sk(i+1,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,1,1)*sk(i,1,j+1,2)*wlms(j+1,2,1)+
     .         ox1lms(i,1,2)*sk(i+1,2,j+1,2)*wlms(j+1,2,1)

      ok1(i+l2,1,j,1)=ox1lms(i,1,1)*sk(i+l2,1,j,1)*wlms(j,1,1)+
     .         ox1lms(i,1,2)*sk(i+l2+1,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,1,1)*sk(i+l2,1,j+1,2)*wlms(j+1,2,1)+
     .         ox1lms(i,1,2)*sk(i+l2+1,2,j+1,2)*wlms(j+1,2,1)

        endif
       endif
      endif

C ... 2-2 block

      if(i .eq. 1 .and. j .gt. 1) then
      ok1(i,2,j,2)=ox1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .           ox1lms(i,2,2)*sk(i,2,j-1,1)*wlms(j-1,1,2)

      ok1(i+l2,2,j,2)=ox1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .           ox1lms(i,2,2)*sk(i+l2,2,j-1,1)*wlms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .gt. 1) then
      ok1(i,2,j,2)=ox1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .           ox1lms(i,2,1)*sk(i-1,1,j,2)*wlms(j,2,2)

      ok1(i+l2,2,j,2)=ox1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .           ox1lms(i,2,1)*sk(i+l2-1,1,j,2)*wlms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. 1) then
      ok1(i,2,j,2)=ox1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)

      ok1(i+l2,2,j,2)=ox1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)

      else

      ok1(i,2,j,2)=ox1lms(i,2,2)*sk(i,2,j,2)*wlms(j,2,2)+
     .         ox1lms(i,2,1)*sk(i-1,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,2,2)*sk(i,2,j-1,1)*wlms(j-1,1,2)+
     .         ox1lms(i,2,1)*sk(i-1,1,j-1,1)*wlms(j-1,1,2)

      ok1(i+l2,2,j,2)=ox1lms(i,2,2)*sk(i+l2,2,j,2)*wlms(j,2,2)+
     .         ox1lms(i,2,1)*sk(i+l2-1,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,2,2)*sk(i+l2,2,j-1,1)*wlms(j-1,1,2)+
     .         ox1lms(i,2,1)*sk(i+l2-1,1,j-1,1)*wlms(j-1,1,2)

        endif
       endif
      endif

C ... 1-2 block
      if(i .eq. ldim .and. j .gt. 1) then
      ok1(i,1,j,2)=ox1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,1)*sk(i,1,j-1,1)*wlms(j-1,1,2)

      ok1(i+l2,1,j,2)=ox1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,1)*sk(i+l2,1,j-1,1)*wlms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .lt. ldim) then
      ok1(i,1,j,2)=ox1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,2)*sk(i+1,2,j,2)*wlms(j,2,2)

      ok1(i+l2,1,j,2)=ox1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,2)*sk(i+l2+1,2,j,2)*wlms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. ldim) then
      ok1(i,1,j,2)=ox1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)

      ok1(i+l2,1,j,2)=ox1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)

      else

      ok1(i,1,j,2)=ox1lms(i,1,1)*sk(i,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,1)*sk(i,1,j-1,1)*wlms(j-1,1,2)+
     .         ox1lms(i,1,2)*sk(i+1,2,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,2)*sk(i+1,2,j-1,1)*wlms(j-1,1,2)

      ok1(i+l2,1,j,2)=ox1lms(i,1,1)*sk(i+l2,1,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,1)*sk(i+l2,1,j-1,1)*wlms(j-1,1,2)+
     .         ox1lms(i,1,2)*sk(i+l2+1,2,j,2)*wlms(j,2,2)+
     .         ox1lms(i,1,2)*sk(i+l2+1,2,j-1,1)*wlms(j-1,1,2)

        endif
       endif
      endif

C ... 2-1 block
      if(i .eq. 1 .and. j .lt. ldim) then
      ok1(i,2,j,1)=ox1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,2)*sk(i,2,j+1,2)*wlms(j+1,2,1)

      ok1(i+l2,2,j,1)=ox1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,2)*sk(i+l2,2,j+1,2)*wlms(j+1,2,1)

      else

      if(j .eq. ldim .and. i .gt. 1) then
      ok1(i,2,j,1)=ox1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,1)*sk(i-1,1,j,1)*wlms(j,1,1)

      ok1(i+l2,2,j,1)=ox1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,1)*sk(i+l2-1,1,j,1)*wlms(j,1,1)

      else

        if(j .eq. ldim .and. i .eq. 1) then
      ok1(i,2,j,1)=ox1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)

      ok1(i+l2,2,j,1)=ox1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)

        else

      ok1(i,2,j,1)=ox1lms(i,2,2)*sk(i,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,1)*sk(i-1,1,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,2)*sk(i,2,j+1,2)*wlms(j+1,2,1)+
     .         ox1lms(i,2,1)*sk(i-1,1,j+1,2)*wlms(j+1,2,1)

      ok1(i+l2,2,j,1)=ox1lms(i,2,2)*sk(i+l2,2,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,1)*sk(i+l2-1,1,j,1)*wlms(j,1,1)+
     .         ox1lms(i,2,2)*sk(i+l2,2,j+1,2)*wlms(j+1,2,1)+
     .         ox1lms(i,2,1)*sk(i+l2-1,1,j+1,2)*wlms(j+1,2,1)

        endif
       endif
      endif

       enddo
      enddo

C ... W*S*ox2

C ... 1-1 block
      do i = 1, ldim
       do j = 1, ldim

      if(i .eq. ldim .and. j .lt. ldim) then

      ok2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox2lms(j,1,1)+
     .        wlms(i,1,1)*sk(i,1,j+1,2)*ox2lms(j+1,2,1)

      ok2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox2lms(j,1,1)+
     .        wlms(i,1,1)*sk(i+l2,1,j+1,2)*ox2lms(j+1,2,1)


      else

      if(j .eq. ldim .and. i .lt. ldim) then
      ok2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox2lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+1,2,j,1)*ox2lms(j,1,1)

      ok2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox2lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+l2+1,2,j,1)*ox2lms(j,1,1)

      else

      if(i .eq. ldim .and. j .eq. ldim) then
      ok2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox2lms(j,1,1)

      ok2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox2lms(j,1,1)

      else
      ok2(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox2lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+1,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,1,1)*sk(i,1,j+1,2)*ox2lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+1,2,j+1,2)*ox2lms(j+1,2,1)

      ok2(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox2lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,1,1)*sk(i+l2,1,j+1,2)*ox2lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j+1,2)*ox2lms(j+1,2,1)

        endif
       endif
      endif

C ... 2-2 block

      if(i .eq. 1 .and. j .gt. 1) then
      ok2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox2lms(j,2,2)+
     .           wlms(i,2,2)*sk(i,2,j-1,1)*ox2lms(j-1,1,2)

      ok2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox2lms(j,2,2)+
     .           wlms(i,2,2)*sk(i+l2,2,j-1,1)*ox2lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .gt. 1) then
      ok2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox2lms(j,2,2)+
     .           wlms(i,2,1)*sk(i-1,1,j,2)*ox2lms(j,2,2)

      ok2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox2lms(j,2,2)+
     .           wlms(i,2,1)*sk(i+l2-1,1,j,2)*ox2lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. 1) then
      ok2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox2lms(j,2,2)

      ok2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox2lms(j,2,2)

      else

      ok2(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox2lms(j,2,2)+
     .         wlms(i,2,1)*sk(i-1,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,2,2)*sk(i,2,j-1,1)*ox2lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i-1,1,j-1,1)*ox2lms(j-1,1,2)

      ok2(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox2lms(j,2,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,2,2)*sk(i+l2,2,j-1,1)*ox2lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j-1,1)*ox2lms(j-1,1,2)

        endif
       endif
      endif

C ... 1-2 block
      if(i .eq. ldim .and. j .gt. 1) then
      ok2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*ox2lms(j-1,1,2)

      ok2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*ox2lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .lt. ldim) then
      ok2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*ox2lms(j,2,2)

      ok2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*ox2lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. ldim) then
      ok2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox2lms(j,2,2)

      ok2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox2lms(j,2,2)

      else

      ok2(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*ox2lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j-1,1)*ox2lms(j-1,1,2)

      ok2(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*ox2lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*ox2lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j-1,1)*ox2lms(j-1,1,2)

        endif
       endif
      endif

C ... 2-1 block
      if(i .eq. 1 .and. j .lt. ldim) then
      ok2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*ox2lms(j+1,2,1)

      ok2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*ox2lms(j+1,2,1)

      else

      if(j .eq. ldim .and. i .gt. 1) then
      ok2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*ox2lms(j,1,1)

      ok2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*ox2lms(j,1,1)

      else

        if(j .eq. ldim .and. i .eq. 1) then
      ok2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox2lms(j,1,1)

      ok2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox2lms(j,1,1)

        else

      ok2(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*ox2lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i-1,1,j+1,2)*ox2lms(j+1,2,1)

      ok2(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*ox2lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*ox2lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j+1,2)*ox2lms(j+1,2,1)

        endif
       endif
      endif

       enddo
      enddo

C ... W*S*ox3

C ... 1-1 block
      do i = 1, ldim
       do j = 1, ldim

      if(i .eq. ldim .and. j .lt. ldim) then

      ok3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox3lms(j,1,1)+
     .        wlms(i,1,1)*sk(i,1,j+1,2)*ox3lms(j+1,2,1)

      ok3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox3lms(j,1,1)+
     .        wlms(i,1,1)*sk(i+l2,1,j+1,2)*ox3lms(j+1,2,1)


      else

      if(j .eq. ldim .and. i .lt. ldim) then
      ok3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox3lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+1,2,j,1)*ox3lms(j,1,1)

      ok3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox3lms(j,1,1)+
     .            wlms(i,1,2)*sk(i+l2+1,2,j,1)*ox3lms(j,1,1)

      else

      if(i .eq. ldim .and. j .eq. ldim) then
      ok3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox3lms(j,1,1)

      ok3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox3lms(j,1,1)

      else
      ok3(i,1,j,1)=wlms(i,1,1)*sk(i,1,j,1)*ox3lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+1,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,1,1)*sk(i,1,j+1,2)*ox3lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+1,2,j+1,2)*ox3lms(j+1,2,1)

      ok3(i+l2,1,j,1)=wlms(i,1,1)*sk(i+l2,1,j,1)*ox3lms(j,1,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,1,1)*sk(i+l2,1,j+1,2)*ox3lms(j+1,2,1)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j+1,2)*ox3lms(j+1,2,1)

        endif
       endif
      endif

C ... 2-2 block

      if(i .eq. 1 .and. j .gt. 1) then
      ok3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox3lms(j,2,2)+
     .           wlms(i,2,2)*sk(i,2,j-1,1)*ox3lms(j-1,1,2)

      ok3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox3lms(j,2,2)+
     .           wlms(i,2,2)*sk(i+l2,2,j-1,1)*ox3lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .gt. 1) then
      ok3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox3lms(j,2,2)+
     .           wlms(i,2,1)*sk(i-1,1,j,2)*ox3lms(j,2,2)

      ok3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox3lms(j,2,2)+
     .           wlms(i,2,1)*sk(i+l2-1,1,j,2)*ox3lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. 1) then
      ok3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox3lms(j,2,2)

      ok3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox3lms(j,2,2)

      else

      ok3(i,2,j,2)=wlms(i,2,2)*sk(i,2,j,2)*ox3lms(j,2,2)+
     .         wlms(i,2,1)*sk(i-1,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,2,2)*sk(i,2,j-1,1)*ox3lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i-1,1,j-1,1)*ox3lms(j-1,1,2)

      ok3(i+l2,2,j,2)=wlms(i,2,2)*sk(i+l2,2,j,2)*ox3lms(j,2,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,2,2)*sk(i+l2,2,j-1,1)*ox3lms(j-1,1,2)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j-1,1)*ox3lms(j-1,1,2)

        endif
       endif
      endif

C ... 1-2 block
      if(i .eq. ldim .and. j .gt. 1) then
      ok3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*ox3lms(j-1,1,2)

      ok3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*ox3lms(j-1,1,2)

      else

      if(j .eq. 1 .and. i .lt. ldim) then
      ok3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*ox3lms(j,2,2)

      ok3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*ox3lms(j,2,2)

      else

      if(j .eq. 1 .and. i .eq. ldim) then
      ok3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox3lms(j,2,2)

      ok3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox3lms(j,2,2)

      else

      ok3(i,1,j,2)=wlms(i,1,1)*sk(i,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i,1,j-1,1)*ox3lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+1,2,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+1,2,j-1,1)*ox3lms(j-1,1,2)

      ok3(i+l2,1,j,2)=wlms(i,1,1)*sk(i+l2,1,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,1)*sk(i+l2,1,j-1,1)*ox3lms(j-1,1,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j,2)*ox3lms(j,2,2)+
     .         wlms(i,1,2)*sk(i+l2+1,2,j-1,1)*ox3lms(j-1,1,2)

        endif
       endif
      endif

C ... 2-1 block
      if(i .eq. 1 .and. j .lt. ldim) then
      ok3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*ox3lms(j+1,2,1)

      ok3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*ox3lms(j+1,2,1)

      else

      if(j .eq. ldim .and. i .gt. 1) then
      ok3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*ox3lms(j,1,1)

      ok3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*ox3lms(j,1,1)

      else

        if(j .eq. ldim .and. i .eq. 1) then
      ok3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox3lms(j,1,1)

      ok3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox3lms(j,1,1)

        else

      ok3(i,2,j,1)=wlms(i,2,2)*sk(i,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i-1,1,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i,2,j+1,2)*ox3lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i-1,1,j+1,2)*ox3lms(j+1,2,1)

      ok3(i+l2,2,j,1)=wlms(i,2,2)*sk(i+l2,2,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j,1)*ox3lms(j,1,1)+
     .         wlms(i,2,2)*sk(i+l2,2,j+1,2)*ox3lms(j+1,2,1)+
     .         wlms(i,2,1)*sk(i+l2-1,1,j+1,2)*ox3lms(j+1,2,1)

        endif
       endif
      endif

       enddo
      enddo

C ... 3C Overlap matrix: ok = 1 + h^(+)ph

C ... Make ok3*hk = W*S*W^(T)*p*W^(T)*S*W
      do i1 = 1, 2
        do j1 = 1, 2
           do i = 1, ldim
              do j = 1, ldim
                 ok4(i,i1,j,j1) = 0d0
                 ok4(i+l2,i1,j,j1) = 0d0
                 omk1(i,i1,j,j1) = 0d0
                 omk2(i,i1,j,j1) = 0d0
                 omk3(i,i1,j,j1) = 0d0
                 omk4(i,i1,j,j1) = 0d0
              enddo
           enddo
        enddo
      enddo


      do i1 = 1, 2
         do j1 = 1, 2
            do i = 1, ldim
               do j = 1, ldim
                  do k = 1, 2
                  do l = 1, ldim
                  omk1(i,i1,j,j1) = omk1(i,i1,j,j1)
     .                            + ok3(i,i1,l,k)*hk(l,k,j,j1)
                  omk2(i,i1,j,j1) = omk2(i,i1,j,j1)
     .                            + ok3(i+l2,i1,l,k)*hk(l+l2,k,j,j1)
                  omk3(i,i1,j,j1) = omk3(i,i1,j,j1)
     .                            + ok3(i+l2,i1,l,k)*hk(l,k,j,j1)
                  omk4(i,i1,j,j1) = omk4(i,i1,j,j1)
     .                            + ok3(i,i1,l,k)*hk(l+l2,k,j,j1)
                enddo
                  enddo
               enddo
            enddo
        enddo
        enddo


         do i1 = 1, 2
         do j1 = 1, 2
            do i = 1, ldim
               do j = 1, ldim
                  ok4(i,i1,j,j1) = omk1(i,i1,j,j1) - omk2(i,i1,j,j1)
                  ok4(i+l2,i1,j,j1) = omk3(i,i1,j,j1) + omk4(i,i1,j,j1)
               enddo
            enddo
        enddo
        enddo

C ... Add all terms

C ... Make Overlap matrix for gamma-rep, 2C: ok = 1

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



      do i1 = 1, 2
      do j1 = 1, 2
       do i = 1, ldim
       do j = 1, ldim
C ... real part
           ok(i,i1,j,j1) = ok1(i,i1,j,j1) + ok2(i,i1,j,j1) +
     .                     ok4(i,i1,j,j1) + ok(i,i1,j,j1)
C ... imaginary part
            ok(i+l2,i1,j,j1) = ok1(i+l2,i1,j,j1) + ok2(i+l2,i1,j,j1) +
     .                         ok4(i+l2,i1,j,j1) + ok(i+l2,i1,j,j1)
       enddo
       enddo
      enddo
      enddo

C ... Add ox0

      do i = 1, ldim
       do i1 = 1, 2
        ok(i,i1,i,i1) = ok(i,i1,i,i1) + ox0lms(i,i1,i1)
       enddo
      enddo

      do i = 2, ldim
       ok(i-1,1,i,2) = ok(i-1,1,i,2) + ox0lms(i-1,1,2)
       ok(i,2,i-1,1) = ok(i,2,i-1,1) + ox0lms(i,2,1)
      enddo



C ... Continue with the construction of the Hamiltonian.

C ... Make hk3*hk = W*S*W^(T)*Enu*p*W^(T)*S*W
      do i1 = 1, 2
        do j1 = 1, 2
           do i = 1, ldim
              do j = 1, ldim
                 hk4(i,i1,j,j1) = 0d0
                 hk4(i+l2,i1,j,j1) = 0d0
                 mk1(i,i1,j,j1) = 0d0
                 mk2(i,i1,j,j1) = 0d0
                 mk3(i,i1,j,j1) = 0d0
                 mk4(i,i1,j,j1) = 0d0
              enddo
           enddo
        enddo
      enddo


      do i1 = 1, 2
         do j1 = 1, 2
            do i = 1, ldim
               do j = 1, ldim
                  do k = 1, 2
                  do l = 1, ldim
                  mk1(i,i1,j,j1) = mk1(i,i1,j,j1)
     .                           + hk3(i,i1,l,k)*hk(l,k,j,j1)
                  mk2(i,i1,j,j1) = mk2(i,i1,j,j1)
     .                           + hk3(i+l2,i1,l,k)*hk(l+l2,k,j,j1)
                  mk3(i,i1,j,j1) = mk3(i,i1,j,j1)
     .                           + hk3(i+l2,i1,l,k)*hk(l,k,j,j1)
                  mk4(i,i1,j,j1) = mk4(i,i1,j,j1)
     .                           + hk3(i,i1,l,k)*hk(l+l2,k,j,j1)
                enddo
                  enddo
               enddo
            enddo
        enddo
        enddo


         do i1 = 1, 2
         do j1 = 1, 2
            do i = 1, ldim
               do j = 1, ldim
                  hk4(i,i1,j,j1) = mk1(i,i1,j,j1) - mk2(i,i1,j,j1)
                  hk4(i+l2,i1,j,j1) = mk3(i,i1,j,j1) + mk4(i,i1,j,j1)
               enddo
            enddo
        enddo
        enddo

C ... Add all terms

      do i1 = 1, 2
      do j1 = 1, 2
       do i = 1, ldim
       do j = 1, ldim
C ... real part
            hk(i,i1,j,j1) = hk(i,i1,j,j1) + hk1(i,i1,j,j1) +
     .                      hk2(i,i1,j,j1) + hk4(i,i1,j,j1)
C ... imaginary part
            hk(i+l2,i1,j,j1) = hk(i+l2,i1,j,j1) + hk1(i+l2,i1,j,j1) +
     .                      hk2(i+l2,i1,j,j1) + hk4(i+l2,i1,j,j1)
       enddo
       enddo
      enddo
      enddo

C ... Add C and mx0

      do i = 1, ldim
       do i1 = 1, 2
        hk(i,i1,i,i1) = hk(i,i1,i,i1) + mx0lms(i,i1,i1) + clms(i,i1,i1)
       enddo
      enddo

C      print*,'ADDED DIAGONAL C'

      do i = 2, ldim
       hk(i-1,1,i,2) = hk(i-1,1,i,2) + mx0lms(i-1,1,2) + clms(i-1,1,2)
       hk(i,2,i-1,1) = hk(i,2,i-1,1) + mx0lms(i,2,1) + clms(i,2,1)
      enddo

C      print*, 'MADE THE HAMILTONIAN'

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

C       call yprm('hk fully rel.',12,hk,(2*ldim)**2,2*ldim,2*ldim,2*ldim)
C       call yprm('ok fully rel.',12,ok,(2*ldim)**2,2*ldim,2*ldim,2*ldim)

       call tcx('hmfr3c')

       end

