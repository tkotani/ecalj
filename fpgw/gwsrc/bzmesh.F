      SUBROUTINE BZMESH(icase,RBAS,QB,N1,N2,N3,G,NG,IPQ,QP,WGT,NQ,NQMX)
C-Divides the reciprocal lattice into microcells
Cr new version. icase=2 for regular mesh not including gamma.
Cr t.kotani. Jan2004
Cr Original version was from  bz.f of ferdi/gw.
C--------------------------------------------------------------------
Ci Inputs:
Ci  rbas, direct lattice vectors; n1,n2,n3, number of divisions to be
Ci  made along each of the reciprocal lattice vectors to divide the
Ci  reciprocal lattice into microcells; g, symmetry group operations
Ci  output from GRPGEN; ng, number of point operations output from
Ci  GRPGEN; nqmx, maximum number of k-points.
Co Outputs:
Co  ipq, wgt, nq (see remarks);
Co  qb, vectors of first microcell for input to BZINTS.
Cr Remarks:
Cr  The reciprocal lattice is divided into n1*n2*n3 microcells which
Cr  are parallelipipeds with 8 corners. The corners are nodes of the
Cr  k-space mesh in the whole reciprocal lattice unit cell. Some of
Cr  these will be symmetry-related leaving nq irreducible k-points.
Cr  These are returned in qp(3,j) j = 1,nq; for each corner defined
Cr  by the triple (i1,i2,i3), ipq(i1,i2,i3) points to the corresponding
Cr  k-vector in qp. wgt(j) then contains the sampling weight associated
Cr  with the j th irreducible k-point in qp : i.e. the number of times
Cr  it occurs in the whole reciprocal unit cell.
C---------------------------------------------------------------------
c icase==2 is for regular-mesh without gamma points.
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION RB(3,3),QB(3,3),WGT(NQMX),V(3),V1(3),IPQ(N1,N2,N3),
     .   RBAS(3,3),QBAS(3,3),QP(3,NQMX),g(3,3,1)
      real(8):: diff2(3),vv(3),diff(3),det,ginv(3,3),xv(3),xx(3)
      integer(4) :: verbose
      hf=0d0
      if(icase==2) hf=0.5d0
      call dinv33(rbas,1,qbas,vol)
      WRITE(*,"(' BZMESH : ',5X,'RBAS',31X,'QBAS')")
      DO K = 1, 3
        WRITE(*,"(3F10.5,5X,3F10.5)") (RBAS(M,K),M=1,3),(QBAS(M,K),M=1,3)
      enddo
      DO  8  M = 1, 3
        QB(M,1) = QBAS(M,1)/N1
        QB(M,2) = QBAS(M,2)/N2
        QB(M,3) = QBAS(M,3)/N3
        RB(M,1) = RBAS(M,1)*N1
        RB(M,2) = RBAS(M,2)*N2
        RB(M,3) = RBAS(M,3)*N3
    8 continue
      IPQ = 0
      W0  = 2D0/(N1*N2*N3)
      NQ  = 0
      SWGT = 0D0
      ixx=0

c      DO IG = 1, NG
c        x1=.5
c        x2=1.5
c        x3=0
c        V = matmul(QBAS(:,:),(/x1,x2,x3/))
c        CALL GRPOP(V,V1,G,IG)
c        write(6,"(' ig qbas1 rot(qbas1)=', 3f8.3,'  ',3f8.3,'  ',3f8.3)")
c     &   v,v1,matmul(v1,rbas)
c      enddo
c      stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      DO  20  I3 = 1, N3
      DO  20  I2 = 1, N2
      DO  20  I1 = 1, N1
        IF (IPQ(I1,I2,I3) .EQ. 0) THEN
          V = matmul(QB,(/(I1-1+hf),(I2-1+hf),(I3-1+hf)/))
          IWGT = 0
          DO  22  IG = 1, NG
            CALL GRPOP(V,V1,G,IG)
            XV = matmul(V1,rbas)
            call rangedq_xxx(n1*xv(1)-hf,xv1)
            call rangedq_xxx(n2*xv(2)-hf,xv2)
            call rangedq_xxx(n3*xv(3)-hf,xv3)

c            write(6,"(' i1 i2 i3 ig=', 4i4,2x,3d13.6)") i1,i2,i3,ig
c     &      , n1*xv(1), n2*xv(2), n3*xv(3)
c            if(abs(xv1)+abs(xv2)+abs(xv3)>4d-6) then
c              write(6,"('             ', 18x,3d13.6)") xv1, xv2, xv3
c            else
c              write(6,"('             ', 18x,3d13.6,' go')") xv1, xv2, xv3
c            endif

            if(abs(xv1)+abs(xv2)+abs(xv3)>3d-6) cycle
            ! this is for the case not on the regular mesh points.
c
            call rangedq_xxx(xv(1)-hf/n1,xv1)
            call rangedq_xxx(xv(2)-hf/n2,xv2)
            call rangedq_xxx(xv(3)-hf/n3,xv3)
c              write(6,"('             ', 18x,3d13.6,' xx')") xv1, xv2, xv3
            j1 = xv1*n1+1 +1d-6
            j2 = xv2*n2+1 +1d-6
            j3 = xv3*n3+1 +1d-6
            if(verbose()>=100) 
     &        write(6,"(' i1 i2 i3 j1 j2 j3=', 3i4,2x,3i4)") i1,i2,i3, j1,j2,j3
            if(abs(j1-xv1*n1-1)/n1+abs(j2-xv2*n2-1)/n2+abs(j3-xv3*n3-1)/n3>3d-6) then
              print *,'xxx1=',abs(j1-xv1*n1-1)
              print *,'xxx2=',abs(j2-xv2*n2-1)
              print *,'xxx3=',abs(j3-xv3*n3-1)
              stop 'bzmesh: bug! abs(j1-xv1*n1-1)+... >4d-6'
            endif
            IF (IPQ(J1,J2,J3) .EQ. 0) then
              IWGT = IWGT+1
            endif
            IPQ(J1,J2,J3) = NQ+1
            ixx  = ixx+1
   22     continue
          NQ = NQ+1
          QP(:,NQ) = V
          WGT(NQ) = IWGT*W0
          SWGT = SWGT + WGT(NQ)
        ENDIF
   20 CONTINUE
      if(verbose()>=100) print *,' Number of Qibz =',nq

cccccccccccccccccccccccccccc
c      stop 'test end of bzmesh 111 ----------------------------'
cccccccccccccccccccccccccccc

c-------------------------------------
c rotation check
      call dinv33(qbas,0,ginv,det)
      DO  I3 = 1, N3
        DO  I2 = 1, N2
          DO  I1 = 1, N1
            V = matmul(QB,(/(I1-1+hf),(I2-1+hf),(I3-1+hf)/))
c        write(6,"(' ===== qbz=',3f8.3)") matmul(ginv,v)
            do k  = 1,nq
              do ir = 1,ng
c        if(verbose()>104) print *,' grp=',ir !;      print *, grp(:,ir)
                diff = matmul(g(:,:,ir),qp(:,k)) - v
                call rangedq3(matmul(ginv,diff), diff2)
c          if(verbose()>104) write(6,"(' matmul(ginv,diff)=',3f8.3,' ',3f8.3)")
c     &      diff, matmul(ginv,diff)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                if(verbose()>60) then
                  write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5,3d13.5)")
     &i1,i2,i3,v,qp(:,k),k,diff2 ! matmul(qp(:,NQx),rbas)*8, matmul(v, rbas)*8
                endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                if(sum(abs(diff2))< 1d-6) then
                  if(verbose()>50) 
     &       write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5)")
     &       i1,i2,i3,v, qp(:,k),k ! matmul(qp(:,NQx),rbas)*8, matmul(v, rbas)*8
c            irotk(k,ir)= kp
c            nstar(k)   = nstar(k) + 1
                  goto 1022
                endif

              enddo
            enddo
            stop'bzmesh: not find irotk; it may require accurate symmetry.'
 1022       continue
          enddo
        enddo
      enddo
c-----------------------------------------
      if (dabs(swgt-2) .gt. 1.d-9) then
        print *,' swgt=',swgt
        stop 'BZMESH: error in weights'
      endif
      write(*,744) nq,n1*n2*n3,n1,n2,n3
  744 FORMAT(/' BZMESH: ',I4,' irreducible k - points. from ',i5,
     .          ' (',3i3,' )')
      write(6,663)
  663 format(13x,'Qx',8x,'Qy',8x,'Qz',6x,'Multiplicity    Weight')
      do  51  iq = 1, nq
        iwgt = wgt(iq)/w0+ .1
        write(6,661) iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,wgt(iq)
   51 continue
  661 format(i5,2x,3f10.4,i10,f16.6)
cccccccccccccccccccccccccccc
c      stop 'test end of bzmesh222 ----------------------------'
cccccccccccccccccccccccccccc
      END


c----------------------------------------------------





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This was old version from  bz.f of ferdi/gw.
c--------------------------------------------------------------------
      SUBROUTINE BZMESH_old(icase,RBAS,QB,N1,N2,N3,G,NG,IPQ,QP,WGT,NQ,NQMX)
C-Divides the reciprocal lattice into microcells
C--------------------------------------------------------------------
Ci Inputs:
Ci  rbas, direct lattice vectors; n1,n2,n3, number of divisions to be
Ci  made along each of the reciprocal lattice vectors to divide the
Ci  reciprocal lattice into microcells; g, symmetry group operations
Ci  output from GRPGEN; ng, number of point operations output from
Ci  GRPGEN; nqmx, maximum number of k-points.
Co Outputs:
Co  ipq, wgt, nq (see remarks);
Co  qb, vectors of first microcell for input to BZINTS.
Cr Remarks:
Cr  The reciprocal lattice is divided into n1*n2*n3 microcells which
Cr  are parallelipipeds with 8 corners. The corners are nodes of the
Cr  k-space mesh in the whole reciprocal lattice unit cell. Some of
Cr  these will be symmetry-related leaving nq irreducible k-points.
Cr  These are returned in qp(3,j) j = 1,nq; for each corner defined
Cr  by the triple (i1,i2,i3), ipq(i1,i2,i3) points to the corresponding
Cr  k-vector in qp. wgt(j) then contains the sampling weight associated
Cr  with the j th irreducible k-point in qp : i.e. the number of times
Cr  it occurs in the whole reciprocal unit cell.
Cm Memory:
Cm  No significant internal requirement. No access to heap.
C---------------------------------------------------------------------
c icase==2 is for regular-mesh without gamma points.
c
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION RB(3,3),QB(3,3),WGT(NQMX),V(3),V1(3),IPQ(N1,N2,N3),
     .   RBAS(3,3),QBAS(3,3),QP(3,NQMX),g(*)
c      print *,'bzmesh: iprint=',iprint(),n1,n2,n3
      hf=0d0
      if(icase==2) hf=0.5d0
      call dinv33(rbas,1,qbas,vol)
      if (iprint() > 50) then
        WRITE(*,"(' BZMESH : ',5X,'RBAS',31X,'QBAS')")
        DO K = 1, 3
          WRITE(*,"(3F10.5,5X,3F10.5)") (RBAS(M,K),M=1,3),(QBAS(M,K),M=1,3)
        enddo
      endif
      DO  8  M = 1, 3
        QB(M,1) = QBAS(M,1)/N1
        QB(M,2) = QBAS(M,2)/N2
        QB(M,3) = QBAS(M,3)/N3
        RB(M,1) = RBAS(M,1)*N1
        RB(M,2) = RBAS(M,2)*N2
        RB(M,3) = RBAS(M,3)*N3
    8 continue
      IPQ = 0  !fortran 90
      W0 = 2D0/(N1*N2*N3)
      NQ = 0
      SWGT = 0D0
      NN1 = 6*N1-1
      NN2 = 6*N2-1
      NN3 = 6*N3-1
      DO  20  I3 = 1, N3
      DO  20  I2 = 1, N2
      DO  20  I1 = 1, N1
        IF (IPQ(I1,I2,I3) .EQ. 0) THEN
          V(1) = (I1-1+hf)*QB(1,1) + (I2-1+hf)*QB(1,2) + (I3-1+hf)*QB(1,3)
          V(2) = (I1-1+hf)*QB(2,1) + (I2-1+hf)*QB(2,2) + (I3-1+hf)*QB(2,3)
          V(3) = (I1-1+hf)*QB(3,1) + (I2-1+hf)*QB(3,2) + (I3-1+hf)*QB(3,3)
          IWGT = 0
          DO  22  IG = 1, NG
            CALL GRPOP(V,V1,G,IG)
            X1 = V1(1)*RB(1,1) + V1(2)*RB(2,1) + V1(3)*RB(3,1) + .1
            X2 = V1(1)*RB(1,2) + V1(2)*RB(2,2) + V1(3)*RB(3,2) + .1
            X3 = V1(1)*RB(1,3) + V1(2)*RB(2,3) + V1(3)*RB(3,3) + .1
            IF (X1 .GT. 0.D0) X1 = X1+1.D0
            IF (X2 .GT. 0.D0) X2 = X2+1.D0
            IF (X3 .GT. 0.D0) X3 = X3+1.D0
            J1 = X1
            J2 = X2
            J3 = X3
            J1 = MOD(J1+NN1,N1) + 1
            J2 = MOD(J2+NN2,N2) + 1
            J3 = MOD(J3+NN3,N3) + 1
            IF (IPQ(J1,J2,J3) .EQ. 0) IWGT = IWGT+1
C temp: try to find bug ...
            if (j1 .lt. 0) j1 = j1+n1
            if (j2 .lt. 0) j2 = j2+n2
            if (j3 .lt. 0) j3 = j3+n3
            IPQ(J1,J2,J3) = NQ+1
   22     continue
          NQ = NQ+1
          QP(1,NQ) = V(1)
          QP(2,NQ) = V(2)
          QP(3,NQ) = V(3)
          WGT(NQ) = IWGT*W0
          SWGT = SWGT + WGT(NQ)
        ENDIF
   20 CONTINUE
      if (dabs(swgt-2) .gt. 1.d-9) stop 'BZMESH: error in weights'
      if (iprint() .ge. 30) then
cfa
c       do  30  j = 1, 2
c  30   write(lgunit(j),744) nq,n1*n2*n3,n1,n2,n3
        write(*,744) nq,n1*n2*n3,n1,n2,n3
  744   FORMAT(/' BZMESH: ',I4,' irreducible k - points. from ',i5,
     .          ' (',3i3,' )')
      endif
      if (iprint() .gt. 40) then
        write(6,663)
  663   format(13x,'Qx',8x,'Qy',8x,'Qz',6x,'Multiplicity    Weight')
        do  51  iq = 1, nq
          iwgt = wgt(iq)/w0+ .1
          write(6,661) iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,wgt(iq)
   51   continue
  661   format(i5,2x,3f10.4,i10,f16.6)
      endif
      END
      SUBROUTINE GRPOP(V,V1,G,I)
      double precision G(3,3,*),V(3),V1(3)
      integer i
      V1(1) = G(1,1,I)*V(1) + G(1,2,I)*V(2) + G(1,3,I)*V(3)
      V1(2) = G(2,1,I)*V(1) + G(2,2,I)*V(2) + G(2,3,I)*V(3)
      V1(3) = G(3,1,I)*V(1) + G(3,2,I)*V(2) + G(3,3,I)*V(3)
      RETURN
      END
