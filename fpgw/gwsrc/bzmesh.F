      SUBROUTINE BZMESH(RBAS,QB,N1,N2,N3,G,NG,IPQ,QP,WGT,NQ,NQMX,nadd,qbzshift) !icase,
!!-Divides the reciprocal lattice into microcells
!! icase=1 only.
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

      IMPLICIT double precision (A-H,O-Z)
      DIMENSION RB(3,3),QB(3,3),WGT(NQMX),V(3),V1(3),IPQ(N1+nadd,N2+nadd,N3+nadd),
     .   RBAS(3,3),QBAS(3,3),QP(3,NQMX),g(3,3,ng)
      real(8):: diff2(3),vv(3),diff(3),det,ginv(3,3),xv(3),xx(3),vi(3),qbzshift(3),xvv(3)
      integer(4) :: verbose
      integer:: nadd,nnn
      integer:: ix=0,iii,jj(3)
c      hf=0d0
c      if(icase==2) hf=0.5d0
      call minv33(rbas,qbas)
      WRITE(*,"(' BZMESH : ',5X,'RBAS',31X,'QBAS')")
      DO K = 1, 3
        WRITE(*,"(3F10.5,5X,3F10.5)") (RBAS(M,K),M=1,3),(QBAS(M,K),M=1,3)
      enddo
      QB(:,1) = QBAS(:,1)/N1
      QB(:,2) = QBAS(:,2)/N2
      QB(:,3) = QBAS(:,3)/N3
      RB(:,1) = RBAS(:,1)*N1
      RB(:,2) = RBAS(:,2)*N2
      RB(:,3) = RBAS(:,3)*N3
      W0  = 2D0/(N1*N2*N3)
      nnn = (N1+nadd)*(N2+nadd)*(N3+nadd)

!! Gamma point first ---
      if(nadd==1) w0= 2d0/nnn
      NQ = 1
      QP(:,1) = 0d0
      WGT(1) = W0
      SWGT = 1d0*W0
      IPQ = 0
      do I3 = 1, N3+nadd
      do I2 = 1, N2+nadd
      do I1 = 1, N1+nadd
        V = matmul(QB,(/(I1-1+qbzshift(1)),(I2-1+qbzshift(2)),(I3-1+qbzshift(3))/))
        if(sum(abs(v))<1d-6) then
          J1=I1
          J2=I2
          J3=I3
        endif
      enddo
      enddo
      enddo
      IPQ(J1,J2,J3) = NQ

!! All other points
      iii=1
      print *,'w0=',w0
      print *,'qbzshift=',qbzshift
      DO  20  I3 = 1, N3+nadd
      DO  20  I2 = 1, N2+nadd
      DO  20  I1 = 1, N1+nadd
        V = matmul(QB,(/(I1-1+qbzshift(1)),(I2-1+qbzshift(2)),(I3-1+qbzshift(3))/))
c        if(sum(abs(v))<1d-6) then
          print *,'0000',i1,i2,i3,v
c        endif
        IF (IPQ(I1,I2,I3) .EQ. 0) THEN
          IWGT = 0
          NQ = NQ+1
          do 22 IG = 1, NG
            v1 = matmul(g(:,:,ig),v(:))
            XV = matmul(V1,rb)
            jj(:)  = nint(xv(:)+1-qbzshift(:))
            xvv(:) = xv(:)+1-qbzshift(:) - jj
            if(nadd==0) then !periodicity. Use mod(x+|x|*n1, n1) so that x+|x|*n1 >=0.
              jj(1)= mod(jj(1)-1+n1*abs(jj(1)-1),n1)+1
              jj(2)= mod(jj(2)-1+n2*abs(jj(2)-1),n2)+1
              jj(3)= mod(jj(3)-1+n3*abs(jj(3)-1),n3)+1
            endif  
cc            write(6,"('aaaaaaaab111 ',i3,3f10.5,3i3,x,3f10.5)") ig,xv,jj(:),v1
c            call rangedq_xxx(xv(1)+1-qbzshift(1),xv1)
c            call rangedq_xxx(xv(2)+1-qbzshift(2),xv2)
c            call rangedq_xxx(xv(3)+1-qbzshift(3),xv3)
            if(sum(abs(xvv(:)))>3d-6) then
              if(nadd==0) then
                print *,' bzmesh: mapping i1i2i3 to: ',i1,i2,i3,xvv(:)
                call rx( 'bzmesh: bug! not mapped to regular mesh point')
              else
                cycle 
              endif
            endif
            j1=jj(1)
            j2=jj(2)
            j3=jj(3)
            if(j1<1 .or. j1>n1+nadd) cycle
            if(j2<1 .or. j2>n2+nadd) cycle
            if(j3<1 .or. j3>n3+nadd) cycle
            if(IPQ(J1,J2,J3) == 0) then
              IWGT = IWGT+1
            endif  
ccccccccccccccccccc
c            IF (IPQ(J1,J2,J3) .EQ. 0) then
c              write(6,"('aaab ',i3,f10.5,3i3,' xxx')") ig,sum(abs(xvv(:))),j1,j2,j3
c            else
c              write(6,"('aaab ',i3,f10.5,3i3)") ig,sum(abs(xvv(:))),j1,j2,j3
c            endif
ccccccccccccccc
            IPQ(J1,J2,J3) = NQ
   22     continue
          QP(:,NQ) = V
          WGT(NQ) = IWGT*W0
          SWGT = SWGT + WGT(NQ)
c          III = III + IWGT
          print *,'sssss',i1,i2,i3,swgt,iwgt,iii
        ENDIF
ccccccccccccccccccc
c        if(i2==i3.and.i1==1) then
c        if( i1==n1+nadd.or.i2==n2+nadd.or.i3==n3+nadd) then
c          VI = matmul(QB,(/(I1-1+hf),(I2-1+hf),(I3-1+hf)/))
c          write(6,"('qqqq: q=',3f11.4,2x,3f11.4,i5)") VI,V,IPQ(I1,I2,I3)
c        endif
cccccccccccccccccccc
   20 CONTINUE
      DO   I3 = 1, N3+nadd
      DO   I2 = 1, N2+nadd
      DO   I1 = 1, N1+nadd
        write(6,*)i1,i2,i3,ipq(i1,i2,i3)
      enddo
      enddo
      enddo
c      if(verbose()>=100) 
      write(6,"(' Number of Qibz =',2i6)") nq,nnn

!! Rotate
c ------------------------------------
      call minv33(qbas,ginv)
      DO  I3 = 1, N3+nadd
        DO  I2 = 1, N2+nadd
          DO  I1 = 1, N1+nadd
            V = matmul(QB,(/(I1-1),(I2-1),(I3-1)/))
c        write(6,"(' ===== qbz=',3f8.3)") matmul(ginv,v)
            do k  = 1,nq
              do ir = 1,ng
c        if(verbose()>104) print *,' grp=',ir !;      print *, grp(:,ir)
                diff = matmul(g(:,:,ir),qp(:,k)) - v
                call rangedq3(matmul(ginv,diff), diff2)
c          if(verbose()>104) write(6,"(' matmul(ginv,diff)=',3f8.3,' ',3f8.3)")
c     &      diff, matmul(ginv,diff)
                if(verbose()>60) then
                  write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5,3d13.5)")
     &             i1,i2,i3,v,qp(:,k),k,diff2 ! matmul(qp(:,NQx),rbas)*8, matmul(v, rbas)*8
                endif
                if(sum(abs(diff2))< 1d-6) then
                  if(verbose()>50) 
     &             write(6,"(' i1i2i3= ',3i3,' v qp k=',2x,3f7.3,2x,3f7.3,i5)")
     &             i1,i2,i3,v, qp(:,k),k ! matmul(qp(:,NQx),rbas)*8, matmul(v, rbas)*8
c            irotk(k,ir)= kp
c            nstar(k)   = nstar(k) + 1
                  goto 1022
                endif
              enddo
            enddo
            call rx( 'bzmesh: not find irotk; it may require accurate symmetry.')
 1022       continue
          enddo
        enddo
      enddo
c-----------------------------------------
      if (dabs(swgt-2d0)> 1.d-9) then
        print *,' swgt=',swgt,swgt/w0*2d0/nnn
        call rx( 'BZMESH: error in weights')
      endif
      write(*,744) nq,(n1+nadd)*(n2+nadd)*(n3+nadd),n1+nadd,n2+nadd,n3+nadd
  744 FORMAT(/' BZMESH: ',I4,' irreducible k - points. from ',i5,' (',3i3,' )')
      write(6,663)
  663 format(13x,'Qx',8x,'Qy',8x,'Qz',6x,'Multiplicity    Weight')
      do  51  iq = 1, nq
        iwgt = wgt(iq)/w0+ .1
        write(6,661) iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,wgt(iq)
   51 continue
  661 format(i5,2x,3f10.4,i10,f16.6)
      END
