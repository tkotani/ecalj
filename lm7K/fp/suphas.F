      subroutine suphas(q,p,ng,iv,n1,n2,n3,qlat,cosgp,singp)
C- Makes exp(-i p* (q+G)) for a list of reciprocal lattice vectors
C ----------------------------------------------------------------------
Ci Inputs
Ci   q     :Bloch wave number to be added to G
Ci   p     :position
Ci   ng    :number of G-vectors
Ci   iv    :G-vectors, as integer multiples of qlat
Ci   n1    :maximum first component of iv
Ci   n2    :maximum second component of iv
Ci   n3    :maximum third component of iv
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Co Outputs
Ci   cosgp :cos(-i p (q+G))
Ci   singp :sin(-i p (q+G))
Cr Remarks
Cr   iv may be calculated from setup suphs0.
Cu Updates
Cu   19 May 00 Adapted from nfp su_phs.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ng,iv(ng,3),n1,n2,n3
      double precision p(3),cosgp(ng),singp(ng),qlat(3,3),q(3)
C ... Local parameters
      integer nmx,i,k
      parameter ( nmx=1000 )
      double precision scalp,tpi
      double complex cph1(-nmx:nmx),cph2(-nmx:nmx),cph3(-nmx:nmx),
     .   cph0,cc

      call tcn('suphas')
      tpi = 8d0*datan(1d0)
      if (max0(n1,n2,n3) .gt. nmx) call rx('suphas: nmx exceeded')

C ... Done straightforwardly
C      do  i = 1, ng
C        scalp=-alat*(p(1)*g(i,1)+p(2)*g(i,2)+p(3)*g(i,3))
C        cosgp(i)=dcos(scalp)
C        singp(i)=dsin(scalp)
C      enddo

C ... Set up phase factors to multiply together as products of
C     phases in the three dimensions.
      scalp = -tpi*(q(1)*p(1)+q(2)*p(2)+q(3)*p(3))
      cph0 = dcmplx(dcos(scalp),dsin(scalp))

      scalp = -tpi*(p(1)*qlat(1,1)+p(2)*qlat(2,1)+p(3)*qlat(3,1))
      cc = dcmplx(dcos(scalp),dsin(scalp))
      cph1(0) = (1d0,0d0)
      do  k = 1, n1
        cph1(k) = cc*cph1(k-1)
        cph1(-k) = dconjg(cc)*cph1(-k+1)
      enddo

      scalp = -tpi*(p(1)*qlat(1,2)+p(2)*qlat(2,2)+p(3)*qlat(3,2))
      cc = dcmplx(dcos(scalp),dsin(scalp))
      cph2(0) = (1d0,0d0)
      do  k = 1, n2
        cph2(k) = cc*cph2(k-1)
        cph2(-k) = dconjg(cc)*cph2(-k+1)
      enddo

      scalp = -tpi*(p(1)*qlat(1,3)+p(2)*qlat(2,3)+p(3)*qlat(3,3))
      cc = dcmplx(dcos(scalp),dsin(scalp))
      cph3(0) = (1d0,0d0)
      do  k = 1, n3
        cph3(k) = cc*cph3(k-1)
        cph3(-k) = dconjg(cc)*cph3(-k+1)
      enddo

      do  i = 1, ng
        cc = cph1(iv(i,1))*cph2(iv(i,2))*cph3(iv(i,3))*cph0
        cosgp(i) = dble(cc)
        singp(i) = dimag(cc)
      enddo

      call tcx('suphas')

      end

      subroutine suphs0(plat,ng,gv,iv)
C- Setup for suphas: write each gv as linear combination of qlat.
Co iv: gv as multiples of qlat 
C     implicit none
C ... Passed parameters
      integer ng,iv(ng,3)
      double precision gv(ng,3),plat(3,3)
C ... Local parameters
      integer i
      double precision s1,s2,s3

      do  i = 1, ng
        s1 = gv(i,1)*plat(1,1)+gv(i,2)*plat(2,1)+gv(i,3)*plat(3,1)
        s2 = gv(i,1)*plat(1,2)+gv(i,2)*plat(2,2)+gv(i,3)*plat(3,2)
        s3 = gv(i,1)*plat(1,3)+gv(i,2)*plat(2,3)+gv(i,3)*plat(3,3)
        iv(i,1) = nint(s1)
        iv(i,2) = nint(s2)
        iv(i,3) = nint(s3)

C        write (6,890) i,gv(i,1),gv(i,2),gv(i,3),s1,s2,s3,
C     .     iv(i,1),iv(i,2),iv(i,3)
C  890   format(i4,' gv',3f10.4,'  si',3f7.3,'  iv',3i4)
      enddo
      end

