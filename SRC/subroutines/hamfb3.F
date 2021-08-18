      subroutine hamfb3k(qin,iq,nk1,nk2,nk3,
c    i k1,k2,k3,ipq,napw_in,ndimh,ldima,ldimb,qb,plat,qlat,ldim,ifac,igstar,
     i k1,k2,k3,ipq,napw_in,ndimh,ldima,ldimb,qb,ldim,ifac,igstar,
     i hq,iaf, !iaf at jun2015
     o gfbz)
!! Generate gfbz (full BZ) from hq (1st BZ).
Ci q iq: q and index for q
Ci   nk1,nk2,nk3:  no. divisions for the qp in 3 recip. latt. vecs
Ci   k1,k2,k3: leading dimensions of gfbz
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci-temp  igstar:contains group operation needed to rotated a k-point
Ci          to its irreducible one (bzmesh.f) !this is just for check.
Ci   ndimh: dim of self energy hq at qin
Ci   ldima :dimensions gfbz; also the number of rows in gfbz to fill.
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   ldimb :dimensions gfbz; also the number of columns in gfbz to fill
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   hq    : Self-energy (or non-local potential for this iq
Ci read date in m_hamindex .
Co Outputs
Co   gfbz  : For those qp in star iq, hq stored
C ----------------------------------------------------------------------
      implicit none
      integer:: nk1,nk2,nk3,k1,k2,k3,ipq(*),igstar(0:*),ndimh,ldima,ldimb,napw_in,debugmode
      real(8)::    qin(3),qb(3,3) !,plat(3,3),qlat(3,3)
      complex(8):: hq(ndimh,ndimh),gfbz(k1,k2,k3,ldima,ldimb)
      integer:: i,i1,i2,i3,ig,iq,iq1,is,j,jj1,jj2,jj3,js,k,nl,ierr,ifac(3),j1,j2,ik1,ik2,ik3,isp,ldim,iaf
      real(8):: q1(3),qk
      character(200)::aaa
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .(jj2*ifac(2)-1)*qb(k,2) +
     .(jj3*ifac(3)-1)*qb(k,3)
      call tcn('hamfb3k')
      if(debugmode()>0) print *, 'hamfb3k: start...'
      iq1 = 0
      do  i3 = 1, nk3
        do  i2 = 1, nk2
          do  i1 = 1, nk1
            iq1 = iq1+1
            if(debugmode()>0) print *,'iq iq1 ipq(iq1)',iq,iq1,ipq(iq1)
!! ipq(iq1) ist gives a table to point irreducible point.
!!          q1,iq1 is target on regular mesh <--- qin,iq is irreducible points; this mapping is by rotsig.
!!          iq=ipq(iq1) shows iq for given iq1.
            if (ipq(iq1) .ne. iq) cycle !this must make things efficient
            q1(1) = qk(1,i1,i2,i3)
            q1(2) = qk(2,i1,i2,i3)
            q1(3) = qk(3,i1,i2,i3)
            if(debugmode()>0) print *,q1
c        ig = igstar(iq1)
c        write(6,"(a,i4)") 'symops',ig
c        write(6,"(3f10.5)")symops(1,:,ig)
c        write(6,"(3f10.5)")symops(2,:,ig)
c        write(6,"(3f10.5)")symops(3,:,ig)
            if(debugmode()>0) write(6,"(a,3f13.5)")' input          qin = ', qin
            if(debugmode()>0) write(6,"(a,3f13.5)")' target a        q1 = ', q1 ! qin = g(:,:,ig)^{-1}*q1
c        write(6,"(a,2i3, 3f13.5)")' target b   qp = ',iq1,ig, matmul(symops(:,:,ig),qin)-q1
c       write(6,"(a,3f13.5)")
            call rotsig(qin,q1,ndimh,napw_in,ldim,hq,gfbz(i1,i2,i3,:,:),ierr,iaf)
            if(ierr/=0) write(aaa,"(' qin=',3f13.6,' q1=',3f13.6)") qin,q1
            if(ierr/=0) call rx('hamfb3: rotsig do not map qin to q1;'//trim(aaa))
          enddo
        enddo
      enddo
      if(debugmode()>0) print *, 'hamfb3k: end...'
      call tcx('hamfb3k')
c      write(6,"(a,i3,3f13.5)")'end of hamfb3k: ierr qin=',ierr,qin
c      if(sum(abs(qin-[0d0,0d0,0.03295d0]))<1d-3.and.ierr==1) stop 'test end xxxxxxxxxxxxxxxxxxxxxxxx'
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$       isp=1
c$$$       open(3013,file="sfz_gfbz")
c$$$       do j1=1,ldima
c$$$       do j2=1,ldimb
c$$$        do ik1=1,k1
c$$$        do ik2=1,k2
c$$$        do ik3=1,k3
c$$$          write(3013,"(6i4,2d13.5)")ik1,ik2,ik3,j1,j2,isp,gfbz(ik1,ik2,ik3,j1,j2)
c$$$c           write(3012,"(6i4,2d13.5)")ik1,ik2,ik3,j1,j2,isp,sfz(ik1,ik2,ik3,j1,j2)
c$$$        enddo
c$$$        enddo
c$$$        enddo
c$$$       enddo
c$$$       enddo
c$$$       close(3013)
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
      end

