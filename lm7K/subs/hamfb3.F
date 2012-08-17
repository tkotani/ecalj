      subroutine hamfb3k(qin,iq,nk1,nk2,nk3,
c    i k1,k2,k3,ipq,napw_in,ndimh,ldima,ldimb,qb,plat,qlat,ldim,ifac,igstar,
     i k1,k2,k3,ipq,napw_in,ndimh,ldima,ldimb,qb,ldim,ifac,igstar,
     i hq,
     o gfbz)
c      use m_hamindex
C- Generate gfbz (full BZ) from hq (1st BZ).
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
      real(8)::    hq(ndimh,ndimh),qin(3),qb(3,3) !,plat(3,3),qlat(3,3)
      complex(8):: gfbz(k1,k2,k3,ldima,ldimb)
      integer:: i,i1,i2,i3,ig,iq,iq1,is,j,jj1,jj2,jj3,js,k,nl,ierr,ifac(3),j1,j2,ik1,ik2,ik3,isp,ldim
      real(8):: q1(3),qk
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
            if (ipq(iq1) .ne. iq) cycle
            q1(1) = qk(1,i1,i2,i3)
            q1(2) = qk(2,i1,i2,i3)
            q1(3) = qk(3,i1,i2,i3)
            if(debugmode()>0) print *,q1
c        ig = igstar(iq1)
c        write(6,"(a,i4)") 'symops',ig
c        write(6,"(3f10.5)")symops(1,:,ig)
c        write(6,"(3f10.5)")symops(2,:,ig)
c        write(6,"(3f10.5)")symops(3,:,ig)
        write(6,"(a,3f13.5)")' input          qin = ', qin
        write(6,"(a,3f13.5)")' target a        q1 = ', q1   ! qin = g(:,:,ig)^{-1}*q1
c        write(6,"(a,2i3, 3f13.5)")' target b   qp = ',iq1,ig, matmul(symops(:,:,ig),qin)-q1
c       write(6,"(a,3f13.5)")
c            call rotsig(qin,q1,ndimh,napw_in,plat,qlat,ldim,hq,gfbz(i1,i2,i3,:,:),ierr)
            call rotsig(qin,q1,ndimh,napw_in,ldim,hq,gfbz(i1,i2,i3,:,:),ierr)
          enddo
        enddo
      enddo
      if(debugmode()>0) print *, 'hamfb3k: end...'
      call tcx('hamfb3k')

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


      subroutine hamfb3(nbas,nl,offH,iprmb,opt,pos,iq,nk1,nk2,nk3,
     .k1,k2,k3,ipq,istab,g,ag,igstar,ifac,lidim,ldima,ldimb,nspc,
     .qb,hq,hw,hw2,gfbz)
C- Kernel called by hamfbz: contr. to entire gf in full BZ from 1 qp.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   opt   :specifies rotation options (see roth.f)
Ci         :1s digit is not used.  (calls roth using 1s digit = 0)
Ci         10s digit
Ci         :0 use real rotation mat for real harmonics
Ci         :1 use complex rotation matrix for cubic harmonics
Ci         :Add 4 if phase convention phi = q * [(g R_j + a) - R_i]
Ci         :  should be scaled by -1
Ci        100s digit distinguishes how complex arithmetic is handled
Ci             for INPUT h,h0 (kcplx mode)
Ci           0: h,h0 have real, imaginary separated
Ci              h = h(ldha,ldhb,2), with h(*,*,1..2) = real..imag
Ci           1: h,h0 are in complex*16 format (see Bugs)
Ci              h = h(2,ldha,ldhb), with s(1,*) = real, s(2,*) = imag
Ci           2: h,h0 have real, imaginary separated by columns
Ci              h = h(ldha,2,ldhb), with h(*,1..2,*) = real..imag
Ci       1000s digit specifies what part of g is to be extracted
Ci         :0 row dimension consists only of lower block
Ci         :1 same as zero
Ci         :2 not allowed yet
Ci         :3 not allowed yet
Ci         :4 row dimension consists lower+intermediate block
Ci   pos   :basis vectors
Ci   iq    :rotate gf to all q in star of iq
Ci   nk1,nk2,nk3:  no. divisions for the qp in 3 recip. latt. vecs
Ci   k1,k2,k3: leading dimensions of gfbz
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   igstar:contains group operation needed to rotated a k-point
Ci          to its irreducible one (bzmesh.f)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   lidim :number of lower+intermediate orbitals, defining dim. of hq
Ci   ldima :dimensions gfbz; also the number of rows in gfbz to fill.
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   ldimb :dimensions gfbz; also the number of columns in gfbz to fill
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   hq    :Green's function or hamiltonian for this iq
Ci   hw    :work array of same dimension as hq
Ci   hw2   :same work array as hw (used internally with diff. dim)
Co Outputs
Co   gfbz  :for those qp in star iq, hq stored
Co         :NB: gfbz is ALWAYS returned in complex*16 mode, regardless
Co         :    of kxplx mode for input h,h0.
Cu Updates
Cu   20 Jun 02  First cut
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer ldima,lidim,ldimb,nspc,opt
      integer nkap0,n0H
      parameter (nkap0=3,n0H=5)
      integer nbas,nk1,nk2,nk3,k1,k2,k3,istab(nbas,1),
     .ipq(1),igstar(0:*),ifac(3),offH(n0H,nkap0,nbas),iprmb(ldima)
      double precision hq(lidim,nspc,lidim,nspc)
      double precision hw(lidim,nspc,lidim,nspc)
      double precision hw2(lidim,2,nspc,lidim,nspc)
      double precision pos(3,*),g(3,3,*),ag(3,*),qb(3,3)
      double complex gfbz(k1,k2,k3,ldima,nspc,ldimb,nspc)
C Local variables
      integer:: i , i1 , i2 , i3 , ig , iq , iq1 , is , j , jj1 , jj2 
     ., jj3 , js , k , nl , kcplx , optrot , lidimx , offi
      complex(8) ,allocatable :: wk_zv(:)

      double precision q1(3),qk

      integer:: roth, i_dummy_roth
C heap:
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .(jj2*ifac(2)-1)*qb(k,2) +
     .(jj3*ifac(3)-1)*qb(k,3)

      call tcn('hamfb3')

C     optrot = opt, with 1s digit stripped
      optrot = opt + 0 - mod(opt,10)
      lidimx = lidim*nspc
      kcplx = mod(opt/100,10)
C     rotate to real here, so roth doesn't do it repeatedly
      if (kcplx .eq. 1) then
C       call zprm('hq before ztoyy',2,hq,lidimx,lidimx,lidimx)
        call ztoyy(hq,lidimx,lidimx,lidimx,lidimx,1,2)
C       roth works with kcplx=2 mode
        optrot = optrot + 200-100
      endif

C     Temporary memory for roth
      allocate(wk_zv(lidimx**2))


C --- Copy gf for each qp in star of qp(iq) ---
      iq1 = 0
      do  12  i3 = 1, nk3
      do  12  i2 = 1, nk2
      do  12  i1 = 1, nk1

        iq1 = iq1+1
C   ... skip this qp unless it is related to iq
        if (ipq(iq1) .ne. iq) goto 12

C   ... Make g by rotation of g(iq): symop relating iq1 to iq
        ig = igstar(iq1)

C   ... q into which h is rotated
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

        if (nspc .ne. 1) call rx('hamfb3: check offsets for nspc')
C   ... Copy hq into hw to preserve gii; use hw
        call dcopy(2*lidimx**2,hq,1,hw,1)
        offi = lidimx**2
        if (kcplx .ne. 0) offi = lidimx

C   ... Rotate and copy
        do  14  is = 1, nspc
        do  14  js = 1, nspc
          if (kcplx .eq. 0) then
            if (ig .ne. 1) then
Ckino roth is integer function,                     call roth ( optrot , nl , nbas , pos , 0 , offh , iprmb , istab
              i_dummy_roth=roth( optrot , nl , nbas , pos , 0 , offh , iprmb , istab 
     .              ( 1 , ig ) , g ( 1 , 1 , ig ) , ag ( 1 , ig ) , q1 , lidim , 
     .              lidim , wk_zv , hw ( 1 , is , 1 , js ) )

            endif

            do  16  j  = 1, ldimb
            do  16  i  = 1, ldima
              gfbz(i1,i2,i3,i,is,j,js) =
     .                dcmplx(hw(i,is,j,js),hw(i+offi,is,j,js))
   16       continue
          else
            if (ig .ne. 1) then
Ckino roth is integer function,                     call roth ( optrot , nl , nbas , pos , 0 , offh , iprmb , istab
              i_dummy_roth=roth( optrot , nl , nbas , pos , 0 , offh , iprmb , istab 
     .              ( 1 , ig ) , g ( 1 , 1 , ig ) , ag ( 1 , ig ) , q1 , lidim , 
     .              lidim , wk_zv , hw2 ( 1 , 1 , is , 1 , js ) )

            endif

            do  18  j  = 1, ldimb
            do  18  i  = 1, ldima
              gfbz(i1,i2,i3,i,is,j,js) =
     .                dcmplx(hw2(i,1,is,j,js),hw2(i,2,is,j,js))
   18       continue
          endif

C         print 357,iq,i1,i2,i3,gfbz(i1,i2,i3,2,1,2,1)
C 357     format(' hamfb3: iq=',i4,' filling i1,i2,i3=',3i4,2f12.5)

   14   continue
   12 continue

      if (allocated(wk_zv)) deallocate(wk_zv)


C     Restore hq to complex*16 storage mode
      if (kcplx .eq. 1) then
        call ztoyy(hq,lidimx,lidimx,lidimx,lidimx,2,1)
      endif

      call tcx('hamfb3')
C     call yprm('gf',3,gfbz,0,k1*k2*k3,k1*k2*k3,ldima*ldimb)


c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$       open(3013,file="sfz_gfbz")
c$$$       do j1=1,ldima
c$$$       do j2=1,ldimb
c$$$        do ik1=1,k1
c$$$        do ik2=1,k2
c$$$        do ik3=1,k3
c$$$          write(3013,"(6i4,2d13.5)")ik1,ik2,ik3,j1,j2,isp,gfbz(ik1,ik2,ik3,j1,1,j2,1)
c$$$c           write(3012,"(6i4,2d13.5)")ik1,ik2,ik3,j1,j2,isp,sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$        enddo
c$$$        enddo
c$$$        enddo
c$$$       enddo
c$$$       enddo
c$$$       close(3013)
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc

      end

