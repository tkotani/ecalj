      subroutine addbkgsm(smrho,k1,k2,k3,nsp,qbg,vol,fac)

C- add uniform background to smooth charge density
Ci Inputs: smrho: smooth density
Ci         qbg: background charge
Ci         k1,k2,k3, dimensions of smrho mesh
Ci         nsp: number of spins
Ci         vol: volume per cell
Ci         fac: fac*qbg/vol//nsp is added to
C-------------------------------------------
C     implicit none
      integer j1,j2,j3,k1,k2,k3,nsp,is
      double precision qbg,rhobg,vol,fac
      double complex smrho(k1,k2,k3,2)
      rhobg=qbg/vol
      do is=1,nsp
        do j1=1,k1
          do j2=1,k2
            do j3=1,k3
              smrho(j1,j2,j3,is)=smrho(j1,j2,j3,is)+rhobg*fac/nsp
            enddo
          enddo
        enddo
      enddo
      return
      end subroutine addbkgsm

      subroutine adbkql ( sv_p_orhoat , nbas , nsp , qbg , vol , fac 
     . , sspec , ssite )


      use m_struc_def  !Cgetarg

C- Add uniform bkg charge density to local smooth rho
Ci orhoat: pointers to local density in spheres
Ci nbas: number of atoms in basis
Ci qbg: background charge
Ci sspec: species structure
Ci nsp: spins
Ci vol: vol of cell
Ci fac: fac * backg density is added
Cu Updates
Cu   01 Jul 05 Zero-radius sites skipped over
C----------------------------------------
C     implicit none
      integer nrmx,nlmx,nlml,lmxl,nbas
      parameter (nrmx=1501,nlmx=64)
      integer:: nsp
      type(s_rv1) :: sv_p_orhoat(3,nbas)

      real(8):: qbg , fac
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)

      integer ib,nr,is
      double precision rhobkg,vol,a,rmt,rofi(nrmx)
      rhobkg = fac*qbg/vol
      do  ib = 1, nbas

        is=ssite(ib)%spec


        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt


        lmxl=sspec(is)%lmxl

        if (lmxl .eq. -1) goto 10
        nlml=(lmxl+1)**2
C      nrml=nr*nlml
        call rxx(nr .gt. nrmx,  'addbkgloc: increase nrmx')
        call rxx(nlml .gt. nlmx,'addbkgloc: increase nlmx')
        call radmsh(rmt,a,nr,rofi)
        call addbkgl ( sv_p_orhoat( 1 , ib )%v , sv_p_orhoat( 2 , ib )%v 
     .   , rhobkg , nr , nsp , rofi , nlml )

   10   continue
      enddo
      end subroutine adbkql

      subroutine addbkgl(rho1,rho2,rhobkg,nr,nsp,rofi,nlml)

C adds uniform background to local smooth density at this site for
C l=0 component (ilm=1)
C for each spin
C     implicit none
      integer nsp,is,nr,nlml,i
      double precision rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rofi(nr)
      double precision rhobkg,pi,srfpi
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4*pi)
C     y0 = 1d0/srfpi
      do is = 1, nsp
        do i = 1, nr
          rho1(i,1,is) = rho1(i,1,is)+srfpi*rofi(i)**2*rhobkg/nsp
          rho2(i,1,is) = rho2(i,1,is)+srfpi*rofi(i)**2*rhobkg/nsp
        enddo
      enddo
      end subroutine addbkgl




