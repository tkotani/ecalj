      subroutine zhevx(n,lh,h,s,lov,lx,nmx,emx,nev,wk,linv,e,lz,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   n:    order of h and s
Ci   lh:   leading dimension of h and s
Ci   h:    hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix, (used only if lov is true)
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci         (not used if LAPACK zhegv is invoked)
Ci   wk:   work array of length at least 11n
Ci         NB: If LAPACK version is used, and eigenvectors are sought
Ci         wk should be dimensioned (n*nmx*2)
Ci   lov:  0 no overlap matrix
Ci         1 overlap matrix, return evecs of nonorthogonal H
Ci   lx:   if T, calls routines to exploit unit stride lengths (risc)
Ci         Not used if LAPACK zhegv is invoked.
Ci   linv: if T, using inverse iteration
Ci         Not used if LAPACK zhegv is invoked.
Ci   lz:   leading dimension of z
Co Outputs:
Co   e:    eigenvalues
Co   nev:  number of eigenvectors found
Co   z:    eigenvectors (1..nev)  (declared as z(n,*)
Co   s:    has been decomposed into and LL+ decomposition.
Co         You can call zhev2 to scale a vector by L
Cr Remarks:
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cr   Aborts on exit
Cp Procedures used:
Cp   (lapack)  zhegv
Cp   (eispack) htribk, htridx, imtql2, tqlrat
Cu Updates
Cu   24 Feb 07 Bug fix when nmx=0
Cu   17 May 03 Adapted from zhev, intended to supersede zhev.
Cu   14 Aug 02 Added zheev when lov is F; new zhev2.
Cu   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
C ----------------------------------------------------------------
      implicit none
      logical linv,lx
      integer lov,lh,lz,nn
      integer n,nev,nmx
      double precision h(n,n),s(*),e(n),wk(*),emx
      integer ier,lwork
      character jobz
      complex(8)::z(n,*)
      call tcn('zhevx')
      if(lh/=n) call rx('zhevx: we assume lh=n now')
      if (nmx .le. 0) then
        jobz = 'N'
        lwork = 4*n
        if (lov .gt. 0) then
          call zhegv(1,jobz,'U',n,h,lh,s,lh,e,wk(1+3*n),lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zhegv cannot find all evals')
        else
          call zheev(jobz,'U',n,h,lh,e,wk(1+3*n),lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zheev cannot find all evals')
        endif
        nev = 0
      else
        jobz = 'V'
        lwork = n*nmx
        if (lov .gt. 0) then
          call zhegv(1,jobz,'U',n,h,lh,s,lh,e,z,lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zhegv cannot find all evals')
        else
          call zheev(jobz,'U',n,h,lh,e,z,lwork,wk(1),ier)
          call rxx(ier.ne.0,'zhevx: zheev cannot find all evals')
        endif
        nn=min(n,nmx)
        z(1:n,1:nn)=h(1:n,1:nn)
c        call zmcpy('N',h,lh,1,z,lz,1,n,min(n,nmx))
C       call zprm('evecs',2,z,lz,n,nmx)
        nev = min(n,nmx)
      endif
  100 call tcx('zhevx')
      end

      subroutine zhev(n,h,s,lov,lx,nmx,emx,nev,wk,linv,ltime,e,z)
C- Eigenvalues and/or some eigenvectors of a Hermitian matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   n:    dimension of h
Ci   h,n:  hermitian matrix, dimensioned h(n,n)
Ci   s:    hermitian overlap matrix, (used only if lov is true)
Ci   nmx:  maximum number of eigenvectors to be found
Ci   emx:  eigenvalue limit for eigenvectors to be found
Ci   wk:   work array of length at least 11n
Ci   lov:  if T, non-orthogonal
Ci   lx:   if T, calls routines to exploit unit stride lengths (risc)
Ci         Not used if LAPACK zhegv is invoked.
Ci   linv: if T, using inverse iteration
Ci         Not used if LAPACK zhegv is invoked.
Co Outputs:
Co   e:    eigenvalues
Co   nev:  number of eigenvectors found
Co   z:    eigenvectors (1..nev)  (declared as z(n,*)
Co   s:    has been decomposed into and LL+ decomposition.
Co         You can call zhev2 to scale a vector by L
Cr Remarks:
Cr   z must be at least of dimension z(n,n), even though nev<n.
Cr   h and s are destroyed on exit.
Cr   Aborts on exit
Cp Procedures used:
Cp   (lapack)  zhegv
Cp   (eispack) htribk, htridx, imtql2, tqlrat
Cu Updates
Cu   14 Aug 02 Added zheev when lov is F; new zhev2.
Cu   21 Jan 02 Added code to invoke LAPACK zhegv in place of diagno
C ----------------------------------------------------------------
      implicit none
      logical linv,lx
      integer n,nev,nmx,ltime
      double precision h(*),s(*),e(n),wk(*),z(*),emx
      logical lov
      integer ier,lwork
      character jobz
      call tcn('zhev')
      if (nmx .le. 0) then
        jobz = 'N'
        lwork = 4*n
        if (lov) then
          call zhegv(1,jobz,'U',n,h,n,s,n,e,wk(1+3*n),lwork,wk(1),ier)
        else
          call zheev(jobz,'U',n,h,n,e,wk(1+3*n),lwork,wk(1),ier)
        endif
        nev = 0
      else
        jobz = 'V'
        lwork = n*min(n,nmx)
        if (lov) then
          call zhegv(1,jobz,'U',n,h,n,s,n,e,z,lwork,wk(1),ier)
        else
          call zheev(jobz,'U',n,h,n,e,z,lwork,wk(1),ier)
        endif
        call zcopy(n*min(n,nmx),h,1,z,1)
        nev = min(n,nmx)
      endif
      call rxx(ier.ne.0,'zhev: zhegv cannot find all evals')
  100 call tcx('zhev')
      end

!!-----------------------------------
      subroutine zhev_tk2(n,h,s,nmx,nev, e,z)
!!== Eigenvalues and/or some eigenvectors of a Hermitian matrix (weighted for first nlmto basis).==
!! ----------------------------------------------------------------
!! Inputs:
!!c   nlmto:dimension of MTO space of 1:nlmto i respected when diagonalization.
!!     n:    dimension of h
!!   h,n:  hermitian matrix, dimensioned h(n,n)
!!   s:    hermitian overlap matrix, 
!!   nmx:  requested number of eigenvectors to be found (and eigenvalues). If nmx>n, nmx is taken to be n.
!!         if nmx=0, nev=n (see NOTE below).
!!   ipr :print switch
!!   ifig, savez,getz:dummy
!! Outputs:
!!   e:    eigenvalues
!!   nev:  number of eigenvectors (=nmx) or (=n if nmx=0)
!!   z:    eigenvectors (1..nev)  (declared as z(n,*)
!!   h and s are destroyed on exit.
!!   july2012takao
!! NOTE: this can be called in the loop of ikp,isp loop. Then data are appended to a file ifig.
!!   If nmx==0, all eigenvalues are returned but without eigenfunctions.
!! -----------------------------------------------------------------------
      implicit none
      integer n,nev,nmx,ltime,ngv,ncut
c      logical ipr 
      complex(8) :: h(n,n),s(n,n),z(n,*)
      complex(8),allocatable:: work(:)
      integer:: i,j,ni,ik,ik2,k
      real(8):: epsovl,epsx,e(n),emx,eee,fac,vldummy,vudummy,abstol
      real(8),allocatable:: rwork(:) 
      integer:: ier,lwork,ix,ifi,ifig,nmx0,ifail(n),nevx
      integer,allocatable:: iwork(:)
      character*1:: jobz
      logical ::nexist
      integer,save:: lworksave=0
      call tcn('zhev_tk2')
!! === diagonalize ===
      if(nmx==0) then
         jobz='N'
         nev=n
      else
        jobz = 'V'
        nev=nmx
      endif
      abstol=1d-10 ! OK?
      lwork=max(1,2*n,lworksave) !OK? efficient?
      allocate(work(lwork),rwork(7*n),iwork(5*n))
      call zhegvx(1,jobz,'I','U',n,h,n,s,n,vldummy,vudummy,1,nev,abstol,nevx,e,z,n,work,lwork,rwork,iwork,ifail,ier)
      lworksave= WORK(1)  !this is optimum lwork
c      print *,'nev nevx n=',nev,nevx,n
      call rxx(nev/=nevx,'zhev_tk2: nev /=nevx something wrong. ')
      call rxx(ier.ne.0, 'zhev_tk2: zhegvx cannot find all eigen.')
      deallocate(work,iwork,rwork)
      call tcx('zhev_tk2')
      end

!!=====================================================================================
      subroutine zhev_tk4(n,h,s,nmx,nev, e,z, epsovl)
!!== Eigenvalues and/or some eigenvectors of a Hermitian matrix (weighted for first nlmto basis).==
!! ----------------------------------------------------------------
!! Inputs:
!!   nlmto:dimension of MTO space of 1:nlmto i respected when diagonalization.
!!     n:    dimension of h
!!   h,n:  hermitian matrix, dimensioned h(n,n)
!!   s:    hermitian overlap matrix, 
!!   nmx:  requested number of eigenvectors to be found (and eigenvalues). If nmx>n, nmx is taken to be n.
!!         if nmx=0, nev=n (see NOTE below).
!!   ipr :print switch
!!   ifig,savez,getz: dummy
!!   epsovl: cutoff to remove Hlbert space.
!! Outputs:
!!   e:    eigenvalues
!!   nev:  number of eigenvectors (=(min(nm,nmx)) or (=nm if nmx=0)
!!   z:    eigenvectors (1..nev)  (declared as z(n,*)
!!   h and s are destroyed on exit.
!!   july2012takao
!! nm is the matrix dimension of the reduced space by epsovl
!!   If nmx==0, all eigenvalues are returned but without eigenfunctions.
!!
!! Essentially similar with zhevo
!! -----------------------------------------------------------------------
      implicit none
      integer n,nev,nmx,ltime,ngv,ncut
c      logical ipr 
      complex(8) :: h(n,n),s(n,n),z(n,nmx)
      complex(8),allocatable:: work(:)
      integer:: i,j,ni,ik,ik2,k
      real(8):: epsovl,epsx,e(n),emx,eee,fac,vldummy,vudummy,abstol
      real(8),allocatable:: rwork(:) 
      integer:: ier,lwork,ix,ifi,ifig,nmx0,ifail(n),nevx
      integer,allocatable:: iwork(:)
      character*1:: jobz
      logical ::nexist
      integer,save:: lworksave=0

      complex(8),allocatable::ii(:,:)
      real(8):: eo(n) 
      complex(8),allocatable:: omat(:,:),wk11(:),
     & zz(:,:),hh(:,:),hhm(:,:),znm(:,:)
      integer:: nevl,nm,nmout,nevout
      logical:: debug=.false.
      
      if(epsovl< 1d-14) then
         call zhev_tk2(n,h,s,nmx,nev, e,z) 
         return
      endif
      call tcn('zhev_tk4')
!! ... eigenvalue of ovarlap matrix
      allocate(omat(n,n))
      omat = s !reserved
      jobz = 'V'
      lwork = n*n
      allocate(wk11(lwork),rwork(max(1,3*n-2)))
      call zheev(jobz,'U',n,omat,n,eo,wk11,lwork,rwork,ier)
      deallocate(wk11,rwork)
      if(debug) then
        write(6,*)'zhev_tk4: ovlmat='
        do ix=0,n,5
          write(6,"(5(i5,d10.2))") (i,eo(i),i=ix+1,min(ix+5,n))
        enddo  
      endif
      do ix= 1,n
        if(eo(ix)>epsovl) then
          ni = ix               ! take i = ni...n
          exit
        endif
      enddo
      nm = n-ni+1            ! this is the dimension.
      nevl=nm
      allocate(zz(n,nm))     ! zz is the projection matrix  
      do ix=ni,n
        zz(:,ix-ni+1) = omat(:,ix)/sqrt(eo(ix))
      enddo
      if(debug) then
      write(6,*)' reduced by OVEPS: n--> nm=',nm
      endif
!! Hamiltonian  <zz|H|zz>
      allocate(hh(nm,nm),hhm(nm,n))
c This failed in ifort when hm >600 or so.-->maybe need ulimit -s unlimited.
c      hh = matmul(dconjg(transpose(zz)),matmul(h,zz))
c In anyway, blas will be better.
      call zgemm('C','N',nm,n,n,(1d0,0d0),zz,n,h,n,(0d0,0d0),hhm,nm)
      call zgemm('N','N',nm,nm,n,(1d0,0d0),hhm,nm,zz,n,(0d0,0d0),hh,nm)
      deallocate(hhm)
      if(nmx==0) then
         jobz='N'
         nev=nm
      else
        jobz = 'V'
        nev = min(nmx,nm) 
      endif
!! note nev: number of output eigenvalues (and eigenfunctions when jobz=V).
      abstol= 1d-10 ! OK?
      lwork = max(1,2*nm,lworksave) !OK? efficient?
      allocate(work(lwork),rwork(7*nm),iwork(5*nm),znm(nm,max(1,nev)))
      call zheevx(jobz,'I','U',nm,hh,nm,vldummy,vudummy,1,nev,abstol,nevout,e,znm,nm,work,lwork,rwork,iwork,ifail,ier)
      lworksave= WORK(1)  !this is optimum lwork right?
      call rxx(nev/=nevout,'zhev_tk4: nev /=nevout something wrong. ')
      call rxx(ier.ne.0, 'zhev_tk4: zheev for hh cause error.')
      deallocate(work,iwork,rwork)
      z=1d99
      do i=1,min(nmx,nm)
        do j=1,n
          z(j,i) = sum(zz(j,:)*znm(:,i)) !this is eigenfunction for original problem.
        enddo
      enddo
      deallocate(znm)
      deallocate(zz)
      if(.false.) then !! === diagonalize === (this part is in zhev_tk2), Kept here for debug purpose
        if(nmx==0) then
          jobz='N'
          nev=n
        else
          jobz = 'V'
          nev=nmx
        endif
        abstol=1d-10            ! OK?
        lwork=max(1,2*n,lworksave) !OK? efficient?
        allocate(work(lwork),rwork(7*n),iwork(5*n))
        call zhegvx(1,jobz,'I','U',n,h,n,s,n,vldummy,vudummy,1,nev,abstol,nevx,e,z,n,work,lwork,rwork,iwork,ifail,ier)
        lworksave= WORK(1)      !this is optimum lwork
        call rxx(nev/=nevx,'zhev_tk4: nev /=nevx something wrong. ')
        call rxx(ier.ne.0, 'zhev_tk4: zhegvx cannot find all eigen.')
        deallocate(work,iwork,rwork) 
      endif 
      call tcx('zhev_tk4')
      end
