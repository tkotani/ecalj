      subroutine psi2b_v3(nctot,ncc,nt0,ntp0,iclass,phase,
     i     cphik,               !rbk,cbk,rhbk,chbk,
     i     cphikq,              !rbkq,cbkq,rhbkq,chbkq,
     i     ppb,                 !pdb,dpb,ddb,
     i     nlnmv,nlnmc,mdim,imdim,iatomp,
     d     mdimx,nlmto,nbloch,nlnmx,
     d     natom,nclass,
     o     zpsi2b)
c     originaly 92.03.17 by Ferdi.
c     takao modified at Apr 2002(v2) Feb2006(v3).
c     
c     calculates <psi(k',t') | psi(k,t) B(R,i)>
c     for all R
c     psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
c     B(R,i)   = Bloch orthonormal product basis for atom R
c     psi(k,t) is stored after nctot
c     
c     nt0        = no. t
c     ntp0       = no. t'
c     coskt,sinkt= exp(ik.T)
c     cphik  b(k)
c     cphikq b(k')
c     
c     ppb        = <phi(RLn) phi(RL'n') B(R,i)>
c     
c     nlnmv      = number of l,n,m for valence
c     nlnmc      = number of n,l,m for core states
c     mdim       = number of optimal product basis functions for each class
c     nctot      = total no. allowed core states
c     nbloch     = total no. optimal product basis
c     nlnmx      = maximum number of l,n,m
c     
c     zpsi2b     =  the matrix elements
C----------------------------------------------------------
      implicit none
      integer(4):: nlmto,nctot,ncc,nt0,ntp0,nlnmx,natom,nbloch,
     &     ia,ic,nc,nv,nc1,ias,iap,icp,i,mdimx,nclass,itp,jp,ib
      complex(8):: cphik(nlmto,*),cphikq(nlmto,ntp0),phase(natom)
      complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc)
      integer(4),allocatable::iasx(:)
      integer(4) :: nzwork1,nzwork2,
     i     nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),
     i     imdim(natom),iatomp(natom)
      real(8)::  ppb(nlnmx,nlnmx,mdimx,nclass)
      complex(8),allocatable :: zz(:,:), zwork(:,:),zppb(:,:)
      complex(8) :: alpha,beta
      logical :: debug=.false.
      if(debug) write(*,*) "============== DEBUG: PASS psi2b_v3 ============ "
!! zwork is used as work array for ppb(:,:,i,ic) and for zpsi2b(ib,:,:).
      allocate( zppb(nlnmx,nlnmx),zz(nlnmx,ntp0), zwork(nt0,ntp0) ) !Work areas
      beta=0d0   
      alpha=1d0
      if(mdimx /= maxval(mdim) ) stop 'psi2b_v3: wrong mdimx'
      if( sum(mdim(iclass(1:natom)))/= nbloch ) stop 'psi2b_v3: wrong nbloch'
      allocate(iasx(natom))
      ias = 1
      do ia = 1,natom
         iasx(ia) = ias
         ias = ias + nlnmv(iclass(ia))
      enddo
      if(ias-1/=nlmto) stop ' psi2b_v3:sum(nlnmv)/= nlmto'
!! sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)>
!! for a given i, for all L'n' and t'
!! bkq is complex but <... > is real
      do ia = 1,natom
         ic   = iclass(ia) 
         nc   = nlnmc(ic)
         nv   = nlnmv(ic)
         nc1  = nc + 1
         if (nc+ nlnmv(ic) > nlnmx) stop 'psi2b_v3: nlnmx exceeded'
         ias  = iasx(ia)
         iap  = iatomp(ia)
         icp  = iclass(iap)
         do i = 1,mdim(icp)     ! loop over optimal product basis
            ib = imdim(iap)-1+i
            zppb(1:nv,1:nv) = ppb(nc+1:nc+nv,nc+1:nc+nv,i,icp)
            call zgemm('T','N',nv,ntp0,nv, alpha, zppb,nlnmx, cphikq(ias,1), nlmto,  beta,
     &           zz,  nlnmx )
            do itp = 1,ntp0
               do jp = 1,nv
                  zz(jp,itp)= dconjg(zz(jp,itp) )
               enddo
            enddo
            call zgemm('T','N', nt0,ntp0,nv, phase(ia), cphik(ias,1),nlmto, zz,nlnmx, beta,
     &           zwork, nt0 )
            zpsi2b(ib,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)=zwork
         end do                 !end of optimal product basis-loop
      end do                    !end of atom-loop
      deallocate(zz,zppb,zwork,iasx)
      if(debug) write(*,*) "=== CHECK: psi2b ", sum(zpsi2b)
      end subroutine psi2b_v3

c------------------------------------------------------------------------------------
      subroutine psicb_v3 (nctot,ncc,nt0,ntp0,  iclass, phase, 
     i     cphik,
     i     cphikq, 
     i     ppb,
     i     nlnmv,nlnmc,mdim,
     i     imdim,iatomp,
     i     mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i     icore,ncore,nl,nnc,  
     o     zpsi2b)
C     - Calculates <psi (k+q,t) |core(k,t) B(R,ibloch)> ---------------------
C     also   <core(k+q,t) |psi (k,t) B(R,ibloch)>
C     r B(R,i)   = Mixed basis.
C     r core(k,t)= core states
C     r <psi(k+q,t') | core(k,t) B(R,i)> = S[RLn]  cphik(RLn,k+q,t')^*  * ppb
C     r
C     i nt0        = no. k states
C     i ntp0       = no. kq states
C     i cphik coeefficients of MT part or valence eigenfunction.
C     i icore      = index for core states
C     i ncore      = no. core states in each class
C     i ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>
C     i nlnmv      = number of l,n,m for valence
C     i nlnmc      = number of l,n,m for core states
C     i mdim       = number of optimal product basis functions for each class
C     i nbloch     = total no. optimal product basis
C     i nlnmx      = maximum number of l,n,m
C     o zpsi2b     =  the matrix elements
C     r coskt,sinkt= exp(ik.T)
C     h           Ferdi 92.03.17       : original for ASA
C     h           takao at Apr2002(v2) : mixed basis version 2
C     h                    Feb2006(v3) : Add case with ncc/=0
C     r
C-------------------------------------------------------------------------
      implicit none
      integer(4):: nl,nlmto,nbloch,nclass,natom,nctot,ncc,nt0,ntp0,mdimx,nlnmx,nnc,
     &     ib,ias,ics,ia,ic,nc,nv,nc1,iap,icp,i,itp,it,icr
      integer(4)::  icore(nl*nl*nnc,nclass),ncore(nclass),
     i     nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom)
     i     ,imdim(natom),iatomp(natom)
      real(8)::   ppb(nlnmx,nlnmx,mdimx,nclass) !  coskt(natom),sinkt(natom)
      complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc),
     &     cphikq(nlmto,*), cphik(nlmto,*),phase(natom)
      integer(4):: verbose,ixx

      write(*,*) "============== DEBUG: PASS psicb_v3 ============ "

      if(ncc/=0 .and. ncc/=nctot) then
         stop "psicb_v3: ncc/=0 and ncc/=ncctot"
      endif

      if(sum(ncore(iclass(1:natom)))  /= nctot) 
     &     stop "psicb_v3:sum(ncore) wrong"

      zpsi2b = 0d0
      ib         = 0
      ias        = 1
      ics        = 0
      do      ia = 1,natom
         ic         = iclass(ia)
         nc         = nlnmc(ic)
         nv         = nlnmv(ic)
         nc1        = nc + 1
         iap        = iatomp(ia)
         icp        = iclass(iap)
         ib         = imdim(iap)-1
         do       i = 1,mdim(icp)
            ib         = ib + 1

            do     itp = 1,ntp0
               do      it = 1,ncore(ic)
                  icr      = icore(it,ic) ! write(6,*),it,ic,icore(it,ic)
                  zpsi2b(ib,ics+it,ncc+itp) = phase(ia)* 
     &                 dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))
               enddo            !end of t'(unoccupied)-loop
            enddo               !end of t(occupied)-loop

            if(verbose()>90) write(6,*)"psicb_v3: ia i=",ia,i
            if(ncc==0) cycle

            do     itp = 1,ncore(ic)
               do      it = 1,nt0
                  icr      = icore(itp,ic) ! write(6,*),it,ic,icore(it,ic)
                  if(verbose()>90) write(6,*)"psicb_v3: aaa itp it ic icr=",itp,it,ic,icr
cccccccccccccccccccccccccc
                  if(.false.) then
                     do ixx = 0,nv-1
                        write(6,"(5i5, 5d13.6)") ixx, nc1,icr,i,icp,cphik(ias+ixx,itp)
                        write(6,"(i5,33x,5d13.6)") ixx, ppb(nc1+ixx,icr,i,icp)
                     enddo
                  endif
ccccccccccccccccccccccccc
                  zpsi2b(ib,nctot+it,ics+itp) = dconjg(phase(ia))* 
     &                 sum(cphik(ias:ias+nv-1,it)*ppb(nc1:nc+nv,icr,i,icp))

                  if(verbose()>90) write(6,*)"psicb_v3: bbb itp it=",itp,it
               enddo            !end of t'(unoccupied)-loop
            enddo               !end of t(occupied)-loop
            if(verbose()>90) write(6,*)"psicb_v3: xxx ia i=",ia,i
         enddo                  !end of optimal product basis-loop
c     end of atom-loop
         ias   = ias + nlnmv(ic)
         ics   = ics + ncore(ic)
      enddo

      write(*,*) "=== CHECK: psicb ", sum(zpsi2b)

      end subroutine psicb_v3
