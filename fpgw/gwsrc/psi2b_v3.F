      subroutine psi2b_v3(nctot,ncc,nt0,ntp0,iclass,phase,
     i                   cphik, !rbk,cbk,rhbk,chbk,
     i                   cphikq,!rbkq,cbkq,rhbkq,chbkq,
     i                   ppb, !pdb,dpb,ddb,
     i                   nlnmv,nlnmc,mdim,imdim,iatomp,
     d                   mdimx,nlmto,nbloch,nlnmx,
     d                   natom,nclass,
     o zpsi2b)
c originaly 92.03.17 by Ferdi.
c takao modified at Apr 2002(v2) Feb2006(v3).
c
c calculates <psi(k',t') | psi(k,t) B(R,i)>
c for all R
c psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
c B(R,i)   = Bloch orthonormal product basis for atom R
c psi(k,t) is stored after nctot
c
c nt0        = no. t
c ntp0       = no. t'
c coskt,sinkt= exp(ik.T)
c cphik  b(k)
c cphikq b(k')
c
c ppb        = <phi(RLn) phi(RL'n') B(R,i)>
c
c nlnmv      = number of l,n,m for valence
c nlnmc      = number of n,l,m for core states
c mdim       = number of optimal product basis functions for each class
c nctot      = total no. allowed core states
c nbloch     = total no. optimal product basis
c nlnmx      = maximum number of l,n,m
c
c zpsi2b     =  the matrix elements
C----------------------------------------------------------
      implicit none
      integer(4):: nlmto,nctot,ncc,nt0,ntp0,nlnmx,natom,nbloch,
     & ia,ic,nc,nv,nc1,ias,iap,icp,i,mdimx,nclass,itp,jp,ib
      complex(8):: cphik(nlmto,*),cphikq(nlmto,ntp0),phase(natom)
      complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc)
      integer(4),allocatable::iasx(:)
      integer(4) :: nzwork1,nzwork2,
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),
     i          imdim(natom),iatomp(natom)
      real(8)::  ppb(nlnmx,nlnmx,mdimx,nclass)
      complex(8),allocatable :: zz(:,:), zwork(:,:),zppb(:,:)
      complex(8) :: alpha,beta

c zwork is used as work array for ppb(:,:,i,ic) and for zpsi2b(ib,:,:).
      allocate( zppb(nlnmx,nlnmx),zz(nlnmx,ntp0), zwork(nt0,ntp0) ) !Work areas
      beta=0d0  ; alpha=1d0
c      allocate(rr(nlnmx,ntp),cc(nlnmx,ntp))

c check dimensions
      if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
      if( sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
      allocate(iasx(natom))
      ias = 1
      do ia = 1,natom
        iasx(ia) = ias
        ias = ias + nlnmv(iclass(ia))
      enddo
Cstop2rx 2013.08.09 kino      if(ias-1/=nlmto) stop ' psi2b_v3:sum(nlnmv)/= nlmto'
      if(ias-1/=nlmto) call rx( ' psi2b_v3:sum(nlnmv)/= nlmto')

c loop over atoms
      do  ia = 1,natom
        ic   = iclass(ia)
        nc   = nlnmc(ic)
        nv   = nlnmv(ic)
        nc1  = nc + 1
Cstop2rx 2013.08.09 kino        if (nc+ nlnmv(ic) > nlnmx) stop 'psi2b_v3: nlnmx exceeded'
        if (nc+ nlnmv(ic) > nlnmx) call rx( 'psi2b_v3: nlnmx exceeded')
c        phase= dcmplx(coskt(ia),sinkt(ia))
        ias  = iasx(ia)
        iap  = iatomp(ia)
        icp  = iclass(iap)
        do   i = 1,mdim(icp) ! loop over optimal product basis
c---------------------------------------------------
cc sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)>
cc for a given i, for all L'n' and t'
cc bkq is complex but < > is real
c1      do     itp = 1,ntp0
c1      do      jp = 1,nlnmv(ic)
c1      zz (jp,itp)=dconjg(
c1     &    sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,nc+jp,i,icp)) )
c1      end do
c1      end do
c
c2      zz(1:nv,1:ntp0) =dconjg(
c2     & matmul(  transpose(ppb(nc1:nc+nv,nc1:nc+nv,i,icp))
c2     &         ,cphikq(ias:ias+nv-1,1:ntp0)) )
c
c3        call dgemm('T','N',nv,ntp0,nv,
c3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
c3     &          dreal(cphikq(ias:ias+nv-1,1:ntp0)), nv,
c3     &          0d0,
c3     &   rr, nlnmx )
c3        call dgemm('T','N',nv,ntp0,nv,
c3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
c3     &          dimag(cphikq(ias:ias+nv-1,1:ntp0)), nv,
c3     &          0d0,
c3     &   cc, nlnmx )
          zppb(1:nv,1:nv) = ppb(nc+1:nc+nv,nc+1:nc+nv,i,icp)
          call zgemm('T','N',nv,ntp0,nv,
     &   alpha, zppb,nlnmx, cphikq(ias,1), nlmto,  beta,
     &   zz,  nlnmx )
          do itp = 1,ntp0
            do jp = 1,nv
              zz(jp,itp)= dconjg(zz(jp,itp) )
            enddo
          enddo
c----------------------------------------------------
c <psi(k+q,t') | psi(k,t) B(i)>
c1      do      it = 1,nt0
c1      do     itp = 1,ntp0
c1       zpsi2b(ib,nctot+it,itp)=
c1     &   phase * sum( zz(1:nv,itp)*cphik(ias:ias+nv-1,it) )
c1c end of t'(unoccupied)-loop
c1      end do
c1c end of t(occupied)-loop
c1      end do
c3        call zgemm('T','N',nt0,ntp0,nv,
c3     &   phase, cphik(ias:ias+nv-1,1:nt0),  nv,
c3     &          dcmplx(rr(1:nv,1:ntp0),-cc(1:nv,1:ntp0)),  nv,
c3     &          0d0,
c3    &   zpsi2b(imdim(iap)-1+i,nctot+1:nctot+nt0,1:ntp0), nt0)
          call zgemm('T','N', nt0,ntp0,nv,
     &   phase(ia), cphik(ias,1),nlmto, zz,nlnmx, beta,
     &   zwork, nt0 )
          ib = imdim(iap)-1+i
          zpsi2b(ib,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)=zwork
c------------------------------------------------------
        end do !end of optimal product basis-loop
      end do !end of atom-loop
c      deallocate(rr,cc,iasx)
      deallocate(zz,zppb,zwork,iasx)
      end


c------------------------------------------------------------------------------------
      subroutine psicb_v3 (nctot,ncc,nt0,ntp0,  iclass, phase, 
     i                   cphik,
     i                   cphikq, 
     i                   ppb,
     i                   nlnmv,nlnmc,mdim,
     i                   imdim,iatomp,
     i                   mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     i                   icore,ncore,nl,nnc,  
     o zpsi2b)
C- Calculates <psi (k+q,t) |core(k,t) B(R,ibloch)> ---------------------
C      also   <core(k+q,t) |psi (k,t) B(R,ibloch)>
Cr B(R,i)   = Mixed basis.
Cr core(k,t)= core states
Cr <psi(k+q,t') | core(k,t) B(R,i)> = S[RLn]  cphik(RLn,k+q,t')^*  * ppb
Cr
Ci nt0        = no. k states
Ci ntp0       = no. kq states
Ci cphik coeefficients of MT part or valence eigenfunction.
Ci icore      = index for core states
Ci ncore      = no. core states in each class
Ci ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>
Ci nlnmv      = number of l,n,m for valence
Ci nlnmc      = number of l,n,m for core states
Ci mdim       = number of optimal product basis functions for each class
Ci nbloch     = total no. optimal product basis
Ci nlnmx      = maximum number of l,n,m
Co zpsi2b     =  the matrix elements
Cr coskt,sinkt= exp(ik.T)
Ch           Ferdi 92.03.17       : original for ASA
Ch           takao at Apr2002(v2) : mixed basis version 2
Ch                    Feb2006(v3) : Add case with ncc/=0
Cr
C-------------------------------------------------------------------------
      implicit none
      integer(4):: nl,nlmto,nbloch,nclass,natom,nctot,ncc,nt0,ntp0,mdimx,nlnmx,nnc,
     & ib,ias,ics,ia,ic,nc,nv,nc1,iap,icp,i,itp,it,icr
      integer(4)::  icore(nl*nl*nnc,nclass),ncore(nclass),
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom)
     i          ,imdim(natom),iatomp(natom)
      real(8)::   ppb(nlnmx,nlnmx,mdimx,nclass) !  coskt(natom),sinkt(natom)
      complex(8):: zpsi2b(nbloch,nt0+nctot,ntp0+ncc),
     &             cphikq(nlmto,*), cphik(nlmto,*),phase(natom)
      integer(4):: verbose,ixx
      if(ncc/=0 .and. ncc/=nctot) then
Cstop2rx 2013.08.09 kino        stop "psicb_v3: ncc/=0 and ncc/=ncctot"
        call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
      endif
      if(sum(ncore(iclass(1:natom)))  /= nctot) 
Cstop2rx 2013.08.09 kino     &   stop "psicb_v3:sum(ncore) wrong"
     &   call rx( "psicb_v3:sum(ncore) wrong")
c      write(6,*)" psicb_v3: ncc nctot=",ncc,nctot,sum(ncore(iclass(1:natom)))
c      stop 'qqqqqqqqqqqqqq'
      zpsi2b = 0d0
c loop over atoms
      ib         = 0
      ias        = 1
      ics        = 0
      do      ia = 1,natom
        ic         = iclass(ia)
        nc         = nlnmc(ic)
        nv         = nlnmv(ic)
        nc1        = nc + 1
c       phase  =  dcmplx(coskt(ia),sinkt(ia))
c loop over optimal product basis
        iap        = iatomp(ia)
        icp        = iclass(iap)
        ib         = imdim(iap)-1
        do       i = 1,mdim(icp)
          ib         = ib + 1
          do     itp = 1,ntp0
            do      it = 1,ncore(ic)
              icr      = icore(it,ic) ! write(6,*),it,ic,icore(it,ic)
              zpsi2b(ib,ics+it,ncc+itp) = phase(ia)* 
     &    dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))
            enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop

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
     &    sum(cphik(ias:ias+nv-1,it)*ppb(nc1:nc+nv,icr,i,icp))
              if(verbose()>90) write(6,*)"psicb_v3: bbb itp it=",itp,it
            enddo  !end of t'(unoccupied)-loop
          enddo  !end of t(occupied)-loop
          if(verbose()>90) write(6,*)"psicb_v3: xxx ia i=",ia,i
        enddo !end of optimal product basis-loop
c end of atom-loop
        ias   = ias + nlnmv(ic)
        ics   = ics + ncore(ic)
      enddo
      end
