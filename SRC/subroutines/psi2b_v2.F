      subroutine psi2b_v2(nt0,ntp0,iclass,coskt,sinkt,
     i                   cphik, !rbk,cbk,rhbk,chbk,
     i                   cphikq,!rbkq,cbkq,rhbkq,chbkq,
     i                   ppb, !pdb,dpb,ddb,
     i                   nlnmv,nlnmc,mdim,nctot,imdim,iatomp,
     d                   mdimx,nlmto,nbloch,nlnmx,noccxv,nt,ntp,
     d                   natom,nclass,
     o zpsi2b)
c originaly 92.03.17 by Ferdi.
c takao modified at Apr 2002
c calculates <psi(k',t') | psi(k,t) B(R,i)>
c for all R
c psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
c B(R,i)   = Bloch orthonormal product basis for atom R
c psi(k,t) is stored after nctot
c
c nt0        = no. t
c ntp0       = no. t'
c coskt,sinkt= exp(ik.T)
c cphik b(k)
c cphikq b(k')
c
c ppb        = <phi(RLn) phi(RL'n') B(R,i)>
c
c ddb        = <phidot(RLn) phidot(RL'n') B(R,i)>, s. ppbl.f
c nlnmv      = number of l,n,m for valence
c nlnmc      = number of n,l,m for core states
c mdim       = number of optimal product basis functions for each class
c nctot      = total no. allowed core states
c nbloch     = total no. optimal product basis
c nlnmx      = maximum number of l,n,m
c noccxv     = max. no. occupied valence states
c nt         = maximum number of occupied states
c ntp        = maximum number of unoccupied states
c
c zpsi2b     =  the matrix elements
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
      integer::nt0,ntp0,natom,nlmto,nclass,nlnmc,mdim,nctot
      complex(8):: cphik(nlmto,noccxv),cphikq(nlmto,ntp0)
     &            ,zpsi2b(nbloch,nt,ntp),phase
      dimension 
     i          ppb(nlnmx,nlnmx,mdimx,nclass),
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),
     i          coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)

      integer(4),allocatable::iasx(:)
      integer(4) :: nzppb1,nzppb2
      complex(8),allocatable :: zz(:,:), zppb(:,:)
      complex(8) :: alpha,beta

c zppb is used as work array for ppb(:,:,i,ic) and for zpsi2b(ib,:,:).
      nzppb1=max(nt0,nlnmx)
      nzppb2=max(ntp0,nlnmx)
      allocate( zz(nlnmx,ntp) )
      allocate( zppb(nzppb1,nzppb2) )
      beta=0d0  ; alpha=1d0
c      allocate(rr(nlnmx,ntp),cc(nlnmx,ntp))

c check dimensions
      if (ntp0 .gt. ntp) call rx( 'psi2bc: ntp exceeded')
      if (mdimx .ne. maxim(mdim,nclass)) call rx( 'psi2bc: wrong mdimx')
      if (nctot+nt0 .gt. nt) call rx( 'psi2bc: nt exceeded')
      if (nt0 .gt. noccxv) call rx( 'psi2bc: noccxv exceeded')
      if ( sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v2: wrong nbloch')
      allocate(iasx(natom))
      ias = 1
      do ia = 1,natom
        iasx(ia) = ias
        ias = ias + nlnmv(iclass(ia))
      enddo
Cstop2rx 2013.08.09 kino      if(ias-1/=nlmto) stop ' psi2b_v2:sum(nlnmv)/= nlmto'
      if(ias-1/=nlmto) call rx( ' psi2b_v2:sum(nlnmv)/= nlmto')

c loop over atoms
      do  ia = 1,natom
        ic   = iclass(ia)
        nc   = nlnmc(ic)
        nv   = nlnmv(ic)
        nc1  = nc + 1
Cstop2rx 2013.08.09 kino        if (nc+ nlnmv(ic) .gt. nlnmx) stop 'psi2b_v2: nlnmx exceeded'
        if (nc+ nlnmv(ic) .gt. nlnmx) call rx( 'psi2b_v2: nlnmx exceeded')
        phase= dcmplx(coskt(ia),sinkt(ia))
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
     &   alpha, zppb,nzppb1, cphikq(ias,1), nlmto,  beta,
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
     &  phase, cphik(ias,1),nlmto, zz,nlnmx, beta,
     &  zppb, nzppb1 )
          ib = imdim(iap)-1+i
          zpsi2b(ib,nctot+1:nctot+nt0,1:ntp0)=zppb(1:nt0,1:ntp0)
c------------------------------------------------------
        end do !end of optimal product basis-loop
      end do !end of atom-loop
c      deallocate(rr,cc,iasx)
      deallocate(zz,zppb,iasx)
      end

c------------------------------------------------------------------------------------
      subroutine psicb_v2 (icore,ncore,ntp0,iclass,coskt,sinkt,
     i                   cphikq, !rbkq,cbkq,rhbkq,chbkq,
     i                   ppb,!pdb,dpb,ddb
     i                   nlnmv,nlnmc,mdim,
     i                   imdim,iatomp,
     d                   mdimx,nlmto,nbloch,nlnmx,nt,ntp,natom,nclass,
     d                   nl,nnc,
     o zpsi2b)! rpsi2b,cpsi2b)
c written by Ferdi  92.03.17
c takao modified at Apr 2002
c
c calculates <psi(k+q,t') | core(k,t) B(R,i)>
c for all R
c psi(k,t) = S[RLn] b(RLn,k,t)*X(RLn,k)
c          = S[RLn] b(RLn,k,t)*Phi(RLn,k) + hb(RLn,k,t)*Phidot(RLn,k)
c core(k,t)= core states
c B(R,i)   = Bloch orthonormal product basis for atom R
c
c <psi(k+q,t') | core(k,t) B(R,i)>
c = S[RLn]  b(RLn,k+q,t')^* <Phi(RLn)    |core(k,t) B(R,i)>
c + S[RLn] hb(RLn,k+q,t')^* <Phidot(RLn) |core(k,t) B(R,i)>

c ntp0       = no. unoccupied states
c coskt,sinkt= exp(ik.T)
c cphikq  = real and imaginary part of b(k+q).
c            coefficients of eigenfunctions for argumentationwaves in each MT
c
c icore      = index for core states
c ncore      = no. core states in each class
c ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>
c
c nlnmv      = number of l,n,m for valence
c nlnmc      = number of l,n,m for core states
c mdim       = number of optimal product basis functions for each class
c nbloch     = total no. optimal product basis
c nlnmx      = maximum number of l,n,m
c nt         = maximum number of occupied states
c ntp        = maximum number of unoccupied states
c
c zpsi2b     =  the matrix elements
c
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
      complex(8):: cphikq(nlmto,ntp0),zpsi2b(nbloch,nt,ntp),phase
      dimension
c                rbkq(nlmto,ntp0),cbkq(nlmto,ntp0),
c     i          rhbkq(nlmto,ntp0),chbkq(nlmto,ntp0),
     i          icore(nl*nl*nnc,nclass),ncore(nclass),
     i          ppb(nlnmx,nlnmx,mdimx,nclass),
c     i          pdb(nlnmx,nlnmx,mdimx,nclass),
c     i          dpb(nlnmx,nlnmx,mdimx,nclass),
c     i          ddb(nlnmx,nlnmx,mdimx,nclass),
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),
     i          coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)
c      dimension rpsi2b(nbloch,nt,ntp),
c     o          cpsi2b(nbloch,nt,ntp)
c initialise matrix elements
c      call dinit   (rpsi2b,nbloch*nt*ntp)
c      call dinit   (cpsi2b,nbloch*nt*ntp)

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
        phase  =  dcmplx(coskt(ia),sinkt(ia))
c loop over optimal product basis
        iap        = iatomp(ia)
        icp        = iclass(iap)
        ib         = imdim(iap)-1
        do       i = 1,mdim(icp)
          ib         = ib + 1

c S[Ln] bkq(Ln,t')^(*) * <Phi(Ln) core(L'n') B(i)> etc.
c for a given i, for all L'n' and t'
c bkq is complex but < > is real
          do     itp = 1,ntp0
            do      it = 1,ncore(ic)
              icr        = icore(it,ic)

cccccccccccccc
c      write(6,*),it,ic,icore(it,ic)
cccccccccccccc
c
c      rs1        = vdv ( rbkq(ias,itp),ppb(nc1,icr,i,icp),nv)
c     .           + vdv (rhbkq(ias,itp),dpb(nc1,icr,i,icp),nv)
c      cs1        = vdv ( cbkq(ias,itp),ppb(nc1,icr,i,icp),nv)
c     .           + vdv (chbkq(ias,itp),dpb(nc1,icr,i,icp),nv)
c
c      rpsi2b(ib,ics+it,itp) = rs1*coskt(ia) + cs1*sinkt(ia)
c      cpsi2b(ib,ics+it,itp) = rs1*sinkt(ia) - cs1*coskt(ia)
cccccccccccccccccccccccccccccccc
c      if(abs(rs1)>1d8.or.abs(cs1)>1d8) then
c        print *,'  psicb2*:'
c        print *, nc1,icr,i,icp, ppb(nc1,icr,i,icp)
c        print *,    dpb(nc1,icr,i,icp)
c       stop
c      endif
cccccccccccccccccccccccccccccccc

              zpsi2b(ib,ics+it,itp) = phase* 
     & dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))

c end of t'(unoccupied)-loop
            end do
c end of t(occupied)-loop
          end do

c end of optimal product basis-loop
        end do

c end of atom-loop
        ias        = ias + nlnmv(ic)
        ics        = ics + ncore(ic)
      end do

      return
      end

















