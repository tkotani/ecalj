      subroutine calcpv2 (zw,rv,cv,
c     w                wpv,ovlp, !wdiag,iwdiag,
c     w                evec,eval,
     d                nbloch,
     o                rpv,cpv,rw,cw )
C- calculate v^{1/2}Pv^{1/2} for given P and v --------
Ci zw : P, the polarization function. <\tilde{M}_i | P | \tilde{M}_j >
Ci rv+ i cv: v, the Coulomb matrix    <M_i |v |M_j>
Ci nbloch  = dimension of vc,P
Co rw  + i cw: (1-Pv)(1-vP)
Co rpv +i cpv:  Pv
c--- History ---------
C tkotani modifyed from Miayke's calcpv Nov2004
ctm 01.07.24
c pv is (1-Pv)(1-vP), not v^{1/2}Pv^{1/2}
c
ctm 01.07.01
c 00.07.03
c---------------------------------------------------------
      implicit double precision (a-h,o-z)

      dimension   rv(nbloch,nbloch),cv(nbloch,nbloch),
     w            rpv(nbloch,nbloch),cpv(nbloch,nbloch),
     w            rw(nbloch,nbloch),cw(nbloch,nbloch)
c     w            wpv(nbloch,nbloch,2),ovlp(nbloch,nbloch,2),
c     w            wdiag(11*nbloch),iwdiag(nbloch),
c     w            evec(nbloch,nbloch,2),eval(nbloch)
      complex(8) zw(nbloch,nbloch)
      real(8),allocatable:: rp(:,:),cp(:,:),
     w            rw1(:,:),cw1(:,:),rw2(:,:),cw2(:,:)

      logical :: debug = .false.
      integer(4):: verbose
c-----------------------------------------
      if(verbose()>100) debug=.true.
      nbloch2    = nbloch*nbloch
      pi         = 4.d0*datan(1.d0)
      allocate ( rp(nbloch,nbloch), cp(nbloch,nbloch),
     &          rw1(nbloch,nbloch),cw1(nbloch,nbloch),
     &          rw2(nbloch,nbloch),cw2(nbloch,nbloch))

C... Polarization function is decompsed into rp + i cp
      rp = dreal(zw(1:nbloch,1:nbloch))
      cp = dimag(zw(1:nbloch,1:nbloch))

C... symmetrize the matrix
      call hermat  (rv,cv,nbloch)
      call hermat  (rp,cp,nbloch)

C... calculate Pv
      if (debug) write(*,*)'debug> mmulc vP in'
      call mmulc   (rp,cp,nbloch,
     i              rv,cv,nbloch,
     i              nbloch,nbloch,nbloch,nbloch,
     o              rpv,cpv )
      rw1       = - rpv
      cw1       = - cpv
      do i      = 1,nbloch
        rw1(i,i)  = 1d0 + rw1(i,i)
      end do

C... calculate vP
      if (debug) write(*,*)'debug> mmulc vP in'
      call mmulc   (rv,cv,nbloch,
     i              rp,cp,nbloch,
     i              nbloch,nbloch,nbloch,nbloch,
     o              rw2,cw2 )
      rw2       = - rw2
      cw2       = - cw2
      do i      = 1,nbloch
        rw2(i,i)  = 1d0 + rw2(i,i)
      end do

C... calculate (1-Pv)(1-vP)
      if (debug) write(*,*)'debug> mmulc vP in'
      call mmulc   (rw1,cw1,nbloch,
     i              rw2,cw2,nbloch,
     i              nbloch,nbloch,nbloch,nbloch,
     o              rw,cw )

c      if (debug) write(*,*)'debug> mmulc out'
c      if (debug) then
c        write(*,*)'*** P'
c        call chhermat (rp,cp,nbloch)
c        write(*,*)'*** v'
c        call chhermat (rv,cv,nbloch)
c        write(*,*)'*** Pv'
c        call chhermat (rpv,cpv,nbloch)
c        pause 'calcpv: check hermite'
c      endif

      deallocate (rp,cp,rw1,cw1,rw2,cw2)
      return
      end
c---------------------------------------------------------------------
      subroutine diagno0(nbloch,wpv,ovlp,evec,eval)
      implicit none
      integer(4):: nbloch,nmx,nev
c      call diagno(nbloch,wpv,ovlp,wdiag,iwdiag,evec,eval)
      complex(8),allocatable:: ovlpc(:,:),wpvc(:,:),evecc(:,:)
      real(8)::emx
      real(8):: wpv(nbloch,nbloch,2),ovlp(nbloch,nbloch,2)
     &   ,evec(nbloch,nbloch,2),eval(nbloch)
      allocate( ovlpc(nbloch,nbloch),wpvc(nbloch,nbloch),
     &          evecc(nbloch,nbloch))
      ovlpc= dcmplx(ovlp(:,:,1), ovlp(:,:,2))
      wpvc = dcmplx( wpv(:,:,1), wpv (:,:,2))
      nev  = nbloch
      nmx  = nbloch
c      print *,' goto diagcv-----------'
      call diagcv(ovlpc,wpvc, evecc, nbloch, eval, nmx, 1d99, nev)
      evec(:,:,1) =dreal(evecc)
      evec(:,:,2) =dimag(evecc)
      deallocate(wpvc,ovlpc,evecc)
      end

c-------------------------------------------
      subroutine getwk(ip,wibz,wqt,nqbz, nqibz,nstibz, nq0i,
     & wk)
      implicit none
      integer(4):: ip,nqibz,nq0i
      integer(4):: nqbz,nstibz(nqibz)
      real(8):: wibz(nqibz),wqt(nq0i),wk
      real(8)   :: wgtq0p
      integer(4):: bzcase
C- weight for k-point sampling
      if(bzcase()==2)then
        if(ip <=nqibz) then
          wk = wibz(ip)*.5d0
          if(nstibz(ip)/=0) wk = wibz(ip)*.5d0 *(1d0-wgtq0p()/nstibz(ip))
        elseif(ip>nqibz) then
          wk = wqt(ip-nqibz)* 1d0/dble(nqbz) * wgtq0p()
        endif
      else
        if(ip <= nqibz) then
          wk = wibz(ip)*.5d0 ! 0.5 for the normalization of wibz
        else
          wk = wqt(ip-nqibz)* 1d0/dble(nqbz)
        endif
        if(abs(wibz(1)-2d0/dble(nqbz)) >1d-10) then
          print *,' wibz(1) nqbz=',wibz(1),nqbz
          print *,' sum wibz=',sum(wibz)
          print *, ' ecorq2: Bug Stop! this may be a bug?  wibz(1) /= 1/dble(2*nqbz)'
          stop     ' ecorq2: Bug Stop! this may be a bug?  wibz(1) /= 1/dble(2*nqbz)'
        endif
      endif
      if(abs(sum(wibz)-2d0)>1d-10) then
        print *,' sum(wibz)=',sum(wibz)
        print *,' ecorq2: Bug Stop! this may be a bug? abs(sum(wibz)-2d0) /=0 '
        stop    ' ecorq2: Bug Stop! this may be a bug? abs(sum(wibz)-2d0) /=0 '
      endif
      end

c----------------------------------------------------------------
      subroutine  ecorq2 (zv,zw, nbloch,  iq,iw,ieceig,
     i    erpaqw, trpvqw, trlogqw)
C- Contribution to Ec from given (q(ip),iw(ix)), Tr(log(1-Pv) +Pv) -----------
c Takao modified from ecorq2. Work arrays are embedded.
c from rpaq.f, July 24, 2001
c diagonalize (1-Pv)(1-vP) instead of Pv
c
c from GW/lw/hecor.F, July 05, 2001
c TM
c
c 00.07.03
c calculates the correlated part of the total energy
c Erpa  = 1/4pi Int[iw=-inf:inf] Tr{log(1-Pv)+Pv}
c       = 1/2pi Int[iw=   0:inf] Tr{log(1-Pv)+Pv}
c Erpa is negative
c
c nbloch  = ngb total number of Bloch basis functions
      implicit real*8(a-h,o-z)
      real(8)::  ecqw
      real(8),allocatable :: rpv(:,:),cpv(:,:),
     &                       rw(:,:),cw(:,:),rv(:,:),cv(:,:),
     &                       wpv(:,:,:),ovlp(:,:,:), !wdiag(:),
     &                       evec(:,:,:),eval(:),eval2(:)
      complex(8) zw(nbloch,nbloch), zv(nbloch,nbloch)
      data      nsngl,ndble/4,8/
      logical :: debug = .false.
      integer(4)::verbose
c-----------------------------------------
      pi = 4d0*datan(1d0)
      if(verbose()>100) debug=.true.
      ngb= nbloch
      allocate( 
     &          rv(ngb,ngb),  cv(ngb,ngb),
     &          rpv(ngb,ngb), cpv(ngb,ngb),
     &          rw(ngb,ngb), cw(ngb,ngb),
     &          wpv(ngb,ngb,2), ovlp(ngb,ngb,2),
!     &          wdiag(11*ngb), iwdiag(ngb),
     &          evec(ngb,ngb,2), eval(ngb), eval2(ngb) )
      rv = dreal(zv)
      cv = dimag(zv)

c--- ix loop is now out of this routine.
c'      do      ix = 1,niw
      call calcpv2 (zw,rv,cv,
     d              nbloch,
     o              rpv,cpv,rw,cw )
      if(debug) write(*,*)' calcpv2_out sumcheck =', sum(abs(rpv)),sum(abs(cpv))

C--- diagonalize Pv
      do j = 1,nbloch
        do i = 1,nbloch
          wpv(i,j,1) = .5d0*(rpv(i,j)+rpv(j,i))
          wpv(i,j,2) = .5d0*(cpv(i,j)-cpv(j,i))
        enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccc
c eigenvlaue test1
c        wpv(:,:,1) = dreal(zv)
c        wpv(:,:,2) = dimag(zv)
c eigenvalue test2
c        wpv(:,:,1) = dreal(zw)
c        wpv(:,:,2) = dimag(zw)
cccccccccccccccccccccccccccccccccccccccccccccccc

      ovlp = 0d0
      do i = 1,nbloch
        ovlp(i,i,1)= 1d0
      enddo
      evec   = 0d0
      eval   = 0d0
      if (debug) write(*,*)'debug> diagno Pv in'
      call diagno0(nbloch,wpv,ovlp,evec,eval)

!      wdiag  = 0.d0
!      iwdiag = 0
ccccccccccccccccccccc
c      call diagno(nbloch,wpv,ovlp,wdiag,iwdiag,evec,eval)
c      do i = 1,nbloch
c         write(6,"(i5,256d10.2)")i,(wpv(i,j,1),j=1,10)
c         write(6,"(i5,256d10.2)")i,(wpv(i,j,2),j=1,10)
c      enddo
ccccccccccccccccccccccc
c      do i=1,nbloch
c      do j=1,nbloch
c         wpv(i,j,1:2)=0d0
c         ovlp(i,j,1:2)=0d0
c         if(i==j) wpv(i,j,1)=id0
c         if(i==j) ovlp(i,j,1)=1d0
c      enddo
c      enddo
c      evec=0d0
c      eval=0d0
c      stop 'xxxxxxxxxxxxxxxxxxxx test end xxxxxxx'
ccccccccccccccccccccc

C--- diagonalize (1-Pv)(1-vP)
      do       j = 1,nbloch
        do       i = 1,nbloch
          wpv(i,j,1) = rw(i,j)
          wpv(i,j,2) = cw(i,j)
          ovlp(i,j,1)= 0d0
          ovlp(i,j,2)= 0d0
        enddo
      enddo
      do       i = 1,nbloch
        ovlp(i,i,1)= 1.d0
      enddo
!      wdiag  = 0.d0 !      iwdiag = 0
      evec   = 0d0
      eval2  = 0d0
      if (debug) write(*,*)'debug> diagno (1-Pv)(1-vP) in'
      call diagno0(nbloch,wpv,ovlp,evec,eval2)

      trlogqw = 0.5d0*sum(dlog(eval2(1:nbloch)))
      trpvqw  = sum (eval(1:nbloch))
      erpaqw  = sum (.5d0*dlog(eval2(1:nbloch)) + eval(1:nbloch) )

! check write
      if(ieceig>0) then
        close(ieceig)
        open (ieceig, file='rpa_eigen.chk',access='append')
        write(ieceig,*)
        write(ieceig,
     &  "('--- iq iw: ',2i6,' Eigen of  Pv+vP (1-Pv)(1-vP) ')") iq,iw
        do i = 1,nbloch
          write(ieceig,"(i5,d17.6,d17.6)")i,eval(i),eval2(i)
        enddo
      endif
      deallocate(rpv,cpv,rw,cw, wpv,ovlp,evec,eval,eval2)
      return
      end




c---------------------------------------------------------------------
      subroutine sqrtmat (rmat,cmat,
     w                    wmat,ovlp,wdiag,iwdiag,
     w                    evec,eval,rw1,cw1,
     d                    ldim,
     o rmat2,cmat2)

      implicit double precision (a-h,o-z)

      dimension   rmat(ldim,ldim),cmat(ldim,ldim),
     .            wmat(ldim,ldim,2),ovlp(ldim,ldim,2),
     .            wdiag(11*ldim),iwdiag(ldim),
     .            evec(ldim,ldim,2),eval(ldim),
     .            rw1(ldim,ldim),cw1(ldim,ldim),
     .            rmat2(ldim,ldim),cmat2(ldim,ldim)

      data tol /1.d-8/

      logical :: debug = .false.

c-----------------------------------------
      do       j  = 1,ldim
        do       i  = 1,ldim
          wmat(i,j,1) = rmat(i,j)
          wmat(i,j,2) = cmat(i,j)
          ovlp(i,j,1) = 0.d0
          ovlp(i,j,2) = 0.d0
        enddo
      enddo
      do        i = 1,ldim
        ovlp(i,i,1) = 1.d0
      enddo

      wdiag       = 0.d0
      iwdiag      = 0
      evec        = 0.d0
      eval        = 0.d0

cccccccccccccccccccccccccccccccccc
c      call diagno(ldim,wmat,ovlp,wdiag,iwdiag,evec,eval)
      call diagno0(ldim,wmat,ovlp,evec,eval)
ccccccccccccccccccccccccccccc

c v^{1/2}
      do        i = 1,ldim
        eval(i)     = dsqrt(eval(i))
      enddo

c U{-1}
      do       j  = 1,ldim
        do       i  = 1,ldim
          wmat(i,j,1) =  evec(j,i,1)
          wmat(i,j,2) = -evec(j,i,2)
        enddo
      enddo

c      if (debug) then
c debug> check if evec is unitary
c      do       i = 1,ldim
c      do       j = 1,ldim
c      rtmp       = 0.d0
c      ctmp       = 0.d0
c      do       k = 1,ldim
c      rtmp       = rtmp
c     .     + evec(i,k,1)*wmat(k,j,1) - evec(i,k,2)*wmat(k,j,2)
c      ctmp       = ctmp
c     .     + evec(i,k,1)*wmat(k,j,2) + evec(i,k,2)*wmat(k,j,1)
c      enddo
c      if ((dabs(rtmp).gt.tol) .or. (dabs(ctmp).gt.tol))
c     .write(*,*)i,j,rtmp,ctmp
c      enddo
c      enddo
c      pause 'sqrtmat: unitary'
c      endif

c v^{1/2} U^{-1}
      do       j  = 1,ldim
        do       i  = 1,ldim
          rw1(i,j)    =  eval(i) * wmat(i,j,1)
          cw1(i,j)    =  eval(i) * wmat(i,j,2)
        enddo
      enddo

c U v^{1/2} U^{-1}
      do        j = 1,ldim
        do        i = 1,ldim
          rtmp        = 0.d0
          ctmp        = 0.d0
          do        k = 1,ldim
            rtmp        = rtmp 
     .            + wmat(k,i,1)*rw1(k,j) + wmat(k,i,2)*cw1(k,j)
            ctmp        = ctmp 
     .            + wmat(k,i,1)*cw1(k,j) - wmat(k,i,2)*rw1(k,j)
          end do
          rmat2(i,j)  = rtmp
          cmat2(i,j)  = ctmp
        enddo
      enddo

c      if (debug) then
cdebug> check mat2*mat2 = mat
c      do        i = 1,ldim
c      write(*,*)'i=',i
c      do        j = 1,ldim
c      wmat(i,j,1) = 0.d0
c      wmat(i,j,2) = 0.d0
c      do        k = 1,ldim
c      wmat(i,j,1) = wmat(i,j,1)
c     .             + rmat2(i,k)*rmat2(k,j) - cmat2(i,k)*cmat2(k,j)
c      wmat(i,j,2) = wmat(i,j,2)
c     .             + rmat2(i,k)*cmat2(k,j) + cmat2(i,k)*rmat2(k,j)
c      enddo
c      rdiff       = wmat(i,j,1) - rmat(i,j)
c      cdiff       = wmat(i,j,2) - cmat(i,j)
c      if ((dabs(rdiff).gt.tol) .or. (dabs(cdiff).gt.tol))
c     .write(*,*)i,j,rdiff,cdiff
c      enddo
c      enddo
c      pause 'check mat2 out'
c      endif

      return
      end
c---------------------------------------------------------------------
      subroutine hermat  (rmat,cmat,ldim)

      implicit double precision (a-h,o-z)
      real(8)   rmat(ldim,ldim),cmat(ldim,ldim)

      do      i = 1,ldim
        do      j = i,ldim
          w1        = 0.5d0 * (rmat(i,j) + rmat(j,i))
          w2        = 0.5d0 * (cmat(i,j) - cmat(j,i))
          rmat(i,j) =  w1
          cmat(i,j) =  w2
          rmat(j,i) =  w1
          cmat(j,i) = -w2
        end do
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine chhermat (rmat,cmat,ldim)

      implicit double precision (a-h,o-z)

      dimension   rmat(ldim,ldim),cmat(ldim,ldim)
      data tol /1.d-8/

      do i = 1,ldim
        do j = 1,ldim
          w1   = rmat(i,j) - rmat(j,i)
          w2   = cmat(i,j) + cmat(j,i)
          if ((dabs(w1).gt.tol).or.(dabs(w2).gt.tol)) then
            write(*,*)i,j,w1,w2
            pause 'chhermat: non-hermite!'
          endif
        end do
      end do

      return
      end
c---------------------------------------------------------------------
