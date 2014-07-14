      subroutine init_unkg(is,qbz,ginv,ef,lein,
     i                     iko_ix,iko_fx,iko_i,iko_f,
     i                     iki_i,iki_f,
     d                     nwf,nband,nqbz,
     o                     amnk,cnk)

      implicit real*8(a-h,o-z)
      
      logical :: lein
c      complex(8),allocatable :: psig(:,:,:)
      complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz),
     &              cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
      real(8) :: qbz(3,nqbz),ginv(3,3),
     &           ovlp(iko_ix:iko_fx)
      integer(4) :: iko_i(nqbz),iko_f(nqbz),
     &              iki_i(nqbz),iki_f(nqbz)


c initialize
      cnk = (0d0,0d0)
      amnk = (0d0,0d0)

c read psig(it,iwf,iqbz) = < psi(it,iqbz) | g(iwf) >
      call readpsig(is,iko_ix,iko_fx,
     d              nqbz,nwf,
c     o              psig)
     o              amnk)

      call amnk2unk(amnk,
     i              iko_ix,iko_fx,iko_i,iko_f,
     d              nwf,nqbz,
     o              cnk)

c inner energy window
      if (lein)
     >call init_iew(iko_ix,iko_fx,iko_i,iko_f,
     i              iki_i,iki_f,
     d              nwf,nband,nqbz,
     u              cnk)

      return
      end
c-----------------------------------------------------------------------
      subroutine pick_nwf(ovlp,iti,itf,nwf,isort)

      implicit real*8(a-h,o-z)
      
      real(8) :: ovlp(iti:itf),otmp(iti:itf)
      integer(4) :: isort(nwf),istate(iti:itf)

c initial
      do it = iti,itf
         istate(it) = it
         otmp(it) = ovlp(it)
      enddo

c sorting
      do it1 = iti,itf-1
      do it2 = it1+1,itf
         if (ovlp(it1) .lt. ovlp(it2)) then
            tmp = ovlp(it2)
            ovlp(it2) = ovlp(it1)
            ovlp(it1) = tmp
            itmp = istate(it2)
            istate(it2) = istate(it1)
            istate(it1) = itmp
         endif
      enddo
      enddo

c sort check
      do i1 = iti,itf-1
         it1 = istate(i1)
         it2 = istate(i1+1)
         if (otmp(it1) .lt. otmp(it2)) stop 'pick_nwf: sort error'
      enddo

c pick largest nwf states
      do it = 1,nwf
         itmp = iti - 1 + it
         isort(it) = istate(itmp)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine read_cnq0(ifhoev,is,qwf0,qbz,ginv,ef,
     i                    itq,
     d                    nwf,nband,nqbz,
     o                    cnq0)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: cks(:,:),hks(:,:),oks(:,:)
      complex(8) :: cnq0(nband,nwf)
      real(8),allocatable :: eval(:)
      real(8) :: qwf0(3),qbz(3,nqbz),q(3),ginv(3,3)
      real(8) :: rydberg
      integer(4) :: itq(nwf)

      iq0 = iqindx(qwf0,ginv,qbz,nqbz)
      cnq0 = (0d0,0d0)

c open
      if (is.eq.1) then 
         ifhoev = iopen('HOEV.UP',0,0,0)
      elseif (is.eq.2) then
         ifhoev = iopen('HOEV.DN',0,0,0)
      else
         stop 'read_cnq0: iopen error'
      endif

c read
      read(ifhoev)ndimh,nqtot
      if (ndimh .ne. nband) stop 'read_cnq0: nband error'
      if (nqtot .lt. nqbz) stop 'read_cnq0: nqbz error'

      allocate(hks(nband,nband),oks(nband,nband),
     &         cks(nband,nband),eval(nband))

      do iq = 1,nqbz
         read(ifhoev)iq2,q(1:3)
         read(ifhoev)hks(1:nband,1:nband)
         read(ifhoev)oks(1:nband,1:nband)
         read(ifhoev)cks(1:nband,1:nband)
         read(ifhoev)eval(1:nband)

         iq3 = iqindx(q,ginv,qbz,nqbz)
         if (iq3 .ne. iq) stop 'read_cnq0: iqindx error'
         if (iq3 .eq. iq0) then
            do it = 1,nwf
               cnq0(:,it) = cks(:,itq(it))
c               ev = (eval(itq(it)) - ef) * rydberg()
c               write(*,*)'iwf,nwf,ev',ev
            enddo
            goto 99
         endif
      enddo
      stop 'read_cnq0: cannot find q0'
 99   continue

      deallocate(hks,oks,cks,eval)

c close
      if (is.eq.1) then 
         ifi = iclose('HOEV.UP')
      else
         ifi = iclose('HOEV.DN')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_amnk(ifhoev,is,qwf0,qbz,ginv,
     i                    cnq0,
     i                    iko_ix,iko_fx,iko_i,iko_f,
     d                    nwf,nband,nqbz,
     o                    amnk)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: cks(:,:),hks(:,:),oks(:,:),
     &                          wmat(:,:)
      complex(8) :: cnq0(nband,nwf),
     &              amnk(iko_ix:iko_fx,nwf,nqbz)
      real(8),allocatable :: eval(:)
      real(8) :: qwf0(3),qbz(3,nqbz),q(3),ginv(3,3)
      integer(4) :: iko_i(nqbz),iko_f(nqbz)


c open
      if (is.eq.1) then 
         ifhoev = iopen('HOEV.UP',0,0,0)
      elseif (is.eq.2) then
         ifhoev = iopen('HOEV.DN',0,0,0)
      else
         stop 'get_amnk: iopen error'
      endif

c read
      read(ifhoev)ndimh,nqtot
      if (ndimh .ne. nband) stop 'get_amnk: nlmto error'
      if (nqtot .lt. nqbz) stop 'get_amnk: nqbz error'

      allocate(hks(nband,nband),oks(nband,nband),
     &         cks(nband,nband),eval(nband),
     &         wmat(nband,nwf))

c initialize
      amnk = (0d0,0d0)

      do iq = 1,nqbz
         read(ifhoev)iq2,q(1:3)
         read(ifhoev)hks(1:nband,1:nband)
         read(ifhoev)oks(1:nband,1:nband)
         read(ifhoev)cks(1:nband,1:nband)
         read(ifhoev)eval(1:nband)

         iq3 = iqindx(q, ginv,qbz,nqbz)
         if (iq3 .ne. iq) stop 'get_amnk: iqindx error'


c wmat = cnq0 * oks
         wmat = (0d0,0d0)
         do in = 1,nwf
         do ij = 1,nband
            do ii = 1,nband
               wmat(ij,in) = wmat(ij,in) 
     &                     + cnq0(ii,in) * oks(ij,ii)
            enddo
         enddo
         enddo

c amnk = cks^{*} * wmat
         do im = iko_i(iq),iko_f(iq)
         do in = 1,nwf
            do ij = 1,nband
               amnk(im,in,iq) = amnk(im,in,iq) + 
     &             dconjg(cks(ij,im)) * wmat(ij,in)
            enddo
         enddo
         enddo

      enddo

      deallocate(hks,oks,cks,eval,wmat)

c close
      if (is.eq.1) then 
         ifi = iclose('HOEV.UP')
      else
         ifi = iclose('HOEV.DN')
      endif

 999  format(i5,3f16.8)

      return
      end
c-----------------------------------------------------------------------
      subroutine amnk2unk(amnk,
     i                    iko_ix,iko_fx,iko_i,iko_f,
     d                    nwf,nqbz,
     o                    cnk0)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-3)
!      complex(8),allocatable :: evecc(:,:),smat(:,:),wmat(:,:),
!     &                          wmat2(:,:),wmat3(:,:)
!      real (8),allocatable :: eval(:)
      complex(8),allocatable :: aa(:,:),cc(:,:),zz(:,:),vv(:,:)
      real (8),allocatable :: dd(:)
      complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz),
     &              cnk0(iko_ix:iko_fx,nwf,nqbz),ctmp
      integer(4) :: iko_i(nqbz),iko_f(nqbz)


c singular value decomposition, 061003
      nks = iko_fx - iko_ix + 1
      allocate (aa(nks,nks),zz(nks,nks),vv(nwf,nwf),dd(nwf))
      cnk0 = (0d0,0d0)
      do iq = 1,nqbz
         aa(1:nks,1:nwf) = amnk(iko_ix:iko_fx,1:nwf,iq)
         call zgesvdmn(nks,nwf,aa,dd,zz,vv)
         do ij = iko_i(iq),iko_f(iq)
            jj = ij - iko_ix + 1
         do ii = 1,nwf
            do ik = 1,nwf
              cnk0(ij,ii,iq) = cnk0(ij,ii,iq) + zz(jj,ik)*vv(ik,ii)
            enddo
         enddo
         enddo
      enddo ! iq
      deallocate (aa,zz,vv,dd)

      return
      end
c-----------------------------------------------------------------------
      subroutine init_iew(iko_ix,iko_fx,iko_i,iko_f,
     i                    iki_i,iki_f,
     d                    nwf,nband,nqbz,
     u                    cnk)

      implicit real*8(a-h,o-z)


      complex(8),allocatable :: cnk2(:,:),vnk(:,:),mat(:,:),evec(:,:)
      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
      real(8),allocatable :: eval(:)
      integer(4),allocatable :: it(:)
      integer(4) :: iko_i(nqbz),iko_f(nqbz),
     &              iki_i(nqbz),iki_f(nqbz)

      allocate(cnk2(iko_ix:iko_fx,nwf),
     &          vnk(iko_ix:iko_fx,nwf))
      nox = iko_fx - iko_ix + 1

      do iq = 1,nqbz

c no innner energy window
      if (iki_i(iq) .eq.0) then
         if (iki_f(iq) .ne. -1) stop 'init_iew: iki_f error'

c inner energy window
      else
         nout = iko_f(iq) - iko_i(iq) + 1
         nin  = iki_f(iq) - iki_i(iq) + 1
         if (nin .lt. 1) stop 'init_iew: nin error'
         cnk2 = (0d0,0d0)

c Nin(k) states in the inner energy window
         do il = 1, nin
            ij = iki_i(iq) - 1 + il
            cnk2(ij,il) = (1d0,0d0)
         enddo

c Nwf - Nin(k) states out of Nout states
! |vnk>
         vnk(:,:) = cnk(:,:,iq)
         vnk(iki_i(iq):iki_f(iq),:) = 0d0

! diagonalization
         allocate(mat(nox,nox),evec(nox,nox),eval(nox))
         mat = 0d0
         do ii = 1,nox
         do ij = 1,nox
            do il = 1,nwf
               mat(ii,ij) = mat(ii,ij)
     &            + vnk(ii+iko_ix-1,il) * dconjg(vnk(ij+iko_ix-1,il))
            enddo
         enddo
         enddo
         call chk_hm(mat,nox)
         call diag_hm(mat,nox,eval,evec)

! cnk2
         do il = 1, nwf-nin
         do ij = iko_ix, iko_fx
            cnk2(ij,il+nin) = evec(ij-iko_ix+1,nox-il+1)
         enddo
         enddo

         deallocate(mat,evec,eval)
         call chk_on(cnk2(iko_ix:iko_fx,:),nox,nwf)

c new cnk(:,:,q)
         cnk(:,:,iq) = cnk2(:,:)

c end of if (iki_i == 0)
      endif
c end of iq-loop
      enddo

      deallocate(cnk2,vnk)

      return
      end
c-----------------------------------------------------------------------
      subroutine   getupu(isc,
     i                    uumat,cnk,
     i                    lein,alpha_in,iq,ikbidx,
     i                    iko_ix,iko_fx,
     i                    iko_i,iko_f,
     i                    iki_i,iki_f,
     i                    ikbo_i,ikbo_f,
     i                    ikbi_i,ikbi_f,
     d                    nwf,nbb,nqbz,
     u                    upu)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: wmat(:,:),wmat2(:,:)
      complex(8) :: upu(iko_ix:iko_fx,iko_ix:iko_fx,nbb), 
     .            uumat(iko_ix:iko_fx,iko_ix:iko_fx,nbb),
     .              cnk(iko_ix:iko_fx,nwf,nqbz),
     .              ctmp
      integer (4) :: ikbidx(nbb),
     .               ikbo_i(nbb),ikbo_f(nbb),
     .               ikbi_i(nbb),ikbi_f(nbb)
      logical :: lein


      if (isc .eq. 1) then
         alpha = 1d0
         upu(:,:,:) = (0d0,0d0)
      else
         alpha = alpha_in
         upu(:,:,:) = upu(:,:,:) * (1d0-alpha)
      endif

c      if (iki_i .eq. 0) then
c         nin = 0
c      else
c         nin = iki_f - iki_i + 1
c      endif
      nin = iki_f - iki_i + 1

      if (nin .ge. nwf) then
         stop 'getupu: Nin >= Nwf'
      endif


      allocate(wmat(iko_ix:iko_fx,iko_ix:iko_fx), 
     &        wmat2(iko_ix:iko_fx,iko_ix:iko_fx))

      do ibb = 1,nbb
         iqb = ikbidx(ibb)


c wmat = cnk * cnk^{*}
         wmat = (0d0,0d0)
         do inp  = ikbo_i(ibb),ikbo_f(ibb)
         do imp  = ikbo_i(ibb),ikbo_f(ibb)
            do il  = nin+1,nwf
                 wmat(inp,imp) = wmat(inp,imp)
     &         + cnk(imp,il,iqb)*dconjg(cnk(inp,il,iqb))
            enddo
         enddo
         enddo

c wmat2 = uumat^{*} * wmat
         wmat2 = (0d0,0d0)
         do in   = iko_i,iko_f
         do imp  = ikbo_i(ibb),ikbo_f(ibb)
            do inp  = ikbo_i(ibb),ikbo_f(ibb)
                 wmat2(imp,in) = wmat2(imp,in)
     &        +  dconjg(uumat(in,inp,ibb)) * wmat(inp,imp)
            enddo
         enddo
         enddo

c upu = upu + uumat * wmat2 * alpha
         do im   = iko_i,iko_f
         do in   = iko_i,iko_f
            do imp  = ikbo_i(ibb),ikbo_f(ibb)
                 upu(im,in,ibb) = upu(im,in,ibb) 
     &        +  uumat(im,imp,ibb) * wmat2(imp,in) * alpha
            enddo
         enddo
         enddo

      enddo

      deallocate(wmat,wmat2)

      return
      end
c-----------------------------------------------------------------------
      subroutine dimz(lein,iko_i,iko_f,iki_i,iki_f,
     o                ndz,nin)

      implicit real*8(a-h,o-z)

      logical :: lein


      nout = iko_f - iko_i + 1
      nin  = iki_f - iki_i + 1
      ndz  = nout - nin

check 
      if (iki_i .eq. 0) then
         if (iki_f .ne. -1) stop 'dimz: iki_f error'
      else
         if (iko_i .gt. iki_i) stop 'dimz: ik_i error'
         if (iko_f .lt. iki_f) stop 'dimz: ik_f error'
         if (iki_i .gt. iki_f) stop 'dimz: iki error'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine getzmn(upu,wbb,lein,
     i                  iko_ix,iko_fx,
     i                  iko_i,iko_f,
     i                  iki_i,iki_f,
     d                  nwf,nbb,nqbz,ndz,
     o                  zmn)

      implicit real*8(a-h,o-z)

      complex(8) :: upu(iko_ix:iko_fx,iko_ix:iko_fx,nbb), 
     &              zmn(ndz,ndz)
      real(8) :: wbb(nbb)
      integer(4) :: it(ndz)
      logical :: lein

      zmn = (0d0,0d0)

c no inner energy window
      if (iki_i.eq.0) then
         no = iko_f - iko_i + 1
         if (no .ne. ndz) stop 'getzmn: ndz error'
         do ibb = 1, nbb
          zmn(1:ndz,1:ndz) = zmn(1:ndz,1:ndz)
     &                 +  wbb(ibb) * upu(iko_i:iko_f,iko_i:iko_f,ibb)
         enddo

c inner energy window
      else

         j = 0
         do i = iko_i,iki_i-1
            j = j + 1
            it(j) = i
         enddo
         do i = iki_f+1,iko_f
            j = j + 1
            it(j) = i
         enddo
         if (j.ne.ndz) stop 'getzmn: ndz error'

         do im = 1,ndz
         do in = 1,ndz
         do ibb = 1, nbb
            zmn(im,in) = zmn(im,in)
     &                 + wbb(ibb) * upu(it(im),it(in),ibb)
         enddo
         enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_hm(zmat,ndim)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-4)
      complex(8):: zmat(ndim,ndim)

      do i = 1,ndim
      do j = 1,ndim
         dr = dabs(dreal(zmat(i,j)) - dreal(zmat(j,i)))
         di = dabs(dimag(zmat(i,j)) + dimag(zmat(j,i)))
         dc = abs(zmat(i,j))
         if (dr .gt. eps) stop 'chk_hm: real part error'
         if (di .gt. eps) stop 'chk_hm: imag part error'
c         if (dr/dc .gt. eps) stop 'chk_hm: real part error'
c         if (di/dc .gt. eps) stop 'chk_hm: imag part error'
c         if (dr/dc .gt. eps) then 
c             write(*,*)zmat(i,j),i,j
c             write(*,*)zmat(j,i),i,j
c             stop 'chk_hm: real part error'
c         endif
c         if (di/dc .gt. eps) then
c             write(*,*)zmat(i,j),i,j
c             write(*,*)zmat(j,i),i,j
c             stop 'chk_hm: imag part error'
c         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_um(zmat,ndim)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-4)
      complex(8):: zmat(ndim,ndim),cij

      do i = 1,ndim
      do j = 1,ndim
         cij = (0d0,0d0)
         do k = 1,ndim
            cij = cij + dconjg(zmat(k,i))*zmat(k,j)
         enddo
         if (i .eq. j) cij = cij - 1d0
         rij = abs(cij)
         if (rij .gt. eps) stop 'chk_um: error'
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_on(zmat,n1,n2)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-4)
      complex(8):: zmat(n1,n2),cij

      do i = 1,n2
      do j = 1,n2
         cij = sum(dconjg(zmat(:,i))*zmat(:,j))
         if (i .eq. j) cij = cij - 1d0
         rij = abs(cij)
         if (rij .gt. eps) stop 'chk_on: error'
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine diag_hm(zmat,ndim,eval,evecc)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
      complex(8):: zmat(ndim,ndim),evecc(ndim,ndim)
      real(8):: eval(ndim)

      allocate(zmat2(ndim,ndim),ovlpc(ndim,ndim))

      nev  = ndim
      nmx  = ndim

      zmat2 = zmat

      ovlpc = (0d0,0d0)
      do i=1,ndim
         ovlpc(i,i) = (1d0,0d0)
      enddo

      evecc = (0d0,0d0)
      eval = 0d0

      call diagcv(ovlpc,zmat2, evecc, ndim, eval, nmx, 1d99, nev)

      deallocate(zmat2,ovlpc)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_eval(wbb,evz,nbb,ndz)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-4)
      real(8):: wbb(nbb),evz(ndz)


c check order
      do i = 1,ndz-1
         e1 = evz(i)
         e2 = evz(i+1)
         if (e1 .gt. e2) stop 'chk_eval: order is wrong'
      enddo

c check 
      ws = sum(wbb)
      if (evz(ndz) .gt. ws) then 
         write(*,*) 'chk_eval: eval is too large'
         write(*,*)'sum(wbb) =',ws
         do i = ndz,1,-1
            write(*,*)'i,evz(i)=',i,evz(i)
         enddo
         stop
      endif

      return
      end
c-----------------------------------------------------------------------
       subroutine new_cnk(cnk,evecc,iq,
     i                    iko_ix,iko_fx,
     i                    iko_i,iko_f,
     i                    iki_i,iki_f,
     d                    nwf,ndz,
     o                    cnk2)

      implicit real*8(a-h,o-z)

      complex(8) :: cnk(iko_ix:iko_fx,nwf),
     i              cnk2(iko_ix:iko_fx,nwf),
     i              evecc(ndz,ndz)
      integer(4) :: it(ndz)

c initialize
      cnk2(:,:) = (0d0,0d0)

c no inner energy window
      if (iki_i .eq. 0) then
         nout = iko_f - iko_i + 1
         if (nout .ne. ndz) stop 'new_cnk: ndz error'
         if (nwf .gt. ndz) stop 'new_cnk: nwf error'

         il2 = ndz + 1
         do il = 1,nwf
            il2 = il2 - 1
            cnk2(iko_i:iko_f,il) = evecc(1:ndz,il2)
         enddo

c inner energy window
      else
         nout = iko_f - iko_i + 1
         nin  = iki_f - iki_i + 1
         if (nout-nin .ne. ndz) stop 'new_cnk: ndz error'
         if (nwf-nin .gt. ndz) stop 'new_cnk: nwf error'
         cnk2(:,1:nin) = cnk(:,1:nin)

! pick nwf-nin states with the largest eigenvalues
         j = 0
         do i = iko_i,iki_i-1
            j = j + 1
            it(j) = i
         enddo
         do i = iki_f+1,iko_f
            j = j + 1
            it(j) = i
         enddo

         il2 = ndz + 1
         do il = nin+1,nwf
            il2 = il2 - 1
            do in = 1,ndz
               cnk2(it(in),il) = evecc(in,il2)
            enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine  get_omgik(wbb,evz,
     i                      iko_i,iko_f,
     i                      iki_i,iki_f,
     d                      nbb,nwf,ndz,
     o                      omgik)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-4)
      real(8):: wbb(nbb),evz(ndz)

      nin = iki_f - iki_i + 1
      nn  = nwf - nin

check
      if (iki_i .eq. 0) then
         if (iki_f .ne. -1) stop 'get_omgik: iki_f error'
      else
         if (nin .le. 0) stop 'get_omgik: nin error'
      endif

c
      if (ndz .lt. nn) stop 'get_omgik: ndz error'
      esum = 0d0
      do i = 1,nn
         j = ndz + 1 - i
         esum = esum + evz(j)
      enddo
c
      omgik = dble(nn)*sum(wbb) - esum

      return
      end
c-----------------------------------------------------------------------
      subroutine zgesvdmn(ngb1,ngb2,zzz, SS,UU,VT)
      implicit none
      integer(4)::lwork,info,ngb1,ngb2,i
      complex(8):: zzz(ngb1,ngb2),UU(ngb1,ngb1),VT(ngb2,ngb2)
      real(8):: ss(ngb2)
      real(8),allocatable:: rwork(:)
      complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
      lwork=4*ngb1
      allocate(work(LWORK),rwork(5*ngb1)) 
      call zgesvd('A','A',ngb1,ngb2,zzz,ngb1,SS,UU,ngb1,VT,ngb2,work,lwork,rwork,info) 
      deallocate(work,rwork)
      end
!--------------------------------------------------------------------------
      subroutine zgesvdnn(ngb,zzz, SS,UU,VT)
C--- SVD zzz= matmul(UU,matmul(SS,VT)) ------------
c$$$C--- SVD of chipm test !SVD procedure is not required to calculate <m|chi|m>
c$$$      lwork=4*ngb
c$$$      allocate(work(LWORK),rwork(5*ngb),zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb),VTT(ngb,ngb),ss0(ngb))
c$$$      zzz= matmul(transpose(conjg(ovlpi)), matmul(zxq(1:ngb,1:ngb,iw), ovlpi) )
c$$$      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS0,UU,ngb,VT,ngb,work,lwork,rwork,info) 
c$$$      write(6,*)' -------SVD: Oinv*chipm*Oinv ----------------'
c$$$      aaax = 0d0
c$$$      do i=1,ngb
c$$$        addx= sum(svec(1:nbloch)*uu(1:nbloch,i)) *ss0(i)* sum(VT(i,1:nbloch)*svec(1:nbloch)) 
c$$$        write(6,"(' SVD OcO: eig_k <m|chi|m>_k=',i4,2x, d13.5,2x,2d14.6)")i,SS0(i),addx
c$$$        if(i<25) aaax= aaax+ addx
c$$$      enddo 
c$$$      aaax= mmnorm**2/aaax
c$$$      deallocate(work,rwork,zzz,uu,vt,vtt)
c$$$      deallocate(ovlpi)
      implicit none
      integer(4)::lwork,info,ngb,i
      complex(8):: zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb)
      real(8):: ss(ngb)
      real(8),allocatable:: rwork(:)
      complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
      lwork=4*ngb
      allocate(zw0bk(ngb,ngb))
      allocate(work(LWORK),rwork(5*ngb)) !,VTT(ngb,ngb))
      zw0bk = zzz
!      print *,' zgesvdnn: singular value decomp '
      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS,UU,ngb,VT,ngb,work,lwork,rwork,info) 
!      do i=1,ngb
!         write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
!         vtt(i,:)=ss(i)*vt(i,:)
!      enddo
!      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
!     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
!      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk))) 
!     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
!      deallocate(vtt)
      end
!--------------------------------------------------------------------------
