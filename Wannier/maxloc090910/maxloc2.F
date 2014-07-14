      subroutine diag_unk(is,qbz,
     i                    iko_ix,iko_fx,iko_i,iko_f,
     d                    nband,nwf,nqbz,
     u                    cnk,
     o                    eunk)

      use m_readeigen

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-3)
      complex(8),allocatable :: evecc(:,:),ham(:,:),cnk2(:,:)
      real (8),allocatable :: eval(:),eks(:)
      complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz),
     &              cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
      real (8) :: qbz(3,nqbz),q(3),eunk(nwf,nqbz)              
      integer(4) :: iko_i(nqbz),iko_f(nqbz)


      allocate (ham(nwf,nwf),eks(nband),
     &          evecc(nwf,nwf),eval(nwf),
     &          cnk2(iko_ix:iko_fx,nwf))

      do iq = 1,nqbz

c read eigenvalues
      q(:) = qbz(:,iq)
      call readeval (q,is,eks)

c construct H
      ham = (0d0,0d0)
      do ii = 1,nwf
      do ij = 1,nwf
         do ik = iko_i(iq),iko_f(iq)
            ham(ii,ij) = ham(ii,ij) + 
     &       dconjg(cnk(ik,ii,iq))*cnk(ik,ij,iq)*eks(ik)
         enddo
      enddo
      enddo

c diagonalize H
      call chk_hm(ham,nwf)
      call diag_hm(ham,nwf,eval,evecc)


c cnk(n,l) = S[m] evecc(m,l)*c(n,m)
      cnk2 = (0d0,0d0)
      do il = 1,nwf
      do im = 1,nwf
         do in = iko_i(iq),iko_f(iq)
            cnk2(in,il) = cnk2(in,il) + 
     &         evecc(im,il) * cnk(in,im,iq)
         enddo
      enddo
      enddo
 
      cnk(:,:,iq) = cnk2
      eunk(:,iq) = eval(:)

c end of iq-loop
      enddo

      deallocate (ham,eks,evecc,eval,cnk2)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_eunk(is,qbz,eunk,ef,
     d                    nqbz,nband,nwf)

      use m_readeigen

      implicit real*8(a-h,o-z)

      real(8) :: qbz(3,nqbz),eunk(nwf,nqbz)
      real(8),allocatable :: eks(:)

      allocate(eks(nband))

      do iq = 1,nqbz
         call readeval (qbz(:,iq),is,eks)

         write(*,*)'Diag energy',iq
         do iband = 1,nwf
            eev = (eunk(iband,iq)-ef)*rydberg()
            write(*,*)iband,eunk(iband,iq),eev
         enddo

         write(*,*)'KS energy  ',iq
         do iband = 1,nband
            eev = (eks(iband)-ef)*rydberg()
            write(*,*)iband,eks(iband),eev
         enddo
      enddo

      deallocate(eks)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_cnk(cnk,
     i                   iko_ix,iko_fx,iko_i,iko_f,
     d                   nband,nwf,nqbz)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-4)
      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
      integer (4) :: iko_i(nqbz),iko_f(nqbz)

      do iq = 1,nqbz

         do im = 1,nwf
         do in = 1,nwf
            ctmp = (0d0,0d0)
            do ii = iko_i(iq),iko_f(iq)
               ctmp = ctmp + dconjg(cnk(ii,im,iq))*cnk(ii,in,iq)
            enddo
            if (in .eq. im) ctmp = ctmp - 1d0
            a = dabs(dreal(ctmp))
            b = dabs(dimag(ctmp))
            if (a.gt.eps) then
               write(*,*)'chk_cnk: real error,iq,im,in, Re'
               write(*,*)iq,im,in,a
            endif
            if (b.gt.eps) then
               write(*,*)'chk_cnk: imag error,iq,im,in, Im'
               write(*,*)iq,im,in,b
            endif
         enddo
         enddo

      enddo

      stop 'chk_cnk: end'

      return
      end
c-----------------------------------------------------------------------
      subroutine init_mmn(cnk,uumat,ikbidx,
     i             iko_ix,iko_fx,iko_i,iko_f,ikbo_i,ikbo_f,
     d             nwf,nqbz,nbb,
     o             mmn)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: wmat(:,:)
      complex(8):: cnk(iko_ix:iko_fx,nwf,nqbz),
     c             uumat(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz),
     c             mmn(nwf,nwf,nbb,nqbz)
      integer(4):: ikbidx(nbb,nqbz),
     c             iko_i(nqbz),iko_f(nqbz),
     c             ikbo_i(nbb,nqbz),ikbo_f(nbb,nqbz)

      allocate(wmat(iko_ix:iko_fx,nwf))

      mmn = (0d0,0d0)
      do iq = 1,nqbz
      do ibb = 1,nbb
         iqb = ikbidx(ibb,iq)

c wmat = cnk * uumat
         wmat = (0d0,0d0)
         do imp = iko_i(iq),iko_f(iq)
         do in = 1,nwf
            do inp = ikbo_i(ibb,iq),ikbo_f(ibb,iq)
               wmat(imp,in) = wmat(imp,in)
     &       + cnk(inp,in,iqb) * uumat(imp,inp,ibb,iq)
            enddo
         enddo
         enddo

c mmn = cnk^{*} * wmat
         do im = 1,nwf
         do in = 1,nwf
            do imp = iko_i(iq),iko_f(iq)
               mmn(im,in,ibb,iq) = mmn(im,in,ibb,iq)
     &       + dconjg(cnk(imp,im,iq)) * wmat(imp,in)
            enddo
         enddo
         enddo

      enddo
      enddo

      deallocate(wmat)

      return
      end
c-----------------------------------------------------------------------
      subroutine init_Umnk(amnk,cnk,
     i                     iko_ix,iko_fx,iko_i,iko_f,
     d                     nwf,nqbz,
     o                     umnk)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-3)
!      complex(8),allocatable :: amn(:,:),
!     &                          evecc(:,:),smat(:,:),wmat(:,:)
!      real (8),allocatable :: eval(:)
      complex(8),allocatable :: amn(:,:),cc(:,:),zz(:,:),vv(:,:)
      real (8),allocatable :: dd(:)
      complex(8) :: umnk(nwf,nwf,nqbz),
     &              amnk(iko_ix:iko_fx,nwf,nqbz),
     &              cnk(iko_ix:iko_fx,nwf,nqbz),ctmp
      integer(4) :: iko_i(nqbz),iko_f(nqbz)


      allocate (amn(nwf,nwf),zz(nwf,nwf),vv(nwf,nwf),dd(nwf))
      umnk = (0d0,0d0)

      do iq = 1,nqbz

c construct A
      amn = (0d0,0d0)
      do in = 1,nwf
      do im = 1,nwf
         do imp = iko_i(iq),iko_f(iq)
            amn(im,in) = amn(im,in) + 
     &       dconjg(cnk(imp,im,iq))*amnk(imp,in,iq)
         enddo
      enddo
      enddo

c singular value decomposition: amn = zz*dd*vv
      call zgesvdmn(nwf,nwf,amn,dd,zz,vv)

c U(m,n) = ( A S^{-1/2} )_mn = (zz*vv)_mn
      do im = 1,nwf
      do in = 1,nwf
         do ik = 1,nwf
           umnk(im,in,iq) = umnk(im,in,iq) + zz(im,ik)*vv(ik,in)
         enddo
      enddo
      enddo

      enddo ! iq

      deallocate(amn,zz,vv,dd)

      return
      end
c-----------------------------------------------------------------------
      subroutine updt_mmn(umnk,mmn0,ikbidx,
     d                    nwf,nqbz,nbb,
     u                    mmn)

      implicit real*8(a-h,o-z)

      parameter (eps = 1d-3)
      complex(8),allocatable :: wmat(:,:),wmat2(:,:)
      complex(8) :: umnk(nwf,nwf,nqbz),
     &              mmn(nwf,nwf,nbb,nqbz),
     &              mmn0(nwf,nwf,nbb,nqbz)
      integer (4) :: ikbidx(nbb,nqbz)

      allocate (wmat(nwf,nwf),wmat2(nwf,nwf))
      wmat = (0d0,0d0)

      do iq = 1,nqbz
      do ibb = 1,nbb
         iqb = ikbidx(ibb,iq)

c wmat = M * U
         wmat = (0d0,0d0)       
         do ii = 1,nwf
         do in = 1,nwf
            do ij = 1,nwf
               wmat(ii,in) = wmat(ii,in) + 
     &             mmn0(ii,ij,ibb,iq)*umnk(ij,in,iqb)
            enddo
         enddo
         enddo

c wmat2 = U^{-1} * wmat
         wmat2 = (0d0,0d0)       
         do im = 1,nwf
         do in = 1,nwf
            do ii = 1,nwf
               wmat2(im,in) = wmat2(im,in) + 
     &         dconjg(umnk(ii,im,iq)) * wmat(ii,in)
            enddo
         enddo
         enddo

c M = wmat2
         mmn(:,:,ibb,iq) = wmat2(:,:)

c end of ibb-loop
      enddo
c end of iq-loop
      enddo

      deallocate (wmat,wmat2)

      return
      end
c-----------------------------------------------------------------------
       subroutine get_rn(mmn,bb,wbb,wbz,
     d                   nwf,nqbz,nbb,
     o                   rn)

      implicit real*8(a-h,o-z)

      complex(8) :: mmn(nwf,nwf,nbb,nqbz)
      real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz),rn(3,nwf)

      rn = 0d0

      do in = 1,nwf
         do ibb = 1,nbb
            rtmp = 0d0
            do iq = 1,nqbz
                rtmp = rtmp
     &        + dimag(log(mmn(in,in,ibb,iq)))*wbz(iq)
            enddo
            rn(:,in) = rn(:,in) + wbb(ibb)*bb(:,ibb)*rtmp
c            write(*,*)in,ibb,rtmp
         enddo
      enddo

      rn(:,:) = - rn(:,:)

      return
      end
c-----------------------------------------------------------------------
      subroutine getrmn(mmn,
     d                  nwf,
     o                  rmn)

c [1] eq.(45)
      implicit real*8(a-h,o-z)

      complex(8) :: mmn(nwf,nwf),rmn(nwf,nwf)

      rmn = (0d0,0d0)

      do im = 1,nwf
      do in = 1,nwf
         rmn(im,in) = mmn(im,in) * dconjg(mmn(in,in))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine getamn(bmn,
     d                  nwf,
     o                  amn)

      implicit real*8(a-h,o-z)

      complex(8) :: bmn(nwf,nwf),amn(nwf,nwf)

      amn = (0d0,0d0)

      do im = 1,nwf
      do in = 1,nwf
         amn(im,in) = (bmn(im,in) - dconjg(bmn(in,im))) * 0.5d0
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine getsmn(bmn,
     d                  nwf,
     o                  smn)

      implicit real*8(a-h,o-z)

      complex(8) :: bmn(nwf,nwf),smn(nwf,nwf),ci

      smn = (0d0,0d0)
      ci = (0d0,1d0)

      do im = 1,nwf
      do in = 1,nwf
         smn(im,in) = -(bmn(im,in) + dconjg(bmn(in,im))) * ci * 0.5d0
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
       subroutine gettmn(rn,mmn,bb,
     d                  nwf,
     o                  qn,tmn)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: rtmn(:,:)
      complex(8) :: mmn(nwf,nwf),tmn(nwf,nwf)
c      real(8),allocatable :: qn(:)
      real(8) :: qn(nwf)
      real(8) :: rn(3,nwf),bb(3)

      allocate (rtmn(nwf,nwf))

c qn ([1] eq.47)
      do in = 1,nwf
         qn(in) = dimag(log(mmn(in,in))) + sum(bb(:)*rn(:,in))
      enddo

c R~mn ([1] eq.48)
      do im = 1,nwf
      do in = 1,nwf
         rtmn(im,in) = mmn(im,in) / mmn(in,in)
      enddo
      enddo

c Tmn = R~mn * qn
      do im = 1,nwf
      do in = 1,nwf
         tmn(im,in) = rtmn(im,in) * qn(in)
      enddo
      enddo

      deallocate(rtmn)

      return
      end
c-----------------------------------------------------------------------
      subroutine updt_uk(wmn,
     d                   nwf,
     u                   umn)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: iwmn(:,:),ewmn(:,:),evecc(:,:),
     &                          zmn(:,:),evalc(:)
      complex(8) :: wmn(nwf,nwf),umn(nwf,nwf),ci
      real(8),allocatable :: eval(:)

      allocate (iwmn(nwf,nwf),ewmn(nwf,nwf),
     &          evecc(nwf,nwf),zmn(nwf,nwf),evalc(nwf),
     &          eval(nwf))

      ci = (0d0,1d0)

c diagonalize W
c notice! W is anti-Hermitian, not Hermitian
      iwmn(:,:) = ci * wmn(:,:)

      call chk_hm(iwmn,nwf)
      call diag_hm(iwmn,nwf,eval,evecc)

c eW = exp(W) = exp(-i * iW)
      do ii = 1, nwf
c         write(*,*)'ev',ii,eval(ii)
          evalc(ii) = exp(-ci*eval(ii))
      enddo

      ewmn = (0d0,0d0)
      do im = 1,nwf
      do in = 1,nwf
         do il = 1,nwf
            ewmn(im,in) = ewmn(im,in) + 
     &       evecc(im,il) * evalc(il) * dconjg(evecc(in,il))
         enddo
      enddo
      enddo

      call chk_um(ewmn,nwf)

c U = U * exp(W)
      zmn = (0d0,0d0)
      do im = 1,nwf
      do in = 1,nwf
         do il = 1,nwf
            zmn(im,in) = zmn(im,in) + umn(im,il) * ewmn(il,in)
         enddo
      enddo
      enddo

      umn = zmn

      deallocate(iwmn,ewmn,evecc,zmn,evalc,eval)

      return
      end
c-----------------------------------------------------------------------
      subroutine getOmg(mmn,rn,bb,wbb,wbz,
     d                  nwf,nqbz,nbb,
     o                  omgi,omgd,omgod,omgdod,omgidod)

      implicit real*8(a-h,o-z)

      complex(8) :: mmn(nwf,nwf,nbb,nqbz)
      real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz),rn(3,nwf)

c omgi ([1] eq.34) 
      omgi = 0d0
      do iq = 1,nqbz
      do ibb = 1,nbb
         a2mmn = 0d0
         do im = 1,nwf
         do in = 1,nwf
            a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
         enddo
         enddo
         omgi = omgi + wbz(iq)*wbb(ibb)*(dble(nwf)-a2mmn)
      enddo
      enddo

c omgod ([1] eq.35) 
      omgod = 0d0
      do iq = 1,nqbz
      do ibb = 1,nbb
         a2mmn = 0d0
         do im = 1,nwf
         do in = 1,nwf
            if (im.ne.in)
     &      a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
         enddo
         enddo
         omgod = omgod + wbz(iq)*wbb(ibb)*a2mmn
      enddo
      enddo

c omgd ([1] eq.36) 
      omgd = 0d0
      do ibb = 1,nbb
         rtmp = 0d0
         do in = 1,nwf
            brn = sum(bb(:,ibb)*rn(:,in))
            do iq = 1,nqbz
               rtmp = rtmp + 
     &         wbz(iq)*(-dimag(log(mmn(in,in,ibb,iq)))-brn)**2
            enddo
         enddo
         omgd = omgd + wbb(ibb)*rtmp
      enddo

c sum up
      omgdod = omgd + omgod
      omgidod = omgi + omgdod

      return
      end
c-----------------------------------------------------------------------
      subroutine writeOmg(is,mmn,rn,bb,wbb,wbz,tpia,
     d                   nwf,nqbz,nbb)

      implicit real*8(a-h,o-z)

      complex(8) :: mmn(nwf,nwf,nbb,nqbz)
      real(8) :: bb(3,nbb),wbb(nbb),wbz(nqbz),rn(3,nwf),omg(nwf)


c omgi ([1] eq.34) 
      omgi = 0d0
      do iq = 1,nqbz
      do ibb = 1,nbb
         a2mmn = 0d0
         do im = 1,nwf
         do in = 1,nwf
            a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
         enddo
         enddo
         omgi = omgi + wbz(iq)*wbb(ibb)*(dble(nwf)-a2mmn)
      enddo
      enddo
c
      omg = 0d0
      do im = 1,nwf
        do iq = 1,nqbz
        do ibb = 1,nbb
         a2mmn = 0d0
         do in = 1,nwf
            a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
         enddo
         omg(im) = omg(im) + wbz(iq)*wbb(ibb)*(1d0-a2mmn)
        enddo
        enddo
      enddo

c omgod ([1] eq.35) 
      omgod = 0d0
      do iq = 1,nqbz
      do ibb = 1,nbb
         a2mmn = 0d0
         do im = 1,nwf
         do in = 1,nwf
            if (im.ne.in)
     &      a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
         enddo
         enddo
         omgod = omgod + wbz(iq)*wbb(ibb)*a2mmn
      enddo
      enddo
c
      do im = 1,nwf
        do iq = 1,nqbz
        do ibb = 1,nbb
         a2mmn = 0d0
         do in = 1,nwf
            if (im.ne.in)
     &      a2mmn = a2mmn + abs(mmn(im,in,ibb,iq))**2
         enddo
         omg(im) = omg(im) + wbz(iq)*wbb(ibb)*a2mmn
        enddo
        enddo
      enddo

c omgd ([1] eq.36) 
      omgd = 0d0
      do ibb = 1,nbb
         rtmp = 0d0
         do in = 1,nwf
            brn = sum(bb(:,ibb)*rn(:,in))
            do iq = 1,nqbz
               rtmp = rtmp + 
     &         wbz(iq)*(-dimag(log(mmn(in,in,ibb,iq)))-brn)**2
            enddo
         enddo
         omgd = omgd + wbb(ibb)*rtmp
      enddo
c
      do in = 1,nwf
        do ibb = 1,nbb
         rtmp = 0d0
            brn = sum(bb(:,ibb)*rn(:,in))
            do iq = 1,nqbz
               rtmp = rtmp + 
     &         wbz(iq)*(-dimag(log(mmn(in,in,ibb,iq)))-brn)**2
            enddo
         omg(in) = omg(in) + wbb(ibb)*rtmp
        enddo
      enddo

c multiply (a/2pi)^2
      omgd = omgd / tpia**2
      omgod = omgod / tpia**2
      omgi = omgi / tpia**2
      omgdod = omgdod / tpia**2
      omg = omg / tpia**2

c sum up
      omgdod = omgd + omgod
      omgidod = omgi + omgdod

c output
      if (is.eq.1) then
        ifbnd = iopen('spread.up',1,-1,0)
      else
        ifbnd = iopen('spread.dn',1,-1,0)
      endif

      write(ifbnd,*)'Total   omega(a.u.)     omega(ang.^2)'
      write(ifbnd,"(5x,2f16.8)")omgidod,omgidod*0.5292d0**2
      write(ifbnd,*)'band    omega(a.u.)     omega(ang.^2)'
      do in = 1,nwf
         write(ifbnd,"(i5,2f16.8)")in,omg(in),omg(in)*0.5292d0**2
      enddo
      omgall = sum(omg(:))
c      write(ifbnd,"(2x,a3,2f16.8)")'all',omgall,omgall*0.5292d0**2

      write(ifbnd,*)'(a/2pi)^2 corrected'
      if (is.eq.1) then
         isx = iclose('spread.up')
      else
         isx = iclose('spread.dn')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine writermn(is,mmn,bb,wbb,qbz,qbz0,wbz,rt,
     d                   nwf,nqbz,nbb,n1,n2,n3)

      implicit real*8(a-h,o-z)

      complex(8) :: mmn(nwf,nwf,nbb,nqbz),ci,czero,cikr,ceikr
     &          ,ctmp1,ctmp2
     &          ,rmn(3,nwf,nwf,nqbz),amn(3,nwf,nwf,nqbz)
     &          ,c1,c2,c3
      real(8) :: bb(3,nbb),wbb(nbb),qbz(3,nqbz),qbz0(3,nqbz),wbz(nqbz),
     &          rt(3,nqbz),dmn,rtmp(3)

      pi = 4d0*atan(1d0)
      ci = (0d0,1d0)
      czero = (0d0,0d0)

c Berry connecrion (in the Wannier gauge)
      amn = czero
      do im = 1,nwf
      do in = 1,nwf
         dmn = 0d0
         if (im.eq.in) dmn = 1d0
         do iq = 1,nqbz
         do ibb = 1,nbb
c            ctmp1 = (mmn(im,in,ibb,iq) - dmn)*wbb(ibb)
            r1 = dimag(0.5d0*(mmn(im,in,ibb,iq)+mmn(in,im,ibb,iq))-dmn)
            r2 = dreal(0.5d0*(mmn(im,in,ibb,iq)-mmn(in,im,ibb,iq)))
            ctmp1 = dcmplx(r2,r1)
            ctmp1 = log(ctmp1+1d0)*wbb(ibb)
            do ix = 1,3 ! x,y,z
               amn(ix,im,in,iq) = amn(ix,im,in,iq)
     &                          + ci*imag(ctmp1) * bb(ix,ibb)
            enddo ! ix
         enddo ! ibb
         enddo ! iq
      enddo ! in
      enddo ! im
      amn = amn * ci

c <0m | r | Rn>
      rmn = czero
      do ir = 1,nqbz
      do iq = 1,nqbz
c         rk = sum(rt(:,ir)*qbz(:,iq))
         rk = sum(rt(:,ir)*qbz0(:,iq))
         cikr = -ci * 2d0 * pi * rk
         ceikr = exp(cikr) * wbz(iq)
         do in = 1,nwf
         do im = 1,nwf
         do ix = 1,3
            rmn(ix,im,in,ir) = rmn(ix,im,in,ir) +
     &      ceikr * amn(ix,im,in,iq)
         enddo ! ix
         enddo ! im
         enddo ! in
      enddo ! iq
      enddo ! ir

c output
      if (is.eq.1) then
        ifrmn = iopen('rmn.up',1,-1,0)
      else
        ifrmn = iopen('rmn.dn',1,-1,0)
      endif

      write(ifrmn,*)'*** nwf,nsite'
      write(ifrmn,*)nwf,nqbz
      write(ifrmn,*)'*** rsite'
      write(ifrmn,*)rt
      write(ifrmn,*)'*** rmn'
      write(ifrmn,*)rmn
      write(ifrmn,*)'*** amn'
      write(ifrmn,*)amn

      if (is.eq.1) then
         isx = iclose('rmn.up')
      else
         isx = iclose('rmn.dn')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine writemmn(is,mmn,bb,wbb,qbz,wbz,rt,
     d                   nwf,nqbz,nbb,n1,n2,n3)

      implicit real*8(a-h,o-z)

      complex(8) :: mmn(nwf,nwf,nbb,nqbz),ci,czero,cikr,ceikr
     &          ,ctmp1,ctmp2
     &          ,rmn(3,nwf,nwf,nqbz),amn(3,nwf,nwf,nqbz)
     &          ,c1,c2,c3
      real(8) :: bb(3,nbb),wbb(nbb),qbz(3,nqbz),wbz(nqbz),
     &          rt(3,nqbz),dmn

      pi = 4d0*atan(1d0)
      ci = (0d0,1d0)
      czero = (0d0,0d0)

c output
      if (is.eq.1) then
        ifrmn = iopen('mmn.up',1,-1,0)
      else
        ifrmn = iopen('mmn.dn',1,-1,0)
      endif

      write(ifrmn,*)'*** nwf,nsite,nb'
      write(ifrmn,*)nwf,nqbz,nbb
      write(ifrmn,*)'*** mmn'
      write(ifrmn,*)mmn
      write(ifrmn,*)'*** bb'
      write(ifrmn,*)bb
      write(ifrmn,*)'*** wbb'
      write(ifrmn,*)wbb
      write(ifrmn,*)'*** qbz'
      write(ifrmn,*)qbz
      write(ifrmn,*)'*** wbz'
      write(ifrmn,*)wbz

      if (is.eq.1) then
         isx = iclose('mmn.up')
      else
         isx = iclose('mmn.dn')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine wmaxloc(ifmlw,ifmlwe,
     i                   qbz,umnk,cnk,eunk,
     i                   iko_ix,iko_fx,iko_i,iko_f,
     d                   nwf,nqbz,nband,nlmto)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: dnk(:,:,:)
      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),
     &              umnk(nwf,nwf,nqbz)
      real(8) :: qbz(3,nqbz),eunk(nwf,nqbz)
      integer(4) :: iko_i(nqbz),iko_f(nqbz)

      allocate(dnk(iko_ix:iko_fx,nwf,nqbz))

      write(ifmlw)nqbz,nwf,iko_ix,iko_fx
c      write(ifmlwe)nqbz,nwf,iko_ix,iko_fx

      dnk = (0d0,0d0)
      do iq = 1,nqbz

         do imp = iko_i(iq),iko_f(iq)
         do in = 1,nwf
            do im = 1,nwf
               dnk(imp,in,iq) = dnk(imp,in,iq) 
     &                 + umnk(im,in,iq) * cnk(imp,im,iq)
            enddo
         enddo
         enddo

         write(ifmlw)iq,qbz(1:3,iq)
         write(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iq)
c         write(ifmlwe)iq,qbz(1:3,iq)
c         write(ifmlwe)eunk(1:nwf,iq)

      enddo

c      call chk_cnkweight(qbz,iko_ix,iko_fx,dnk,
c     &     nqbz,nwf,nband,nlmto)

      deallocate(dnk)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_dnk(is,eunk,qbz,
     i                   umnk,cnk,
     i                   iko_ix,iko_fx,iko_i,iko_f,
     d                   nband,nwf,nqbz)

      use m_readeigen

      implicit real*8(a-h,o-z)

      parameter (eps=1d-3)
      complex(8),allocatable :: dnk(:,:)
      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),
     &              umnk(nwf,nwf,nqbz),
     &              ctmp,ham(nwf,nwf),evecc(nwf,nwf)
      real(8) :: eunk(nwf,nqbz),qbz(3,nqbz),eks(nband),eval(nwf)
      integer(4) :: iko_i(nqbz),iko_f(nqbz)

      allocate(dnk(iko_ix:iko_fx,nwf))

      do iq = 1,nqbz

         call readeval (qbz(:,iq),is,eks)

         dnk = (0d0,0d0)
         do imp = iko_i(iq),iko_f(iq)
         do in = 1,nwf
            do im = 1,nwf
               dnk(imp,in) = dnk(imp,in) 
     &                 + umnk(im,in,iq) * cnk(imp,im,iq)
            enddo
         enddo
         enddo


         ham = 0d0
         do ii=1,nwf
         do ij=1,nwf
            ctmp = 0d0
            do ik = iko_i(iq),iko_f(iq)
               ctmp = ctmp + eks(ik) *
     &                  dconjg(dnk(ik,ii))*dnk(ik,ij)
            enddo
            ham(ii,ij) = ctmp
            write(97,"(3i5,2f12.6)")iq,ii,ij,dreal(ctmp),dimag(ctmp)
         enddo
         enddo
         call diag_hm(ham,nwf,eval,evecc)
         do ii = 1,nwf
            write(98,"(2i5,2f12.6)")iq,ii,eval(ii),eunk(ii,iq)
         enddo

      enddo

      deallocate(dnk)

      return
      end
c-----------------------------------------------------------------------
