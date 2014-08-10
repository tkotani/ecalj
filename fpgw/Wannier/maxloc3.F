      subroutine get_nqbze(nqbz,nqbze)
      implicit real*8(a-h,o-z)

      logical :: lq0p

      inquire(file='Q0P',exist=lq0p)

      if (lq0p) then
         ifq0p = iopen('Q0P',1,0,0)
         read (ifq0p,"(i5)") nq0i
         ifq0P = iclose('Q0P')         
         nqbze = nqbz*(nq0i+1)
      else
         nqbze = nqbz
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_nqbze2(nqbz,nqbze)
      implicit real*8(a-h,o-z)

      logical :: lq0p

      inquire(file='Q0P',exist=lq0p)

      if (lq0p) then
         ifq0p = iopen('Q0P',1,0,0)
         read (ifq0p,"(i5)") nq0i
         ifq0P = iclose('Q0P')         
         nqbze = nqbz*(nq0i+2)
      else
         nqbze = nqbz*2
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_qbze(qbz,nqbz,
     o                    qbze,nqbze)
      implicit real*8(a-h,o-z)

      logical :: lq0p
      real(8), allocatable :: wqt(:),q0i(:,:)
      real(8) :: qbz(3,nqbz),qbze(3,nqbze)

      inquire(file='Q0P',exist=lq0p)

      if (lq0p) then
         ifq0p = iopen('Q0P',1,0,0)
         read (ifq0p,"(i5)") nq0i
         allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
         do iq=1,nq0i
           read (ifq0p, * ) wqt(iq),q0i(1:3,iq)
         enddo
         ifq0P = iclose('Q0P')         

c         nqbze = nqbz*(nq0i+1)
c         allocate(qbze(3,nqbze))
         qbze(:,1:nqbz) = qbz(:,1:nqbz)
         do iq = 1,nq0i
           ini = nqbz*(1 + iq -1)
           do ix=1,nqbz
             qbze (:,ini+ix)   = q0i(:,iq) + qbz(:,ix)
           enddo
         enddo
         deallocate(wqt,q0i)
      else
c         nqbze = nqbz
c         allocate(qbze(3,nqbze))
         qbze = qbz
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_qbze2(qbz,qbzs,nqbz,
     o                    qbze,nqbze)
      implicit real*8(a-h,o-z)

      logical :: lq0p
      real(8), allocatable :: wqt(:),q0i(:,:)
      real(8) :: qbz(3,nqbz),qbzs(3,nqbz),qbze(3,nqbze)

      inquire(file='Q0P',exist=lq0p)

      if (lq0p) then
         ifq0p = iopen('Q0P',1,0,0)
         read (ifq0p,"(i5)") nq0i
         allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
         do iq=1,nq0i
           read (ifq0p, * ) wqt(iq),q0i(1:3,iq)
         enddo
         ifq0P = iclose('Q0P')         

c         nqbze = nqbz*(nq0i+2)
c         allocate(qbze(3,nqbze))
         qbze(:,1:nqbz) = qbz(:,1:nqbz)
         qbze(:,nqbz+1:2*nqbz) = qbzs(:,1:nqbz)
         do iq = 1,nq0i
           ini = nqbz*(2 + iq -1)
           do ix=1,nqbz
             qbze (:,ini+ix)   = q0i(:,iq) + qbz(:,ix)
           enddo
         enddo
         deallocate(wqt,q0i)
      else
c         nqbze = nqbz*2
c         allocate(qbze(3,nqbze))
         qbze(:,1:nqbz) = qbz(:,1:nqbz)
         qbze(:,nqbz+1:2*nqbz) = qbzs(:,1:nqbz)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rot_hmnk(umnk,eunk,
     d                    nwf,nqbz,
     o                    hrotk)

c see Ref.[2] eq.24

      implicit real*8(a-h,o-z)

      complex(8) :: umnk(nwf,nwf,nqbz),hrotk(nwf,nwf,nqbz),ctmp
      real(8) :: eunk(nwf,nqbz)

      hrotk = (0d0,0d0)

      do iq = 1,nqbz
      do im = 1,nwf
      do in = 1,nwf
         do ii = 1,nwf
            hrotk(im,in,iq) = hrotk(im,in,iq) + 
     &      dconjg(umnk(ii,im,iq))*eunk(ii,iq)*umnk(ii,in,iq)
         enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_hrotr_ws(hrotk,qbz,wbz,
     i                  rws,irws,drws,
     d                  nwf,nqbz,nrws,
     o                  hrotr)

c see Ref.[2] eq.25

      implicit real*8(a-h,o-z)

      complex(8) :: hrotk(nwf,nwf,nqbz),hrotr(nwf,nwf,nrws),
     &              ci,cikr,ceikr,ctmp
      real(8) :: qbz(3,nqbz),wbz(nqbz),q(3),r(3)
     &          ,rws(3,nrws),drws(nrws)
      integer(4) :: irws(nrws)
 
      pi = 4d0* atan(1d0)
      ci = (0d0,1d0)

      hrotr = (0d0,0d0)

      ir = 0
      do ir = 1,nrws
      do iq = 1,nqbz
         rk = 2d0*pi*sum(rws(:,ir)*qbz(:,iq))
         ceikr = exp(-ci*rk)
         do im = 1,nwf
         do in = 1,nwf
            hrotr(im,in,ir) = hrotr(im,in,ir) + 
     &      ceikr * hrotk(im,in,iq) / dble(nqbz)
         enddo ! in
         enddo ! im
      enddo ! iq
      enddo ! ir

      return
      end
c-----------------------------------------------------------------------
      subroutine get_hrotkp_ws(hrotr,rws,drws,irws,q,
     d                     nwf,nqbz,nrws,
     o                     hrotkp)

c see Ref.[2] eq.26

      implicit real*8(a-h,o-z)

      complex(8) :: hrotr(nwf,nwf,nrws),hrotkp(nwf,nwf),
     &              ci,cikr,ceikr,ctmp
      real(8) :: q(3),rws(3,nrws),drws(nrws)
      integer(4) :: irws(nrws)

      pi = 4d0* atan(1d0)
      ci = (0d0,1d0)

      hrotkp = (0d0,0d0)

      do ir = 1,nrws
         rk = sum(rws(:,ir)*q(:))
         cikr = ci * 2d0 * pi * rk
         ceikr = exp(cikr) / dble(irws(ir))
         do im = 1,nwf
         do in = 1,nwf
            hrotkp(im,in) = hrotkp(im,in) + 
     &      ceikr * hrotr(im,in,ir)
         enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_hrotkp_tb_ws(rcut,plat,alat,
     i                     hrotr,rws,drws,irws,q,
     d                     nwf,nqbz,nrws,
     o                     hrotkp)

c truncate long-range part of hrotr
c from get_hrotkp_ws
c see Ref.[2] eq.26

      implicit real*8(a-h,o-z)

      parameter (delta = 1d-3)
      complex(8) :: hrotr(nwf,nwf,nrws),hrotkp(nwf,nwf),
     &              ci,cikr,ceikr,ctmp
      real(8) :: q(3),rws(3,nrws),drws(nrws),plat(3,3)
      integer(4) :: irws(nrws)

ctm
c      rc = rcut + delta
      rc = rcut

      pi = 4d0* atan(1d0)
      ci = (0d0,1d0)

      hrotkp = (0d0,0d0)

      do ir = 1,nrws
         ceikr = (0d0,0d0)
         rtmp = alat*dsqrt(sum(rws(:,ir)**2))
         if (rtmp.le.rc) then
            rk = sum(rws(:,ir)*q(:))
            cikr = ci * 2d0 * pi * rk
            ceikr = exp(cikr) / dble(irws(ir))
         endif  
         do im = 1,nwf
         do in = 1,nwf
            hrotkp(im,in) = hrotkp(im,in) + 
     &      ceikr * hrotr(im,in,ir)
         enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wmaxloc_diag(ifmlw,ifmlwe,
     i                 iq,q,umnk,cnk,eunk,evecc,eval,
     i                 iko_ix,iko_fx,iko_i,iko_f,
     d                 nwf,nqbz)

      implicit real*8(a-h,o-z)

      complex(8),allocatable :: dnk(:,:)
      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),
     &              umnk(nwf,nwf,nqbz),evecc(nwf,nwf)
      real(8) :: q(3),eunk(nwf,nqbz),eval(nwf)
      integer(4) :: iko_i(nqbz),iko_f(nqbz)


      write(ifmlw)iq,q(1:3)
      write(ifmlw)evecc(1:nwf,1:nwf)

      write(ifmlwe)iq,q(1:3)
      write(ifmlwe)eval(1:nwf)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_diag(q,ham,nwf,eval,evec)
      implicit real*8(a-h,o-z)

      complex(8) :: ham(nwf,nwf),evec(nwf,nwf),zm1(nwf,nwf),zm2(nwf,nwf)
      real(8) :: eval(nwf),q(3)


      write(7100,*)'***'
      write(7100,"(3f12.6)")q

c zm1(i1,iwf2)   = S[i2] ham(i1,i2) * evec(i2,iwf2)
      zm1 = 0d0
      do iwf2 = 1,nwf
      do i1   = 1,nwf
         do i2 = 1,nwf
            zm1(i1,iwf2) = zm1(i1,iwf2) + ham(i1,i2) * evec(i2,iwf2)
         enddo
      enddo
      enddo

c zm2(iwf1,iwf2) = S[i1] conjg(evec(i1,iwf1)) * zm1(i1,iwf2)
      zm2 = 0d0
      do iwf1 = 1,nwf
      do iwf2 = 1,nwf
         do i1 = 1,nwf
            zm2(iwf1,iwf2) = zm2(iwf1,iwf2)
     &              + dconjg(evec(i1,iwf1)) * zm1(i1,iwf2)
         enddo
      enddo
      enddo

c output
      do iwf1 = 1,nwf
      do iwf2 = 1,nwf
         tmp = 0d0
         if (iwf1.eq.iwf2) tmp = eval(iwf1)
         write(7100,"(2i5,3f12.6)")iwf1,iwf2,
     &        dreal(zm2(iwf1,iwf2)),dimag(zm2(iwf1,iwf2)),tmp
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_umnk(q,ham,nwf,eval,umn)
      implicit real*8(a-h,o-z)

      complex(8) :: ham(nwf,nwf),evec(nwf,nwf),zm1(nwf,nwf),zm2(nwf,nwf),
     &              umn(nwf,nwf)
      real(8) :: eval(nwf),q(3)


      do i1 = 1,nwf
      do i2 = 1,nwf
         evec(i1,i2) = dconjg(umn(i2,i1))
      enddo
      enddo

      write(7000,*)'***'
      write(7000,"(3f12.6)")q

c zm1(i1,iwf2)   = S[i2] ham(i1,i2) * evec(i2,iwf2)
      zm1 = 0d0
      do iwf2 = 1,nwf
      do i1   = 1,nwf
         do i2 = 1,nwf
            zm1(i1,iwf2) = zm1(i1,iwf2) + ham(i1,i2) * evec(i2,iwf2)
         enddo
      enddo
      enddo

c zm2(iwf1,iwf2) = S[i1] conjg(evec(i1,iwf1)) * zm1(i1,iwf2)
      zm2 = 0d0
      do iwf1 = 1,nwf
      do iwf2 = 1,nwf
         do i1 = 1,nwf
            zm2(iwf1,iwf2) = zm2(iwf1,iwf2)
     &              + dconjg(evec(i1,iwf1)) * zm1(i1,iwf2)
         enddo
      enddo
      enddo

c output
      do iwf1 = 1,nwf
      do iwf2 = 1,nwf
         tmp = 0d0
         if (iwf1.eq.iwf2) tmp = eval(iwf1)
         write(7000,"(2i5,3f12.6)")iwf1,iwf2,
     &        dreal(zm2(iwf1,iwf2)),dimag(zm2(iwf1,iwf2)),tmp
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_umn_evec(q,umn,evec,eval,nwf)
      implicit real*8(a-h,o-z)

      complex(8) :: evec(nwf,nwf),zm1(nwf,nwf),zm2(nwf,nwf),
     &              umn(nwf,nwf)
      real(8) :: eval(nwf),q(3)


      do i1 = 1,nwf
      do i2 = 1,nwf
         zm1(i1,i2) = dconjg(umn(i2,i1))
      enddo
      enddo

      zm2 = zm1 - evec

c output
      write(7000,*)'***'
      write(7000,"(3f12.6)")q
      do iwf = 1,nwf
         write(7000,*)iwf,eval(iwf)
      enddo
      do iwf1 = 1,nwf
      do iwf2 = 1,nwf
         write(7000,"(2i5,3f12.6)")iwf1,iwf2,
     &        dreal(zm2(iwf1,iwf2)),dimag(zm2(iwf1,iwf2))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine writeham(ifi,is,ef,alat,plat,pos,qbz,wbz,rws,irws,hrotk,
     &              nspin,natom,nwf,nqbz,nrws)
      implicit none

      integer(4) :: nspin,nwf,natom,nqbz,nrws,is,irws(nrws)
      double precision :: ef,alat,plat(3,3),pos(3,natom),qbz(3,nqbz),
     &             wbz(nqbz),rws(3,nrws)
      complex(8) :: hrotk(nwf,nwf,nqbz)

      integer(4) :: iopen,iclose,ifi

      if (is.eq.1) then
         ifi = iopen('HMLWF',1,-1,0)
         write(ifi,*)nspin,natom,nwf,nqbz,nrws,ef
         write(ifi,*)alat
         write(ifi,*)plat
         write(ifi,*)pos
         write(ifi,*)qbz
         write(ifi,*)wbz
         write(ifi,*)rws
         write(ifi,*)irws
      endif

      write(ifi,*)hrotk

      if (is.eq.nspin) ifi = iclose('HMLWF')

      return
      end
c-----------------------------------------------------------------------
