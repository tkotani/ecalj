!========================================================================
      program htbplot
!! Tight binding band plot.
!! hrotr is the main input
c References
c [1] N. Marzari and D.Vanderbilt, PRB56,12847(1997)
c [2] I. Souza, N. Marzari and D.Vanderbilt, PRB65,035109(2002)
c      use m_readqg
c      use m_readeigen
c      use m_read_bzdata,ngrp2=>ngrp
c      use m_genallcf_v3
c      use keyvalue
      implicit none
      integer(4):: nlinex,ntmp
      parameter (nlinex=100)
      integer(4)::nline,np(nlinex)
      real(8):: qi(3,nlinex),qf(3,nlinex)

      complex(8),allocatable:: hrotr(:,:,:), hrotkp(:,:),evecc(:,:)
      integer:: iftb

      real(8),allocatable :: rws(:,:),drws(:),eval(:)
      integer,allocatable:: irws(:)
      integer:: nrws

      real(8)::qold(3)
      real(8),allocatable:: xq(:),eval3(:,:),pos(:,:),q(:,:)
      integer::npin 
      real(8):: qiin(3),qfin(3)

      integer, allocatable:: 
     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:),iposwf(:)
      real(8):: rydberg=13.6058d0, alat,ef,plat(3,3),rcut,rcuta
      integer:: ifile_handle,iband,ifh,n1,n2,n3,natom,i,iq,j,nq,nwf,is,if99
c---------------------------------------
c      hartree=2d0*rydberg()
c      iii=verbose()
c      write(6,*)' verbose=',iii
c      pi   = 4d0*datan(1d0)
c      tpia = 2d0*pi/alat
c      call dinv33(plat,1,xxx,vol)
c      voltot = dabs(vol)*(alat**3)
      write(6,*)' start htbplot: read range cutoff rcut(ang)=?'
      read(5,*) rcuta
      write(5,*)' rcut (ang)=', rcuta
      rcut=rcuta/0.529177
!! QPNT data, 
      write(*,*)'---Read k points for bands from SYML ---'
      if99=ifile_handle()
      open(if99,file='SYML',status='old') !symmetry line
      nline=0
      do i = 1,nlinex
        read(if99,*,err=551,end=552)npin,qiin,qfin
        if (npin==0) exit
        nline = nline+1
        np(nline)=npin
        qi(1:3,nline)=qiin  !cartesian (2pi/alat unit).
        qf(1:3,nline)=qfin
 551    continue
      enddo
 552  continue
      if (nline>nlinex) stop 'hmaxloc: too many lines in SYML'
      close(if99)
      nq = 0
      do i = 1,nline
        nq = nq + np(i)
      enddo
      allocate(q(3,nq),xq(nq))
      iq = 0
      xq=0d0
      qold=q(:,1)
      do i = 1,nline
      do j = 0,np(i)-1
         iq = iq + 1
         q(:,iq) = qi(:,i) + (qf(:,i)-qi(:,i))*dble(j)/dble(np(i)-1)
         if(iq>1) then
            xq(iq)= xq(iq-1) + dsqrt( sum((q(:,iq)-qold)**2) )
         endif
         qold=q(:,iq)
      enddo
      enddo

!! --- Required data set ---------------------
!! q(:,1:nq) 
!! ef: fermi energy or VBM
!! alat: unit for primitive vector, atomic positions (in a.u.)
!! plat: primitive vector
!! rcut: real-space cutoff for tb parameters.
!! nwf:  # of Wannier function
!! nrws: # of R-R' pair
!! rws(3,i): R-R' vector, i=1,nrws
!! irws(i) : degeneracy, i=1,nrws
!! hrotr(nwf,nwf,nrws) = Hmn(R) = <0m|H |Rn> (m=1,nwf; n=1,nwf, R=1,nrws)
!! 
!! natom: number of atoms in the primitive cell.
!! ibaswf(nwf): atomic position for the Wannier orbital.
!! pos(3,natom): atomic poisition in the cell.
      print *,' --- Reading HrotRS ---'
      is=1
      ifh=ifile_handle()
      if(is==1) open(ifh,file='HrotRS.up',form='unformatted')
      if(is==2) open(ifh,file='HrotRS.dn',form='unformatted')
      read(ifh)alat,plat,natom
      allocate(pos(3,natom))
      read(ifh)pos
      read(ifh)ef 
      read(ifh)nwf,nrws,n1,n2,n3
      print *,' nwf,nrws,n1,n2,n3=',nwf,nrws,n1,n2,n3
      allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),iposwf(nwf))
      allocate(hrotr(nwf,nwf,nrws)) !real space Hamiltonian in Wannier funciton basis
      read(ifh) irws,rws,hrotr, iposwf
      close(ifh)
      allocate(eval3(nwf,nq),evecc(nwf,nwf),eval(nwf)) 

!! hrotkp= Hrot_mn(q)  -- Tight-binding Hamiltonian at q(1:3,iq) ---
!! Conversion from hrotr --> hrotkp by FT.
      allocate(hrotkp(nwf,nwf))
      do iq = 1,nq
        write(6,"(a,i5,3f11.5)")' --- get hrotkp and diagonalize --- iq=',iq,q(:,iq)
        call get_hrotkp_tb_ws2(rcut,alat,
     i         hrotr,rws,irws,q(:,iq), iposwf,pos,natom,
     d         nwf,nrws,
     o         hrotkp) 
        call diag_hm(hrotkp,nwf,eval,evecc) ! diagonalization 
        eval3(1:nwf,iq)=eval
      enddo

!! Write band data
      iftb=ifile_handle()
      if(is==1) open(iftb,file='bnds.tb.up')         
      if(is==2) open(iftb,file='bnds.tb.dn')         
      write(iftb,*)nq
      write(iftb,*)nwf
      write(iftb,*)ef,rcuta,' ! ef, rcut(angstrom)'
      do iband = 1,nwf
        do iq = 1,nq
          write(iftb,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' !eee! x eval-ef(ev) iband' )") 
     &     iq,q(1:3,iq),  xq(iq),(eval3(iband,iq)-ef)*rydberg,iband
        enddo
        write(iftb,*)
      enddo
      stop 'OK! end of htpplot'
      end

!-------------------------------------------------------------------------
      subroutine get_hrotkp_tb_ws2(rcut,alat,
     i                     hrotr,rws,irws,q,  iposwf,pos,natom,
     d                     nwf,nrws,
     o                     hrotkp)
c truncate long-range part of hrotr
c from get_hrotkp_ws
c see Ref.[2] eq.26  <0m | H| Rn>
c      implicit real*8(a-h,o-z)
      implicit none
      real(8),parameter::delta = 1d-3
      integer:: nrws,nwf
      complex(8) :: hrotr(nwf,nwf,nrws),hrotkp(nwf,nwf),
     &              ci,cikr,ceikr,ctmp
      real(8) :: q(3),rws(3,nrws),alat
      integer(4) :: irws(nrws)
      integer:: iposwf(nwf),natom,im,in,inp,imp,ir
      real(8):: pos(3,natom),rc,pi,rcut,rk,rtmp
c      print *,'get_hrotkp_tb_ws2: rcut nef nrws=',rcut,nwf,nrws
      rc = rcut
      pi = 4d0* atan(1d0)
      ci = (0d0,1d0)
      hrotkp = (0d0,0d0)
      do ir = 1,nrws
         ceikr = (0d0,0d0)
c         rtmp = alat*dsqrt(sum(rws(:,ir)**2))
c         if (rtmp.le.rc) then
            rk = sum(rws(:,ir)*q(:))
            cikr = ci * 2d0 * pi * rk
            ceikr = exp(cikr) / dble(irws(ir))
c         endif  
         do im = 1,nwf
           imp= iposwf(im)
         do in = 1,nwf
           inp= iposwf(in)
!! Rn -0m
           rtmp = alat*dsqrt( sum( (rws(:,ir)+pos(:,inp)-pos(:,imp))**2 ) )
           if(rtmp<rc) then
             hrotkp(im,in) = hrotkp(im,in) +  ceikr * hrotr(im,in,ir)
           endif 
         enddo
         enddo
      enddo
      return
      end

!-------------------------------------------------------------------------
      subroutine diag_hm(zmat,ndim,eval,evecc)
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
      complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
      complex(8):: zmat(ndim,ndim),evecc(ndim,ndim)
      real(8):: eval(ndim),wk(ndim,11)
      integer iwk(ndim)
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
!      call diagno(ndim,zmat2,ovlpc,wk,iwk,evecc,eval)
      call diagcv(ovlpc,zmat2, evecc, ndim, eval, nmx, 1d99, nev)
      deallocate(zmat2,ovlpc)
      return
      end
!-------------------------------------------------------------------------
      subroutine diagcv(s,h,t,n,evl,nmx,emx,nev)
C  diagonalizes and returns the nev lowest eigenstates
C  eigenvecs are returned for i.lt.nmx. emx is dummy
      implicit none
      integer n,nmx,nev
      double precision emx,evl(n),abstol,dlamch
      complex*16 s(n,n),h(n,n),t(n,n)
      logical lx,lnv
      integer i,ipr,iprint,j,oww
      complex(8),allocatable:: work(:),s_(:,:),h_(:,:)
      integer(4),allocatable:: iwork(:),ifail(:)
      real(8),allocatable:: rwork(:)
      real(8):: vl,vu  !bugfix (whis was not decleared) jan2013
      integer:: lwork,info
c
c      print *,'diagcv in diagcv2:'
      abstol= 2d0*DLAMCH('S')
c      print *,' xxx abstol=',abstol
c---find optimum work size
      allocate( work(1))
      LWORK = -1
      call ZHEGVX( 1, 'V', 'I', 'U', n, h, n, s, n,
     $                   VL, VU, 1, nmx, abstol, nev, evl, t, n, 
     $                   work,LWORK, RWORK, IWORK, IFAIL, INFO )
      lwork = work(1) + 1
      deallocate(work)

c      LWORK = max(1,2*n-1)    ! This caused a  problem---why?
c      LWORK = max(1,2*n-1)+1  ! This caused no problem---why?

      allocate( work(LWORK),rwork(7*n),iwork(5*n),ifail(n))
c      print *,' goto zhegvx'
      call ZHEGVX( 1, 'V', 'I', 'U', n, h, n, s, n,
     $                   VL, VU, 1, nmx, abstol, nev, evl, t, n, 
     $                   work,LWORK, RWORK, IWORK, IFAIL, INFO )
      if(INFO/=0) then
        print *, 'See http://www.netlib.org/cgi-bin/netlibget.pl',
     &   '/lapack/complex16/zhegvx.f'
        print *, 'ZHEGVX error info=',info
        print *, 'ZHEGVX error info=',ifail
      endif
      deallocate(work,rwork,iwork,ifail)
      return
      end

      integer function ifile_handle()
!! find open file handle
      implicit none
      integer:: i
      logical:: nexist
      do i=5001,9999
         inquire(unit=i,opened=nexist)
         if(.not.nexist) then
            ifile_handle=i
            return
         endif
      enddo
      stop 'ifile_handle: we did not find open file hundle'
      end
