      subroutine checkeq(i,j)
Cstop2rx 2013.08.09 kino      if(i/=j) stop 
      if(i/=j) call rx( ''//
     & " checkeq in hsfp0: dim of WVR and WVI not compatible")
      end

c--------------------------------------------------------------------
      subroutine wsex2 (ifsex,q,itq,eq,sexq,
     d                 ntq)
c 92.02.28
c 92.04.20 from wsec, write SEx or vxc

c ifsex = unit file
c q     = q-vector
c itq(ntq) = state label
c eq    = eigenvalues at q
c sexq  = SEx(q,t)
c ntq   = no. states t

      implicit real*8 (a-h,o-z)
      dimension q(3),itq(ntq),eq(ntq),sexq(ntq)
      data itwo/2/

      write (ifsex,*)' '
      write (ifsex,*)' ***'
      write (ifsex,6600)ntq,itwo,(q(i),i=1,3)
      do       i = 1,ntq
        write (ifsex,6100)itq(i),eq(i),sexq(i)
ctakao for debug
        write (6,'(i5,2d15.6)')itq(i),eq(i),sexq(i)
c
      end do
      write (ifsex,*)
     .'-------------------------------------------------------------'
c formats
 6100 format (1x,i4,2f10.4)
 6600 format (1x,2i4,3f10.4,'  q-vector')
      end

c--------------------------------------------------------------------
c      subroutine wsec2 (ifsec,is,q,itq,eqx,deltaw,nbloch,rsec,csec,
c     d                 nw,ntq,nlmto)
c
c 92.02.28
c write SEc(qt,w)
c
c ifsec = unit file
c is    = spin label
c q     = q-vector
c itq(ntq) = state label
c eq    = eigenvalues at q
c deltaw= energy mesh in SEc(qt,w)
c nbloch= no. optimal product basis
c rsec,csec(nw,ntq) = real and imaginary part of SEc(qt,w)
c nw    = no. w in SEc(qt,w)
c ntq   = no. states t
c nlmto = no. LMTO-basis
c
c      implicit real*8 (a-h,o-z)
c      dimension q(3),itq(ntq),eqx(ntq),rsec(nw,ntq),csec(nw,ntq)
c      real(8) ::  rydberg,hartree
c      data ifive,isix/5,6/
c     hartree=2d0*rydberg()
c
c     write (ifsec,*)' '
c     write (ifsec,*) 'Self-energy correlated SEc(qt,w)'
c     write (ifsec,6500) is,nbloch,deltaw
c
c      write (ifsec,*)' '
c      write (ifsec,6400) q(1),q(2),q(3)
c      write (ifsec,*)
c     .'states,e(qt), Re(SEc(qt,w)) eV and Z'
c      write (ifsec,*)' ***'
c      write (ifsec,6600)ntq,isix,(q(i),i=1,3)
c      do       i = 1,ntq
c      slope      = (rsec(3,i) - rsec(1,i))/(2.d0*deltaw)
c      slope      = 1.d0/(1.d0-slope)
c      elda       =  eqx(i)  !rydberg()*(eq(itq(i))-ef)
c      rse1       = hartree*rsec(1,i)
c      rse2       = hartree*rsec(2,i)
c      rse3       = hartree*rsec(3,i)
c      write (ifsec,6100)itq(i),elda,rse1,rse2,rse3,slope
c      end do
c
c      write (ifsec,*)' '
c      write (ifsec,*)
c     .'states,e(qt), Im(SEc(qt,w)) eV '
c      write (ifsec,*)' ***'
c      write (ifsec,6600)ntq,ifive,(q(i),i=1,3)
c      do       i = 1,ntq
c      elda       =  eqx(i) ! rydberg()*(eq(itq(i))-ef)
c      cse1       = hartree*csec(1,i)
c      cse2       = hartree*csec(2,i)
c      cse3       = hartree*csec(3,i)
c      write (ifsec,6100)itq(i),elda,cse1,cse2,cse3
c      end do
c      write (ifsec,*)
c     .'-------------------------------------------------------------'
c
c formats
c 6100 format (1x,i4,5d13.5)
c 6400 format (' q-vector =',3f10.4)
c 6500 format (' spin =',i2,'   nbloch =',i4,'  deltaw =',f7.4,' a.u.' )
c 6600 format (1x,2i4,3f10.4)
c
c      return
c      end
c--------------------------------------------------------------------
      subroutine rsexx2 (nspin, itq, q, ntq,nq,ginv, symgg,ng, !crystal symmetry operation. july2012takao
     o                 vxco)
      implicit real*8 (a-h,o-z)
      dimension vxco(ntq,nq,nspin),q(3,nq),itq(ntq) !itq is not dependent on q, right?
      real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
      logical ::nocore,lfind
      real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3),qr(3),symgg(3,3,ng),sym(3,3)

c      if (nspin == 2) stop 'nspin=2 require some change!'
c
      inquire(file='NoCore',exist=nocore)
      if(nocore) then
        ifvxcfp  = iopen('VXCFPV',0,0,0)
        write(6,*)' OPEN VXCFPV = valence only Vxc'
      else
        write(6,*)' OPEN VXCFP '
        ifvxcfp  = iopen('VXCFP',0,0,0)
      endif
      read(ifvxcfp) ldim,nqbz
      write(6,*)' rsexx ldim,nqbz',ldim,nqbz
      allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
      do ikp = 1,nqbz
c        read(ifvxcfp) qqq(1:3,ikp),(vxcfpx(i,ikp),i=1,ldim)
        read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
cccccccccccccccccccccccccccccccccc
        write(6,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
cccccccccccccccccccccccccccccccccc
      enddo
      if(nocore) then
        isx = iclose('VXCFPV')
      else
        isx = iclose('VXCFP')
      endif
c
      do iq=1,nq
        do ikp=1,nqbz
        do ig = 1,ng
          sym = symgg(:,:,ig)
c          qr= qqq(1:3,ikp)
          qr=matmul(sym,qqq(1:3,ikp))

cccccccc
c          write(6,"('ikp qqq q=',i5,3d12.4,3x,3d12.4)") ikp,qqq(1:3,ikp),q(1:3,iq)
cccccccc

c         if(sum( (qqq(1:3,ikp)-q(1:3,iq))**2) <1d-10) then
c            ikpx=ikp
c            goto 100
c          endif
!above block to find out ikpx is replaced (even when translational equivalent points)
          lfind=.false.
          if(sum( (qr-q(1:3,iq))**2) <tolq) then
            lfind=.true.
          else
            call rangedq( matmul(ginv,q(1:3,iq)-qr), qx)
            if(sum(abs(qx))< tolq) lfind=.true.
          endif
          if(lfind) then
            ikpx=ikp
            goto 100
          endif
c
        enddo
        enddo
Cstop2rx 2013.08.09 kino        stop ' rsexx: not find ikp'
        call rx( ' rsexx: not find ikp')
 100    continue
        vxco(1:ntq,iq,1:nspin)=rydberg()*vxcfpx(itq(1:ntq),ikpx,1:nspin)
      enddo
      end


c------------------------------------------------------------------------
      subroutine rsexx (nspin, itq, q, ntq,nq,ginv, 
     o                 vxco)
      implicit real*8 (a-h,o-z)
      dimension vxco(ntq,nq,nspin),q(3,nq),itq(ntq) !itq is not dependent on q, right?
      real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
      logical ::nocore,lfind
      real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3)

c      if (nspin == 2) stop 'nspin=2 require some change!'
c
      inquire(file='NoCore',exist=nocore)
      if(nocore) then
        ifvxcfp  = iopen('VXCFPV',0,0,0)
        write(6,*)' OPEN VXCFPV = valence only Vxc'
      else
        write(6,*)' OPEN VXCFP '
        ifvxcfp  = iopen('VXCFP',0,0,0)
      endif
      read(ifvxcfp) ldim,nqbz
      write(6,*)' rsexx ldim,nqbz',ldim,nqbz
      allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
      do ikp = 1,nqbz
c        read(ifvxcfp) qqq(1:3,ikp),(vxcfpx(i,ikp),i=1,ldim)
        read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
cccccccccccccccccccccccccccccccccc
        write(6,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
cccccccccccccccccccccccccccccccccc
      enddo
      if(nocore) then
        isx = iclose('VXCFPV')
      else
        isx = iclose('VXCFP')
      endif
c
      do iq=1,nq
        do ikp=1,nqbz
cccccccc
c          write(6,"('ikp qqq q=',i5,3d12.4,3x,3d12.4)") ikp,qqq(1:3,ikp),q(1:3,iq)
cccccccc

c         if(sum( (qqq(1:3,ikp)-q(1:3,iq))**2) <1d-10) then
c            ikpx=ikp
c            goto 100
c          endif
!above block to find out ikpx is replaced (even when translational equivalent points)
          lfind=.false.
          if(sum( (qqq(1:3,ikp)-q(1:3,iq))**2) <tolq) then
            lfind=.true.
          else
            call rangedq( matmul(ginv,q(1:3,iq)-qqq(:,ikp)), qx)
            if(sum(abs(qx))< tolq) lfind=.true.
          endif
          if(lfind) then
            ikpx=ikp
            goto 100
          endif
c
        enddo
Cstop2rx 2013.08.09 kino        stop ' rsexx: not find ikp'
        call rx( ' rsexx: not find ikp')
 100    continue
        vxco(1:ntq,iq,1:nspin)=rydberg()*vxcfpx(itq(1:ntq),ikpx,1:nspin)
      enddo
      end

c--------------------------------------------------------------------
      double precision function egex (q,ef)

c 92.02.22
c electron gas bare exchange
c SEx(q) = -(kf/pi) { 1 + ((kf^2-q^2)/2q*kf) ln|q+kf/q-kf|
c ef     = kf^2

c q      = magnitude of q-vector
c ef     = fermi energy in Rydberg

c egex   = exchange energy in Hartree

      implicit real*8 (a-h,o-z)

      pi         = 4.d0*datan(1.d0)
      qf         = dsqrt(ef)
c
      if(q==0d0) then
        fac       =  1d0
      else
        fac       =  ((qf*qf-q*q)/(2.d0*q*qf) )
     .                             *dlog(dabs( (q+qf)/(q-qf) ))
      endif
      egex       = -(qf/pi)*(1d0 + fac)
      return
      end

      subroutine winfo(ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz,ef
     &    ,deltaw,alat,esmr)
      implicit none
      integer(4) :: ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz
      real(8) :: ef,deltaw,alat,esmr
      write (ifi,*)' ***'
      write (ifi,6700) nspin,nq,ntq
      write (ifi,6501) is,nbloch,ngp,ngc,nqbz,nqibz,ef
     &  ,deltaw,alat,ef,esmr
 6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4
     &        ,'  nqbz =',i6,'  nqibz =',
     .        i6,'   ef=', f10.4,' Rydberg'
     &        ,/,d23.16,' <= deltaw(Hartree)'
     &        ,/,d23.16,' <= alat'
     &        ,/,d23.16,' <= ef '
     &        ,/,d23.16,' <= esmr')
 6700 format (1x,3i4,'  nspin  nq  ntq')
      end

c------------------------------------------------------------------
      real(8) function lowesteigen(nspin,nband,qbz,nqbz)
      use m_readeigen, only :readeval
c lowest eband -----------------------------------------
      implicit none
      integer(4):: nband,nqbz,nspin,iq,is
      real(8)   :: ekt(nband,nqbz,nspin),qbz(3,nqbz)
      do is = 1,nspin
        do iq = 1,nqbz
          call readeval(qbz(:,iq),is, ekt(1,iq,is))
        enddo
      enddo
      lowesteigen = minval(ekt)
      end
c
c------------------------------------------------------------------
      real(8) function lowesteb(ifev,ispin, nband,nqbz)
c lowest eband -----------------------------------------
      implicit real*8 (a-h,o-z)
      integer(4)::   ifev(2)
      real(8):: ekt(nband,nqbz,ispin)
      do is = 1,ispin
        call rwdd   (ifev(is),   !read eigenvalues
     d              nband,nqbz, 
     o              ekt(1,1,is))
      enddo
      lowesteb = minval(ekt)
      end
c

c-----------------------------------------------------
      subroutine q0iwgt2(symops,ngrp,wqt,q0i,nq0i, 
     o                  wgt0)
c Get weight for each k-points near 0.
c wgt0(irreducible-k, irotation)
      implicit none
      integer(4) :: ixx,ix,i, ngrp,ig,nq0i
      real(8)     :: q0(3,6),q0i(3,6),symops(3,3,ngrp),sym(3,3),
     &  sym2(3,3),qt(3), q0in(3,ngrp*nq0i), wgt0(nq0i,ngrp), wqt(nq0i)
c      write(6,*)' q0iwgt:'
      wgt0 = 0d0
C      ixx = 0
      do i = 1,nq0i
c       ixx = 0
        qt = q0i(:,i)
c equivalence check
        do ig = 1,ngrp
          sym = symops(:,:,ig)

ctakao comment out once more to keep symmetry !
c          write(66,'(2i3,3f12.5)')i,ig,matmul(sym,qt)
cccccccccccccccccccccccccccccccccccccccc
c          if (ig > 1) then
c          do ix = 1,ig-1
cC            if(sum(abs(q0in(:,ix)-matmul(sym,qt)))<1d-10) then
c             sym2=symops(:,:,ix)
c             if(sum(abs(  matmul(sym2,qt)-
c     .                    matmul(sym,qt)  ))<1d-8) then
c               goto 1111
c             endif
c          enddo
c          endif
ccccccccccccccccccccccccccccccccccccccc
c          ixx = ixx+1
c          q0in (:,ixx) = matmul(sym,qt)
          wgt0 (i,ig)  = wqt(i) /dble(ngrp)
cccccccccccc
          write(6,'(2i3, f12.5)')i,ig, wgt0(i,ig)
cccccccccccc
c 1111     continue
        enddo !ig
c        do ig = 1, ngrp
c         wgt0 (i,ig) = wgt0 (i,ig)/dble(ixx)
c        enddo !ig
      enddo !i
ccccccccccccccccccccc
c     do ig=1,ngrp
c     write(6,*)' sym ig  =',ig
c     write(6,'(3f12.6)') symops(1,1:3,ig)
c     write(6,'(3f12.6)') symops(2,1:3,ig)
c     write(6,'(3f12.6)') symops(3,1:3,ig)
c      enddo
ccccccccccccccccccccc
      end




c This is taken from Ferdi's rw.f
c-----------------------------------------------------------------
      subroutine readxx(ifil)
      character*72 rchar
      integer(4):: n=1000
      do 10 i = 1,n
        read(ifil,'(a72)')rchar
c      write(6,*) rchar
c        j       = 0
c        call rmvbl (rchar,72,j)
c        rchar      = rchar(j+1:72)
        rchar=trim(adjustl(rchar))
        if(rchar(1:5) .eq. '*****') return
        if(rchar(1:5) .eq. '#####') return
   10 continue
Cstop2rx 2013.08.09 kino      stop 'readx: cannot find the string (gwsrc/wse.f/readxx)'
      call rx( 'readx: cannot find the string (gwsrc/wse.f/readxx)')
      end
c-------------------------------------------------------------------


      subroutine winfo2(ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz,ef
     &    ,ef2,deltaw,alat,esmr,esmr2)
      implicit none
      integer(4) :: ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz
      real(8) :: ef,ef2,deltaw,alat,esmr,esmr2
      write (ifi,*)' ***'
      write (ifi,6700) nspin,nq,ntq
      write (ifi,6501) is,nbloch,ngp,ngc,nqbz,nqibz,ef
     &  ,deltaw,alat,ef,ef2,esmr,esmr2
 6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4
     &        ,'  nqbz =',i6,'  nqibz =',
     .        i6,'   ef=', f10.4,' Rydberg'
     &        ,/,d23.16,' <= deltaw(Hartree)'
     &        ,/,d23.16,' <= alat'
     &        ,/,2d23.15,' <= ef '
     &        ,/,2d23.15,' <= esmr')
 6700 format (1x,3i4,'  nspin  nq  ntq')
      end


!! -------- S.F.Jan06 --------------------------------------------------
      subroutine q0iwgt3(allq0i,symops,ngrp,wqt,q0i,nq0i,           
     o                  wgt0)
!! Get weight for each k-points near 0.
!! wgt0(irreducible-k, irotation)
      implicit none
      logical allq0i !if true --> all rot. incl., even equivivalent
      !if false -> include only nonequivivalent rot.
      integer(4) :: ixx,ix,i, ngrp,ig,nq0i
      real(8)     :: q0(3,6),q0i(3,6),symops(3,3,ngrp),sym(3,3),
     &  sym2(3,3),qt(3), q0in(3,ngrp*nq0i), wgt0(nq0i,ngrp), wqt(nq0i)
      wgt0 = 1d0
      do i = 1,nq0i
        ixx = 0
        qt = q0i(:,i)
c equivalence check
        do ig = 1,ngrp
          sym = symops(:,:,ig)
          if (allq0i) then
            wgt0 (i,ig)  = wqt(i) /dble(ngrp)
            if(wgt0(i,ig)/=0d0) write(6,'(a,2i3,f12.5)')'allq0i=T ',i,ig, wgt0(i,ig)
          else
            do ix = 1,ig-1
              sym2=symops(:,:,ix)
              if(sum(abs(  matmul(sym2,qt)- matmul(sym,qt)  ))<1d-6) then
                wgt0(i,ig) = 0d0
                goto 1111
              endif
            enddo
            ixx = ixx+1
          endif
 1111     continue
        enddo !ig
        if (.not. allq0i) then
          do ig = 1, ngrp
            wgt0(i,ig) = wgt0(i,ig)*wqt(i)/dble(ixx)
            if(wgt0(i,ig)/=0d0) write(6,'(a,2i3, f12.5)')'allq0i=F ',i,ig, wgt0(i,ig)
          enddo !ig
        endif
      enddo !i
c      if (abs(sum(wgt0)-1d0) .gt. 1d-8) stop 'q0iwgt3: wrong wgt0'
c      if (abs(sum(wqt)-1d0) .gt. 1d-8)  stop 'q0iwgt3: wrong wqt'
      end


      integer*4 function ivsumxxx(ia,n)
      integer(4) ::ia(n)
      ivsumxxx=0
      do i=1,n
        if(ia(i)/=0) ivsumxxx=ivsumxxx+1
      enddo
      end
