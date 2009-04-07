      subroutine suq0x(alat,plat,alp,qbz,nnn,nx0,xn,q0x,wt0)
C- Generate q0x - special qp to substitute for Gamma point
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   alp   :parameter in auxilliary function
Ci   qbz   :qp (full BZ)
Ci   nnn   :number of qp in the BZ
Ci   nx0   :1 or 2 for number of special qp q0x
Ci   xn    :ratio of weights is xn:1-xn for nx0=2 case
Co Outputs
Co   q0x   :6 * nx0 special qp to substitute for q(0)
Co   wt0   :weighting for each of the q0x
Cr Remarks
Cr   Linear response has 1/q**2 singularity at gamma
Cr   q0x is specially constructed from auxilliary function
Cr   to get around it.
Cr
Cr   wtt=1d0/(n1q*n2q*n3q) is the weight for each q-point.
Cr   Determine q0x so that wtrue  = wsumau + wtt*auxfun(q0x) .
Cr
Cr  nx0 = number of q points near q=0.  nx0=1 or 2 now.
Cr  TK tested nx0=2 case for si222, but only slightly different
Cr  (0.01eV) dSE
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nnn,nx0
      double precision alat,plat(3,3),alp,qbz(3,nnn),xn
      double precision q0x(3,6*nx0),wt0(6*nx0)
C ... Local parameters
      integer,allocatable :: ngc(:),kv(:,:),ngvect(:,:,:)
      integer iq,n1,n2,n3,ngmx,ngcmx,i,k,stdo,lgunit,iprint
      double precision q(3),QpGcut,xx,qlat(3,3),vol,q2oq1,volinv
      double precision funa(nnn),pi,wtt,wsumau,wtrue,auxf0,qini,qc,
     .  snorm,qa,qb,fa,fb,fder,auxfun,auxfun6,auxfun6xn
C     double precision fc,qa0

C     call tcn('suq0x')
      call dinv33(plat,1,qlat,vol)
      vol = alat**3*abs(vol)
      pi = 3.1415926535897932d0
      stdo = lgunit(1)

C ... Get the q-dependent G vectors, and their number
C     Energy cutoff, in a.u. : alp * QpGcut**2 = 25.  (exp(- alp * QpGcut**2))
      QpGcut = sqrt(25d0/alp)
      allocate(ngc(nnn))
      call pshpr(0)
      do  iq  = 1, nnn
        q  = qbz(1:3,iq)
        n1 = 0
        n2 = 0
        n3 = 0
        call gvlst2(alat,plat,q,n1,n2,n3,0d0,QpGcut,0,100,0,ngc(iq),xx,
     .    xx,xx,xx)
      enddo
      ngcmx = maxval(ngc)
      n1 = n1*2
      n2 = n3*2
      n3 = n3*2
      allocate(ngvect(3,ngcmx,nnn),kv(3,n1*n2*n3*8))
      do  iq  = 1, nnn
        q = qbz(1:3,iq)
        ngmx = ngc(iq)
        call gvlst2(alat,plat,q,n1,n2,n3,0d0,QpGcut,0,2,ngmx,ngc(iq),kv,
     .    xx,xx,ngvect(1,1,iq))
      enddo
      call poppr
      deallocate(kv)

C ... Weight each of 6 points equally
      if (nx0 .eq. 1) then
       wt0 = 1d0/6d0
      elseif (nx0 .eq. 2) then
        wt0(1:6)  = xn/6d0
        wt0(7:12) = (1d0-xn)/6d0
        q2oq1 = sqrt(xn/(xn-1d0))
      else
        call rxi(' suq0x: unsupported input, nx0 = ',nx0)
      endif

C ... Make initial estimate qini
      volinv  = (2*pi)**3/vol
      funa(1) = 0d0
      do  iq = 2, nnn
      funa(iq) = auxfun(qbz(1,iq),alp,alat,qlat,ngc(iq),ngvect(1,1,iq))
      enddo
      wtt = 1d0/dble(nnn)
      wsumau = sum(funa)*wtt
      wtrue  = 4*pi/volinv*sqrt(pi/alp)/2d0
      auxf0  = (wtrue - wsumau)/wtt
      qini   = sqrt(1/auxf0)*alat/(2*pi)
      if (iprint() .ge. 80)
     .  call awrit3(' suq0x : auxf(2..nq) = %,6d  auxf(true) = %,6d'//
     .  '  starting q = %,6d',' ',80,stdo,wsumau,wtrue,qini)

C --- Find q such that auxf0 = auxfun(q) ---
      snorm = sqrt( sum ((qlat(:,1)**2)) )
      q0x(:,1)=  qlat(:,1)/snorm
      q0x(:,2)= -qlat(:,1)/snorm
      snorm = sqrt( sum ((qlat(:,2)**2)) )
      q0x(:,3)=  qlat(:,2)/snorm
      q0x(:,4)= -qlat(:,2)/snorm
      snorm = sqrt( sum ((qlat(:,3)**2)) )
      q0x(:,5)=  qlat(:,3)/snorm
      q0x(:,6)= -qlat(:,3)/snorm
      if (nx0 .eq. 2) q0x(:,7:12) = q0x(:,1:6)

C --- find qa for auxf0 = aufun(qa), Newtonian minimization ---
      iq = 1
      qa = qini
      qb = qini + 0.001d0
      do  i = 1, 100
        k = i
        if (nx0 .eq. 1) then
          fa = auxfun6(qa,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq))
          fb = auxfun6(qb,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq))
        else
          fa = auxfun6xn(qa,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq),xn)
          fb = auxfun6xn(qb,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq),xn)
        endif
C        write(stdo,*) 'qa fa-auxf0 =',qa,fa-auxf0
C        write(stdo,*) 'qb fb-auxf0 =',qb,fb-auxf0
        fder = (fb-fa)/(qb-qa)
        qc = qa + (auxf0-fa)/fder
        qb = qa
        qa = qc
        if (abs(qa-qb) .lt. 1d-10) exit
      enddo

      if (k .eq. 100) call rx('search for q failed to converge')

C --- refine using bisection (probably not needed) ---
C      qa0 = qa
C      qb = qa0 + 1d-9
C      qa = qa0 - 1d-9
C      do  i = 1, 100
C        k = i
C        if (nx0 .eq. 1) then
C          fa = auxfun6(qa,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq))
C          fb = auxfun6(qb,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq))
C          qc = 0.5d0*(qa+qb)
C          fc = auxfun6(qc,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq))
C        else
C          fa = auxfun6xn(qa,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq),xn)
C          fb = auxfun6xn(qb,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq),xn)
C          qc = 0.5d0*(qa+qb)
C          fc = auxfun6xn(qc,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq),xn)
C        endif
C
C        if ( (auxf0-fa)*(auxf0-fb) .gt. 0d0 ) then
C          write(stdo,*) i,fa,fb,auxf0
C          write(stdo,*) qa,qb
C          call rx('suq0x: something wrong (auxf0-fa)*(auxf0-fb) >0d0')
C        elseif ((auxf0-fa)*(auxf0-fc) .lt. 0d0 ) then
C          qb = qc
C        elseif ((auxf0-fb)*(auxf0-fc) .lt. 0d0 ) then
C          qa = qc
C        endif
C        if (abs(fa-fb) .lt. 1d-12) exit
C      enddo
C      if (k .eq. 100) call rx('search for q failed to converge')

C --- q0x from qa; printout ---
      if (nx0 .eq. 1) then
        fa = wsumau +
     .       auxfun6(qa,q0x,alp,alat,qlat,ngc(iq),ngvect(1,1,iq))*wtt
        q0x(:,1:6) =  qa * q0x(:,1:6)
        if (iprint() .ge. 20) then
          call awrit3(' suq0x : %i special gamma-points  auxf = %,6d'//
     .      ' (error = %;2,2g)',' ',80,stdo,6*nx0,fa,wtrue-fa)
          call awrit2('%9fqa = %,6d  q0x(1) =%3:1,8;8d',' ',80,stdo,qa,
     .      q0x)
        endif
      else
        q2oq1 = sqrt(xn/(xn-1d0))
        qb    = q2oq1*qa
        fa = wsumau + auxfun6xn(qa,q0x, alp,alat,qlat,ngc(iq),
     .                          ngvect(1,1,iq),xn)*wtt
        q0x(:,1:6)  = qa * q0x(:,1:6)
        q0x(:,7:12) = qb * q0x(:,7:12)
        if (iprint() .ge. 20) then
          call awrit4(' suq0x: qa = %,6d  qb = %,6d  auxf = %,6d'
     .      //' (error = %;2,2g)',' ',80,stdo,qa,qb,fa,wtrue-fa)
          call awrit1('        q0x(1) =%3:1;12F',' ',80,stdo,q0x)
          call awrit1('        q0x(7) =%3:1;12F',' ',80,stdo,q0x(1,7))
        endif
      endif
   99 continue
      deallocate(ngc)
      deallocate(ngvect)
C     call tcx('suq0x')
      end
      double precision function auxfun(q,alp,alat,qlat,ngc,ngvect)
C- Auxilliary function for 1 qp
C     implicit none
      integer ig,ngc,ngvect(3,ngc)
      double precision alat,q(3),tpiba,qg(3),qlat(3,3),qg2,pi,alp
      parameter (pi=3.1415926535897932d0)
      tpiba  = 2*pi/alat
      auxfun = 0d0
      do  ig = 1, ngc
        qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvect(1:3,ig)))
        qg2     = sum(qg(1:3)**2)
        auxfun  = auxfun + exp(-alp*qg2)/qg2
      enddo
      end

      double precision function auxfun6(qa,q0x,alp,alat,qlat,ngc,ngvect)
C- Auxilliary function averaged over 6 qp
C     implicit none
      integer i,ngc,ngvect(3,ngc)
      double precision qa,alat,qlat(3,3),alp,auxfun,q0x(3,6)

      auxfun6 = 0d0
      do  i = 1, 6
        auxfun6 = auxfun6 +
     .            auxfun( (/qa*q0x(1,i), qa*q0x(2,i), qa*q0x(3,i)/),
     .            alp,alat,qlat,ngc,ngvect)
      enddo
      auxfun6 = auxfun6/6d0
      end

      double precision function auxfun6xn(qa,q0x,alp,alat,qlat,ngc,
     .  ngvect,xn)
C- Auxilliary function for nx0=2 case averaged over 6 qp
C     implicit none
      integer i,ngc,ngvect(3,ngc)
      double precision xn,wtq1,wtq2,q2oq1,qb,alp,auxfun,q0x(3,6)
      double precision qa,alat,qlat(3,3)
      wtq1  = xn/6d0
      wtq2 = (1d0-xn)/6d0
      q2oq1 = sqrt(xn/(xn-1d0))
      qb = q2oq1*qa
      auxfun6xn = 0d0
      do  i = 1, 6
        auxfun6xn = auxfun6xn
     .    + wtq1*auxfun( (/qa*q0x(1,i),qa*q0x(2,i),qa*q0x(3,i)/),
     .    alp,alat,qlat,ngc,ngvect)
     .    + wtq2*auxfun( (/qb*q0x(1,i),qb*q0x(2,i),qb*q0x(3,i)/),
     .    alp,alat,qlat,ngc,ngvect)
      enddo
      end

      subroutine q0irre(q0,wt0,nx06,g,ng, q0i,nq0i,wt)
C- Find inequivalent q0 points
C ----------------------------------------------------------------------
Ci Inputs
Ci   q0    :6 * nx0 special qp to substitute for q(0)
Ci   wt0   :weighting for each of the q0
Ci   nx06  :6 * nx0
Ci   g     :point group operations
Ci   ng    :number of point group operations
Co Outputs
Co   q0i   :irreducible points of the q0
Co   nq0i  :number of q0i
Co   wt    :weightings for points
Cu Updates
Cu   16 Feb 01 routine copied from Kotani's
C ----------------------------------------------------------------------
C     implicit none
      integer ng,nx06,nq0i
      real(8) :: q0(3,nx06),q0i(3,nx06),g(3,3,ng),wt0(nx06),wt(nx06)
C ... Local parameters
      integer ixx,i,ix,ig,stdo,iprint,lgunit
      double precision sym(3,3),qt(3)

      stdo = lgunit(1)
      wt = 1d99
      ixx = 0
      do  i = 1, nx06
        qt = q0(:,i)
        do ix = 1,ixx
        do ig = 1,ng
          sym = g(:,:,ig)
          if(sum(abs(q0i(:,ix)-matmul(sym,qt)))<1d-6) then
            wt(ix) = wt(ix)+wt0(i)
            goto 990
          endif
        enddo
        enddo

        ixx = ixx+1
        q0i(:,ixx) = qt
        wt   (ixx) = wt0(i)
  990   continue
      enddo
      nq0i = ixx

      if (iprint() .lt. 20) return
      call awrit1(' q0irre: %i irreducible special gamma-point(s)',
     .  ' ',120,stdo,nq0i)
      if (iprint() .le. 40) return
      call awrit1('  iq%20fqp%21fweight',' ',120,stdo,nq0i)
      write(stdo,"(i4,f12.7,2x,3f12.7)")(i,q0i(1:3,i),wt(i),i=1,nq0i)

      end

