c$$$      subroutine rotspu(mode,nbas,nl,eula,neul,u)
c$$$C- Sets up rotation matrices for all sites
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   mode  :1s digit
c$$$Ci           1 if to reverse the sense of rotation
c$$$Ci   nbas  :size of basis
c$$$Ci   nl    :(global maximum l) + 1
c$$$Ci   eula  :Euler angles for noncollinear spins
c$$$Ci   neul  :1  if Euler angles are l-independent
c$$$Ci          nl if Euler angles are l-dependent
c$$$Ci          nl*nl if Euler angles are lm-dependent
c$$$Co Outputs
c$$$Co   u     :spinor rotation matrices
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$      integer mode,nbas,neul,nl
c$$$      double complex u(2,2,nl*nl,nbas)
c$$$      double precision eula(nbas,neul,3)
c$$$C Local variables
c$$$      double precision xxc,xxs,alpha,gamma
c$$$      integer ib,il,im,ilm
c$$$
c$$$C --- Spinor rotation matrices for all sites ---
c$$$      do  102  ib = 1, nbas
c$$$C ... Assume initially Euler angles are not l-dependent
c$$$        xxc = dcos(eula(ib,1,2)/2)
c$$$        xxs = dsin(eula(ib,1,2)/2)
c$$$        alpha = eula(ib,1,1)
c$$$        gamma = eula(ib,1,3)
c$$$        ilm = 0
c$$$      do  101  il = 1, nl
c$$$C   ... If euler angles are l dependent
c$$$        if (neul .eq. nl) then
c$$$          xxc = dcos(eula(ib,il,2)/2)
c$$$          xxs = dsin(eula(ib,il,2)/2)
c$$$          alpha = eula(ib,il,1)
c$$$          gamma = eula(ib,il,3)
c$$$        endif
c$$$      do  10  im = -il+1, il-1
c$$$        ilm = ilm+1
c$$$C   ... If euler angles are lm dependent
c$$$        if (neul .eq. nl*nl) then
c$$$          xxc = dcos(eula(ib,ilm,2)/2)
c$$$          xxs = dsin(eula(ib,ilm,2)/2)
c$$$          alpha = eula(ib,ilm,1)
c$$$          gamma = eula(ib,ilm,3)
c$$$        endif
c$$$        if (mode .eq. 1) then
c$$$          u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
c$$$          u(1,2,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
c$$$          u(2,1,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
c$$$          u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
c$$$        else
c$$$          u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
c$$$          u(2,1,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
c$$$          u(1,2,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
c$$$          u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
c$$$        endif
c$$$C        print *, ilm,ib
c$$$C        if (ilm .eq. 1 .and. iprint() .ge. 100)
c$$$C     .    call zprm('u in rotspu',2,u(1,1,ilm,ib),2,2,2)
c$$$   10 continue
c$$$ 101  continue
c$$$ 102  continue
c$$$      end
c$$$
