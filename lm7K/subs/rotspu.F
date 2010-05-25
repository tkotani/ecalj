      subroutine rotspu(mode,nbas,nl,eula,neul,u)
C- Sets up rotation matrices for all sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci           1 if to reverse the sense of rotation
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1  if Euler angles are l-independent
Ci          nl if Euler angles are l-dependent
Ci          nl*nl if Euler angles are lm-dependent
Co Outputs
Co   u     :spinor rotation matrices
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nbas,neul,nl
      double complex u(2,2,nl*nl,nbas)
      double precision eula(nbas,neul,3)
C Local variables
      double precision xxc,xxs,alpha,gamma
      integer ib,il,im,ilm

C --- Spinor rotation matrices for all sites ---
      do  10  ib = 1, nbas
C ... Assume initially Euler angles are not l-dependent
        xxc = dcos(eula(ib,1,2)/2)
        xxs = dsin(eula(ib,1,2)/2)
        alpha = eula(ib,1,1)
        gamma = eula(ib,1,3)
        ilm = 0
        do  10  il = 1, nl
C   ... If euler angles are l dependent
          if (neul .eq. nl) then
            xxc = dcos(eula(ib,il,2)/2)
            xxs = dsin(eula(ib,il,2)/2)
            alpha = eula(ib,il,1)
            gamma = eula(ib,il,3)
          endif
          do  10  im = -il+1, il-1
            ilm = ilm+1
C   ... If euler angles are lm dependent
            if (neul .eq. nl*nl) then
              xxc = dcos(eula(ib,ilm,2)/2)
              xxs = dsin(eula(ib,ilm,2)/2)
              alpha = eula(ib,ilm,1)
              gamma = eula(ib,ilm,3)
            endif
            if (mode .eq. 1) then
              u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
              u(1,2,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
              u(2,1,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
              u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
            else
              u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
              u(2,1,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
              u(1,2,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
              u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
            endif
C        print *, ilm,ib
C        if (ilm .eq. 1 .and. iprint() .ge. 100)
C     .    call zprm('u in rotspu',2,u(1,1,ilm,ib),2,2,2)
   10 continue
      end

