      subroutine rlxstp(natrlx,indrx_iv,xyzfrz,pdim)
      use m_lmfinit,only:  ctrl_nbas,ctrl_nitmv,ctrl_mdprm,ctrl_defm,ctrl_lfrce,ssite=>v_ssite
C- Set up variables for relaxation
C ----------------------------------------------------------------------
Ci Inputs:
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nitmv mdprm defm ltb lfrce
Ci     Stored:
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: relax
Ci     Stored:
Co Outputs:
Co   indrx_iv(1,i) points to the ith relaxing component and
Co   indrx_iv(2,i) points to the corresponding site
Co   natrlx      # of relaxing degrees of freedom for atoms
Co   xyzfrz(i)   T means all the ith components are frozen (T on input)
Co   pdim:       dimension of the work array p, needed in relax
Cu Updates
Cu   21 Mar 06 mdprm(1)>100 signifies shear relaxation
C ----------------------------------------------------------------------
      implicit none
      logical xyzfrz(3)
      integer indrx_iv(2,*),nvar,natrlx,nitrlx,pdim
      double precision mdprm(6),defm(6)
      integer nbas,i,j,k,iprint,i1mach,ifrlx(3),igets,lrlx !,ltb
      logical force,mdxx
      nbas  = ctrl_nbas
      nitrlx= ctrl_nitmv
      mdprm =  ctrl_mdprm !call dcopy(size(ctrl_mdprm),ctrl_mdprm,1,mdprm,1)
      defm  = ctrl_defm   !call dcopy(size(ctrl_defm),ctrl_defm,1,defm,1)
      force = int(ctrl_lfrce) .gt. 0
      if (.not. force .or. nint(mdprm(1)) .eq. 0) goto 9299
      mdxx = nint(mdprm(1)) .le. 3
      lrlx = mod(nint(mdprm(1)),100)
C --- Set relaxation variables ---
      j = 0
      if (mdxx) then
        xyzfrz = .false.
        goto 9299
      elseif (force .and. mdprm(1) .ge. 100) then
        do  i = 1, 6
          if (defm(i) .eq. 1) then
            j = j+1
            indrx_iv(1,j) = i
          endif
        enddo
      elseif (force) then
        do  i = 1, nbas
          call icopy(size(ssite(i)%relax),ssite(i)%relax,1,ifrlx,1)
          do  k = 1, 3
            if (ifrlx(k) .eq. 1) then
              j = j + 1
              indrx_iv(1,j) = k
              indrx_iv(2,j) = i
              xyzfrz(k) = .false.
            endif
          enddo
        enddo
      endif
      natrlx = j
      if (natrlx .eq. 0) goto 9299
      pdim = 0
      if (.not. mdxx) then
        if (lrlx .eq. 4) pdim = natrlx*7
        if (lrlx .eq. 5) pdim = natrlx*(7+natrlx)
        if (lrlx .eq. 6) pdim = natrlx*(12+2*natrlx)
      endif

C --- Printout ---
      if (iprint() .ge. 30) then
        if (lrlx .eq. 4) then
          call info(0,1,0,' RLXSTP: Molecular statics (conjugate gradients) ..',0,0)
        elseif (lrlx .eq. 5) then
          call info(0,1,0, ' RLXSTP: Molecular statics (Fletcher-Powell) ..',0,0)
        else
          call info(0,1,0, ' RLXSTP: Molecular statics (Broyden) ..',0,0)
        endif
        call info2(0,0,0, '         relaxing %i variables, %i iterations',natrlx,nitrlx)
        call awrit4('         x-tol=%d, g-tol=%d, step=%d (pdim=%i)',
     .  ' ',120,i1mach(2),mdprm(3),mdprm(4),mdprm(5),pdim)
      endif
 9299 continue
      end subroutine rlxstp
