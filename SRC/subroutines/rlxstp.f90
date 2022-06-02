subroutine rlxstp(natrlx,indrx_iv,xyzfrz,pdim)
  use m_lmfinit,only:  ctrl_nbas,ctrl_nitmv,ctrl_mdprm,ctrl_defm,ctrl_lfrce,ssite=>v_ssite
  !- Set up variables for relaxation
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   sctrl :struct for program flow parameters; see routine uctrl
  !i     Elts read: nbas nitmv mdprm defm ltb lfrce
  !i     Stored:
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: relax
  !i     Stored:
  !o Outputs:
  !o   indrx_iv(1,i) points to the ith relaxing component and
  !o   indrx_iv(2,i) points to the corresponding site
  !o   natrlx      # of relaxing degrees of freedom for atoms
  !o   xyzfrz(i)   T means all the ith components are frozen (T on input)
  !o   pdim:       dimension of the work array p, needed in relax
  !u Updates
  !u   21 Mar 06 mdprm(1)>100 signifies shear relaxation
  ! ----------------------------------------------------------------------
  implicit none
  logical :: xyzfrz(3)
  integer :: indrx_iv(2,*),nvar,natrlx,nitrlx,pdim
  double precision :: mdprm(6),defm(6)
  integer :: nbas,i,j,k,iprint,i1mach,ifrlx(3),igets,lrlx !,ltb
  logical :: force,mdxx
  nbas  = ctrl_nbas
  nitrlx= ctrl_nitmv
  mdprm =  ctrl_mdprm !call dcopy(size(ctrl_mdprm),ctrl_mdprm,1,mdprm,1)
  defm  = ctrl_defm   !call dcopy(size(ctrl_defm),ctrl_defm,1,defm,1)
  force = int(ctrl_lfrce) .gt. 0
  if ( .NOT. force .OR. nint(mdprm(1)) == 0) goto 9299
  mdxx = nint(mdprm(1)) .le. 3
  lrlx = mod(nint(mdprm(1)),100)
  ! --- Set relaxation variables ---
  j = 0
  if (mdxx) then
     xyzfrz = .false.
     goto 9299
  elseif (force .AND. mdprm(1) >= 100) then
     do  i = 1, 6
        if (defm(i) == 1) then
           j = j+1
           indrx_iv(1,j) = i
        endif
     enddo
  elseif (force) then
     do  i = 1, nbas
        call icopy(size(ssite(i)%relax),ssite(i)%relax,1,ifrlx,1)
        do  k = 1, 3
           if (ifrlx(k) == 1) then
              j = j + 1
              indrx_iv(1,j) = k
              indrx_iv(2,j) = i
              xyzfrz(k) = .false.
           endif
        enddo
     enddo
  endif
  natrlx = j
  if (natrlx == 0) goto 9299
  pdim = 0
  if ( .NOT. mdxx) then
     if (lrlx == 4) pdim = natrlx*7
     if (lrlx == 5) pdim = natrlx*(7+natrlx)
     if (lrlx == 6) pdim = natrlx*(12+2*natrlx)
  endif

  ! --- Printout ---
  if (iprint() >= 30) then
     if (lrlx == 4) then
        call info(0,1,0,' RLXSTP: Molecular statics (conjugate gradients) ..',0,0)
     elseif (lrlx == 5) then
        call info(0,1,0, ' RLXSTP: Molecular statics (Fletcher-Powell) ..',0,0)
     else
        call info(0,1,0, ' RLXSTP: Molecular statics (Broyden) ..',0,0)
     endif
     call info2(0,0,0, '         relaxing %i variables, %i iterations',natrlx,nitrlx)
     write(stdo,"('        x-tol=',d13.5,' g-tol=',d13.5,' step=',d13.5,' pdim=',d13.5)") &
          mdprm(3),mdprm(4),mdprm(5),pdim
  endif
9299 continue
end subroutine rlxstp
