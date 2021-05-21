      subroutine rlxstp(ssite,natrlx,nvar,indrlx,xyzfrz,pdim)!sctrl,
      use m_struc_def           !Cgetarg
      use m_lmfinit,only:  ctrl_nbas,ctrl_nitmv,ctrl_mdprm,ctrl_defm,ctrl_lfrce
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
Co   indrlx(1,i) points to the ith relaxing component and
Co   indrlx(2,i) points to the corresponding site
Co   natrlx      # of relaxing degrees of freedom for atoms
Co   nvar:       # of variables to relax (natrlx)
Co   xyzfrz(i)   T means all the ith components are frozen (T on input)
Co   pdim:       dimension of the work array p, needed in relax
Cr Remarks:
Cr   At version 6.15 first attempt to restore volume and shear relaxations
Cr   At version 6.10 volume and shear relaxations have been removed;
Cr   hence nvar=natrlx.
Cu Updates
Cu   21 Mar 06 mdprm(1)>100 signifies shear relaxation
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      logical xyzfrz(3)
      integer indrlx(2,*),nvar,natrlx,nitrlx,pdim
c      type(s_ctrl)::sctrl
      type(s_site)::ssite(*)

C Local variables
      double precision mdprm(6),defm(6)
C     double precision autime,autmp
      integer nbas,i,j,k,iprint,i1mach,ifrlx(3),igets,lrlx !,ltb
      logical force,md,bittst
C     character*1 comp(3)
C     data comp /'x','y','z'/
C     data autime/0.048377d0/
C     data autmp/0.6333328d-5/


      nbas=ctrl_nbas
      nitrlx=ctrl_nitmv
c       i_copy_size=size(ctrl_mdprm)
c       call dcopy(i_copy_size,ctrl_mdprm,1,mdprm,1)
c       i_copy_size=size(ctrl_defm)
c       call dcopy(i_copy_size,ctrl_defm,1,defm,1)
      call dcopy(size(ctrl_mdprm),ctrl_mdprm,1,mdprm,1)
      call dcopy(size(ctrl_defm),ctrl_defm,1,defm,1)
c       ltb=ctrl_ltb

      nvar = 0
c       force = bittst ( ltb , 16 ) .or. int(ctrl_lfrce) .gt. 0
      force = int(ctrl_lfrce) .gt. 0

      if (.not. force .or. nint(mdprm(1)) .eq. 0) return
      md = nint(mdprm(1)) .le. 3
      lrlx = mod(nint(mdprm(1)),100)

C --- Set relaxation variables ---
      j = 0
      if (md) then
        xyzfrz(1) = .false.
        xyzfrz(2) = .false.
        xyzfrz(3) = .false.
        return
      elseif (force .and. mdprm(1) .ge. 100) then
        do  i = 1, 6
          if (defm(i) .eq. 1) then
            j = j+1
            indrlx(1,j) = i
          endif
        enddo
      elseif (force) then
        do  i = 1, nbas
c           i_copy_size=size(ssite(i)%relax)
c           call i8icopy(i_copy_size,ssite(i)%relax,1,ifrlx,1)
          call icopy(size(ssite(i)%relax),ssite(i)%relax,1,ifrlx,1)

          do  k = 1, 3
            if (ifrlx(k) .eq. 1) then
              j = j + 1
              indrlx(1,j) = k
              indrlx(2,j) = i
              xyzfrz(k) = .false.
            endif
          enddo
        enddo
      endif
      natrlx = j
      nvar = natrlx
      if (nvar .eq. 0) return
      pdim = 0
      if (.not. md) then
        if (lrlx .eq. 4) pdim = nvar*7
        if (lrlx .eq. 5) pdim = nvar*(7+nvar)
        if (lrlx .eq. 6) pdim = nvar*(12+2*nvar)
      endif

C --- Printout ---
      if (iprint() .ge. 30) then
        if (lrlx .eq. 4) then
          call info(0,1,0,
     .    ' RLXSTP: Molecular statics (conjugate gradients) ..',0,0)
        elseif (lrlx .eq. 5) then
          call info(0,1,0,
     .    ' RLXSTP: Molecular statics (Fletcher-Powell) ..',0,0)
        else
          call info(0,1,0,
     .    ' RLXSTP: Molecular statics (Broyden) ..',0,0)
        endif
        call info2(0,0,0,
     .  '         relaxing %i variables, %i iterations',nvar,nitrlx)
        call awrit4('         x-tol=%d, g-tol=%d, step=%d (pdim=%i)',
     .  ' ',120,i1mach(2),mdprm(3),mdprm(4),mdprm(5),pdim)
      endif

      end subroutine rlxstp


