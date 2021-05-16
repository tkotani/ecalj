      subroutine relax(ssite,sspec,it,indrlx,natrlx,nvar,f, !,sctrl prgnam,
     .p,w,nelts,delta,bas,icom)

      use m_struc_def  !Cgetarg
      use m_lmfinit,only:ctrl_nbas,ctrl_nitmv,ctrl_mdprm

C- Relax atomic positions and volume using variable metric algorithm
C ----------------------------------------------------------------------
Ci Inputs:
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nitmv mdprm ltb
Ci     Stored:
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec relax
Ci     Duplicate: spec relax relax
Ci     Stored:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:
Ci     Stored:    name
Ci   it:        iteration number
Ci   indrlx(1,i) points to the ith relaxing component and
Ci   indrlx(2,i) points to the corresponding site (see rlxstp)
Ci   natrlx:    # of relaxing degrees of freedom for atoms or shear
Ci   nvar:      # of relaxing degrees of freedom for atoms or shear
Ci   f:         forces
Ci   p:         work array for gradzr (dimensioned in rlxstp)
Ci   w:         work array for the hessian
Ci   delta,nelts used for printout only, if U[L] turned on (TBE)
Cio Inputs/Outputs
Cio  bas        atom positions or shear coordinates; se Remarks
Cio             New estimate on output
Co Outputs
Co   icom       1 if relaxation is complete
Co             -1 if relaxation is not complete and current gradient
Co                is larger than prior case.
Co                p(*,nm) holds positions for prior case
Co              0 otherwise
Cl Local variables
Cl   wkg:       work array used by gradzr.
Cl              It should not be modified between iterations
Cr Remarks:
Cr   At version 6.10 volume and shear relaxations have been removed;
Cr   hence nvar=natrlx.
Cr   When ctrl->mdprm(1) > 100, bas are shear coords
Cu Updates
Cu   21 Mar 06 mdprm(1)>100 signifies shear relaxation
Cu   08 Mar 06 Some improvements to gradzr; repackaged MPI
Cu   15 Feb 02 (ATP) Added MPI parallelization
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
c      character prgnam*8
      integer it,nit,natrlx,nvar,nelts,icom,procid,master,mpipid
      integer indrlx(2,natrlx)
      real(8):: f(3,*) , w(nvar,nvar) , p(nvar,6) , delta(nelts,*) 
     ., bas(3,*)
c      type(s_ctrl)::sctrl
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C Local Variables
      integer i,j,ipr,ifi,fopna,ix,lgunit,nbas,ltb,ifrlx(3),
     .ir,iprint,isw,rdm,lrlx,is,idamax,stdo,stdl,nd,ns,nkill
      parameter (nd=4,ns=6)
      logical rdhess,bittst,lpp,cmdopt,a2bin,lshr !U,
      double precision mdprm(6),step,xtol,gtol,xtoll,grfac,wkg(28),ddot,
     .xv(10)
      equivalence (step,mdprm(5)), (xtol,mdprm(3)), (gtol,mdprm(4))
      character clablj*8,dumstr*6,strn*128
      save ir,wkg
      data wkg /28*0d0/

      stdo = lgunit(1)
      stdl = lgunit(2)
      master = 0
      procid = mpipid(1)
      if (procid .eq. master) call pshpr(iprint()+30)
      call tcn('relax')

C --- Setup ---
      nbas=  ctrl_nbas
      nit =  ctrl_nitmv
      mdprm= ctrl_mdprm
      nkill  = nint(mdprm(6))
      lrlx   = nint(mdprm(1))
      lshr   = lrlx .gt. 100
      lrlx   = mod(lrlx,100)
      rdhess = nint(mdprm(2)) .eq. 1
      call getpr(ipr)
C --- Make vector of positions and gradients ---
      j = 1
      if (.not. lshr) then
        do  10  i = 1, natrlx
          p(j,1) = bas(indrlx(1,i),indrlx(2,i))
          p(j,2) = -f(indrlx(1,i),indrlx(2,i))
          j = j + 1
   10   continue
      else
        call dcopy(natrlx,bas,1,p,1)
        call dcopy(natrlx,f,1,p(1,2),1)
      endif
C --- Initialization ---
      if (it .eq. 1) then
        ir = 0
        call dcopy(nvar*nvar,0d0,0,w,1)
        call dcopy(nvar,1d0,0,w,nvar+1)
C   ... Read Hessian from disc
        if (rdhess) then
          if (procid .eq. master) then
            ifi = fopna('HSSN',-1,4)
            rewind ifi
C     ... Hessian written in rdm-compatible format
            if (rdm(ifi,2,nvar*nvar,' ',w,nvar,nvar) .lt. 0) then
C       ... Attempt to read Hessian failed.  Make it unity
              call info0(-20,0,0,
     .        ' RELAX: (warning) failed to read Hessian; set to unity')
              call dcopy(nvar*nvar,0d0,0,w,1)
              call dcopy(nvar,1d0,0,w,nvar+1)
            else
              call info0(-20,0,0,' RELAX: Hessian read from disc')
            endif
            call fclose(ifi)
          endif
          call mpibc1(w,nvar*nvar,4,.false.,'relax','hssn')
        else
          call info0(-30,0,0,' RELAX: no Hessian read from disc')
        endif
      endif

C --- Relax ---
      grfac = 1.2d0
      if (cmdopt('-grfac=',6,0,strn)) then
        j = 6
        if (.not. a2bin(strn,grfac,4,0,' ',j,len(strn))) then
          print *, 'RELAX: Ignored command line value of grfac'
          grfac = 1.2d0
        endif
      endif
      xtoll = abs(xtol)
      if (cmdopt('-xtoll=',6,0,strn)) then
        j = 6
        if (.not. a2bin(strn,xtoll,4,0,' ',j,len(strn))) then
          print *, 'RELAX: Ignored command line value of xtoll'
          xtoll = abs(xtol)
        endif
      endif
C ... the user sets xtol _and_ gtol.
      if (lrlx .eq. 4) isw = 00021 + 40
      if (lrlx .eq. 5) isw = 00121 + 40
      if (lrlx .eq. 6) isw = 00221 + 00
      if (lrlx.eq.4 .or. lrlx.eq.5 .or. lrlx.eq.6) then
        if (xtol .eq. 0 .and. gtol .eq. 0)
     .  call rx('RELAX: both xtol and gtol are zero')
        if (gtol .eq. 0) isw = isw-10
        if (xtol .eq. 0) isw = isw-20
      endif
      call gradzr(nvar,p,w,xtoll,step,xtol,gtol,grfac,wkg,isw,ir)
      if (ir .gt. 0) then
        call rx1('RELAX: gradzr returned ir=%i ... aborting',ir)
      endif
      if ((lrlx.eq.4 .or. lrlx.eq.5) .and. ipr .ge. 20) then
        if (procid .eq. master) then
          call awrit4(' fp rlx ln %d:  %1j'//
     .    '%-1j%?#n==-0#cnvg##'//
     .    '%-1j%?#n==-1#new ln##'//
     .    '%-1j%?#n==-2#brack ln##'//
     .    '%-1j%?#n==-3#brack ln##'//
     .    '%-1j%?#n==-4#extrap ln ##'//
     .    '%-1j%?#n==-5#extrap ln ##'//
     .    '%-1j%?#n==-6#trouble!!##'//
     .    '  dxmx=%;7F  |g|=%1;3g',' ',80,stdl,wkg(19),ir,
     .    wkg(1)*p(idamax(nvar,p(1,nd),1),nd),
     .    dsqrt(ddot(nvar,p(1,2),1,p(1,2),1)))
          if (ipr .ge. 20) call awrit4(' RELAX line %d:  %1j'//
     .    '%-1j%?#n==-0#converged to tolerance##'//
     .    '%-1j%?#n==-1#new line minimization##'//
     .    '%-1j%?#n==-2#bracketed root this line##'//
     .    '%-1j%?#n==-3#bracketed root this line##'//
     .    '%-1j%?#n==-4#extrapolated along this line##'//
     .    '%-1j%?#n==-5#extrapolated along this line##'//
     .    '%-1j%?#n==-6#is in trouble##'//
     .    ';  max shift=%;7F  |g|=%1;3g',' ',80,stdo,wkg(19),ir,
     .    wkg(1)*p(idamax(nvar,p(1,nd),1),nd),
     .    dsqrt(ddot(nvar,p(1,2),1,p(1,2),1)))
        endif
      elseif (lrlx .eq. 6) then
        if (procid .eq. master) then
          call awrit3(' fp rlx Br '//
     .    '%?#n==0#cnvg#%-1jit %i#  dxmx=%;7F  |g|=%1;3g',' ',
     .    80,stdl,-ir,p(idamax(nvar,p(1,ns),1),ns),
     .    dsqrt(ddot(nvar,p(1,2),1,p(1,2),1)))
          if (ipr .ge. 20) call awrit4(' RELAX:  %1j'//
     .    '%-1j%?#n==-0#converged to tolerance##'//
     .    '%-1j%?#n<0#completed Broyden iter %i#%j#'//
     .    ';  max shift=%;7F  |g|=%1;3g',' ',80,stdo,
     .    ir,-ir,p(idamax(nvar,p(1,ns),1),ns),
     .    dsqrt(ddot(nvar,p(1,2),1,p(1,2),1)))
        endif
      endif
      if (ir .eq. 0) icom = 1
C --- Printout  ---
      if (ipr .ge. 40) then
        call info0(-40,0,0,'        Gradients:')
        print 100, (p(i,2), i = 1, nvar)
        if (lrlx .ne. 4 .and. (it .ne. 1 .or. xtol .ne. 0)) then
          call info0(-40,0,0,'      Diagonal inverse Hessian:')
          print 100, (w(i,i), i = 1,nvar)
        endif
      endif
  100 format(10f8.3)
C --- Update atom positions ---
      call prelx1(1,1,lshr,natrlx,nvar,indrlx,p,bas)
C --- Write Hessian to disc ---
      if (rdhess .and. (icom .eq. 1 .or. it .eq. nit) .or. .true.) then
        if (procid .eq. master) then
          ifi = fopna('HSSN',-1,4)
C   ... Hessian written in rdm-compatible format
          rewind ifi
          write(ifi) nvar,nvar,11
          write(ifi) w
          call fclose(ifi)
        endif
      endif
      if (procid .eq. master) call poppr
      call getpr(ipr)
C --- Periodically remove hessian ---
      if (nkill .gt. 0) then
        if (mod(it,nkill) .eq. 0) then
          call info(-20,0,0,
     .    '   ...  resetting hessian : iter=%i nkill=%i',it,nkill)
          ifi = fopna('HSSN',-1,4)
          call dfclos(ifi)
          call dcopy(nvar*nvar,0d0,0,w,1)
          call dcopy(nvar,1d0,0,w,nvar+1)
        endif
      endif
C --- Printout ---
      ifi = lgunit(2)
      lpp = (icom .eq. 1 .or. it .eq. nit) .and. ipr .ge. 30 .or.ipr .ge. 31
      if (natrlx .gt. 0 .and. lpp .and. lshr) then
        call dpzero(xv,6)
        call grdep2(1,nvar,indrlx,bas,xv)
        call info2(-20,0,0,' Update shear%N   PDEF=%6;8,4D%N STRAIN=%6;8,4D',bas,xv)
      elseif (natrlx .gt. 0 .and. lpp) then
        print 120
        print *,' Site   Class                      Position(relaxed)'
        do  70  j = 1, nbas
          is=ssite(j)%spec
          clablj = sspec(is)%name
          call icopy(size(ssite(j)%relax),ssite(j)%relax,1,ifrlx,1)
          write (stdo,130) j,clablj,(bas(ix,j),ifrlx(ix).eq.1,ix=1,3)
   70   continue
      endif

C ... Write new atom positions to LOG file, for use in CTRL file
      if (ipr .ge. 20 .and. .not. lshr) then
        dumstr = 'SITE'
        do  80  j = 1, nbas
          call icopy(size(ssite(j)%relax),ssite(j)%relax,1,ifrlx,1)
          if (j .eq. 2) dumstr = ' '
          is=ssite(j)%spec
          clablj=sspec(is)%name
          call icopy(size(ssite(j)%relax),ssite(j)%relax,1,ifrlx,1)
          
          call awrit2(dumstr//'ATOM='//clablj//
     .    ' POS=%3;12,7D  RELAX=%3,2i',' ',80,lgunit(2), bas(1,j),ifrlx)
   80   continue
      endif
      if (icom .eq. 0 .and. wkg(28) .lt. 0) icom = -1
      flush(stdo)
      flush(ifi)
  120 format(/' Updated atom positions:')
  130 format(i4,6x,a4,3x,3(f14.8,'(',l1,')'))
  160 format(10x,'DELTA=',3(f13.8),:,4(/31x,3(f13.8),:))
      call tcx('relax')
      end subroutine relax

      subroutine prelx1(mode,nm,lshr,natrlx,nvar,indrlx,p,bas)
C- Copy vector of variables to be minimized from/to bas
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :>0 copy p(*,mode) into bas(*)
Ci         :<0 copy bas(*) into p(*,mode)
Ci   nm    :which index to p that should be copied
Ci   lshr  :F p correspond to basis positions
Ci         :T p correspond to simple vector
Ci   natrlx:number of variables to relax
Ci   indrlx:permutation table mappng p -> bas
Ci   bas   :basis vectors
Ci   p     :vector of variables undergoing minimization process
Co Outputs
Cr Remarks
Cr
Cu Updates
Cu   09 Mar 06 Extracted from relax.f for use by other routines
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nm,natrlx,nvar,indrlx(2,natrlx)
      double precision p(nvar,nm),bas(3,*)
      logical lshr
C ... Local parameters
      integer i,j,iat,ix

C --- Update atom positions, or vice-versa ---
      if (.not. lshr) then
        j = 1
        do  60  i = 1, natrlx
          ix = indrlx(1,i)
          iat = indrlx(2,i)
          if (mode .gt. 0) then
            bas(ix,iat) = p(j,nm)
          else
            p(j,nm) = bas(ix,iat)
          endif
          j = j + 1
   60   continue
      else
        if (mode .gt. 0) then
          call dcopy(natrlx,p(1,nm),1,bas,1)
        else
          call dcopy(natrlx,bas,1,p(1,nm),1)
        endif

      endif
      end subroutine prelx1

      subroutine grdep2(i1,i2,indrlx,dstprm,dist)

C- Build up actual distortion from distortion parms
Cr Distortion is added into dist
C     implicit none
      integer i1,i2,indrlx(i2)
      double precision dstprm(i2),dist(6)
      integer i,ip

      do  i = i1, i2
        ip = indrlx(i)
C       1 + 2 + 3
        if (ip .eq. 1) then
          dist(1) = dist(1) + dstprm(i)
          dist(2) = dist(2) + dstprm(i)
          dist(3) = dist(3) + dstprm(i)
C       1 - 2
        else if (ip .eq. 2) then
          dist(1) = dist(1) + dstprm(i)
          dist(2) = dist(2) - dstprm(i)
C       1-3
        else if (ip .eq. 3) then
          dist(1) = dist(1) + dstprm(i)
          dist(3) = dist(3) - dstprm(i)
C       4, 5, 6
        else
          dist(ip) = dist(ip) + dstprm(i)
        endif
      enddo

      end subroutine grdep2


