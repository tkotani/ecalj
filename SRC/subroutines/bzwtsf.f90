!!== FSMOMMETHOD=1 == June2011 takao
!! This is for molecules.
!! 1st step : Set initial bias magnetic field under assuming all eigenvalues are discrete.
!!            Search elumo1,ehomo1 (for up spin) and elumo2,ehomo2 (for down spin) below.
!! 2nd step : refine the bias field for given temperature.
!!
!!  (Takao think this bzwtsf2 may need to be modified for solids).
subroutine bzwtsf2(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
     fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,lswtk, &
     swtk,efermi,sumev,wtkb,qval,lwtkb,lfill,vnow)
  use m_lmfinit,only: stdo
  use m_ftox

  !- BZ integration for fermi level, band sum and qp weights, fixed-spin
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbmx  :leading dimension of eb
  !i   nevx  :leading dimension of wtkb and max number of evals calculated
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
  !i   n1..n3:number of divisions for the k-point mesh
  !i   nkp   :number of inequivalent k-points (bzmesh.f)
  !i   ntet  :number of inequivalent tetrahedra (tetirr.f)
  !i   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
  !i         :corners of tetrahedron;
  !i         :idtet(0,i) number of tetrahedra of the i'th kind
  !i   zval  :valence charge
  !i   fmom  :fixed spin moment.  If zero, no constraint is applied.
  !i   metal :T => metal, F => nonmetal
  !i   tetra :T => tetrahedron integration
  !i   norder,npts,width,rnge: parameters for sampling integr. (maknos)
  !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
  !i   eb    :energy bands; alias eband
  !i   eb    :energy bands
  !i   lswtk :Flags indicating whether 'spin weights' swtk are available
  !i   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
  !i         :where z are eigenvectors, sigma is the Pauli spin matrix
  !i         :Supplies information about spin moment in noncoll case.
  !i         :Used when lswtk is set
  ! oxxx  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
  ! o        :0 weights are not available; no moment calculation
  ! o        :if 1, weights were generated with no constraint
  ! o        :In this case, print moment, and if fmom ne 0 remake weights
  ! o        :with constraint; set to lwtkb=2 on output.
  ! o        :if 2, weights were generated with constrained global moment
  ! o        :if -1, same as 1
  !o Outputs
  !o   efermi:Fermi energy
  !o   sumev :sum of eigenvalues
  !o   wtkb  :integration weights (not generated for nonmetal case)
  !o   qval  :qval(1) = total charge; qval(2) = magnetic moment
  !u Updates
  !u   12 Jul 08 change arg list in bzwts -- now returns entropy term
  !u   02 Jan 06 return qval (valence charge and moment)
  !u   22 Sep 01 Adapted from bzwts.
  ! ----------------------------------------------------------------------
  implicit none
  ! Passed parameters
  logical :: metal,tetra
  integer :: nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet, &
       idtet(5,ntet),lswtk,lwtkb
  double precision :: zval,fmom,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp), &
       wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efermi,sumev,qval(2)
  ! Local variables
  integer :: ikp,ib,ipr,itmax,iter,iprint
  double precision :: amom,dosef(2),vhold(12),vnow,dvcap,dv,ef0,ent
  parameter (dvcap=.2d0,itmax=50)

  logical:: agreemom
  real(8),parameter::    NULLR =-99999
  integer:: nmom1,nmom2
  real(8):: ehomo1,ehomo2,elumo1,elumo2
  real(8),allocatable:: ebs(:,:,:)
  real(8):: ele1,ele2
  integer:: itermx
  logical:: quitvnow,lfill

  ! --- Fermi level without spin constraint ---
  call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
       metal,tetra,norder,npts,width,rnge,wtkp,eb,efermi, &
       sumev,wtkb,dosef,qval,ent,lfill)
  if (nsp == 1) return

  call getpr(ipr)
  ! angenglob      stdo = nglob('stdo')
  !      stdo = globalvariables%stdo
  !     stdl = nglob('stdl')
  ! --- Make and print out magnetic moment ---
  if ((lswtk == 1 .OR. nspc == 1) .AND. metal) then
     call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
     if (ipr >= 20) write(stdo,922) amom
922  format(9x,'Mag. moment:',f15.6)
     qval(2) = amom
  else
     write(stdo,*) 'spin weights not available ... no spin moment calculated'
     return
  endif
  !! --- Setup for fixed-spin moment method ---
  if (fmom==NULLR .OR. lwtkb == 0) return
  !      if (fmom==NULLR .or. lwtkb .eq. 0) call rx('bzwtsf2: check logic of lwtkb')

  call tcn('bzwtsf2')

  ele1 = (zval+fmom)/2d0
  ele2 = (zval-fmom)/2d0
  nmom1 = ele1+1d-8
  nmom2 = ele2+1d-8
  elumo1= minval(eb(nmom1+1,1,:))
  elumo2= minval(eb(nmom2+1,2,:))
  if(nmom1/=0) then
     ehomo1= maxval(eb(nmom1,1,:))
  else
     ehomo1= elumo1 - 0.5d0
  endif
  if(nmom2/=0) then
     ehomo2= maxval(eb(nmom2,2,:))
  else
     ehomo2= elumo2 - 0.5d0
  endif
  vnow = (ehomo1+elumo1)/2d0 -(ehomo2+elumo2)/2d0
  efermi= ((ehomo1+elumo1)/2d0 +(ehomo2+elumo2)/2d0)/2d0 !/2d0 bug fix Jun26,2014 this only affects to the message.
  write(stdo,"('bzwtsf2: zval fmom nmon1 nmom2=',2f12.8,2x,2i3)")zval,fmom,nmom1,nmom2
  write(stdo,"('bzwtsf2: HOMOup LUMOup Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo1,elumo1,elumo1-ehomo1
  write(stdo,"('bzwtsf2: HOMOdn LUMOdn Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo2,elumo2,elumo2-ehomo2
  write(stdo,"('bzwtsf2: Set Bias initial cond. -Vup+Vdn=',f20.15)")vnow

  !!= takao interted a block taken from original version of bzwtsf.F June-2 2011.=
  vhold= 0d0
  ef0  = efermi
  write(stdo,*)' Seek potential shift for fixed-spin mom ...'

  !!== do loop for new guess at potential shift ==
  ! bisection method takao
  itermx=100
  quitvnow=.false.
  do 10 iter=1,itermx
     !! Potential shift
     allocate(ebs(nevx,2,nkp))
     if (nspc == 2) then
        ebs = eb + vnow/2*swtk
     else
        ebs(:,1,:) = eb(:,1,:) - vnow/2d0
        ebs(:,2,:) = eb(:,2,:) + vnow/2d0
     endif
     !! Fermi level with dv shift
     if( .NOT. quitvnow) call pshpr(ipr-50)
     call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
          metal,tetra,norder,npts,width,rnge,wtkp,ebs,efermi, &
          sumev,wtkb,dosef,qval,ent,lfill)
     if (iprint()>= 20) then
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk, &
             amom)
        write(stdo,922) amom
     endif
     deallocate(ebs)
     if ( .NOT. quitvnow) call poppr
     if(quitvnow) exit
     !!=== Magnetic moment ===
     call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
     if(ipr>=41)write(stdo,ftox)' -Vup+Vdn=',ftof(vnow,8),'yields ef=',ftof(efermi), &
          'amom',ftof(amom),'when seeking mom=',ftof(fmom)
     ! takao for molecule Dec1 2010
     agreemom= abs(amom-fmom) < 1d-6 ! 1d-6 on June-2 2011
     if(iprint()>60) print *,'ttttt amom fmom=',amom,fmom,agreemom
     call dvdos(vnow,amom,dosef,vhold,fmom,dvcap,dv)
     if(agreemom) vhold(12)=0
     quitvnow=.false.
     if (abs(dv) < 1d-6 .OR. agreemom) then
        if (vhold(12) == -2 .OR. vhold(12) == -3 .OR. &
             vhold(12) ==  0 .OR. vhold(12) ==  1) then
           if (ipr >= 10) &
                write(stdo,ftox)' BZWTSF2: potential shift bracketed.', &
                'Unconstrained efermi=',ftof(ef0), &
                'constraint fmom=',ftof(fmom),'actual mmom=',amom, &
                'ef=',efermi,'-Vup+Vdn=',vnow
           quitvnow=.true.
        endif
     else if (iter == itmax) then
        if (ipr >= 10) &
             write(stdo,ftox)' BZWTSF2: failed to converge potential shift', &
             'after',iter,'iterations.', &
             'constraint fmom=',ftof(fmom),'actual amom=',ftof(amom), &
             'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vnow)
        quitvnow=.true.
     endif
10 enddo
  ele1 = (zval+amom)/2d0
  ele2 = (zval-amom)/2d0
  sumev=sumev+ ele1*(vnow/2d0) - ele2*(vnow/2d0) !sumev correction takao
  if(iprint()>20) write(stdo,ftox)' bzwtsf2(METHOD=1): Set Bias field -Vup+Vdn=',ftof(vnow,8)
  if (lswtk == 1 .AND. lwtkb == 1) then
     call rx('bzwtsf2:111 tk think not used here')
  elseif (lswtk == 1 .AND. lwtkb == 2) then
     call rx('bzwtsf2:222 tk think not used here')
  endif
  call tcx('bzwtsf2')
end subroutine bzwtsf2


!!== FSMOMMETHOD=0 ogiginal version(modified version. fmom=0 is allowed.)==
subroutine bzwtsf(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
     fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,lswtk, &
     swtk,efermi,sumev,wtkb,qval,lwtkb,lfill,vnow)
  use m_lmfinit,only: stdo
  use m_ftox

  !- BZ integration for fermi level, band sum and qp weights, fixed-spin
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbmx  :leading dimension of eb
  !i   nevx  :leading dimension of wtkb and max number of evals calculated
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
  !i   n1..n3:number of divisions for the k-point mesh
  !i   nkp   :number of inequivalent k-points (bzmesh.f)
  !i   ntet  :number of inequivalent tetrahedra (tetirr.f)
  !i   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
  !i         :corners of tetrahedron;
  !i         :idtet(0,i) number of tetrahedra of the i'th kind
  !i   zval  :valence charge
  !i   fmom  :fixed spin moment.  Even zero is allowed.
  !i   metal :T => metal, F => nonmetal
  !i   tetra :T => tetrahedron integration
  !i   norder,npts,width,rnge: parameters for sampling integr. (maknos)
  !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
  !i   eb    :energy bands; alias eband
  !i   eb    :energy bands
  !i   lswtk :Flags indicating whether 'spin weights' swtk are available
  !i   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
  !i         :where z are eigenvectors, sigma is the Pauli spin matrix
  !i         :Supplies information about spin moment in noncoll case.
  !i         :Used when lswtk is set
  ! oxxx  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
  ! o        :0 weights are not available; no moment calculation
  ! o        :if 1, weights were generated with no constraint
  ! o        :In this case, print moment, and if fmom ne 0 remake weights
  ! o        :with constraint; set to lwtkb=2 on output.
  ! o        :if 2, weights were generated with constrained global moment
  ! o        :if -1, same as 1
  !o Outputs
  !o   efermi:Fermi energy
  !o   sumev :sum of eigenvalues
  !o   wtkb  :integration weights (not generated for nonmetal case)
  !o   qval  :qval(1) = total charge; qval(2) = magnetic moment
  !u Updates
  !u   12 Jul 08 change arg list in bzwts -- now returns entropy term
  !u   02 Jan 06 return qval (valence charge and moment)
  !u   22 Sep 01 Adapted from bzwts.
  ! ----------------------------------------------------------------------
  implicit none
  ! Passed parameters
  logical :: metal,tetra
  integer :: nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet, &
       idtet(5,ntet),lswtk,lwtkb
  double precision :: zval,fmom,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp), &
       wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efermi,sumev,qval(2)
  ! Local variables
  integer :: ikp,ib,ipr,itmax,iter,iprint
  double precision :: amom,dosef(2),vhold(12),vnow,dvcap,dv,ef0,ent
  parameter (dvcap=.2d0,itmax=50)

  real(8):: ele1,ele2
  integer:: itermx
  logical:: agreemom
  real(8),parameter::    NULLR =-99999
  real(8),allocatable:: ebs(:,:,:)
  logical:: quitvnow,lfill

  !!== Fermi level without spin constraint ==
  call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
       metal,tetra,norder,npts,width,rnge,wtkp,eb,efermi, &
       sumev,wtkb,dosef,qval,ent,lfill)
  if (nsp == 1) return

  call getpr(ipr)
  !      stdo = globalvariables%stdo

  !!== Make and print out magnetic moment ==
  if ((lswtk == 1 .OR. nspc == 1) .AND. metal) then
     call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
     if (ipr >= 20) write(stdo,922) amom
922  format(9x,'Mag. moment:',f15.6)
     qval(2) = amom
  else
     write(stdo,*)'spin weights not available ... no spin moment calculated'
     return
  endif

  !!== Setup for fixed-spin moment method ==
  !      if (fmom .eq. 0 .or. lwtkb .eq. 0) return
  if (fmom==NULLR .OR. lwtkb == 0) return
  !      if (fmom==NULLR .or. lwtkb .eq. 0) call rx('bzwtsf:check logic of lwtkb')
  call tcn('bzwtsf')
  call dpzero(vhold,12)
  vnow = 0
  ef0 = efermi
  write(stdo,*)' Seek potential shift for fixed-spin mom ...'

  !!== do loop for new guess at potential shift ==
  ! bisection method takao
  itermx=100
  do 10 iter=1,itermx
     !!=== Magnetic moment ===
     call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
     if(ipr>=41) write(stdo,ftox)' -Vup+Vdn=',ftof(vnow,8),'yields ', &
          'ef=',ftof(efermi),'amom=',ftof(amom),'when seeking',ftof(fmom)
     ! takao for molecule Dec1 2010
     agreemom= abs(amom-fmom) < 1d-3
     if(iprint()>60) print *,'ttttt amom fmom=',amom,fmom,agreemom
     !      if(original_dvdos) then
     call dvdos(vnow,amom,dosef,vhold,fmom,dvcap,dv)
     if(agreemom) vhold(12)=0
     !      if (abs(dv) .lt. 1d-6) then
     quitvnow=.false.
     if (abs(dv) < 1d-6 .OR. agreemom) then
        !       A root was found
        if (vhold(12) == -2 .OR. vhold(12) == -3 .OR. &
             vhold(12) ==  0 .OR. vhold(12) ==  1) then
           if (ipr >= 10) write(stdo,ftox)' BZWTSF: potential shift bracketed.', &
                'Unconstrained efermi=',ftof(ef0), &
                'constraint fmom=',ftof(fmom),'actual mmom=',ftof(amom), &
                'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vnow,8)
           quitvnow=.true.
        endif
     else if (iter == itmax) then
        if(ipr>=10)then
           write(stdo,ftox)' BZWTSF: failed to converge potential shift after',iter,'iterations.'
           write(stdo,ftox)' constraint fmom=',ftof(fmom),'actual amom=',ftof(amom), &
                'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vnow,8)
        endif
        quitvnow=.true.
     endif

     !! Potential shift
     allocate(ebs(nevx,2,nkp))
     if (nspc == 2) then
        ebs = eb + vnow/2*swtk
     else
        ebs(:,1,:) = eb(:,1,:) - vnow/2
        ebs(:,2,:) = eb(:,2,:) + vnow/2
     endif

     !! Fermi level with dv shift
     if( .NOT. quitvnow) call pshpr(ipr-50)
     call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
          metal,tetra,norder,npts,width,rnge,wtkp,ebs,efermi, &
          sumev,wtkb,dosef,qval,ent,lfill)
     if (iprint()>= 20) then
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk, &
             amom)
        write(stdo,922) amom
     endif
     deallocate(ebs)
     if ( .NOT. quitvnow) call poppr
     if(quitvnow) exit
10 enddo

  ele1 = (zval+fmom)/2d0
  ele2 = (zval-fmom)/2d0
  sumev=sumev+ ele1*(vnow/2d0) - ele2*(vnow/2d0) !sumev correction takao
  if(iprint()>20) write(stdo,"(' bzwtsf: Set Bias field -Vup+Vdn=',f20.15)")vnow
  if (lswtk == 1 .AND. lwtkb == 1) then
     !        lwtkb = 2
     call rx('bzwtsf:111 tk think not used here')
  elseif (lswtk == 1 .AND. lwtkb == 2) then
     !        lwtkb = 1
     call rx('bzwtsf:222 tk think not used here')
  endif
  call tcx('bzwtsf')
end subroutine bzwtsf


subroutine bzwtsm(lswtk,nkp,nsp,nevx,wtkb,swtk,amom)
  use m_ftox
  !- Determine the magnetic moment, collinear or noncollinear case
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lswtk :if true, swtk is used.  Otherwise, collinear case assumed:
  !i         :swtk(*,1,*) = 1  and swtk(*,2,*) = -1
  !i   nkp   :number of irreducible k-points (bzmesh.f)
  !i   nevx  :Maximum number of bands
  !i   wtkb  :band weights
  !i   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
  !i         :where z are eigenvectors, sigma is the Pauli spin matrix
  !i         :Used when lswtk is set
  !o Outputs
  !o   amom  :magnetic moment
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   09 Jun 07
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  logical :: lswtk
  integer :: nkp,nevx,nsp
  double precision :: wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),amom
  ! ... Local parameters
  integer :: ikp,k
  double precision :: dsum

  if (nsp == 1) return

  !      if (lswtk .eq. 1 .and. nspc .eq. 2) then
  if (lswtk) then
     amom = 0
     do  ikp = 1, nkp
        do  k = 1, nevx
           amom = amom + wtkb(k,1,ikp)*swtk(k,1,ikp) &
                + wtkb(k,2,ikp)*swtk(k,2,ikp)
        enddo
     enddo
  else
     amom = 0
     do  ikp = 1, nkp
        amom = amom + dsum(nevx,wtkb(1,1,ikp),1) - &
             dsum(nevx,wtkb(1,2,ikp),1)
     enddo
  endif
end subroutine bzwtsm

