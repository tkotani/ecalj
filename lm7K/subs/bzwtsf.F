!!== FSMOMMETHOD=1 == June2011 takao
!! This is for molecules.
!! 1st step : Set initial bias magnetic field under assuming all eigenvalues are discrete.
!!            Search elumo1,ehomo1 (for up spin) and elumo2,ehomo2 (for down spin) below.
!! 2nd step : refine the bias field for given temperature.
!!
!!  (Takao think this bzwtsf2 may need to be modified for solids).
      subroutine bzwtsf2(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,lswtk,
     .swtk,efermi,sumev,wtkb,qval,lwtkb,lfill,vnow)
      use m_globalvariables
C- BZ integration for fermi level, band sum and qp weights, fixed-spin
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbmx  :leading dimension of eb
Ci   nevx  :leading dimension of wtkb and max number of evals calculated
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   n1..n3:number of divisions for the k-point mesh
Ci   nkp   :number of inequivalent k-points (bzmesh.f)
Ci   ntet  :number of inequivalent tetrahedra (tetirr.f)
Ci   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
Ci         :corners of tetrahedron;
Ci         :idtet(0,i) number of tetrahedra of the i'th kind
Ci   zval  :valence charge
Ci   fmom  :fixed spin moment.  If zero, no constraint is applied.
Ci   metal :T => metal, F => nonmetal
Ci   tetra :T => tetrahedron integration
Ci   norder,npts,width,rnge: parameters for sampling integr. (maknos)
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   eb    :energy bands; alias eband
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   eb    :energy bands
Ci   efmax :eigenvalue limit for eigenvectors calculated in diagno
Ci   lswtk :Flags indicating whether 'spin weights' swtk are available
Ci   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
Ci         :where z are eigenvectors, sigma is the Pauli spin matrix
Ci         :Supplies information about spin moment in noncoll case.
Ci         :Used when lswtk is set
Cio  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
Cio        :0 weights are not available; no moment calculation
Cio        :if 1, weights were generated with no constraint
Cio        :In this case, print moment, and if fmom ne 0 remake weights
Cio        :with constraint; set to lwtkb=2 on output.
Cio        :if 2, weights were generated with constrained global moment
Cio        :if -1, same as 1
Co Outputs
Co   efermi:Fermi energy
Co   sumev :sum of eigenvalues
Co   wtkb  :integration weights (not generated for nonmetal case)
Co   qval  :qval(1) = total charge; qval(2) = magnetic moment
Cu Updates
Cu   12 Jul 08 change arg list in bzwts -- now returns entropy term
Cu   02 Jan 06 return qval (valence charge and moment)
Cu   22 Sep 01 Adapted from bzwts.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical metal,tetra
      integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,
     .idtet(5,ntet),lswtk,lwtkb
      double precision zval,fmom,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),
     .wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efmax,efermi,sumev,qval(2)
C Local variables
      integer ikp,ib,stdo,nglob,ipr,itmax,iter,iprint
      double precision amom,dosef(2),vhold(12),vnow,dvcap,dv,ef0,ent
      parameter (dvcap=.2d0,itmax=50)

      logical:: agreemom
      real(8),parameter::    NULLR =-99999
      integer:: nmom1,nmom2
      real(8):: ehomo1,ehomo2,elumo1,elumo2
      real(8),allocatable:: ebs(:,:,:)
      real(8):: ele1,ele2
      integer:: itermx
      logical:: quitvnow,lfill

C --- Fermi level without spin constraint ---
      call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,efermi,
     .sumev,wtkb,dosef,qval,ent,lfill)
      if (nsp .eq. 1) return

      call getpr(ipr)
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
C     stdl = nglob('stdl')
C --- Make and print out magnetic moment ---
      if ((lswtk .eq. 1 .or. nspc .eq. 1) .and. metal) then
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
        if (ipr .ge. 20) write(stdo,922) amom
  922   format(9x,'Mag. moment:',f15.6)
        qval(2) = amom
      else
        call info0(20,0,0,
     .  '%9fspin weights not available ... no spin moment calculated')
        return
      endif
!! --- Setup for fixed-spin moment method ---
!!      if (fmom .eq. 0 .or. lwtkb .eq. 0) return
      if (fmom==NULLR .or. lwtkb .eq. 0) return

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
c      vnow=0d0
      call info0(41,1,0,' Seek potential shift for fixed-spin mom ...')

!!== do loop for new guess at potential shift ==
!!bisection method takao
      itermx=100
      quitvnow=.false.
      do 10 iter=1,itermx
!! Potential shift
        allocate(ebs(nevx,2,nkp))
        if (nspc .eq. 2) then
          ebs = eb + vnow/2*swtk
        else
          ebs(:,1,:) = eb(:,1,:) - vnow/2d0
          ebs(:,2,:) = eb(:,2,:) + vnow/2d0
        endif
!! Fermi level with dv shift
        if(.not.quitvnow) call pshpr(ipr-50)
        call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .  metal,tetra,norder,npts,width,rnge,wtkp,ebs,efmax,efermi,
     .  sumev,wtkb,dosef,qval,ent,lfill)
        if (iprint()>= 20) then
          call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,
     .  amom)
          write(stdo,922) amom
        endif
        deallocate(ebs)
        if (.not.quitvnow) call poppr
        if(quitvnow) exit

!!=== Magnetic moment ===
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
        if (ipr .ge. 41) call awrit5(' -Vup+Vdn=%;10,6D  yields '//
     . 'ef=%;10,6D  amom=%;10,6D;  seeking %;10,6D',' ',160,
     . stdo,vnow,efermi,amom,fmom,vnow)
!!takao for molecule Dec1 2010
        agreemom= abs(amom-fmom) < 1d-6 ! 1d-6 on June-2 2011
        if(iprint()>60) print *,'ttttt amom fmom=',amom,fmom,agreemom
c      if(original_dvdos) then
        call dvdos(vnow,amom,dosef,vhold,fmom,dvcap,dv)
        if(agreemom) vhold(12)=0
c      if (abs(dv) .lt. 1d-6) then
        quitvnow=.false.
        if (abs(dv).lt.1d-6 .or. agreemom) then
C       A root was found
          if (vhold(12) .eq. -2 .or. vhold(12) .eq. -3 .or.
     .  vhold(12) .eq.  0 .or. vhold(12) .eq.  1) then
            if (ipr .ge. 10)
     .    call awrit5('%N BZWTSF2: potential shift bracketed.'//
     .    '  Unconstrained efermi=%,6;6d.'//
     .    '%N constraint fmom=%,6;6d  actual mmom=%,6;6d'//
     .    '  ef=%,6;6d  -Vup+Vdn=%,6;6d',' ',160,stdo,ef0,
     .    fmom,amom,efermi,vnow)
            quitvnow=.true.
          endif
        else if (iter .eq. itmax) then
          if (ipr .ge. 10)
     .  call awrit5('%N BZWTSF2: failed to converge potential shift'//
     .  ' after %i iterations.'//
     .  '%N constraint fmom=%,6;6d  actual amom=%,6;6d'//
     .  '  ef=%,6;6d  -Vup+Vdn=%,6;6d',' ',160,stdo,iter,fmom,amom,efermi,
     .  vnow)
          quitvnow=.true.
        endif
 10   continue
C--------------------------------------

cccccccccccccccccccccccccccccccccc
!! weitht modification by hand for Vanadium atom.
c      wtkb(8,1,1)=0d0                                                            
c      wtkb(9,1,1)=1d0                                                            
c      wtkb(10,1,1)=1d0                                                           
c      do ie=1,20                                                                 
c       write(*,"('vvvvvv up:',i3,f13.6,' ',2f13.4)")ie,eb(ie,1,1),wtkb(ie,1,1)   
c      enddo                                                                      
c      do ie=1,20                                                                 
c       write(*,"('vvvvvv dn:',i3,f13.6,' ',2f13.4)")ie,eb(ie,2,1),wtkb(ie,2,1)   
c      enddo                                                                      
ccccccccccccccccccccccccccccccccccc              
c
      ele1 = (zval+amom)/2d0
      ele2 = (zval-amom)/2d0
      sumev=sumev+ ele1*(vnow/2d0) - ele2*(vnow/2d0) !sumev correction takao
      if(iprint()>20) write(stdo,"(' bzwtsf2(METHOD=1): Set Bias field -Vup+Vdn=',f20.15)")vnow
      if (lswtk.eq.1 .and. lwtkb .eq. 1) then
        lwtkb = 2
      elseif (lswtk.eq.1 .and. lwtkb .eq. 2) then
        lwtkb = 1
      endif
      call tcx('bzwtsf2')
      end subroutine bzwtsf2


!!== FSMOMMETHOD=0 ogiginal version(modified version. fmom=0 is allowed.)==
      subroutine bzwtsf(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,lswtk,
     .swtk,efermi,sumev,wtkb,qval,lwtkb,lfill,vnow)
      use m_globalvariables
C- BZ integration for fermi level, band sum and qp weights, fixed-spin
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbmx  :leading dimension of eb
Ci   nevx  :leading dimension of wtkb and max number of evals calculated
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   n1..n3:number of divisions for the k-point mesh
Ci   nkp   :number of inequivalent k-points (bzmesh.f)
Ci   ntet  :number of inequivalent tetrahedra (tetirr.f)
Ci   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
Ci         :corners of tetrahedron;
Ci         :idtet(0,i) number of tetrahedra of the i'th kind
Ci   zval  :valence charge
Ci   fmom  :fixed spin moment.  Even zero is allowed.
Ci   metal :T => metal, F => nonmetal
Ci   tetra :T => tetrahedron integration
Ci   norder,npts,width,rnge: parameters for sampling integr. (maknos)
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   eb    :energy bands; alias eband
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   eb    :energy bands
Ci   efmax :eigenvalue limit for eigenvectors calculated in diagno
Ci   lswtk :Flags indicating whether 'spin weights' swtk are available
Ci   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
Ci         :where z are eigenvectors, sigma is the Pauli spin matrix
Ci         :Supplies information about spin moment in noncoll case.
Ci         :Used when lswtk is set
Cio  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
Cio        :0 weights are not available; no moment calculation
Cio        :if 1, weights were generated with no constraint
Cio        :In this case, print moment, and if fmom ne 0 remake weights
Cio        :with constraint; set to lwtkb=2 on output.
Cio        :if 2, weights were generated with constrained global moment
Cio        :if -1, same as 1
Co Outputs
Co   efermi:Fermi energy
Co   sumev :sum of eigenvalues
Co   wtkb  :integration weights (not generated for nonmetal case)
Co   qval  :qval(1) = total charge; qval(2) = magnetic moment
Cu Updates
Cu   12 Jul 08 change arg list in bzwts -- now returns entropy term
Cu   02 Jan 06 return qval (valence charge and moment)
Cu   22 Sep 01 Adapted from bzwts.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical metal,tetra
      integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,
     .idtet(5,ntet),lswtk,lwtkb
      double precision zval,fmom,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),
     .wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efmax,efermi,sumev,qval(2)
C Local variables
      integer ikp,ib,stdo,nglob,ipr,itmax,iter,iprint
      double precision amom,dosef(2),vhold(12),vnow,dvcap,dv,ef0,ent
      parameter (dvcap=.2d0,itmax=50)

      real(8):: ele1,ele2
      integer:: itermx
      logical:: agreemom
      real(8),parameter::    NULLR =-99999
      real(8),allocatable:: ebs(:,:,:)
      logical:: quitvnow,lfill

!!== Fermi level without spin constraint ==
      call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,efermi,
     .sumev,wtkb,dosef,qval,ent,lfill)
      if (nsp .eq. 1) return

      call getpr(ipr)
      stdo = globalvariables%stdo

!!== Make and print out magnetic moment ==
      if ((lswtk .eq. 1 .or. nspc .eq. 1) .and. metal) then
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
        if (ipr .ge. 20) write(stdo,922) amom
  922   format(9x,'Mag. moment:',f15.6)
        qval(2) = amom
      else
        call info0(20,0,0,
     .  '%9fspin weights not available ... no spin moment calculated')
        return
      endif

!!== Setup for fixed-spin moment method ==
c      if (fmom .eq. 0 .or. lwtkb .eq. 0) return
      if (fmom==NULLR .or. lwtkb .eq. 0) return
      call tcn('bzwtsf')
      call dpzero(vhold,12)
      vnow = 0
      ef0 = efermi
      call info0(41,1,0,' Seek potential shift for fixed-spin mom ...')

!!== do loop for new guess at potential shift ==
!!bisection method takao
      itermx=100
      do 10 iter=1,itermx
!!=== Magnetic moment ===
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
        if (ipr .ge. 41) call awrit5(' -Vup+Vdn=%;10,6D  yields '//
     . 'ef=%;10,6D  amom=%;10,6D;  seeking %;10,6D',' ',160,
     . stdo,vnow,efermi,amom,fmom,vnow)
!!takao for molecule Dec1 2010
        agreemom= abs(amom-fmom) < 1d-3
        if(iprint()>60) print *,'ttttt amom fmom=',amom,fmom,agreemom
c      if(original_dvdos) then
        call dvdos(vnow,amom,dosef,vhold,fmom,dvcap,dv)
        if(agreemom) vhold(12)=0
c      if (abs(dv) .lt. 1d-6) then
        quitvnow=.false.
        if (abs(dv).lt.1d-6 .or. agreemom) then
C       A root was found
          if (vhold(12) .eq. -2 .or. vhold(12) .eq. -3 .or.
     .  vhold(12) .eq.  0 .or. vhold(12) .eq.  1) then
            if (ipr .ge. 10)
     .    call awrit5('%N BZWTSF: potential shift bracketed.'//
     .    '  Unconstrained efermi=%,6;6d.'//
     .    '%N constraint fmom=%,6;6d  actual mmom=%,6;6d'//
     .    '  ef=%,6;6d  -Vup+Vdn=%,6;6d',' ',160,stdo,ef0,
     .    fmom,amom,efermi,vnow)
            quitvnow=.true.
          endif
        else if (iter .eq. itmax) then
          if (ipr .ge. 10)
     .  call awrit5('%N BZWTSF: failed to converge potential shift'//
     .  ' after %i iterations.'//
     .  '%N constraint fmom=%,6;6d  actual amom=%,6;6d'//
     .  '  ef=%,6;6d  -Vup+Vdn=%,6;6d',' ',160,stdo,iter,fmom,amom,efermi,
     .  vnow)
          quitvnow=.true.
        endif

!! Potential shift
        allocate(ebs(nevx,2,nkp))
        if (nspc .eq. 2) then
          ebs = eb + vnow/2*swtk
        else
          ebs(:,1,:) = eb(:,1,:) - vnow/2
          ebs(:,2,:) = eb(:,2,:) + vnow/2
        endif

!! Fermi level with dv shift
        if(.not.quitvnow) call pshpr(ipr-50)
        call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .  metal,tetra,norder,npts,width,rnge,wtkp,ebs,efmax,efermi,
     .  sumev,wtkb,dosef,qval,ent,lfill)
        if (iprint()>= 20) then
          call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,
     .  amom)
          write(stdo,922) amom
        endif
        deallocate(ebs)
        if (.not.quitvnow) call poppr
        if(quitvnow) exit
 10   continue
c
      ele1 = (zval+fmom)/2d0
      ele2 = (zval-fmom)/2d0
      sumev=sumev+ ele1*(vnow/2d0) - ele2*(vnow/2d0) !sumev correction takao
      if(iprint()>20) write(stdo,"(' bzwtsf: Set Bias field -Vup+Vdn=',f20.15)")vnow
      if (lswtk.eq.1 .and. lwtkb .eq. 1) then
        lwtkb = 2
      elseif (lswtk.eq.1 .and. lwtkb .eq. 2) then
        lwtkb = 1
      endif
      call tcx('bzwtsf')
      end


!!== FSMOMMETHOD unsed now==
      subroutine bzwtsf3(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,lswtk,
     .swtk,efermi,sumev,wtkb,qval,lwtkb)
      use m_globalvariables
!!t.kotani from bzwtsf. now fs
C- BZ integration for fermi level, band sum and qp weights, fixed-spin
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbmx  :leading dimension of eb
Ci   nevx  :leading dimension of wtkb and max number of evals calculated
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   n1..n3:number of divisions for the k-point mesh
Ci   nkp   :number of inequivalent k-points (bzmesh.f)
Ci   ntet  :number of inequivalent tetrahedra (tetirr.f)
Ci   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
Ci         :corners of tetrahedron;
Ci         :idtet(0,i) number of tetrahedra of the i'th kind
Ci   zval  :valence charge
Ci   fmom  :fixed spin moment.  Even for zero, it works.
Ci   metal :T => metal, F => nonmetal
Ci   tetra :T => tetrahedron integration
Ci   norder,npts,width,rnge: parameters for sampling integr. (maknos)
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   eb    :energy bands; alias eband
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   eb    :energy bands
Ci   efmax :eigenvalue limit for eigenvectors calculated in diagno
Ci   lswtk :Flags indicating whether 'spin weights' swtk are available
Ci   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
Ci         :where z are eigenvectors, sigma is the Pauli spin matrix
Ci         :Supplies information about spin moment in noncoll case.
Ci         :Used when lswtk is set
Cio  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
Cio        :0 weights are not available; no moment calculation
Cio        :if 1, weights were generated with no constraint
Cio        :In this case, print moment, and if fmom ne 0 remake weights
Cio        :with constraint; set to lwtkb=2 on output.
Cio        :if 2, weights were generated with constrained global moment
Cio        :if -1, same as 1
Co Outputs
Co   efermi:Fermi energy
Co   sumev :sum of eigenvalues
Co   wtkb  :integration weights (not generated for nonmetal case)
Co   qval  :qval(1) = total charge; qval(2) = magnetic moment
Cu Updates
Cu   12 Jul 08 change arg list in bzwts -- now returns entropy term
Cu   02 Jan 06 return qval (valence charge and moment)
Cu   22 Sep 01 Adapted from bzwts.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical metal,tetra
      integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,
     .idtet(5,ntet),lswtk,lwtkb
      double precision zval,fmom,fmomx,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),
     .wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efmax,efermi,sumev,qval(2)
C Local variables
      integer ikp,ib,stdo,nglob,ipr,iter,iprint
      double precision amom,dosef(2),vhold(12),vnow,dvcap,dv,ef0,ent
      integer,parameter:: itmax=150

c      integer:: itermx
      logical:: agreemom
      real(8),parameter::    NULLR =-99999
      real(8),allocatable:: ebs(:,:,:)
      logical:: quitvnow,lfill
      real(8):: vnowmin,vnowmax,widthx,vnow1,vnow2
      integer:: ix
      call tcn('bzwtsf3')
      call getpr(ipr)
      stdo = globalvariables%stdo

!!== Fermi level without spin constraint ==
      call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,efermi,
     .sumev,wtkb,dosef,qval,ent,lfill)
      if (nsp .eq. 1) return

!!== Make and print out magnetic moment ==
      if ((lswtk .eq. 1 .or. nspc .eq. 1) .and. metal) then
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
        if (ipr .ge. 20) write(stdo,922) amom
  922   format(9x,'Mag. moment:',f15.6)
        qval(2) = amom
      else
        call info0(20,0,0,
     .  '%9fspin weights not available ... no spin moment calculated')
        return
      endif

!!== Setup for fixed-spin moment method ==
c      if (fmom .eq. 0 .or. lwtkb .eq. 0) return
      if (fmom==NULLR .or. lwtkb .eq. 0) return
      call dpzero(vhold,12)
      vnow = 0d0
      ef0 = efermi
      call info0(41,1,0,' Seek potential shift for fixed-spin mom ...')
c      dv1=-5d0!min
c      dv2=5d0 !max
c      dv0=0.1
      if(amom>fmom) then
        vnowmax=vnow
        vnowmin=-5d0
      else
        vnowmax=5d0
        vnowmin=vnow
      endif

!!== do loop for new guess at potential shift ==
!!bisection method takao
c      bzwtsf_bisection=.true.
      allocate(ebs(nevx,2,nkp))
      do 20 ix=1,3
        if(ix==1) fmomx=fmom+1d0 !search vnow1 for fsmom to fsmom+1
        if(ix==2) fmomx=fmom     !search vnow2 for fsmom-1 to fsmom
        quitvnow=.false.
        do 10 iter=1,itmax
          vnow=(vnowmax+vnowmin)/2d0 !trial vnow bisection method. new mode if this exist
!!=== Potential shift ===
          if (nspc .eq. 2) then
            ebs = eb + vnow/2*swtk
          else
            ebs(:,1,:) = eb(:,1,:) - vnow/2d0
            ebs(:,2,:) = eb(:,2,:) + vnow/2d0
          endif
!!=== Magnetic moment ===
          if(.not.quitvnow) call pshpr(ipr-50)
c        widthx=width
c        if(.not.quitvnow) widthx=0.01d0
          call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .   metal,tetra,norder,npts,width,rnge,wtkp,ebs,efmax,efermi,
     .   sumev,wtkb,dosef,qval,ent,lfill)
          call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
          if(.not.quitvnow) call poppr()
c       if (ipr .ge. 41) call awrit5(' dv=%;10,6D  yields '//
c     . 'ef=%;10,6D  amom=%;10,6D;  seeking %;10,6D',' ',160,
c     . stdo,vnow,efermi,amom,fmom,vnow)
          if(ix==3.or.iter==itmax) goto 11
          if(amom>fmomx) then
            vnowmax=vnow
          else
            vnowmin=vnow
          endif
!!takao for molecule Dec1 2010
          if(iprint()>60) print *,'ttttt iter fmomx amom=',iter,fmomx,amom,vnowmax,vnowmin
          if(amom>fmomx-1d-12.and.amom<fmomx+1d-12) quitvnow=.true.
          if(quitvnow) goto 11
 10     continue
        call awrit5('%N BZWTSF: failed to converge potential shift'//
     .  ' after %i iterations.'//
     .  '%N constraint fmom=%,6;6d  actual amom=%,6;6d'//
     .  '  ef=%,6;6d  -Vup+Vdn=%,6;6d',' ',160,stdo,iter,fmom,amom,efermi,
     .  vnow)
 11     continue
        if(ix==1) then
          vnow1  = vnow
          vnowmin= vnow-5d0
          vnowmax= vnow+0.1d0
        elseif(ix==2) then
          vnow2  = vnow
          vnowmax= vnow1
          vnowmin= vnow2
        endif
        if(iprint()>60) print *,' qqqqq ix=',ix, 'vnowmin vnowmax amom=',vnowmin,vnowmax,amom
 20   continue
      deallocate(ebs)
      if (lswtk.eq.1 .and. lwtkb .eq. 1) then
        lwtkb = 2
      elseif (lswtk.eq.1 .and. lwtkb .eq. 2) then
        lwtkb = 1
      endif
      call tcx('bzwtsf3')
      end






      subroutine bzwtsm(lswtk,nkp,nsp,nevx,wtkb,swtk,amom)
C- Determine the magnetic moment, collinear or noncollinear case
C ----------------------------------------------------------------------
Ci Inputs
Ci   lswtk :if true, swtk is used.  Otherwise, collinear case assumed:
Ci         :swtk(*,1,*) = 1  and swtk(*,2,*) = -1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nevx  :Maximum number of bands
Ci   wtkb  :band weights
Ci   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
Ci         :where z are eigenvectors, sigma is the Pauli spin matrix
Ci         :Used when lswtk is set
Co Outputs
Co   amom  :magnetic moment
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   09 Jun 07
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lswtk
      integer nkp,nevx,nsp
      double precision wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),amom
C ... Local parameters
      integer ikp,k
      double precision dsum

      if (nsp .eq. 1) return

C      if (lswtk .eq. 1 .and. nspc .eq. 2) then
      if (lswtk) then
        amom = 0
        do  ikp = 1, nkp
          do  k = 1, nevx
            amom = amom + wtkb(k,1,ikp)*swtk(k,1,ikp)
     .      + wtkb(k,2,ikp)*swtk(k,2,ikp)
          enddo
        enddo
      else
        amom = 0
        do  ikp = 1, nkp
          amom = amom + dsum(nevx,wtkb(1,1,ikp),1) -
     .    dsum(nevx,wtkb(1,2,ikp),1)
        enddo
      endif
      end

