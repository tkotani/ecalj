      subroutine bzwtsf(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,lswtk,
     .swtk,efermi,sumev,wtkb,qval,lwtkb)
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
C     implicit none
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

#ifdef KINODEBUG
      write(*,*) 'kino bzwtsf in'
#endif
C --- Fermi level without spin constraint ---
      call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,efermi,
     .sumev,wtkb,dosef,qval,ent)
      if (nsp .eq. 1) return

      call getpr(ipr)
      stdo = nglob('stdo')
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

C --- Setup for fixed-spin moment method ---
      if (fmom .eq. 0 .or. lwtkb .eq. 0) return
      call dpzero(vhold,12)
      vnow = 0
      ef0 = efermi
      call info0(41,1,0,' Seek potential shift for fixed-spin mom ...')
      iter = 0

C --- Re-entry point for new guess at potential shift ---
C     iterate to itmax times.  iter<0 => last iteration
   10 continue
      iter = iter + 1

C ... Magnetic moment
      call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
      if (ipr .ge. 41) call awrit5(' dv=%;10,6D  yields '//
     .'ef=%;10,6D  amom=%;10,6D;  seeking %;10,6D',' ',160,
     .stdo,vnow,efermi,amom,fmom,vnow)
      call dvdos(vnow,amom,dosef,vhold,fmom,dvcap,dv)
      if (abs(dv) .lt. 1d-6) then
C       A root was found
        if (vhold(12) .eq. -2 .or. vhold(12) .eq. -3 .or.
     .  vhold(12) .eq.  0 .or. vhold(12) .eq.  1) then
          if (ipr .ge. 10)
     .    call awrit5('%N BZWTSF: potential shift bracketed.'//
     .    '  Unconstrained efermi=%,6;6d.'//
     .    '%N constraint fmom=%,6;6d  actual mmom=%,6;6d'//
     .    '  ef=%,6;6d  dv=%,6;6d',' ',160,stdo,ef0,
     .    fmom,amom,efermi,vnow)
          iter = -iter
        endif
      else if (iter .eq. itmax) then
        if (ipr .ge. 10)
     .  call awrit5('%N BZWTSF: failed to converge potential shift'//
     .  ' after %i iterations.'//
     .  '%N constraint fmom=%,6;6d  actual amom=%,6;6d'//
     .  '  ef=%,6;6d  dv=%,6;6d',' ',160,stdo,iter,fmom,amom,efermi,
     .  vnow)
        iter = -iter
      endif

C ... Potential shift
      do  ikp = 1, nkp
        do  ib = 1, nevx
          if (nspc .eq. 2) then
            eb(ib,1,ikp) = eb(ib,1,ikp) + vnow/2*swtk(ib,1,ikp)
            eb(ib,2,ikp) = eb(ib,2,ikp) + vnow/2*swtk(ib,2,ikp)
          else
            eb(ib,1,ikp) = eb(ib,1,ikp) - vnow/2
            eb(ib,2,ikp) = eb(ib,2,ikp) + vnow/2
          endif
        enddo
      enddo

C ... Fermi level with dv shift
      if (iter .gt. 0) call pshpr(ipr-50)
      call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,efermi,
     .sumev,wtkb,dosef,qval,ent)
      if (iprint() .ge. 20) then
        call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,
     .  amom)
        write(stdo,922) amom
      endif
      if (iter .gt. 0) call poppr

      do  ikp = 1, nkp
        do  ib = 1, nevx
          if (nspc .eq. 2) then
            eb(ib,1,ikp) = eb(ib,1,ikp) - vnow/2*swtk(ib,1,ikp)
            eb(ib,2,ikp) = eb(ib,2,ikp) - vnow/2*swtk(ib,2,ikp)
          else
            eb(ib,1,ikp) = eb(ib,1,ikp) + vnow/2
            eb(ib,2,ikp) = eb(ib,2,ikp) - vnow/2
          endif
        enddo
      enddo

      if (iter .gt. 0) goto 10
      if (lswtk.eq.1 .and. lwtkb .eq. 1) then
        lwtkb = 2
      elseif (lswtk.eq.1 .and. lwtkb .eq. 2) then
        lwtkb = 1
      endif

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

