      subroutine magtrq(nbas,nl,nclass,ipc,sdmod,sdprm,ssite,sspec,
     .  dclabl,pp,rhos,nrhos,ehf,eula,neul,frc,aamom)
C- Magnetic torque and update of the Euler angles
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   nclass:number of inequivalent classes
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   sdmod :specifies kind of spin statics or dynamics; see Remarks
Ci         :1s digit
Ci           0: Euler angles updated by zeroing out off-diagonal
Ci              parts of the spin-density matrix
Ci           1: Euler angles updated by force * sdprm(2)
Ci           2: Simple spin dynamics w/; see remarks
Ci           3: Spin dynamics with friction terms
Ci         :10's digit affects mixing of eula (not used here)
Ci         :100's digit
Ci            1: Nose thermostat (see Remarks)
Ci         :1000's digit
Ci            1: Reserved by lm to prevent updating of atomic pp's
Ci   sdprm :parameters governing dynamics (see Remarks)
Ci         :(1) scale factor amplifying magnetic forces
Ci         :(2) time step tau for Landau dynamics
Ci         :(3) reference energy etot0 for Nose thermostat
Ci              (not used unless fscal is zero)
Ci         :(4) maximum allowed change in angle (not used unless nonneg)
Ci         :(5) etol: set tau=0 this iter if etot-ehf>etol
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec class relax
Ci     Stored:    *
Ci     Passed to: *
Ci   dclabl:class name, packed as a real number
Ci   pp    :potential parameters (atomsr.f)
Ci   rhos  : spin density matrix
Ci         : Note: rhos should be hermitian in spin space, but may not
Ci         : be owing to energy integration errors in the complex
Ci         : plane. Use a symmetrized form to minimize errors.
Ci   neul  :1 if Euler angles are l- and m-independent,
Ci         :nl if Euler are l-dependent and m-independent
Ci         :nl**2 if Euler are l- and m-dependent
Ci   ehf
Co Outputs
Co   eula  :updated, if sdmod<>3
Co   frc   :if sdmod=2
Co   aamom :local magnetic moments
Cl Local variables
Cl   ifrlx :suppress updating any Euler angles for any ifrlx nonzero
Cl   ila   :1 if Euler angles are l-dependent, otherwise 0
Cr Remarks
Cr   Spin statics: (sdmod=1):  Direction vector e (=z in loc. coord)
Cr     updated by fscal * J, where J = "force"
Cr   Spin dynamics: (sdmod=2):  Direction vector rotated by
Cr     z = tau*gyro/M (z cross J) + Nose thermostat
Cr     tau, M, gyro:  time step, moment amplitude, gyromagnetic ratio.
Cr     Nose thermostat accounts for coupling to the lattice, and adds
Cr     'statics' rotation, whose magnitude magtrq estimates to change
Cr     the internal energy to the target value etot0.  It estimates
Cr     the rotation angle from dE = 1/2 sum f_j th_j, and sets each
Cr     th_j in proportion to the f_j.  The th_j estimated by the above
Cr     is additionally scaled by sdprm(1).
Cr   Formulas for moments:
Cr            M_x =  2 Re(rho+-) = (rho+-)+(rho-+)
Cr            M_y =  2 Im(rho+-) = ((rho+-)-(rho-+))/i
Cr            M_z =  (rho++)-(rho--)
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be quite hermitian; e.g. when rhos is generated
Cr   by a Green's function technique.
Cb Bugs
Cr   No attempt to symmetrize updated Euler angles by class.
Cu Updates
Cu   21 Apr 04 Additions for an m-dependent spin-density matrix
Cu   04 Sep 03 Fix improper weighting of spin density matrix for
Cu             sdmod=0
Cu   27 Jan 03 Some alterations : sdmod is now consistent with docs
Cu   22 May 02 Bug fix to handle high and neglected orbitals
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer nbas,nl,nrhos,neul,nclass,ipc(*),sdmod
      double precision pp(6,nl,2,nclass),ssite(1),sspec(1),
     .  eula(nbas,neul,3),frc(nbas,3),aamom(nbas),ehf,
     .  sdprm(5),rhos(2,0:2,nrhos,2,2,nclass),dclabl(nclass)

      call rx('magtrq implemented in nc package')
      end

