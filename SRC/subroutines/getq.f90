module m_getqvc ! Gets the charge and related parameters within an atom. Returns principal quantum numbers from pnu, estimating those unknown
  public getqvc,config
  private
contains
  subroutine getqvc(nsp,nl,lmx,z,pnu,qnu,ncmx,nvmx,kcor,lcor,qcor, qc,qt,dq,ec,ev)! Gets the charge and related parameters within an atom
    use m_ftox
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nl    :(global maximum l) + 1
    !i   lmx   :lmx(j) = maximum l for atom j
    !i   ncmx: 0 or maximum allowed number of core orbitals (see remarks)
    !i   nvmx: 0 or maximum allowed number of valence orbitals (see remarks
    !i   kcor  :(partial core occupation) p.q.n for occupation
    !i   lcor  :(partial core occupation) l quantum for occupation
    !i   qcor  :(partial core occupation) core charge and moment
    !o Outputs
    !o   qc:    core electronic charge
    !o   qt:    total charge (nuclear + electronic) within sphere
    !o   dq:    difference between spin up and spin down charge
    !o   ec:    starting guess for core eigenvalues (see remarks)
    !o   ev:    starting guess for valence eigenvalues (see remarks)
    !r Remarks
    !r   Generates charges and (optionally) initializes ec and ev
    !r   ec is initialized if ncmx > 0 and is otherwise unused
    !r   ev is initialized if nvmx > 0 and is otherwise unused
    !u Updates
    !u   26 Jan 06 bug fix for core hole, case lmx < lcor
    !r   (ATP) adapted from getq for one sphere, enables core hole
    ! ----------------------------------------------------------------
    implicit none
    integer :: nl,nsp,lmx,ncmx,nvmx,kcor,lcor
    double precision :: pnu(nl,nsp),qnu(3,nl,nsp),qcor(2), qc,qt,dq,z
    real(8),optional:: ec(*),ev(*)
    integer :: l,k,konf,konfig,konfg(0:10),isp,lmaxc,ncore,nval
    integer :: iprint
    double precision :: ecore0,eval0,deg
    character(8):: xt
    parameter (ecore0=-5.d0, eval0=-.5d0)
    call config(pnu(1,1),lmx,z,konfg,lmaxc)
    if (kcor > 0) then
       lmaxc = max(lmaxc,lcor)
       konfg(lcor) = max(konfg(lcor),kcor+1)
    endif
    qc = 0
    qt = 0
    dq = 0
    ncore = 0
    nval  = 0
    do  l = 0, lmaxc
       do  isp = 1, nsp
          deg = 2*(2*l+1)
          konfig = konfg(l)
          if (l <= lmx) konfig = pnu(l+1,isp)
          if (konfig <= 0) then
             !            if (iprint() .eq. 0) call setpr(10)
             call rx('GETQVC: bad l pnu='//trim(xt(l))//' '//ftof(pnu(l+1,isp)))
          endif
          do  konf = l+1, konfig-1
             if (ncmx /= 0) then
                ncore = ncore+1
                ec(ncore) = ecore0
             endif
             if (konf == kcor .AND. l == lcor) then
                deg = deg + qcor(1)
                if (nsp > 1) dq = dq + qcor(2)/nsp
             endif
             qc = qc + deg/nsp
          enddo
       enddo
    enddo
    do  k = 1, lmx+1
       if (nsp > 1) dq = dq + qnu(1,k,1) - qnu(1,k,2)
       do  isp = 1, nsp
          if (nvmx /= 0) then
             nval = nval+1
             ev(nval) = eval0
          endif
          qt = qt + qnu(1,k,isp)
       enddo
    enddo
    if (ncore > ncmx .OR. nval > nvmx) call rx('GETQVC: too many orbitals')
    qt = qc + qt - z
  end subroutine getqvc
  subroutine getq(nsp,nl,lmx,nc,z,pnu,qnu,ics,sspec,qc,qt,dq)  !- Gets the charge and related parameters for a set of atoms
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nl    :(global maximum l) + 1
    !i   lmx   :lmx(j) = maximum l for atom j
    !i   nc    :number of classes
    !i   z     :nuclear charge
    !i   pnu   :pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i         :Integer part used to get princ. quant. number for counting
    !i   qnu   :energy-weighted moments of the sphere charges
    !i   ics   :species table: class ic belongs to species ics(ic)
    !i         :Set to zero to suppress information for core hole
    !i   sspec :struct for species-specific information; see routine uspec
    !i         :Not used unless ics(1)>0
    !i     Elts read: kcor,lcor,qcor
    !o Outputs
    !o   qc:    core electronic charge
    !o   qt:    total charge (nuclear + electronic) within sphere
    !o   dq:    difference between spin up and spin down charge
    !u Updates
    !u   04 Jan 06 Redesign to enable core hole
    !u   17 Feb 03 Changed dq to correct convention (q+ - q-)
    ! ----------------------------------------------------------------
    implicit none
    integer :: nl,nsp,nc,lmx(nc),ics(nc)
    double precision :: pnu(nl,nsp,nc),qnu(3,nl,nsp,nc), qc(nc),qt(nc),dq(nc),z(nc),sspec(1)
    integer :: ic,kcor,lcor,is
    double precision :: qcor(2)
    do  ic = 1, nc
       if (ics(1) == 0) then
          kcor = 0
          lcor = 0
       else
          is = ics(ic)
          call gtpcor(is,kcor,lcor,qcor)
       endif
       call getqvc(nsp,nl,lmx(ic),z(ic),pnu(1,1,ic),qnu(1,1,1,ic),0,0, kcor,lcor,qcor,qc(ic),qt(ic),dq(ic))
    enddo
  end subroutine getq
  subroutine config(pnu,lmax,z,konfig,lmaxc) !- Returns principal quantum numbers from pnu, estimating those unknown
    !i Inputs
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   lmax  :maximum l for which pnu is supplied
    !i   z     :nuclear charge
    !o Outputs
    !o   konfig:estimated principal quant. no. of valence state for l>lmax
    !o         :Core orbitals are specified by:
    !o         :  1, 2, ..., konf(0)-1 for s
    !o         :  2, 3, ..., konf(1)-1 for p
    !o         :  3, 4, ..., konf(2)-1 for d, and so on.
    !o   lmaxc :largest l for which there is a core state
    !r Remarks
    !u Updates
    !u   08 Feb 01 generate config for l=0..8
    implicit none
    integer :: lmax,konfig(0:8),lmaxc
    double precision :: pnu(0:*),z
    integer :: l,LI,NA,K,RB,CS,FR,CU,AG,AU,HF
    parameter (LI=3, NA=11, K=19, RB=37, CS=55, FR=87, CU=29, AG=47, AU=79, HF=72)
    ! --- Calculate lmaxc ---
    lmaxc = lmax
    if (z >= LI) lmaxc = max(lmax,0)
    if (z >= NA) lmaxc = max(lmax,1)
    if (z >= CU) lmaxc = max(lmax,2)
    if (z >= HF) lmaxc = max(lmax,3)
    ! --- Estimate konfig ---
    do  10  l = 0, 8
       konfig(l) = l+1
10  enddo
    if (z >= LI) konfig(0) = 2
    if (z >= NA) konfig(0) = 3
    if (z >=  K) konfig(0) = 4
    if (z >= RB) konfig(0) = 5
    if (z >= CS) konfig(0) = 6
    if (z >= FR) konfig(0) = 7
    konfig(1) = max(konfig(0),2)
    if (z >= CU) konfig(2) = 4
    if (z >= AG) konfig(2) = 5
    if (z >= AU) konfig(2) = 6
    if (z >= HF) konfig(3) = 5
    ! --- Override konfig with given pnu ---
    do  20  l = 0, lmax
       konfig(l) = pnu(l)
20  enddo
  end subroutine config
end module m_getqvc

