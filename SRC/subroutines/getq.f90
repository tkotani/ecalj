
module m_getqvc
contains
  subroutine getqvc(nsp,nl,lmx,z,pnu,qnu,ncmx,nvmx, &
       kcor,lcor,qcor,qc,qt,dq,ec,ev)
    !- Gets the charge and related parameters within an atom
    ! ----------------------------------------------------------------
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
    !     implicit none
    ! ... Passed parameters
    integer :: nl,nsp,lmx,ncmx,nvmx,kcor,lcor
    double precision :: pnu(nl,nsp),qnu(3,nl,nsp),qcor(2), &
         qc,qt,dq,z
    real(8),optional:: ec(*),ev(*)
    ! ... Local parameters
    integer :: l,k,konf,konfig,konfg(0:10),isp,lmaxc,ncore,nval
    integer :: iprint
    double precision :: ecore0,eval0,deg
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
             call fexit2(-1,1,' Exit -1 GETQVC: '// &
                  'bad pnu(l=%i) (%,1d)',l,pnu(l+1,isp))
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

    if (ncore > ncmx .OR. nval > nvmx) &
         call rx('GETQVC: too many orbitals')

    qt = qc + qt - z
  end subroutine getqvc

end module m_getqvc

subroutine getq(nsp,nl,lmx,nc,z,pnu,qnu,ics,sspec,qc,qt,dq)
  !- Gets the charge and related parameters for a set of atoms
  ! ----------------------------------------------------------------
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
  use m_getqvc
  !     implicit none
  ! ... Passed parameters
  integer :: nl,nsp,nc,lmx(nc),ics(nc)
  double precision :: pnu(nl,nsp,nc),qnu(3,nl,nsp,nc), &
       qc(nc),qt(nc),dq(nc),z(nc),sspec(1)
  ! ... Local parameters
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
     call getqvc(nsp,nl,lmx(ic),z(ic),pnu(1,1,ic),qnu(1,1,1,ic),0,0, &
          kcor,lcor,qcor,qc(ic),qt(ic),dq(ic))
  enddo

end subroutine getq
