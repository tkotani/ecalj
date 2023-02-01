! FCPP#define F90 1
subroutine rnatm(pl,ql,n0,irchan,lmax,z,a,b,rofi,ev, &
     nr,rc,nsp,v,rho,plplus,qlplus)
  use m_lmfinit,only: stdo
  use m_ftox

  !- Renormalise charge or potential for a free atom
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   pl    :boundary conditions.  If Dl = log. deriv. at rmax,
  !i         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
  !i   ql    :valence charges for each l channel
  !i   n0    :dimensioning parameter
  !i   irchan:irchan(l+1) => suppress renormalization of that l
  !i   lmax  :maximum l for a given site
  !i   z     :nuclear charge
  !i   a     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   rofi  :radial mesh points
  !i   ev    :valence eigenvalues
  !i   nr    :number of radial mesh points
  !i   rc    :Fermi radius for renormalization cutoff
  !i         :rc>0 => charge renormalization
  !i         :rc<0 => potential renormalization
  !i         :rc(2) is (optional) width.  If zero, choose rc/8.
  !i   nsp   :2 for spin-polarized case, otherwise 1
  ! o Inputs/Outputs
  ! o  v     :On input, spherical potential (atomsr.f)
  ! o        :if rc<0, renormalized on output  NOT IMPLEMENTED
  !o Outputs
  !o   rho   :renormalized on output
  !l Local variables
  !l   ltop  :max l for which charge is nonzero
  !r Remarks
  !r
  !u Updates
  !u   01 Feb 06 Adapted from mol/rnatm.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nr,nsp,n0,lmax,irchan(1+lmax)
  double precision :: z,a,b,rc(2),pl(n0,2),ql(n0,2),qqq
  double precision :: rho(nr,nsp),v(nr,nsp),rofi(nr),ev(25)
  ! ... Local parameters
  logical :: pot
  integer :: l,isp,lp1,konfig,ir,iprint,ltop
  double precision :: sum1(10,2),sum2(10,2),tol,r,fac,sum,rmax,rcw
  double precision :: rh1(nr),rh2(nr),rwgt(nr),ddot
  external dpcopy,info2,makrvl,radwgt,rseq,rx
  integer:: iz
  real(8):: plplus(0:lmax,nsp),qlplus(0:lmax,nsp)

  !     call prrmsh('input rho',rofi,rho,nr,nr,nsp)

  if (rc(1) == 0) return
  ! FCPP#if ! (F90 | AUTO_ARRAY)
  ! FCPP      if (nr .gt. nrmax) call rx('RNATM: nr > nrmax')
  ! FCPP#endif

  ! --- Setup ---
  ! angenglob      stdo = nglob('stdo')
  !      stdo = globalvariables%stdo
  rmax = rofi(nr)
  call radwgt(rmax,a,nr,rwgt)
  rcw = rc(2)
  if (rcw == 0) rcw = abs(rc(1)/8)
  tol = 1d-10
  if (lmax > 8) call rx('RNATM: lmax > 8')
  if(iprint()>20) write(stdo,ftox)&
       ' RNATM: renormalize sphere density rc=',ftof(rc),'rcw=',ftof(rcw)
  pot = rc(1) .lt. 0
  if (pot) call rx('RNATM not set up for pot now')

  ltop = -1
  do   isp = 1, nsp
     if (nsp == 2) write(stdo,ftox)' Spin=',isp
     do l = 0, lmax
        do iz= 0,1
           lp1 = l+1
           !! bug fix around here jun2012
           !! but not tested yet since this routine is used only when 'SPEC_ATOM_RCFA'.
           if(iz==1) qqq= ql(lp1,isp)
           if(iz==0) qqq= qlplus(l,isp)
           if (qqq < 1d-6) goto 20
           if(iz==1) konfig = pl(lp1,isp)
           if(iz==0) konfig = plplus(l,isp)
           sum1(lp1,isp) = 0d0
           sum2(lp1,isp) = 0d0
           call makrvl(z,l,a,b,nr,rofi,konfig,v(1,isp),ev,tol,rh1)
           !       call intrho(rh1,rofi,a,b,nr,sum1(lp1,isp))
           sum1(lp1,isp) = qqq * ddot(nr,rwgt,1,rh1,1)
           if (irchan(lp1) == 0) then
              do  ir = 1, nr
                 r = rofi(ir)
                 fac = 1d0/(dexp((r-dabs(rc(1)))/rcw) + 1)
                 rh2(ir) = rh1(ir)*fac
              enddo
           else
              call dpcopy(rh1,rh2,1,nr,1d0)
           endif
           !       call intrho(rh2,rofi,a,b,nr,sum2(lp1,isp))
           sum2(lp1,isp) = qqq* ddot(nr,rwgt,1,rh2,1)
           fac = sum1(lp1,isp)/sum2(lp1,isp)
           do  ir = 1, nr
              rho(ir,isp) = rho(ir,isp) + qqq* (fac*rh2(ir)-rh1(ir))
           enddo
           ltop = max(ltop,l)
20         continue
        enddo
     enddo

     ! --- Printout ---
     if (iprint() >= 20) then
        sum = ddot(nr,rwgt,1,rho(1,isp),1)
        write(stdo,333) (sum1(lp1,isp), lp1=1,ltop+1)
        write(stdo,334) (sum2(lp1,isp), lp1=1,ltop+1)
        write(stdo,335) sum
333     format(' Ql before scaling:',8f9.5)
334     format(' Ql after  scaling:',8f9.5)
335     format(' valence q after renormalisation:',f9.5)
     endif
  enddo

end subroutine rnatm
!      subroutine intrho(rho,rofi,a,b,nr,sum)
!      implicit none
!      integer ir,nr
!      double precision rho(1),rofi(1),a,b,sum,r
!      sum = 0d0
!      do  1  ir = 2, nr
!        r = rofi(ir)
!        sum = sum + (mod(ir+1,2)+1)*(r+b)*rho(ir)
!    1 continue
!      r = rofi(nr)
!      sum = sum + .5d0*(r+b)*rho(nr)
!      sum = 2d0*sum*a/3d0
!      end

subroutine makrvl(z,l,a,b,nr,rofi,konfig,v,ev,tol,rho)
  use m_lmfinit,only:lrel_g=>lrel


  !- Makes contribution to charge density from a given spin and l,
  !  for free atom and given spherical potential (Adapted from nwrofp)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   z     :complex energy
  !i   z     :nuclear charge
  !i   l     :charge for particular l
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !i   konfig:core configuration
  !i   v     :spherical potential (atomsr.f)
  !i   ev
  !i   tol
  !o Outputs
  !o   rho  :partial density for given l
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   01 Feb 06
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ir,ival,konfig,l,nn,nr,nre,jr
  double precision :: v(nr),rofi(nr),rho(*),ev(*)
  double precision :: a,b,z,tol
  ! ... Local parameters
  integer :: lrel, nrx
  parameter( nrx=1501 )
  double precision :: eb1,eb2,eval,slo,sum,val,gfac,fllp1,c,tmc,r
  double precision :: g(2*nrx)
  eb1 = -50.d0
  eb2 = 15d0
  ival = 0
  nn = konfig-(l+1)
  ival = ival+1
  eval = ev(ival)
  val = 1.d-30
  slo = -val
  lrel = lrel_g
  if (lrel == 0) then
     call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v,g, &
          sum,a,b,rofi,nr,nre)
     do  77  ir = 1, nre
        rho(ir)= a*(rofi(ir)+b)*g(ir)**2
77   enddo
  else
     call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v,g, &
          sum,a,b,rofi,nr,nre)
     fllp1 = l*(l+1)
     c = 274.074d0
     rho(1) = 0
     do  78  ir = 2, nre
        jr = ir+nr
        r = rofi(ir)
        tmc = c - (v(ir) - 2*z/r - eval)/c
        gfac = 1 + fllp1/(tmc*r)**2
        rho(ir) = gfac*g(ir)**2 + g(jr)**2
78   enddo
  endif
  do  79  ir = nre+1, nr
     rho(ir)= 0d0
79 enddo
  !      print *, 'makrvl: l=',l
  !      call prrmsh('psi**2 in makrvl',rofi,rho,nre,nre,1)
end subroutine makrvl

