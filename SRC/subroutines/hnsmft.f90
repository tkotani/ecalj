subroutine hnsmft(rofi,rho,nr,qout,a,b,nrmt,e,qcst,xi,f,nf,err)! Fit the charge density tail to nf functions
  use m_lgunit,only:stdo
  !i Inputs
  !i   rofi  :radial mesh points
  !i   rho   :spherical valence charge density times 4*pi*r*r
  !i   nr    :number of radial mesh points
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nrmt  :number of points between 0..rmt
  !i   e     :vector of fitting function energies (printout only)
  !i   qcst  :vector of integrals of fitting functions from rmt...infty
  !i   xi    :vector of functions
  !i   nf    :number of functions to fit
  !o  Outputs
  !o   qout  :true charge in tail outside rmt
  !o   f     :fit coefficients (see remarks)
  !o   err   :rms error in fit
  !r  Remarks
  !r   The input charge density is defined on the mesh of points
  !r   1..nr, though only points nrmt..nr are included in the fitting.
  !r
  !r   rho is fit to the form rho(r) = sum_i=1,nf f_i H_i
  !r   where H_i is input as a function tabulated on the rofi mesh
  !r
  !r   Fits to the constraint of total charge conservation.
  !-----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: nr,nrmt,nf
  double precision :: rofi(nr),rho(nr),e(nf),f(nf+1),xi(nr,nf),qcst(nf), &
       qout
  ! Local variables
  integer :: i,ip,ipr,ir,ire,irs,is,j,nfmax,info,lpr
  parameter (nfmax = 10)
  integer :: ipiv(nfmax)
  double precision :: r,rhot,rhotol,err,a,b,dif,fpi,fit,qtf,qtr,qof,qor
  double precision :: s(nfmax,nfmax),wk(2*nfmax),wt
  real(8):: ff(nf)
  ! --- Setup ---
  if (nf+1 > nfmax) call rx('HNSMFT: nf gt nfmax')
  fpi = 16*datan(1d0)
  rhotol = 1d-12
  call getpr(ipr)
  !      stdo = lgunit(1)

  ! --- Find odd mesh points for integrations of the tail density ---
  irs = nrmt
  ire = irs + 3
  do  1  i = irs+5, nr-2, 2
     ire = i
     if (dabs(rho(i)) < rhotol) goto 2
1 enddo
2 continue

  ! --- Calculate charge outside the MT sphere; Simpson's rule ---
  qout = 0d0
  do  3  ir = irs+1, ire-1
     r = rofi(ir)
     qout = qout + (mod(ir+1,2)+1)*rho(ir)*(r+b)
3 enddo
  r = rofi(irs)
  qout = qout + .5d0*rho(irs)*(r + b)
  r = rofi(ire)
  qout = qout + .5d0*rho(ire)*(r + b)
  qout = 2d0*qout*a/3d0

  ! --- Calculate F-integrals (rho,H_i) and S-integrals (H_i,H_j) ---
  is = 0
  call dpzero(s,nfmax**2)
  do  5  i = 1, nf
     call fint(xi(nrmt,i),nrmt,rho,a,b,rofi,irs,ire,f(i))
     !       call awrit2('hnsmft f(%i) = %g',' ',80,stdo,i,f(i))
     do  4  j = 1, i
        is = is + 1
        call sint(xi(nrmt,i),xi(nrmt,j),nrmt, &
             a,b,rofi,irs,ire,s(i,j))
        s(j,i) = s(i,j)
        !         call awrit3('hnsmft s(%i,%i) = %g',' ',80,stdo,i,j,s(is))
4    enddo
5 enddo

  ! --- Add constraint that charge outsite rmt matches rho ---
  call dcopy(nf,qcst,1,s(1,nf+1),1)
  call dcopy(nf,qcst,1,s(nf+1,1),nfmax)
  s(nf+1,nf+1) = 0
  f(nf+1) = qout
  ! ccccccccccccc
  !      print *,'qqqq qcst=',qcst
  ! ccccccccccccc
  ! terative fitting.


  ! --- Least-squares fit ---
  !     call prmx('s',s,nfmax,nf+1,nf+1)
  !     call prmx('f',f,nfmax,nf,1)

  call dsytrf('L',nf+1,s,nfmax,ipiv,wk,2*nfmax,info)
  if (info /= 0) call rx('hnsmft: normal matrix not pos def')
  call dsytrs('L',nf+1,1,s,nfmax,ipiv,f,nf+1,info)

  ! akao simple inversion
  !      call dsytrf('L',nf,s,nfmax,ipiv,wk,2*nfmax,info)
  !      if (info .ne. 0) call rx('hnsmft: normal matrix not pos def')
  !      call dsytrs('L',nf,1,s,nfmax,ipiv,f,nf,info)

  ! akao exi(nf) first,...
  !      ff=f(1:nf)
  !      f(nf)=ff(nf)/s(nf,nf)
  !      do i=nf-1,1
  !        f(i)=(ff(i)-sum(s(i,i+1:nf)*ff(1:i+1)))/s(i,i)
  !      enddo


  ! --- Printout ---
  if (ipr > 30) then
     write(stdo,*)" ---E:energies of smHankels."// &
          " C:fitting coeeficient for core tail. ---"
     !        call awrit2(' E:%n:4;8F',' ',80,stdo,nf,e)
     !        call awrit2(' C:%n:4;8F',' ',80,stdo,nf,f)
     write(stdo,"(' E:',10f12.5)")(e(ip), ip=1,nf)
     write(stdo,"(' C:',10f12.5)")(f(ip), ip=1,nf)
  endif
  !     print 100, ire-irs+1,rofi(irs),rofi(ire),qout
  !     print 101, 'E',(e(ip), ip=1,nf)
  !     print 101, 'C',(f(ip), ip=1,nf)
  !  100 format(' HNSMFT:',i4,' points in interval',f9.5,f10.5,
  !     .  ';  q=',f10.6)
  ! 101 format(5x,a1,':',10f12.5)

  err = 0
  qof  = 0
  qor  = 0
  qtf  = 0
  qtr  = 0
  if (ipr > 30) write(stdo,200)
200 format(7x,' r ',9x,'rho',9x,'fit',9x,'diff')
  do  7  ir = 2, ire
     r = rofi(ir)
     rhot = rho(ir)/fpi/r/r
     fit = 0
     do  8  ip = 1, nf
        fit = fit + f(ip)*xi(ir,ip)
8    enddo
     fit = fit/dsqrt(fpi)
     dif = rhot - fit
     wt  = (mod(ir+1,2)+1)*fpi*r*r*(r+b)*2*a/3d0
     qtf = qtf + wt*fit
     qtr = qtr + wt*rhot
     if (ir >= irs) then
        if (ir == irs .OR. ir == ire) wt = wt*.5d0/(mod(ir+1,2)+1)
        err = err + wt*dif**2
        qof = qof + wt*fit
        qor = qor + wt*rhot
     endif
     !       rho(ir) = fpi*r*r*fit
     if (ipr <= 30 .OR. ir < irs .AND. ipr < 50) goto 7
     lpr = 0
     if (ipr >= 80 .OR. ir == irs) lpr = 1
     if (ipr >= 40 .AND. ir == 2) lpr = 1
     if (dabs(rhot) > 1d-6  .AND. (mod(ir-irs,10) == 0)) lpr = 1
     if (lpr >= 1) write(stdo,300) r,rhot,fit,dif
7 enddo
  err = err/(rofi(ire)-rofi(irs))
  err = dsqrt(err)
300 format(4f12.6)
  if (ipr >= 30) then
     write(stdo,400) qof,err
400  format('    q(fit):',f13.6,4x,'rms diff:',f11.6)
     write (stdo,500) qof,qtf-qof,qtf,qor,qtr-qor,qtr
500  format(4x,'fit: r>rmt',f10.6,'   r<rmt',f10.6,'   qtot',f10.6, &
          /4x,'rho: r>rmt',f10.6,'   r<rmt',f10.6,'   qtot',f10.6)
  endif
end subroutine hnsmft
subroutine sint(xi,xj,nrmt,a,b,rofi,irs,ire,sum)
  implicit none
  integer :: ire,irs,ir,nrmt,jr
  double precision :: rofi(1),r,a,b,sum,xi(1),xj(1)
  sum = 0d0
  do  1  ir = irs+1, ire-1
     jr = ir+1-nrmt
     r = rofi(ir)
     sum = sum + (mod(ir+1,2)+1)*(r+b)*r**2*xi(jr)*xj(jr)
1 enddo
  r = rofi(irs)
  sum = sum + .5d0*(r+b)*r**2*xi(irs+1-nrmt)*xj(irs+1-nrmt)
  r = rofi(ire)
  sum = sum + .5d0*(r+b)*r**2*xi(ire+1-nrmt)*xj(ire+1-nrmt)
  sum = 2d0*sum*a/3d0
end subroutine sint
subroutine fint(xi,nrmt,rho,a,b,rofi,irs,ire,sum)
  !     implicit none
  integer :: ire,irs,ir,nrmt
  double precision :: rho(1),rofi(1),r,a,b,sqfpi,sum,xi(1)
  sum = 0d0
  sqfpi = dsqrt(16*datan(1d0))
  do  1  ir = irs+1, ire-1
     r = rofi(ir)
     sum = sum + (mod(ir+1,2)+1)*rho(ir)*(r+b)*xi(ir+1-nrmt)
1 enddo
  r = rofi(irs)
  sum = sum + .5d0*rho(irs)*(r+b)*xi(irs+1-nrmt)
  r = rofi(ire)
  sum = sum + .5d0*rho(ire)*(r+b)*xi(ire+1-nrmt)
  sum = 2d0*sum*a/3d0/sqfpi
end subroutine fint

