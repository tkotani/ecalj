subroutine splwts(nqp,nband,nbmx,nsp,wgts,evl,n,w,efermi, &
     metal,sumev,bndwts,wtot,entrpy,dosef,cv)
  use m_lmfinit,only: stdo
  !- make sampling weights for integrals under the Fermi surface
  !-----------------------------------------------------------------------
  !i  Input
  !i    nqp : number of q-points; nband : number of bands
  !i    wgts: band weights
  !i    evl : energy eigenvalues
  !i    n   : n>0 Methfessel-Paxton polynomial order
  !i        : n<0 sampling done with Fermi-Dirac statistics
  !i    w   : n>0 gaussian width in Methfessel-Paxton integration (Ry)
  !i        : n<0 Temperature for Fermi distribution (Ry)
  !i    nbmx : first dimension of evl ;
  !i    metal : if F, weights unity below E_f and zero above.
  !i    efermi : Fermi energy
  !o  Output
  !o    bndwts : band and E_F - dependent k-point weights for integration
  !o    wtot   : sum of all weights (charge)
  !o    entrpy : electron entropy
  !o    dosef  : DOS at Fermi energy
  !o    cv     : electronic specific heat
  !o             (only evaluated with Fermi-Dirac statistics)
  !r  Remarks
  !r    sum of occupied eigenvalues = sum_n,k  w_nk E_nk
  !r    w_nk are generalised occupation numbers;
  !r    see Needs et al. Phys Rev B 33 (1986) 3778, eqs 1 & 2.
  !u Updates
  !u   16 Jul 08 returns entropy as TS for all n
  !u   04 Aug 07 Generates dos(efermi), cv(T=w) for F-D statistics
  !u   02 May 07 (MvS) prints entropy to stdout
  !u   21 Jun 06 (ATP) generates entrpy as output
  !u   17 Jan 05 Output wtot
  !-----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nqp,nband,nbmx,nsp,n,ix,isplwts,i_copy_size
  logical :: metal
  double precision :: wgts(nqp),evl(nbmx,nsp,nqp),w,efermi,sumev, &
       bndwts(nband,nsp,nqp),wtot,entrpy,dosef,cv
  ! ... Local parameters
  integer :: iqp,iband,isp,iprint,i1mach
  double precision :: e,s,d,wt,x,xx,dsdt,tdsdt
  logical :: fractional
  sumev = 0d0
  wtot = 0d0
  entrpy = 0d0
  dosef = 0
  tdsdt = 0
  fractional=.false.
  do  3  iqp = 1, nqp
     do  2  iband = 1, nband
        do  1  isp = 1, nsp
           e = evl(iband,isp,iqp)
           if (metal) then

              !             debugging: check derivative numerically
              !              x = (efermi - e) / (w+1d-7)
              !              call delstp(n,x,d,s,sp)
              !              x = (efermi - e) / (w-1d-7)
              !              call delstp(n,x,d,s,sm)
              !              dsdt1 = (sp-sm)/2d-7

              x = (efermi - e) / w
              !             call delstp(n,x,d,s,xx)
              call delstd(n,x,d,s,xx,dsdt)
              if (abs(x) < 36) then
                 dsdt = -(efermi-e)/w**2 * dsdt
              else
                 dsdt = 0
              endif

              !C             debugging: compare analytical, numerical derivative
              !              if (abs(x) .lt. 30) then
              !                print 222, x,dsdt1,dsdt,dsdt-dsdt1
              !  222           format(3f14.8,1pe12.3)
              !              endif
           else
              s = 1d0
              if (e <= efermi) s = 0d0
              xx = 0
              d = 0
           endif
           wt = abs(wgts(iqp)) * (1d0 - s) / nsp
           bndwts(iband,isp,iqp) = wt
           if(0.1d0<wt .AND. wt<0.9d0) then
              fractional=.true.
           endif
           dosef = dosef + d*abs(wgts(iqp))/w/nsp
           wtot = wtot + wt
           sumev = sumev + e * wt
           entrpy = entrpy + xx  * abs(wgts(iqp)) / nsp
           tdsdt  = tdsdt + dsdt * abs(wgts(iqp)) / nsp
1       enddo
2    enddo
3 enddo
  tdsdt = tdsdt*w
  entrpy = entrpy*w
  if (n < 0) cv = tdsdt

  ! ... Print out band weights, if only 1 kp
  if (iprint() > 30 .AND. nqp == 1) then
     isplwts=1093
     open(isplwts,file="BandWeight.dat")
     do  isp = 1, nsp
        call info2(30,0,0,' SPLWTS: band weights .. '// &
             '%?;(n==2);Spin %i;;%N       eval      weight',nsp,isp)
        if(fractional) then
           write(isplwts,*)"! Fractional occupation (criterion 0.1<wgt<0.9)"
           write(stdo,*)   "! Fractional occupation (criterion 0.1<wgt<0.9)"
        endif
        ix=0
        do  iband = 1, nband
           if(bndwts(iband,isp,1)<1d-7) ix=ix+1
           if(ix==10) then
              write (stdo,"('     ... ')")
              exit
           endif
           write (stdo,20) iband,evl(iband,isp,1),bndwts(iband,isp,1)
           write (isplwts,20) iband,evl(iband,isp,1),bndwts(iband,isp,1)
20         format (4x,i5,2f10.6)
        enddo
     enddo
     close(isplwts)
  endif
  ! ... Print out various k-point integrations
  if (iprint() >= 10) then
     if (n >= 0) then
        !          call awrit6(' N=%i, W=%d, E_F=%d, sumev=%d, entropy term:'
        !     .    //' %d, %d electrons',' ',256,i1mach(2),
        !     .    n,w,efermi,sumev,entrpy,wtot)
        write(stdo,"(a,i5,5d13.5)")'N W E_F sumev TS Nele=', &
             n,w,efermi,sumev,entrpy,wtot
     else
        !          call awrit5(' T=%dK, E_F=%d, sumev=%d, TS=%;3g,'
        !     .          //' %d electrons',' ',256,i1mach(2),
        write(stdo,"(a,5d13.5)")'T(K) E_F sumev TS Nele=', &
             0.1579d6*w,efermi,sumev,entrpy,wtot
        write(stdo,"(a,2d13.5)")'Entropy S, specific heat TdS/dT=',entrpy/w,tdsdt
        write(stdo,"(a)")'Fermi-Dirac;sampling'
     endif
     !      call info5(10,0,0,' SPLWTS: Fermi energy:%;6d;'//
     !     .  '  band energy=%;6d;  %;6d electrons  DOS(E_f)=%;4g',
     !     .    efermi,sumev,wtot,dosef,0)
  endif
end subroutine splwts

subroutine delstd(n,x,d,s,e,ep)
  !- Returns generalised delta and step functions (Methfessel & Paxton)
  !-----------------------------------------------------------------------
  !i  Inputs
  !i    n  : order of approximant; see Remarks
  !i       : n>=0 returns Methfessel-Paxton broadening
  !i       : n<0  returns Fermi-Dirac broadening
  !i    x  : (efermi - e) / width
  !i       : width should be gaussian width (n>=0)
  !i       : or temperature (n<0)
  !o  Outputs
  !o    D_n (x): smeared delta-function
  !o    S_n (x): smeared heaviside function
  !o    e_n (x): entropy
  !o    ep     : de/dx (Fermi-Dirac case only)
  !r  Remarks
  !r    For Methfessel-Paxton (generalized gaussian) broadening
  !r    (see Phys Rev B40, 3616 (1989))
  !r      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
  !r      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
  !r      where H is a Hermite polynomial and
  !r      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
  !r    For Fermi-Dirac broadening
  !r      s = 1/(exp(x)+1)   (fermi function)
  !r      d = ds/dx = exp(x)*s*s
  !r      e = -( s*log(s) + (1-s)*log(1-s) )
  !r     ep = log(s) - log(1-s)
  !u Updates
  !u   04 Aug 07 extended delstp to returns ep in Fermi-dirac case
  !u   23 May 00 extended to handle Fermi-Dirac broadening
  !-----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n
  double precision :: x,d,s,e,ep
  ! ... Local parameters
  integer :: i,k
  double precision :: a,h1,h2,h3,s0,ex2,derfc,srpi
  !      intrinsic dsqrt,datan,dexp
  external derfc

  srpi = dsqrt(4d0*datan(1d0))

  ! ... Fermi-Dirac broadening
  if (n < 0) then
     if (x < -36d0) goto 91
     if (x >  36d0) goto 92
     s = 1d0/(dexp(x)+1d0)
     d = dexp(x)*s*s
     e = -( s*dlog(s) + (1-s)*dlog(1-s) )
     ep = (dlog(s) - dlog(1-s)) * s**2 * exp(x)
     return
  endif

  ! ... Methfessel-Paxton broadening
  if (x < -6d0) goto 91
  if (x >  6d0) goto 92
  ex2 = dexp(-x*x)
  s0 = .5d0 * derfc(x)
  a = 1d0/srpi
  k = 0
  h1 = 1d0
  h2 = 2d0 * x
  s = 0d0
  d = a
  do  1  i = 1, n
     a = -a / ( 4d0*i )
     k = k+1
     h3 = h1
     h1 = h2
     h2 = 2*x*h2 - 2*k*h3
     s = s + a*h1
     k = k+1
     h3 = h1
     h1 = h2
     h2 = 2*x*h2 - 2*k*h3
     d = d + a*h1
1 enddo
  d = d * ex2
  s = s0 + s*ex2
  e = 0.5d0*a*h1*ex2
  ep = 0
  return

  ! ... Branch for very small or very large x
91 s = 1d0
  e = 0d0
  d = 0d0
  ep = 0d0
  return

92 s = 0d0
  e = 0d0
  d = 0d0
  ep = 0d0
  return

end subroutine delstd

