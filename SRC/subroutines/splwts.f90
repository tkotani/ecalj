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
