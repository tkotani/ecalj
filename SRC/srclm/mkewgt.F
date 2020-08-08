      subroutine mkewgt(lmet,wgt,qval,ndimh,evl,ef0,def,esmear,
     .numq,nevec,ewgt,sumev,sumqv)
C- State-dependent weights for sampling BZ integration.
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmet  :0 assume insulator
Ci         :1 save eigenvectors to disk
Ci         :2 read weights from file, if they are needed
Ci         :3 always make two band passes; weights never needed a priori
Ci         :4 BZ integration with 3-point scheme
Ci   wgt   :symmetry weight for this qp
Ci   qval  :valence charge; number of occupied states
Ci   ndimh :dimensions evl
Ci   evl   :eigenvalues
Ci   ef0   :trial Fermi level
Ci   def   :uncertainty in Fermi level; quantities are accumulated for
Ci         :ef0-def, ef0, ef0+def
Ci   esmear:Parameter that describes gaussian broadening.
Ci         :Integer part >0 for for generalized gaussian broadening
Ci         :and is the the Methfessel-Paxton integration order
Ci         :Fractional part is the broadening width.
Ci         :Integer part <0 => Fermi-Dirac broadening used
Ci         :Fractional part is the temperature
Ci         :(see delstp.f)
Ci         :integer part above 100's digit is stripped.
Ci   numq  :number of trial Fermi levels
Co Outputs
Co   nevec :number of states with nonzero weights
Co   ewgt  :weights
Co   sumev :sumev(*,1..numq) are sumev(1..2) for numq Fermi levels
Co         :numq = 1 or 3, in which case Ef = ef0-def, ef0, ef0+def
Co         :sumev(1,1..numq) = sum of eigenvalues
Co         :sumev(2,1..numq) = entropy
Co   sumqv :sum of charges for numq Fermi levels
Cr Remarks
Cr   Weights are made for generalized Gaussian broadening, where a
Cr   delta function is expressed as a hermite polynomial*gaussian,
Cr   or for Fermi-Dirac broadening (see delstp.f)
Cu Updates
Cu   17 Jan 05 Extension of esmear to Fermi distribution
Cu   31 May 00 Adapted from nfp mkewgt
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmet,ndimh,nevec,numq
      double precision def,ef0,esmear,qval,wgt
      double precision evl(ndimh),ewgt(numq,1),sumev(2,3),sumqv(3)
C ... Local parameters
      integer stdo,lgunit,i,i1,i2,ie,iprint,iq,k,nord
      double precision dn,eff,fevec,fn,sn,width,wtop,wx,x

      stdo = lgunit(1)
      if (qval/2 .gt. ndimh) call fexit2(-1,111,'%N Exit -1 : MKEWGT: '
     .//'Basis with ndimh=%i is insufficient to carry q/2=%d states',
     .ndimh,qval/2)

C --- Nonmetal case ---
      if (lmet .eq. 0) then
        if (numq .ne. 1) call rxi('mkewgt: nonmetal but numq=',numq)
        fevec = qval/2d0
        nevec = fevec + 0.500001d0
        wtop = fevec-(nevec-1)
        if (dabs(wtop-1d0) .gt. 1d-6) call awrit1(' mkewgt (warning):'//
     .  ' uneven occupation for nonmet, top wgt=%d',' ',80,stdo,wtop)
        do  i = 1, nevec
          ewgt(1,i) = 1d0
          if (i .eq. nevec) ewgt(1,i) = wtop
          sumev(1,1) = sumev(1,1) + wgt*ewgt(1,i)*evl(i)
          sumev(2,1) = 0d0
          sumqv(1) = sumqv(1) + wgt*ewgt(1,i)
        enddo

C --- Metal case: require three fermi levels ---
      else
        if (numq .ne. 3) call rxi('mkewgt: metal but numq=',numq)
        width = dabs(esmear) - int(dabs(esmear))
        nord = dsign(1d0,esmear) * mod(int(dabs(esmear)),100)
        nevec = 0
        do  ie = 1, ndimh
          wx = 0d0
          do  k = 1, 3
            eff = ef0 + (k-2)*def
            x = (evl(ie)-eff)/width
            call delstp(nord,x,dn,fn,sn)
            ewgt(k,ie) = fn
            sumqv(k) = sumqv(k) + wgt*ewgt(k,ie)
            sumev(1,k) = sumev(1,k) + wgt*ewgt(k,ie)*evl(ie)
            sumev(2,k) = sumev(2,k) - wgt*width*sn
            wx = dmax1(wx,dabs(ewgt(k,ie)))
          enddo
          if (wx .lt. 1d-8) goto 80
C          write(stdo,268) ie,evl(ie),wx,(ewgt(k,ie),k=1,numq)
C  268     format(i5,f10.6,d13.2,8f12.6)
          nevec = ie
        enddo
  80    continue
      endif

!! this print out is confusing because we make truncation of sigm now in fpgw/exec/hqpe_sc
c$$$C --- Printout weights near frontier ---
c$$$      if (iprint() .gt. 40) then
c$$$C       write(stdo,200) lmet,numq,nevec
c$$$C 200   format(' mkewgt: lmet=',i2,'   numq=',i2,'  nevec=',i4)
c$$$        i2 = nevec
c$$$        i1 = max0(nevec-5,1)
c$$$        write(stdo,501) (i,i=i1,i2)
c$$$        write(stdo,502) (evl(i),i=i1,i2)
c$$$        do  iq = 1, numq
c$$$          write(stdo,503) iq,(ewgt(iq,i),i=i1,i2)
c$$$        enddo
c$$$  501   format(' state:',7i11)
c$$$  502   format(' evl:  ',7f11.6)
c$$$  503   format(' w',i1,':   ',7f11.6)
C        do  ie = 1, nevec
C          write(stdo,210) ie,evl(ie),(ewgt(iq,ie),iq=1,numq)
C  210     format('  state',i3,'   evl',f10.6,'   wgts:',3f12.6)
C        enddo
c$$$      endif
      end

