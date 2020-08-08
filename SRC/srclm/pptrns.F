      subroutine pptrns(iopt,nl,ipc,nclass,nsp,alpha,nbas,pp,oold)
C- Transform the set of potential parameters into another repsn
C ----------------------------------------------------------------
Ci Inputs
Ci   iopt: 1s digit
Ci         0: to passed alpha
Ci         1: to alpha = gamma.
Ci         2: to alpha=0
Ci         3: to alpha = gamma(spin1)
Ci         4: to alpha = gamma(spin2)
Ci         5: to alpha = (gamma(spin1) + gamma(spin2)/2
Ci         10s digit
Ci         1: set p(gamma) to zero
Ci         NB: sign of iopt<0 used to flag returning new alpha in alpha
Ci   nl,nclass,nsp,eny
Ci   alpha,ipc,nbas (needed only for iopt=0)
Ci   alph,c,sqdel,gam,palph (old)
Co Outputs
Co   pp are transformed to new representation (alph,c,sqdel,gam,palph)
Co   oold -- overlap in old alpha representation
Co   alpha   (iopt<0)
Cr Remarks
Cr   pp(1) : enu
Cr   pp(2) : calpha
Cr   pp(3) : srdel = sqrt(delta) but with proper sign (that of phi-).
Cr   pp(4) : palpha
Cr   pp(5) : gamma, or Q in Varenna notes
Cr   pp(6) : alpha, or qbar in Varenna notes
Cr  Transformations use the following (see Varenna, p88)
Cr    (1) srdnew/srdold = 1 + (cold-enu)*(alpnew-alpold)/srdold**2
Cr  where srdnew,srdold are sqrt(delta) for new and old representations
Cr    (2) 1 / o^a =  delta^a / (alpha - gamma) - (C^a - Enu)
Cr  (calculated by an external function);
Cr  change in palpha = change in oalpha**2 since
Cr    (3)  p = o**2 + p^gam
Cb Bugs
Cb   iopt=-1 ALWAYS returns alpha(spin 2), though alpha differs for
Cb   spins 1,2
Cu Updates
Cu   18 Apr 05 New 10's digit iopt
C ----------------------------------------------------------------
C     implicit none
      integer iopt,nl,nclass,nsp,nbas,ipc(1)
      double precision alpha(0:nl**2-1,nbas)
      double precision pp(6,nl,nsp,nclass),oold(nl,nsp,nclass)
      integer isp,jc,il,iclbas,jb,m,ib,nl2
      double precision xx,enu,gamma,pgam
      double precision cold,srdold,alpold,cnew,srdnew,alpnew,pold,pnew
      double precision oalpha
      external oalpha,iclbas

      do  10  jc = nclass, 1, -1
        jb = iclbas(jc,ipc)
        do  12  isp = 1, nsp
        do  12  il = 1, nl

c        goto (1,2,3,4,5,6), mod(iabs(iopt),10)+1
          select case(mod(iabs(iopt),10)+1)
          case(1)
          alpnew = alpha((il-1)**2,jb)
          case(2)
          alpnew = pp(5,il,isp,jc)
          case(3)
          alpnew = 0
          case(4)
          alpnew = pp(5,il,1,jc)
          case(5)
          alpnew = pp(5,il,nsp,jc)
          case(6)
          alpnew = (pp(5,il,1,jc) + pp(5,il,nsp,jc))/2
          end select

          if (iopt .lt. 0) then
            do  9  m = 1, 2*il-1
              alpha((il-1)**2+m-1,jb) = alpnew
    9       continue
          endif

C --- Calculate potential parameters in new representation from old ---
          enu = pp(1,il,isp,jc)
          cold = pp(2,il,isp,jc)
          gamma = pp(5,il,isp,jc)
          alpold = pp(6,il,isp,jc)
          pold = pp(4,il,isp,jc)
          srdold = pp(3,il,isp,jc)

C   ... delta=0 => no potential parameters for this l channel
          if (alpnew .eq. alpold) goto 12
          if (srdold .eq. 0) goto 12

          xx = 1 + (cold-enu)*(alpnew-alpold)/srdold**2
          srdnew = srdold*xx
          cnew = enu + (cold-enu)*xx

          oold(il,isp,jc) = oalpha(enu,cold,srdold**2,alpold,gamma)
          pgam = pold - oold(il,isp,jc)**2
          if (iabs(iopt) .ge. 10) pgam = 0
C        pnew = pold - oold(il,isp,jc)**2 +
C     .         oalpha(enu,cnew,srdnew**2,alpnew,gamma)**2
          pnew = pgam + oalpha(enu,cnew,srdnew**2,alpnew,gamma)**2

          pp(2,il,isp,jc) = cnew
          pp(3,il,isp,jc) = srdnew
          pp(6,il,isp,jc) = alpnew
          pp(4,il,isp,jc) = pnew
   12   continue
   10 continue

C --- If alpha is returned, copy alpha to all ib ---
      if (iopt .lt. 0) then
        nl2 = nl*nl
        do  20  ib = 1, nbas
          jb = iclbas(ipc(ib),ipc)
          call dpscop(alpha,alpha,nl2,nl2*(jb-1)+1,nl2*(ib-1)+1,1d0)
   20   continue
      endif

      end
      double precision function oalpha(enu,c,delta,alpha,gamma)
C- Calculate overlap in alpha representation from pp's
C ----------------------------------------------------------------
Ci Inputs
Ci   enu,c,delta,alpha,gamma
Co Outputs
Co   oalpha
Cr Remarks
Cr   Varenna p 88, Eq 91 has:
Cr   1 / o^a =  (C^gam - Enu) - delta^gam / (gamma - alpha)
Cr   More generally, it can be written:
Cr   1 / o^a =  delta^a / (alpha - gamma) - (C^a - Enu)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      double precision enu,c,delta,alpha,gamma
C Local parameters
      double precision xx

      xx = (alpha-gamma)/delta
      oalpha = xx/(1 - xx*(c-enu))

      end

