subroutine pptrns(iopt,nl,ipc,nclass,nsp,alpha,nbas,pp,oold)
  !- Transform the set of potential parameters into another repsn
  ! ----------------------------------------------------------------
  !i Inputs
  !i   iopt: 1s digit
  !i         0: to passed alpha
  !i         1: to alpha = gamma.
  !i         2: to alpha=0
  !i         3: to alpha = gamma(spin1)
  !i         4: to alpha = gamma(spin2)
  !i         5: to alpha = (gamma(spin1) + gamma(spin2)/2
  !i         10s digit
  !i         1: set p(gamma) to zero
  !i         NB: sign of iopt<0 used to flag returning new alpha in alpha
  !i   nl,nclass,nsp,eny
  !i   alpha,ipc,nbas (needed only for iopt=0)
  !i   alph,c,sqdel,gam,palph (old)
  !o Outputs
  !o   pp are transformed to new representation (alph,c,sqdel,gam,palph)
  !o   oold -- overlap in old alpha representation
  !o   alpha   (iopt<0)
  !r Remarks
  !r   pp(1) : enu
  !r   pp(2) : calpha
  !r   pp(3) : srdel = sqrt(delta) but with proper sign (that of phi-).
  !r   pp(4) : palpha
  !r   pp(5) : gamma, or Q in Varenna notes
  !r   pp(6) : alpha, or qbar in Varenna notes
  !r  Transformations use the following (see Varenna, p88)
  !r    (1) srdnew/srdold = 1 + (cold-enu)*(alpnew-alpold)/srdold**2
  !r  where srdnew,srdold are sqrt(delta) for new and old representations
  !r    (2) 1 / o^a =  delta^a / (alpha - gamma) - (C^a - Enu)
  !r  (calculated by an external function);
  !r  change in palpha = change in oalpha**2 since
  !r    (3)  p = o**2 + p^gam
  !b Bugs
  !b   iopt=-1 ALWAYS returns alpha(spin 2), though alpha differs for
  !b   spins 1,2
  !u Updates
  !u   18 Apr 05 New 10's digit iopt
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: iopt,nl,nclass,nsp,nbas,ipc(1)
  double precision :: alpha(0:nl**2-1,nbas)
  double precision :: pp(6,nl,nsp,nclass),oold(nl,nsp,nclass)
  integer :: isp,jc,il,iclbas,jb,m,ib,nl2
  double precision :: xx,enu,gamma,pgam
  double precision :: cold,srdold,alpold,cnew,srdnew,alpnew,pold,pnew
  double precision :: oalpha
  external oalpha,iclbas

  do  10  jc = nclass, 1, -1
     jb = iclbas(jc,ipc)
     do  121  isp = 1, nsp
        do  12  il = 1, nl

           !        goto (1,2,3,4,5,6), mod(iabs(iopt),10)+1
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

           if (iopt < 0) then
              do  9  m = 1, 2*il-1
                 alpha((il-1)**2+m-1,jb) = alpnew
9             enddo
           endif

           ! --- Calculate potential parameters in new representation from old ---
           enu = pp(1,il,isp,jc)
           cold = pp(2,il,isp,jc)
           gamma = pp(5,il,isp,jc)
           alpold = pp(6,il,isp,jc)
           pold = pp(4,il,isp,jc)
           srdold = pp(3,il,isp,jc)

           !   ... delta=0 => no potential parameters for this l channel
           if (alpnew == alpold) goto 12
           if (srdold == 0) goto 12

           xx = 1 + (cold-enu)*(alpnew-alpold)/srdold**2
           srdnew = srdold*xx
           cnew = enu + (cold-enu)*xx

           oold(il,isp,jc) = oalpha(enu,cold,srdold**2,alpold,gamma)
           pgam = pold - oold(il,isp,jc)**2
           if (iabs(iopt) >= 10) pgam = 0
           !        pnew = pold - oold(il,isp,jc)**2 +
           !     .         oalpha(enu,cnew,srdnew**2,alpnew,gamma)**2
           pnew = pgam + oalpha(enu,cnew,srdnew**2,alpnew,gamma)**2

           pp(2,il,isp,jc) = cnew
           pp(3,il,isp,jc) = srdnew
           pp(6,il,isp,jc) = alpnew
           pp(4,il,isp,jc) = pnew
12      enddo
121  enddo
10 enddo

  ! --- If alpha is returned, copy alpha to all ib ---
  if (iopt < 0) then
     nl2 = nl*nl
     do  20  ib = 1, nbas
        jb = iclbas(ipc(ib),ipc)
        call dpscop(alpha,alpha,nl2,nl2*(jb-1)+1,nl2*(ib-1)+1,1d0)
20   enddo
  endif

end subroutine pptrns

real(8) function oalpha(enu,c,delta,alpha,gamma)
  !- Calculate overlap in alpha representation from pp's
  ! ----------------------------------------------------------------
  !i Inputs
  !i   enu,c,delta,alpha,gamma
  !o Outputs
  !o   oalpha
  !r Remarks
  !r   Varenna p 88, Eq 91 has:
  !r   1 / o^a =  (C^gam - Enu) - delta^gam / (gamma - alpha)
  !r   More generally, it can be written:
  !r   1 / o^a =  delta^a / (alpha - gamma) - (C^a - Enu)
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  double precision :: enu,c,delta,alpha,gamma
  ! Local parameters
  double precision :: xx
  xx = (alpha-gamma)/delta
  oalpha = xx/(1 - xx*(c-enu))
END function oalpha

