      subroutine scalsr(iax,ldot,ndimL,nl2,npr,trad,tral,sc,sdotc)
C- Scales real-space screened S and Sdot
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci         :nlma = iax(9,1) and for channel k, nlmb = iax(9,k)
Ci   ldot  :T: calculate energy derivatives
Ci   ndimL :leading dimension of sc and sdotc
Ci   nl2   :second dimension of tral,trad
Ci   npr   :number of neighbors around each basis atom
Ci   trad  :(kappa*avw)^2-derivative of tral
Ci   tral  :transformation matrix for head and tail functions (mktra2.f)
Cio Inputs/Outputs:
Cio  sc    :On input, unscaled structure constants
Cio        :sc = sc(ndimL,nlma)
Cio        :On output sc is overwritten (scaled) by:
Cio        :
Cio        : sc_i,k -> 1/t3_i  sc_i,k  1/t3_k det(tral_k) + t1/t3 delta_i,k
Cio        :
Cio  sdotc :On input, energy derivative of unscaled sc
Cio        :sdotc = sdotc(ndimL,nlma)
Cio        :On output sdotc is overwritten (scaled) by:
Cio        :
Cio        : sdotc_i,k -> -1/t3_i  sdotc_i,k  1/t3_k  det(tral_k)
Cio        :              -  sc_i,k  td3_k/t3_k
Cio        :              -  td3_i/t3_i sc_i,k
Cio        :              +  (td1 - td3 t1/t3)/t3 delta_i,k
Cio        :
Cr Remarks:
Cr  scalsr  was adapted from the Stuttgart LMTO scals written by R. Tank
Cr  To make identification with OKA, "Developing the MTO formalism",
Cr  note his equation 52
Cr    B^a = 1/j(ka) [B(E) - k cot(alpha)]^-1 1/j(ka) + a partial{j(ka)}
Cr  is equivalent to the present with the following substitutions
Cr          code                                 book
Cr       t3 = -2/avw*J(ka)                       j(ka)
Cr       alpha^-1 = J(ka)/N(ka)                  k cot(alpha^OKA)
Cr       t1/t3 = (2*a^2/w J'(ka))/(2/w J(ka))    a^2 j'(ka)/j(ka)
Cr             = a^2 J'(ka)/ J(ka)
Cr       ...
Cr       NB: (This code sticks to old definition of alpha
Cr            while book's alpha^OKA refers to phase shift)
Cr       NB: W{N,J} = w/2 -> J'*N-J*N' = w/2/a^2 ; whereas W{j,n}=1/k
Cr       ? should be making S^a, not B^a?
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      logical ldot
      integer npr,niax,nl2,ndimL
      parameter (niax=10)
      integer iax(niax,npr)
      double precision sc(ndimL,npr),sdotc(ndimL,npr),
     .                 tral(4,nl2,*),trad(4,nl2,*)
C Local variables:
      integer ii,ilm,ipr,klm,nlmb,ia,ib,nlma
      double precision dt,scli,sclk,sd,sdd

      ia = iax(1,1)
      nlma = iax(9,1)

      ii = 1
      do  30  ipr = 1, npr
        ib = iax(2,ipr)
        nlmb = iax(9,ipr)
        do  32  ilm = 1, nlmb
          scli = 1d0/tral(3,ilm,ib)
          do  36  klm = 1, nlma
            sclk = 1d0/tral(3,klm,ia)
            dt = tral(1,klm,ia)*tral(4,klm,ia) -
     .           tral(2,klm,ia)*tral(3,klm,ia)
            sc(ii,klm) = scli*sc(ii,klm)*sclk*dt
            if (ldot) then
              sdotc(ii,klm) = -scli*sdotc(ii,klm)*sclk*dt
              sdotc(ii,klm) = sdotc(ii,klm)
     .                        -sc(ii,klm)*trad(3,klm,ia)/tral(3,klm,ia)
     .                        -trad(3,ilm,ib)/tral(3,ilm,ib)*sc(ii,klm)

C              sdotc(ii,klm) = -scli*sdotc(ii,klm)*sclk*dt
C     .                        -sc(ii,klm)*trad(3,klm,ia)/tral(3,klm,ia)
C     .                        -trad(3,ilm,ib)/tral(3,ilm,ib)*sc(ii,klm)
C              print *, -scli*sdotc(ii,klm)*sclk*dt
C              print *, sdotc(1,1)
              if (klm .eq. ii) then
                sd = tral(1,ilm,ia)/tral(3,ilm,ia)
                sdd = (trad(1,ilm,ib) -sd*trad(3,ilm,ib))/tral(3,ilm,ib)
                sdotc(klm,klm) = sdotc(klm,klm) + sdd
              endif
            endif
            if (klm .eq. ii) then
              sd = tral(1,ilm,ia)/tral(3,ilm,ia)
              sc(klm,klm) = sc(klm,klm) + sd
            endif
   36     continue
          ii = ii+1
   32   continue
   30 continue

      end

