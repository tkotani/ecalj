      subroutine mkalpv(alpha,adot,iax,iop,nla,nbas,nkap,npr,nd,diag)
C- Makes diagonal part vector for cluster centered at ib
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alpha :vector, e.g. screening parms, ordered by ibas
Ci   adot  :like alpha, but energy derivative, used for iop=2
Ci   iax   :structure-constant information relating to s
Ci   iop   0: return alpha
Ci         1: return alpha^-1
Ci         2: return (alpha^-1)-dot = -adot * alpha^-2
Ci   nla   :leading dimension of alpha,adot
Ci   nkap  :Number kappas in alpha,adot
Ci   npr   :number of neighbors around each atom
Ci   nd    :leading dimenion of diag
Co Outputs:
Co   diag  :alpha ordered as in the cluster
Cu Updates
Cu   03 Aug 06 Extended to 2-kappa case.  New argument list
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer niax,iop,npr,nla,nd,nkap,nbas
      parameter (niax=10)
      integer iax(niax,npr)
      double precision diag(nd,nkap,nkap),
     .  adot(nla,nbas,nkap,nkap),alpha(nla,nbas,nkap,nkap)
C ... Local parameters
      integer ilm,ipr,ib,iv,nlmi,j,k
      double precision det

      call isanrg(iop,0,2,'MKALPV','mode',.true.)
      iv = 0
      do  10  ipr = 1, npr
        ib = iax(2,ipr)
        nlmi = iax(9,ipr)
        do  12  ilm = 1, nlmi
          iv = iv+1
          if (iop .eq. 0) then
            do  j = 1, nkap
            do  k = 1, nkap
              diag(iv,j,k) = alpha(ilm,ib,j,k)
            enddo
            enddo
          elseif (iop .eq. 1 .and. nkap .eq. 1) then
            diag(iv,1,1) = 1/alpha(ilm,ib,1,1)
          elseif (iop .eq. 1 .and. nkap .eq. 2) then
            det = alpha(ilm,ib,1,1)*alpha(ilm,ib,2,2) -
     .            alpha(ilm,ib,1,2)*alpha(ilm,ib,2,1)
            diag(iv,1,1) = alpha(ilm,ib,2,2)/det
            diag(iv,2,2) = alpha(ilm,ib,1,1)/det
            diag(iv,1,2) = -alpha(ilm,ib,1,2)/det
            diag(iv,2,1) = -alpha(ilm,ib,2,1)/det
          elseif (nkap .eq. 1) then
            diag(iv,1,1) = adot(ilm,ib,1,1)/alpha(ilm,ib,1,1)**2
          else
            call rx(' MKALPV: mode not implemented')
          endif

   12   continue
   10 continue

      end

