      subroutine hstr(pv,le,strx,drstrx,nlmf,nlm,nlf,hl,nlmst,cg,
     .                indxcg,jcg,vol)
C- Make structure constants from reduced strux at energy zero
C ----------------------------------------------------------------------
Ci Inputs: logical pv,le; nlmf,nlm,nlf,hl,nlmst,cg,indxcg,jcg,vol
Ci   
Co Outputs: strx,drstrx
Co   
Cr Remarks B_LL'(tau) = 4\pi \sum_L" (-1)^l' H_l"(tau) Y_L"(tau)
Cr         HY are reduced structure constants (see rcnsl0, soldh)
Cr         If pv is set, returns tau*d/dtau B in drstrx
C ----------------------------------------------------------------------
C     implicit none
      integer nlf,nlm,nlmf,nlmst,indxcg(1),jcg(1),le
      double precision cg(1),strx(nlf,nlm),drstrx(nlf,nlm),hl(1),vol
      integer lmxx,icg,icg1,icg2,ii,ilm,indx,ipow,klm,l,lk,ll,llm,lm,lp,
     .  lmax,mlm
      parameter (lmxx=28)
      double precision sig(0:lmxx),fourpi,fpibv,sum,sumr
      logical pv

      call tcn('hstr0: make strux')
      fourpi = 16d0*datan(1d0)
      lmax = ll(nlmf) + ll(nlm)
      if (lmax .gt. 12) call rx('change dimensions in hstr0''')
C --- (-1)^l ---
      sig(0) = 1d0
      do  l = 1, lmax
        sig(l) = -sig(l-1)
      enddo
C --- add together Gaunt-sums ---
      do  mlm = 1, nlmf
        lm = ll(mlm)
        do  klm = 1, nlm
          lk = ll(klm)
          sum = 0d0
          sumr = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2+min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  icg = icg1, icg2
            llm = jcg(icg)
            lp = ll(llm)
            ipow = (lm + lk - lp)/2
            if (ipow .eq. 0) then
              sum = sum + cg(icg)*hl(llm)
              if (pv) then
                sumr = sumr - (lp+1) * cg(icg)*hl(llm)
              endif
            endif
          enddo
          strx(mlm,klm) = sum*fourpi*sig(lk)
          if (pv) then
            drstrx(mlm,klm) = sumr*fourpi*sig(lk)
          endif
        enddo
      enddo
C --- the following includes extra p terms 'implicitly' ---
      if (le .eq. 0) then
        fpibv = fourpi/vol
        do  ilm = 2, min0(4,nlm,nlmf)
          strx(ilm,ilm) = strx(ilm,ilm) - fpibv
          if (pv) then
            drstrx(ilm,ilm) = drstrx(ilm,ilm) - fpibv
          endif
        enddo
      endif
      call tcx('hstr0: make strux')
      end
