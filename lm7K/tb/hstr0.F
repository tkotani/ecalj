      subroutine hstr0(le,strx,nlmf,nlm,nlf,hl,nlmst,cg,indxcg,jcg,vol)
C- Make structure constants from reduced strux at energy zero
C ----------------------------------------------------------------------
Ci Inputs:
Ci   
Co Outputs:
Co   
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
      integer nlf,nlm,nlmf,nlmst,indxcg(1),jcg(1),le
      double precision cg(1),strx(nlf,nlm),hl(1),vol
      integer lmxx,icg,icg1,icg2,ii,ilm,indx,ipow,klm,l,lk,ll,llm,lm,
     .  lmax,mlm
      parameter (lmxx=28)
      double precision sig(0:lmxx),fourpi,fpibv,sum

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
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2+min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  icg = icg1, icg2
            llm = jcg(icg)
            ipow = (lm + lk - ll(llm))/2
            if (ipow .eq. 0) then
              sum = sum + cg(icg)*hl(llm)
            endif
          enddo
          strx(mlm,klm) = sum*fourpi*sig(lk)
        enddo
      enddo
C --- the following includes extra p terms 'implicitly' ---
      if (le .eq. 0) then
        fpibv = fourpi/vol
        do  ilm = 2, min0(4,nlm,nlmf)
          strx(ilm,ilm) = strx(ilm,ilm) - fpibv
        enddo
      endif
      end
