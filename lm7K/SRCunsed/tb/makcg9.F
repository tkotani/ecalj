      subroutine makcg9(indxcg,jcg,cg,cg9,ak)
C- Make 9x9x25 matrix of "Clebsch Gordon" (actually Gaunt) coefficients 
C ----------------------------------------------------------------------
Ci Inputs: indxcg,jcg,cg
Ci   
Co Outputs: cg9
Co   
Cr Remarks
Cr  From demonstration program written by Michael.
Cu Updates
Cu 26 may 05 (ATP) For TB+U makes 
Cu A^k(LL'L''L''')=4pi/2k+1 \sum_{p=-k}^{k} C_{LL'''K} C_{L'L''K}
Cu for k=0,2,4
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer jcg(6500),indxcg(1300)
      double precision cg(6500),cg9(9,9,25),ak(9,9,9,9,3)
C Local Variables
      integer i,k,p,ii,lmax,nlm,ilm1,ilm2,ilm3,ilm4,ilmk,indx,icg,icg1,
     .        icg2,iprint
      double precision fpi,datan,f

      fpi = 16d0*datan(1d0)
      lmax = 2
      nlm = (lmax+1)**2
      call dpzero(cg9,9*9*25)

c ... loop over the first two indices
      do  ilm1 = 1, nlm
        do ilm2 = 1, nlm

c ...     Find top and bottom limits icg1,icg2 in the lists.
c         The min and max functions are used because cases like 
c         (1,3) and (3,1) are only tabulated once. 
c         What you normally want is
c            indx=(ilm1*(ilm1-1))/2+ilm2
c         where it is known that ilm1 is larger or equal to ilm2.
c         The max, min stuff always gets the proper icg1 and icg2,
c         no matter which one of ilm1 and ilm2 is larger.

          ii = max0(ilm1,ilm2)
          indx = (ii*(ii-1))/2 + min0(ilm1,ilm2)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1

c ...     loop over the relevant part of list, get out
c         the coefficient and the third index.
c         Note that ilm3 will run to higher values then either
c         ilm1 and ilm2 (twice the lmax value). 
c         If you want ilm3 to stay below nlm also, you need 
c         extra if..else  statements.

          do  icg = icg1, icg2
            ilm3 = jcg(icg)
            cg9(ilm1,ilm2,ilm3) = cg(icg)
          enddo
        enddo
      enddo
      call xxxxin(cg9)

      call dcopy(19683,0d0,0,ak,1)
      do k = 0, 4, 2
        i = k/2 + 1
        f = fpi/(2*k + 1)
        do  ilm1 = 1, nlm
          do  ilm2 = 1, nlm
            do  ilm3 = 1, nlm
              do  ilm4 = 1, nlm
                do  p = -k, k
                  ilmk = (k**2 + 1) + (p+k)
                  ak(ilm1,ilm2,ilm3,ilm4,i) = ak(ilm1,ilm2,ilm3,ilm4,i)
     .              + f * cg9(ilm1,ilm4,ilmk) * cg9(ilm2,ilm3,ilmk)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      if (iprint() .lt. 120) return
      print *, 'MAKCG9: Gaunt coefficients ...'
      do  ilm3 = 1, 25
        write (*,10) ((cg9(ilm1,ilm2,ilm3),ilm2=1,9),ilm1=1,9)
        print *
   10   format (9f10.6)
      enddo

      end
      subroutine xxxxin(cg9)
C- Shuffle 9X9X25 matrix from MSM lm indices to TBE lm indices
C ----------------------------------------------------------------------
Ci Inputs:
Ci   cg9 : 9X9X25 matrix of Clebsch Gordons
Co Outputs:
Co   cg9 reordered
Cr Remarks
Cr   Michael's structure constants are ordered according to the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     y     z     x    xy    yz 3z^2-r^2 zx   x^2-y^2
Cr   while the TBE programs use the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     x     y     z    xy    yz    zx  x^2-y^2  3z^2-r^2 
Cr   This widget rearranges the matrix made by makcg9 into
Cr   the TBE order.
Cr   The l>2 ordering is unchanged
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      double precision cg9(9,9,25)
C Local Variables
      double precision wk(9,9,25)
      integer i, j, k, ind(9)
      data ind /1,4,2,3,5,6,8,9,7/

      do  1  i = 1, 9
      do  1  j = 1, 9
      do  1  k = 1, 9
        wk(i,j,k) = cg9(ind(i),ind(j),ind(k))
    1 continue 
      do  2  i = 1, 9
      do  2  j = 1, 9
      do  2  k = 10, 25
        wk(i,j,k) = cg9(ind(i),ind(j),k)
    2 continue 
      call dcopy(9*9*25,wk,1,cg9,1)

      end
