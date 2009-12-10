      subroutine asadmu(nbas,nsp,lmaxu,nl,idu,nqpp,
     .  qpp,ldim,indxsh,dmatu)
C- Extract LDA+U output density-matrix from qpp matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmaxu :dimensioning parameter for U matrix
Ci   nl    :(global maximum l) + 1
Ci   idu   :idu(l+1)=1 => this l has a nonlocal U matrix
Ci   nqpp  :dimensions qpp; should be nl**2(nl**2+1)/2
Co   qpp   : matrix contains wave function products (makwts.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Co Output
Co   dmatu: output density-matrix for lda+u blocks
Cu Updates
Cu   08 Nov 07 (J. Xu) first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,lmaxu,nl,nqpp
      integer ldim,indxsh(ldim),idu(4,nbas)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      double complex qpp(nqpp,4,nsp,nbas)
C ... Local parameters
      integer isp,ib,l1,l2,m1,m2,ilm1,ilm2,iblu,iqpp,j1,j2,id

      do  isp = 1, nsp
      id = 0
      iblu = 1
      do  ib = 1, nbas
        ilm1 = 0
        iqpp = 0
C       m1 and m2 zero => no LDA+U block in ib
        call imxmn(4,idu(1,ib),1,m1,m2)
        if (m1.ne.0 .or. m2.ne.0) then
          do  l1 = 0, nl-1
            do  m1 = -l1, l1
              ilm1 = ilm1 + 1
              j1 = indxsh(id+ilm1)
              if (j1 .le. ldim) then
                ilm2 = 0
                do  l2 = 0, nl-1
                  do  m2 = -l2, l2
                    ilm2 = ilm2 + 1
                    j2 = indxsh(id+ilm2)
                    if (j2 .le. ldim) then
                      if (ilm2 .le. ilm1) then
                        iqpp = iqpp + 1
                        if (idu(l1+1,ib) .ne. 0 .and. l2 .eq. l1) then
                          dmatu(m1,m2,isp,iblu) = dmatu(m1,m2,isp,iblu)
     .                                          + qpp(iqpp,1,isp,ib)
                          if (ilm2 .ne. ilm1) then
                            dmatu(m2,m1,isp,iblu) =dmatu(m2,m1,isp,iblu)
     .                        + dconjg(qpp(iqpp,1,isp,ib))
                          endif
                        endif
                      endif
                    endif
                  enddo
                enddo
              endif
            enddo
            if (idu(l1+1,ib) .ne. 0) iblu = iblu + 1
          enddo
        end if
        id = id + nl**2
      end do
      end do

      end

