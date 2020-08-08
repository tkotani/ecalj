      subroutine maknrs(nbas,ips,nspec,nrspec)
C- Counts species and makes array nrspec = no. species per site
C     implicit none
      integer nbas,ips(nbas),nspec,nrspec(*)
      integer i1,ib,is

      call imxmn(nbas,ips,1,i1,nspec)
      call iinit(nrspec,nspec)
      do  ib = 1, nbas
        is = ips(ib)
        nrspec(is) = nrspec(is) + 1
      enddo

      end

