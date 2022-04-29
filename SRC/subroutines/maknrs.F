c$$$      subroutine maknrs(nbas,ips,nspec,nrspec)
c$$$C- Counts species and makes array nrspec = no. species per site
c$$$C     implicit none
c$$$      integer nbas,ips(nbas),nspec,nrspec(*)
c$$$      integer i1,ib,is
c$$$
c$$$      call imxmn(nbas,ips,1,i1,nspec)
c$$$      call iinit(nrspec,nspec)
c$$$      do  ib = 1, nbas
c$$$        is = ips(ib)
c$$$        nrspec(is) = nrspec(is) + 1
c$$$      enddo
c$$$
c$$$      end
c$$$
