      subroutine ncutcorrect(ncut,n,gv,ng)
C- correct ncuti because of the denereracy of gv.
c ncuti give by lmfp-suham-sugcut(1,..) lmfp-bndfp-suham2-sugcut(2,..) (1 is for normal MTO, 2 is for lo)
c are only at Gamma point; thus symmetry can not be kept well for other k points.
Ci  n, gv(ng,3),ng
Cio ncut
Cr tkotani Apr2009
C----------------------
      implicit none
      integer:: i,n,ni,ng,ig, ncut(n) !,ncut0(n)
      real(8):: gv(ng,3),gg, eps=1d-6
c      ncut=100000
c      return
c      ncut0=ncut
      do i=1,n
        ni = ncut(i)
        if(ni>ng) then
           ncut(i)=ng           !takao Sep2010 correct?
           cycle
        endif
        gg = sum(gv(ni,:)**2)
        do ig= ni+1,ng
          if( abs( sum(gv(ig,:)**2)/gg -1d0)>eps) then
            ncut(i)=ig-1
            goto 99
          endif
        enddo
        ncut(i)=ng
 99     continue
      enddo
cccccccccccccccccccccc
c      write(6,"(a,100i6)")' hhh0: ng ncut=',ng,ncut0
c      write(6,"(a,100i6)")' hhh1: ng ncut=',ng,ncut
cccccccccccccccccccccc
      end
