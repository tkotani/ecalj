subroutine ncutcorrect(ncut,n,gv,ng)
  !- correct ncuti because of the denereracy of gv.
  ! ncuti give by lmfp-suham-sugcut(1,..) lmfp-bndfp-suham2-sugcut(2,..) (1 is for normal MTO, 2 is for lo)
  ! are only at Gamma point; thus symmetry can not be kept well for other k points.
  !i  n, gv(ng,3),ng
  ! o ncut
  !r tkotani Apr2009
  !----------------------
  implicit none
  integer:: i,n,ni,ng,ig, ncut(n) !,ncut0(n)
  real(8):: gv(ng,3),gg, eps=1d-6
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
99   continue
  enddo
end subroutine ncutcorrect
