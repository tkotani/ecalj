subroutine rotcg(lmxax,symops,ng,cgr)
  !! -- Rotated CG coefficients.
  !! output: cgr
  !! NOTE: key is the line
  !!      cgr(lm1,lm2,lm,ig) = sum(cgn(lx-l:lx+l,lm1,lm2)*dlmm(-l:l,m,l,ig))
  !! Here, dlmm(m,m') is the rotation matrix  of angular momentum space lm
  !! for given symops(3,3,ig).
  use m_lldata,only: ll
  implicit none
  integer(4) :: lmxax, ng, nlmxa, &
       lnjcg, lnxcg, &
       ilma,la,ilmb,lh,ii,indx,icg1,icg2,icg, &
       ig,lm1,lm2,lm,l,m,md,lmd,lmxcg,ilm ,lx
  real(8) :: &
       cgr((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2,ng), &
       symops(9,ng)        ,sumr
  real(8),allocatable:: cg(:),dlmm(:,:,:,:),cgn(:,:,:)
  integer(4),allocatable :: jcg(:),indxcg(:)
  integer:: ll1,ll2,lxm
  ! --- CG coefficienets. <LM3|lm1 lm2>
  ! inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
  ! Injcg = indxcg(inxcg) to indxcg(inxcg)-1
  ! cg(inxcg)  : = <lm3|lm1 lm2>
  ! jcg(lnjcg) : = lm3
  !c      write(6,*)' rotcg:'
  !      do ig=1,ng
  !      write(6,*)' transposed symope ig  =',ig
  !     write(6,'(3f12.6)') symops(1:3,ig)
  !     write(6,'(3f12.6)') symops(4:6,ig)
  !     write(6,'(3f12.6)') symops(7:9,ig)
  !      enddo
  ll2 = (2*lmxax+1)**2
  ll1 = (lmxax+1)**2
  allocate( cgn(ll2,ll1,ll1) )
  cgn = 0d0
  nlmxa = (lmxax+1)**2
  lmxcg = lmxax
  if (lmxcg <= 6) then
     lnjcg = 6500
     lnxcg = 1300
  else if (lmxcg <= 8) then
     lnjcg = 22700
     lnxcg = 3400
  else if (lmxcg <= 10) then
     lnjcg = 62200
     lnxcg = 7400
  else
     write(6,*) 'rotcg: cannot handle lmxcg=',lmxcg
     call rx( 'rotcg: cannot handle lmxcg')
  endif
  allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
  call scg(lmxcg,cg,indxcg,jcg)
  !----
  do ilma = 1, nlmxa
     la = ll(ilma)
     do ilmb = 1, nlmxa
        lh = ll(ilmb)
        ii = max0(ilma,ilmb)
        indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1)-1
        do icg = icg1, icg2
           ilm  = jcg(icg)
           cgn(ilm, ilma,ilmb)  = cg(icg) ! ilm is move to 1st argument.!
        enddo
     enddo
  enddo

  !! --- no Rotation case
  ig=1
  if(ng==1 .AND. sum(abs(symops(:,ig)-[1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0]))<1d-12) then
     do lm =  1, (2*lmxax+1)**2
        l = ll(lm)
        m = lm - l**2 -l -1
        lxm = l**2 +l +1 +m
        ig  = 1
        do lm2 = 1, nlmxa
           do lm1 = 1, nlmxa
              cgr(lm1,lm2,lm,1) = cgn(lxm,lm1,lm2)
           enddo
        enddo
     enddo
     deallocate( cg,cgn,jcg,indxcg)
     return
  endif

  !! --- Rotation matrix
  allocate(dlmm( -2*lmxax:2*lmxax, -2*lmxax:2*lmxax, 0:2*lmxax,ng))
  call rotdlmm(symops, ng, 2*lmxax+1,dlmm)
  !! --- Rotated CG
  do lm =  1, (2*lmxax+1)**2
     l = ll(lm)
     m = lm - l**2 -l -1
     lx = l**2 +l +1
     do ig  = 1, ng
        do lm2 = 1, nlmxa
           do lm1 = 1, nlmxa
              ! ilm is move to 1st argument.!
              cgr(lm1,lm2,lm,ig) = sum(cgn(lx-l:lx+l,lm1,lm2)*dlmm(-l:l,m,l,ig))
              !        sumr = 0d0
              !        do md = -l,l
              !          lmd = l**2 +l +1 + md
              !          sumr = sumr + cgn(lm1,lm2,lmd)*dlmm(md,m,l,ig)
              !        enddo
              !        cgr(lm1,lm2,lm,ig) = sumr
           enddo
        enddo
     enddo
  enddo
  deallocate( cg,dlmm,cgn,jcg,indxcg)
  !c      write(6,*)' rotcg end:'
end subroutine rotcg
