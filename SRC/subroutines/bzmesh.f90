subroutine bzmesh(plat,qb,ifac,n1,n2,n3,lshft,g,ng,ipq,qp,wgt,nq,nqmx)! Divides the reciprocal lattice into microcells
  use m_lgunit,only:stdo
  use m_ext,only: sname
  !i  plat     :primitive lattice vectors
  !i  n1,n2,n3 :no. divisions for the 3 recip. latt. vecs; (see Remarks)
  !i  g,ng     :symmetry point group operations, and number
  !i  nqmx     : number of k-points exceeds this maximum
  !i  lshft(1:3)  :logical switch, for each recip. latt. vec. 
  !i           T:center the mesh straddling the origin --- no point at (0,0,0).
  !o Outputs:
  !o   ipq   :ipq(i1,i2,i3) points to irreducible qp corresponding to
  !o         :qp labelled (i1,i2,i3), where i1=1..n1, i2=1..n2, i3=1..n3
  !    ipq(i1,i2,i3) marks which irreducible qp to which each point in the
  !    full BZ belongs: point (i1,i2,i3) is equivalent to irreducible point ipq(i1,i2,i3).
  !o   wgt   :weight assigned to this qp 
  !    wgt(j) contains the sampling weight associated with qp(j), i.e.
  !    2/(n1*n2*n3) * no. points equivalent to qp(j); factor 2 for spin.
  !o   nq    :number of irreducible qp
  !o   qb    :vectors of first microcell for input to BZINTS (see bzmsh0)
  !o   qp
  !    Some of the qp will be symmetry-related, leaving nq irreducible k-points, which are returned in qp(1..3,j), j=1,nq.
  !o   ifac(i)    1 if not shifted, 2 if shifted for i'th axis; see qb.
  !o   rb:        a microcell in the Wigner Seitz cell.
  !o   qb:        a microcell in the Brillouin zone
  !o              This, together with ifac provides information how to generate the actual q-point from a triplet of integers
  !o              specifying a point on the mesh.  Given a triplet (j_1,j_2,j_3) of ipq, the i'th component q_i is
  !o                 q_i(j_1,j_2,j_3) = sum_n (j'_n)*qb(i,n),  with      j'_n = j_n*ifac(n)-1
  !   is(i)      0 mesh points centered at the origin for i'th axis
  !              1 mesh points off-centered from the origin for i'th axis
  !NOTE:
  !  The reciprocal lattice is divided into n1*n2*n3 microcells which are parallelipipeds with 8 corners. The corners are nodes of the
  !  k-space mesh in the whole reciprocal lattice unit cell. Thus, for i1=1..n1, i2=1..n2, i3=1..n3 the qp(i1,i2,i3) are
  !  q_k = (i1*ifac(1)-1)*qb(k,1) +  (i2*ifac(2)-1)*qb(k,2) +  (i3*ifac(3)-1)*qb(k,3),    where ifac is 1 or 2; see bzmsh0.
  implicit none
  logical :: lshft(3)
  integer :: n1,n2,n3,nqmx,ng,nq,ipq(n1,n2,n3)
  double precision :: qb(3,3),wgt(nqmx),plat(3,3),qp(3,nqmx),g(3,3,*)
  integer :: i1,i2,i3,ifac(3),ig,igcnt,ii,ii1,ii2,ii3,ipr,iq,is(3),iwgt,jj(3),lgunit,m1,m2,m3,ndmx,nnn(3),mmm(3),ifi
  double precision :: w0,swgt,v(3),v1(3),rb(3,3),xx(3)
  character(1) :: chr(0:2)
  real(8):: tolq
  call getpr(ipr)
  bzmesh0: block
    integer:: k,m,iprint,mvec(3)
    real(8):: vol,qlat(3,3)
    call dinv33(plat,1,qlat,vol)
    is   = [(merge(1,0,lshft(m)),m=1,3)]
    ifac = [(merge(2,1,lshft(m)),m=1,3)]
    mvec = [n1,n2,n3]*ifac
    do  m = 1, 3
       qb(m,:) = qlat(m,:)/mvec
       rb(m,:) = plat(m,:)*mvec
    enddo
    if(iprint() > 80) then
       write(stdo,"(' BZMESH : ',5X,'Plat',31X,'Qlat')")
       do k = 1, 3
          write(stdo,"(3f10.5,5x,3f10.5)") (plat(m,k),m=1,3),(qlat(m,k),m=1,3)
       enddo
    endif
  endblock bzmesh0
  ipq=0 
  nq = 0
  swgt = 0d0
  igcnt = 0
  w0 = 2d0/(n1*n2*n3)
  nnn = [n1,n2,n3]
  mmm = 6*nnn*ifac
  i3loop: do 23 i3=0,n3-1 !For each of (n1*n2*n3) qp, find irreducible set ---
     do 22      i2=0,n2-1
        do 21   i1=0,n1-1   
           if(ipq(i1+1,i2+1,i3+1) == 0) then ! Add qp to list if not flagged as symmetry-related to a prior
              v(:) = matmul(qb(:,:),[i1*ifac(1)+is(1), i2*ifac(2)+is(2), i3*ifac(3)+is(3)])
              iwgt = 0
              v1=v
              igloop: do ig = 1, max(ng,1)
                 if (ng > 0) v1=matmul(g(:,:,ig),v) 
                 xx = matmul(v1(:),rb(:,:))-is
                 jj = nint(xx)
                 if(sum(abs(xx-jj)) > tolq()) then
                    write(stdo,"(a,3f9.4,' ',3f9.4)") ' qp mapped to is not on k-mesh',v,v1
                    write(stdo,"(a,3f9.4,' ',3i5)")   '             x j=',xx,jj(1),jj(2),jj(3)
                    open(newunit=ifi,file='bzmesh.'//trim(sname)//'.err')
                    write(ifi,*)'BZMESH: symops incompatible with this mesh'
                    close(ifi)
                    call rx('BZMESH: symops incompatible with this mesh')
                 endif
                 if(any(lshft(:).AND.mod(abs(jj(:)),2)== 1)) cycle ! discard if shifted off mesh
                 if(lshft(1)) jj(1) = jj(1)/2 
                 if(lshft(2)) jj(2) = jj(2)/2
                 if(lshft(3)) jj(3) = jj(3)/2
                 jj = mod(jj+2*mmm,nnn) + 1 !Ensure (jj(1),jj(2),jj(3)) in first quadrant of Q
                 call rxx(any(jj<=0),'neg j in bzmesh')
                 if(ipq(jj(1),jj(2),jj(3)) == 0) then
                    ipq(jj(1),jj(2),jj(3)) = nq+1
                    iwgt = iwgt+1
                    igcnt = igcnt+1
                 endif
              enddo igloop
              nq = nq+1
              qp(:,nq) = v
              wgt(nq) = iwgt*w0
              swgt = swgt + abs(wgt(nq))
           endif
21      enddo
22   enddo
23 enddo i3loop
  if(igcnt/=n1*n2*n3 ) call rx('bug in bzmesh')
  if(abs(swgt-2)>1d-9) call rx1('BZMESH: QP weights sum to ',swgt)
  if(ipr>=20)write(stdo,"(a,i3,i5,a,3i4,a,3l)") " BZMESH: ngrp nq ",ng,nq," QP from ",n1,n2,n3," shift=",lshft
  if(ipr>=50) then
     chr(2) = ' '
     chr(0) = '*'
     write(stdo,"(14x,'Qx',10x,'Qy',10x,'Qz',6x,'Multiplicity    Weight')")
     do iq = 1, nq
        ii = 1+dsign(1d0,wgt(iq))
        iwgt = abs(wgt(iq)/w0) + .1d0
        write(stdo,"(i5,2x,3f12.6,i10,1x,a,f14.6)") iq,qp(1,iq),qp(2,iq),qp(3,iq),iwgt,chr(ii),abs(wgt(iq))
     enddo
  endif
end subroutine bzmesh
