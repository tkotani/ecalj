module m_gaugm !Generic routine to make augmentation matrices
  use m_lmfinit,only: n0,nab
  public gaugm
  private
contains
  subroutine gaugm(nr,nsp,lso,rofi,rwgt,lmxa,lmxl,nlml,vsm,gpotb,gpot0,hab,vab,sab,sodb,qum,vum,&
       lmaxu,vumm,lldau,idu,lmux, &
       nf1,nf1s,lmx1,lx1,nlx1,f1,x1,v1,d1, &
       nf2,nf2s,lmx2,lx2,nlx2,f2,x2,v2,d2, &
       sig,tau,ppi,hsozz,hsopm) !OUTPUT: augmentation integrals. 
    use m_lmfinit,only: cg=>rv_a_ocg,jcg=>iv_a_ojcg,indxcg=>iv_a_oidxcg
    use m_lgunit,only: stdo
    implicit none
    integer :: nr,nsp,lmxa,lmxl,nlml,lmux,lso, & 
         nf1,nf1s,lmx1,nlx1,lx1(nf1),nf2,nf2s,lmx2,nlx2,lx2(nf2)
    integer :: lmaxu,lldau,idu(4)
    complex(8):: vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,3,3,2,0:lmaxu)
    real(8) :: rofi(nr),rwgt(nr),vsm(nr,nlml,nsp), &
         qum(0:lmxa,0:lmxa,0:lmxl,3,3,nsp),gpotb(1),gpot0(1), &
         hab(3,3,n0,nsp),vab(3,3,n0,nsp),sab(3,3,n0,nsp), &
         f1(nr,0:lmx1,nf1s),v1(0:lmx1,nf1s),d1(0:lmx1,nf1s), &
         f2(nr,0:lmx2,nf2s),v2(0:lmx2,nf2s),d2(0:lmx2,nf2s), &
         x1(nr,0:lmx1,nf1s),x2(nr,0:lmx2,nf2s), &
         sig(nf1,nf2,0:lmux,nsp),tau(nf1,nf2,0:lmux,nsp),vum(0:lmxa,0:lmxa,nlml,3,3,nsp)
    double precision :: sodb(3,3,n0,nsp,2) !Spin-Orbit related
    integer :: ilm1,ilm2,ix
    integer :: i1,i2,ilm,l,ll,nlm,nlm1,nlm2,i,iprint,ii
    real(8),parameter:: pi= 4d0*datan(1d0),y0= 1d0/dsqrt(4d0*pi)
    real(8):: vsms(nr),ppi0(nf1,nf2,0:lmux), &
         qm(nf1*nf2*(lmx1+1)*(lmx2+1)*(lmxl+1)), ssum((lmx1+1)*(lmx2+1)*nlml)
    complex(8),allocatable:: hso(:,:,:,:,:,:),ppiz(:,:,:,:,:)
    complex(8):: ppi  (nf1,nf2,nlx1,nlx2,nsp)
    complex(8):: hsozz(nf1,nf2,nlx1,nlx2,nsp)
    complex(8):: hsopm(nf1,nf2,nlx1,nlx2,nsp) !offdiag parts for <isp|ispo> (see sodb(*,2))
    real(8)::   ppir(nf1,nf2,nlx1,nlx2,nsp),mmm
    allocate(ppiz(nf1,nf2,nlx1,nlx2,nsp))
    if(lso/=0)  hsozz=0d0
    if(lso/=0)  hsopm=0d0
    ppir=0d0
    sig=0d0 
    tau=0d0 
    ppi0=0d0 
    do i = 1, nsp
       vsms=y0*vsm(1:nr,1,i)!Spherical part of the smooth potential
       ! --- Make sig, tau, ppi0 = spherical part of ppi ---
       call pvagm(nr,rofi,rwgt,vsms,&
            f1,f2,x1,x2,&
            nf1,nf1s,lmx1,lx1,nlx1,v1,d1,& !f1~f2~ part of sig and corresponding tau,ppi0
            nf2,nf2s,lmx2,lx2,nlx2,v2,d2,&
            i,lso,lmux, &
            hab(1,1,1,i),vab(1,1,1,i),sab(1,1,1,i),sodb(1,1,1,i,1),sodb(1,1,1,i,2), &
            sig(1,1,0,i),tau(1,1,0,i),ppi0,hsozz(1,1,1,1,i),hsopm(1,1,1,1,i)) 
       call ppiNS(nr,nlml,vsm(1,1,i),rwgt,cg,jcg,indxcg,& ! ppi from Non-Spherical potential ---
            nf1,nf1s,lmx1,lx1,f1, & !smooth integral f1^ (-vsm) f2^ for nonspherical part of vsm
            nf2,nf2s,lmx2,lx2,f2,ssum,nlx1,nlx2,& 
            v1,d1,v2,d2,lmxa,vum(0,0,1,1,1,i), & 
            ppir(1,1,1,1,i)) ! integral f1~ vtrue f2~ for nonspherical part of vtrue
       call getqm(nf1,nf1s,lmx1,lx1,f1,v1,d1,nf2,nf2s,lmx2,lx2,f2,v2,d2,&
            nr,rofi,rwgt,lmxa,qum(0,0,0,1,1,i),lmxl, &
            qm)! Moments qm = (f1~*f2~ - f1^*f2^) r^m Y_m ---
       call ppiall(nf1,lmx1,lx1,nf2,lmx2,lx2,lmxl,nlml,cg,jcg, & 
            indxcg,qm,gpotb,gpot0,lmux,ppi0,nlx1,nlx2,& !qm * (gpot0-gpotb) added to ppi
            ppir(1,1,1,1,i)) !Assemble ppi from ppi0, NS part and multipole qm contr
       do  i1 = 1, nf1 ! --- Add tau into ppi ---
          do  i2 = 1, nf2
             do  ilm = 1, min0((lx1(i1)+1)**2,(lx2(i2)+1)**2)
                ppir(i1,i2,ilm,ilm,i) = ppir(i1,i2,ilm,ilm,i) + tau(i1,i2,ll(ilm),i)
             enddo
          enddo
       enddo
    enddo
    if(lldau /= 0) then !LDA+U contribution.
       ppiz=0d0
       do i=1,nsp
          call ppildau(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2,lmaxu,vumm,nlx1,nlx2,&
               ppiz(1,1,1,1,i),i,idu)
       enddo
       ppi=ppir+ppiz
    else
       ppi = ppir
    endif
    deallocate(ppiz)
  end subroutine gaugm
  subroutine pvagm(nr,rofi,rwgt,vsms,& ! sig, tau, ppi (spherical part), and SO part.
       f1,f2,x1,x2,&
       nf1,nf1s,lmx1,lx1,nlx1,v1,d1,&
       nf2,nf2s,lmx2,lx2,nlx2,v2,d2,&
       isp,lso,lmux, &
       hab,vab,sab,sodb,sondb,&
       sig,tau,ppi,hsozz,hsopm) !output
    implicit none
    integer :: lmux,lmx1,lmx2,nf1,nf2,nf1s,nf2s,lx1(nf1),lx2(nf2), nlx1,nlx2,isp,nr,ii
    double precision :: hab(3,3,0:n0-1),sab(3,3,0:n0-1),vab(3,3,0:n0-1), &
         v1(0:lmx1,nf1),d1(0:lmx1,nf1), &
         v2(0:lmx2,nf2),d2(0:lmx2,nf2), &
         sig(nf1,nf2,0:lmux),tau(nf1,nf2,0:lmux),ppi(nf1,nf2,0:lmux)
    double precision :: sodb(3,3,0:n0-1),sondb(3,3,0:n0-1)
    double precision :: a1,a2,vd1(2),vd2(2)
    complex(8):: hsozz(nf1,nf2,nlx1,nlx2),hsopm(nf1,nf2,nlx1,nlx2),somatpm,somatzz
    integer :: i1,i2,lmax1,lmax2,lmax,l
    integer :: m1,m2,l1,l2,lso
    complex(8),parameter:: img=(0d0,1d0)
    real(8)::sim,tum,xbc,rwgt(nr),vsms(nr),rofi(nr),&
         f1(nr,0:lmx1,nf1),x1(nr,0:lmx1,nf1), &
         f2(nr,0:lmx2,nf2),x2(nr,0:lmx2,nf2) 
    real(8):: sorad
    sig=0d0 
    tau=0d0 
    ppi=0d0 
    do  i1 = 1, nf1
       do  i2 = 1, nf2
          do  l = 0, min0(lx1(i1),lx2(i2))
             sig(i1,i2,l) = sum([(rwgt(ii)*f1(ii,l,i1)*f2(ii,l,i2),ii=2,nr)])
             ppi(i1,i2,l) = sum([(rwgt(ii)*f1(ii,l,i1)*f2(ii,l,i2)*vsms(ii),ii=2,nr)])
             sim  =  sum([(rwgt(ii)*f1(ii,l,i1)*f2(ii,l,i2)/rofi(ii)**2,ii=2,nr)])
             tum  =  sum([(rwgt(ii)*x1(ii,l,i1)*x2(ii,l,i2),ii=2,nr)])
             xbc = f1(nr,l,i1) * x2(nr,l,i2)
             tau(i1,i2,l) = tum + l*(l+1)*sim - xbc
             !       ... hzz,vzz,szz
             if (i1 > nf1s .AND. i2 > nf2s) then
                sig(i1,i2,l) = -sig(i1,i2,l) + sab(3,3,l)
                tau(i1,i2,l) = -tau(i1,i2,l) + (hab(3,3,l)-vab(3,3,l))
                ppi(i1,i2,l) = -ppi(i1,i2,l) + vab(3,3,l)
                !       ... hzu,vzu,szs
             elseif (i1 > nf1s) then
                sig(i1,i2,l)= -sig(i1,i2,l) + sum(sab(3,1:2,l)*[v2(l,i2),d2(l,i2)])
                tau(i1,i2,l)= -tau(i1,i2,l) + sum((hab(3,1:2,l)-vab(3,1:2,l))*[v2(l,i2),d2(l,i2)])
                ppi(i1,i2,l) = -ppi(i1,i2,l)+ sum(vab(3,1:2,l)*[v2(l,i2),d2(l,i2)])
                !       ... huz,vuz,ssz
             elseif (i2 > nf2s) then
                sig(i1,i2,l) = -sig(i1,i2,l) +sum([v1(l,i1),d1(l,i1)]*sab(1:2,3,l))
                tau(i1,i2,l) = -tau(i1,i2,l) +sum([v1(l,i1),d1(l,i1)]*(hab(1:2,3,l)-vab(1:2,3,l)))
                ppi(i1,i2,l) = -ppi(i1,i2,l) +sum([v1(l,i1),d1(l,i1)]*vab(1:2,3,l))
             else
                vd1= [v1(l,i1),d1(l,i1)]
                vd2= [v2(l,i2),d2(l,i2)]
                sig(i1,i2,l)= -sig(i1,i2,l) + sum(vd1*matmul(sab(1:2,1:2,l),vd2))
                tau(i1,i2,l)= -tau(i1,i2,l) + sum(vd1*matmul(hab(1:2,1:2,l)-vab(1:2,1:2,l),vd2))
                ppi(i1,i2,l)= -ppi(i1,i2,l) + sum(vd1*matmul(vab(1:2,1:2,l),vd2))
             endif
          enddo
       enddo
    enddo
    
    if(lso==0) return  !SO part from now on ===================
    hsopm=0d0  
    hsozz=0d0
    i1loop2: do  i1 = 1, nf1
       i2loop2: do  i2 = 1, nf2
          l1 = 0
          l2 = 0
          lloop2: do  l = 0, min0(lx1(i1),lx2(i2))
             mloop2: do  m1 = -l, l               !l1 = l1 + 1
                l1 = l**2 + l+1+ m1     !l1 + 1   !if (m1 >= (-l+1)) l2 = l2 - (2*l + 1)
                do  m2 = -l, l
                   l2 = l**2 + l+1 + m2 !l2 + 1    !l2 = l2 + 1
                   vd1= [v1(l,i1),d1(l,i1)]
                   vd2= [v2(l,i2),d2(l,i2)]
                   if ( m1+m2==0.AND.m2/=0.and.l1/=l2) then ! ... LzSz part
                      if (i1 > nf1s .AND. i2 > nf2s) then !     ... hso_zz
                         sorad = sodb(3,3,l)
                      elseif (i1 > nf1s) then                      !             ... hso_zu
                         sorad = sum(sodb(3,1:2,l)*[v2(l,i2),d2(l,i2)])
                      elseif (i2 > nf2s) then                      !             ... hso_uz
                         sorad = sum([v1(l,i1),d1(l,i1)]*sodb(1:2,3,l))
                      else
                         sorad = sum(vd1*matmul(sodb(1:2,1:2,l),vd2))
                      endif
                      if (l1 < l2 .AND. m1 == -m2) somatzz =  abs(m1)*img 
                      if (l1 > l2 .AND. m1 == -m2) somatzz = -abs(m1)*img
                      hsozz(i1,i2,l1,l2) = somatzz * sorad
                   endif
                   if (lso == 1) then ! LxSx+LySy part
                      if (i1 > nf1s .AND. i2 > nf2s) then
                         sorad = sondb(3,3,l)
                      elseif (i1 > nf1s) then                         !         ... hso_zu
                         sorad = sum(sondb(3,1:2,l)*[v2(l,i2),d2(l,i2)])
                      elseif (i2 > nf2s) then                         !         ... hso_uz
                         sorad = sum([v1(l,i1),d1(l,i1)]*sondb(1:2,3,l))
                      else
                         sorad = sum(vd1*matmul(sondb(1:2,1:2,l),vd2))
                      endif
                      call mksomat(l,m1,m2,isp,somatpm) ! LxSx+LySy lm part
                      hsopm(i1,i2,l1,l2)= somatpm* sorad
                   endif
                enddo
             enddo mloop2
          enddo lloop2
       enddo i2loop2
    enddo i1loop2
  end subroutine pvagm
  subroutine ppiNS(nr,nlml,vsm,rwgt,cg,jcg,indxcg, &
       nf1,nf1s,lmx1,lx1,f1,&
       nf2,nf2s,lmx2,lx2,f2,ssum,nlx1,nlx2,& 
       v1,d1,v2,d2, lmxa,vum,&
       ppi) !- Add non-spherical constribution to ppi.
    implicit none
    integer :: lmx1,lmx2,nf1,nf2,nf1s,nf2s,nlml,nlx1,nlx2,nr
    integer :: lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
    double precision :: vsm(nr,nlml),rwgt(nr), &
         ppi(nf1,nf2,nlx1,nlx2),f1(nr,0:lmx1,nf1),f2(nr,0:lmx2,nf2), &
         cg(1),ssum(0:lmx1,0:lmx2,nlml)
    integer :: i1,i2,icg,ilm1,ilm2,ix,l1,l2,ll,mlm,nlm1,nlm2,i
    double precision :: sam
    integer :: lmxa 
    double precision :: vum(0:lmxa,0:lmxa,nlml,3,3), &
         v1(0:lmx1,nf1),d1(0:lmx1,nf1), &
         v2(0:lmx2,nf2),d2(0:lmx2,nf2),vd1(2),vd2(2)
    double precision :: add=1d99
    ppi=0d0 
    do  i1 = 1, nf1 ! Sum over CG coefficients, make radial integrals as needed
       do  i2 = 1, nf2
          nlm1 = (lx1(i1)+1)**2
          nlm2 = (lx2(i2)+1)**2
          ssum(0:lx1(i1),0:lx2(i2),1:nlml) = 2d10 ! --- integral not yet calculated
          do  ilm1 = 1, nlm1
             l1 = ll(ilm1)
             do  ilm2 = 1, nlm2
                l2 = ll(ilm2)
                ix = max0(ilm1,ilm2)
                ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
                do icg = indxcg(ix),indxcg(ix+1)-1
                   mlm = jcg(icg)
                   if (mlm > 1 .AND. mlm <= nlml) then
                      if (ssum(l1,l2,mlm) > 1d10) then
                         ssum(l1,l2,mlm) = sum([(rwgt(i)*vsm(i,mlm)*f1(i,l1,i1)*f2(i,l2,i2),i=2,nr)])
                      endif
                      ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) - cg(icg)*ssum(l1,l2,mlm)
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    do  i1 = 1, nf1
       do  i2 = 1, nf2
          nlm1 = (lx1(i1)+1)**2
          nlm2 = (lx2(i2)+1)**2
          do  ilm1 = 1, nlm1
             l1 = ll(ilm1)
             do  ilm2 = 1, nlm2
                l2 = ll(ilm2)
                ix = max0(ilm1,ilm2)
                ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
                do  icg = indxcg(ix),indxcg(ix+1)-1
                   mlm = jcg(icg)
                   vd1=[v1(l1,i1),d1(l1,i1)]
                   vd2=[v2(l2,i2),d2(l2,i2)]
                   if (mlm > 1 .AND. mlm <= nlml) then
                      if (i1 > nf1s .AND. i2 > nf2s) then!    <g | V | g>
                         add = vum(l1,l2,mlm,3,3)
                      elseif (i1 > nf1s) then  !              <g | V | (u,s)>
                         add = sum( vum(l1,l2,mlm,3,1:2) * [v2(l2,i2),d2(l2,i2)])
                      elseif (i2 > nf2s) then !               <(u,s) | V | g>
                         add = sum([v1(l1,i1),d1(l1,i1)] * vum(l1,l2,mlm,1:2,3))
                      else
                         add= sum(vd1*matmul(vum(l1,l2,mlm,1:2,1:2),vd2))
                      endif
                      ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) + cg(icg)*add
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ppiNS
  subroutine getqm(nf1,nf1s,lmx1,lx1,f1,v1,d1,nf2,nf2s,lmx2,lx2,f2, &
       v2,d2,nr,rofi,rwgt,lmxa,qum,lmxl,qm) !output = qm = integrals (f1~*f2~ - f1*f2) r**l
    implicit none
    integer :: nf1,nf2,nf1s,nf2s,lmx1,lmx2,lmxa,lmxl,nr
    integer :: lx1(nf1),lx2(nf2)
    double precision :: rofi(nr),rwgt(nr), &
         qm(nf1,nf2,0:lmx1,0:lmx2,0:lmxl), &
         f1(nr,0:lmx1,nf1),v1(0:lmx1,nf1s),d1(0:lmx1,nf1s), &
         f2(nr,0:lmx2,nf2),v2(0:lmx2,nf2s),d2(0:lmx2,nf2s), &
         qum(0:lmxa,0:lmxa,0:lmxl,3,3),vd1(2),vd2(2)
    integer :: i1,i2,l1,l2,lm,i
    double precision :: ssum,sam
    qm=0d0 
    do  i1 = 1, nf1
       do  i2 = 1, nf2
          do  l1 = 0, lx1(i1)
             do  l2 = 0, lx2(i2)
                do  lm = 0, lmxl
                   ssum = sum([(rwgt(i)*f1(i,l1,i1)*f2(i,l2,i2)*rofi(i)**lm,i=2,nr)])
                   vd1=[v1(l1,i1),d1(l1,i1)]
                   vd2=[v2(l2,i2),d2(l2,i2)]
                   if (i1 > nf1s .AND. i2 > nf2s) then ! Both f1~ and f2~ are local orbitals
                      sam = qum(l1,l2,lm,3,3) 
                   elseif (i1 > nf1s) then!  f1~ is local, f2~ is linear combination of (u,s)
                      sam = sum([v2(l2,i2),d2(l2,i2)] * qum(l2,l1,lm,1:2,3))
                   elseif (i2 > nf2s) then !f1~ is linear combination of (u,s), f2~ is local
                      sam = sum([v1(l1,i1),d1(l1,i1)] * qum(l1,l2,lm,1:2,3))
                   else
                      sam = sum(vd1*matmul(qum(l1,l2,lm,1:2,1:2),vd2))
                   endif
                   if(sam==0.AND.ssum/=0)call rx('gaugm: inconsistent treatment of local orbitals')
                   qm(i1,i2,l1,l2,lm) = sam-ssum
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine getqm
  subroutine ppiall(nf1,lmx1,lx1,nf2,lmx2,lx2,lmxl,nlml, &
       cg,jcg,indxcg,qm,gpotb,gpot0,lmux,ppi0,nlx1,nlx2,&
       ppi)!- Assemble ppi0,qm ppi to ppi
    implicit none
    integer :: lmux,lmx1,lmx2,lmxl,nf1,nf2,nlml,nlx1,nlx2
    integer :: lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
    integer :: i1,i2,nlm1,nlm2,ilm1,ilm2,l1,l2,ll,ix,icg,mlm,lm
    real(8):: ppi0(nf1,nf2,0:lmux),ppi(nf1,nf2,nlx1,nlx2), &
         cg(1),gpotb(nlml),gpot0(nlml),qm(nf1,nf2,0:lmx1,0:lmx2,0:lmxl)
             ! ... Add terms from moments of f~g~ - fg and ppi0 to ppi
    do  i1 = 1, nf1
       do  i2 = 1, nf2
          nlm1 = (lx1(i1)+1)**2
          nlm2 = (lx2(i2)+1)**2
          do  ilm1 = 1, nlm1
             l1 = ll(ilm1)
             if (ilm1 <= nlm2) ppi(i1,i2,ilm1,ilm1) = ppi(i1,i2,ilm1,ilm1) + ppi0(i1,i2,l1)
             do  ilm2 = 1, nlm2
                l2 = ll(ilm2)
                ix = max0(ilm1,ilm2)
                ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
                do  icg = indxcg(ix),indxcg(ix+1)-1
                   mlm = jcg(icg)
                   if (mlm <= nlml) then
                      lm = ll(mlm)
                      ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) + &
                           cg(icg)*qm(i1,i2,l1,l2,lm)*(gpot0(mlm)-gpotb(mlm))
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ppiall
  subroutine ppildau(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2, &
       lmaxu,vumm,nlx1,nlx2,ppiz,isp,idu) ! ppiz from non-local part of pot (LDA+U)
    implicit none
    integer :: lmx1,lmx2,nf1,nf1s,nf2,nf2s,nlx1,nlx2,lmaxu, lx1(nf1),lx2(nf2),isp,idu(4),&
         i1,i2,ilm1,ilm2,l1,l2,m1,m2
    real(8):: &
         v1(0:lmx1,nf1),d1(0:lmx1,nf1), &
         v2(0:lmx2,nf2),d2(0:lmx2,nf2),vd1(2),vd2(2)
    complex(8):: vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,3,3,2,0:lmaxu),add,ppiz(nf1,nf2,nlx1,nlx2)
    do  i1 = 1, nf1
       do  i2 = 1, nf2 ! ... Matrix elements of vumm constructed from (u,s)
          ilm1 = 0
          do  l1 = 0, min(lx1(i1),lmaxu)
             do  m1 = -l1, l1
                ilm1 = ilm1+1
                ilm2 = 0
                do  l2 = 0, min(lx2(i2),lmaxu)
                   do  m2 = -l2, l2
                      ilm2 = ilm2+1
                      if (idu(l1+1) /= 0 .AND. l1 == l2) then
                         if (i1 <= nf1s .AND. i2 <= nf2s) then!.. (u,s)V(u,s)
                            vd1=[v1(l1,i1),d1(l1,i1)]
                            vd2=[v2(l2,i2),d2(l2,i2)]
                            add = sum(vd1*matmul(vumm(m1,m2,1:2,1:2,isp,l1),vd2))
                         elseif (i1 > nf1s .AND. i2 > nf2s) then! ... zVz
                            add = vumm(m1,m2,3,3,isp,l1) 
                         elseif (i1 > nf1s) then!           ... zV(u,s)
                            add = sum( vumm(m1,m2,3,1:2,isp,l1)*[v2(l2,i2),d2(l2,i2)])
                         elseif (i2 > nf2s) then!           ... (u,s)Vz
                            add = sum([v1(l1,i1),d1(l1,i1)]* vumm(m1,m2,1:2,3,isp,l1))
                         endif
                         ppiz(i1,i2,ilm1,ilm2) = ppiz(i1,i2,ilm1,ilm2) + add
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine ppildau

  subroutine mksomat(l,m1,m2,isp,somatpm) ! L-(for isp=1),L+(for isp=2) matrix
    integer:: l,m1,m2,isp
    real(8):: a1,a2
    complex(8):: somatpm,img=(0d0,1d0)
    somatpm=0d0
    a1 = dsqrt(dble((l-abs(m2))*(l+abs(m2)+1)))
    a2 = dsqrt(dble((l+abs(m2))*(l-abs(m2)+1)))
    if (isp == 1) then ! <l,m1|L-|l,m2> for spin (1,2) block
       if (abs(m2) > 1 .AND. (abs(m2)+1) <= l) then !             Case A
          if (m2 > 0) then
             if ((abs(m2)-1) == m1) somatpm=     (-1)**(2*m2-1)*a2*0.5d0
             if ((abs(m2)-1) ==-m1) somatpm= img*(-1)**(2*m2-1)*a2*0.5d0
             if ((abs(m2)+1) == m1) somatpm=     a1*0.5d0
             if ((abs(m2)+1) ==-m1) somatpm=-img*a1*0.5d0
          else
             if ((abs(m2)-1) == m1) somatpm=-img*(-1)**(2*m2-1)*a2*0.5d0
             if ((abs(m2)-1) ==-m1) somatpm=     (-1)**(2*m2-1)*a2*0.5d0
             if ((abs(m2)+1) == m1) somatpm= img*a1*0.5d0
             if ((abs(m2)+1) ==-m1) somatpm=     a1*0.5d0
          endif
       elseif (abs(m2) > 1 .AND. (abs(m2)+1) > l) then !             Case B
          if (m2 > 0) then
             if ((abs(m2)-1) == m1) somatpm =      (-1)**(2*m2+1)*a2*0.5d0
             if ((abs(m2)-1) ==-m1) somatpm=  img*(-1)**(2*m2+1)*a2*0.5d0
          endif
          if (m2 < 0) then
             if ((abs(m2)-1) == m1) somatpm =  -img*(-1)**(2*m2-1)*a2*0.5d0
             if ((abs(m2)-1) ==-m1) somatpm=     (-1)**(2*m2-1)*a2*0.5d0
          endif
       elseif (abs(m2) == 1 .AND. (abs(m2)+1) <= l) then!             Case C
          if (m2 > 0) then
             if (m1 == 0) somatpm =       (-1)**m2*a2*dsqrt(0.5d0)
             if ((abs(m2)+1) == m1) somatpm =      a1*0.5d0
             if ((abs(m2)+1) == -m1) somatpm= -img*a1*0.5d0
          endif
          if (m2 < 0) then
             if (m1 == 0) somatpm =           -img*(-1)**m2*a2*dsqrt(0.5d0)
             if ((abs(m2)+1) == m1) somatpm =  img*a1*0.5d0
             if ((abs(m2)+1) == -m1) somatpm=     a1*0.5d0
          endif
       elseif (abs(m2) == 1 .AND. (abs(m2)+1) > l) then!             Case D
          if (m2 > 0) then
             if (m1 == 0) somatpm =       (-1)**m2*a2*dsqrt(0.5d0)
          endif
          if (m2 < 0) then
             if (m1 == 0) somatpm =  -img*(-1)**m2*a2*dsqrt(0.5d0)
          endif
       elseif (abs(m2) == 0) then!             Case m=0
          if (m1 == 1)  somatpm =       a1*dsqrt(0.5d0)
          if (m1 == -1) somatpm =    -img*a1*dsqrt(0.5d0)
       else
          call rx('mksomat this can not occur 111')
       endif
    else! <l,m1|L+|l,m2> for spin (2,1) block
       if (abs(m2) > 1 .AND. (abs(m2)+1) <= l) then!               Case A
          if (m2 > 0) then
             if ((abs(m2)-1) == m1) somatpm =   a2*0.5d0
             if ((abs(m2)-1) == -m1) somatpm =  -img*a2*0.5d0
             if ((abs(m2)+1) == m1) somatpm =  (-1)**(2*m2+1)*a1*0.5d0
             if ((abs(m2)+1) == -m1) somatpm =   img*(-1)**(2*m2+1)*a1*0.5d0
          else
             if ((abs(m2)-1) == m1) somatpm =  img*a2*0.5d0
             if ((abs(m2)-1) == -m1) somatpm =  a2*0.5d0
             if ((abs(m2)+1) == m1) somatpm =  -img*a1*(-1)**(2*m2+1)*0.5d0
             if ((abs(m2)+1) == -m1) somatpm =      (-1)**(2*m2+1)*a1*0.5d0
          endif
       elseif (abs(m2) > 1 .AND. (abs(m2)+1) > l) then !               Case B
          if (m2 > 0) then
             if ((abs(m2)-1) == m1) somatpm =  a2*0.5d0
             if ((abs(m2)-1) == -m1) somatpm=     -img*a2*0.5d0
          endif
          if (m2 < 0) then
             if ((abs(m2)-1) == m1) somatpm =  img*a2*0.5d0
             if ((abs(m2)-1) == -m1) somatpm= a2*0.5d0
          endif
       elseif (abs(m2) == 1 .AND. (abs(m2)+1) <= l) then!               Case C
          if (m2 > 0) then
             if (m1 == 0) somatpm =                a2*dsqrt(0.5d0)
             if ((abs(m2)+1) == m1) somatpm =   (-1)**(2*m2+1)*a1*0.5d0
             if ((abs(m2)+1) == -m1) somatpm=  img*(-1)**(2*m2+1)*a1*0.5d0
          endif
          if (m2 < 0) then
             if (m1 == 0) somatpm =              img*a2*dsqrt(0.5d0)
             if ((abs(m2)+1) == m1) somatpm =   -img*(-1)**(2*m2+1)*a1*0.5d0
             if ((abs(m2)+1) == -m1) somatpm =       (-1)**(2*m2+1)*a1*0.5d0
          endif
       elseif (abs(m2) == 1 .AND. (abs(m2)+1) > l) then!               Case D
          if (m2 > 0) then
             if (m1 == 0) somatpm =         a2*dsqrt(0.5d0)
          endif
          if (m2 < 0) then
             if (m1 == 0) somatpm =     img*a2*dsqrt(0.5d0)
          endif
       elseif (abs(m2) == 0) then !               Case m=0
          if (m1 == 1) somatpm =  -a1*dsqrt(0.5d0)
          if (m1 == -1) somatpm=  -img*a1*dsqrt(0.5d0)
       else
          call rx('mksomat this can not occur 222')
       endif
    end if
  end subroutine mksomat
end module m_gaugm

    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nf1   :number of functions of first kind for each l
    !i   lmx1  :dimensions f1
    !i   lx1   :l-cutoffs for each of the nf1 functions
    !i   nf2   :number of functions of second kind for each l
    !i   lmx2  :dimensions f2
    !i   lx2   :l-cutoffs for each of the nf2 functions
    !i   lmxl  :l-cutoff for density, potential on the radial mesh
    !i   nlml  :(lmxl+1)*(lmxl+1)
    !i   cg    :Clebsch Gordon coefficients
    !i   jcg   :L q.n. for the C.G. coefficients
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   qm    :Moments of the products f1~*f2~ - f1*f2
    !i         :For local orbitals, the term f1*f2 is absent
    !i   gpotb :integrals of compensating gaussians * local smooth estat
    !i         :pot calculated from the compensated smoothed local density
    !i   gpot0 :integrals of local gaussians * phi0~ (smves.f)
    !i         :phi0~ is the estatic potential of the interstitial
    !i         :smooth density including compensating local charges.
    !i   lmux  :l-cutoff for sigma,tau, and spherical part of ppi
    !i   ppi0  :contribution to ppi from spherical part of potential
    !i   nlx1  :dimensions ppi
    !i   nlx2  :dimensions ppi
    !i   ppi   :nonspherical parts only of potential integrals (paug1,paug2)
    !o Outputs
    !o   ppi   :augmentation potential integrals assembled
    !r Remarks
    !r  The terms proportional to the spherical part of V are added to
    !r  ppi (generated in pvagm1, pvagm1c, pvagm2), and also the term
    !r     qm (gpot0-gpotb)
    !r  which are the terms proportional to Qkk'LL'M of Eqs. 28,29 in
    !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
    !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
    !r      ed. (Springer-Verlag, Berlin) 2000.
    !r   However, see comments in augmat.f about indices k,k'.
    !u Updates
    !u   14 Sep 01 Extended to local orbitals, which have no smooth part
    !u   17 May 00 Adapted from nfp paug1.f
    ! ----------------------------------------------------------------------

  ! subroutine paug2(nr,nlml,v2,rwgt,cg,jcg,indxcg, &
  !      nf1,nf1s,lmx1,lx1,f1,nf2,nf2s,lmx2,lx2,f2,ssum,nlx1,nlx2,ppi)
  !   !- Put in ppi constribution from smooth pot, smooth wave functions
  !   ! ----------------------------------------------------------------------
  !   !i Inputs
  !   !i   nr    :number of radial mesh points
  !   !i   nlml  :L-cutoff for charge density on radial mesh
  !   !i   v2    :smooth potential, seen by unaugmented functions
  !   !i   rwgt  :radial mesh weights
  !   !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !   !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !   !i   indxcg:index for Clebsch Gordon coefficients
  !   !i   nf1   :number of functions of first kind for each l
  !   !i   nf1s  :number of functions of first kind for each l, for which
  !   !i         :there is a smooth part to be subtracted
  !   !i   lmx1  :dimensions f1,sum
  !   !i   lx1   :l-cutoffs for each of the nf1 functions
  !   !i   f1    :`bra' radial functions tabulated numerically; see Outputs
  !   !i   nf2   :number of functions of second kind for each l
  !   !i   nf2s  :number of functions of second kind for each l, for which
  !   !i         :a smooth part is to be subtracted.
  !   !i   lmx1  :dimensions f2,sum
  !   !i   lx2   :l-cutoffs for each of the nf2 functions
  !   !i   f2    :`ket' radial functions tabulated numerically; see Outputs
  !   !i   sum   :work array holding integrations
  !   !i   nlx1  :dimensions ppi
  !   !i   nlx2  :dimensions ppi
  !   !o Outputs
  !   !o   ppi   :<f1^ | V2~ | f2^> subtracted for each f1,f2 pair
  !   !r Remarks
  !   !r    Makes the 2nd half of the first term in Eq. 29, Springer book.
  !   !r    But see Remarks in augmat.f:  there are three flavors of this
  !   !r    contribution to pi:
  !   !r         P V2~ P      H V2~ P      H V2~ H
  !   !r    where V2~ is the one-center repsn'f of the smooth potential.
  !   !u Updates
  !   !u   24 Aug 01 Extended to local orbitals, which have no smooth part
  !   ! ----------------------------------------------------------------------
  !   implicit none
  !   integer :: lmx1,lmx2,nf1,nf2,nf1s,nf2s,nlml,nlx1,nlx2,nr
  !   integer :: lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
  !   double precision :: v2(nr,nlml),rwgt(nr), &
  !        ppi(nf1,nf2,nlx1,nlx2),f1(nr,0:lmx1,nf1s),f2(nr,0:lmx2,nf2s), &
  !        cg(1),ssum(0:lmx1,0:lmx2,nlml)
  !   integer :: i1,i2,icg,ilm1,ilm2,ix,l1,l2,ll,mlm,nlm1,nlm2,i
  !   double precision :: sam
  !   ppi=0d0 
  !   do  i1 = 1, nf1s ! Sum over CG coefficients, make radial integrals as needed
  !      do  i2 = 1, nf2s
  !         nlm1 = (lx1(i1)+1)**2
  !         nlm2 = (lx2(i2)+1)**2
  !         ssum(0:lx1(i1),0:lx2(i2),1:nlml) = 2d10 ! --- integral not yet calculated
  !         do  ilm1 = 1, nlm1
  !            l1 = ll(ilm1)
  !            do  ilm2 = 1, nlm2
  !               l2 = ll(ilm2)
  !               ix = max0(ilm1,ilm2)
  !               ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
  !               do icg = indxcg(ix),indxcg(ix+1)-1
  !                  mlm = jcg(icg)
  !                  if (mlm > 1 .AND. mlm <= nlml) then
  !                     if (ssum(l1,l2,mlm) > 1d10) then
  !                        ssum(l1,l2,mlm) = sum([(rwgt(i)*v2(i,mlm)*f1(i,l1,i1)*f2(i,l2,i2),i=2,nr)])
  !                     endif
  !                     ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) - cg(icg)*ssum(l1,l2,mlm)
  !                  endif
  !               enddo
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  ! end subroutine paug2
  ! subroutine paug1(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2, &
  !      lmxa,nlml,cg,jcg,indxcg,vum,nlx1,nlx2,ppi)
  !   !- Add to ppi constribution from true pot, true wave functions
  !   ! ----------------------------------------------------------------------
  !   !i Inputs
  !   !i   nf1   :number of functions of first kind for each l
  !   !i   nf1s  :number of functions of first kind for each l, for which
  !   !i         :there is a smooth part to be subtracted (which also
  !   !i         :corresponds to the functions which connect to envelope
  !   !i         :functions)
  !   !i   lmx1  :dimensions f1
  !   !i   lx1   :l-cutoffs for each of the nf1 functions
  !   !i   v1    :values of f1 at rofi(nr)
  !   !i   d1    :slopes of f1 at rofi(nr)
  !   !i   nf2   :number of functions of second kind for each l
  !   !i   nf2s  :number of functions of second kind for each l, for which
  !   !i         :a smooth part is to be subtracted (which also
  !   !i         :corresponds to the functions which connect to envelope
  !   !i         :functions)
  !   !i   lmx2  :dimensions f2
  !   !i   lx2   :l-cutoffs for each of the nf2 functions
  !   !i   v1    :values of f1 at rofi(nr)
  !   !i   d1    :slopes of f1 at rofi(nr)
  !   !i   lmxa  :augmentation l-cutoff
  !   !i   nlml  :L-cutoff for charge density on radial mesh
  !   !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !   !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !   !i   indxcg:index for Clebsch Gordon coefficients
  !   !i   vum   :integrals ((u,s,gz) * v1 * (u,s,gz)),
  !   !i         :decomposed by M where v1 = sum_M v1_M Y_M
  !   !i         :vum(l1,l2,M,i,j) = (ui_l1 v1_M uj_l2)
  !   !i   nlx1  :dimensions ppi
  !   !i   nlx2  :dimensions ppi
  !   !o Outputs
  !   !o   ppi   :partial matrix element of potential; see Remarks
  !   !r Remarks
  !   !r    Makes the first half of the first term in Eq. 29, Springer book.
  !   !r    But see Remarks in augmat.f:  there are three flavors of this
  !   !r    contribution to pi:
  !   !r         P~ V1 P~      H~ V1 P~      H~ V1 H~
  !   !r    where V1 is true potential.
  !   !u Updates
  !   !u   27 Aug 01 Extended to local orbitals.  Altered argument list.
  !   ! ----------------------------------------------------------------------
  !   implicit none
  !   integer :: lmx1,lmx2,lmxa,nf1,nf1s,nf2,nf2s,nlml,nlx1,nlx2
  !   integer :: lx1(nf1),lx2(nf2),jcg(1),indxcg(1)
  !   double precision :: vum(0:lmxa,0:lmxa,nlml,3,3), &
  !        v1(0:lmx1,nf1),d1(0:lmx1,nf1), &
  !        v2(0:lmx2,nf2),d2(0:lmx2,nf2), ppi(nf1,nf2,nlx1,nlx2),cg(1),vd1(2),vd2(2)
  !   integer :: i1,i2,icg,ilm1,ilm2,ix,l1,l2,ll,mlm,nlm1,nlm2
  !   double precision :: add=1d99
  !   ! ... Combine with CG coefficents
  !   do  i1 = 1, nf1s
  !      do  i2 = 1, nf2s
  !         nlm1 = (lx1(i1)+1)**2
  !         nlm2 = (lx2(i2)+1)**2
  !         do  ilm1 = 1, nlm1
  !            l1 = ll(ilm1)
  !            do  ilm2 = 1, nlm2
  !               l2 = ll(ilm2)
  !               ix = max0(ilm1,ilm2)
  !               ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
  !               do  icg = indxcg(ix),indxcg(ix+1)-1
  !                  mlm = jcg(icg)
  !                  if (mlm > 1 .AND. mlm <= nlml) then
  !                     vd1=[v1(l1,i1),d1(l1,i1)]
  !                     vd2=[v2(l2,i2),d2(l2,i2)]
  !                     ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) &
  !                          + cg(icg)*sum(vd1*matmul(vum(l1,l2,mlm,1:2,1:2),vd2))
  !                  endif
  !               enddo
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   ! --- Matrix elements involving local orbitals ---
  !   if (nf1s >= nf1 .AND. nf2s >= nf2) return
  !   do  i1 = 1, nf1
  !      do  i2 = 1, nf2
  !         if (i1 > nf1s .OR. i2 > nf2s) then
  !            nlm1 = (lx1(i1)+1)**2
  !            nlm2 = (lx2(i2)+1)**2
  !            do  ilm1 = 1, nlm1
  !               l1 = ll(ilm1)
  !               do  ilm2 = 1, nlm2
  !                  l2 = ll(ilm2)
  !                  ix = max0(ilm1,ilm2)
  !                  ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
  !                  do  icg = indxcg(ix),indxcg(ix+1)-1
  !                     mlm = jcg(icg)
  !                     if (mlm > 1 .AND. mlm <= nlml) then
  !                        if (i1 > nf1s .AND. i2 > nf2s) then!    <g | V | g>
  !                           add = vum(l1,l2,mlm,3,3)
  !                        elseif (i1 > nf1s) then  !              <g | V | (u,s)>
  !                           add = sum( vum(l1,l2,mlm,3,1:2) * [v2(l2,i2),d2(l2,i2)])
  !                        elseif (i2 > nf2s) then !               <(u,s) | V | g>
  !                           add = sum([v1(l1,i1),d1(l1,i1)] * vum(l1,l2,mlm,1:2,3))
  !                        endif
  !                        ppi(i1,i2,ilm1,ilm2) = ppi(i1,i2,ilm1,ilm2) + cg(icg)*add
  !                     endif
  !                  enddo
  !               enddo
  !            enddo
  !         endif
  !      enddo
  !   enddo
  ! end subroutine paug1
!   subroutine pvaglc(&
!        nf1,nf1s,lmx1,lx1,nlx1,v1,d1, &
!        nf2,nf2s,lmx2,lx2,nlx2,v2,d2,&
!        isp,lso, &
!        hab,vab,sab,sodb,sondb, &
!        lmux,sig,tau,ppi,hsozz,hsopm)
!     !- Augmentation part of sig and tau and ppi, loc. orbitals (spher. pot)
!     ! ----------------------------------------------------------------------
!     !i Inputs
!     !i   nf1   :number of 'bra' function types for each l
!     !i   nf1s  :number of 'bra' function types for each l
!     !i         :which are formed by linear combination of (u,s)
!     !i   lmx1  :dimensions f1
!     !i   lx1   :l-cutoffs for each of the nf1 functions
!     !i   nlx1  :dimensions hsozz,hsopm
!     !i   v1    :values of f1 at rofi(nr) (not multiplied by r)
!     !i   d1    :slopes of f1 at rofi(nr) (not multiplied by r)
!     !i   isp   :(used only when adding to SO hamiltonian)
!     !i         :isp=1 => make <l,m|L-|l,m'> of hsopm
!     !i         :isp=2 => make <l,m|L+|l,m'> of hsopm
!     !i   nf2   :number of 'ket' function types for each l
!     !i   nf2s  :number of 'ket' function types for each l
!     !i         :which are formed by linear combination of (u,s)
!     !i   lmx2  :dimensions f2
!     !i   lx2   :l-cutoffs for each of the nf2 functions
!     !i   nlx2  :dimensions hsozz,hsopm
!     !i   v2    :values of f2 at rofi(nr) (not multiplied by r)
!     !i   d2    :slopes of f2 at rofi(nr) (not multiplied by r)
!     !i   lso   :1 include L.S coupling; 2 include LzSz part only
!     !i   hab   :hamiltonian matrix elements of (u,s); see Remarks
!     !i         :Leading dimension of hab corresponds to these entries:
!     !i         :(1)=<ul|h|ul>   (5)=<ul|h|gz>   (9)=<gz|h|sl>
!     !i         :(2)=<ul|h|sl>   (6)=<sl|h|gz>
!     !i         :(3)=<sl|h|ul>   (7)=<gz|h|gz>
!     !i         :(4)=<sl|h|sl>   (8)=<gz|h|ul>
!     !i   vab   :potential matrix elements of (u,s)
!     !i   sab   :overlap matrix elements of (u,s)
!     !i   sodb
!     !i   sondb
!     !i   lmux  :l-cutoff for sig,tau,ppi
!     !i         :Usually min (largest lmx1, largest lmx2)
!     ! o Inputs/Outputs
!     ! o  sig   :On input, sig is smoothed counterpart of overlap
!     ! o        :  matrix <f1^ f2^> (see pvagm2)
!     ! o        :On output, sig is overwritten by
!     ! o        :  <f1~ f2~> - <f1^ f2^>
!     ! o        :  where the first term is calculated from sab
!     ! o  tau   :On input, tau is smoothed counterpart of k.e.
!     ! o        :  matrix <f1^ -nabla f2^> (see pvagm2)
!     ! o        :On output, tau is overwritten by
!     ! o        :  <f1~ -nabla f2~> - <f1^ -nabla f2^>
!     ! o        :  where the first term is calculated from hab-vab
!     ! o  ppi   :On input, ppi is smoothed counterpart of pot
!     ! o        :  matrix <f1^ v2~ f2^> (see pvagm2)
!     ! o        :On output, ppi is overwritten by
!     ! o        :  <f1~ v1(l=0) f2~> - <f1^ v2~(l=0) f2^>
!     ! o        :  where the first term is calculated from vab
!     ! o  hsozz :
!     ! o  hsopm :
!     !r Remarks
!     !r   See Remarks for routine gaugm, above.
!     !r
!     !r   Symmetrization of kinetic energy.  The kinetic energy matrix
!     !r   should always be symmetric.  In general the
!     !r   difference T_fg - T_gf between orbitals f and g is
!     !r    T_fg - T_gf = -[r^2 (fg' - f'g)]_rmax = -W{f,g;rmax}
!     !r   which is nonzero.  However, the kinetic energy matrix always
!     !r   takes of two forms:
!     !r      1) it consists of a difference T_fg = <f T g> - <f^ T g^>
!     !r         between true and smoothed functions where f matches f^
!     !r         and g matches g^ at rmax.  In this case, T_fg - T_gf = 0
!     !r         because W{f,g;rmax} = W{f^,g^;rmax}
!     !r      2) either f or g is a local orbital with value,slope=0.
!     !r         Then -W{f,g;rmax} = 0.
!     !r
!     !r   See also Remarks for routine gaugm, above.
!     !u Updates
!     !u   26 Sep 05 (A. Chantis) added matrix elements for SO coupling
!     ! ----------------------------------------------------------------------
!     implicit none
!     integer :: lmux,lmx1,lmx2,nf1,nf2,nf1s,nf2s,lx1(nf1),lx2(nf2),isp,nlx1,nlx2
!     double precision :: hab(3,3,0:n0-1),sab(3,3,0:n0-1),vab(3,3,0:n0-1), &
!          v1(0:lmx1,nf1),d1(0:lmx1,nf1), &
!          v2(0:lmx2,nf2),d2(0:lmx2,nf2), &
!          sig(nf1,nf2,0:lmux),tau(nf1,nf2,0:lmux),ppi(nf1,nf2,0:lmux)
!     !     Spin-Orbit related
!     double precision :: sodb(3,3,0:n0-1),sondb(3,3,0:n0-1)
!     complex(8):: hsozz(nf1,nf2,nlx1,nlx2),hsopm(nf1,nf2,nlx1,nlx2),img=(0d0,1d0),somatzz,somatpm
!     integer :: i1,i2,lmax1,lmax2,lmax,l
!     integer :: m1,m2,l1,l2,lso
!     real(8):: tmp(nf1,nf2,nlx1,nlx2),tmp1(nf1,nf2,nlx1,nlx2),a1,a2
! !    print *,'nnnnnnffffffff222',nf1,nf2,nf1s,nf2s
!     do  i1 = 1, nf1
!        do  i2 = 1, nf2
!           if (i1 > nf1s .OR. i2 > nf2s) then
!              lmax1 = lx1(i1)
!              lmax2 = lx2(i2)
!              lmax = min0(lmax1,lmax2)
!              l1 = 0
!              l2 = 0
!              do  l = 0, lmax
!                 !       ... hzz,vzz,szz
!                 if (i1 > nf1s .AND. i2 > nf2s) then
!                    sig(i1,i2,l) = -sig(i1,i2,l) + sab(3,3,l)
!                    tau(i1,i2,l) = -tau(i1,i2,l) + (hab(3,3,l)-vab(3,3,l))
!                    ppi(i1,i2,l) = -ppi(i1,i2,l) + vab(3,3,l)
!                    !       ... hzu,vzu,szs
!                 elseif (i1 > nf1s) then
!                    sig(i1,i2,l)= -sig(i1,i2,l)  + sum(sab(3,1:2,l)*[v2(l,i2),d2(l,i2)])
!                    tau(i1,i2,l)= -tau(i1,i2,l)  + sum((hab(3,1:2,l)-vab(3,1:2,l))*[v2(l,i2),d2(l,i2)])
!                    ppi(i1,i2,l) = -ppi(i1,i2,l) + sum(vab(3,1:2,l)*[v2(l,i2),d2(l,i2)])
!                    !       ... huz,vuz,ssz
!                 elseif (i2 > nf2s) then
!                    sig(i1,i2,l) = -sig(i1,i2,l) + sum([v1(l,i1),d1(l,i1)]*sab(1:2,3,l))
!                    tau(i1,i2,l) = -tau(i1,i2,l) + sum([v1(l,i1),d1(l,i1)]*(hab(1:2,3,l)-vab(1:2,3,l)))
!                    ppi(i1,i2,l) = -ppi(i1,i2,l) + sum([v1(l,i1),d1(l,i1)]*vab(1:2,3,l))
!                 endif
!              enddo
!           endif
!        enddo
!     enddo
!     ! if (lso == 0) return
!     ! if (lso /= 0) then
!     !    tmp=0d0  !call dpzero(tmp,   nf1*nf2*nlx1*nlx2)
!     !    tmp1=0d0 !call dpzero(tmp1,  nf1*nf2*nlx1*nlx2)
!     ! endif
!     ! !       ... Spin-Orbit matrix elements in real harmonics.
!     ! i1loop: do  i1 = 1, nf1
!     !    i2loop: do  i2 = 1, nf2
!     !       if (i1 > nf1s .OR. i2 > nf2s) then
!     !          continue
!     !       else
!     !          cycle
!     !       endif
!     !       lmax1 = lx1(i1)
!     !       lmax2 = lx2(i2)
!     !       lmax = min0(lmax1,lmax2)
!     !       l1 = 0
!     !       l2 = 0
!     !       lloop: do  l = 0, lmax
!     !          mloop: do  m1 = -l, l
!     !             l1 = l1 + 1
!     !             if (m1 >= (-l+1)) l2 = l2 - (2*l + 1)
!     !             do  m2 = -l, l
!     !                l2 = l2 + 1
                   
!     !                if (i1 > nf1s .AND. i2 > nf2s) then !     ... hso_zz
!     !                   tmp1(i1,i2,l1,l2) = sodb(3,3,l)
!     !                elseif (i1 > nf1s) then                      !             ... hso_zu
!     !                   tmp1(i1,i2,l1,l2) = sum(sodb(3,1:2,l)*[v2(l,i2),d2(l,i2)])
!     !                elseif (i2 > nf2s) then                      !             ... hso_uz
!     !                   tmp1(i1,i2,l1,l2) = sum([v1(l,i1),d1(l,i1)]*sodb(1:2,3,l))
!     !                endif
!     !                if ( m1+m2==0.AND.m2/=0.and.l1/=l2) then ! ... LzSz part
!     !                   if (l1 < l2 .AND. m1 == -m2) somatzz =  abs(m1)*img 
!     !                   if (l1 > l2 .AND. m1 == -m2) somatzz = -abs(m1)*img
!     !                   hsozz(i1,i2,l1,l2) = somatzz * tmp1(i1,i2,l1,l2)
!     !                endif
                   
!     !                if (lso == 1) then ! LxSx+LySy part
!     !                   if (i1 > nf1s .AND. i2 > nf2s) then
!     !                      tmp(i1,i2,l1,l2) = sondb(3,3,l)
!     !                   elseif (i1 > nf1s) then                         !         ... hso_zu
!     !                      tmp(i1,i2,l1,l2) = sum(sondb(3,1:2,l)*[v2(l,i2),d2(l,i2)])
!     !                   elseif (i2 > nf2s) then                         !         ... hso_uz
!     !                      tmp(i1,i2,l1,l2) = sum([v1(l,i1),d1(l,i1)]*sondb(1:2,3,l))
!     !                   endif
!     !                   call mksomat(l,m1,m2,isp,somatpm) ! LxSx+LySy lm part
!     !                   hsopm(i1,i2,l1,l2)= somatpm*tmp(i1,i2,l1,l2)
!     !                endif
!     !             enddo
!     !          enddo mloop
!     !       enddo lloop
!     !    enddo i2loop
!     ! enddo i1loop
!   end subroutine pvaglc
!   subroutine pvagm2(nf1,lmx1,lx1,f1,x1,nf2,lmx2,lx2,f2,x2, &
!        nr,rofi,rwgt,vsms,lmux,sig,tau,ppi)
!     !- Smooth part of sig, tau, ppi (spherical part of local smooth pot)
!     ! ----------------------------------------------------------------------
!     !i Inputs
!     !i   nf1   :number of 'bra' function types for each l
!     !i   lmx1  :dimensions f1
!     !i   lx1   :l-cutoffs for each of the nf1 functions
!     !i   f1    :'bra' unaugmented envelope fn (radial part * r).
!     !i         :Must be zero in channels that have no envelope functions.
!     !i         :See Remarks
!     !i   x1    :(radial derivative of 'bra' functions) * r = r*d(f1/r)/dr
!     !i   nf2   :number of 'ket' function types for each l
!     !i   lmx2  :dimensions f2
!     !i   lx2   :l-cutoffs for each of the nf2 functions
!     !i   f2    :'ket' unaugmented envelope fn (radial part . r).
!     !i         :Must be zero in channels that have no envelope functions
!     !i         :See Remarks.
!     !i   x2    :(radial derivative of 'ket' functions) * r = r*d(f2/r)/dr
!     !i   nr    :number of radial mesh points
!     !i   rofi  :radial mesh points
!     !i   rwgt  :radial mesh weights
!     !i   vsms  :spherical smooth potential V2~
!     !i   lmux  :maximum l for which to calculate sig,tau,ppi.
!     !i         :Usually min (largest lmx1, largest lmx2)
!     !o Outputs
!     !o   sig   :overlap matrix <f1 f2>, with f1,f2 on a radial mesh
!     !o   tau   :kinetic energy matrix <f1 -nabla f2>
!     !o         := <grad f1 grad f2> + l(l+1)< f1 r^-2 f2> + surface term
!     !o   ppi   :potential matrix integral <f1 vsms f2>, spherical part of vsms
!     !r Remarks
!     !r   This routine computes the matrix elements of smoothed functions
!     !r     <f1^ f2^>,  -<grad f1^ grad f2^>,   <f1^ (V2~_l=0) f2^>
!     !r   which correspond to the second term in Eq. 21 for overlap,
!     !    and the second half of the first term Eq. 29 of the
!     !r   Springer book chapter.  (pvagm2 computes only ppi matrix element
!     !r   for the spherical part of V2~).  Note that there are three
!     !r   flavors of sig,tau,ppi as described in the  Remarks in augmat.f:
!     !r        P op P     H op P      H op H
!     !r   with op = one of (1, -nabla, or V2~) and V2~ is the
!     !r   one-center repsn'f of the smooth potential.
!     !r   This routine makes one of these three; which one depends on
!     !r   the functions f1^,f2^ passed to pvagm2.
!     !r
!     !r   sig_kL,k'L' for k=1..nf1 and k'=1..nf2 is diagonal in LL' and
!     !r   depends only on l.  Only sig(nf1,nf2,0..l) is stored.  Ditto for
!     !r   tau and ppi (for spherical part of potential) treated here.
!     !r
!     !r   Formula for kinetic energy.  If f1=r*phi1, x1=r*phi1'  and
!     !r   f2=r*phi2,x2=r*phi2', the kinetic energy in the sphere to radius
!     !r   R for channel l, excluding angular part, is
!     !r     T_fg =  int_0^R (phi1) (-1/r d^2/dr^2) (r*phi2) r^2 dr
!     !r          = -int_0^R (r*phi1) (d^2/dr^2) (r*phi2) dr
!     !r          = -[(r phi1) d(r*phi2)/dr]_R
!     !r            +int_0^R d(r*phi1)/dr * d(r*phi2)/dr
!     !r          = -[r^2 phi1 dphi2/dr]_R + int_0^R r*(dphi1/dr)*r*(dphi2/dr)
!     !r          = -[f1*x2]_R + int_0^R r*x1 r*x2 dr
!     !r     The fourth step follows after some simple algebra.
!     !u Updates
!     !u   20 Jul 04 Added treatment for extended local orbitals.
!     !u             Envelopes f1,f2 must be zero for all channels that
!     !u             have no smooth counterparts to subtract.
!     ! ----------------------------------------------------------------------
!     implicit none
!     integer :: lmux,lmx1,lmx2,nf1,nf2,nr,lx1(nf1),lx2(nf2),ii
!     double precision :: rofi(nr),rwgt(nr),vsms(nr), &
!          f1(nr,0:lmx1,nf1),x1(nr,0:lmx1,nf1), &
!          f2(nr,0:lmx2,nf2),x2(nr,0:lmx2,nf2), &
!          ppi(nf1,nf2,0:lmux),sig(nf1,nf2,0:lmux),tau(nf1,nf2,0:lmux)
!     integer :: i1,i2,lmax1,lmax2,lmax,l,i
!     double precision :: ssum,sim,tum,vum,xbc
!     sig=0d0 !call dpzero(sig, nf1*nf2*(lmux+1))
!     tau=0d0 !call dpzero(tau, nf1*nf2*(lmux+1))
!     ppi=0d0 !call dpzero(ppi, nf1*nf2*(lmux+1))
!     do  i1 = 1, nf1
!        do  i2 = 1, nf2
!           do  l = 0, min0(lx1(i1),lx2(i2)) !lmax
!              sig(i1,i2,l) = sum([(rwgt(ii)*f1(ii,l,i1)*f2(ii,l,i2),ii=2,nr)])
!              ppi(i1,i2,l) = sum([(rwgt(ii)*f1(ii,l,i1)*f2(ii,l,i2)*vsms(ii),ii=2,nr)])
!              sim  =  sum([(rwgt(ii)*f1(ii,l,i1)*f2(ii,l,i2)/rofi(ii)**2,ii=2,nr)])
!              tum  =  sum([(rwgt(ii)*x1(ii,l,i1)*x2(ii,l,i2),ii=2,nr)])
!              xbc = f1(nr,l,i1) * x2(nr,l,i2)
!              tau(i1,i2,l) = tum + l*(l+1)*sim - xbc
!           enddo
!        enddo
!     enddo
!   end subroutine pvagm2

    !- Generic routine to make augmentation matrices
    !o Outputs
    !o   sig   :augmentation overlap integrals; see Remarks.
    !o   tau   :augmentation kinetic energy integrals
    !o   ppi   :augmentation potential integrals
    !o         :NB: tau is added to pi, so ppi = (kinetic energy + potential)
    !o         :ppi is returned complex (lcplxp=1)
    !o         :In the noncollinear case:
    !o         :ppi(:,:,:,:,isp,1,1) = spin-diagonal part of potential ! last argment=2 means impart or work area?, a little confusing.
    !o         :ppi(:,:,:,:,1,2,1) = up-down   (12) block of potential
    !o         :ppi(:,:,:,:,2,2,1) = down-down (21) block of potential
    !o         :In the SO case: the SO hamiltonian is added to ppi
    !o         :and also stored in ppi(:,:,:,:,:,:,2)
    !o    hsozz(:,:,:,:,isp) = SO diagonal  part of potential
    !o    hsopm(:,:,:,:,1) = SO L-S+ up-down   (12) block of potential
    !o    hsopm(:,:,:,:,2) = SO L+S- down-down (21) block of potential
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   rwgt  :radial mesh weights
    !i   lmxa  :augmentation L-cutoff
    !i   lmxl  :l-cutoff for density, potential on the radial mesh
    !i   nlml  :(lmxl+1)*(lmxl+1)
    !i   vsm   :smooth local potential
    !i   gpotb :integrals of compensating gaussians * local smooth estat
    !i         :pot calculated from the compensated smoothed local density
    !i   gpot0 :integrals of local gaussians * phi0~ (smves.f)
    !i         :phi0~ is the estatic potential of the interstitial
    !i         :smooth density including compensating local charges.
    !i   hab   :matrix elements of the hamiltonian for the true wave
    !i         :functions and true (spherical part of) potential,
    !i         :where w.f. are expressed in forms of (ul,sl,gz) as
    !i         :described in potpus.f
    !i   vab   :corresponding matrix elements of the (spherical part of) V
    !i   sab   :corresponding matrix elements of the overlap
    !i   sodb  :corresponding matrix elements of SO = 2/(c^2) dV/dr*(1/r)
    !i   qum   :moments (u,s) * (u,s) * r**l (momusl.f)
    !i   vum   :integrals ((u,s,gz) * (true potential) * (u,s,gz))
    !i         :for the full (nonspherical) true potential (momusl.f)
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg)
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   ...    Two sets of functions f1 and f2 are input, given as:
    !i   nf1   :number of 'bra' function types for each l
    !i   nf1s  :number of 'bra' function types for each l
    !i         :which are formed by linear combination of (u,s)
    !i   lmx1  :dimensions f1
    !i   lx1   :l-cutoffs for each of the nf1 functions
    !i   f1    :'bra' unaugmented envelope fn (radial part . r).
    !i         :Must be zero in channels that have no envelope functions.
    !i         :See Remarks
    !i   x1    :radial part of grad(f1/r) times r
    !i   v1    :values of f1 at rofi(nr) (not multiplied by r)
    !i   d1    :slopes of f1 at rofi(nr) (not multiplied by r)
    !i   nf2   :number of 'ket' function types for each l
    !i   nf2s  :number of 'ket' function types for each l
    !i         :which are formed by linear combination of (u,s)
    !i   lmx2  :dimensions f2
    !i   lx2   :l-cutoffs for each of the nf2 functions
    !i   f2    :'ket' unaugmented envelope fn (radial part . r).
    !i         :Must be zero in channels that have no envelope functions
    !i         :See Remarks.
    !i   x2    :radial part of grad(f2/r) times r
    !i   v2    :values of f2 at rofi(nr) (not multiplied by r)
    !i   d2    :slopes of f2 at rofi(nr) (not multiplied by r)
    !i   lmux  :l-cutoff for sigma,tau, and spherical part of ppi
    !i   nlx1  :dimensions ppi
    !i   nlx2  :dimensions ppi
    !i   lso   :if 2 calculate LzSz matrix elements, if 1 calculate LS
    !i  lcplxp=1 only now :0 if ppi is real; 1 if ppi is complex
    !i   ...   The following are associated with LDA+U
    !i   lmaxu :dimensioning parameter for U matrix
    !i   vumm  : orbital dependent pot matrix in (u,s) representation
    !i   lldau :lldau=0 => no U on this site otherwise
    !i         :U on this site
    !i   idu   :idu(l+1)=1 => this l has a nonlocal U matrix
    !l Local variables
    !l   ppi0  :contribution to ppi from spherical part of potential
    !l   qm    :multipole moments; see Remarks
    !xxxxl         :On exit, hso is added to ppi; see ppi
    !r Remarks
    !r   This subroutine assembles the various terms in the computation of
    !r   sigma, tau, pi that comprise the local (augmented) part of the
    !r   hamiltonian.  See also Remarks in routine augmat.f.  Equations
    !r   references are found in methods paper; see book by Springer:
    !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
    !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
    !r      ed. (Springer-Verlag, Berlin) 2000.
    !r
    !r   Augmentation integrals sig, tau, ppi are constructed as a
    !r   difference between augmented function products and unaugmented ones.
    !r
    !r     sig = integral f1~*f2~ - f1^*f2^.
    !r     tau = integral f1~ (-nabla) f2~ - f1^ (-nabla) f2^.
    !r         = integral -1* (grad f1~ grad f2~ - grad f1^ grad f2^)
    !r           since  f1~ grad f2~ = f1^ grad f2^ at the MT boundary.
    !r           NB: Integration by parts also works for true local orbitals
    !r           (without smooth parts ^) because val, slo = 0 at rmt.
    !r     ppi = integral f1~ vtrue f2~ - f1^ vsm f2^.
    !r
    !r   f1^, f2^ are one-center expansions to the envelope function.
    !r
    !r   The augmented part of sig and tau and ppi (spherical potential only)
    !r   is obtained from linear combinations of sab, hab, and vab.
    !r   The smoothed part is calculated directly by numerical integration.
    !r
    !r   hab,vab,sab are matrix elements of h, spherical-v, and 1 of (u,s)
    !r   where u and s are linear combinations of and phi,phidot defined as:
    !r   u has val=1, slo=1 at rmax, s has val=0, slo=1.  Matrix elements
    !r   of the potential are divided into a spherical part and the
    !r   remainder: the (dominant) spherical term is computed with the
    !r   relativistic wave function including the small component, while
    !r   the residual nonspherical part does not.
    !r
    !r   There two additional contributions to ppi.  One involves integrals
    !r
    !r   QkkLL'M integral g_M (V0~-V2~), where the integrals are passed as
    !r   gpot0 and gpotb.   QkkLL'M (cf the Springer book chapter Eq 27)
    !r   is evaluated in getqm (QkkLL'M is called qm here).
    !r      QkkLL'M gpot0 corresponds to Eq. 28;
    !r      QkkLL'M gpotb corresponds to the last term in Eq. 29.
    !r
    !r   QkkLL'M  = integrals ((u,s) * 1 * (u,s)) for valence functions,
    !r   and corresponding integrals for local orbitals, is obtained from
    !r   input qum.
    !r
    !r   The second term involves matrix elements of v-v(l=0) which are
    !r   obtained in the augmented term from linear combinations of input
    !r      vum = integrals ((u,s) * vtrue * (u,s))
    !r
    !r   Notes on local orbitals:
    !r   They occur either as true local orbitals (val,slo=0 at rmt) or
    !r   as extended orbitals, where a tail is attached at the MT sphere.
    !r   In the latter case, a smooth part may be subtracted, to
    !r   be included later as part of matrix elements of the interstitial
    !r   potential.  Caller must ensure that f1,f2 for all orbitals with
    !r   no smooth part are zero; and for extended local orbitals where
    !r   a smooth part is to be subtracted, that they match continuously
    !r   and differentiably to the smooth functions f1,f2 at rmt.
    !r
    !r   Notes on Spin-Orbit Coupling
    !r
    !r   The (spin-diagonal) LzSz inside augmentation sphere
    !r   has matrix elements in real harmonics:
    !r    <l,m|Lz|l,m'> = +i|m| delta_m,-m if m>0
    !r                    -i|m| delta_m,-m if m<0
    !r   where delta the Kronecker delta. The Hamiltonian is ordered
    !r   -m,...,0,...,m.  The p block for example looks like:
    !r          m:   -1       0     1
    !r         -1     0       0   i<|so|>
    !r          0     0       0      0
    !r          1 -i<|so|>    0      0
    !r
    !r   The (spin-off-diagonal) LxSx+LySy inside augmentation sphere.
    !r   has these matrix elements:
    !r     Spin up-down block is 1/2<l,m|L-|l,m'>; L- = Lx - i Ly
    !r     Spin down-up block is 1/2<l,m|L+|l,m'>; L+ = Lx + i Ly
    !r     The p and blocks have the following form: <p,m|L-|p,m'> =
    !r               -1                    0               1
    !r      -1        0           -i*a1*<|so|>/sqrt(2)     0
    !r       0 i*a2*<|so|>/sqrt(2)         0       -a2*<|so|>/sqrt(2)
    !r       1        0              a1*<|so|>/sqrt(2)     0
    !r     where a1 = sqrt{(l-m)*(l+m+1)} and a2 = sqrt{(l+m)*(l-m+1)}
    !r     <p,m|L+|p,m'> is the hermitian conjugate of <p,m|L-|p,m'>.
    !r
    !r   Notes on LDA+U
    !r   Call to ppildau adds nonlocal potential part to ppi matrix
    !r
    !u Updates
    !u   27 Nov 05 LDA+U => complex potential
    !u   09 Nov 05 (wrl) Convert dmat to complex form
    !u   26 Sep 05 (A. Chantis) Bug fix: local orbitals in conjunction w/ SO
    !u   29 Jun 05 (MvS) redesign to keep hso separated
    !u   27 Apr 05 LDA+U (Lambrecht)
    !u   03 Feb 05 (A. Chantis) calculate matrix elements of L.S
    !u             inside augmentation sphere
    !u    1 Sep 04 Adapted to handle complex ppi; so folded into ppi
    !u   12 Aug 04 Added treatment for extended local orbitals.
    !u             Envelopes f1,f2 must be zero for all channels that
    !u             have no smooth counterparts to subtract.
    !u   20 Jul 04 bug fix in pvagm3
    !u   29 Jun 04 (A. Chantis) calculate matrix elements of LzSz
    !u             inside augmentation sphere
    !u   14 Sep 01 Extended to local orbitals.  Altered argument list.
    !u   17 May 00 Adapted from nfp gaugm.f
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nf1   :number of 'bra' function types for each l
    !i   nf1s  :number of 'bra' function types for each l
    !i         :which are formed by linear combination of (u,s)
    !i   lmx1  :dimensions f1
    !i   lx1   :l-cutoffs for each of the nf1 functions
    !i   f1    :'bra' unaugmented envelope fn (radial part . r).
    !i         :Must be zero in channels that have no envelope functions.
    !i         :See Remarks
    !i   v1    :values of f1 at rofi(nr) (not multiplied by r)
    !i   d1    :slopes of f1 at rofi(nr) (not multiplied by r)
    !i   nf2   :number of 'ket' function types for each l
    !i   nf2s  :number of 'ket' function types for each l
    !i         :which are formed by linear combination of (u,s)
    !i   lmx2  :dimensions f2
    !i   lx2   :l-cutoffs for each of the nf2 functions
    !i   f2    :'ket' unaugmented envelope fn (radial part . r).
    !i         :Must be zero in channels that have no envelope functions
    !i         :See Remarks.
    !i   v2    :values of f2 at rofi(nr) (not multiplied by r)
    !i   d2    :slopes of f2 at rofi(nr) (not multiplied by r)
    !i   isp   :current spin index; used for SO coupling
    !i   lso   :if nonzero calculate hsozz = LzSz matrix elements
    !i         :if 1 calculate also  hsopm
    !i   hab   :hamiltonian matrix elements of (u,s); see Remarks
    !i   vab   :potential matrix elements of (u,s)
    !i   sab   :overlap matrix elements of (u,s)
    !i   sodb  :corresponding matrix elements of SO = 2/(c^2) dV/dr*(1/r)
    !i         :see comments for routine gaugm
    !i   sondb :corresponding matrix elements of SO = 2/(c^2) dV/dr*(1/r)
    !i         :for off-diagonal part; see comments for routine gaugm
    !i   lmux  :l-cutoff for sig,tau,ppi
    !i         :Usually min (largest lmx1, largest lmx2)
    ! o Inputs/Outputs
    ! o  sig   :On input, sig is smoothed counterpart of overlap
    ! o        :  matrix <f1^ f2^> (see pvagm2)
    ! o        :On output, sig is overwritten by
    ! o        :  <f1~ f2~> - <f1^ f2^>
    ! o        :  where the first term is calculated from sab
    ! o  tau   :On input, tau is smoothed counterpart of k.e.
    ! o        :  matrix <f1^ -nabla f2^> (see pvagm2)
    ! o        :On output, tau is overwritten by
    ! o        :  <f1~ -nabla f2~> - <f1^ -nabla f2^>
    ! o        :  where the first term is calculated from hab-vab
    ! o  ppi   :On input, ppi is smoothed counterpart of pot
    ! o        :  matrix <f1^ v2~ f2^> (see pvagm2)
    ! o        :On output, ppi is overwritten by
    ! o        :  <f1~ v1(l=0) f2~> - <f1^ v2~(l=0) f2^>
    ! o        :  where the first term is calculated from vab
    !o Outputs
    !o  hsozz  :spin-orbit contribution to ppi, spin diagonal block
    !o         :(LzSz part)
    !o         :isp=1: makes up-up block
    !o         :isp=2: makes down-down block
    !o  hsopm  :spin-orbit contribution to ppi, spin off-diagonal block
    !o         :(LxSx + LySy part)
    !o         :isp=1: makes up-down block
    !o         :isp=2: makes down-up block
    !u Updates
    !u   03 Feb 05 (A. Chantis) calculate matrix elements of L.S
    !r Remarks
    !r   See Remarks for routine gaugm, above.
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nf1   :number of 'bra' function types for each l
    !i   nf1s  :number of 'bra' function types for each l
    !i         :which are formed by linear combination of (u,s)
    !i   lmx1  :dimensions f1
    !i   lx1   :l-cutoffs for each of the nf1 functions
    !i   f1    :'bra' unaugmented envelope fn (radial part . r).
    !i         :Must be zero in channels that have no envelope functions.
    !i         :See Remarks
    !i   v1    :values of f1 at rofi(nr) (not multiplied by r), for
    !i         :functions 1..nf1s
    !i   d1    :slopes of f1 at rofi(nr) (not multiplied by r), for
    !i         :functions 1..nf1s
    !i   nf2   :number of 'ket' function types for each l
    !i   nf2s  :number of 'ket' function types for each l
    !i         :which are formed by linear combination of (u,s)
    !i   lmx2  :dimensions f2
    !i   lx2   :l-cutoffs for each of the nf2 functions
    !i   f2    :'ket' unaugmented envelope fn (radial part . r).
    !i         :Must be zero in channels that have no envelope functions
    !i         :See Remarks.
    !i   v2    :values of f2 at rofi(nr) (not multiplied by r), for
    !i         :functions 1..nf2s
    !i   d2    :slopes of f2 at rofi(nr) (not multiplied by r), for
    !i         :functions 1..nf2s
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   rwgt  :radial mesh weights
    !i   lmxa  :augmentation l-cutoff
    !i   qum   :integrals (u,s) * (u,s) * r**l (see momusl)
    !i         :qum(l1,l2,l,i,j) = int (ui_l1 uj_l2) * r**l
    !i         :u and s are augmented functions defined as:
    !i         :u has val=1, slo=1 at rmax, s has val=0, slo=1
    !i   lmxl  :l-cutoff for density, potential on the radial mesh
    !o Outputs
    !o   qm    :integrals (f1~*f2~ - f1*f2) r**l
    !r Remarks
    !r   The qm correspond to Qkk'LL'M of Eq. 27 in
    !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
    !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
    !r      ed. (Springer-Verlag, Berlin) 2000.
    !r   However, see comments in augmat.f about indices k,k'.
    !r
    !r   f1~ and f2~ are not supplied explicitly, but their contribution
    !r   to qm is supplied through qum.
    !r
    !r   If f1~ (f2~) is a valence orbital, it is the particular linear
    !r   combination of u,s that has value v1 (v2) and slope d1 (d2)
    !r   at rofi(nr); see qum(*,1..3), i.e. f1~ = V1*u1 + D1*s1
    !r
    !r   If f1~ (f2~) is a local orbital, the function is given directly;
    !r   no linear combination of functions is needed.
    !r
    !r   Define the following:
    !r   Q12uu = qum(l1,l2,lm,1,1)
    !r   Q12ss = qum(l1,l2,lm,2,2)
    !r   Q12zz = qum(l1,l2,lm,3,3)
    !r   Q12us = qum(l1,l2,lm,1,2) = Q21su;   Q12su = qum(l1,l2,lm,2,1)
    !r   Q12uz = qum(l1,l2,lm,1,3) = Q21zu;   Q12zu = qum(l1,l2,lm,3,1)
    !r   Q12sz = qum(l1,l2,lm,2,3) = Q21zs;   Q12zs = qum(l1,l1,lm,3,2)
    !r
    !r   Then
    !r   f1~ f2~ = V1*V2*Q12uu + V1*D2*Q12us + D1*V2*Q12su + D1*D2*Q12ss
    !r             if f1~ and f2~ are both valence orbitals;
    !r   f1~ f2~ = V1*Q12uz + D1*Q12sz
    !r             if f1~ = valence, f2~ = local;
    !r   f1~ f2~ = V2*Q12zu + D2*Q12zs = V2*Q21uz + D2*Q21sz
    !r             if f1~ = local, f2~ = valence;
    !r   f1~ f2~ = Q12zz
    !r             if f1~ = local, f2~ = local.
    !r
    !r   In principle, we must distinguish types of extended local orbitals.
    !r   Some types have no corresponding smooth part, e.g. ones with no
    !r   extension outside the MT region, while others do.  This distinction
    !r   is automatically incorporated if f1 and f2 are zero in the channels
    !r   where there no smooth part is to be subtracted.
    !r
    !r   The smooth contribution is computed by numerical integration.
    !u Updates
    !u   20 Jul 04 Added treatment for extended local orbitals.
    !u             Envelopes f1,f2 must be zero for all channels that
    !u             have no smooth counterparts to subtract.
    !u   20 Jul 04 bug fix: improper indices in qum for local orbitals
    !u   14 Sep 01 Added treatment for local orbitals.  Altered argument list.
    ! ----------------------------------------------------------------------


    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nr    :number of radial mesh points
    !i   nlml  :L-cutoff for charge density on radial mesh
    !i   vsm    :smooth potential, seen by unaugmented functions
    !i   rwgt  :radial mesh weights
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   nf1   :number of functions of first kind for each l
    !i   nf1s  :number of functions of first kind for each l, for which
    !i         :there is a smooth part to be subtracted
    !i   lmx1  :dimensions f1,sum
    !i   lx1   :l-cutoffs for each of the nf1 functions
    !i   f1    :`bra' radial functions tabulated numerically; see Outputs
    !i   nf2   :number of functions of second kind for each l
    !i   nf2s  :number of functions of second kind for each l, for which
    !i         :a smooth part is to be subtracted.
    !i   lmx1  :dimensions f2,sum
    !i   lx2   :l-cutoffs for each of the nf2 functions
    !i   f2    :`ket' radial functions tabulated numerically; see Outputs
    !i   sum   :work array holding integrations
    !i   nlx1  :dimensions ppi
    !i   nlx2  :dimensions ppi
    !o Outputs
    !o   ppi   :<f1^ | Vsm~ | f2^> subtracted for each f1,f2 pair
    !r Remarks
    !r    Makes the 2nd half of the first term in Eq. 29, Springer book.
    !r    But see Remarks in augmat.f:  there are three flavors of this
    !r    contribution to pi:
    !r         P Vsm~ P      H Vsm~ P      H Vsm~ H
    !r    where Vsm~ is the one-center repsn'f of the smooth potential.
    !u Updates
    !u   24 Aug 01 Extended to local orbitals, which have no smooth part
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nf1   :number of functions of first kind for each l
    !i   nf1s  :number of functions of first kind for each l, for which
    !i         :there is a smooth part to be subtracted (which also
    !i         :corresponds to the functions which connect to envelope
    !i         :functions)
    !i   lmx1  :dimensions f1
    !i   lx1   :l-cutoffs for each of the nf1 functions
    !i   v1    :values of f1 at rofi(nr)
    !i   d1    :slopes of f1 at rofi(nr)
    !i   nf2   :number of functions of second kind for each l
    !i   nf2s  :number of functions of second kind for each l, for which
    !i         :a smooth part is to be subtracted (which also
    !i         :corresponds to the functions which connect to envelope
    !i         :functions)
    !i   lmx2  :dimensions f2
    !i   lx2   :l-cutoffs for each of the nf2 functions
    !i   v1    :values of f1 at rofi(nr)
    !i   d1    :slopes of f1 at rofi(nr)
    !i   lmaxu : used to dimension vumm
    !i   vumm  :matrix elements of non-local potential
    !i         :vumm(m1,m2,1) = <u_i| vorb(m1,m2) |u_j>
    !i   nlx1  :dimensions ppiz
    !i   nlx2  :dimensions ppiz
    !i   isp   : spin we're working on
    !i   idu   : idu(l)=1 => this l has a U
    !o Outputs
    !o   ppiz  :partial matrix element of potential; see Remarks
    !r Remarks
    !u Updates
    !u   09 Nov 05 Convert ppi to complex form
    !u   08 Jun 05 (MvS) extended to local orbitals
    !u   27 Apr 05 (Lambrecht) LDA+U
    ! ----------------------------------------------------------------------
