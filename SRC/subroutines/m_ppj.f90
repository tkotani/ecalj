module m_ppj !------------ Get ppj --------------------
  use m_lmfinit,only: lmxax,ndima,nsp,nbas,pos,alat
  use m_rdata1,only: nrc,aac,bbc
  use m_ftox
  use m_lgunit,only:stdo
  use m_hamindex0,only: iclass=>iclasst
  use m_rdata1,only:rdata1init,nradmx,nnc,nrad,nindx_r,lindx_r,iord,nvmax,nrc,mindx,&
       gval_n,gcore_n,aac,bbc,gval_orth,zzpi,nrmxe=>nrmx
  complex(8),allocatable,public:: ppj(:,:,:)
contains
  subroutine m_ppj_init()
    implicit none
    integer,allocatable:: m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)
    integer:: lxx,nl,nrx,ldim2,ifoc,ix,ixx,nn
    allocate(ppj(ndima,ndima,nsp),source=(0d0,0d0)) ! ppj: ovalap matrix within MT. ibb=0 is for Gramschmidt
    lxx= 2*lmxax
    nl =   lmxax+1
    nrx= maxval(nrc)
    ldim2 = ndima
    open(newunit=ifoc,file='@MNLA_CPHI')
    read(ifoc,*)
    allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
    do ix =1,ldim2
      read(ifoc,*) m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx !m,m,l,ibas index 
      if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
    enddo
    close(ifoc)
    write(stdo,ftox)' --- read @MNLA_CPHI'
    nn=maxval(n_indx)
    ppjmat: block
      use m_ll,only: ll
      use m_gint,only: gintpp
      real(8),allocatable:: rofi(:)
      real(8):: cg(nl**2,nl**2,(2*nl-1)**2),symope(3,3),cy((lxx+1)**2),yl((lxx+1)**2),absdq,absqg2,absqg,r2s,ylk,dq(3)
      real(8):: ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nsp,nbas), rprodx(nrx,0:lxx), phij(0:lxx),psij(0:lxx),rphiphi(nrx)
      real(8),parameter::  pi = 4d0*atan(1d0),fpi =    4d0*pi
      complex(8),parameter:: img=(0d0,1d0)
      complex(8)::phaseatom
      integer:: lm1,lm2,m1,m2,l1,l2,lm3,l3,ia1,ia2,ibas,ibas1,ibas2,ir,lx,ngrpx=1,irad1,irad2,irad,n1x,n2x,ic,isp
      cg=0d0
      symope=reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],shape=[3,3])
      call rotcg(nl-1,symope,ngrpx,cg) !CG coefficient
      dq     =  [1d-10,0d0,0d0]
      absdq  = sqrt(sum(dq**2))
      absqg2 = (2*pi/alat)**2 *sum(dq**2)
      absqg  = sqrt(absqg2)
      call sylmnc(cy,lxx)
      call sylm(dq/absdq,yl,lxx,r2s) !spherical factor Y(dq)
      ppbrd=0d0
      ibasloop0: do ibas = 1,nbas ! radial integral  ppbrd = <phi phi j_l>
        ic = iclass(ibas)
        allocate(rofi(nrc(ic)),source=[(bbc(ic)*(exp((ir-1)*aac(ic))-1d0), ir=1,nrc(ic))])
        rprodx = 0d0
        do ir =2,nrc(ic)
          call bessl(absqg2*rofi(ir)**2,lxx,phij,psij) ! phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
          rprodx(ir,0:lxx) = [(rofi(ir)* phij(lx)* (absqg*rofi(ir))**lx,lx=0,lxx)] != r \times j_l(|dq|r) !bessel function for the expansion of exp(i(q1-q2) r)
        enddo
        ispinloop00: do concurrent(isp = 1:nsp, irad1=1:nrad(ibas),irad2=1:nrad(ibas))
          l1= lindx_r(irad1,ibas); n1x= nindx_r(irad1,ibas); l2= lindx_r(irad2,ibas); n2x= nindx_r(irad2,ibas)
          rphiphi(1:nrc(ic)) = [0d0,(gval_orth(ir,l1,n1x,isp,ic)*gval_orth(ir,l2,n2x,isp,ic)/rofi(ir),ir=2,nrc(ic))] ! phi = u = r \phi
          do lx = abs(l1-l2), l1+l2
            ppbrd(l1,n1x,l2,n2x,lx, isp,ibas)= gintpp(rprodx(1,lx), rphiphi,aac(ic),bbc(ic),nrc(ic)) 
          enddo
        enddo ispinloop00
        deallocate(rofi)
      enddo ibasloop0
      ispinloop02: do concurrent(isp=1:nsp,ia1=1:ndima,ia2=1:ndima)
        ibas1= ibas_indx(ia1); ibas2= ibas_indx(ia2)
        if(ibas2/=ibas1) cycle
        l1=l_indx(ia1); m1=m_indx(ia1); n1x=n_indx(ia1); lm1= l1**2+l1+1+ m1
        l2=l_indx(ia2); m2=m_indx(ia2); n2x=n_indx(ia2); lm2= l2**2+l2+1+ m2 !+ nc_max(l2,ibas2)
        phaseatom = exp( img* 2d0*pi*sum(dq*pos(:,ibas1)) ) !correct?
        do lm3= (l1-l2)**2+1, (l1+l2+1)**2 ! l3 takes |l1-l2|,...l1+l2
          l3 = ll(lm3)
          ylk= cy(lm3)*yl(lm3)
          ppj(ia1,ia2,isp) = ppj(ia1,ia2,isp) + ppbrd(l1,n1x,l2,n2x,l3,isp,ibas1) *cg(lm1,lm2, lm3) * fpi * img**l3* phaseatom * ylk
          ! cg(lm1,lm2,lm3)= \int Y_lm3(\hat(r)) Y_lm2(\hat(r)) Y_lm1(\hat(r)) \frac{d \Omega}{4\pi}.
          ! This is based on inverse expansion. See the book of angular momentum book of Rose.Eq.3.8.
        enddo
      enddo ispinloop02
      write(stdo,ftox)' --- end of ppjmat',sum(abs(ppbrd)),sum(abs(ppj(:,:,:)))
    endblock ppjmat
  end subroutine m_ppj_init
end module m_ppj
