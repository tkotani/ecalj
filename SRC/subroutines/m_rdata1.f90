!>refine mesh for GW and index sets for sugw
module m_rdata1
  use m_nvfortran,only:findloc
  integer,public:: nradmx,nrmx,nnc
  integer,allocatable,public:: nrad(:),nindx_r(:,:),lindx_r(:,:), iord(:,:,:,:),nvmax(:,:), nrc(:),mindx(:)
  real(8),allocatable,public:: gval_n (:,:,:,:,:), gcore_n  (:,:,:,:), aac(:), bbc(:) 
  real(8),allocatable,public:: gval_orth (:,:,:,:,:),zzpi(:,:,:,:,:)
  public rdata1init
  private
contains
  subroutine rdata1init(ncores,ndima,ncoremx,konf0,gval,gcore)
    use m_lmfinit,only: lmxax,nbas,ispec,nsp,lmxa,rmt,nspec,nr,spec_a,nrmx0=>nrmx
    use m_hamindex0,only: nclass,iclass=>iclasst,nindx,lindx,ibasindx,nphimx
    use m_lgunit,only:stdo
    use m_nvfortran,only:findloc
    implicit none
    integer::nxx,lxx,ibxx,ix,irr,mmm,l,ic,iorb,lx,nx,ispecc(nclass),icore,n1,n2,l1,l2,ierr,i,i1,ir,irad1,irad2,nmax
    integer,allocatable:: nrnn(:),nncx(:,:)
    real(8),allocatable:: aann(:),eb(:), zzp(:,:,:,:,:)
    real(8)::delrset,aa_n,bb_n, bb(nspec),delr
    real(8),allocatable:: rofi_n(:) ,rofi(:), gx_in(:,:,:,:,:)
    real(8),allocatable :: ovv(:,:,:,:,:)
    complex(8),allocatable::  zzpx(:,:)
    integer:: nr_n,isp,ibas,i2,is,mx,nn,ndima,ncores(nspec),ncoremx,konf0(0:lmxax,nclass)
    real(8):: gval(nrmx0,0:lmxax,nphimx,nsp,nclass),gcore(nrmx0,ncoremx,nsp,nclass)
    nxx= -999; lxx= -999; ibxx=-999
    allocate(aac(nclass), bbc(nclass),nrc(nclass))
    allocate(nrad(nbas), nindx_r(ndima,nbas),lindx_r(ndima,nbas),mindx(ndima))
    nrad = 0
    do ix = 1,ndima
       if(ix>1 .AND. nxx == nindx(ix) .AND. lxx == lindx(ix) .AND. ibxx== ibasindx(ix)) then
          mmm=mmm+1
       else
          lxx = lindx(ix)
          ibxx= ibasindx(ix)
          nxx = nindx(ix)
          nrad(ibxx)  = nrad(ibxx)+1
          irr = nrad(ibxx)
          nindx_r(irr,ibxx) = nxx
          lindx_r(irr,ibxx) = lxx
          mmm = -lxx
       endif
       mindx(ix) = mmm
    enddo
    nradmx= maxval(nrad)
    allocate(nvmax(0:lmxax,nclass),source=0)
    do ix =1,ndima
       l  = lindx(ix)
       ic = iclass(ibasindx(ix))
       nvmax(l,ic) = max(nindx(ix),nvmax(l,ic))      !if( nxx> nvmax(l,ic) ) nvmax(l,ic) = nxx
    enddo
    nmax = maxval(nvmax)
    allocate( iord(-lmxax:lmxax,nmax,0:lmxax,nbas))
    iorb=0
    ix=0
    do ibas = 1,nbas
       ic = iclass(ibas)
       is = ispec(ibas)
       do lx = 0,lmxa(is)
          do nx = 1,nvmax(lx,ic)
             iorb=iorb+1
             do mx = -lx,lx
                ix = ix+1
                iord(mx,nx,lx,ibas)=ix !        write(ifoc,"(10i6)")mx,nx,lx,ibas,ix,iorb
             enddo
          enddo
       enddo
    enddo
    if(ix/=ndima) call rx( 'rdata4gw:ix/=ndima')
    ! Refining mesh for GW. Get new  nrmx, gval, gcore, aa(ic),bb(ic),nr(ic)
    !   For given two conditions; a. dr/dI (delrset() in switches.F), b. keeping dr/dI at r=0 (= a*b),
    !   we can deternie required nr(ic), a(ic), b(ic).
!    write(stdo,ftox)'rmeshrefine:'
    do ibas=1,nbas !iclass is crytalographically equivalent atoms. 
       ispecc(iclass(ibas)) = ispec(ibas) !note iclass(ibas) ispec(ibas). The same iclass should have the same ispec.
    enddo
    delr = delrset()       !delr is dr/di at rmax in a.u.
!    write(stdo,*)' meshrefine : delr nclass=',delr,nclass
    allocate(nrnn(nclass),aann(nclass)) !,bbnn(nclass))
    do ic = 1,nclass
       is= ispecc(ic)
       bb(is) = rmt(is)/(exp(spec_a(is)*(nr(is)-1))-1d0)
       aa_n= (delr-spec_a(is)*bb(is))/rmt(is) 
       bb_n= spec_a(is)*bb(is)/aa_n
!       write(stdo,"(' ic aa bb=',i5,2d13.6,i5,' rmt aa_n bb_n=',3d13.6)") ic,spec_a(is),bb(is),nr(is),rmt(is),aa_n,bb_n
       ir = findloc([(bb_n*( exp(aa_n*(i-1))-1d0 ) >rmt(is),i=1,100000)],value=.true.,dim=1)
       nrnn(ic) = (ir/2)*2 + 1 !odd for simpson integral later on.
       aann(ic) = aa_n
    enddo
    nrmx = maxval(nrnn(1:nclass))
!    write(stdo,*) ' New nrmx =',nrmx
    allocate( gval_n(nrmx,0:lmxax, nphimx, nsp,nclass), gcore_n(nrmx, ncoremx, nsp,nclass) )
    do ic = 1,nclass
       !     if(minval(abs(iclass-ic))/=0) cycle !jan2008
       is=ispecc(ic)
!       write(stdo,"('  input  nr a b =',i5,3d13.6)") nr(is),spec_a(is),bb(is)
       allocate(rofi,source = [(bb(is)*(exp((ir-1)*spec_a(is))-1d0),ir=1,nr(is))])
       nr_n = nrnn(ic)
       aa_n = aann(ic)
       bb_n = rmt(is)/(exp(aa_n*(nr_n-1))-1d0)
       allocate(rofi_n,source=[(bb_n*(exp((ir-1)*aa_n)-1d0),ir=1,nr_n)])
       do isp = 1, nsp
          do lx = 0,lmxa(is)
             do nx = 1,nvmax(lx,ic)
                call rrefine(rofi,nr(is),rofi_n,nrnn(ic), gval(1, lx,nx,isp,ic), gval_n(1, lx, nx, isp,ic) )
             enddo
          enddo
          do icore = 1,ncores(is) !  write(stdo,ftox)' ggg2 gcore: icore2 isp ic=',icore, isp,ic, gcore(100,icore,isp,ic)
             call rrefine(rofi,nr(is),rofi_n,nrnn(ic), gcore(1,icore,isp,ic), gcore_n(1,icore, isp,ic) )
          enddo
       enddo
       aac(ic) = aa_n
       bbc(ic) = bb_n
       nrc(ic) = nr_n
       deallocate(rofi,rofi_n)
!       write(stdo,"(' output  nr a b =',i5,3d13.6)") nrc(ic),aac(ic),bbc(ic)
    enddo   ! scaled gval to avoid degeneracy of overalp OrthoNormalized
    allocate( gx_in(nrmx,0:lmxax,nphimx,nsp,nclass), gval_orth(nrmx,0:lmxax,nphimx,nsp,nclass) )
    allocate( zzp(nmax,nmax,0:lmxax,nclass,nsp),zzpi(nmax,nmax,0:lmxax,nclass,nsp),eb(nmax))
    allocate( nncx(0:lmxax,nbas),source=0)
    do ibas= 1,nbas
       ic = iclass(ibas)
       is = ispec(ibas)
       do l = 0,lmxa(is)
          nncx(l,ibas) = konf0(l,ic) -1 -(l+1) +1
       enddo
    enddo
    nnc = maxval(nncx)
    ! zzpi --- This section it to keep the numerical stability when we have degeneracy. (mainly in the case of orthnormalized input)
    !  ovv(1:nm,1:nm,l1,ic,isp) ---> zzp (1:nm,1:nm,l1,ic,isp), where nm = nvmax(l1,ic) \sum_j ovv(i,j)*zzp0(j,k)= e_k zzp0(i,k)
    !     zzp(i,j) = zzp0(i,k)  /sqrt(e_k)  !   zzpi = inverse of zzp
    allocate(ovv(nphimx,nphimx,0:lmxax, nclass,nsp) )
    ibasloop: do 1010 ibas = 1,nbas
       ic    = iclass(ibas)
       is    = ispec(ibas)
!       write(stdo,*); write(stdo,*)' ### ibas ic =',ibas,ic
!       write(stdo,"(4i4,2d14.6)")nrc(ic),lmxa(is),nsp,ncores(is),aac(ic),bbc(ic)
       do isp = 1, nsp ! ... Get overlap matrix ovv of radial functions.
          do l1  = 0,lmxa(is)
             do nn = 1, nvmax(l1,ic)
                gx_in (1:nrc(ic),l1,nn,isp,ic)=gval_n(1:nrc(ic),l1,nn,isp,ic)*sqrt(1d0+0.1d0*nn) !sqrt(1d0+0.1d0*nn) is to avoid degeneracy. See zzp
             enddo
          enddo
          do irad1 = 1,nrad(ibas)
             do irad2 = 1,nrad(ibas)
                l1 = lindx_r (irad1,ibas); n1 = nindx_r (irad1,ibas)
                l2 = lindx_r (irad2,ibas); n2 = nindx_r (irad2,ibas)
                if(l1/=l2) cycle
                call gintxx( gx_in(1,l1,n1,isp,ic),gx_in(1,l1,n2,isp,ic),aac(ic),bbc(ic), nrc(ic), ovv(n1,n2,l1,ic,isp) )
             enddo
          enddo
          do l1  = 0,lmxa(is)
             n1 = nvmax(l1,ic) 
             call rss(n1, ovv(1:n1,1:n1,l1,ic,isp), eb, zzp(1:n1,1:n1,l1,ic,isp), ierr) !!Get zzp : eigenfunctions of ovv=<gx_in|gx_in>
!             write(stdo,"(' eb=',10f12.6)") eb(1:n1)
             if(ierr/=0) call rx(' rdata4gw: error in rs ')
             forall(i2=1:n1) zzp(:,i2,l1,ic,isp)=zzp(:,i2,l1,ic,isp)/eb(i2)**.5
             allocate(zzpx(1:n1,1:n1))
             zzpx=zzp(1:n1,1:n1,l1,ic,isp) !Get zzpi : inverse of zzp
             call matcinv(n1,zzpx)
             zzpi(1:n1,1:n1,l1,ic,isp) = dreal(zzpx) !connect gx_in and gval_orth
             deallocate(zzpx)
             forall(ir=1:nrc(ic)) gval_orth(ir,l1,1:n1,isp,ic)= matmul(gx_in(ir,l1,1:n1,isp,ic),zzp(1:n1,1:n1,l1,ic,isp))
          enddo
       enddo
1010 enddo ibasloop
  endsubroutine rdata1init
endmodule m_rdata1
subroutine rrefine(rofio,nro,rofin,nrn,go, gn )
  implicit none 
  intent(in):: rofio,nro,rofin,nrn,go
  intent(out):: gn
  integer:: nro,nrn,ir
  real(8):: polinta,rofio(nro),rofin(nrn),go(nro),gn(nrn)
  gn = [(polinta(rofin(ir), rofio,go,nro),ir=1,nrn)]
endsubroutine rrefine
