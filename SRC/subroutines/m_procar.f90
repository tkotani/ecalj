module m_procar !for wriging PROCAR (VASP format)
  use m_lmfinit,only: nlmax,nsp,nbas,stdo,sspec=>v_sspec,ispec,nlmax,nspc,n0
  use m_suham,only: ndhamx=>ham_ndhamx,ndham=>ham_ndham,nspx=>ham_nspx
  use m_igv2x,only: igv2x,napw
  use m_mkpot,only: sab_rv 
  use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize,procid,xmpbnd2
  use m_qplist,only: nkp,xdatt,qplist

  public m_procar_init, m_procar_closeprocar, m_procar_writepdos, dwgtall,nchanp

  private
  integer::  nchanp=25      !total of s,p,d,f
  real(8),allocatable,protected:: dwgtall(:,:,:,:,:)
  logical,private:: isp1init=.true.,isp2init=.true.,init=.true.
  integer,private:: iprocar1,iprocar2
  logical,private:: cmdopt0,fullmesh,debug,procaron
contains
  subroutine m_procar_closeprocar()
    logical:: nexist
    close(iprocar1)
    inquire(iprocar2,exist=nexist)
    if(nexist) close(iprocar2)
  end subroutine m_procar_closeprocar
  subroutine m_procar_init(iq,isp,ef0,evl,ndimh,jsp,qp,nev,evec,ndimhx,nmx)
    use m_makusq,only: makusq
    implicit none
    complex(8):: evec(ndimhx,nmx)
    character*1000::ccc
    real(8):: ef0
    complex(8):: auasaz(3)
    real(8):: s11,s22,s33,s12,s13,s23,dwgt(100),dwgtt(100),xdat,qold(3),qp(3)
    complex(8),allocatable:: auspp(:,:,:,:,:)
    integer:: iq,isp,jspini,jspend,jspp,iprocar,iband,is,ilm,ndimh,nspc,jsp,ib,nev,i,m,l,ndimhx,nmx
    real(8):: rydberg=13.6058d0,evl(ndhamx,nspx)
    logical:: cmdopt0
    fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface')
    debug = cmdopt0('--debugbndfp')
    PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
    if(isp1init .AND. isp==1) then
       open(newunit=iprocar1,file='PROCAR.UP.'//trim(strprocid))
       isp1init=.false.
    endif
    if(isp2init .AND. isp==2) then
       open(newunit=iprocar2,file='PROCAR.DN.'//trim(strprocid))
       isp2init=.false.
    endif
    if(procaron .AND. fullmesh .AND. init) then
       allocate(dwgtall(nchanp,nbas,ndhamx,nsp,nkp))
       dwgtall=0d0
       init=.false.
    endif
    allocate( auspp(nlmax,ndhamx,3,nsp,nbas) )
    auspp = 0d0
    call makusq(nbas,[-999], nev,jsp,1,qp,evec, auspp ) !for (u,s,gz) !correct 2023-jan mode=0
    jspini=isp
    jspend=isp
    jspp=isp
    if(jspp==1) iprocar=iprocar1
    if(jspp==2) iprocar=iprocar2
    if(debug) write(stdo,*) 'm_procar',iprocar1,iprocar2,jspp,iprocar,ef0,nlmax,ndham,nspc,nsp,nbas
    ccc="ion        s       py       pz       px      dxy      dyz      dz2      dxz   dx2-y2"// &
         "      f-3      f-2      f-1       f0       f1       f2       f3"// &
         "                                                                                 tot"
    write(iprocar,*)
    write(iprocar,*)
    write(iprocar,"('k-point ',i4,' :    ',3f13.8,'     weight = -------  : x =',f15.8)")iq,qp,xdatt(iq)
    write(iprocar,*)
    do iband = 1, nev
       write(iprocar,*)
       write(iprocar,"('band ',i3,' # energy ',f13.8,' # occ. -----' )")iband,(evl(iband,jsp)-ef0)*rydberg
       write(iprocar,*)
       dwgtt=0d0
       do ib = 1, nbas
          is  = ispec(ib)!ssite(ib)%spec
          ilm = 0
          dwgt=0d0
          do  l = 0, sspec(is)%lmxa
             do  m = -l, l
                ilm = ilm+1 !ilm,ib --> evec(ix,
                auasaz = auspp(ilm,iband,1:3,jspp,ib)
                !as = auspp(ilm,iband,2,jspp,ib)
                !az = auspp(ilm,iband,3,jspp,ib)
                !Note au,as,az are coefficients for phi1*Ylm phi2*Ylm phi3*Ylm.
                ! If --ylmc, Ylm(complex) is assumed.
                !  u=phi1: linear combination of phi,phidot (val=1 slo=0) at MP
                !  s=phi2: linear combination of phi,phidot (val=0 slo=1) at MT
                !  gz=phi3:  LO  (phi + phidot) with val=slo=0 at MT is zero
                pdosc: Block
                  complex(8):: img=(0d0,1d0)
                  integer:: ilmm
                  real(8),parameter:: dsq=1d0/2d0**0.5d0
                  logical:: cmdopt0
                  if(cmdopt0('--ylmc') .AND. m/=0) then !feb 2022 !based on spherical harmonics
                     ilmm = ilm-2*m
                     if(m>0) then
                        auasaz =[&
                             dsq*(-1)**m*(auspp(ilm,iband,1,jspp,ib)-img*auspp(ilmm,iband,1,jspp,ib)),&
                             dsq*(-1)**m*(auspp(ilm,iband,2,jspp,ib)-img*auspp(ilmm,iband,2,jspp,ib)),&
                             dsq*(-1)**m*(auspp(ilm,iband,3,jspp,ib)-img*auspp(ilmm,iband,3,jspp,ib))]
                     elseif(m<0) then
                        auasaz = [dsq*(img*auspp(ilm,iband,1,jspp,ib) + auspp(ilmm,iband,1,jspp,ib)),&
                             dsq*(img*auspp(ilm,iband,2,jspp,ib) + auspp(ilmm,iband,2,jspp,ib)),&
                             dsq*(img*auspp(ilm,iband,3,jspp,ib) + auspp(ilmm,iband,3,jspp,ib))]
                     endif
                  endif
                EndBlock pdosc
                dwgt(ilm)= sum( dconjg(auasaz) & ! auasaz is for phi,phidot,pz(val=slo=0)
                     *matmul( sab_rv(:,:,l+1+n0*(jspp-1)+n0*nsp*(ib-1)),auasaz)) !bugfix 2023-4-28 based on suzuki's report for cDyN. ! sab(3,3,l+1,isp,ib)
!bug before 2023-4-28           *matmul( sab_rv(:,:,l+1+n0*(ib-1)+n0*nbas*(jspp-1)),auasaz)) 
             enddo
          enddo
          dwgtt = dwgtt + dwgt(1:ilm)
          if(ib==1)  write(iprocar,"(a)") trim(ccc)
          write(iprocar,"(i3,100(x,f8.5))")ib,(dwgt(i),i=1,nchanp),sum(dwgt)
          if(ib==nbas) write(iprocar,"('tot',100(x,f8.5))")(dwgtt(i),i=1,nchanp),sum(dwgtt)
          if(fullmesh) dwgtall(1:nchanp,ib,iband,jspp,iq) = dwgt(1:nchanp)
       enddo
    enddo
    deallocate( auspp )
  end subroutine m_procar_init
  !!--------------------------------------------------------
  subroutine m_procar_writepdos(evlall,nevmin,ef0,kpproc)
    use m_mkqp,only: nkabc=> bz_nabc
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
    use m_ext,only:sname
    use m_tetirr,only: tetirr
    real(8) evlall(:,:,:)
!    integer:: nev_(:)
    logical:: cmdopt0
    integer:: kpproc(*)
    integer,allocatable:: ipqe(:,:,:),idtete(:,:)
    integer:: ik1,ik2,ik3,iq,nkk1,nkk2,nkk3,ntete,ifip,nevmin
    real(8):: qx(3),ef0!,emin,emax
    nkk1=nkabc(1)
    nkk2=nkabc(2)
    nkk3=nkabc(3)
    !! pdos mode (--mkprocar and --fullmesh). ===
    if(debug) print *,'mmmm procid sum dwgt check=',procid,sum(dwgtall)
    call xmpbnd2(kpproc,nbas*nchanp*ndhamx,nkp,nspx,dwgtall)
    if(master_mpi) then
       allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
       iq=0
       do ik3 = 1, nkk3
          do ik2 = 1, nkk2
             do ik1 = 1, nkk1
                iq = iq+1
                ipqe(ik1,ik2,ik3)=iq
                qx = matmul(qlat, [dble(ik1-1)/nkk1,dble(ik2-1)/nkk2, dble(ik3-1)/nkk3])
                if(abs(sum(qx-qplist(:,iq)))>1d-6)call rx("bndfp: qx/=qplist something strang")
             enddo
          enddo
       enddo
       call tetirr(qlat, nkk1,nkk2,nkk3, ipqe, ntete,idtete) !tetrahedron
       write(*,*)" ntete 6*nkk1*nkk2*nkk3 nkp=",ntete,6*nkk1*nkk2*nkk3,nkp
       write(stdo,"(' pdosmode: ndhamx nsp nspx =',4i7)") ndhamx, nsp, nspx,nevmin 
       write(stdo,"(' pdosmode: nchanp nbas ef0=',2i7,3f12.5)") nchanp,nbas,ef0
       open(newunit=ifip,file='pdosdata.'//trim(sname),form='unformatted')
       write(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp 
       write(ifip) idtete,ipqe 
       write(ifip) evlall  
       write(ifip) dwgtall 
       write(ifip) ef0
       close(ifip)
       if(nkp/=nkk1*nkk2*nkk3) call rx('pdosmode but nkp/=nkk1*nkk2*nkk3')
       if( cmdopt0('--tetraw')) then
          open(newunit=ifip,file='tetradata.dat',form='unformatted')
          write(ifip) ndhamx,nkp,ntete
          write(ifip) idtete
          close(ifip)
          open(newunit=ifip,file='qplistf.dat')
          do iq=1,nkp
             write(ifip,"(3f15.8)") qplist(:,iq)
          enddo
          close(ifip)
       endif
       deallocate(idtete,ipqe) 
    endif
  end subroutine m_procar_writepdos
end module m_procar
