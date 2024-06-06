module m_rdovfa
  use m_lmfinit,only: z_i=>z,nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,lmxb_i=>lmxb,lmxl_i=>lmxl,spec_a
  use m_lmfinit,only: kmxt_i=>kmxt,lfoca_i=>lfoca,rfoca_i=>rfoca,rsmv_i=>rsmv
  use m_ll,only:ll
  public rdovfa
contains
  subroutine rdovfa() !- Read atm files and overlap free atom densities.
    use m_density,only: zv_a_osmrho=>osmrho,sv_p_orhoat=>orhoat,v1pot,v0pot,eferm !Outputs. allocated

    use m_supot,only: ng=>lat_ng,rv_a_ogv,iv_a_okv,rv_a_ogv,n1,n2,n3
    use m_lmfinit,only:alat=>lat_alat,nsp,nbas,nspec,ispec,qbg=>zbak,slabl,v0fix
    use m_lattic,only: vol=>lat_vol
    use m_struc_def,only: s_rv1,s_rv2
    use m_ext,only: sname
    use m_lgunit,only:stdo,stdl
    use m_ftox
    use m_MPItk,only: master_mpi,comm
    use m_fatom,only:sspec,mpibc1_s_spec
    !i Inputs
    !i   nbas  :size of basis
    !i   nspec :number of species
    !i   sspec :struct containing species-specific information
    !i   qbg  :constant background charge  qcore+qval-qbg = \sum_i Zi
    !o Outputs
    !o   orhoat: vector of offsets containing site density, in standard
    !o           3-component form (true rho, smoothed rho, core rho)
    !o   smrho :smoothed interstitial density, complex (for computational convenience)
    ! ----------------------------------------------------------------------
    implicit none
    integer :: nrmx, n0,i_spec
    parameter ( nrmx=1501, n0=10 )
    integer:: nxi(nspec)
    type(s_rv1) :: rv_a_orhofa(nspec)
    type(s_rv2) :: rv_a_ov0a(nspec)
    real(8):: rsmfa(nspec),exi(n0,nspec), hfc(n0,2,nspec),hfct(n0,2,nspec),&
         a,rmt,z,rfoc,z0,rmt0,a0,qc,ccof, &
         ceh,stc,ztot,ctot,corm,ssum,fac,sum1,sum2,sqloc,dq,smom, slmom,qcor(2)
    character(8) :: spid(nspec),spidr
    integer:: ipr,iprint,iofa,kcor,lcor, i,ifi,is, nr,lfoc,nr0,i1,nch,ib,igetss,lmxl,nlml,ierr 
    real(8) ,allocatable :: rwgt_rv(:)
    complex(8) ,allocatable :: cv_zv(:)
    character msg*23, strn*120
    logical :: lfail, l_dummy_isanrg,isanrg,mlog
    include 'mpif.h'
    call tcn('rdovfa')
    ipr   = iprint()
    msg   = '         File mismatch:'
    if(ipr>=10) write(stdo,"(/'rdovfa: read and overlap free-atom densities',' (mesh density) ...')")
    hfc=0d0
    exi=0d0
    hfc=0d0
    hfct=0d0
    if (master_mpi) open(newunit=ifi,file='atm.'//trim(sname))  !! Read free-atom density for all species ---
    isloop: do  10  is = 1, nspec
       spid(is)=slabl(is)
       a= spec_a(is)
       nr=nr_i(is)
       allocate(rv_a_orhofa(is)%v(nr*nsp),   source=0d0)
       if(allocated(sspec(is)%rv_a_orhoc)) deallocate(sspec(is)%rv_a_orhoc)
       allocate(sspec(is)%rv_a_orhoc(nr*nsp),source=0d0)
       allocate(rv_a_ov0a(is)%v(nr,nsp),     source=0d0)
       rmt=rmt_i(is)
       z=z_i(is)
       lfoc=lfoca_i(is)
       rfoc=rfoca_i(is)
       if (master_mpi) then
          if (z == 0 .AND. rmt == 0) then
             nxi(is) = 0
             rsmfa(is) = 0
             z0=0
             rmt0=0
             a0=0
             nr0=0
             qc=0
             ccof=0
             ceh=0
             stc=0
             if (allocated(rv_a_ov0a(is)%v)) deallocate(rv_a_ov0a(is)%v)
             if (allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
          else
             nr0=nrmx 
             lfail = .false.
             lfail = iofa ( spidr,n0,nxi ( is ),exi ( 1,is ),hfc &
                  ( 1,1,is ),hfct ( 1,1,is ),rsmfa ( is ),z0,rmt0 &
                 ,a0,nr0,qc,ccof,ceh,stc,rv_a_orhofa( is )%v,sspec &
                 ( is ) %rv_a_orhoc,rv_a_ov0a ( is ) %v,ifi,'read' )    < 0
             if(lfail) call rx('you did lmfa? Readin freeatom')
          endif
       endif
       call mpi_barrier(comm,ierr)
       call mpibc1(nr0,1,2,mlog,'rdovfa','nr0')
       call mpibc1(nxi(is),1,2,mlog,'rdovfa','nxi')
       call mpibc1(exi(1,is),nxi(is),4,mlog,'rdovfa','exi')
       call mpibc1(hfc(1,1,is),nsp*n0,4,mlog,'rdovfa','hfc')
       call mpibc1(hfct(1,1,is),nsp*n0,4,mlog,'rdovfa','hfct')
       call mpibc1(rsmfa(is),1,4,mlog,'rdovfa','rsmfa')
       call mpibc1(a0,1,4,mlog,'rdovfa','a0')
       call mpibc1(rv_a_orhofa( is )%v,nr0 * nsp,4,mlog,'rdovfa' ,'rhofa' )
       call mpibc1(sspec(is)%rv_a_orhoc,nr0*nsp,4,mlog,'rdovfa','rhoca')
       call mpibc1( rv_a_ov0a( is )%v,nr0 * nsp,4,mlog,'rdovfa', 'v0a' )
       if(master_mpi.and. ipr >= 30 .AND. rmt0 /= 0) write(stdo,400) trim(spid(is)),spidr,rmt0,nr0,a0
400       format(' rdovfa: expected ',a,',',T27,' read ',a, ' with rmt=',f8.4,'  mesh',i6,f7.3)
       if(nr <= 0)   nr = nr0
       if(a <= 1d-6) a = a0
       if(z == 0 .AND. rmt == 0) then
          a = 0
          nr = 0
       endif
       sspec(is)%qc=qc
       sspec(is)%rsmfa=rsmfa(is)
       sspec(is)%ctail=ccof
       sspec(is)%etail=ceh
       sspec(is)%stc=stc
       sspec(is)%nxi=nxi(is)
       sspec(is)%exi=exi(:,is)
       sspec(is)%chfa=hfc(:,:,is)
10  enddo isloop
    do i_spec=1,nspec ! Re-broadcast entire species structure, and arrays used below
       call mpibc1_s_spec(sspec(i_spec))!,'rdovfa_sspec')
    enddo
    if (master_mpi) close(ifi)
    ! --- Define arrays for local densities rho1,rho2,rhoc and v0,v1 ---
    ztot = 0d0
    ctot = 0d0
    corm = 0d0
    if(allocated(v0pot))deallocate(v0pot,v1pot)!this may cause mem leak?(v0pot%v is not deallocated).
    allocate(v1pot(nbas),v0pot(nbas))
    if(allocated(sv_p_orhoat)) deallocate(sv_p_orhoat)
    allocate(sv_p_orhoat(3,nbas))
    ibloop: do  20  ib = 1, nbas
       is = ispec(ib) 
       a=   spec_a(is)
       nr=  nr_i(is)
       rmt= rmt_i(is)
       lmxl=lmxl_i(is)
       z=z_i(is)
       qc=sspec(is)%qc
       lfoc=lfoca_i(is)
       nlml = (lmxl+1)**2
       allocate(sv_p_orhoat(1,ib)%v(nr*nlml*nsp))
       allocate(sv_p_orhoat(2,ib)%v(nr*nlml*nsp))
       allocate(sv_p_orhoat(3,ib)%v(nr*nsp))
       if (nsp == 2 .AND. lmxl > -1) then
          allocate(rwgt_rv(nr))
          call radwgt ( rmt,a,nr,rwgt_rv )
          call radsum ( nr,nr,1,nsp,rwgt_rv,sspec(is)%rv_a_orhoc,ssum )
          call radsum ( nr,nr,1,1,rwgt_rv,sspec(is)%rv_a_orhoc,sum1 )
          sum2 = ssum - sum1
          call gtpcor(is,kcor,lcor,qcor)
          if(dabs(qcor(2)-(sum1-sum2)) > 0.01d0.and.ipr>=10) &
               write(stdo,ftox)' (warning) core moment mismatch spec ',is,'input file=',ftof(qcor(2)),'atom file=',ftof(sum1-sum2)
          corm = corm + qcor(2)
          deallocate(rwgt_rv)
       endif
       if (lmxl > -1) then
          allocate(v0pot(ib)%v(nr,nsp))
          allocate(v1pot(ib)%v(nr,nsp))
          v0pot(ib)%v=rv_a_ov0a( is )%v
          v1pot(ib)%v=rv_a_ov0a( is )%v
          sv_p_orhoat(3,ib)%v(1:nr*nsp) = sspec(is)%rv_a_orhoc(1:nr*nsp)
          if (lfoc == 0) then
             allocate(rwgt_rv(nr))
             call radwgt ( rmt,a,nr,rwgt_rv )
             call radsum ( nr,nr,1,nsp,rwgt_rv,sv_p_orhoat( 3,ib )%v,ssum )
             fac = merge(qc/ssum,1d0,dabs(ssum) > 1d-7)
             if(ipr>=40) write(stdo,"(' scale foca=0 core species',i2,': qc,sum,scale=', 3f12.6,f12.6)") is,qc,ssum,fac
             sv_p_orhoat( 3,ib )%v(1:nr*nsp)=fac*sv_p_orhoat(3,ib)%v(1:nr*nsp)
             deallocate(rwgt_rv)
          endif
       endif
       ztot = ztot+z
       ctot = ctot+qc
20  enddo ibloop
    v0wrireblock:block
      integer:: ir,isp
      character(8):: charext
      if(v0fix.and.master_mpi) then
         do ib=1,nbas
            is = ispec(ib) 
            nr = nr_i(is)
            do ir=1,nr
               v0pot(ib)%v(ir,:)= sum(v0pot(ib)%v(ir,1:nsp))/nsp
            enddo
            write(stdo,*)' v0fix=T: writing v0pot',ib,nr
            open(newunit=ifi,file='v0pot.'//trim(charext(ib)),form='unformatted')
            write(ifi) v0pot(ib)%v
            close(ifi)
         enddo
      endif
    endblock v0wrireblock
    if(allocated(zv_a_osmrho)) deallocate(zv_a_osmrho)
    allocate(zv_a_osmrho(n1*n2*n3,nsp),source=(0d0,0d0))
    eferm=0d0
    allocate(cv_zv(ng*nsp))
    call ovlpfa( nbas,nxi,n0,exi,hfc,rsmfa, ng,rv_a_ogv,cv_zv ) !Overlap smooth hankels to get smooth interstitial density ---
    call gvputf( ng,nsp,iv_a_okv,n1,n2,n3,cv_zv,zv_a_osmrho )
    deallocate(cv_zv)
    call fftz3(zv_a_osmrho,n1,n2,n3,n1,n2,n3,nsp, 0,1 )! ... FFT to real-space mesh
    zv_a_osmrho=zv_a_osmrho-qbg/vol/nsp   !Add compensating uniform density to compensate background
    sum1 = dreal(sum(zv_a_osmrho(:,:)))*vol/(n1*n2*n3)
    if(nsp==2) smom = 2d0*dreal(sum(zv_a_osmrho(:,1)))*vol/(n1*n2*n3) - sum1
    call ovlocr(nbas,n0,nxi,exi,hfc,rsmfa,rv_a_orhofa,sv_p_orhoat,sqloc,slmom) !Set up local densities using rmt from atm file ---
    call adbkql ( sv_p_orhoat,nbas,nsp,qbg,vol,- 1d0 ) ! Add compensating uniform electron density to compensate background
    if(abs(qbg)/=0d0.and. ipr>=10) write(stdo,ftox) ' Uniform density added to neutralize background q=',ftof(qbg)
    dq = sum1+sqloc+ctot-ztot+qbg !charge
    if (nsp == 1.and.ipr>=10) then
       write(stdo,895) sum1,sqloc,sum1+sqloc,ctot,-ztot,qbg,dq
895    format(/  ' Smooth charge on mesh:    ',f16.6 &
            /    ' Sum of local charges:     ',f16.6 &
            /    ' Total valence charge:     ',f16.6 &
            /    ' Sum of core charges:      ',f16.6 &
            /    ' Sum of nuclear charges:   ',f16.6 &
            /    ' Homogeneous background:   ',f16.6 &
            /    ' Deviation from neutrality:',f16.6)
       write (stdl,710) sum1+sqloc,sum1,sqloc,qbg,dq
710    format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,'   bg',f10.6,'  dQ',f10.6)
    elseif(ipr>=10) then
       write(stdo,896) sum1,smom,sqloc,slmom,sum1+sqloc,smom+slmom,ctot,corm,-ztot,qbg,dq
896    format(/  ' Smooth charge on mesh:    ',f16.6,4x,'moment', f12.6, &
            /    ' Sum of local charges:     ',f16.6,4x,'moments',f11.6, &
            /    ' Total valence charge:     ',f16.6,4x,'moment', f12.6, &
            /    ' Sum of core charges:      ',f16.6,4x,'moment', f12.6, &
            /    ' Sum of nuclear charges:   ',f16.6 &
            /    ' Homogeneous background:   ',f16.6 &
            /    ' Deviation from neutrality:',f16.6)
       write (stdl,711) sum1+sqloc,sum1,sqloc,qbg,smom+slmom
711    format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,'   bg',f11.6,' mm',f11.6)
    endif
    if(dabs(dq)>1d-4.AND.ipr>0) write(stdo,"(' rdovfa (warning) overlapped density not neutral'//', dq= ',d13.5)") dq
    do is=1,nspec
       if(allocated(rv_a_ov0a(is)%v)) deallocate(rv_a_ov0a(is)%v)
       if(allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
    enddo
    call tcx('rdovfa')
  end subroutine rdovfa
  subroutine ovlocr(nbas,nxi0,nxi,exi,hfc,rsmfa,rv_a_orhofa, sv_p_orhoat,sqloc, slmom )!- Makes the site densities for overlapped free atoms.
    use m_lmfinit,only: nsp,ispec
    use m_struc_def
    use m_lgunit,only:stdo
    use m_smhankel,only: hxpbl
    use m_lattic,only: rv_a_opos
    use m_smhankel,only: hxpos
!    use m_hansr,only:corprm
    !i   nbas  :size of basis
    !i   nxi   :number of Hankels
    !i   nxi0  :leading dimension of hfc
    !i   exi   :smoothed Hankel energies; see Remarks
    !i   hfc   :coefficients to smoothed Hankels
    !i   rsmfa :Hankel smoothing radius
    !i   orhofa:free-atom density, by species
    !o Outputs
    !   orhoat :local density, given by true and smooth densities
    !    sqloc :sum of local charges (integral over rho1-rho2)
    !    slmom :sum of local magnetic moments
    implicit none
    integer:: kmxv=15 !k-cutoff for 1-center projection of free-atom rho. ! kmxv = cutoff to expand smoothed potential.    hardwired for now
    integer:: nbas,nxi(*),nxi0
    type(s_rv1) :: sv_p_orhoat(3,nbas)
    type(s_rv1) :: rv_a_orhofa(nbas)
    real(8):: rsmfa(*),exi(nxi0,*),hfc(nxi0,2,*),sqloc,slmom
    integer:: ib,ipr,iprint,is,jb,je,js,lfoca,lmxl,nlmh,nlml,nr,i
    real(8) :: ceh,cofg,cofh,eh,qcorg,qcorh,qsc,qcsm,qloc,rfoca,rmt,rsmh,rsmv,z,amom, a,p1(3), p2(3)
    real(8),allocatable:: acof(:,:,:),b0(:,:),rofi(:),rwgt(:)
    complex(8),allocatable:: b(:,:)
    integer:: ibini,ibend
    call tcn('ovlocr')
    ipr  = iprint()
    sqloc = 0
    slmom = 0
    if(ipr>=30) write(stdo,300)
300 format(/' Free atom and overlapped crystal site charges:'/'   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local')
    ibini= 1
    ibend= nbas
    do  ib = ibini,ibend
       is=ispec(ib) 
       lmxl=lmxl_i(is)
       if (lmxl == -1) cycle
       p1(:)=rv_a_opos(:,ib) 
       rsmv=rsmv_i(is)
       nlml = (lmxl+1)**2
       allocate(acof(0:kmxv,nlml,nsp),source=0d0)
       allocate(b(0:kmxv,nlml),source=(0d0,0d0))
       a= spec_a(is)
       nr=nr_i(is)
       rmt=rmt_i(is)
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoca,rfoca, z)
       qcsm = qcorg+qcorh
       allocate(rofi(nr),rwgt(nr))
       call radmsh(rmt,a,nr,rofi)
       call radwgt(rmt,a,nr,rwgt)
       do  jb = 1, nbas !add up tail expansion
          js=ispec(jb) 
          p2(:)=rv_a_opos(:,jb) 
          do  je = 1, nxi(js)
             rsmh = rsmfa(js)
             eh   = exi(je,js)
             nlmh = 1
             call hxpbl(p2,p1,[0d0,0d0,0d0],[rsmh],rsmv,[eh],kmxv,nlmh,nlml,kmxv,nlml, b) 
             allocate(b0(0:kmxv,nlmh),source=0d0)
             if (ib == jb) call hxpos([rsmh],rsmv,[eh],kmxv,nlmh,kmxv,b0)
             do  i = 1, nsp
                call p1ovlc(kmxv,nlml,hfc(je,i,js),b,b0,acof(0,1,i))
             enddo
             deallocate(b0)
          enddo
       enddo
       call p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,acof,rofi &
           ,rwgt,nxi0,nxi(is),exi(1,is),hfc(1,1,is),rsmfa(is),rv_a_orhofa(is)%v,sv_p_orhoat(3,ib)%v &
           ,lfoca,qcsm,qloc,amom,sv_p_orhoat(1,ib)%v,sv_p_orhoat( 2,ib )%v )
       sqloc = sqloc + qloc
       slmom = slmom + amom
       deallocate(rofi,rwgt,acof,b)
    enddo
    call tcx('ovlocr')
  end subroutine ovlocr
  subroutine p1ovlc(kmxv,nlml,hfc,b,b0,a)!- Adds contribution to P_kl expansion of density from one basis function
    !i Inputs
    !i   nlml  :density expanded to nlml
    !i   kmxv  :k-cutoff for P_kl expansion
    !i   hfc   :coefficient to basis function
    !i   b     :P_kl expansion of density from one basis function
    !i   b0    :P_kl expansion of on-site density
    !o Outputs
    !o   a     :cumulative P_kl expansion of density for this site
    use m_ftox
    implicit none
    integer :: nlml,kmxv
    real(8) :: a(0:kmxv,nlml),b0(0:kmxv,1),hfc
    double complex b(0:kmxv,nlml)
    integer :: k,ilm
!    write(6,ftox) 'ppppppp111111111',hfc,sum(b),sum(b0(:,1))
    a(:,:) = a(:,:) + hfc*b(:,:)
    a(:,1) = a(:,1) - hfc*b0(:,1)
  end subroutine p1ovlc
  subroutine p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,acof,rofi,rwgt,nxi0,nxi,exi,hfc,rsmfa,rhofa,rhoc,lfoca,qcsm,qloc,amom,rho1,rho2) !Assemble local density from P_kl expansion for one site
    use m_lgunit,only:stdo
    use m_hansmr,only: hansmr,hansmronly
    use m_hansr,only:  hansr
    !i   ib    :site for which to assemble local density
    !i   nsp   :number of spin channels
    !i   rsmv  :smoothing radius for P_kl expansion
    !i   kmxv  :k-cutoff for P_kl expansion
    !i   nr    :number of radial mesh points
    !i   nlml  :L-cutoff for P_kl expansion
    !i   acof  :coefficients to P_kl expansion
    !i   rofi  :radial mesh points for tabulation on a radial mesh
    !i   rwgt  :radial mesh weights for integration on a radial mesh
    !o   rhohd :work array (holds on-site smoothed head density)
    !i   nxi0  :leading dimension of hfc
    !i   nxi   :number of smoothed Hankel energies in head expansion
    !i   exi   :smoothed Hankel energies in head expansion
    !i   hfc   :coefficients to Hankel energies in head expansion
    !i   rsmfa :Hankel smoothing radius in head expansion
    !i   rhofa :head free-atom density
    !i   rhoc  :core density --- used to integrate core charge
    !i   lfoca :switch specifying treatment of core density.
    !i          0 => val,slo = 0 at sphere boundary
    !i          1 => core tails included explicitly with valence
    !i          2 => tails included perturbatively
    !i   qcsm  :smoothed core density
    !o Outputs
    !i   rho1  :local true density, tabulated on a radial mesh
    !i   rho2  :local smoothed density, tabulated on a radial mesh
    !o   qloc  :sphere charge
    !o   amom  :sphere magnetic moment
    implicit none
    integer :: nr,nxi0,ib,nsp,kmxv,nlml,nxi,lfoca, i,ie,ilm,ipr,iprint,k,l,lmax,lmxl,isp
    real(8) :: qcsm,qloc,rhofa(nr,nsp),rho1(nr,nlml,nsp), &
         rho2(nr,nlml,nsp),rofi(nr),rwgt(nr),rhohd(nr,nsp),exi(*), &
         hfc(nxi0,nsp),rhoc(nr,nsp),acof(0:kmxv,nlml,nsp),rsmv,rsmfa,amom
    real(8) :: asm,gam,pi,qall,qexa,qin,qlc,qnum,qout,qsmo,qut, &
         r,rl,rmt,srfpi,ssum,sumfa,sumhd,sumsm,sumtr,y0, xi(0:10),x0(0:2),ddot !pkl(0:kmx,0:lmx)
    real(8),allocatable:: pkl(:,:)
    ipr   = iprint()
    pi    = 4d0*datan(1d0)
    srfpi = dsqrt(4*pi)
    y0    = 1d0/srfpi
    lmxl  = ll(nlml)
    allocate(pkl(0:kmxv,0:lmxl))
    qnum = 0d0
    qexa = 0d0
    qsmo = 0d0
    qut = 0d0
    rhohd=0d0
    asm = 1d0/rsmfa
    lmax = 0
    do  ie = 1, nxi
       ssum = 0d0
       do  i = 1, nr
          r = rofi(i)
          call hansmr(r,exi(ie),asm,xi,lmax)
          ssum = ssum + srfpi*rwgt(i)*xi(0)*r*r
          do  isp = 1, nsp
             rhohd(i,isp) = rhohd(i,isp) + srfpi*hfc(ie,isp)*xi(0)*r*r !Assemble smooth on-site head density in rhohd ---
          enddo
       enddo
       gam = 0.25d0*rsmfa**2
       qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
       rmt = rofi(nr)
       call hansmr(rmt,0d0,1/rsmfa,x0,1)
       call hansmr(rmt,exi(ie),1/rsmfa,xi,1)
       qout = srfpi/exi(ie)*(-dexp(rsmfa**2/4*exi(ie)) - rmt**3*(xi(1)-dexp(rsmfa**2/4*exi(ie))*x0(1)))
       qin = qall-qout
       do  isp = 1, nsp
          qnum = qnum + hfc(ie,isp)*ssum
          qexa = qexa + hfc(ie,isp)*qin
          qsmo = qsmo + hfc(ie,isp)*qall
          qut  = qut  + hfc(ie,isp)*qout
       enddo
    enddo
    rho2=0d0 
    do  i = 1, nr
       r = rofi(i)
       call radpkl(r,rsmv,kmxv,lmxl,kmxv,pkl)
       do  ilm = 1, nlml
          l = ll(ilm)
          rl = 0d0
          if( r > 0d0 ) rl = r**l
          do k = 0, kmxv
             rho2(i,ilm,:) = rho2(i,ilm,:) + acof(k,ilm,:)*pkl(k,l)*r*r*rl !Assemble overlapped tail density in rho2 ---
          enddo
       enddo
    enddo
    rho1=rho2
!!!!!!!!!!!!!!!!!!!!
!       block
!         use m_density,only: orhoat
!         write(6,*)'rrrrrrdensity111222xxx',sum(rho1),sum(abs(rho1)),sum(abs(acof)),sum(abs(pkl))
!       endblock
!!!!!!!!!!!!!!!!!!!!    
    do   i = 1, nr ! ... Make the true density in rho1, smooth density in rho2     !call dpcopy(rho2,rho1,1,nr*nlml*nsp,1d0)
       rho1(i,1,:) = rho1(i,1,:) + y0*rhofa(i,:)
       rho2(i,1,:) = rho2(i,1,:) + y0*rhohd(i,:)
    enddo
    sumfa = sum([(rwgt(i)*sum(rhofa(i,:)),i=1,nr)])
    sumhd = sum([(rwgt(i)*sum(rhohd(i,:)),i=1,nr)])
    sumsm = sum([(rwgt(i)*sum(rho2(i,1,:)),i=1,nr)])*srfpi
    sumtr = sum([(rwgt(i)*sum(rho1(i,1,:)),i=1,nr)])*srfpi
    qlc   = sum([(rwgt(i)*sum(rhoc(i,:)),i=1,nr)])
    qloc = sumtr-sumsm
    amom = -srfpi* sum([(rwgt*(rho1(:,1,isp)-rho2(:,1,isp))*(3-2*isp),isp=1,nsp)])
    !amom = -srfpi* (ddot(nr,rwgt,1,rho1(1,1,nsp),1)-ddot(nr,rwgt,1,rho1(1,1,1),1)-ddot(nr,rwgt,1,rho2(1,1,nsp),1)+ddot(nr,rwgt,1,rho2(1,1,1),1))
    if (lfoca == 0) qloc = qloc + qlc - qcsm
    if(ipr>=30)            write(stdo,"(i5,6f12.6)") ib,sumfa,sumhd,sumtr,sumsm,qloc
    if(ipr>=30.and.nsp==2) write(stdo,"(' amom',6f12.6)") ddot(nr,rwgt,1,rhofa,1)-ddot(nr,rwgt,1,rhofa(1,2),1), &
            ddot(nr,rwgt,1,rhohd,1)-ddot(nr,rwgt,1,rhohd(1,2),1), &
            srfpi*(ddot(nr,rwgt,1,rho1,1)-ddot(nr,rwgt,1,rho1(1,1,2),1)), &
            srfpi*(ddot(nr,rwgt,1,rho2,1)-ddot(nr,rwgt,1,rho2(1,1,2),1)), amom
  end subroutine p2ovlc
  subroutine adbkql( sv_p_orhoat,nbas,nsp,qbg,vol,fac )!- Add uniform bkg charge density to local smooth rho
    use m_struc_def
    use m_lmfinit,only: ispec
    !i orhoat: pointers to local density in spheres
    !i nbas: number of atoms in basis
    !i qbg: background charge
    !i nsp: spins
    !i vol: vol of cell
    !i fac: fac * backg density is added
    implicit none
    integer :: nrmx,nlmx,nlml,lmxl,nbas, nsp,ib,nr,is
    parameter (nrmx=1501)
    type(s_rv1) :: sv_p_orhoat(3,nbas)
    real(8):: qbg,fac,rhobkg,vol,a,rmt,rofi(nrmx)
    rhobkg = fac*qbg/vol
    do  ib = 1, nbas
       is=ispec(ib)
       lmxl=lmxl_i(is)
       if (lmxl == -1) cycle
       a=spec_a(is)
       nr=nr_i(is)
       rmt=rmt_i(is)
       nlml=(lmxl+1)**2
       call rxx(nr   > nrmx,'addbkgloc: increase nrmx')
       call radmsh(rmt,a,nr,rofi)
       call addbkgl(sv_p_orhoat(1,ib )%v,sv_p_orhoat(2,ib)%v, rhobkg,nr,nsp,rofi,nlml )
    enddo
  end subroutine adbkql
  subroutine addbkgl(rho1,rho2,rhobkg,nr,nsp,rofi,nlml) ! adds uniform background to local smooth density at this site for l=0 component (ilm=1) 
    implicit none
    integer :: nsp,is,nr,nlml,i
    real(8):: rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rofi(nr),rhobkg
    real(8),parameter:: pi = 4d0*datan(1d0),srfpi = dsqrt(4*pi)
    do is = 1, nsp
       rho1(:,1,is) = rho1(:,1,is)+srfpi*rofi(:)**2*rhobkg/nsp
       rho2(:,1,is) = rho2(:,1,is)+srfpi*rofi(:)**2*rhobkg/nsp
    enddo
  end subroutine addbkgl
  subroutine ovlpfa(nbas,nxi,nxi0,exi,hfc,rsmfa,ng,  gv,cv)!- Set up Fourier coeffs to overlap the smooth part of FA densities.
    use m_lmfinit,only:lat_alat,nsp,ispec
    use m_lattic,only: lat_vol,rv_a_opos
    use m_lgunit,only:stdo
    use m_ftox
    use m_struc_def
    !i   nxi   :number of Hankels
    !i   nxi0  :leading dimension of hfc
    !i   exi   :smoothed Hankel energies; see Remarks
    !i   hfc   :coefficients to smoothed Hankels
    !i   rsmfa :Hankel smoothing radius
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (glist.f)
    !o Outputs
    !o   cv    :Fourier coefficients
    implicit none
    integer :: nbas,nxi(*),nxi0,ng
    real(8)::gv(ng,3),rsmfa(*),exi(nxi0,*),hfc(nxi0,2,*),v(3),alat,vol,tpiba,ssum(2),pos(3),sam(2),e,cof,rsm,gam,v2
    integer :: ipr,iprint,ib,is,nx,i,ixi,ig,isp
    complex(8):: img=(0d0,1d0), phase,cv(ng,nsp)
    real(8),parameter:: pi   = 4d0*datan(1d0),  y0   = 1d0/dsqrt(4*pi)
    call tcn('ovlpfa')
    ipr  = iprint()
    if(ipr>=30) write(stdo,*)' ovlpfa: overlap smooth part of FA densities'
    alat=lat_alat
    vol=lat_vol
    tpiba = 2*pi/alat
    cv=0d0
    ssum = 0d0
    do ib=1,nbas
       is=ispec(ib) 
       pos=rv_a_opos(:,ib) 
       nx = nxi(is)
       sam = 0
       do  isp = 1, nsp
          do  ixi = 1, nx
             e = exi(ixi,is)
             cof = hfc(ixi,isp,is)
             rsm = rsmfa(is)
             gam = 0.25d0*rsm**2
             sam(isp) = sam(isp) - cof*y0*4d0*pi*dexp(gam*e)/e
             do ig = 1, ng              !       ... Loop over reciprocal lattice vectors
                v  = gv(ig,:)*tpiba
                v2 = sum(v**2)
                cv(ig,isp) = cv(ig,isp) + -4d0*pi*dexp(gam*(e-v2))/(e-v2)* cof*exp(-img*alat*sum(pos*v))*y0/vol
             enddo
          enddo
          ssum(isp) = ssum(isp) + sam(isp)
       enddo
       if(ipr>30.AND.nx>0)write(stdo,ftox)' site',ib,'spec',is,'pos',ftof(pos,4),'Qsmooth',sam(1)+sam(2),'mom', sam(1)-sam(2)
       if(ipr>= 40.and.nx>0) then
          write(stdo,"(2x,a7,16f11.3)") 'energy:',(exi(i,is),i=1,nx)
          write(stdo,"(2x,a7,16f11.3)") 'coeff:',(hfc(i,1,is),i=1,nx)
          if(nsp==2) write(stdo,"(2x,a7,16f11.3)") 'spin2:',(hfc(i,2,is),i=1,nx)
          write(stdo,'(1x)')
       endif
    enddo
    if(ipr>30) write(stdo,ftox)' total smooth Q = ',sum(ssum)  !   write(stdo,ftox)' FT(0,0,0)=',(cv(1,1)+cv(1,nsp))*vol/(3-nsp))
    if(ipr>30.and.nsp==2) write(stdo,ftox)' total moment=',ssum(1)-ssum(2)
    call tcx('ovlpfa')
  end subroutine ovlpfa
end module m_rdovfa
