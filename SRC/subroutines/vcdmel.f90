! Valence-core dipole matrix elements. We may need when we use this branch.
module m_vcdmel
  use m_lmfinit,only: z_i=>z,nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,spec_a
  public vcdmel
  private
contains  
  subroutine vcdmel(nlmax,ndimh,nbandmx,nspx,nq,nsp,nspc,ef,evl,aus,nsite,isite,iclsl,iclsn,dosw)! Valence-core dipole matrix elements
    use m_lmfinit,only: rv_a_ocg,iv_a_ojcg,iv_a_oidxcg,ispec,n0,lmxax
    use m_mkqp,only: iv_a_oidtet ,bz_nabc, bz_ntet
    use m_struc_def
    use m_ext,only: sname     !extention for file
    use m_lmfinit,only: bz_ndos,bz_dosmax,slabl
    use m_elocp,only: rsmlss=>rsml,ehlss=>ehl
    use m_density,only: v0pot,pnuall,pnzall
    use m_makusp,only: makusp
!    use m_igv2x,only: nbandmx
    !i   nlmax :first dimension of aus; largest augmentation (l+1)^2
    !i   ndham :second dimension of aus, at least as large as ndimh
    !i   ndimh :number of eigenvalues
    !i   nq    :number of k-points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
    !i   ef    :Fermi energy
    !i   evl   :energy bands at the nq k-points
    !i   aus   :values and slopes of eigenstates at MT sphere surfaces (makusq)
    !i   nsite,isite,iclsl,iclsn see suclst
    !o Outputs:
    !o   weights for each channel output in iomoms style
    !NOTE: lmxa=-1 -> no augmentation
    ! ----------------------------------------------------------------------
    implicit none
    integer:: nlmax,ndimh,nq,nsp,nspc,nsite, isite(nsite),iclsl(nsite),iclsn(nsite),nbandmx
    real(8):: ef,evl(nbandmx,nspx,nq)
    integer:: ifi,isp,ib,is,lcls,ncls,i,j,nr,lmxa,iq,nlma,igets,igetss,i1mach,nfstg,nchan
    integer:: nbandx,nspx,npts,ifid,ikp,ichib,ifdos,ie,ild, lh(10)
    real(8),allocatable :: rofi_rv(:), ul_rv(:), sl_rv(:), gz_rv(:),ruu_rv(:), rus_rv(:),rss_rv(:), g_rv(:), s_rv(:,:,:)
    real(8),pointer:: pnu(:,:),pnz(:,:)
    real(8):: a,rmt,z,xx,rsml(n0),ehl(n0), ume(0:lmxax,nsp,nsite),sme(0:lmxax,nsp,nsite),emin,dosw(2),emax,del
    character clabl*8
    logical:: lidos
    real(8),allocatable:: wk_rv(:),dos_rv(:,:,:)
    complex(8):: aus(nlmax,nbandmx,3,nsp,nsite,nq)
    call tcn ('vcdmel')
    rsml=0d0
    ehl=0d0 
    do  i = 1, nsite
       ib = isite(i)
       ncls = iclsn(i)
       lcls = iclsl(i)
       is = ispec(ib) 
       pnu=>pnuall(:,:,ib)
       pnz=>pnzall(:,:,ib)
       clabl=slabl(is) 
       a= spec_a(is)
       nr=nr_i(is)
       rmt=rmt_i(is)
       z=z_i(is)
       lmxa=lmxa_i(is)
       if (lmxa > lmxax) call rxi('vcdmel needs lmxax ',lmxa)
       if (lmxa == -1) cycle 
       allocate(rofi_rv(nr))
       call radmsh( rmt,a,nr,rofi_rv )
       allocate( ul_rv(nr*(lmxa+1)*nsp)) !   --- Augmented wave functions u,s
       allocate( sl_rv(nr*(lmxa+1)*nsp))
       allocate( gz_rv(nr*(lmxa+1)*nsp))
       allocate(ruu_rv(nr*(lmxa+1)*2*nsp))
       allocate(rus_rv(nr*(lmxa+1)*2*nsp))
       allocate(rss_rv(nr*(lmxa+1)*2*nsp))
       rsml= rsmlss(:,is)
       ehl = ehlss(:,is)
       call makusp(n0,z,nsp,rmt,lmxa,v0pot(ib)%v,a,nr,pnu,pnz,rsml,ehl,ul_rv,sl_rv,gz_rv, ruu_rv,rus_rv,rss_rv)
       write(6,"('CLS atom: ib name n l=',i5,' ',a,2i5)") ib, trim(clabl),ncls,lcls !Matrix elements of u,s with core
       allocate(g_rv(nr*2))
       call pvcdm1(ncls,lcls,g_rv,z,lmxa,v0pot(ib)%v,a,nr,rofi_rv,ul_rv,sl_rv,nsp,lmxax,ume(0,1,i ),sme ( 0,1,i ))
       deallocate(g_rv,rss_rv,rus_rv,ruu_rv,gz_rv,sl_rv,ul_rv,rofi_rv)
    enddo
    allocate( s_rv(3*nsite*2*ndimh,nsp,nq) ) ! Open CLS weights file and write first line
    nfstg = 11
    nchan = 3*nsite !this means only half of 3*nsite*2. S
    s_rv=0d0
    do   iq = 1, nq !For each qp, make <nk|x,y,z|core> at each site and save to disk in
       do  isp = 1, nsp
          do  i = 1, nsite
             lcls = iclsl(i)
             ib   = isite(i)
             is   = ispec(ib) 
             lmxa = lmxa_i(is)
             nlma = (lmxa+1)**2
             if(lmxa>-1) call pvcdm2(i,nsite,nbandmx,ndimh,nlma,nlmax,aus(1,1,1,isp,i,iq),ume ( 0,isp,i ),sme ( 0,isp,i ),lcls ,&
                     rv_a_ocg,iv_a_ojcg,iv_a_oidxcg, s_rv(1,isp,iq))
          enddo
          s_rv(:,isp,iq)= 1d2*s_rv(:,isp,iq) !Scale weights arbitrarily by 100 for plotting etc ..
       enddo
    enddo
!    nbandx = ndimh*nspc       !q-dependent only for no PW case.
!    nspx = nsp / nspc
    lidos=.false.
    npts= bz_ndos
    allocate(wk_rv(npts))
    emin = dosw(1)
    emax = bz_dosmax + ef
    allocate(dos_rv(npts,nsp,nchan))
    call dostet ( nbandmx,nsp,nspx,nchan,bz_nabc(1),bz_nabc(2),bz_nabc(3), bz_ntet,iv_a_oidtet, evl &
        ,s_rv(1:nchan*nbandmx,1:nsp,1:nq), npts, emin, emax, lidos,wk_rv,dos_rv )
    del = 0d0
    open(newunit=ifdos,file='dos-vcdmel.'//trim(sname))
    write(ifdos,"(2f10.5,3i5,2f10.5,i5)") emin,emax,npts,nchan,nsp,ef,del
    do ie=1,npts
       write(ifdos,"(15f14.6)") emin+(emax-emin)/(npts-1d0)*(ie-1),((dos_rv(ie,isp,ild),ild=1,nchan),isp=1,nsp)
    enddo
    if (allocated(s_rv)) deallocate(s_rv)
    call tcx ('vcdmel')
  end subroutine vcdmel
  subroutine pvcdm1(ncls,lcls,gcore,z,lmxa,v,a,nr,rofi,ul,sl,nsp,lmxax,ume,sme) !- Radial matrix elements < (u,s) | r | core >
    use m_lgunit,only:stdo
    use m_rseq,only: rseq
    use m_ftox
    implicit none
    integer::ncls,lcls,lmxa,nr,nsp,lmxax,lx, nodes,l,nre,isp,ir,i1mach,iprint
    real(8):: a,z,gcore(nr,2),rofi(nr),v(nr,nsp),ul(nr,0:lmxa,nsp),sl(nr,0:lmxa,nsp),ume(0:lmxax,nsp),sme(0:lmxax,nsp),&
         e1,e2,slo,val,rmax,b,ecore,tol,yyy,dlml,slo1,r,wgt,uc,sc,ecor0,sum
    do  isp = 1, nsp
       if (nsp .eq. 2.and.iprint()>30) write(stdo,ftox)' Spin',isp,'..'
       !   --- gcore <- core level wave function * r ---
       tol = 1.d-8
       e1 = -2.5d0*z*z - 5
       e2 = 20.d0
       val = 1.d-30
       slo = -val
       l = lcls
       rmax = rofi(nr)
       b = rmax/(dexp(a*nr-a)-1.d0)
       nodes = ncls - (l+1)
       ecore = (e1+e2)/2
       call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,v(1,isp),gcore,sum, a,b,rofi,nr,nre)
       ecor0 = ecore
       !   ... Correct core energy by using hankel bc's
       yyy = ecore - v(nr,isp) + 2*z/rmax
       if(nre .eq. nr .and. yyy .lt. 0.d0) then
          dlml = -1.d0-dsqrt(-yyy)*rmax
          do  lx = 1, l
             dlml = -yyy*rmax*rmax/dlml - (2*lx+1)
          enddo
          slo1 = val*(dlml+l+1)/rmax
          call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo1,v(1,isp),gcore, sum,a,b,rofi,nr,nre)
       endif
       write(6,"('vcdmel: ecor0=',f15.8,' ecore=',f15.8)")ecor0,ecore
       write(6,"('(not including electrostatic potential shift)')")
       ! --- Matrix elements < (u,s) | r | core > ---
       print 332
332    format( '   l',3x,'<u|core>',5x,'<s|core>',4x,'<u|r|core>',2x,'<s|r|core>')
       do  l = 0, lmxa
          ume(l,isp) = 0
          sme(l,isp) = 0
          uc = 0
          sc = 0
          do  ir = 2, nre-1
             r = rofi(ir)
             wgt = (mod(ir+1,2)+1) * (r+b)
             uc     = uc + wgt * ul(ir,l,isp) * gcore(ir,1)
             sc     = sc + wgt * sl(ir,l,isp) * gcore(ir,1)
             ume(l,isp) = ume(l,isp) + wgt * ul(ir,l,isp) * r * gcore(ir,1)
             sme(l,isp) = sme(l,isp) + wgt * sl(ir,l,isp) * r * gcore(ir,1)
          enddo
          ir = nre
          r = rofi(ir)
          wgt = .5d0 * (r+b)
          uc     = uc + wgt * ul(ir,l,isp) * gcore(ir,1)
          sc     = sc + wgt * sl(ir,l,isp) * gcore(ir,1)
          ume(l,isp) = ume(l,isp) + wgt * ul(ir,l,isp) * r * gcore(ir,1)
          sme(l,isp) = sme(l,isp) + wgt * sl(ir,l,isp) * r * gcore(ir,1)
          uc = uc*2d0*a/3d0
          sc = sc*2d0*a/3d0
          ume(l,isp) = ume(l,isp)*2d0*a/3d0
          sme(l,isp) = sme(l,isp)*2d0*a/3d0
          print 335, l,uc,sc,ume(l,isp),sme(l,isp)
335       format(i4,4f12.6)
       enddo
    enddo
  end subroutine pvcdm1
  subroutine pvcdm2(isite,nsite,nbandmx,ndimh,nlma,nlmax,aus,ume,sme,lcls,cg,jcg,indxcg,s) !Kernel called by vcmdel
    use m_ll,only:ll 
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   indxcg:index for !lebsch Gordon coefficients
    !o Outputs
    !o   s     :Matrix elements
    implicit none
    integer isite,lcls,nbandmx,ndimh,nlma,nlmax,nsite,indxcg(*),jcg(*)
    double precision cg(*),ume(0:*),sme(0:*),s(3,nsite,ndimh,2)
    double complex aus(nlmax,nbandmx,2)
    integer kk(4),mlm,lm,klm,ii,indx,icg1,icg2,icg,llm,ib,k
    double complex cxx     !     Transposes (y,z,x) to (x,y,z)
    data kk /0,2,3,1/
    do  11  mlm = 1, nlma !Loop over lm of (u,s)
       lm = ll(mlm)
       !  Selection rule would be handled by CG anyway:
       if(lm .eq. lcls-1 .or. lm .eq. lcls+1) then
          do 14  klm = 2, 4
             ii = max0(mlm,klm)
             indx = (ii*(ii-1))/2 + min0(mlm,klm)
             icg1 = indxcg(indx)
             icg2 = indxcg(indx+1) - 1
             do  15  icg = icg1, icg2
                llm  = jcg(icg)
                if (ll(llm) .eq. lcls) then !lm of core
                   do  10  ib = 1, ndimh
                      cxx = cg(icg)* (dconjg(aus(mlm,ib,1))*ume(lm) + dconjg(aus(mlm,ib,2))*sme(lm))
                      s(kk(klm),isite,ib,1) = s(kk(klm),isite,ib,1) + dble(cxx)
                      s(kk(klm),isite,ib,2) = s(kk(klm),isite,ib,2) + dimag(cxx)
10                 enddo
                endif
15           enddo
14        enddo
       endif
11  enddo
    do  k = 1, 3
       do  ib = 1, ndimh
          s(k,isite,ib,1) = s(k,isite,ib,1)*s(k,isite,ib,1) + s(k,isite,ib,2)*s(k,isite,ib,2)
       enddo
    enddo
  end subroutine pvcdm2
  subroutine dostet(nbmx,nsp,nspx,nchan,n1,n2,n3,ntet,idtet,eband,doswt,npts,emin,emax,lidos,wk,zos)! Density of states to third order by tetrahedron method
    !i Inputs:
    !i   nbmx, first dim. of eband;
    !i   nsp=1 spin degenerate, =2 non-deg;
    !i   nspx: 1 for spin up and down coupled, otherwise nsp
    !i   nbandx, no. of bands;
    !i   nchan, no. of DOS channels; n1,n2,n3;
    !i   ntet, idtet, o/p from tetirr
    !i   eband, work array for bands; doswt, work array for weights;
    !i   npts, no of points in energy range: [emin, emax];
    !i   lidos :F zos = dos
    !i         :T zos = energy integral of dos
    !o Outputs: zos : DOS (or integrated DOS for lidos=T) for each spin and nchan
    implicit none
    integer :: nchan,nsp,nspx,nbmx,npts,ntet,idtet(0:4,*),n1,n2,n3,isp,ib,i,itet,ichan,iq1234(4),nspc,jsp,ksp
    real(8) :: eband(nbmx,nspx,*),emin,emax,wk(npts),zos(npts,nsp,nchan),doswt(nchan,nbmx,nsp,*),bin,eigen(4),v,wt,ebot,dmin1
    logical :: lidos
    if (npts <= 1 .OR. npts <= 2 .AND. .NOT. lidos) call rx1('dostet: npts(=%i) too small for DOS : require npts>2',npts)
    nspc = nsp / nspx
    zos=0d0
    bin = (emax - emin) / (npts - 1)
    v = ( 3d0  -  nsp ) / ( n1 * n2 * n3 * 6d0 )
    tetrahdronloop: do  5  itet = 1, ntet
       iq1234 = idtet(1:4,itet)
       do  4  isp = 1, nspx !Loop over spins and sum over bands ---
          do  3  ib = 1, nbmx !nbandx
             eigen = eband(ib,isp,iq1234)
             if (minval(eigen) > emax) cycle
             do  jsp = 1, nspc
                ksp = max(jsp,isp) !ksp is isp for uncoupled spins, and jsp for coupled spins
                do  ichan = 1, nchan !       ... Accumulate no. states assuming constant wt from this tet
                   wt = sum( doswt(ichan,ib,ksp,iq1234) ) * idtet(0,itet) * v / 4d0
                   call slinz(wt,eigen,emin,emax,zos(1,ksp,ichan),npts)
                enddo
             enddo
3         enddo
4      enddo
5   enddo tetrahdronloop
    if (lidos) return
    do 111  isp  = 1, nsp !DOS from finite difference of NOS ---
       do ichan = 1, nchan
          zos(2:,  isp,ichan) = [((zos(i+1,isp,ichan) - zos(i-1,isp,ichan)) /bin/2d0, i=2,npts-1)]
          zos(1,   isp,ichan) = zos(2,isp,ichan)
          zos(npts,isp,ichan) = zos(npts-1,isp,ichan)
       enddo
111 enddo
  end subroutine dostet
end module m_vcdmel
