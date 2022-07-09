subroutine ovlocr(nbas,nxi0,nxi,exi,hfc,rsmfa,rv_a_orhofa, sv_p_orhoat , sqloc, slmom )
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,nsp,ispec,sspec=>v_sspec
  use m_struc_def
  use m_lgunit,only:stdo
  use m_smhankel,only: hxpbl
  use m_lattic,only: rv_a_opos
  !- Makes the site densities for overlapped free atoms.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct containing site-specific information
  !i   sspec :struct containing species-specific information
  !i   slat  :struct containing information about the lattice
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
  !r Remarks
  !u Updates
  !u   xxx 12 May 07 parallelized (MPI)
  !u   01 Jul 05 Zero-radius empty spheres skip as having no local part
  !u   14 Jun 00 spin polarized
  !u   24 Apr 00 Adapted from nfp ovlocr.f
  ! ----------------------------------------------------------------------
  implicit none
  integer:: kmxv=15 !Hardwired taken from original code.
  
  integer:: nbas , nxi(1) , nxi0
  type(s_rv1) :: sv_p_orhoat(3,nbas)
  type(s_rv1) :: rv_a_orhofa(nbas)
  real(8):: rsmfa(1) , exi(nxi0,1) , hfc(nxi0,2,1) , sqloc , slmom
!  type(s_site)::ssite(*)
!  type(s_spec)::sspec(*)
  integer:: ib , ipr , iprint , is , jb , je , js , lfoca &
       , lmxl , nlmh , nlml , nr , i
  double precision :: ceh,cofg,cofh,eh,qcorg,qcorh,qsc,qcsm,qloc,rfoca,rmt,rsmh,rsmv,z,amom
  double precision :: a,p1(3), p2(3),q(3) !,b0(ktop0+1),acof((ktop0+1),nlmx,2)
  real(8),allocatable:: acof(:,:,:),b0(:,:),rofi(:),rwgt(:)
  complex(8),allocatable:: b(:,:)
  data q /0d0,0d0,0d0/
  integer:: ibini,ibend

  call tcn('ovlocr')
  ipr  = iprint()
  sqloc = 0
  slmom = 0
  if (ipr >= 30) write (stdo,300)
300 format(/' Free atom and overlapped crystal site charges:' &
       /'   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local')
  ibini= 1
  ibend= nbas
  do  ib = ibini,ibend
     is=ispec(ib) !ssite(ib)%spec
     p1(:)=rv_a_opos(:,ib) !ssite(ib)%pos(:)
     lmxl=sspec(is)%lmxl
     !kmxv=sspec(is)%kmxv
     rsmv=sspec(is)%rsmv
     nlml = (lmxl+1)**2
     allocate(acof(0:kmxv,nlml,nsp),b(0:kmxv,nlml))
     acof=0d0
     b=0d0
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoca,rfoca, z)
     qcsm = qcorg+qcorh
     if (lmxl == -1) goto 10
     allocate(rofi(nr),rwgt(nr))
     call radmsh(rmt,a,nr,rofi)
     call radwgt(rmt,a,nr,rwgt)
     !   ... Loop over other sites, add up tail expansion
     do  jb = 1, nbas
        js=ispec(jb) !ssite(jb)%spec
        p2(:)=rv_a_opos(:,jb) !ssite(jb)%pos(:)
        do  je = 1, nxi(js)
           rsmh = rsmfa(js)
           eh   = exi(je,js)
           nlmh = 1
           call hxpbl ( p2 , p1 , q , [rsmh], rsmv , [eh] , kmxv , nlmh , nlml &
                , kmxv , nlml , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy ,  b ) 
           allocate(b0(0:kmxv,nlmh))
           b0=0d0
           if (ib == jb) call hxpos([rsmh],rsmv,[eh],kmxv,nlmh,kmxv,b0)
           do  i = 1, nsp
              call p1ovlc(kmxv,nlml,hfc(je,i,js),b,b0,acof(0,1,i))
           enddo
           deallocate(b0)
        enddo
     enddo
     call p2ovlc ( ib , nsp , rsmv , kmxv , nr , nlml , acof , rofi &
          , rwgt , nxi0 , nxi ( is ) , exi ( 1 , is ) , hfc ( 1 &
          , 1 , is ) , rsmfa ( is ) , rv_a_orhofa ( is ) %v , sv_p_orhoat( 3 , ib )%v &
          , lfoca , qcsm , qloc , amom , sv_p_orhoat( 1 , ib )%v , sv_p_orhoat( 2 , ib )%v )
     sqloc = sqloc + qloc
     slmom = slmom + amom
     deallocate(rofi,rwgt)
10   continue
     deallocate(acof,b)
  enddo
  call tcx('ovlocr')
end subroutine ovlocr


subroutine p1ovlc(kmxv,nlml,hfc,b,b0,a)
  !- Adds contribution to P_kl expansion of density from one basis function
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nlml  :density expanded to nlml
  !i   kmxv  :k-cutoff for P_kl expansion
  !i   hfc   :coefficient to basis function
  !i   b     :P_kl expansion of density from one basis function
  !i   b0    :P_kl expansion of on-site density
  !o Outputs
  !o   a     :cumulative P_kl expansion of density for this site
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nlml,kmxv
  double precision :: a(0:kmxv,nlml),b0(0:kmxv,1),hfc
  double complex b(0:kmxv,nlml)
  integer :: k,ilm
  do  10  k = 0, kmxv
     do  12  ilm = 1, nlml
        a(k,ilm) = a(k,ilm) + hfc*dble(b(k,ilm))
12   enddo
     a(k,1) = a(k,1) - hfc*b0(k,1)
10 enddo
end subroutine p1ovlc


subroutine p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,acof,rofi,rwgt, &
     nxi0,nxi,exi,hfc,rsmfa,rhofa,rhoc,lfoca,qcsm,qloc,amom,rho1,rho2)
  use m_lgunit,only:stdo
  use m_hansr,only: hansmr
  !- Assemble local density from P_kl expansion for one site
  ! ----------------------------------------------------------------------
  !i Inputs
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
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nr,nxi0,ib,nsp,kmxv,nlml,nxi,lfoca
  double precision :: qcsm,qloc,rhofa(nr,nsp),rho1(nr,nlml,nsp), &
       rho2(nr,nlml,nsp),rofi(nr),rwgt(nr),rhohd(nr,nsp),exi(1), &
       hfc(nxi0,nsp),rhoc(nr,nsp),acof(0:kmxv,nlml,nsp),rsmv,rsmfa,amom
  ! ... Local parameters
  !      integer stdo !kmx,lmx,
  !      parameter (kmx=20, lmx=6)
  integer :: i,ie,ilm,ipr,iprint,k,l,ll,lmax,lmxl,isp
  double precision :: asm,gam,pi,qall,qexa,qin,qlc,qnum,qout,qsmo,qut, &
       r,rl,rmt,srfpi,sum,sumfa,sumhd,sumsm,sumtr,y0, &
       xi(0:10),x0(0:2),ddot !pkl(0:kmx,0:lmx)
  real(8),allocatable:: pkl(:,:)

  ipr   = iprint()
  !      stdo  = lgunit(1)
  pi    = 4d0*datan(1d0)
  srfpi = dsqrt(4*pi)
  y0    = 1d0/srfpi
  lmxl  = ll(nlml)
  allocate(pkl(0:kmxv,0:lmxl))
  !      if (lmxl .gt. lmx) call rxi('ovlocr: increase lmx, need',lmxl)

  !     do  ilm = 1, nlml
  !       do  k = 0, kmxv
  !         if (dabs(acof(k,ilm,1)).gt.1d-6)
  !    .      write(stdo,780) ilm,k,acof(k,ilm,1),acof(k,ilm,nsp)
  ! 780     format('ilm,k',2i5,2f14.8)
  !       enddo
  !     enddo

  ! --- Assemble smooth on-site head density in rhohd ---
  qnum = 0d0
  qexa = 0d0
  qsmo = 0d0
  qut = 0d0
  call dpzero(rhohd, nr*nsp)
  asm = 1d0/rsmfa
  lmax = 0
  do  ie = 1, nxi
     sum = 0d0
     do  i = 1, nr
        r = rofi(i)
        call hansmr(r,exi(ie),asm,xi,lmax)
        sum = sum + srfpi*rwgt(i)*xi(0)*r*r
        do  isp = 1, nsp
           rhohd(i,isp) = rhohd(i,isp) + srfpi*hfc(ie,isp)*xi(0)*r*r
        enddo
     enddo
     gam = 0.25d0*rsmfa**2
     qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
     rmt = rofi(nr)
     call hansmr(rmt,0d0,1/rsmfa,x0,1)
     call hansmr(rmt,exi(ie),1/rsmfa,xi,1)
     qout = srfpi/exi(ie)*(-dexp(rsmfa**2/4*exi(ie)) &
          - rmt**3*(xi(1)-dexp(rsmfa**2/4*exi(ie))*x0(1)))
     qin = qall-qout
     do  isp = 1, nsp
        qnum = qnum + hfc(ie,isp)*sum
        qexa = qexa + hfc(ie,isp)*qin
        qsmo = qsmo + hfc(ie,isp)*qall
        qut  = qut  + hfc(ie,isp)*qout
     enddo
  enddo

  !|      write(stdo,917) qnum,qexa,qsmo,qut
  !|  917 format('summed smooth charge:  num',f14.8,'   exact',f14.8
  !|     .   /' total smooth q',f14.8,'  outside',f14.8)

  ! --- Assemble overlapped tail density in rho2 ---
  !      if (kmxv .gt. kmx) call rx('ovlocr: increase kmx')
  call dpzero(rho2,  nr*nlml*nsp)
  do  i = 1, nr
     r = rofi(i)
     call radpkl(r,rsmv,kmxv,lmxl,kmxv,pkl)
     do  isp = 1, nsp
        do  ilm = 1, nlml
           l = ll(ilm)
           rl = 0.d0
           if ( r > 0.d0 ) rl = r**l
           do  k = 0, kmxv
              rho2(i,ilm,isp) = rho2(i,ilm,isp) + &
                   acof(k,ilm,isp)*pkl(k,l)*r*r*rl
           enddo
        enddo
     enddo
  enddo

  ! ... Make the true density in rho1, smooth density in rho2
  call dpcopy(rho2,rho1,1,nr*nlml*nsp,1d0)
  do   isp = 1, nsp
     do   i = 1, nr
        rho1(i,1,isp) = rho1(i,1,isp) + y0*rhofa(i,isp)
        rho2(i,1,isp) = rho2(i,1,isp) + y0*rhohd(i,isp)
     enddo
  enddo

  ! ... Do some integrals
  sumfa = 0d0
  sumsm = 0d0
  sumhd = 0d0
  sumtr = 0d0
  qlc = 0d0
  do   isp = 1, nsp
     do   i = 1, nr
        sumfa = sumfa + rwgt(i)*rhofa(i,isp)
        sumhd = sumhd + rwgt(i)*rhohd(i,isp)
        sumtr = sumtr + rwgt(i)*rho1(i,1,isp)
        sumsm = sumsm + rwgt(i)*rho2(i,1,isp)
        qlc = qlc + rwgt(i)*rhoc(i,isp)
     enddo
  enddo
  sumsm = sumsm*srfpi
  sumtr = sumtr*srfpi
  qloc = sumtr-sumsm
  amom = -srfpi* &
       (ddot(nr,rwgt,1,rho1(1,1,nsp),1)-ddot(nr,rwgt,1,rho1(1,1,1),1) &
       -ddot(nr,rwgt,1,rho2(1,1,nsp),1)+ddot(nr,rwgt,1,rho2(1,1,1),1))

  if (lfoca == 0) qloc = qloc + qlc - qcsm
  if (ipr >= 30) then
     write(stdo,810) ib,sumfa,sumhd,sumtr,sumsm,qloc
     if (nsp == 2) write(stdo,811) &
          ddot(nr,rwgt,1,rhofa,1)-ddot(nr,rwgt,1,rhofa(1,2),1), &
          ddot(nr,rwgt,1,rhohd,1)-ddot(nr,rwgt,1,rhohd(1,2),1), &
          srfpi*(ddot(nr,rwgt,1,rho1,1)-ddot(nr,rwgt,1,rho1(1,1,2),1)), &
          srfpi*(ddot(nr,rwgt,1,rho2,1)-ddot(nr,rwgt,1,rho2(1,1,2),1)), &
          amom

  endif
810 format(i5,6f12.6)
811 format(' amom',6f12.6)

end subroutine p2ovlc


