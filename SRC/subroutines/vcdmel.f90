module m_vcdmel
  public vcdmel
  private
contains  
subroutine vcdmel(nl,nlmax,ndham,ndimh,& !- Valence-core dipole matrix elements
     nq,nsp,nspc,ef,evl,aus,nsite,isite,iclsl,iclsn,dosw)
  use m_lmfinit,only: rv_a_ocg , iv_a_ojcg , iv_a_oidxcg,ispec,sspec=>v_sspec
  use m_mkqp,only: iv_a_oidtet ,bz_nabc, bz_ntet
  use m_struc_def
  use m_ext,only: sname     !extention for file
  use m_lmfinit,only: bz_ndos,bz_dosmax,slabl
  use m_elocp,only: rsmlss=>rsml,ehlss=>ehl
  use m_density,only: v0pot,pnuall,pnzall
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   nlmax :first dimension of aus; largest augmentation (l+1)^2
  !i   ndham :second dimension of aus, at least as large as ndimh
  !i   ndimh :number of eigenvalues
  !i   nq    :number of k-points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
  !i   ef    :Fermi energy
  !i   evl   :energy bands at the nq k-points
  !i   aus   :values and slopes of eigenstates at MT sphere surfaces
  !i          (makusq)
  !i   nsite,isite,iclsl,iclsn see suclst
  !o Outputs:
  !o   weights for each channel output in iomoms style
  !r Remarks
  !u Updates
  !u   08 Jul 08 Dimension aus separately from ndimh
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   19 Sep 03 (ATP) Bug fixes
  !u   28 Mar 01 (MvS) rearrangement of indices to accommodate new makusq.
  !u   20 Mar 01 (ATP) extended to handle multiple core levels
  !u   20 Feb 01 Written by ATP
  ! ----------------------------------------------------------------------
  implicit none
  integer nlmax,ndham,ndimh,nq,nsp,nspc,nsite
  integer isite(nsite),iclsl(nsite),iclsn(nsite)
  real(8):: ef , evl(ndham,nsp,nq)
!  type(s_site)::ssite(nsite)
!  type(s_spec)::sspec(*)
  double complex aus(nlmax,ndham,3,nsp,nsite,nq)
  integer n0,lmxax
  parameter (n0=10,lmxax=10)
  integer ifi,isp,ib,is,lcls,ncls,nl,i,j,nr,lmxa,iq,nlma,igets,igetss,i1mach,nfstg,nchan
  integer lh(10)
  real(8) ,allocatable :: rofi_rv(:)
  real(8) ,allocatable :: ul_rv(:)
  real(8) ,allocatable :: sl_rv(:)
  real(8) ,allocatable :: gz_rv(:)
  real(8) ,allocatable :: ruu_rv(:)
  real(8) ,allocatable :: rus_rv(:)
  real(8) ,allocatable :: rss_rv(:)
  real(8) ,allocatable :: g_rv(:)
  real(8) ,allocatable :: s_rv(:,:,:)
  real(8),pointer:: pnu(:,:),pnz(:,:)!pnu(n0,2),pnz(n0,2)
  real(8):: a,rmt,z,xx,rsml(n0),ehl(n0)
  double precision ume(0:lmxax,nsp,nsite),sme(0:lmxax,nsp,nsite)
  character clabl*8
  integer:: nbandx,nspx,npts,ifid,ikp,ichib,ifile_handle,ifdos,ie,ild
  logical:: lidos
  real(8),allocatable:: wk_rv(:),dos_rv(:,:,:)
  real(8):: emin,dosw(2),emax,del
  call tcn ('vcdmel')
  rsml=0d0
  ehl=0d0 
  do  i = 1, nsite
     ib = isite(i)
     ncls = iclsn(i)
     lcls = iclsl(i)
     is = ispec(ib) !ssite(ib)%spec
     pnu=>pnuall(:,:,ib)
     pnz=>pnzall(:,:,ib)
!     pnu=ssite(ib)%pnu
!     pnz=ssite(ib)%pz
     clabl=slabl(is) !sspec(is)%name
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     z=sspec(is)%z
     lmxa=sspec(is)%lmxa
     if (lmxa .gt. lmxax) call rxi('vcdmel needs lmxax ',lmxa)
     if (lmxa .eq. -1) cycle !  if (lmxa .eq. -1) goto 10
     allocate(rofi_rv(nr))
     call radmsh ( rmt , a , nr , rofi_rv )
     !   --- Augmented wave functions u,s
     allocate(ul_rv(nr*(lmxa+1)*nsp))
     allocate(sl_rv(nr*(lmxa+1)*nsp))
     allocate(gz_rv(nr*(lmxa+1)*nsp))
     allocate(ruu_rv(nr*(lmxa+1)*2*nsp))
     allocate(rus_rv(nr*(lmxa+1)*2*nsp))
     allocate(rss_rv(nr*(lmxa+1)*2*nsp))
     rsml= rsmlss(:,is)
     ehl = ehlss(:,is)
     call makusp ( n0 , z , nsp , rmt , lmxa , v0pot(ib)%v , a , nr , &
          xx , xx , pnu , pnz , rsml , ehl , ul_rv , sl_rv , gz_rv, ruu_rv , rus_rv , rss_rv )
     !   --- Matrix elements of u,s with core
     write(6,"('CLS atom: ib name n l=',i5,' ',a,2i5)") ib, trim(clabl),ncls,lcls
     allocate(g_rv(nr*2))
     call pvcdm1 ( ncls , lcls , g_rv , z , lmxa , v0pot(ib)%v , a &
          , nr , rofi_rv , ul_rv , sl_rv , nsp , lmxax , ume(0,1,i ) , sme ( 0 , 1 , i ) )
     deallocate(g_rv,rss_rv,rus_rv,ruu_rv,gz_rv,sl_rv,ul_rv,rofi_rv)
  enddo
  ! --- Open CLS weights file and write first line
  allocate( s_rv(3*nsite*2*ndimh,nsp,nq) )
  nfstg = 11
  nchan = 3*nsite !this means only half of 3*nsite*2. S
  ! --- For each qp, make <nk|x,y,z|core> at each site and save to disk in
  s_rv=0d0
  do   iq = 1, nq
     do  isp = 1, nsp
        do  i = 1, nsite
           lcls = iclsl(i)
           ib   = isite(i)
           is   = ispec(ib) !ssite(ib)%spec
           lmxa = sspec(is)%lmxa
           nlma = (lmxa+1)**2
           if (lmxa .gt. -1) then
              call pvcdm2(i,nsite,ndham,ndimh,nlma,nlmax,aus(1,1,1,isp,i,iq),&
                   ume ( 0 , isp , i ) , sme ( 0 , isp , i ) , lcls ,&
                   rv_a_ocg , iv_a_ojcg , iv_a_oidxcg, s_rv(1,isp,iq))
           endif
        enddo
        ! --- Scale weights arbitrarily by 100 for plotting etc ..
        s_rv(:,isp,iq)= 1d2*s_rv(:,isp,iq)
        !call dscal ( 3 * ndimh * nsite * 2, 1d2 , s_rv(1,isp,iq) , 1 )
     enddo
  enddo
  nbandx = ndimh*nspc       !q-dependent only for no PW case.
  nspx = nsp / nspc
  lidos=.false.
  npts= bz_ndos
  allocate(wk_rv(npts))
  emin = dosw(1)
  emax = bz_dosmax + ef
  allocate(dos_rv(npts,nsp,nchan))
  call dostet ( ndham , nsp , nspx , nbandx , nchan , &
       bz_nabc(1),bz_nabc(2),bz_nabc(3), bz_ntet , iv_a_oidtet, evl &
       , s_rv(1:nchan*nbandx,1:nsp,1:nq), npts, emin, emax, lidos , wk_rv , dos_rv )
  del = 0d0
  open(newunit=ifdos,file='dos-vcdmel.'//trim(sname))
  write(ifdos,"(2f10.5,3i5,2f10.5,i5)") emin,emax,npts,nchan,nsp,ef,del
  do ie=1,npts
     write(ifdos,"(15f14.6)") &
          emin+(emax-emin)/(npts-1d0)*(ie-1),((dos_rv(ie,isp,ild),ild=1,nchan),isp=1,nsp)
  enddo
  if (allocated(s_rv)) deallocate(s_rv)
  call tcx ('vcdmel')
end subroutine vcdmel
subroutine pvcdm1(ncls,lcls,gcore,z,lmxa,v,a,nr,rofi,ul,sl,nsp,lmxax,ume,sme)
  !- Radial matrix elements < (u,s) | r | core >
  implicit none
  integer ncls,lcls,lmxa,nr,nsp,lmxax
  double precision a,z,gcore(nr,2),rofi(nr),v(nr,nsp),&
       ul(nr,0:lmxa,nsp),sl(nr,0:lmxa,nsp),ume(0:lmxax,nsp),sme(0:lmxax,nsp)
  integer nodes,l,nre,isp,ll,ir,i1mach
  double precision e1,e2,slo,val,rmax,b,ecore,tol,yyy,dlml,slo1,r,wgt,uc,sc,ecor0,sum
  do  isp = 1, nsp
     if (nsp .eq. 2) then
        call info2(30,0,0,' Spin %i ..',isp,0)
     endif
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
        do  ll = 1, l
           dlml = -yyy*rmax*rmax/dlml - (2*ll+1)
        enddo
        slo1 = val*(dlml+l+1)/rmax
        call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo1,v(1,isp),gcore, sum,a,b,rofi,nr,nre)
     endif
     write(6,"('vcdmel: ecor0=',f15.8,' ecore=',f15.8)")ecor0,ecore
     write(6,"('(not including electrostatic potential shift)')")
     ! --- Matrix elements < (u,s) | r | core > ---
     print 332
332  format( '   l',3x,'<u|core>',5x,'<s|core>',4x,'<u|r|core>',2x,'<s|r|core>')
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
335     format(i4,4f12.6)
     enddo
  enddo
end subroutine pvcdm1
subroutine pvcdm2(isite,nsite,ndham,ndimh,nlma,nlmax,aus,ume,sme,lcls,cg,jcg,indxcg,s)
  !- Kernel called by vcmdel
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   isite :
  !i   nsite :
  !i   ndimh :
  !i   nlma  :
  !i   nlmax :
  !i   aus   :
  !i   ume   :
  !i   sme   :
  !i   lcls  :
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   indxcg:index for !lebsch Gordon coefficients
  !o Outputs
  !o   s     :Matrix elements
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  integer isite,lcls,ndham,ndimh,nlma,nlmax,nsite,indxcg(*),jcg(*)
  double precision cg(*),ume(0:*),sme(0:*),s(3,nsite,ndimh,2)
  double complex aus(nlmax,ndham,2)
  integer kk(4),mlm,lm,ll,klm,ii,indx,icg1,icg2,icg,llm,ib,k
  double complex cxx
  !     Transposes (y,z,x) to (x,y,z)
  data kk /0,2,3,1/
  ! ... Loop over lm of (u,s)
  do  11  mlm = 1, nlma
     lm = ll(mlm)
     !       Selection rule would be handled by CG anyway:
     if (lm .eq. lcls-1 .or. lm .eq. lcls+1) then
        do  14  klm = 2, 4
           ii = max0(mlm,klm)
           indx = (ii*(ii-1))/2 + min0(mlm,klm)
           icg1 = indxcg(indx)
           icg2 = indxcg(indx+1) - 1
           do  15  icg = icg1, icg2
              !             lm of core
              llm  = jcg(icg)
              if (ll(llm) .eq. lcls) then
                 do  10  ib = 1, ndimh
                    cxx =  cg(icg)* (dconjg(aus(mlm,ib,1))*ume(lm) + &
                         dconjg(aus(mlm,ib,2))*sme(lm))
                    s(kk(klm),isite,ib,1) = s(kk(klm),isite,ib,1) + dble(cxx)
                    s(kk(klm),isite,ib,2) = s(kk(klm),isite,ib,2) + dimag(cxx)
10               enddo
              endif
15         enddo
14      enddo
     endif
11 enddo
  do  k = 1, 3
     do  ib = 1, ndimh
        s(k,isite,ib,1) = s(k,isite,ib,1)*s(k,isite,ib,1) + s(k,isite,ib,2)*s(k,isite,ib,2)
     enddo
  enddo
end subroutine pvcdm2
end module m_vcdmel
