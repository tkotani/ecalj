      subroutine ovlocr ( nbas , ssite , sspec ,  nxi0 , nxi 
     .  , exi , hfc , rsmfa , rv_a_orhofa , sv_p_orhoat , sqloc, slmom )
      use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,nsp
      use m_struc_def
      use m_lgunit,only:stdo
C- Makes the site densities for overlapped free atoms.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Ci   slat  :struct containing information about the lattice
Ci   nxi   :number of Hankels
Ci   nxi0  :leading dimension of hfc
Ci   exi   :smoothed Hankel energies; see Remarks
Ci   hfc   :coefficients to smoothed Hankels
Ci   rsmfa :Hankel smoothing radius
Ci   orhofa:free-atom density, by species
Co Outputs
C   orhoat :local density, given by true and smooth densities
C    sqloc :sum of local charges (integral over rho1-rho2)
C    slmom :sum of local magnetic moments
Cr Remarks
Cu Updates
Cu   12 May 07 parallelized (MPI)
Cu   01 Jul 05 Zero-radius empty spheres skip as having no local part
Cu   14 Jun 00 spin polarized
Cu   24 Apr 00 Adapted from nfp ovlocr.f
C ----------------------------------------------------------------------
      implicit none
      integer:: nbas , nxi(1) , nxi0
      type(s_rv1) :: sv_p_orhoat(3,nbas)
      type(s_rv1) :: rv_a_orhofa(nbas)
      real(8):: rsmfa(1) , exi(nxi0,1) , hfc(nxi0,2,1) , sqloc , slmom
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
c      parameter(nrmx=1501,nlmx=64,ktop0=25)
      integer:: ib , ipr , iprint , is , jb , je , js , kmxv , lfoca 
     .    , lmxl , nlmh , nlml , nr , i
      double precision ceh,cofg,cofh,eh,qcorg,qcorh,qsc,qcsm,qloc,rfoca,
     .rmt,rsmh,rsmv,z,amom
      double precision a,p1(3),
     .p2(3),q(3) !,b0(ktop0+1),acof((ktop0+1),nlmx,2)
      real(8),allocatable:: acof(:,:,:),b0(:,:),rofi(:),rwgt(:)
      complex(8),allocatable:: b(:,:)
c$$$#if MPI | MPIK
c$$$C      include "mpif.h"
c$$$      integer, dimension(:),allocatable :: ibproc
c$$$      integer ierr,pid
c$$$      integer procid, master, mpipid, numprocs
c$$$      logical mlog,cmdopt
c$$$      character strn*120
c$$$#endif
      data q /0d0,0d0,0d0/
      integer:: ibini,ibend

      call tcn('ovlocr')
      ipr  = iprint()
c      stdo = lgunit(1)
c      nsp  = globalvariables%nsp
c$$$#if MPI | MPIK
c$$$      procid = mpipid(1)
c$$$      numprocs = mpipid(0)
c$$$      master = 0
c$$$      mlog = cmdopt('--mlog',6,0,strn)
c$$$#endif
C --- Loop over sites where density is made ---
      sqloc = 0
      slmom = 0
c$$$#if MPI | MPIK
c$$$      call info0(31,1,0,' ovlocr: make sphere densities from '//
c$$$     .'overlapping FA densities (parallel)')
c$$$      allocate (ibproc(0:numprocs), stat=ierr)
c$$$      call pshpr(ipr-10)
c$$$      call dstrbp(nbas,numprocs,1,ibproc(0))
c$$$      call poppr
c$$$      ipr = 0
c$$$      call pshpr(10)
c$$$      ibini= ibproc(procid)
c$$$      ibend= ibproc(procid+1)-1
c$$$#else
      if (ipr .ge. 30) write (stdo,300)
  300 format(/' Free atom and overlapped crystal site charges:'
     ./'   ib    true(FA)    smooth(FA)  true(OV)',
     .'    smooth(OV)    local')
      ibini= 1
      ibend= nbas
c$$$#endif
      do  ib = ibini,ibend
        is=ssite(ib)%spec
        p1(:)=ssite(ib)%pos(:)
        lmxl=sspec(is)%lmxl
        kmxv=sspec(is)%kmxv
        rsmv=sspec(is)%rsmv
        nlml = (lmxl+1)**2
        allocate(acof(0:kmxv,nlml,nsp),b(0:kmxv,nlml))
        acof=0d0
        b=0d0
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoca,rfoca, z)
        qcsm = qcorg+qcorh
        if (lmxl .eq. -1) goto 10
        allocate(rofi(nr),rwgt(nr))
        call radmsh(rmt,a,nr,rofi)
        call radwgt(rmt,a,nr,rwgt)
C   ... Loop over other sites, add up tail expansion
        do  jb = 1, nbas
          js=ssite(jb)%spec
          p2(:)=ssite(jb)%pos(:)
          do  je = 1, nxi(js)
            rsmh = rsmfa(js)
            eh   = exi(je,js)
            nlmh = 1
            call hxpbl ( p2 , p1 , q , rsmh , rsmv , eh , kmxv , nlmh , nlml
     .          , kmxv , nlml , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy 
     .          ,  b ) !slat ,

            allocate(b0(0:kmxv,nlmh))
            b0=0d0
c            call dpzero(b0,(kmxv+1)*nlmh)
            if (ib .eq. jb) call hxpos(rsmh,rsmv,eh,kmxv,nlmh,kmxv,b0)
            do  i = 1, nsp
              call p1ovlc(kmxv,nlml,hfc(je,i,js),b,b0,acof(0,1,i))
            enddo
            deallocate(b0)
          enddo
        enddo
        call p2ovlc ( ib , nsp , rsmv , kmxv , nr , nlml , acof , rofi
     .    , rwgt , nxi0 , nxi ( is ) , exi ( 1 , is ) , hfc ( 1 
     .    , 1 , is ) , rsmfa ( is ) , rv_a_orhofa ( is ) %v , sv_p_orhoat( 3 , ib )%v 
     .    , lfoca , qcsm , qloc , amom , sv_p_orhoat( 1 , ib )%v , sv_p_orhoat( 2 , ib )%v 
     .    )
        sqloc = sqloc + qloc
        slmom = slmom + amom
        deallocate(rofi,rwgt)
   10   continue
        deallocate(acof,b)
      enddo

c$$$#if MPI | MPIK
c$$$      call poppr
c$$$
c$$$C ... Combine sphere densities from separate threads
c$$$      do  pid = 0, numprocs-1
c$$$        do  ib = ibproc(pid), ibproc(pid+1)-1
c$$$
c$$$          is=ssite(ib)%spec
c$$$
c$$$
c$$$          nr=sspec(is)%nr
c$$$          lmxl=sspec(is)%lmxl
c$$$
c$$$          nlml = (lmxl+1)**2
c$$$          call mpibc3 ( sv_p_orhoat( 1 , ib )%v , nr * nlml * nsp , 4 , 
c$$$     .     pid , 0 , 'ovlocr' , 'rhoat(1)' )
c$$$
c$$$          call mpibc3 ( sv_p_orhoat( 2 , ib )%v , nr * nlml * nsp , 4 , 
c$$$     .     pid , 0 , 'ovlocr' , 'rhoat(1)' )
c$$$
c$$$          call mpibc3 ( sv_p_orhoat( 3 , ib )%v , nr * nsp , 4 , pid , 
c$$$     .     0 , 'ovlocr' , 'rhoat(1)' )
c$$$
c$$$C          call MPI_BCAST(w(orhoat(1,ib)),nr*nlml*nsp,
c$$$C     .      MPI_DOUBLE_PRECISION,pid,MPI_COMM_WORLD,ierr)
c$$$C          call MPI_BCAST(w(orhoat(2,ib)),nr*nlml*nsp,
c$$$C     .      MPI_DOUBLE_PRECISION,pid,MPI_COMM_WORLD,ierr)
c$$$C          call MPI_BCAST(w(orhoat(3,ib)),nr*nsp,
c$$$C     .      MPI_DOUBLE_PRECISION,pid,MPI_COMM_WORLD,ierr)
c$$$        enddo
c$$$      enddo
c$$$C ... Combine sum-of-sphere charges from separate threads
c$$$      call mpibc2_real(sqloc,1,'ovlocr_sqloc')
c$$$      call mpibc2_real(slmom,1,'ovlocr_slmom')
c$$$
c$$$C     Debugging printout
c$$$C      print 345, procid, (dval(w(orhoat(1,ib)),445), ib=1,20)
c$$$C      print 345, procid, (dval(w(orhoat(2,ib)),445), ib=1,20)
c$$$C      print 345, procid, (dval(w(orhoat(3,ib)),445), ib=1,20)
c$$$C  345 format(i4,20f10.5)
c$$$C     print *, procid, iprint()
c$$$
c$$$      deallocate(ibproc, stat=ierr)
c$$$C     call rx0('done')
c$$$#endif
c$$$
      call tcx('ovlocr')
      end subroutine ovlocr


      subroutine p1ovlc(kmxv,nlml,hfc,b,b0,a)

C- Adds contribution to P_kl expansion of density from one basis function
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlml  :density expanded to nlml
Ci   kmxv  :k-cutoff for P_kl expansion
Ci   hfc   :coefficient to basis function
Ci   b     :P_kl expansion of density from one basis function
Ci   b0    :P_kl expansion of on-site density
Co Outputs
Co   a     :cumulative P_kl expansion of density for this site
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlml,kmxv
      double precision a(0:kmxv,nlml),b0(0:kmxv,1),hfc
      double complex b(0:kmxv,nlml)
C ... Local parameters
      integer k,ilm

      do  10  k = 0, kmxv
        do  12  ilm = 1, nlml
          a(k,ilm) = a(k,ilm) + hfc*dble(b(k,ilm))
   12   continue
        a(k,1) = a(k,1) - hfc*b0(k,1)
   10 continue
      end subroutine p1ovlc


      subroutine p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,acof,rofi,rwgt, !rhohd,
     .nxi0,nxi,exi,hfc,rsmfa,rhofa,rhoc,lfoca,qcsm,qloc,amom,rho1,rho2)
      use m_lgunit,only:stdo
C- Assemble local density from P_kl expansion for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ib    :site for which to assemble local density
Ci   nsp   :number of spin channels
Ci   rsmv  :smoothing radius for P_kl expansion
Ci   kmxv  :k-cutoff for P_kl expansion
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for P_kl expansion
Ci   acof  :coefficients to P_kl expansion
Ci   rofi  :radial mesh points for tabulation on a radial mesh
Ci   rwgt  :radial mesh weights for integration on a radial mesh
Co   rhohd :work array (holds on-site smoothed head density)
Ci   nxi0  :leading dimension of hfc
Ci   nxi   :number of smoothed Hankel energies in head expansion
Ci   exi   :smoothed Hankel energies in head expansion
Ci   hfc   :coefficients to Hankel energies in head expansion
Ci   rsmfa :Hankel smoothing radius in head expansion
Ci   rhofa :head free-atom density
Ci   rhoc  :core density --- used to integrate core charge
Ci   lfoca :switch specifying treatment of core density.
Ci          0 => val,slo = 0 at sphere boundary
Ci          1 => core tails included explicitly with valence
Ci          2 => tails included perturbatively
Ci   qcsm  :smoothed core density
Co Outputs
Ci   rho1  :local true density, tabulated on a radial mesh
Ci   rho2  :local smoothed density, tabulated on a radial mesh
Co   qloc  :sphere charge
Co   amom  :sphere magnetic moment
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nxi0,ib,nsp,kmxv,nlml,nxi,lfoca
      double precision qcsm,qloc,rhofa(nr,nsp),rho1(nr,nlml,nsp),
     .rho2(nr,nlml,nsp),rofi(nr),rwgt(nr),rhohd(nr,nsp),exi(1),
     .hfc(nxi0,nsp),rhoc(nr,nsp),acof(0:kmxv,nlml,nsp),rsmv,rsmfa,amom
C ... Local parameters
c      integer stdo !kmx,lmx,
c      parameter (kmx=20, lmx=6)
      integer i,ie,ilm,ipr,iprint,k,l,ll,lmax,lmxl,isp
      double precision asm,gam,pi,qall,qexa,qin,qlc,qnum,qout,qsmo,qut,
     .r,rl,rmt,srfpi,sum,sumfa,sumhd,sumsm,sumtr,y0,
     .xi(0:10),x0(0:2),ddot !pkl(0:kmx,0:lmx)
      real(8),allocatable:: pkl(:,:)

      ipr   = iprint()
c      stdo  = lgunit(1)
      pi    = 4d0*datan(1d0)
      srfpi = dsqrt(4*pi)
      y0    = 1d0/srfpi
      lmxl  = ll(nlml)
      allocate(pkl(0:kmxv,0:lmxl))
c      if (lmxl .gt. lmx) call rxi('ovlocr: increase lmx, need',lmxl)

C     do  ilm = 1, nlml
C       do  k = 0, kmxv
C         if (dabs(acof(k,ilm,1)).gt.1d-6)
C    .      write(stdo,780) ilm,k,acof(k,ilm,1),acof(k,ilm,nsp)
C 780     format('ilm,k',2i5,2f14.8)
C       enddo
C     enddo

C --- Assemble smooth on-site head density in rhohd ---
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
        qout = srfpi/exi(ie)*(-dexp(rsmfa**2/4*exi(ie))
     .  - rmt**3*(xi(1)-dexp(rsmfa**2/4*exi(ie))*x0(1)))
        qin = qall-qout
        do  isp = 1, nsp
          qnum = qnum + hfc(ie,isp)*sum
          qexa = qexa + hfc(ie,isp)*qin
          qsmo = qsmo + hfc(ie,isp)*qall
          qut  = qut  + hfc(ie,isp)*qout
        enddo
      enddo

c|      write(stdo,917) qnum,qexa,qsmo,qut
c|  917 format('summed smooth charge:  num',f14.8,'   exact',f14.8
c|     .   /' total smooth q',f14.8,'  outside',f14.8)

C --- Assemble overlapped tail density in rho2 ---
c      if (kmxv .gt. kmx) call rx('ovlocr: increase kmx')
      call dpzero(rho2,  nr*nlml*nsp)
      do  i = 1, nr
        r = rofi(i)
        call radpkl(r,rsmv,kmxv,lmxl,kmxv,pkl)
        do  isp = 1, nsp
          do  ilm = 1, nlml
            l = ll(ilm)
            rl = 0.d0
            if ( r .gt. 0.d0 ) rl = r**l
            do  k = 0, kmxv
              rho2(i,ilm,isp) = rho2(i,ilm,isp) +
     .        acof(k,ilm,isp)*pkl(k,l)*r*r*rl
            enddo
          enddo
        enddo
      enddo

C ... Make the true density in rho1, smooth density in rho2
      call dpcopy(rho2,rho1,1,nr*nlml*nsp,1d0)
      do   isp = 1, nsp
      do   i = 1, nr
        rho1(i,1,isp) = rho1(i,1,isp) + y0*rhofa(i,isp)
        rho2(i,1,isp) = rho2(i,1,isp) + y0*rhohd(i,isp)
      enddo
      enddo
      
C ... Do some integrals
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
      amom = -srfpi*
     .(ddot(nr,rwgt,1,rho1(1,1,nsp),1)-ddot(nr,rwgt,1,rho1(1,1,1),1)
     .-ddot(nr,rwgt,1,rho2(1,1,nsp),1)+ddot(nr,rwgt,1,rho2(1,1,1),1))

      if (lfoca .eq. 0) qloc = qloc + qlc - qcsm
      if (ipr .ge. 30) then
        write(stdo,810) ib,sumfa,sumhd,sumtr,sumsm,qloc
        if (nsp .eq. 2) write(stdo,811)
     .  ddot(nr,rwgt,1,rhofa,1)-ddot(nr,rwgt,1,rhofa(1,2),1),
     .  ddot(nr,rwgt,1,rhohd,1)-ddot(nr,rwgt,1,rhohd(1,2),1),
     .  srfpi*(ddot(nr,rwgt,1,rho1,1)-ddot(nr,rwgt,1,rho1(1,1,2),1)),
     .  srfpi*(ddot(nr,rwgt,1,rho2,1)-ddot(nr,rwgt,1,rho2(1,1,2),1)),
     .  amom

      endif
  810 format(i5,6f12.6)
  811 format(' amom',6f12.6)

      end subroutine p2ovlc


