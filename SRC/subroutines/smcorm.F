      subroutine smcorm(nbas,ssite,sspec,ng,gv,!,slat
     .cgh1,cgh2,lfoc1,lfoc2)
      use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg


      use m_struc_def           !Cgetarg
      use m_lmfinit,only:lat_alat,lat_vol


C- For foca, add together density of smoothed part of core
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:
Ci     Stored:    *
Ci     Passed to: corprm
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat vol ocy
Ci     Stored:    *
Ci     Passed to: *
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Co Outputs
Co   cgh1  :Portion of smoothed core that is treated directly
Co   cgh2  :Portion of smoothed core that is treated perturbatively
Co   lfoc1 :returned nonzero if any site lfoca is direct (1)
Co   lfoc2 :returned nonzero if any site lfoca is perturbative
Cu Updates
Cu   02 Jul 05  skip sites for which cofh=0
Cu    1 May 00  Adapted from nfp smc_mesh.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,nbas,lfoc1,lfoc2
      real(8):: gv(ng,3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
c      type(s_lat)::slat
      double complex cgh1(ng),cgh2(ng)
C ... Local parameters
      integer:: k0 , nlmx , kmax , ib , is , lfoc , i
      double precision tau(3),v(3),alat,vol,qcorg,qcorh,qsc,cofg,cofh,
     .ceh,rfoc,z
      parameter (k0=3, nlmx=25)
      double complex gkl(0:k0,nlmx)
      alat=lat_alat
      vol=lat_vol
      kmax = 0
C --- Accumulate FT of smooth-Hankel foca heads ---
      call dpzero(cgh1,  2*ng)
      call dpzero(cgh2,  2*ng)
      lfoc1 = 0
      lfoc2 = 0
      do  ib = 1, nbas
        is=ssite(ib)%spec
c        i_copy_size=size(ssite(ib)%pos)
c     call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1)
        tau=ssite(ib)%pos
        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
C       qc = qcorg+qcorh
C        if (iprint() .ge. 50) write(stdo,351) qc,lfoc,qcorg,qcorh
C  351   format(' qc=',f12.6,'   lfoc',i2,'   qcorg,qcorh',2f12.6)
        if (cofh .ne. 0) then
          if (lfoc .eq. 1) then
            lfoc1 = 1
            do  i = 1, ng
              v(1) = gv(i,1)
              v(2) = gv(i,2)
              v(3) = gv(i,3)
              call hklft ( v , rfoc , ceh , tau , alat , kmax , 1 , k0 , rv_a_ocy , gkl )
              cgh1(i) = cgh1(i) + cofh*gkl(0,1)/vol
            enddo
          else if (lfoc .eq. 2) then
            call rx('smcorm: we do not allow now lfoc=2 anymore takao') !Aug2010
            lfoc2 = 1
            do  i = 1, ng
              v(1) = gv(i,1)
              v(2) = gv(i,2)
              v(3) = gv(i,3)
              call hklft ( v , rfoc , ceh , tau , alat , kmax , 1 , k0 , rv_a_ocy , gkl )
              cgh2(i) = cgh2(i) + cofh*gkl(0,1)/vol
            enddo
          endif
        endif
      enddo
      end subroutine smcorm






