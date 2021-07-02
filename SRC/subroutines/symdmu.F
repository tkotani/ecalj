CSFCPP#define F90 1
      subroutine symdmu(dmatu,dmatw,nbas,nsp,lmaxu,sspec,ssite,ng,g,
     .istab,lldau,rms)

      use m_struc_def  !Cgetarg

C- Symmetrize LDA+U density matrix dmatu
C ----------------------------------------------------------------------
Ci Inputs
Ci   dmatu :density matrix for LDA+U
Ci   dmatw :work array of same dimension as dmatu.
Ci         :Returns original dmatu on output.
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmaxu :dimensioning parameter for U matrix
Ci   sspec :struct for species-specific information; see routine uspec
Ci         Elts read: lmxa idu
Ci   ssite :struct for site-specific information; see routine usite
Ci         Elts read: spec
Ci   ng    :number of group operations.  Program does nothing if ng=0
Ci   g     :point group operations
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   istab :site istab(i,ig) is transformed into site i by grp op ig
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci          U on site ib with dmat beginning at dmats(*,lldau(ib))
Cl Local variables
Cl   ofjbl :number of dens mat arrays preceding current one for
Cl         :current site
Cio Inputs/Outputs
Cio dmatu  :density matrix or vorb symmetrized on output.
Cio        :the output dmatu is the sum of the original dmatw
Cio        :and the symmetrized dmatu.  If output dmatu is
Cio        :to be a symmetrized version of the input,
Cio        :dmatw MUST be initially zero on input
Co Outputs
Co   rms   :rms change in dmatu from symmetrization
Cr Notes
Cr   Routine uses real harmonics
Cu Updates
Cu   30 Jan 06 dmats now symmetrized across different sites.
Cu             Unsuccessful attmempt to include spinor rotations
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   30 Apr 05 Lambrecht first created
C--------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,lldau(nbas),ng,nsp,lmaxu,istab(nbas,ng),i_copy_size
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      double complex dmatw(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)

      double precision g(9,*),rms
C ... Local parameters
      integer is,igetss,lmxa,idu(4),m1,m2,ilm1,ilm2,ib,l,isp,m3,m4,ig,
     .iblu,nlibu,jb,jblu,ofjbl,lwarn
      double precision rmat(16,16),r(-3:3,-3:3),ddot
      double complex sdmat(-3:3,-3:3,2,2)
C ... for spinor rotations
      logical cmdopt0
      double precision eula(3),ddet33,xx
CSFCPP#if F90 | AUTO_ARRAY
C      integer i,j
      double complex u(2,2,ng)
C     double complex s11,s12,s21,s22,su(2,2),usu(2,2)
CSFCPP#else
CSFCPP      integer ngmx
CSFCPP      parameter(ngmx=48)
CSFCPP      double complex u(2,2,ngmx)
CSFCPP#endif
      character *40 str

      rms = 0
      if (cmdopt0('--nosymdm') .or. ng .eq. 0) return

C --- Setup for spinor part ---
      lwarn = 0
      do  ig = 1, ng
        xx = ddet33(g(1,ig))
C       Extract to rdmat pure rotation part of g; make Euler angles
        if (xx .gt. 0) then
          call dpcopy(g(1,ig),rmat,1,9,1d0)
        else
          call dpcopy(g(1,ig),rmat,1,9,-1d0)
        endif
        call dvset(eula,1,3,0d0)
        call rm2eua(rmat,eula(1),eula(2),eula(3))
        call rotspu(0,1,1,eula,1,u(1,1,ig))

C       For debugging
C        call asymop(rmat,eula,' ',str)
C        print *, ig, str, sngl(ddet33(g(1,ig))*g(9,ig))
C        print 433, ig, eula, (ddet33(g(1,ig)) .lt. 0)
C  433   format(' ig=',i4,' eula = ', 3f10.5, '  inv=', L1)
C        do  m1 = 1, 2
C          write(6,'(2f10.5,4x,2f10.5)') (dble(u(m1,m2,ig)), m2=1,2)
C        enddo
C        do  m1 = 1, 2
C          write(6,'(2f10.5,4x,2f10.5)') (dimag(u(m1,m2,ig)), m2=1,2)
C        enddo

C   .. for now
        if (dabs(xx*g(9,ig)-1) .gt. 1d-6) lwarn = lwarn+1
      enddo

      if (lwarn .gt. 0) call info(10,0,0,
     .' symdmu  (warning): %i symops rotate z axis',lwarn,0)

C --- For each site density-matrix, do ---
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)


          lmxa=sspec(is)%lmxa
          i_copy_size=size(sspec(is)%idu)
          call icopy(i_copy_size,sspec(is)%idu,1,idu,1)

          ofjbl = -1
          do  l = 0, min(lmxa,3)
            if (idu(l+1) .ne. 0) then
              iblu = iblu+1
              ofjbl = ofjbl+1

C             call zprm('dm spin 2',2,dmatu(-l,-l,2,iblu),7,2*l+1,2*l+1)

C         --- Loop over group operations ---
              do  ig = 1, ng

                jb = istab(ib,ig)
                jblu = lldau(jb) + ofjbl

C               Rotation matrices for spherical harmonics up to f orbitals
                call ylmrtg(16,g(1,ig),rmat)
C               Pick out the one we need
                ilm1 = l**2
                do  m1 = -l, l
                  ilm1 = ilm1+1
                  ilm2 = l**2
                  do  m2 = -l, l
                    ilm2 = ilm2+1
                    r(m1,m2) = rmat(ilm1,ilm2)
                  enddo
                enddo

C               call prmx('rot',r(-l,-l),7,2*l+1,2*l+1)

C           ... Spatial rotation: dmatu(iblu) -> sdmat
                do  isp = 1, nsp
                  do  m1 = -l, l
                    do  m2 = -l, l
                      sdmat(m1,m2,isp,isp) = 0
                      sdmat(m1,m2,isp,3-isp) = 0
                      do  m3 = -l, l
                        do  m4 = -l, l
                          sdmat(m1,m2,isp,isp) = sdmat(m1,m2,isp,isp)
     .                    + r(m1,m3)*dmatu(m3,m4,isp,iblu)*r(m2,m4)
                        enddo
                      enddo
                    enddo
                  enddo
                enddo

C               call zprm('rot dm(sp2)',2,sdmat(-l,-l,2),7,2*l+1,2*l+1)


C           ... Spinor rotation ... give up for now
C               Debugging
C                print 432, ig, eula, (ddet33(g(1,ig)) .lt. 0)
C  432           format(' ig=',i4,' eula = ', 3f10.5, '  inv=', L1,' u:')
C                do  m1 = 1, 2
C                  write(6,'(2f10.5)') (dble(u(m1,m2,ig)), m2=1,2)
C                enddo
C                do  m1 = 1, 2
C                  write(6,'(2f10.5)') (dimag(u(m1,m2,ig)), m2=1,2)
C                enddo
C
C                do  isp = 1, 2
C                print *, 'init, spin',isp,' ig=',ig, ' iblu=',iblu
C                do  m1 = -l, l, 2*l
C                  write(6,'(2f10.5)')
C     .              (dble(dmatu(m1,m2,isp,iblu)), m2=-l, l, 2*l)
C                enddo
C                do   m1 = -l, l, 2*l
C                  write(6,'(2f10.5)')
C     .              (dimag(dmatu(m1,m2,isp,iblu)), m2=-l, l, 2*l)
C                enddo
C                enddo
C                do  isp = 1, 2
C                print *, 'spatial rot spin',isp
C                do  m1 = -l, l, 2*l
C                  write(6,'(2f10.5)')
C     .              (dble(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l)
C                enddo
C                do   m1 = -l, l, 2*l
C                  write(6,'(2f10.5)')
C     .              (dimag(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l)
C                enddo
C                enddo
C
C                do  m1 = -l, l
C                do  m2 = -l, l
C                  if (m1 .eq. l .and. m2 .eq. l) then
C                    print *, 'hi'
C                  endif
C                  s11 = sdmat(m1,m2,1,1)
C                  s12 = sdmat(m1,m2,1,2)
C                  s21 = sdmat(m1,m2,2,1)
C                  s22 = sdmat(m1,m2,2,2)
C                  su(1,1) = s11*u(1,1,ig) + s12*u(2,1,ig)
C                  su(2,1) = s21*u(1,1,ig) + s22*u(2,1,ig)
C                  su(1,2) = s11*u(1,2,ig) + s12*u(2,2,ig)
C                  su(2,2) = s21*u(1,2,ig) + s22*u(2,2,ig)
C                  do  i = 1, 2
C                  do  j = 1, 2
C                    sdmat(m1,m2,i,j) = dconjg(u(1,i,ig))*su(1,j) +
C     .                                 dconjg(u(2,i,ig))*su(2,j)
CC                    usu(i,j) = dconjg(u(1,i,ig))*su(1,j) +
CC     .                         dconjg(u(2,i,ig))*su(2,j)
CC                    sdmat(m1,m2,i,j) = su(i,j)
C                  enddo
C                  enddo
C                enddo
C                enddo
C
C                do  isp = 1, 2
C                  print *, 'spinor rot, isp=',isp,' ig=',ig,'jblu=',jblu
C                do  m1 = -l, l, 2*l
C                  write(6,'(2f10.5:2x,2f10.5)')
C     .              ((dble(sdmat(m1,m2,isp,is)), m2=-l,l,2*l), is=1,2)
CC     .              (dble(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l),
C                enddo
C                do   m1 = -l, l, 2*l
C                  write(6,'(2f10.5:2x,2f10.5)')
C     .              ((dimag(sdmat(m1,m2,isp,is)), m2=-l,l,2*l), is=1,2)
CC    .              (dimag(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l)
C                enddo
C                enddo
C               pause

C           ... Add sdmat/ng into dmat
                do   isp = 1, nsp
                  do  m1 = -l, l
                    do  m2 = -l, l
                      dmatw(m1,m2,isp,jblu) = dmatw(m1,m2,isp,jblu) +
     .                sdmat(m1,m2,isp,isp)/ng
                    enddo
                  enddo
                enddo

              enddo

            endif
          enddo
        endif
      enddo

C     Exchange original for symmetrized dmatu
      nlibu = iblu
      is = nsp*nlibu*(lmaxu*2+1)**2
      call dswap(2*is,dmatw,1,dmatu,1)

C     RMS change in dmatu
      call daxpy(2*is,-1d0,dmatu,1,dmatw,1)
      rms = dsqrt(ddot(2*is,dmatw,1,dmatw,1)/(2*is))
      call daxpy(2*is,1d0,dmatu,1,dmatw,1)

      end subroutine symdmu


