subroutine symdmu(nlibu,dmatu,nbas,nsp,lmaxu,sspec,ssite,ng,g, &
  istab,lldau,rms)
  use m_struc_def
  !- Symmetrize LDA+U density matrix dmatu
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   dmatu :density matrix for LDA+U
  !i   dmatw :work array of same dimension as dmatu.
  !i         :Returns original dmatu on output.
  !i   nbas  :size of basis
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   lmaxu :dimensioning parameter for U matrix
  !i   sspec :struct for species-specific information; see routine uspec
  !i         Elts read: lmxa idu
  !i   ssite :struct for site-specific information; see routine usite
  !i         Elts read: spec
  !i   ng    :number of group operations.  Program does nothing if ng=0
  !i   g     :point group operations
  !i   istab :table of site permutations for each group op (mksym.f,symtbl.f)
  !i   istab :site istab(i,ig) is transformed into site i by grp op ig
  !i   lldau :lldau(ib)=0 => no U on this site otherwise
  !i          U on site ib with dmat beginning at dmats(*,lldau(ib))
  !l Local variables
  !l   ofjbl :number of dens mat arrays preceding current one for
  !l         :current site
  ! o Inputs/Outputs
  ! o dmatu  :density matrix or vorb symmetrized on output.
  ! o        :the output dmatu is the sum of the original dmatw
  ! o        :and the symmetrized dmatu.  If output dmatu is
  ! o        :to be a symmetrized version of the input,
  ! o        :dmatw MUST be initially zero on input
  !o Outputs
  !o   rms   :rms change in dmatu from symmetrization
  !r Notes
  !r   Routine uses real harmonics
  !u Updates
  !u   30 Jan 06 dmats now symmetrized across different sites.
  !u             Unsuccessful attmempt to include spinor rotations
  !u   09 Nov 05 (wrl) Convert dmat to complex form
  !u   30 Apr 05 Lambrecht first created
  !--------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nbas,lldau(nbas),ng,nsp,lmaxu,istab(nbas,ng),i_copy_size
  double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
  double complex dmatw(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
  type(s_spec)::sspec(*)
  type(s_site)::ssite(*)

  double precision :: g(9,*),rms
  ! ... Local parameters
  integer :: is,igetss,lmxa,idu(4),m1,m2,ilm1,ilm2,ib,l,isp,m3,m4,ig, &
       iblu,nlibu,jb,jblu,ofjbl,lwarn
  double precision :: rmat(16,16),r(-3:3,-3:3),ddot
  double complex sdmat(-3:3,-3:3,2,2)
  ! ... for spinor rotations
  logical :: cmdopt0
  real(8):: ddet33,xx
  !      double precision eula(3),ddet33,xx
  !      double complex u(2,2,ng)
  character (40) :: str
  dmatw=0d0
  rms = 0
  if (cmdopt0('--nosymdm') .OR. ng == 0) return

  ! --- Setup for spinor part ---
  lwarn = 0
  do  ig = 1, ng
     xx = ddet33(g(1,ig))
     !       Extract to rdmat pure rotation part of g; make Euler angles
     if (xx > 0) then
        call dpcopy(g(1,ig),rmat,1,9,1d0)
     else
        call dpcopy(g(1,ig),rmat,1,9,-1d0)
     endif
     !        call dvset(eula,1,3,0d0)
     !        call rm2eua(rmat,eula(1),eula(2),eula(3))
     !        call rotspu(0,1,1,eula,1,u(1,1,ig))

     !       For debugging
     !        call asymop(rmat,eula,' ',str)
     !        print *, ig, str, sngl(ddet33(g(1,ig))*g(9,ig))
     !        print 433, ig, eula, (ddet33(g(1,ig)) .lt. 0)
     !  433   format(' ig=',i4,' eula = ', 3f10.5, '  inv=', L1)
     !        do  m1 = 1, 2
     !          write(6,'(2f10.5,4x,2f10.5)') (dble(u(m1,m2,ig)), m2=1,2)
     !        enddo
     !        do  m1 = 1, 2
     !          write(6,'(2f10.5,4x,2f10.5)') (dimag(u(m1,m2,ig)), m2=1,2)
     !        enddo

     !   .. for now
     if (dabs(xx*g(9,ig)-1) > 1d-6) lwarn = lwarn+1
  enddo

  if (lwarn > 0) call info(10,0,0, &
       ' symdmu  (warning): %i symops rotate z axis',lwarn,0)

  ! --- For each site density-matrix, do ---
  iblu = 0
  do  ib = 1, nbas
     if (lldau(ib) /= 0) then
        is = int(ssite(ib)%spec)


        lmxa=sspec(is)%lmxa
        i_copy_size=size(sspec(is)%idu)
        call icopy(i_copy_size,sspec(is)%idu,1,idu,1)

        ofjbl = -1
        do  l = 0, min(lmxa,3)
           if (idu(l+1) /= 0) then
              iblu = iblu+1
              ofjbl = ofjbl+1

              !             call zprm('dm spin 2',2,dmatu(-l,-l,2,iblu),7,2*l+1,2*l+1)

              !         --- Loop over group operations ---
              do  ig = 1, ng

                 jb = istab(ib,ig)
                 jblu = lldau(jb) + ofjbl

                 !               Rotation matrices for spherical harmonics up to f orbitals
                 call ylmrtg(16,g(1,ig),rmat)
                 !               Pick out the one we need
                 ilm1 = l**2
                 do  m1 = -l, l
                    ilm1 = ilm1+1
                    ilm2 = l**2
                    do  m2 = -l, l
                       ilm2 = ilm2+1
                       r(m1,m2) = rmat(ilm1,ilm2)
                    enddo
                 enddo

                 !               call prmx('rot',r(-l,-l),7,2*l+1,2*l+1)

                 !           ... Spatial rotation: dmatu(iblu) -> sdmat
                 do  isp = 1, nsp
                    do  m1 = -l, l
                       do  m2 = -l, l
                          sdmat(m1,m2,isp,isp) = 0
                          sdmat(m1,m2,isp,3-isp) = 0
                          do  m3 = -l, l
                             do  m4 = -l, l
                                sdmat(m1,m2,isp,isp) = sdmat(m1,m2,isp,isp) &
                                     + r(m1,m3)*dmatu(m3,m4,isp,iblu)*r(m2,m4)
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo

                 !               call zprm('rot dm(sp2)',2,sdmat(-l,-l,2),7,2*l+1,2*l+1)


                 !           ... Spinor rotation ... give up for now
                 !               Debugging
                 !                print 432, ig, eula, (ddet33(g(1,ig)) .lt. 0)
                 !  432           format(' ig=',i4,' eula = ', 3f10.5, '  inv=', L1,' u:')
                 !                do  m1 = 1, 2
                 !                  write(6,'(2f10.5)') (dble(u(m1,m2,ig)), m2=1,2)
                 !                enddo
                 !                do  m1 = 1, 2
                 !                  write(6,'(2f10.5)') (dimag(u(m1,m2,ig)), m2=1,2)
                 !                enddo

                 !                do  isp = 1, 2
                 !                print *, 'init, spin',isp,' ig=',ig, ' iblu=',iblu
                 !                do  m1 = -l, l, 2*l
                 !                  write(6,'(2f10.5)')
                 !     .              (dble(dmatu(m1,m2,isp,iblu)), m2=-l, l, 2*l)
                 !                enddo
                 !                do   m1 = -l, l, 2*l
                 !                  write(6,'(2f10.5)')
                 !     .              (dimag(dmatu(m1,m2,isp,iblu)), m2=-l, l, 2*l)
                 !                enddo
                 !                enddo
                 !                do  isp = 1, 2
                 !                print *, 'spatial rot spin',isp
                 !                do  m1 = -l, l, 2*l
                 !                  write(6,'(2f10.5)')
                 !     .              (dble(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l)
                 !                enddo
                 !                do   m1 = -l, l, 2*l
                 !                  write(6,'(2f10.5)')
                 !     .              (dimag(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l)
                 !                enddo
                 !                enddo

                 !                do  m1 = -l, l
                 !                do  m2 = -l, l
                 !                  if (m1 .eq. l .and. m2 .eq. l) then
                 !                    print *, 'hi'
                 !                  endif
                 !                  s11 = sdmat(m1,m2,1,1)
                 !                  s12 = sdmat(m1,m2,1,2)
                 !                  s21 = sdmat(m1,m2,2,1)
                 !                  s22 = sdmat(m1,m2,2,2)
                 !                  su(1,1) = s11*u(1,1,ig) + s12*u(2,1,ig)
                 !                  su(2,1) = s21*u(1,1,ig) + s22*u(2,1,ig)
                 !                  su(1,2) = s11*u(1,2,ig) + s12*u(2,2,ig)
                 !                  su(2,2) = s21*u(1,2,ig) + s22*u(2,2,ig)
                 !                  do  i = 1, 2
                 !                  do  j = 1, 2
                 !                    sdmat(m1,m2,i,j) = dconjg(u(1,i,ig))*su(1,j) +
                 !     .                                 dconjg(u(2,i,ig))*su(2,j)
                 !C                    usu(i,j) = dconjg(u(1,i,ig))*su(1,j) +
                 !C     .                         dconjg(u(2,i,ig))*su(2,j)
                 !C                    sdmat(m1,m2,i,j) = su(i,j)
                 !                  enddo
                 !                  enddo
                 !                enddo
                 !                enddo

                 !                do  isp = 1, 2
                 !                  print *, 'spinor rot, isp=',isp,' ig=',ig,'jblu=',jblu
                 !                do  m1 = -l, l, 2*l
                 !                  write(6,'(2f10.5:2x,2f10.5)')
                 !     .              ((dble(sdmat(m1,m2,isp,is)), m2=-l,l,2*l), is=1,2)
                 !C     .              (dble(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l),
                 !                enddo
                 !                do   m1 = -l, l, 2*l
                 !                  write(6,'(2f10.5:2x,2f10.5)')
                 !     .              ((dimag(sdmat(m1,m2,isp,is)), m2=-l,l,2*l), is=1,2)
                 !C    .              (dimag(sdmat(m1,m2,isp,isp)), m2=-l, l, 2*l)
                 !                enddo
                 !                enddo
                 !               pause

                 !           ... Add sdmat/ng into dmat
                 do   isp = 1, nsp
                    do  m1 = -l, l
                       do  m2 = -l, l
                          dmatw(m1,m2,isp,jblu) = dmatw(m1,m2,isp,jblu) + &
                               sdmat(m1,m2,isp,isp)/ng
                       enddo
                    enddo
                 enddo

              enddo

           endif
        enddo
     endif
  enddo

  !     Exchange original for symmetrized dmatu
  nlibu = iblu
  is = nsp*nlibu*(lmaxu*2+1)**2
  call dswap(2*is,dmatw,1,dmatu,1)

  !     RMS change in dmatu
  call daxpy(2*is,-1d0,dmatu,1,dmatw,1)
  rms = dsqrt(ddot(2*is,dmatw,1,dmatw,1)/(2*is))
  call daxpy(2*is,1d0,dmatu,1,dmatw,1)
end subroutine symdmu


double precision function ddet33(matrix)
  !- Calculates the determinant of a 3X3 matrix
  double precision :: matrix(9),ddot,m1cm2(3)
  call cross(matrix(4),matrix(7),m1cm2)
  ddet33 = ddot(3,matrix(1),1,m1cm2,1)
end function ddet33
