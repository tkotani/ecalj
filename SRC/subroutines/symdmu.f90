subroutine symdmu(nlibu,dmatu,nbas,nsp,lmaxu,ng,g, istab,lldau,rms)!- Symmetrize LDA+U density matrix dmatu
  use m_struc_def
  use m_ftox
  use m_lmfinit,only:idu,ispec,sspec=>v_sspec
  use m_lgunit,only:stdo
  !i Inputs
  !i   dmatu :density matrix for LDA+U
  !i   dmatw :work array of same dimension as dmatu.
  !i         :Returns original dmatu on output.
  !i   nbas  :size of basis
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   lmaxu :dimensioning parameter for U matrix
  !i   sspec :struct for species-specific information; see routine uspec
  !i         Elts read: lmxa idu
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
  integer :: nbas,lldau(nbas),ng,nsp,lmaxu,istab(nbas,ng),i_copy_size
  integer :: is,igetss,lmxa,m1,m2,ilm1,ilm2,ib,l,isp,m3,m4,ig,iblu,nlibu,jb,jblu,ofjbl,lwarn
  real(8):: rmat(16,16),r(-3:3,-3:3),ddot,g(9,*),rms,ddet33,xx
  complex(8):: sdmat(-3:3,-3:3,2,2),&  ! ... for spinor rotations
       dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),&
       dmatw(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
  logical :: cmdopt0
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
     if (dabs(xx*g(9,ig)-1) > 1d-6) lwarn = lwarn+1
  enddo
  if (lwarn > 0) write(stdo,ftox)'symdmu (warning): ',lwarn,'symops rotate z axis'
  ! --- For each site density-matrix, do ---
  iblu = 0
  do  ib = 1, nbas
     if (lldau(ib) /= 0) then
        is = ispec(ib) !ssite(ib)%spec
        lmxa=sspec(is)%lmxa
        ofjbl = -1
        do  l = 0, min(lmxa,3)
           if (idu(l+1,is) /= 0) then
              iblu = iblu+1
              ofjbl = ofjbl+1
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
