!> ESM by Masao Obata, Kanazawa univ. Make ESM electrostatic potential of density given in real space
module m_esmsmves 
  public esmsmves
  private
contains
  subroutine esmsmves(qmom,ng,gv,kv,cv,cg1,cgsum,smrho,qbg,smpot,f,gpot0,hpot0,qsmc,zsum,vrmt)
    use m_lmfinit,only: nsp, alat=>lat_alat,nbas,nlmxlx
    use m_MPItk,only:  master_mpi,mlog
    use m_lattic,only: vol=>lat_vol,plat=>lat_plat
    use m_lgunit,only:stdo
    use m_supot,only:n1,n2,n3
    use m_vesgcm,only: vesgcm
    implicit none
    include "mpif.h"
    integer, intent(in) :: ng, kv(ng,3)
    real(8), intent(in) :: gv(ng,3), qmom(nlmxlx,nbas), qbg
    real(8), intent(out) :: f(3,nbas), gpot0(*), hpot0(nbas), vrmt(nbas)
    real(8), intent(out) :: zsum, qsmc
    complex(8), intent(in)  :: smrho(n1, n2, n3, 2), cgsum(ng)
    complex(8), intent(out) :: smpot(n1, n2, n3), cv(ng), cg1(ng)
    integer, save :: jesm = 0, jtresm
    real(8), save :: sa0, tresm, z1esm, z2esm, z0esm, vesmp, vesmm, eesmp, eesmm
    logical, save :: esm_init = .true.
    logical ::  twrite = .true.
    integer :: ig, i, j, mpipid, ifiese, ierr, ib, is !, i_copy_size, in3, in2, in1!, nsp
    real(8), parameter ::angstr=1d0/0.5291769, ev=2d0/27.2113834, eps=1d-10
    real(8) :: eh, tau(3), vg(3), ddot, gvr, rmt, fac, gvb, g2
    real(8) :: pi, epi, tpiba, tpiba2
    complex(8), allocatable:: vtmp(:,:,:), cwork(:)
    call tcn('esmsmves')
    allocate(vtmp(n1, n2, n3), cwork(ng))
    pi  = 4d0*datan(1d0)
    epi = 8d0*pi
    tpiba = 2*pi/alat
    tpiba2 = tpiba*tpiba
    if( esm_init ) then
       jtresm=0
       tresm=0d0
       z1esm=0d0
       z2esm=0d0
       vesmp=0d0
       vesmm=0d0
       eesmp=0d0
       eesmm=0d0
       if (master_mpi ) then
          open(newunit=ifiese,file='esm_input.dat',status='old',err=201)
          read(ifiese,*,err=201) jesm
          read(ifiese,*,err=201) jtresm,tresm
          read(ifiese,*,err=201) z1esm,z2esm
          read(ifiese,*,err=201) vesmp,vesmm
          read(ifiese,*,err=201) eesmp,eesmm
          close(ifiese)
201       continue
       endif
       call mpibc1_int(jesm, 1,'esmsmves_jesm')
       if(jesm==0) then
          if(master_mpi) write(stdo,"(a)")'  esmsmves: ESM is not turned on, you need esm_input.dat for ESM mode'
          call tcx('esmsmves')
          return
       endif
       if (master_mpi ) then
          tresm=-tresm
          z0esm =  alat*plat(3,3)*0.5d0
          sa0 = alat*alat*abs(plat(1,1)*plat(2,2)-plat(2,1)*plat(1,2))
          write(stdo,*)
          write(stdo,*) '  esmsmves'
          write(stdo, '(3I5,2I10)') n1, n2, n3, ng, ng
          write(stdo,'(3F10.5)') (alat*plat(1,i),i=1,3)
          write(stdo,'(3F10.5)') (alat*plat(2,i),i=1,3)
          write(stdo,'(3F10.5)') (alat*plat(3,i),i=1,3)
          write(stdo,'(A,I3)')'effective screening medium method jesm=',jesm
          write(stdo,'(A,I3,F10.3)')' jtresm,tresm=',jtresm,-tresm
          write(stdo,'(A,3F10.3,A,F10.2)') ' z0esm,z1esm,z2esm=',z0esm,z1esm,z2esm,' sa0=',sa0
          if( jesm == 1 ) then
             write(stdo,*)'system : vaccum(-z)/slab/vaccum(+z) '
             write(stdo,'(A,2F10.5,A,F10.5)')' vesmp,vesmm=',vesmp,vesmm,' vave=',(vesmp+vesmm)*0.5d0
          elseif( jesm == 2 ) then
             write(stdo,*)'system : metal(-z)/slab/metal(+z) '
             write(stdo,'(A,2F10.5,A,F10.5)')' vesmp,vesmm=',vesmp,vesmm,' qbg=',qbg
          elseif( jesm == 3 ) then
             write(stdo,*)'system : vaccum(-z)/slab/metal(+z) '
             write(stdo,'(A,2F10.5)')' qbg=',qbg
          elseif( jesm == 4 ) then
             write(stdo,*)'system : metal(-z)/slab/vaccum(+z) '
             write(stdo,'(A,2F10.5)')' qbg=',qbg
          elseif( jesm == 5 ) then
             write(stdo,*)'system : vacuum(-z)/slab/vacuum(+z) '
             write(stdo,'(A,2F10.5)')' eesmp,eesmm=',eesmp,eesmm
          elseif( jesm == 6 ) then
             write(stdo,*)'system : metal(-z)/slab/metal(+z) '
             write(stdo,'(A,2F10.5)')' vesmp,eesmm=',vesmp,eesmm
          elseif( jesm == 7 ) then
             write(stdo,*)'system : metal(-z)/slab/metal(+z) '
             write(stdo,'(A,2F10.5)')' vesmm,eesmp=',vesmm,eesmp
          elseif( jesm == 10 .OR. jesm == 11 ) then
             write(stdo,*)'same as jesm=0, but employ esm scheme '
             write(stdo,*)'system : vaccum(-z)/slab/vaccum(+z) '
          endif
          write(stdo,*)
       endif
       call mpi_barrier(mpi_comm_world,ierr)
       !        call mpibc1_int (jesm,  1         ,'esmsmves_jesm')
       call mpibc1_int (jtresm,1,'esmsmves_jtresm')
       call mpibc1_real(tresm, 1,'esmsmves_tresm')
       call mpibc1_real(z1esm, 1,'esmsmves_z1esm')
       call mpibc1_real(z2esm, 1,'esmsmves_z2esm')
       call mpibc1_real(vesmp, 1,'esmsmves_vesmp')
       call mpibc1_real(vesmm, 1,'esmsmves_vesmm')
       call mpibc1_real(eesmp, 1,'esmsmves_eesmp')
       call mpibc1_real(eesmm, 1,'esmsmves_eesmm')
       call mpibc1_real(z0esm, 1,'esmsmves_z0esm')
       call mpibc1_real(sa0,   1,'esmsmves_sa0')
       if( jesm == 0 ) then
          if ( master_mpi ) write(stdo,*)'jesm=0 return'
          esm_init = .false.
          call tcx('esmsmves')
          return
       endif
       esm_init = .false.
    endif
    if( jesm == 0 ) then
       if ( master_mpi ) write(stdo,*)'jesm=0 return'
       call tcx('esmsmves')
       return
    endif
    call esmhartp(cgsum, smrho, gv, kv, ng, eh, n1, n2, n3, cv, cg1, cwork, smpot, 'charge_pot.dat')
    cv(1) = 0d0
    do  ig = 2, ng
       g2 = tpiba*tpiba*(gv(ig,1)**2+gv(ig,2)**2+gv(ig,3)**2)
       cv(ig) = epi*cgsum(ig)/g2
    enddo
    call gvputf(ng, 1, kv, n1, n2, n3, cv, vtmp)
    smpot(1:n1,1:n2,1:n3) = smpot(1:n1,1:n2,1:n3) - vtmp(1:n1,1:n2,1:n3)
    call vesgcm(qmom,ng,gv,kv,cv,cg1,cwork,smpot,f,gpot0,hpot0,qsmc,zsum,vrmt)
    if ( master_mpi ) write(stdo,*) "eh=",eh
    call tcx('esmsmves')
    return

  contains
    !-----------------------------------------------------------------------
    subroutine esmhartp(cgsum, smrho, gv, kv, ng, eh,  n1, n2, n3, cv, cg1, cwork, smpot, fname)
      !      INSTALLED DEVELOVED BY Otani&Sugino (2006); Phys.Rev.B73,115407(2006).
      !      OPTION JESM=1,2,3      ORIGINAL OPTIONS (i),(ii),(iii) IN THE PAPER
      use m_lgunit,only:stdo
      implicit none
      integer, intent(in) :: ng, n1, n2, n3
      integer, intent(in) :: kv(ng,3)
      real(8), intent(out) :: eh
      real(8), intent(in) :: gv(ng,3)
      complex(kind(0d0)), intent(in) :: smrho(n1, n2, n3, 2), cgsum(ng)
      complex(kind(0d0)), intent(inout) :: cv(ng), cg1(ng), cwork(ng)
      complex(kind(0d0)), intent(out) :: smpot(n1, n2, n3)
      character(len=*) ,intent(in) :: fname

      complex(kind(0d0)), parameter :: ci=(0d0,1d0)
      complex(kind(0d0)), allocatable :: v(:,:,:), v3(:), xidz(:,:,:)
      complex(kind(0d0)) :: temp
      real(8), allocatable :: xdr(:), ydr(:), xdl(:), ydl(:), g(:), gpas(:,:), v3work(:,:)
      real(8) :: temp1, temp2, temp3, temp4, temp5, clat, delc, z1esm1, z2esm1, z, tg, als, bls, ars, brs
      integer, parameter :: ndat = 5

      integer :: ir, ii, ir1, ir2, ir3, ig1, ig2,ifi

      allocate (v(n1,n2,n3), g(ng), v3(n3), gpas(n1,n2), xidz(n1,n2,4))
      xidz(1:n1,1:n2,1:4) = (0.d0, 0.d0)

      if( twrite ) then
         allocate (xdl(ndat), ydl(ndat), xdr(ndat), ydr(ndat))
         allocate (v3work(n3,5))
      endif

      call gvgetf(ng, 1, kv, n1, n2, n3, smrho(1,1,1,1), cv)

      if ( master_mpi ) then
         write(stdo,'(A,F15.8)') "smrho val =",dble(cv(1))*vol
         write(stdo,'(A,F15.8)') "comepnsating gaussinal val =",dble(cgsum(1))*vol
      endif

      cv(1:ng) = cgsum(1:ng) + cv(1:ng)
      clat=z0esm*2d0
      delc=clat/dble(n3)

      !....TRANSLATION IN Z-COORDINATE
      if( jtresm == 1 ) then
         do ig=2,ng
            cv(ig) = cv(ig)*exp(-ci*tpiba*tresm*delc*gv(ig,3))
         enddo
      else
         do ig=2,ng
            cv(ig) = cv(ig)*exp(-ci*tpiba*tresm*gv(ig,3))
         enddo
      endif

      !....BARE COULOMB PART
      v(1:n1,1:n2,1:n3) = (0.d0,0.d0)
      g(1:ng) = 0d0
      gpas(1:n1,1:n2) = 0d0
      do ig = 2, ng
         g(ig) = gv(ig,1)**2+gv(ig,2)**2+gv(ig,3)**2
         gpas(kv(ig,1),kv(ig,2)) = sqrt(gv(ig,1)**2+gv(ig,2)**2)
         v(kv(ig,1), kv(ig,2), kv(ig,3)) = epi*cv(ig)/(tpiba2*g(ig))
      enddo

      do ir2=1, n2
         do ir1=1, n1
            v3(1:n3) = v(ir1,ir2,1:n3)
            call fftz3(v3, n3, 1, 1, n3, 1, 1, 1, 0, 1)
            smpot(ir1,ir2,1:n3) = v3(1:n3)
         enddo
      enddo

      !-----------------------------------------------------------------------
      ! JESM=1 vacuum/slab/vacuum
      ! Eeffective boundary (+Z) : Z1ESM1=Z0ESM
      !                     (-Z) : Z2ESM1=-Z0ESM
      !-----------------------------------------------------------------------
      if( jesm == 1 ) then

         z1esm1 = z0esm
         z2esm1 = -z0esm

         do ig=2, ng
            temp = epi*exp(ci*tpiba*z1esm1*gv(ig,3))/(tpiba2*g(ig))
            xidz(kv(ig,1), kv(ig,2), 1) = xidz(kv(ig,1), kv(ig,2), 1)+temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 3) = xidz(kv(ig,1), kv(ig,2), 3)+ci*tpiba*gv(ig,3)*temp*cv(ig)
            temp = epi*exp(ci*tpiba*z2esm1*gv(ig,3))/(tpiba2*g(ig))
            xidz(kv(ig,1), kv(ig,2), 2) = xidz(kv(ig,1), kv(ig,2), 2)+temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 4) = xidz(kv(ig,1), kv(ig,2), 4)+ci*tpiba*gv(ig,3)*temp*cv(ig)
         enddo

         !.......G_para = 0
         do ir3=1, n3
            if(ir3 <= n3/2) then
               z=dble(ir3-1)*delc
            else
               z=dble(ir3-n3-1)*delc
            endif

            smpot(1,1,ir3) = smpot(1,1,ir3) &
                 -0.5d0*(xidz(1,1,1)+xidz(1,1,2)) &
                 -0.5d0*(z-z1esm1)*xidz(1,1,3) &
                 -0.5d0*(z-z2esm1)*xidz(1,1,4) &
                 +2d0*pi*(vesmp+vesmm)
         enddo

      elseif( jesm == 2 ) then

         do ig=2, ng
            temp = epi*exp(ci*tpiba*z0esm*gv(ig,3))/(tpiba2*g(ig))
            xidz(kv(ig,1), kv(ig,2), 1) = xidz(kv(ig,1), kv(ig,2), 1)+temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 2) = xidz(kv(ig,1), kv(ig,2), 2)+ci*tpiba*gv(ig,3)*temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 3) = xidz(kv(ig,1), kv(ig,2), 3)+conjg(temp)*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 4) = xidz(kv(ig,1), kv(ig,2), 4)+ci*tpiba*gv(ig,3)*conjg(temp)*cv(ig)
         enddo

         do ir3=1, n3
            if(ir3 <= n3/2) then
               z=dble(ir3-1)*delc
            else
               z=dble(ir3-n3-1)*delc
            endif

            smpot(1,1,ir3) = smpot(1,1,ir3) &
                 +epi*0.5d0*cv(1)*((2d0*z1esm-z0esm)*z0esm-z**2) &
                 +0.5d0*(-(z/z1esm+1d0)*xidz(1,1,1) &
                 +(z/z1esm-1d0)*xidz(1,1,3) &
                 +(-z1esm+z0esm+z*(z0esm/z1esm-1d0))*xidz(1,1,2) &
                 +( z1esm-z0esm+z*(z0esm/z1esm-1d0))*xidz(1,1,4)) &
                 +z/z1esm*(vesmp-vesmm)*0.5D0+(vesmp+vesmm)*0.5D0

            do ig2=1, n2
               do ig1=1, n1
                  if( gpas(ig1,ig2) <= 1d-5 ) cycle
                  tg=tpiba*gpas(ig1,ig2)
                  temp1 = 1d0/(exp((-z+z1esm)*tg)-exp((-z-3d0*z1esm)*tg)) &
                       -1d0/(exp((z+3d0*z1esm)*tg)-exp((z-z1esm)*tg))
                  temp2 = 1d0/(exp((-z+3d0*z1esm)*tg)-exp((-z-z1esm)*tg)) &
                       -1d0/(exp((z+z1esm)*tg)-exp((z-3d0*z1esm)*tg))
                  temp3=exp( (z1esm-z0esm)*tg)
                  temp4=exp(-(z1esm-z0esm)*tg)
                  smpot(ig1,ig2,ir3) = smpot(ig1,ig2,ir3) &
                       -0.5d0*temp1*(xidz(ig1,ig2,1)*(temp3+temp4)+xidz(ig1,ig2,2)*(temp3-temp4)/tg) &
                       +0.5d0*temp2*(xidz(ig1,ig2,3)*(temp3+temp4)+xidz(ig1,ig2,4)*(-temp3+temp4)/tg)
               enddo
            enddo
         enddo

      elseif( jesm == 3 ) then

         do ig=2, ng
            temp = epi*exp(ci*tpiba*z0esm*gv(ig,3))/(tpiba2*g(ig))
            xidz(kv(ig,1), kv(ig,2), 1) = xidz(kv(ig,1), kv(ig,2), 1)+temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 2) = xidz(kv(ig,1), kv(ig,2), 2)+ci*tpiba*gv(ig,3)*temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 3) = xidz(kv(ig,1), kv(ig,2), 3)+conjg(temp)*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 4) = xidz(kv(ig,1), kv(ig,2), 4)+ci*tpiba*gv(ig,3)*conjg(temp)*cv(ig)
         enddo

         do ir3=1, n3
            if( ir3 <= n3/2 ) then
               z = dble(ir3-1)*delc
            else
               z = dble(ir3-n3-1)*delc
            endif

            smpot(1,1,ir3) = smpot(1,1,ir3) &
                 +epi*0.5d0*cv(1)*(4d0*z1esm*z0esm-(z+z0esm)**2) &
                 -xidz(1,1,1)+(z0esm-z1esm)*xidz(1,1,2)+(z1esm-z)*xidz(1,1,4)

            do ig2=1, n2
               do ig1=1, n1
                  if( gpas(ig1,ig2) <= 1d-5 ) cycle
                  tg = tpiba*gpas(ig1,ig2)
                  temp1 = exp((z-z0esm)*tg)
                  temp2 = exp(-(z+z0esm)*tg)
                  temp3 = exp((z-2d0*z1esm+z0esm)*tg)
                  temp4 = exp((z-2d0*z1esm-z0esm)*tg)
                  smpot(ig1,ig2,ir3) = smpot(ig1,ig2,ir3) &
                       -0.5d0*(temp1+temp3)*xidz(ig1,ig2,1) &
                       -0.5d0*(temp2-temp4)*xidz(ig1,ig2,3) &
                       -0.5d0*(temp1-temp3)*xidz(ig1,ig2,2)/tg &
                       +0.5d0*(temp2-temp4)*xidz(ig1,ig2,4)/tg
               enddo
            enddo
         enddo

      elseif( jesm == 4 )then

         do ig=2,ng
            temp=epi*exp(-ci*tpiba*z2esm*gv(ig,3))/(tpiba2*g(ig))
            xidz(kv(ig,1), kv(ig,2), 2) = xidz(kv(ig,1), kv(ig,2), 2)+conjg(temp)*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 3) = xidz(kv(ig,1), kv(ig,2), 3)+temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 4) = xidz(kv(ig,1), kv(ig,2), 4)+conjg(ci*tpiba*gv(ig,3)*temp*cv(ig))
         enddo

         do ir3=1,n3
            if( ir3 <= n3/2 ) then
               z = dble(ir3-1)*delc
            else
               z = dble(ir3-n3-1)*delc
            endif

            smpot(1,1,ir3) = smpot(1,1,ir3) &
                 +epi*0.5d0*cv(1)*(-z2esm+z)*(-z2esm-z+2d0*z1esm) &
                 -xidz(1,1,2)-(-z2esm+z)*xidz(1,1,4)
            do ig2=1, n2
               do ig1=1, n1
                  if( gpas(ig1,ig2) <= 1d-5 ) cycle
                  temp1 = exp(-(z-z2esm)*tpiba*gpas(ig1,ig2))
                  temp2 = exp( (z-z1esm)*tpiba*gpas(ig1,ig2))
                  temp3 = exp(-(z-2d0*z2esm+z1esm)*tpiba*gpas(ig1,ig2))
                  smpot(ig1,ig2,ir3) = smpot(ig1,ig2,ir3) &
                       -xidz(ig1,ig2,2)*temp1 &
                       -0.5d0*(temp2-temp3)*(xidz(ig1,ig2,3)+ &
                       +xidz(ig1,ig2,4)/(tpiba*gpas(ig1,ig2)))
               enddo
            enddo
         enddo

      elseif( jesm == 5 )then

         do ig=2,ng
            temp=epi*exp(ci*tpiba*z1esm*gv(ig,3))/(tpiba2*g(ig))
            xidz(kv(ig,1), kv(ig,2), 3) = xidz(kv(ig,1), kv(ig,2), 3) + ci*tpiba*gv(ig,3)*temp*cv(ig)
            xidz(kv(ig,1), kv(ig,2), 4) = xidz(kv(ig,1), kv(ig,2), 4) + ci*tpiba*gv(ig,3)*conjg(temp)*cv(ig)
         enddo

         do ir3=1,n3
            if( ir3 <= n3/2 ) then
               z = dble(ir3-1)*delc
            else
               z = dble(ir3-n3-1)*delc
            endif

            smpot(1,1,ir3) = smpot(1,1,ir3) &
                 -0.25d0*((z+z1esm)**2/z1esm)*(xidz(1,1,3)-eesmp) &
                 +0.25d0*((z-z1esm)**2/z1esm)*(xidz(1,1,4)-eesmm)

            do ig2=1, n2
               do ig1=1, n1
                  if( gpas(ig1,ig2) <= 1d-5 ) cycle
                  temp1 = cosh((z+z1esm)*tpiba*gpas(ig1,ig2))/(sinh(2d0*z1esm*tpiba*gpas(ig1,ig2))*tpiba*gpas(ig1,ig2))
                  temp2 = cosh((z-z1esm)*tpiba*gpas(ig1,ig2))/(sinh(2d0*z1esm*tpiba*gpas(ig1,ig2))*tpiba*gpas(ig1,ig2))
                  smpot(ig1,ig2,ir3) = smpot(ig1,ig2,ir3) &
                       -(xidz(ig1,ig2,3)-eesmp)*temp1 &
                       +(xidz(ig1,ig2,4)-eesmm)*temp2
               enddo
            enddo
         enddo

      endif
      !.....endif jesm

      call gvputf(ng, 1, kv, n1, n2, n3, cv, v)

      do ir2=1, n2
         do ir1=1, n1
            v3(1:n3) = v(ir1,ir2,1:n3)
            call fftz3(v3, n3, 1, 1, n3, 1, 1, 1, 0, 1)
            v(ir1,ir2,1:n3) = v3(1:n3)
         enddo
      enddo

      !     NOW V:
      !           n1*n2 (G-SPACE)
      !           n3      (R-SPACE)

      temp = (0d0,0d0)
      do ir3=1, n3
         do ig2=1, n2
            do ig1=1, n1
               temp = temp+conjg(v(ig1,ig2,ir3))*smpot(ig1,ig2,ir3)
            enddo
         enddo
      enddo

      eh = dble(temp)*0.5d0*sa0*delc

      do ir3=1, n3
         call fftz3(smpot(1,1,ir3), n1, n2, 1, n1, n2, 1, 1, 0, 1)
      enddo

      !     NOW smpot:
      !           n1*n2 (R-SPACE)
      !           n3      (R-SPACE)


      !.....WRITE CHARGE DISTRIBUTION AND HARTREE POTENTIAL
      if( twrite ) then
         do ir3=1, n3
            call fftz3(v(1,1,ir3), n1, n2, 1, n1, n2, 1, 1, 0, 1)
         enddo

         !.......smooth density and esm potential
         do ir3=1,n3
            temp1 = 0d0
            temp2 = 0d0
            do ir2=1,n2
               do ir1=1,n1
                  temp1 = temp1+dble(v(ir1,ir2,ir3))
                  temp2 = temp2+dble(smpot(ir1,ir2,ir3))
               enddo
            enddo
            v3work(ir3,1) = temp1/dble(n1*n2)
            v3work(ir3,2) = temp2/dble(n1*n2)
         enddo

         !.......spin density and gaussian density
         if ( nsp == 2 ) then
            call gvgetf(ng, 1, kv, n1, n2, n3, smrho(1,1,1,1), cv)

            v(1:n1,1:n2,1:n3) = smrho(1:n1,1:n2,1:n3,2)
            call fftz3(v, n1, n2, n3, n1, n2, n3, 1, 0, -1)
            call gvgetf(ng, 1, kv, n1, n2, n3, v, cwork)
            cwork(1:ng) = cv(1:ng) - 2d0*cwork(1:ng)
         endif

         cg1(1:ng) = cgsum(1:ng)

         if( jtresm == 1 ) then
            cg1(1:ng) = cg1(1:ng)*exp(-ci*tpiba*tresm*delc*gv(1:ng,3))
            cwork(1:ng) = cwork(1:ng)*exp(-ci*tpiba*tresm*delc*gv(1:ng,3))
         else
            cg1(1:ng) = cg1(1:ng)*exp(-ci*tpiba*tresm*gv(1:ng,3))
            cwork(1:ng) = cwork(1:ng)*exp(-ci*tpiba*tresm*gv(1:ng,3))
         endif

         call gvputf(ng, 1, kv, n1, n2, n3, cg1, v)
         call fftz3(v, n1, n2, n3, n1, n2, n3, 1, 0, 1)

         do ir3=1,n3
            temp3 = 0d0
            do ir2=1,n2
               do ir1=1,n1
                  temp3 = temp3+dble(v(ir1,ir2,ir3))
               enddo
            enddo
            v3work(ir3,3) = temp3/dble(n1*n2)
         enddo

         call gvputf(ng, 1, kv, n1, n2, n3, cwork, v)
         call fftz3(v, n1, n2, n3, n1, n2, n3, 1, 0, 1)

         do ir3=1,n3
            temp4 = 0d0
            do ir2=1,n2
               do ir1=1,n1
                  temp4 = temp4+dble(v(ir1,ir2,ir3))
               enddo
            enddo
            v3work(ir3,4) = temp4/dble(n1*n2)
         enddo

         !.......splined esm potential
         v(1:n1,1:n2,1:n3) = smpot(1:n1,1:n2,1:n3)
         call cubic(n1, n2, n3, v)
         do ir3=1,n3
            temp5 = 0d0
            do ir2=1,n2
               do ir1=1,n1
                  temp5 = temp5+dble(v(ir1,ir2,ir3))
               enddo
            enddo
            v3work(ir3,5) = temp5/dble(n1*n2)
         enddo

         if ( master_mpi ) then
            open(newunit=ifi,file=trim(fname))
            write(ifi,'(a)') '# charge and hartree potential (planar average)'
            write(ifi,'(a,f18.8)')'# eh=',eh
            write(ifi,'(a)')'#z-coordinate, smooth density, esm potential, gaussian density, spin density, splined esm potential'

            do ir3=n3/2+1,n3
               temp1 = v3work(ir3,1)
               temp2 = v3work(ir3,2)
               temp3 = v3work(ir3,3)
               temp4 = v3work(ir3,4)
               temp5 = v3work(ir3,5)
               write(ifi,'(E16.7,5E17.8)') dble(ir3-n3-1)*delc,temp1,temp2,temp3,temp4,temp5
               if( jesm == 2 )then
                  if( ir3-n3/2-1 < 5 ) then
                     ii = ir3-n3/2
                     xdl(ii) = dble(ir3-1)*delc
                     ydl(ii) = temp2
                  endif
               endif
            enddo

            do ir3=1,n3/2
               temp1 = v3work(ir3,1)
               temp2 = v3work(ir3,2)
               temp3 = v3work(ir3,3)
               temp4 = v3work(ir3,4)
               temp5 = v3work(ir3,5)
               write(ifi,'(E16.7,5E17.8)') dble(ir3-1)*delc,temp1,temp2,temp3,temp4,temp5
               if( jesm == 2 .OR. jesm == 3 ) then
                  if( n3/2-ir3 < 5) then
                     ii = n3/2-ir3+1
                     xdr(ii) = dble(ir3-1)*delc
                     ydr(ii) = temp2
                  endif
               endif
            enddo

            if( jesm == 2 ) then
               !           least-square fit at left surface
               call lsfit(ndat,xdl,ydl,als,bls)
               write(ifi,'(A)')'# left surface: e = a * z + b'
               write(ifi,'(A,2F18.8)')'# als(au,ev/a)=',als,als*angstr/ev
               write(ifi,'(A,2F18.8)')'# bls(au,ev)  =',bls,bls/ev
            endif

            if( jesm == 2 .OR. jesm == 3 ) then
               !           least-square fit at right surface
               call lsfit(ndat,xdr,ydr,ars,brs)
               write(ifi,'(A)')'# right surface: e = a * z + b'
               write(ifi,'(A,2F18.8)')'# ars(au,ev/a)=',ars,ars*angstr/ev
               write(ifi,'(A,2F18.8)')'# brs(au,ev)  =',brs,brs/ev
            endif
            close(ifi)
         endif
      endif !endif of if(twrite)then

      !.....INTERPOLATE THE POTENTIAL BY SPLINE METHOD
      if (jesm /= 11) then
         call cubic(n1, n2, n3, smpot)
      endif

      !.....BACK TO G-SPACE
      call fftz3(smpot, n1, n2, n3, n1, n2, n3, 1, 0, -1)

      call gvgetf(ng, 1, kv, n1, n2, n3, smpot, cv)
      !.....BACK TRANSLATION IN Z-COORDINATE
      if( jtresm == 1 ) then
         do ig=2,ng
            cv(ig) = cv(ig)*exp(ci*tpiba*tresm*delc*gv(ig,3))
         enddo
      else
         do ig=2,ng
            cv(ig) = cv(ig)*exp(ci*tpiba*tresm*gv(ig,3))
         enddo
      endif

      call gvputf(ng, 1, kv, n1, n2, n3, cv, smpot)

      deallocate (v,xidz)

      if ( twrite ) then
         deallocate (xdl,ydl,xdr,ydr)
      endif

      return
    end subroutine esmhartp
    !-----------------------------------------------------------------------
    subroutine lsfit(ndat,xd,yd,a,b)

      integer, intent(in) :: ndat
      real(8), intent(in) :: xd(ndat),yd(ndat)
      real(8), intent(out) :: a,b
      real(8) :: rm(1:2,1:2),vs(1:2),det

      !     least-square fit with f=a*x+b (by m.tsujikawa)

      rm(1,1) = 0d0
      rm(1,2) = 0d0
      rm(2,1) = 0d0
      rm(2,2) = 0d0
      vs(1) = 0d0
      vs(2) = 0d0

      do i=1,ndat
         rm(1,1) = rm(1,1)+xd(i)**2
         rm(1,2) = rm(1,2)+xd(i)
         vs(1) = vs(1)+yd(i)*xd(i)
         vs(2) = vs(2)+yd(i)
      enddo
      rm(2,1) = rm(1,2)
      rm(2,2) = dble(ndat)

      det = rm(1,1)*rm(2,2)-rm(2,1)*rm(1,2)
      a = (rm(2,2)*vs(1)-rm(1,2)*vs(2))/det
      b = (-rm(2,1)*vs(1)+rm(1,1)*vs(2))/det

      return
    end subroutine lsfit

    !-----------------------------------------------------------------------
    subroutine cubic(n1, n2, n3, v)
      !     COMPUTE : CONNECT THE EDGES OF HARTREE POTENTIAL
      !     V  INPUT  : HARTREE POTENTIAL (R-SPACE)
      !     V  OUTPUT : HARTREE POTENTIAL (R-SPACE)

      implicit none
      integer, intent(in) :: n1,n2,n3
      complex(kind(0d0)), intent(inout) :: v(n1,n2,n3)

      real(8) :: xnrh, slope1, slope2, zl4, zr4, a1, a2, a3, a4
      integer :: ir1, ir2, nrh

      xnrh = n3*0.5d0
      nrh = n3/2

      do ir2=1,n2
         do ir1=1,n1

            zl4 = real(v(ir1,ir2,nrh-3))
            zr4 = real(v(ir1,ir2,nrh+5))

            slope1 = 0.5*(v(ir1,ir2,nrh-2)-v(ir1,ir2,nrh-4))
            slope2 = 0.5*(v(ir1,ir2,nrh+6)-v(ir1,ir2,nrh+4))

            a4 = (zl4-zr4+4d0*(slope1+slope2))/256d0
            a3 = (slope2-slope1)/16d0-3d0*a4*(xnrh+1d0)
            a2 = slope1-2d0*a3*(xnrh-3d0)-3d0*a4*(xnrh-3d0)**2
            a1 = zr4-a2*(xnrh+5d0)-a3*(xnrh+5d0)**2-a4*(xnrh+5d0)**3

            v(ir1,ir2,nrh-2) = a1+a2*(xnrh-2d0)+a3*(xnrh-2d0)**2+a4*(xnrh-2d0)**3
            v(ir1,ir2,nrh-1) = a1+a2*(xnrh-1d0)+a3*(xnrh-1d0)**2+a4*(xnrh-1d0)**3
            v(ir1,ir2,nrh)   = a1+a2*xnrh+a3*xnrh**2+a4*xnrh**3
            v(ir1,ir2,nrh+1) = a1+a2*(xnrh+1d0)+a3*(xnrh+1d0)**2+a4*(xnrh+1d0)**3
            v(ir1,ir2,nrh+2) = a1+a2*(xnrh+2d0)+a3*(xnrh+2d0)**2+a4*(xnrh+2d0)**3
            v(ir1,ir2,nrh+3) = a1+a2*(xnrh+3d0)+a3*(xnrh+3d0)**2+a4*(xnrh+3d0)**3
            v(ir1,ir2,nrh+4) = a1+a2*(xnrh+4d0)+a3*(xnrh+4d0)**2+a4*(xnrh+4d0)**3

         enddo
      enddo

      return
    end subroutine cubic
  end subroutine esmsmves
end module m_esmsmves
