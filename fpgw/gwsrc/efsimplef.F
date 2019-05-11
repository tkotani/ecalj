      subroutine efsimplef2a (ispin,wibz,qibz,ginv,
     i             nband,nqibz,
     i             konfig,z,nl,natom,iclass,nclass, 
     i             valn, legas, esmr, !input only for empty case
     i             qbz,nqbz, !index_qbz, n_index_qbz,
     o             ef)
      use m_readeigen, only: readeval 
      use m_mpi, only: mpi__root 
!!== Calculate efermi for discrete sum. (not for tetrahedron method) ==
!! You need to call init_reaeigen before you call this.
!! user readeval (readeigen.f) to get eigenvalues.
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level
c nband   = no. states
c nqbz    = no. k-points
c valn    = number of valence electron.
c
c -------------------
c      e(iband) < efm : occupation is one
c efm< e(iband) < efp : occupation is wfacef.
c efp< e(iband)       : occupation is zero
c
      implicit none
      integer(4):: is,iq,ispin,nqibz,ik,isig,kpx,nband,ifev(2)
      integer(4):: ieaord(nband*nqibz*ispin),mbytes,mwords,iwksize,
     &             iqibz  !,iindxk
      real(8)   :: ekt(nband, nqibz,ispin), ektx(nband*nqibz*ispin)
      real(8)   :: wgt(nband, nqibz,ispin), wgtx(nband*nqibz*ispin)
      real(8)   :: qbzx(3), qx(3),qbas(3,3),ginv(3,3),wwg
c
      integer(4):: nclass,natom,nl,ncore,l,ia,ic !,indxk(*)
      real(8)   :: wibz(nqibz),valn,ef, z(nclass),qibz(3,nqibz)
      integer(4):: iclass(natom),konfig(0:nl-1,nclass),ierr
c
      integer(4) :: nbnqnsp,ix,ikx,ikini,nne
      real(8)    :: ew1,ew2,ein,valx,enumef_gauss,esmr, efini
     &            ,eee2,wwg2 ,enumef
c      real(8) :: efp,efm,wwgo,wfacef
      logical :: legas,autoew,GaussSmear=.true. !is external

      integer(4):: nqbz,ifile_handle,if8301,if8302
c      integer(4):: n_index_qbz,index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
      real(8)   :: qbz(3,nqbz)
C--------------------------------------------------------------------
      autoew =.false.
      if(GaussSmear) then
        write(6,*)' efsimplef2(gaussian mode):start'
      else
        write(6,*)' efsimplef2(rectangular mode):start'
      endif
      if(esmr<=0d0) autoew=.true.
c total valence charge
      if(legas) then
        write(6,*)' efsimplef2: legas=T use given valn = ',valn
      else
        valn    = 0d0
        do ia   = 1,natom
          ic    = iclass(ia)
          valn  = valn + z(ic)
          write(6,*)' ia z(ic)=',ia, z(ic)
          do    l = 0,nl-1
            write(6,*)' l (konfig(l,ic)-l-1) 2*(2l+1)=',l,(konfig(l,ic)-l-1),( 2*l +1)*2
            valn  = valn - (konfig(l,ic)-l-1) *( 2*l +1)*2
          end do
        end do
      endif
ccccccccccccccccccccccccccccccccccccccccccc
c      valn=29d0
cccccccccccccccccccccccccccccccccccccccccc
        write(6,*)' valn=',valn
c
      do is = 1,ispin
        do iq = 1,nqibz
          call readeval(qibz(:,iq),is, ekt(1,iq,is) )
c        call catch1 (qibz,iq,3,nqibz,qx)
c         call fbz2(qx,  ginv,index_qbz,n_index_qbz,qbz,nqbz,
c     o              qbzx,kpx)
c        call rwdd1( ifev(is),   !read eigenvalues
c     d              kpx,nband,
c     o              ekt(1,iq,is) )
        enddo
      enddo
c
      if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
        write(6,*) 'sum (wibz)=', sum(wibz(1:nqibz))
Cstop2rx 2013.08.09 kino        stop 'efsimplef2: wibzsumerr'
        call rx( 'efsimplef2: wibzsumerr')
      endif
      do is = 1,ispin
        do iq = 1,nqibz
          wgt(1:nband,iq,is) = wibz(iq)
          if(ispin==2) wgt(1:nband,iq,is) = wgt(1:nband,iq,is)/2d0
        enddo
      enddo

c ekt and wgt
      call dcopy ( nband*nqibz*ispin, ekt,1, ektx ,1)
      call dcopy ( nband*nqibz*ispin, wgt,1, wgtx ,1)
      call sortea( ektx,ieaord,nband*nqibz*ispin,isig)
c
c      write(6,*)nband,nqibz,ispin,nband*nqibz*ispin
c      do ik  = 1, nband*nqibz*ispin
c        write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wgtx(ieaord(ik))
c      enddo
c
      if (mpi__root) then
         if8301=ifile_handle()
         if8302=ifile_handle()
         open(if8301,file = "DOSACC.lda")
         open(if8302,file = "DOSACC2.lda")
      endif
      wwg = 0d0
      wwg2= 0d0
      eee2= -1d99
      ikini= 1
      ierr= 1
      nne = nband*nqibz*ispin
      do ik = 1, nne
c---
        if(eee2 +1d-4< ektx(ieaord(ik)).or.ik==nne ) then
! degeneracy check
          if (mpi__root) then 
            if(ik/=1) write(8302,"(2i6,2d23.15)") ikini,ik-1,eee2,wwg2
          endif 
          wwg2 = wgtx(ieaord(ik))
          eee2 = ektx(ieaord(ik))
          ikini =ik
        else
          wwg2= wwg2 + wgtx(ieaord(ik))
        endif
c---
        wwg = wwg + wgtx(ieaord(ik))
        if(wwg<valn+2d0) write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg
c
        if (mpi__root) then 
          write(8301,"(2i6,3d23.15)")
     &    ik,ieaord(ik),ektx(ieaord(ik)),wwg,wgtx(ieaord(ik))
        endif 
c
        if( wwg>valn-1d-8.and.ierr==1 ) then
          write(6,*)
          efini = .5*(ektx(ieaord(ik+1))+ ektx(ieaord(ik)))
          if(autoew) then
Cstop2rx 2013.08.09 kino            if(ik<3) stop ' efsimplef2: ik<3'
            if(ik<3) call rx( ' efsimplef2: ik<3')
            esmr  = ektx(ieaord(ik)) - ektx(ieaord(ik-1))
          endif
          ierr=0
        endif
        if( wwg > valn+1d0) ikx=ik
      enddo
c
Cstop2rx 2013.08.09 kino      if(ierr==1) stop ' efsimplef2: ierr=1 given nval is too large'
      if(ierr==1) call rx( ' efsimplef2: ierr=1 given nval is too large')
c
      nbnqnsp = nband*nqibz*ispin
ccc gaussian
      if(GaussSmear) then
        valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,efini,esmr,nbnqnsp)
      else
        valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,efini,esmr,nbnqnsp)
      endif
      write(6,*) 'valx at efini=',efini,valx
      if(abs(valx-valn)<1d-8) then
        ef=efini
        goto 8891
      endif

      ew1= ektx(ieaord(1))-0.01d0
      ew2= ektx(ieaord(ikx))
      nbnqnsp =nband*nqibz*ispin
      do ix = 1,100
        ein = 0.5d0*(ew1+ew2)
        if(GaussSmear) then
          valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp))
     &          ,ektx(ieaord(1:nbnqnsp)),ein,esmr,nbnqnsp)
        else
          valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,ein,esmr,nbnqnsp)
        endif
        if(valx>valn) ew2=ein
        if(valx<valn) ew1=ein
        if(abs(ew1-ew2)<1d-15) exit
      enddo
      ef = 0.5d0*(ew1+ew2)
c
 8891 continue 
      if (mpi__root) then 
        write(8301,*) " ef=",ef
        close(8301)
        write(8302,*) " ef=",ef
        close(8302)
      endif 

c      write(6,*)' esmr        =',esmr
c      write(6,*)' determined ef =',ef
c----------------------------------
c      wwg = 0d0
c      do ik  = 1, nband*nqibz*ispin
c       wwgo = wwg
c       wwg  = wwg + wgtx(ieaord(ik))
c       if( abs(wwg-valn)<1d-6) then
c         ef = 0.5d0*( ektx(ieaord(ik))+ektx(ieaord(ik+1)) )
c         efp = ef + 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
c         efm = ef - 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efm is just below the fermi
c        elseif(wwg>valn) then
c         ef      = ektx(ieaord(ik))
c         wfacef  = (valn-wwgo)/wgtx(ieaord(ik))
c         efp = ef + 0.5d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
c         efm = ef - 0.5d0*(ektx(ieaord(ik))-ektx(ieaord(ik-1)))  !efm is just below the fermi
c          write(6,*)' determined ef    =',ef
c          write(6,*)'            efp   =',efp
c          write(6,*)'            efm   =',efm
c          write(6,*)'           wfacef =',wfacef
c          return
c        endif
c     enddo
      if(GaussSmear) then
        write(6,*)' efsimplef2(gauss):end'
      else
        write(6,*)' efsimplef2:end'
      endif
      end
c----------------------------------------------------------------------------

      subroutine efsimplef2 (ifev,ispin,wibz,qibz,ginv,
c      subroutine efsimplef(ifev,ispin,wibz,indxk, qibz
c     i             ,n1,n2,n3,qbas,ginv, nband,nqibz,
     i             nband,nqibz,
     i             konfig,z,nl,natom,iclass,nclass, 
     i             valn, legas, esmr, !input only for empty case
c
     i             index_qbz, n_index_qbz,qbz,nqbz,
     o             ef)
c Calculate efermi for discrete sum. (not for tetrahedron method)
c ifev(2) = direct access unit file for eigenvalues
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level
c nband   = no. states
c nqbz    = no. k-points
c valn    = number of valence electron.
c
c -------------------
c      e(iband) < efm : occupation is one
c efm< e(iband) < efp : occupation is wfacef.
c efp< e(iband)       : occupation is zero
c
c
      use m_mpi, only: mpi__root 
      implicit none
      integer(4):: is,iq,ispin,nqibz,ik,isig,kpx,nband,ifev(2)
      integer(4):: ieaord(nband*nqibz*ispin),mbytes,mwords,iwksize,
     &             iqibz  !,iindxk
      real(8)   :: ekt(nband, nqibz,ispin), ektx(nband*nqibz*ispin)
      real(8)   :: wgt(nband, nqibz,ispin), wgtx(nband*nqibz*ispin)
      real(8)   :: qbzx(3), qx(3),qbas(3,3),ginv(3,3),wwg
c
      integer(4):: nclass,natom,nl,ncore,l,ia,ic !,indxk(*)
      real(8)   :: wibz(nqibz),valn,ef, z(nclass),qibz(*)
      integer(4):: iclass(natom),konfig(0:nl-1,nclass),ierr
c
      integer(4) :: nbnqnsp,ix,ikx,ikini,nne
      real(8)    :: ew1,ew2,ein,valx,enumef_gauss,esmr, efini
     &            ,eee2,wwg2 ,enumef
c      real(8) :: efp,efm,wwgo,wfacef
      logical :: legas,autoew,GaussSmear=.true. !is external

      integer(4):: n_index_qbz,nqbz,ifile_handle,if8301,if8302
      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
      real(8)   :: qbz(3,nqbz)
C--------------------------------------------------------------------
      autoew =.false.
      if(GaussSmear) then
        write(6,*)' efsimplef2(gaussian mode):start'
      else
        write(6,*)' efsimplef2(rectangular mode):start'
      endif
      if(esmr<=0d0) autoew=.true.
c total valence charge
      if(legas) then
        write(6,*)' efsimplef2: legas=T use given valn = ',valn
      else
        valn    = 0d0
        do ia   = 1,natom
          ic    = iclass(ia)
          valn  = valn + z(ic)
          do    l = 0,nl-1
            valn  = valn - (konfig(l,ic)-l-1) *( 2*l +1)*2
          end do
        end do
        write(6,*)' valn=',valn
      endif
c
      do is = 1,ispin
        do iq = 1,nqibz
          call catch1 (qibz,iq,3,nqibz,qx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        call fbz  ( qx,n1,n2,n3,qbas,ginv,indxk,
c     o              qbzx,kpx)
          call fbz2(qx,  ginv,index_qbz,n_index_qbz,qbz,nqbz,
     o              qbzx,kpx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          call rwdd1( ifev(is),   !read eigenvalues
     d              kpx,nband,
     o              ekt(1,iq,is) )
        enddo
      enddo
c
      if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
        write(6,*) 'sum (wibz)=', sum(wibz(1:nqibz))
Cstop2rx 2013.08.09 kino        stop 'efsimplef2: wibzsumerr'
        call rx( 'efsimplef2: wibzsumerr')
      endif
      do is = 1,ispin
        do iq = 1,nqibz
          wgt(1:nband,iq,is) = wibz(iq)
          if(ispin==2) wgt(1:nband,iq,is) = wgt(1:nband,iq,is)/2d0
        enddo
      enddo

c ekt and wgt
      call dcopy ( nband*nqibz*ispin, ekt,1, ektx ,1)
      call dcopy ( nband*nqibz*ispin, wgt,1, wgtx ,1)
      call sortea( ektx,ieaord,nband*nqibz*ispin,isig)
c
c      write(6,*)nband,nqibz,ispin,nband*nqibz*ispin
c      do ik  = 1, nband*nqibz*ispin
c        write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wgtx(ieaord(ik))
c      enddo
c
      if (mpi__root) then 
         if8301=ifile_handle()
         if8302=ifile_handle()
         open(if8301,file = "DOSACC.lda")
         open(if8302,file = "DOSACC2.lda")
      endif 
      wwg = 0d0
      wwg2= 0d0
      eee2= -1d99
      ikini= 1
      ierr= 1
      nne = nband*nqibz*ispin
      do ik = 1, nne
c---
        if(eee2 +1d-4< ektx(ieaord(ik)).or.ik==nne ) then
! degeneracy check
          if (mpi__root) then 
            if(ik/=1) write(if8302,"(2i6,2d23.15)") ikini,ik-1,eee2,wwg2
          endif 
          wwg2 = wgtx(ieaord(ik))
          eee2 = ektx(ieaord(ik))
          ikini =ik
        else
          wwg2= wwg2 + wgtx(ieaord(ik))
        endif
c---
        wwg = wwg + wgtx(ieaord(ik))
        if(wwg<valn+2d0) write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg
c
        write(if8301,"(2i6,3d23.15)")
     &  ik,ieaord(ik),ektx(ieaord(ik)),wwg,wgtx(ieaord(ik))
c
        if( wwg>valn-1d-8.and.ierr==1 ) then
          write(6,*)
          efini = .5*(ektx(ieaord(ik+1))+ ektx(ieaord(ik)))
          if(autoew) then
Cstop2rx 2013.08.09 kino            if(ik<3) stop ' efsimplef2: ik<3'
            if(ik<3) call rx( ' efsimplef2: ik<3')
            esmr  = ektx(ieaord(ik)) - ektx(ieaord(ik-1))
          endif
          ierr=0
        endif
        if( wwg > valn+1d0) ikx=ik
      enddo
c
Cstop2rx 2013.08.09 kino      if(ierr==1) stop ' efsimplef2: ierr=1 given nval is too large'
      if(ierr==1) call rx( ' efsimplef2: ierr=1 given nval is too large')
c
      nbnqnsp = nband*nqibz*ispin
ccc gaussian
      if(GaussSmear) then
        valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,efini,esmr,nbnqnsp)
      else
        valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,efini,esmr,nbnqnsp)
      endif
      write(6,*) 'valx at efini=',efini,valx
      if(abs(valx-valn)<1d-8) then
        ef=efini
        goto 8891
      endif

      ew1= ektx(ieaord(1))-0.01d0
      ew2= ektx(ieaord(ikx))
      nbnqnsp =nband*nqibz*ispin
      do ix = 1,100
        ein = 0.5d0*(ew1+ew2)
        if(GaussSmear) then
          valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp))
     &          ,ektx(ieaord(1:nbnqnsp)),ein,esmr,nbnqnsp)
        else
          valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,ein,esmr,nbnqnsp)
        endif
        if(valx>valn) ew2=ein
        if(valx<valn) ew1=ein
        if(abs(ew1-ew2)<1d-15) exit
      enddo
      ef = 0.5d0*(ew1+ew2)
c
 8891 continue 
      if (mpi__root) then 
        write(if8301,*) " ef=",ef
        close(if8301)
        write(if8302,*) " ef=",ef
        close(if8302)
      endif

c      write(6,*)' esmr        =',esmr
c      write(6,*)' determined ef =',ef
c----------------------------------
c      wwg = 0d0
c      do ik  = 1, nband*nqibz*ispin
c       wwgo = wwg
c       wwg  = wwg + wgtx(ieaord(ik))
c       if( abs(wwg-valn)<1d-6) then
c         ef = 0.5d0*( ektx(ieaord(ik))+ektx(ieaord(ik+1)) )
c         efp = ef + 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
c         efm = ef - 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efm is just below the fermi
c        elseif(wwg>valn) then
c         ef      = ektx(ieaord(ik))
c         wfacef  = (valn-wwgo)/wgtx(ieaord(ik))
c         efp = ef + 0.5d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
c         efm = ef - 0.5d0*(ektx(ieaord(ik))-ektx(ieaord(ik-1)))  !efm is just below the fermi
c          write(6,*)' determined ef    =',ef
c          write(6,*)'            efp   =',efp
c          write(6,*)'            efm   =',efm
c          write(6,*)'           wfacef =',wfacef
c          return
c        endif
c     enddo
      if(GaussSmear) then
        write(6,*)' efsimplef2(gauss):end'
      else
        write(6,*)' efsimplef2:end'
      endif
      end

c------------------------------------------------------
      real(8) function enumef_gauss( wgtx,ektx,ein,esmr,nbnqnsp)
      implicit none
      integer(4):: nbnqnsp,ik
      real(8) :: ektx(nbnqnsp),wgtx(nbnqnsp),wwg,
     & derfcx,ein,esmr
      wwg = 0d0
      do ik = 1, nbnqnsp
        wwg= wwg + wgtx(ik)
     &    *0.5d0* derfcx( -(ein-ektx(ik))/sqrt(2d0)/esmr )
      enddo
      enumef_gauss = wwg
c     write(6,*)' ein enumef=', ein, enumef
      end

      real(8) function derfcx(a)
      real(8):: a, derfc,ax, amx=12d0
      ax = a
      if( abs(a)>amx) ax= a/abs(a) * amx
c      write(6,*)' xxx ',ax
c      write(6,*)' yyy ',derfc(ax)
      derfcx=derfc(ax)
      end

c------------
      real(8) function enumef( wgtx,ektx,ein,esmr,nbnqnsp)
      implicit real*8(a-h,o-z)
      real(8) :: ektx(nbnqnsp),wgtx(nbnqnsp)
c     write(6,*) esmr
      wwg = 0d0
      do ik = 1, nbnqnsp
c       write(6,*)'ik=',ik,ektx(ik),wgtx(ik)
        if    (  ektx(ik) + 0.5d0*esmr < ein ) then
          wwg  = wwg + wgtx(ik)
        elseif(  ektx(ik) - 0.5d0*esmr < ein ) then
          wwg  = wwg + wgtx(ik)*(ein- (ektx(ik)-0.5d0*esmr))/esmr
        endif
      enddo
      enumef = wwg
c     write(6,*)' ein enumef=', ein, enumef
      end



c---------------------------------------------------------------
      subroutine findemaxmin(nband,qbz,nqbz,nspin,
     o   emax,emin)
      use m_readeigen, only: readeval
      implicit none
      integer(4) :: nband,nqbz,nspin,isp,kx,i !,ifev(2)
      real(8)::emax,emin,qbz(3,nqbz),eee
      real(8),allocatable:: ekxxx(:,:,:)
      allocate( ekxxx(nband,nqbz,nspin))
      Emax=-1d9
      do isp =1, nspin
        do kx = 1, nqbz
c        call rwdd1 (ifev(isp), kx, nband, ekxxx(1:nband,kx,isp) )
c          write(6,"('readeval in finemaxmin',3f10.6,i8,i3,255f12.8)") qbz(:,kx), kx,isp
          call readeval(qbz(:,kx), isp,ekxxx(1:nband,kx,isp) )
c          write(6,"('readeval in finemaxmin',255f12.8)")ekxxx(1:nband,kx,isp)
          do i=1,nband
            eee= ekxxx(i,kx,isp)
c            print *,i,eee 
            if(eee>Emax.and.eee<1d9) Emax=eee !not eee<1d9 corresponds to 1d20 for padding in lmf2gw.F and sugw.Fago
          enddo
        enddo
      enddo
c      Emax = maxval(ekxxx(1:nband)
      Emin = minval(ekxxx)
      deallocate(ekxxx)
      end
c$$$      subroutine findemaxmin(ifev,nband,nqbz,nspin,
c$$$     o   emax,emin)
c$$$      implicit none
c$$$      integer(4) :: nband,nqbz,nspin,isp,kx,ifev(2)
c$$$      real(8)::emax,emin
c$$$      real(8),allocatable:: ekxxx(:,:,:)
c$$$      allocate( ekxxx(nband,nqbz,nspin))
c$$$      do isp =1, nspin
c$$$      do kx = 1, nqbz
c$$$        call rwdd1 (ifev(isp), kx, nband, ekxxx(1:nband,kx,isp) )
c$$$      enddo
c$$$      enddo
c$$$      Emax = maxval(ekxxx)
c$$$      Emin = minval(ekxxx)
c$$$      deallocate(ekxxx)
c$$$      end
      subroutine readomgc(ifinin,omg_c)
      real(8)::omg_c
      ifinin=iopen('GWIN_V2',1,0,0)
      read(ifinin,*)
      read(ifinin,*)
      read(ifinin,*) blank,omg_c !omg_c is frequency parameter
      ! freq2(iw)=dw(iw-1)+dw**2(iw-1)**2/2/omg_c
      ! quadratic term is essential for energy > omg_c
      is=iclose('GWIN_V2')
      end
