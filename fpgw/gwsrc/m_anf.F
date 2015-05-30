      module m_anf
      implicit none
      logical:: laf !antiferro switch
      integer:: natom
      integer,allocatable:: ibasf(:),ldima(:),iantiferro(:),iclasst(:)
      real(8),allocatable:: pos(:,:),anfvec(:),qlat(:,:),plat(:,:)
      contains

      subroutine anfcond()
      implicit none
      integer:: ifile_handle,ilmto,ildima,ificlass
      character(256):: aaa,keyplat
      real(8)::vecs(3),vece(3),basdiff(3)
      integer:: ibas,lkeyplat,i,ibasx
      character(3)::iaaa
!! read LMTO file
      write(6,*) 'anfcond:'
      ilmto=ifile_handle()
      open(ilmto,file='LMTO')
      do 
        read(ilmto,"(a)",end=1011,err=1011) aaa
        aaa = adjustl(aaa)
c        print *,trim(aaa)
        if(trim(aaa)=='primitive lattice vectors (plat)') then
          allocate(plat(3,3),qlat(3,3))
          do i=1,3
           read(ilmto,*) plat(1:3,i)
          enddo
          call dinv33x (plat,qlat)
        endif  
        if(trim(aaa)=='number of atoms (natom)') then
          read(ilmto,*) natom
          read(ilmto,*)
          allocate(iantiferro(natom),pos(3,natom))
          iantiferro=0
          do ibas = 1,natom
            read(ilmto,*) pos(1:3,ibas)
            write(6,*) pos(1:3,ibas)
          enddo
        endif
        if(aaa(1:6)=='AFcond') then
          read(ilmto,*) iantiferro(1:natom)
        endif
      enddo
 1011 continue
      close(ilmto)
      if(sum(abs(iantiferro))==0) then
        laf=.false. !no AF case
        return
      endif  
!! Antiferro case --------------
      laf=.true.
      write(6,"(a,100i4)") ' Antiferromode=',iantiferro
!! Read ldima
      ildima=ifile_handle()
      open(ildima,file='ldima')
      read(ildima,'(a3)') iaaa
      close(ildima)
! Is this correct description? takao2015may
!     ... June2007 for floating orbitals
!     ldima is generaged by lmfgw; it contains number of MTO including
!     floating orbital, and positions.
      if(iaaa=='***') then
        open(ildima,file='ldima')
        read(ildima,*) iaaa,natom
        deallocate(pos)
        allocate(ldima(natom),pos(1:3,natom) )
        do ibas = 1,natom
          read(ildima,*) ldima(ibas),pos(1:3,ibas)
          write(6,"('ldima pos=',i5,3f10.4)")ldima(ibas),pos(1:3,ibas)
        enddo
        close(ildima)
      else
        open(ildima,file='ldima')
        allocate(ldima(natom))
        do ibas=1,natom
          read(ildima,*) ldima(ibas)
          write(6,*) 'ldima=',ldima(ibas)
        enddo
        close(ildima)
      endif
!! Read CLASS
      ificlass=ifile_handle()
      open (ificlass,file='CLASS')
      allocate(iclasst(natom))
      write(6,*)'  --- Readingin CLASS info ---'
      do ibas = 1,natom
        read(ificlass,*) ibasx, iclasst(ibas)
      enddo
      close(ificlass)

!! Get anfvec and ibasf
      allocate(ibasf(natom))
c          ifianf = 211
c          open(ifianf,file='ANFcond')
c          read(ifianf,*)
c          read(ifianf,*)
c          read(ifianf,*) anfvec(1:3)
!!
      do ibas=1,natom
        if(iantiferro(ibas)==-1) then
          vecs = pos(:,ibas)
          iclasst(ibas)=  999 !overwrite by 999
        elseif(iantiferro(ibas)==1) then
          vece = pos(:,ibas)
          iclasst(ibas)=  999
        endif
      enddo
!!
      allocate(anfvec(3))
      anfvec = vece-vecs
      do ibas=1,natom
        do ibasx=1,natom
          do i=1,3
            basdiff(i)= sum((pos(:,ibas)+anfvec-pos(:,ibasx))*qlat(:,i))
          enddo
c         write(6,"(a,4i4,3f13.6)")' ibas ibasx iclass iclassx basdiff=',ibas,ibasx,iclasst(ibas),iclasst(ibasx),basdiff
          if(sum(abs(basdiff-anint(basdiff)))<1d-6.and.iclasst(ibas)==iclasst(ibasx)) then
            ibasf(ibas)=ibasx
            write(6,"(a,2i5)")' ibas ibasf=',ibas,ibasf(ibas)
            goto 888
          endif
        enddo
        call rx('m_anf: ibasf did not found')
 888    continue
      enddo
      write(6,'(" antiferro translation vector=",3f13.6)') anfvec
      end subroutine anfcond
      end module

