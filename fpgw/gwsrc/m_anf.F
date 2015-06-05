      module m_anf
      implicit none
      logical:: laf !antiferro switch
      integer,allocatable:: ibasf(:)
c      integer:: natom
c     ,ldima(:),iantiferro(:),iclasst(:)
c      real(8),allocatable:: pos(:,:),anfvec(:),qlat(:,:),plat(:,:)
      contains

      subroutine anfcond()
      implicit none
      integer,allocatable:: iantiferro(:)
      integer:: ifile_handle,ilmto,ildima,ificlass
      character(256):: aaa,keyplat
      real(8)::vecs(3),vece(3),basdiff(3)
      integer:: ibas,lkeyplat,i,ibasx,natom
      character(3)::iaaa


!! read LMTO file
      write(6,*) 'anfcond:'
      ilmto=ifile_handle()
      open(ilmto,file='LMTO')
      do 
        read(ilmto,"(a)",end=1011,err=1011) aaa
        aaa = adjustl(aaa)
c$$$c        print *,trim(aaa)
c$$$        if(trim(aaa)=='primitive lattice vectors (plat)') then
c$$$          allocate(plat(3,3),qlat(3,3))
c$$$          do i=1,3
c$$$           read(ilmto,*) plat(1:3,i)
c$$$          enddo
c$$$          call dinv33x (plat,qlat)
c$$$        endif  
        if(trim(aaa)=='number of atoms (natom)') then
          read(ilmto,*) natom
          read(ilmto,*)
          allocate(iantiferro(natom),ibasf(natom))
c$$$          allocate(iantiferro(natom),pos(3,natom))
c$$$          iantiferro=0
c$$$          do ibas = 1,natom
c$$$            read(ilmto,*) pos(1:3,ibas)
c$$$            write(6,*) pos(1:3,ibas)
c$$$          enddo
        endif
        if(aaa(1:6)=='AFcond') then
          read(ilmto,*) iantiferro(1:natom)
          do ibas=1,natom
            do ibasx=ibas+1,natom
              if(abs(iantiferro(ibas))/=0 .and. iantiferro(ibas)+iantiferro(ibasx)==0) then
                ibasf(ibas)=ibasx
                exit
              endif
            enddo
            write(6,"(a,2i5)")' AF pair: ibas ibasf(ibas)=',ibas,ibasf(ibas)
          enddo
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
      end subroutine anfcond
      end module

!! ---- followings are wrong
c$$$!! Read ldima
c$$$      ildima=ifile_handle()
c$$$      open(ildima,file='ldima')
c$$$      read(ildima,'(a3)') iaaa
c$$$      close(ildima)
c$$$! Is this correct description? takao2015may
c$$$!     ... June2007 for floating orbitals
c$$$!     ldima is generaged by lmfgw; it contains number of MTO including
c$$$!     floating orbital, and positions.
c$$$      if(iaaa=='***') then
c$$$        open(ildima,file='ldima')
c$$$        read(ildima,*) iaaa,natom
c$$$        deallocate(pos)
c$$$        allocate(ldima(natom),pos(1:3,natom) )
c$$$        do ibas = 1,natom
c$$$          read(ildima,*) ldima(ibas),pos(1:3,ibas)
c$$$          write(6,"('ldima pos=',i5,3f10.4)")ldima(ibas),pos(1:3,ibas)
c$$$        enddo
c$$$        close(ildima)
c$$$      else
c$$$        open(ildima,file='ldima')
c$$$        allocate(ldima(natom))
c$$$        do ibas=1,natom
c$$$          read(ildima,*) ldima(ibas)
c$$$          write(6,*) 'ldima=',ldima(ibas)
c$$$        enddo
c$$$        close(ildima)
c$$$      endif
c$$$!! Read CLASS
c$$$      ificlass=ifile_handle()
c$$$      open (ificlass,file='CLASS')
c$$$      allocate(iclasst(natom))
c$$$      write(6,*)'  --- Readingin CLASS info ---'
c$$$      do ibas = 1,natom
c$$$        read(ificlass,*) ibasx, iclasst(ibas)
c$$$      enddo
c$$$      close(ificlass)
c$$$
c$$$!! Get anfvec and ibasf
c$$$      allocate(ibasf(natom))
c$$$c          ifianf = 211
c$$$c          open(ifianf,file='ANFcond')
c$$$c          read(ifianf,*)
c$$$c          read(ifianf,*)
c$$$c          read(ifianf,*) anfvec(1:3)
c$$$!!
c$$$      do ibas=1,natom
c$$$        if(iantiferro(ibas)==-1) then
c$$$          vecs = pos(:,ibas)
c$$$          iclasst(ibas)=  999 !overwrite by 999
c$$$        elseif(iantiferro(ibas)==1) then
c$$$          vece = pos(:,ibas)
c$$$          iclasst(ibas)=  999
c$$$        endif
c$$$      enddo
c$$$!!
c$$$      allocate(anfvec(3))
c$$$      anfvec = vece-vecs
c$$$      do ibas=1,natom
c$$$        do ibasx=1,natom
c$$$          do i=1,3
c$$$            basdiff(i)= sum((pos(:,ibas)+anfvec-pos(:,ibasx))*qlat(:,i))
c$$$          enddo
c$$$c         write(6,"(a,4i4,3f13.6)")' ibas ibasx iclass iclassx basdiff=',ibas,ibasx,iclasst(ibas),iclasst(ibasx),basdiff
c$$$          if(sum(abs(basdiff-anint(basdiff)))<1d-6.and.iclasst(ibas)==iclasst(ibasx)) then
c$$$            ibasf(ibas)=ibasx
c$$$            write(6,"(a,2i5)")' ibas ibasf=',ibas,ibasf(ibas)
c$$$            goto 888
c$$$          endif
c$$$        enddo
c$$$        call rx('m_anf: ibasf did not found')
c$$$ 888    continue
c$$$      enddo
c$$$      write(6,'(" antiferro translation vector=",3f13.6)') anfvec
c$$$      end subroutine anfcond
c$$$      end module

