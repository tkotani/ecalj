!> Antiferro condition module. We have line AFcond at the bottom of 'LMTO' file.
!! Currently(feb2016), only laf is used (thus AF symmetry is not used yet for hx0fp0_sc)
!! To access laf, need to call anfcond() in advance.

      module m_anf
      implicit none
      logical,protected:: laf !! - laf: antiferro switch
      integer,allocatable,protected:: ibasf(:) !! - ibasf(ibas) specify AF pair atom.

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
      write(6,*) 'Read AFcond section in LMTO file, call anfcond in m_anf.F:'
      ilmto=ifile_handle()
      open(ilmto,file='LMTO')
      do 
        read(ilmto,"(a)",end=1011,err=1011) aaa
        aaa = adjustl(aaa)
        if(trim(aaa)=='number of atoms (natom)') then
          read(ilmto,*) natom
          read(ilmto,*)
          allocate(iantiferro(natom),ibasf(natom))
        endif
        if(aaa(1:6)=='AFcond') then
          read(ilmto,*) iantiferro(1:natom)
          ibasf=-999
          do ibas=1,natom
            do ibasx=ibas+1,natom
              if(abs(iantiferro(ibas))/=0 .and. iantiferro(ibas)+iantiferro(ibasx)==0) then
                ibasf(ibas)=ibasx
                exit
              endif
            enddo
            if(ibasf(ibas)/=-999) write(6,"(a,2i5)")' AF pair: ibas ibasf(ibas)=',ibas,ibasf(ibas)
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
      if(laf) write(6,"(a,100i4)") ' Antiferromode=',iantiferro
      end subroutine anfcond
      end module

