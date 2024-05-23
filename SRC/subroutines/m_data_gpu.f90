module m_data_gpu
  use m_zmel, only: ppbir
  use m_itq, only: ntq
  use m_read_bzdata,only: nqibz
  use m_zmel, only: ppovlz 
  use m_readVcoud, only: vcoud, ngb

  public :: SetDataGPU, ExitDataGPU, SetDataGPU_inkx, ExitDataGPU_inkx
  private
  logical :: has_ppbir_in_gpu = .false.
  logical :: has_ppovlz_in_gpu = .false.
  logical :: has_vcoud_in_gpu = .false.
contains
  subroutine SetDataGPU(set_ppbir_in_gpu)
    logical, intent(in), optional :: set_ppbir_in_gpu
    if(present(set_ppbir_in_gpu)) then
      if(set_ppbir_in_gpu .and. .not.has_ppbir_in_gpu) then !ppbir is set in mptauof_zmel
        block
          integer :: nlnmx, mdimx, nclass, ng, nspin
          nlnmx = size(ppbir,1)
          mdimx = size(ppbir,3)
          nclass = size(ppbir,4)
          ng = size(ppbir,5)
          nspin = size(ppbir,6)
          !$acc enter data copyin(ppbir(1:nlnmx,1:nlnmx,1:mdimx,1:nclass,1:ng,1:nspin))
          has_ppbir_in_gpu = .true.
        end block
      endif
    endif
  end subroutine
  subroutine ExitDataGPU()
    if(has_ppbir_in_gpu) then
      !$acc exit data delete(ppbir)
      has_ppbir_in_gpu = .false.
    endif
  end subroutine
  subroutine SetDataGPU_inkx(set_ppovlz_in_gpu, set_vcoud_in_gpu)
    logical, intent(in), optional :: set_ppovlz_in_gpu, set_vcoud_in_gpu
    if(present(set_ppovlz_in_gpu)) then
      if(set_ppovlz_in_gpu .and. .not.has_ppovlz_in_gpu) then  !ppovlz is set in setppovlz
        !$acc enter data copyin(ppovlz(1:ngb,1:ngb))
        has_ppovlz_in_gpu = .true.
      endif
    endif
    if(present(set_vcoud_in_gpu)) then
      if(set_vcoud_in_gpu .and. .not.has_vcoud_in_gpu) then !vcoud is set in readvcoud
        !$acc enter data copyin(vcoud(1:ngb))
        has_vcoud_in_gpu = .true.
      endif
    endif
  end subroutine
  subroutine  ExitDataGPU_inkx()
    if(has_ppovlz_in_gpu) then
      !$acc exit data delete(ppovlz)
      has_ppovlz_in_gpu = .false.
    endif
    if(has_vcoud_in_gpu) then
      !$acc exit data delete(vcoud)
      has_vcoud_in_gpu = .false.
    endif
  end subroutine
end module m_data_gpu
