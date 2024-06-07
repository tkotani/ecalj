module m_stopwatch !M.Obata 2024/04/21
  implicit none
  private
  public :: t_sw_zmel, t_sw_x0
  public :: stopwatch, stopwatch_start, stopwatch_init, stopwatch_pause, stopwatch_show, stopwatch_reset
  type stopwatch
    private
    character(len=32) :: task_name = 'task'
    logical :: running = .false.
    real(8) :: elapsed_time = 0d0, stop_time = 0d0, start_time = 0d0
    ! ???? object-oriented style gets the execuation error in intel depending on the case
    ! contains
    !   procedure :: start => stopwatch_start
    !   procedure :: pause => stopwatch_pause
    !   procedure :: show => stopwatch_show
    !   procedure :: setname => stopwatch_init
    !   procedure :: reset => stopwatch_reset
  end type stopwatch

  type(stopwatch) :: t_sw_zmel, t_sw_x0

  ! interface stopwatch
  !   module procedure new_stopwatch
  ! end interface stopwatch

  contains
  ! type(stopwatch) function new_stopwatch(task_name)
  !   character(len=*), intent(in), optional ::task_name 
  !   if(present(task_name)) call new_stopwatch%setname(task_name)
  ! end function new_stopwatch

  subroutine stopwatch_init(self, task_name)
    class(stopwatch), intent(inout) :: self
    character(len=*), intent(in) ::task_name 
    self%task_name = task_name
    self%running = .false.
    self%start_time = 0d0
    self%elapsed_time = 0d0
    self%stop_time = 0d0
  end subroutine stopwatch_init

  subroutine stopwatch_start(self)
    class(stopwatch), intent(inout) :: self
    if (self%running) then
      self%elapsed_time = self%elapsed_time + measure_time() - self%start_time
    endif
    self%start_time = measure_time()
    self%running = .true.
  end subroutine stopwatch_start

  subroutine stopwatch_pause(self)
    class(stopwatch), intent(inout) :: self
    if (self%running) then
      self%stop_time = measure_time()
      self%elapsed_time = self%elapsed_time + self%stop_time - self%start_time
      self%running = .false.
    endif
  end subroutine stopwatch_pause

  subroutine stopwatch_reset(self)
    class(stopwatch), intent(inout) :: self
    self%elapsed_time = 0d0
  end subroutine stopwatch_reset

  subroutine stopwatch_show(self)
    use m_lgunit, only:stdo
    class(stopwatch), intent(inout) :: self
    if (self%running) then
      self%stop_time = measure_time()
      self%elapsed_time = self%elapsed_time + self%stop_time - self%start_time
    end if
    write(stdo,'(X,A,A20,X,F10.4,X,A)') 'Time:', trim(self%task_name), self%elapsed_time, '(sec)'
  end subroutine stopwatch_show

  real(8) function measure_time() result(time)
    !$ use omp_lib
#ifdef __GPU
    use cudafor
    integer :: ierr
    ierr = cudadevicesynchronize() 
#endif 
    call cpu_time(time)
    !$ time = omp_get_wtime()
  end function measure_time
end module m_stopwatch
