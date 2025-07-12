module m_stopwatch !M.Obata 2024/04/21
  implicit none
  public :: stopwatch, stopwatch_start, stopwatch_init, stopwatch_pause, stopwatch_show, stopwatch_reset, &
            stopwatch_lap_time, stopwatch_elapsed_time, stopwatch_begin, stopwatch_show_list
  private
  type stopwatch
    private
    character(len=32) :: task_name = 'task'
    logical :: running = .false.
    real(8) :: elapsed_time = 0d0, stop_time = 0d0, start_time = 0d0, lap_time = 0d0
    ! ???? object-oriented style gets the execution error in intel depending on the case
    ! contains
    !   procedure :: start => stopwatch_start
    !   procedure :: pause => stopwatch_pause
    !   procedure :: show => stopwatch_show
    !   procedure :: setname => stopwatch_init
    !   procedure :: reset => stopwatch_reset
  end type stopwatch

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
    self%lap_time = 0d0
    self%stop_time = 0d0
  end subroutine stopwatch_init

  subroutine stopwatch_start(self)
    class(stopwatch), intent(inout) :: self
    call stopwatch_pause(self)
    self%start_time = measure_time()
    self%running = .true.
  end subroutine stopwatch_start

  subroutine stopwatch_begin(self, task_name)
    class(stopwatch), intent(inout) :: self
    character(len=*), intent(in) :: task_name
    call stopwatch_init(self, task_name)
    call stopwatch_start(self)
  end subroutine stopwatch_begin

  subroutine stopwatch_pause(self)
    class(stopwatch), intent(inout) :: self
    if (self%running) then
      self%stop_time = measure_time()
      self%lap_time = self%stop_time - self%start_time
      self%elapsed_time = self%elapsed_time + self%lap_time
      self%running = .false.
    endif
  end subroutine stopwatch_pause

  subroutine stopwatch_reset(self)
    class(stopwatch), intent(inout) :: self
    self%elapsed_time = 0d0
    self%lap_time = 0d0
    self%running = .false.
  end subroutine stopwatch_reset

  subroutine stopwatch_show(self)
    use m_lgunit, only:stdo
    class(stopwatch), intent(inout) :: self
    real(8) :: elapsed_time
    elapsed_time = stopwatch_elapsed_time(self)
    write(stdo,'(X,A,A20,X,F10.4,X,A)') 'Time:', trim(self%task_name), elapsed_time, '(sec)'
  end subroutine stopwatch_show

  subroutine stopwatch_show_list(self)
    use m_lgunit, only: stdo
    class(stopwatch), intent(inout) :: self(:)
    integer :: i
    real(8) :: elapsed_time
    do i = 1, size(self)
      elapsed_time = stopwatch_elapsed_time(self(i))
      write(stdo,'(X,A,A20,X,F10.4,X,A)') 'Time:', trim(self(i)%task_name), elapsed_time, '(sec)'
    enddo
  end subroutine stopwatch_show_list

  real(8) function stopwatch_elapsed_time(self) result(elapsed_time)
    class(stopwatch), intent(inout) :: self
    real(8) :: lap_time
    elapsed_time = self%elapsed_time
    if(self%running) then
      lap_time = stopwatch_lap_time(self)
      elapsed_time = elapsed_time + lap_time
    endif
  endfunction

  real(8) function stopwatch_lap_time(self) result(lap_time)
    class(stopwatch), intent(inout) :: self
    if (self%running) then
      self%lap_time = measure_time() - self%start_time
    endif
    lap_time = self%lap_time
  endfunction

  real(8) function measure_time() result(time)
#ifdef __GPU
    block
      use cudafor
      integer :: ierr
      ierr = cudadevicesynchronize() 
    end block
#endif 
    block
      !$ use omp_lib
      integer(kind=8) :: time_c, count_rate, count_max
      call system_clock(time_c, count_rate, count_max)
      time = dble(time_c)/count_rate
      ! call cpu_time(time)
      !$ time = omp_get_wtime()
    end block
  end function measure_time
end module m_stopwatch
