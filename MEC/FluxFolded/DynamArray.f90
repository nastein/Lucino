module dynvec_real8
  implicit none
  private
  public :: r8vec_t, r8vec_init, r8vec_push, r8vec_shrink, r8vec_free

  type :: r8vec_t
     real(8), allocatable :: a(:)
     integer :: n = 0
     integer :: cap = 0
  end type

contains

  subroutine r8vec_init(v, cap0)
    type(r8vec_t), intent(inout) :: v
    integer, intent(in), optional :: cap0
    integer :: c
    c = 100000
    if (present(cap0)) c = cap0
    v%n = 0
    v%cap = c
    allocate(v%a(v%cap))
  end subroutine

  subroutine r8vec_push(v, x)
    type(r8vec_t), intent(inout) :: v
    real(8), intent(in) :: x
    real(8), allocatable :: tmp(:)
    integer :: newcap

    if (v%n == v%cap) then
      newcap = max(2*v%cap, v%cap + 100000)   ! grow: double, or +chunk
      allocate(tmp(newcap))
      tmp(1:v%cap) = v%a(1:v%cap)
      call move_alloc(tmp, v%a)
      v%cap = newcap
    end if

    v%n = v%n + 1
    v%a(v%n) = x
  end subroutine

  subroutine r8vec_shrink(v)
    type(r8vec_t), intent(inout) :: v
    real(8), allocatable :: tmp(:)
    if (.not. allocated(v%a)) return
    allocate(tmp(v%n))
    if (v%n > 0) tmp = v%a(1:v%n)
    call move_alloc(tmp, v%a)
    v%cap = v%n
  end subroutine

  subroutine r8vec_free(v)
    type(r8vec_t), intent(inout) :: v
    if (allocated(v%a)) deallocate(v%a)
    v%n = 0; v%cap = 0
  end subroutine

end module