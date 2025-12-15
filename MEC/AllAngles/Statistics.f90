module select_tools
  use mympi
  implicit none
  private

  public :: quickselect_kth
  public :: compute_global_p99_gatherv
  public :: check_threshold
contains

  subroutine quickselect_kth(a, k)
    ! Reorders a so that a(k) is the kth smallest (1-based). a is modified.
    real(8), intent(inout) :: a(:)
    integer, intent(in) :: k
    integer :: left, right, pivotIndex

    left = 1
    right = size(a)

    do
      if (left == right) return
      pivotIndex = (left + right) / 2
      pivotIndex = partition(a, left, right, pivotIndex)

      if (k == pivotIndex) then
        return
      else if (k < pivotIndex) then
        right = pivotIndex - 1
      else
        left = pivotIndex + 1
      end if
    end do
  end subroutine

  integer function partition(a, left, right, pivotIndex) result(p)
    real(8), intent(inout) :: a(:)
    integer, intent(in) :: left, right, pivotIndex
    real(8) :: pivotValue, tmp
    integer :: i, storeIndex

    pivotValue = a(pivotIndex)

    ! move pivot to end
    tmp = a(pivotIndex); a(pivotIndex) = a(right); a(right) = tmp

    storeIndex = left
    do i = left, right-1
      if (a(i) < pivotValue) then
        tmp = a(storeIndex); a(storeIndex) = a(i); a(i) = tmp
        storeIndex = storeIndex + 1
      end if
    end do

    ! move pivot to its final place
    tmp = a(storeIndex); a(storeIndex) = a(right); a(right) = tmp
    p = storeIndex
  end function

  subroutine compute_global_p99_gatherv(comm, w_local, p, w_p)
    implicit none
    integer, intent(in) :: comm
    real(8), intent(in) :: w_local(:)   ! local warmup weights (suggest abs(w))
    real(8), intent(in) :: p            ! e.g. 0.99
    real(8), intent(out) :: w_p         ! global p-quantile estimate (exact order statistic)

    integer :: rank, nranks, ierr, root
    integer :: nloc, i
    integer, allocatable :: counts(:), displs(:)
    integer :: ntot, k
    real(8), allocatable :: w_all(:)

    root = 0
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nranks, ierr)

    nloc = size(w_local)

    if (rank == root) then
      allocate(counts(nranks), displs(nranks))
    else
      allocate(counts(1), displs(1))
    end if

    call MPI_Gather(nloc, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

    if (rank == root) then
      displs(1) = 0
      do i = 2, nranks
        displs(i) = displs(i-1) + counts(i-1)
      end do
      ntot = displs(nranks) + counts(nranks)
      allocate(w_all(ntot))
    else
      ntot = 0
      allocate(w_all(0))
    end if

    call MPI_Gatherv(w_local, nloc, MPI_DOUBLE_PRECISION, &
                     w_all, counts, displs, MPI_DOUBLE_PRECISION, &
                     root, comm, ierr)

    if (rank == root) then
      ! 1-based order statistic index for p-quantile:
      ! k = ceil(p * ntot)
      k = int(ceiling(p * dble(ntot)))
      if (k < 1) k = 1
      if (k > ntot) k = ntot

      call quickselect_kth(w_all, k)
      w_p = w_all(k)
    end if

    call MPI_Bcast(w_p, 1, MPI_DOUBLE_PRECISION, root, comm, ierr)
  end subroutine

  subroutine check_threshold(comm, w_local, wthr, p)
    implicit none
    integer, intent(in) :: comm
    real(8), intent(in) :: w_local(:)
    real(8), intent(in) :: wthr, p

    integer :: rank, ierr
    integer :: nloc, i
    integer(kind=8) :: n_over_loc, n_over, n_tot_loc, n_tot
    real(8) :: max_loc, max_glob
    real(8) :: max_over_loc, max_over_glob

    call MPI_Comm_rank(comm, rank, ierr)

    nloc = size(w_local)
    n_tot_loc = int(nloc, kind=8)

    max_loc = maxval(w_local)

    n_over_loc = 0_8
    max_over_loc = -1.0d0
    do i = 1, nloc
      if (w_local(i) > wthr) then
        n_over_loc = n_over_loc + 1_8
        if (w_local(i) > max_over_loc) max_over_loc = w_local(i)
      end if
    end do

    call MPI_Allreduce(n_tot_loc,  n_tot,       1, MPI_INTEGER8,        MPI_SUM, comm, ierr)
    call MPI_Allreduce(n_over_loc, n_over,      1, MPI_INTEGER8,        MPI_SUM, comm, ierr)
    call MPI_Allreduce(max_loc,    max_glob,    1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    call MPI_Allreduce(max_over_loc, max_over_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

    if (rank == 0) then
      write(*,'(a,1pe14.6)') "w_thr = ", wthr
      write(*,'(a,1pe14.6)') "global max = ", max_glob
      write(*,'(a,i0,a,i0,a,f10.6)') "N_over = ", n_over, " / N_tot = ", n_tot, &
           "  frac = ", dble(n_over)/dble(n_tot)
      write(*,'(a,1pe14.6)') "max(weight > w_thr) = ", max_over_glob
      write(*,'(a,f10.6)') "expected overflow frac ~ ", (1.0d0 - p)
    end if
  end subroutine

end module