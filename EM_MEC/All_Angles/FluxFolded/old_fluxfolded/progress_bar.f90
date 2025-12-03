module progress_bar
   implicit none
   save

   integer*4 :: total_steps     = 0
   integer*4 :: last_step_print = -1
   integer*4 :: start_count     = 0
   integer*4 :: counts_per_sec  = 1
   real*8    :: last_pct_print  = -1.0d0
   real*8    :: min_delta_pct   = 1.0d0   ! only print if we advanced by ≥ this many percent

contains

   subroutine progress_init(total, delta_pct)
      implicit none
      integer*4, intent(in) :: total
      real*8,    intent(in), optional :: delta_pct
      integer*4 :: rate

      total_steps     = max(total, 1)
      last_step_print = -1
      last_pct_print  = -1.0d0

      if (present(delta_pct)) then
         min_delta_pct = max(delta_pct, 0.1d0)
      else
         min_delta_pct = 1.0d0
      end if

      call system_clock(start_count, rate)
      if (rate > 0) then
         counts_per_sec = rate
      else
         counts_per_sec = 1
      end if
   end subroutine progress_init


   subroutine progress_update(current_step)
      implicit none
      integer*4, intent(in) :: current_step
      integer*4 :: now, rate
      real*8 :: pct, elapsed, eta

      if (total_steps <= 0) return

      pct = 100.0d0 * dble(current_step) / dble(total_steps)
      if (pct > 100.0d0) pct = 100.0d0

      ! Only print if we've advanced at least min_delta_pct,
      ! unless we're at the very end.
      if (pct - last_pct_print < min_delta_pct .and. current_step /= total_steps) return

      call system_clock(now, rate)
      if (rate <= 0) then
         elapsed = 0.0d0
         eta     = 0.0d0
      else
         if (counts_per_sec <= 0) counts_per_sec = rate
         elapsed = dble(now - start_count) / dble(counts_per_sec)
         if (pct > 0.0d0) then
            eta = elapsed * (100.0d0/pct - 1.0d0)
         else
            eta = 0.0d0
         end if
      end if

      write(6,'("Progress: ",F6.2,"%", &
     &         " | elapsed ",F8.1," s", &
     &         " | ETA ",F8.1," s")') pct, elapsed, max(eta, 0.0d0)

      last_pct_print  = pct
      last_step_print = current_step
   end subroutine progress_update

end module progress_bar