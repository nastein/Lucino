module two_nucleon_sf
  implicit none
  private
  public :: tn_sf_t, tn_sf_init, tn_sf_eval, tn_sf_norms, tn_sf_get_moms, tn_sf_xpf, tn_sf_FG_normalize

  real(8),parameter :: pi=acos(-1.0d0), hbarc=197.327053d0

  !========================
  ! Types
  !========================
  type tn_sf_t
     integer :: np = 80
     real(8) :: xpf_p = 0.0d0               ! Fermi momentum (same units as p)
     real(8) :: xpf_n = 0.0d0               ! Fermi momentum (same units as p)
     integer :: is_fg = 1                   ! 1 => build FG; 0 => read file
     real(8) :: dp_p_step = 0.0d0           ! proton grid spacing Δp (assumed uniform)
     real(8) :: dp_n_step = 0.0d0           ! neutron grid spacing Δp (assumed uniform)
     real(8), allocatable :: p_p(:)         ! momentum grid protons(size np)
     real(8), allocatable :: p_n(:)         ! momentum grid neutrons(size np)
     real(8), allocatable :: dp_pp(:,:)     
     real(8), allocatable :: dp_np(:,:)     
     real(8), allocatable :: dp_pn(:,:)     
     real(8), allocatable :: dp_nn(:,:)    
     real(8) :: norm_pp, norm_nn, norm_np, norm_pn 
  end type tn_sf_t

contains

  !========================
  ! Initialization
  !========================
  subroutine tn_sf_init(sf, use_fg, xpf_p, xpf_n, np, filename)
    type(tn_sf_t), intent(inout) :: sf
    integer,       intent(in)    :: use_fg
    real(8),       intent(in)    :: xpf_p, xpf_n
    integer,       intent(out)   :: np     
    character(*),  intent(in),   optional :: filename ! required if .not. use_fg

    integer :: i,j
    real(8) :: hp_p,hp_n
    real(8) :: dummy
    character(len=256) :: fname

    sf%xpf_p  = xpf_p
    sf%xpf_n  = xpf_n
    sf%is_fg = use_fg

    if (use_fg.eq.0) then
       ! ---------- FILE MODE ----------
       if (.not. present(filename)) then
          write(*,*) 'tn_sf_init: ERROR: filename is required in file mode.'
          stop 1
       end if

       fname = filename
       open(unit=8, file=fname, status='old', form='formatted', action='read')
       read(8,*) np
       sf%np = np

       allocate(sf%p_p(sf%np),sf%p_n(sf%np), &
       &  sf%dp_pp(sf%np,sf%np), sf%dp_np(sf%np,sf%np), &
       & sf%dp_pn(sf%np,sf%np), sf%dp_nn(sf%np,sf%np))

       do i=1,sf%np
          do j=1,sf%np
             read(8,*) sf%p_p(i),sf%p_p(j),dummy,sf%dp_pp(i,j), sf%dp_np(i,j)

             sf%dp_pn(i,j) =sf%dp_np(i,j)

             sf%p_n(i)=sf%p_p(i)
          end do
       end do
       close(8)

       ! p <- p*hbarc ; dp <- dp / hbarc^6 / (2π)^6
       sf%p_p   = sf%p_p * hbarc
       sf%p_n   = sf%p_n * hbarc

       sf%dp_pp = sf%dp_pp / hbarc**6
       sf%dp_np = sf%dp_np / hbarc**6
       sf%dp_pn = sf%dp_pn / hbarc**6

       sf%dp_pp = sf%dp_pp / (2.0d0*pi)**6
       sf%dp_np = sf%dp_np / (2.0d0*pi)**6
       sf%dp_pn = sf%dp_pn / (2.0d0*pi)**6

       sf%dp_nn = sf%dp_pp
       sf%dp_pn = sf%dp_np

       sf%dp_p_step = sf%p_p(2) - sf%p_p(1)
       sf%dp_n_step = sf%p_n(2) - sf%p_n(1)

    else
       ! ---------- FERMI-GAS MODE ----------       
       allocate(sf%p_p(sf%np),sf%p_n(sf%np), &
       &  sf%dp_pp(sf%np,sf%np), sf%dp_np(sf%np,sf%np), &
       & sf%dp_pn(sf%np,sf%np), sf%dp_nn(sf%np,sf%np))

       !Use default number of bins
       np = sf%np

       hp_p = xpf_p / dble(sf%np)
       hp_n = xpf_n / dble(sf%np)

       sf%dp_p_step = hp_p
       sf%dp_n_step = hp_n

       do i=1,sf%np
          sf%p_p(i) = dble(i-0.5d0)*hp_p
          sf%p_n(i) = dble(i-0.5d0)*hp_n
          do j=1,sf%np
             sf%dp_pp(i,j) = 1.0d0
             sf%dp_np(i,j) = 1.0d0
             sf%dp_pn(i,j) = 1.0d0
             sf%dp_nn(i,j) = 1.0d0
          end do
       end do
    end if

  end subroutine tn_sf_init

  !========================
  ! Get a p1,p2 depending on isospin 
  !========================
  subroutine tn_sf_get_moms(sf, i, j, t1, t2, p1, p2)
    type(tn_sf_t), intent(in)  :: sf
    real(8),      intent(out)  :: p1, p2
    integer,       intent(in) :: i, j, t1, t2

    if(t1.eq.1) then 
      p1 = sf%p_p(i)
    else
      p1 = sf%p_n(i)
    endif

    if(t2.eq.1) then
      p2 = sf%p_p(j)
    else
      p2 = sf%p_n(j)
    endif
  end subroutine tn_sf_get_moms

  !========================
  ! Evaluate the SF depending on isospin 
  !========================
  subroutine tn_sf_eval(sf, i, j, t1, t2, value, norm)
    type(tn_sf_t), intent(in)  :: sf
    real(8),      intent(out)  :: value, norm
    integer,       intent(in) :: i, j, t1, t2

    if(t1.eq.1 .and. t2.eq.1) then
      value = sf%dp_pp(i,j)
      norm = sf%norm_pp
    endif
    if(t1.eq.1 .and. t2.eq.2) then
      value = sf%dp_pn(i,j)
      norm = sf%norm_pn
    endif
    if(t1.eq.2 .and. t2.eq.1) then
      value = sf%dp_np(i,j)
      norm = sf%norm_np
    endif
    if(t1.eq.2 .and. t2.eq.2) then
      value = sf%dp_nn(i,j)
      norm = sf%norm_nn
    endif
  end subroutine tn_sf_eval

  !========================
  ! Current norms (after normalization)
  !========================
  subroutine tn_sf_norms(sf, norm_pp, norm_np, norm_pn, norm_nn)
    type(tn_sf_t), intent(inout)  :: sf
    real(8),       intent(out) :: norm_pp, norm_np, norm_pn, norm_nn
    call compute_norms(sf, norm_pp, norm_np, norm_pn, norm_nn)
    sf%norm_pp = norm_pp
    sf%norm_nn = norm_nn
    sf%norm_pn = norm_pn
    sf%norm_np = norm_np
  end subroutine tn_sf_norms
  
  subroutine tn_sf_xpf(sf, i1, i2, xpf1, xpf2)
    type(tn_sf_t), intent(in)  :: sf
    integer, intent(in) :: i1,i2
    real(8), intent(out) :: xpf1,xpf2
    if(i1.eq.1) then
      xpf1 = sf%xpf_p
    else
      xpf1 = sf%xpf_n
    endif

    if(i2.eq.1) then
      xpf2 = sf%xpf_p
    else
      xpf2 = sf%xpf_n
    endif
  end subroutine tn_sf_xpf

  !========================
  ! Internal: normalization (exactly your recipe)
  !========================
  subroutine tn_sf_FG_normalize(sf)
    type(tn_sf_t), intent(inout) :: sf
    real(8) :: norm_pp, norm_np, norm_pn, norm_nn
    real(8) :: V_p,V_n

    call compute_norms(sf, norm_pp, norm_np, norm_pn, norm_nn)

    ! Target normalization factors:
    !   multiply so that ∑ dp(i,j) p_i^2 p_j^2 (4πΔp)^2 = (4π xpf^3/3)^2 (and same for dp1, dp0)
    V_p = (4.0d0*pi*sf%xpf_p**3/3.0d0)    ! one-body k-space volume
    V_n = (4.0d0*pi*sf%xpf_n**3/3.0d0)    ! one-body k-space volume

    sf%dp_pp = sf%dp_pp / norm_pp * V_p * V_p
    sf%dp_np = sf%dp_np / norm_np * V_n * V_p
    sf%dp_pn = sf%dp_pn / norm_pn * V_p * V_n
    sf%dp_nn = sf%dp_nn / norm_nn * V_n * V_n

    call compute_norms(sf, norm_pp, norm_np, norm_pn, norm_nn)

    sf%norm_pp = norm_pp
    sf%norm_nn = norm_nn
    sf%norm_pn = norm_pn
    sf%norm_np = norm_np

  end subroutine tn_sf_FG_normalize

  !========================
  ! Internal: compute norms with weights
  !========================
  subroutine compute_norms(sf, norm_pp, norm_np, norm_pn, norm_nn)
    type(tn_sf_t), intent(in)  :: sf
    real(8),       intent(out) :: norm_pp, norm_pn, norm_np, norm_nn
    integer :: i,j

    norm_pp = 0.0d0
    norm_np = 0.0d0
    norm_pn = 0.0d0
    norm_nn = 0.0d0

    do i=1,sf%np
       do j=1,sf%np
          norm_pp = norm_pp + sf%dp_pp(i,j)*(4.0d0*pi)**2*sf%p_p(i)**2*sf%p_p(j)**2*sf%dp_p_step*sf%dp_p_step
          norm_np = norm_np + sf%dp_np(i,j)*(4.0d0*pi)**2*sf%p_n(i)**2*sf%p_p(j)**2*sf%dp_n_step*sf%dp_p_step
          norm_pn = norm_pn + sf%dp_pn(i,j)*(4.0d0*pi)**2*sf%p_p(i)**2*sf%p_n(j)**2*sf%dp_p_step*sf%dp_n_step
          norm_nn = norm_nn + sf%dp_nn(i,j)*(4.0d0*pi)**2*sf%p_n(i)**2*sf%p_n(j)**2*sf%dp_n_step*sf%dp_n_step
       end do
    end do
  end subroutine compute_norms

end module two_nucleon_sf
