program test_spinors
  use dirac_matrices
  implicit none
  real*8 :: p(4),pp(4),q(4),qt(4)
  real*8 :: E, mom, mass
  complex*16 :: norm

  mass = .938d0
  mom = .5d0
  E = sqrt(mom**2 + mass**2)

  norm = (0.0d0,0.0d0)
  p(1) = E
  p(2) = 0.0d0
  p(3) = 0.0d0
  p(4) = mom
  !p = (/E, 0.0, 0.0, mom/)
  pp = (/0.0, 0.0, 0.0, 0.0/)
  q = (/0.0,0.0,0.0,0.0/)
  qt = (/0.0,0.0,0.0,0.0/)
  
  call dirac_matrices_in(mass)
  call current_init(p,pp,q,qt)
  call define_spinors()
  call normalization()

end program 
