program test_ff
  implicit none
  real*8 :: Q2,ff1s,ff2s,ff1v,ff2v,ges,gms,gev,gmv,Q2step
  real*8 :: ff1p, ff1n, ff2p, ff2n, hc, ffa, ffsa, ffp
  integer*4 :: ig
   character*50 :: fname

  hc = .197327d0
  Q2step = 2.0d0/1000.0


  fname='kellyFF_fortran.out'
  fname=trim(fname)
  open(unit=7, file=fname)

  Q2 = 0.0d0
  do ig=0,999
    Q2 = Q2 + Q2step
    call nform(Q2/hc/hc,ff1p,ff1n,ff2p,ff2n,ffa,ffp)
    ff1v=ff1p - ff1n
    ff2v=ff2p - ff2n
    write(7,*) Q2, ff1v, ff2v, ffa, ffp
  enddo

end program 