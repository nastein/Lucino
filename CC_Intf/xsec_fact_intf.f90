    module xsec_fact
    implicit none
    real*8, parameter :: xmn=938.272d0,xmp=939.565d0,xme=0.0d0
    real*8, parameter :: hbarc=197.327053d0,G_F = 1.1664e-11*hbarc**2, cb=0.9741699d0,sb=sqrt(1.0d0-cb**2)
    real*8, parameter :: pi= acos(-1.0d0),c0=1.0d0/16.0d0/pi/pi,sw2=0.2224d0
    real*8, save :: xm,xm_in,xm_f, G_Ft

  end module xsec_fact


  subroutine masses_in(xmin, xmout)
    use xsec_fact
    implicit none
    real*8 :: xmin, xmout
    
    xm=(xmp+xmn)/2.0d0
    xm_in=xmin
    xm_f=xmout
    G_Ft=G_F*cb
    
  end subroutine masses_in

    !Computes nucleon cross section 
    subroutine sig_nuN_CCQE(pdel,pot_del,np_del,xq,w,wt,probeP4,outlepP4,p_4, &
        &   pp_4,pspec_4,ppspec_4,spec_iso,enu,contraction)
    use xsec_fact 
    use dirac_matrices 
    implicit none
    real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6
    real*8, parameter :: fstar=2.13d0,eps=1.0d0
    integer*4 :: ilept,i,spec_iso,np_del
    real*8 :: pdel(np_del), pot_del(np_del)
    real*8 :: xq,q2,w,wt,q2t,enu,thetalept
    real*8 :: elept,xklept,ep,ek
    real*8 :: p_4(4),pp_4(4),q_4(4),qt_4(4),xsec
    real*8 :: probeP4(4), outlepP4(4), pspec_4(4), ppspec_4(4)
    real*8 :: k2e_4(4),k1e_4(4)
    complex*16 :: had_del(4,4),had_pi(4,4),had_intf(4,4)
    complex*16 :: lepton_tensor(4,4), contraction
    real*8 :: re_contraction
    real*8 :: ca5,cv3,cv4,cv5
    real*8:: f1p,f1n,f2p,f2n,gep,gen,gmp,gmn
    real*8 :: ff1v,ff2v,ffa,ffp

    q2 = w**2 - xq**2
    q2t = wt**2 - xq**2

    q_4(1)=w
    q_4(2:3)=0.0d0
    q_4(4)=xq

    qt_4(1)=wt
    qt_4(2:3)=0.0d0 
    qt_4(4)=xq

    ! Define pion momenta
    k1e_4(:)=ppspec_4(:)-p_4(:) 
    k2e_4(:)= pp_4(:)-pspec_4(:)

    ! Initialize form factors for delta currents
    cv3=fstar/(1.0d0-q2t/lsq)**2/(1.0d0-q2t/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
    cv4=-1.15/(1.0d0-q2t/lsq)**2/(1.0d0-q2t/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
    cv5=0.48/(1.0d0-q2t/lsq)**2/(1.0d0-q2t/0.776/lsq)*sqrt(3.0d0/2.0d0)
    ca5=1.20d0/(1.0d0-q2t/xma2)**2/(1.0d0-q2t/3.0d0/xma2)*sqrt(3.0d0/2.0d0)

    ! Get the 1 body form factors
    call nform(-q2t/hbarc**2,f1p,f1n,f2p,f2n,gep,gen,gmp,gmn,ffa,ffp)
    ff1v=f1p-f1n
    ff2v=f2p-f2n

    ! Initialize all hadron kinematic variables
    call current_init(w,p_4,pspec_4,pp_4,ppspec_4,qt_4,k1e_4,k2e_4,2)
    call define_spinors()

    ! Initilalize lepton kinematic variables
    ! Compute leptonic tensor
    call lepton_current_init(probeP4, outlepP4)
    call define_lept_spinors()
    call lept_tens(lepton_tensor)

    !Compute 1 body currents
    call det_Ja(ff1v,ff2v,ffa,ffp)
   
    !Compute interference tensor from 2body pion and delta induced currents
    call det_Jpi(gep-gen)
    call det_JaJb_JcJd(cv3,cv4,cv5,ca5,np_del,pdel,pot_del)
    call det_J1Jdel_exc(had_del,-1,spec_iso)
    call det_J1Jpi_exc(had_pi,-1,spec_iso)
 
    !Hadron interference tensor
    had_intf = had_del + had_pi

    ! Compute |M|^2 
    ! Contraction of lepton and hadron tensor
    ! Taking minus sign in here from definition of matrix element
    contraction = contract(-had_intf,lepton_tensor)
    re_contraction = contraction + conjg(contraction)

    end subroutine sig_nuN_CCQE



