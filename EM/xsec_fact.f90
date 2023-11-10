    module xsec_fact
    implicit none
    real*8, parameter :: xmn=938.272d0,xmp=939.565d0,xme=0.0d0
    real*8, parameter :: hbarc=197.327053d0, alpha = 1.0d0/137.0d0, cb=0.9741699d0,sb=sqrt(1.0d0-cb**2)
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
    
  end subroutine masses_in

    !Computes nucleon cross section 
    subroutine sig_eN_EMQE(xq,w,wt,probeP4,outlepP4,p_4, &
        &   pp_4,enu,xsec)
    use xsec_fact 
    use dirac_matrices 
    implicit none
    integer*4 :: ilept,i
    real*8 :: xq,q2,w,wt,q2t,enu,thetalept
    real*8 :: elept,xklept,ep,ek
    real*8 :: sig
    real*8 :: p_4(4),pp_4(4),q_4(4),qt_4(4),xsec
    real*8 :: probeP4(4), outlepP4(4)


    q2 = w**2 - xq**2
    q2t = wt**2 - xq**2

    elept = enu - w

    xklept = sqrt(sum(outlepP4(2:4)**2))

    q_4(1)=w
    q_4(2:3)=0.0d0
    q_4(4)=xq

    qt_4(1)=wt
    qt_4(2:3)=0.0d0 
    qt_4(4)=xq

    sig=0.0d0
    
    call current_init(p_4, pp_4, q_4, qt_4, probeP4, outlepP4)
    call define_spinors()
    call define_lept_spinors()
    call sf_hyp(sig,q2t,q2)
    xsec = sig*xklept*elept
    

    end subroutine sig_eN_EMQE

    !This calculates all relavent form factors and then computes the response tensor and grabs the linear combinations that we need (r_cc, r_cl, etc.)
    subroutine sf_hyp(sig,q2t,q2)

    use xsec_fact
    use dirac_matrices
    implicit none
    real*8:: f1p,f1n,f2p,f2n,q2t,sig_p,sig,sigmott
    real*8 :: ff1v,ff2v,ffa,ffp,q2

    sigmott = (8.0d0*alpha**2)/(q2**2)*2.0d0*pi

    call nform(-q2t/hbarc**2,f1p,f1n,f2p,f2n,ffa,ffp)
   
    call det_Ja(f1p, f2p, ffa, ffp)
    call contract(sig_p)
    sig=sigmott*0.5d0*(sig_p)

  end subroutine sf_hyp


