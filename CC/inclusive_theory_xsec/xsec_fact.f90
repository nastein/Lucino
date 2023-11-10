    module xsec_fact
    implicit none
    real*8, parameter :: xmn=938.272d0,xmp=939.565d0,xme=0.0d0,xmmu=105.658357
    !xmLambda=1115.70d0
    ![hbarc] = MeV * fm 
    ![G_F] = 
    real*8, parameter :: hbarc=197.327053d0,G_F = 1.1664e-11*hbarc**2, cb=0.9741699d0,sb=sqrt(1.0d0-cb**2)
    real*8, parameter :: pi= acos(-1.0d0),c0=1.0d0/16.0d0/pi/pi,sw2=0.2224d0
    integer*4, save :: iproc,ip,im
    real*8, save :: Alept(5),resp(5)
    real*8, save :: xm,xm_in,xm_f, G_Ft

  end module xsec_fact


  subroutine masses_in(iproc_in,ip_in,im_in)
    use xsec_fact
    implicit none
    integer*4 :: iproc_in,ip_in,im_in
    iproc=iproc_in
    ip=ip_in
    im=im_in
    
    if(iproc.ne.1)then 
       if(im.eq.1) then !CC nu scattering
          xm=(xmp+xmn)
          xm_in=xmn
          xm_f=xmp
          G_Ft=G_F*cb
       else if(im.eq.2) then !!CC anti-nu scattering
          xm=(xmp+xmn)
          xm_in=xmp
          xm_f=xmn
          G_Ft=G_F*cb
       endif
    elseif(iproc.eq.1.and.ip.eq.1)then !NC and scattering off of protons
       xm=(xmp+xmp)
       xm_in=xmp
       xm_f=xmp
       G_Ft=G_F*cb
    elseif(iproc.eq.1.and.ip.eq.2)then !NC and scattering off of neutrons
       xm=(xmn+xmn)
       xm_in=xmn
       xm_f=xmn
       G_Ft=G_F*cb
    endif
    
  end subroutine masses_in

    !Computes nucleon cross section 
    subroutine sig_nuN_CCQE(xmlept,xq,w,wt,xk,xp,enu,thetalept,xsec)
    use xsec_fact 
    use dirac_matrices 
    implicit none
    integer*4, parameter :: nphi=20
    integer*4 :: ilept,i
    real*8 :: xq,q2,w,wt,q2t,enu,thetalept
    real*8 :: xk,xp,sig(5)
    real*8 :: sin2,cos2,tan2,cosa,sina2,elept,xklept,xmlept,c1,c2,ep,ek
    real*8 :: xk_x,xk_z, xk_y, sig0
    real*8 :: aux1,c3,xkx2,xkz2,cont, ctot, hphi,phi
    real*8 :: p_4(4),pp_4(4),q_4(4),qt_4(4),qp(4),xsec
    real*8 :: v_cc,v_cl,v_ll,v_t,v_tp,kdotq
    real*8 :: nu !changes sign for anti neutrinos


    q2 = w**2 - xq**2
    !write(6,*)'q2 = ', q2
    q2t = wt**2 - xq**2

    elept = enu - w
    xklept = sqrt(elept**2 - xmlept**2)
    kdotq=(xmlept**2)/2.0d0+enu*w-(w**2-xq**2)/2.0d0
    !write(6,*)'kdotq = ', kdotq

    q_4(1)=w
    q_4(2:3)=0.0d0
    q_4(4)=xq

    qt_4(1)=wt
    qt_4(2:3)=0.0d0 
    qt_4(4)=xq


    !qp is k + k' as apposed to k - k'
    qp(1)=enu+elept
    qp(2)=2.0d0*sqrt(enu**2-kdotq**2/xq**2)
    qp(3)=0.0d0
    qp(4)=2.0d0*kdotq/xq-xq

    !Construct kinematic factors which multiply response tensor
    !Comes from multiplying Lmunu*Wmunu
    v_cc=qp(1)**2-xq**2-xmlept**2
    v_cl=-(-qp(1)*qp(4)+w*xq)
    v_ll=qp(4)**2-xq**2-q2+xmlept**2
    v_t=qp(2)**2/2.0d0-q2+xmlept**2
    v_tp=(qp(1)*xq-w*qp(4))

     !write(6,*) 'v_cc = ', v_cc, ', v_cl = ', v_cl, ', v_ll = ', v_ll, ', v_t = ', v_t, ', v_tp = ', v_tp

   
    ek = sqrt(xk**2 + xm_in**2)
    ep = sqrt(xp**2 + xm_f**2)

    if (xk.ne.0.d0) then 
       cosa = (xp**2-xk**2-xq**2)/2.0d0/xk/xq
       sina2 = 1.0d0-cosa**2
    else
       sig=0.0d0
       return
    endif

    hphi=2.0d0*pi/dble(nphi)
    sig=0.0d0
    sig0=(G_Ft)**2/(2.0d0*pi)*xklept/2.0d0/enu/hbarc**2
      !write(6,*)'sig0 = ', sig0
    !Integration over phi (not needed for QE because azimuthal symmetry)
    do i=1,nphi
       phi=(dble(i)-0.5d0)*hphi
       xk_x=xk*sqrt(sina2)*cos(phi)
       xk_y=xk*sqrt(sina2)*sin(phi)      
       xk_z=xk*cosa
       xkx2 = xk**2*sina2/2.0d0   
       xkz2= (xk*cosa)**2
       p_4(1)=ek
       p_4(2)=xk_x
       p_4(3)=xk_y
       p_4(4)=xk_z
       pp_4(1)=ep
       pp_4(2)=xk_x
       pp_4(3)=xk_y
       pp_4(4)=xk_z+xq
       call current_init(p_4, pp_4, q_4,qt_4)!initialize four momentum
       call define_spinors()!create four spinors
       call sf_hyp(q2t)!calculates form factors and response tensor
       sig = sig+0.5*resp*hphi !factor of 1/2 comes from averaging over spins
    enddo

    !Changes the sign of r_t' for anti-neutrinos
    nu=1.0d0
    if(im.ne.1)then
      nu=-1.0d0
    endif
    xsec=sig0*(v_cc*sig(1)+2.0d0*v_cl*sig(2)+v_ll*sig(3)+ &
         & v_t*sig(4)+2.0d0*nu*v_tp*sig(5))

    end subroutine sig_nuN_CCQE

    !This calculates all relavent form factors and then computes the response tensor and grabs the linear combinations that we need (r_cc, r_cl, etc.)
    subroutine sf_hyp(q2t)

    use xsec_fact
    use dirac_matrices
    implicit none
    real*8:: f1p, f1n, f2p,f2n, q2,pdotq,q2t
    real*8 :: ff1s,ff2s,ff1v,ff2v,ges,gms,gev,gmv,ffa,ffp,ffsa

    call nform(-q2t/hbarc**2,f1p,f1n,f2p,f2n,ffa,ffp)


    if(iproc.ne.1)then
       ff1v=f1p-f1n
       ff2v=f2p-f2n
    elseif(iproc.eq.1.and.ip.eq.1)then
       ff1v=(0.5d0-2.0d0*sw2)*f1p-0.5d0*f1n
       ff2v=(0.5d0-2.0d0*sw2)*f2p-0.5d0*f2n
       ffa=0.5d0*ffa+0.5d0*ffsa
    elseif(iproc.eq.1.and.ip.eq.2)then
       ff1v=(0.5d0-2.0d0*sw2)*f1n-0.5d0*f1p
       ff2v=(0.5d0-2.0d0*sw2)*f2n-0.5d0*f2p
       ffa=-0.5d0*ffa+0.5d0*ffsa
    endif

    !ffp = 0.0d0
    call det_Ja(ff1v, ff2v, ffa, ffp)
    call det_res1b(resp)

  end subroutine sf_hyp


