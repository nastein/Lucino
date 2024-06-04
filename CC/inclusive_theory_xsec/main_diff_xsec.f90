program response_ia
    use mathtool
    use dirac_matrices

    implicit none
    real*8, parameter :: mp=938.00d0,mn=938.00d0,mu=931.494061d0,hbarc=197.327053d0
    real*8, parameter :: xmLambda=1115.70d0,xmSigma_m=1197.449d0,xm_Sigma_z=1192.642d0

    real*8, parameter :: xme=0.0d0,enu_max=5000.0d0,pi=acos(-1.0d0),xmmu=105.658357
    real*8, parameter :: small=1e-12

    integer*4 :: np0,ne0,nA,nZ,np,nc,nw,i,j,k,nskip,iform,ilept,iw,l,nq,iproc,ip,fg,im,ip_in,ip_max
    integer*4 :: ie,nen,ie0
    real*8, allocatable :: p(:),dp(:),cost_d(:),coste(:),dp0(:,:),ee0(:),p0(:),xnp(:),Pke(:,:)
    real*8, allocatable :: w(:),wt(:),sigma(:,:,:),Q2(:),my(:),mqe(:),sigma_tot(:,:),enu_v(:)
    real*8 :: enu,thetalept,cos_theta,jac_c,ep,he
    real*8 :: pmax,hp,wmin,wmax,hw,rdummy,mspec,mnuc,qval,qfm,eps1,eps2,sig,arg,mstar
    real*8 :: cost,cost_rel,pf,dpf,kf,TA,dwfold,ffold,hwfold,hc,delta_w,espec,epf
    real*8 :: q4qe,emu,pmu,delta_cos,arg1,hq,q2max,mlept,costheta_p,res,q2_p,np1,norm,emax
    character * 40 :: nk_fname
    character*50 :: fname,en_char,theta_char,int_char,dec_char
    logical :: explicit_d

! inputs

    read(5,*) enu,thetalept!q2max is given in GeV^2
    read(5,*) ilept
    read(5,*) fg
    read(5,*) iproc
    read(5,*) nZ, nA
    read(5,*) kF
    read(5,*) mstar
    read(5,*) nk_fname
    read(5,*) np0, nE0
    read(5,*) np,nw
    read(5,*) wmax
    read(5,*) wmin
    read(5,*) eps1,eps2

    !Formatting output file names
    write(int_char,'(i3)') int(thetalept)
    write(dec_char,'(i1)') int(mod(thetalept*10.0d0,10.0d0))

    theta_char=trim(int_char)//'p'//trim(dec_char)
    theta_char=adjustl(theta_char)
    theta_char=trim(theta_char)
    write(en_char,'(i4)') int(enu)
    en_char=adjustl(en_char)
    en_char=trim(en_char)

    fname='test_xsec_diff_nu_'//trim(en_char)//'_'//trim(theta_char)//'_CC.out'
    fname=trim(fname)

    if(iproc.ne.1) then !charged current
       ip_in=0
       ip_max=0
    elseif(iproc.eq.1)then !neutral current
       ip_in=1
       ip_max=2
    endif
    write(6,*) 'iproc = ',iproc, ', ip_in = ', ip_in, ', ip_max = ',ip_max
    thetalept=thetalept/180.0d0*pi

! initialize useful constants
    mspec=dble(nA-1)*mu !mass of spectator nucleus
    mnuc=dble(nA)*mu    !mass of initial nucleus
    if(iproc.ne.1)then  !charged current
       if(ilept.eq.1)then !muons
          mlept=xmmu
       else
          mlept=0.0d0
       endif
    elseif(iproc.eq.1)then !neutral current
       mlept=0.0d0
    endif
    
    allocate(mqe(2),mY(2))


! read and interpolate the momentum distribution
    if(1==0) then

      open(unit=10,file='Pke_c12_sep.out',status='unknown',form='formatted')
      read(10,*) np0, nE0
      allocate(p0(np0),dp0(ne0,np0),ee0(nE0),dp(np0))


      do i=1,ne0
         do j=1,np0
            read(10,*) p0(j), ee0(i), dp0(i,j)
         enddo
            read(10,*)
      enddo  
      p0 = p0*hbarc
      dp0=dp0/hbarc**3/(2.0d0*pi)**3
    
      close(10)
      he=(ee0(2)-ee0(1))
      do j=1,np0
         dp(j)=sum(dp0(:,j))*he
      enddo
    
      norm=sum(p0(:)**2*dp(:))*4.d0*pi*(p0(2)-p0(1))
      write(6,*) 'norm nk' ,norm
    
      dp0=dp0/norm
      pmax=p0(np0)
      emax=ee0(ne0)
      write(6,*)'pmax=',pmax,emax


      allocate(p(np))
      hp=p0(np0)/dble(np)
      do i=1,np
         p(i)=(dble(i)-0.5d0)*hp
      enddo
      write(6,*)'n(k) norm=',sum(dp(:)*p(:)**2)*4.0d0*pi*hp
    endif

    if(1==1) then
      open(unit=10,file='pke12_tot.data',status='unknown',form='formatted')
      allocate(p0(np0),dp0(nE0,np0),ee0(nE0),xnp(np0))
      do j=1,np0
         read(10,*) p0(j)
         read(10,'(4(f6.1,2x,e10.3))')(ee0(i),dp0(i,j),i=1,nE0)
         dp0(:,j)=dp0(:,j)/dble(nZ) ! omar c12, carlo
         !dp0(:,j)=dp0(:,j)/dble(nA) ! omar o16 

      enddo
      close(10)
      he=(ee0(2)-ee0(1))
      do j=1,np0
         xnp(j)=sum(dp0(:,j))*he
      enddo
      write(6,*) 'norm nk' ,sum(p0(:)**2*xnp(:))*4.d0*pi*(p0(2)-p0(1))
      allocate(p(np),dp(np))
      hp=p0(np0)/dble(np)
      do i=1,np
         p(i)=(dble(i)-0.5d0)*hp
         if(p(i).gt.kf) then
            dp(i)=0.0d0
         else
            dp(i)=3.0d0/(4.0d0*pi*kf**3)
         endif
      enddo
      write(6,*)'n(k) norm=',sum(dp(:)*p(:)**2)*4.0d0*pi*hp
    endif

! construct omega grid and the form factors
    print *, "Constructing omega grid and form factors"
    allocate(w(nw),wt(nw))

    wmin=0.0d0
    hw=(wmax-wmin)/dble(nw)
    do iw=1,nw
       w(iw)=wmin + (dble(iw))*hw
       !w(iw) = 120.0d0
    end do
    if(fg.eq.1) then
            nE0=1
            he=1.0d0
    endif    

    call dirac_matrices_in(mn)    

       ! compute the response functions in the impulse approximation
    allocate(sigma(2,3,nw))
    sigma(:,:,:)=0.0d0
    do im=1,1 !im goes from 1 to 2
       do ip=ip_in,ip_max !ip goes from ip_in = 0 (1) for CC (NC) to ip_max = 0 (2) for CC (NC), in other words ip = 0 for CC, and ip = 1, 2 for CC
          if(iproc.ne.1) then !CC -> mY = (mp, mn), mqe = (mn, mp) Initial and final state nucleon masses?
             mY(1)=mp 
             mqe(1)=mn
             mY(2)=mn
             mqe(2)=mp
          elseif(iproc.eq.1.and.ip.eq.1) then !NC and scattering off of protons -> mY = (mp, mp), mqe = (mp, mp)
             mY(1)=mp
             mqe(1)=mp
             mY(2)=mp
             mqe(2)=mp
          elseif(iproc.eq.1.and.ip.eq.2) then !NC and scattering off of neutrons -> mY = (mn, mn), mqe = (mn, mn)
             mY(1)=mn
             mqe(1)=mn
             mY(2)=mn
             mqe(2)=mn
          endif
          
          call masses_in(iproc,ip,im) !Sets up couplings and masses of nucleons based on if we have charged current or neutral current
          do iw=1,nw !Loop over omega grid
             emu = enu - w(iw) !Energy of lepton
             pmu = sqrt( emu**2 - mlept**2 ) !Magnitude of momentum of lepton 
             q2_p=2.0d0*Enu*emu*(1.0d0-(pmu/emu)*cos(thetalept)) - mlept**2 !Q^2 = -q^2
             !if(q2_p.lt.100.0d0) cycle !Genie cuts on Q2
             !write(6,*)'Q2 = ', q2_p/1000000.0
             qval=sqrt(q2_p+ w(iw)**2) !$\sqrt{Q^{2} + w^{2}}
             res=0.0d0
             do ie=1,ne0 !Loop over energy of spectral function
                do j=1,np !Loop over momentum of spectral function
                   ep=sqrt(p(j)**2+mqe(im)**2)!On Shell energy of nucleon 
                   if(fg.eq.1) then
                      wt(iw)=w(iw)-20.0d0
                   else
                      wt(iw)=w(iw)-ep-abs(ee0(ie))+mqe(im) !$\tilde{\omega} = \omega - E_{p}^{\rm{onshell}}
                      !What is ee0? 
                   endif
                 !  if(fg.eq.1) then
                  !    costheta_p=((w(iw)+ep)**2-p(j)**2-qval**2-my(im)**2)/(2.0d0*p(j)*qval)
                  ! else
                      costheta_p=((wt(iw)+ep)**2-p(j)**2-qval**2-my(im)**2)/(2.0d0*p(j)*qval)
                   !endif
                   if(abs(costheta_p).gt.1.0d0) cycle
                   if(wt(iw).lt.0.0d0) cycle
                   pf=sqrt(p(j)**2+qval**2+2.0d0*qval*p(j)*costheta_p)
                   !if(dp(j).le.small) cycle
                   !if (pf.lt.kf) cycle
                   epf=sqrt(pf**2+mY(im)**2)
                   call sig_nuN_CCQE(mlept,qval,w(iw),wt(iw),p(j),pf,enu,thetalept,sig)  
                   if(fg.eq.1) then
                      np1=dp(j)
                   else
                      call interpolint(p0,dp0(ie,:),np0,p(j),np1,3)
                   endif
                   res=res+p(j)**2*np1*sig*epf/(p(j)*qval)*he
                enddo
             enddo
             sigma(im,ip+1,iw)=res*hp*1.e15*dble(nZ)
             write(6,*) w(iw),im,ip,sigma(im,ip+1,iw)
          enddo
       enddo
    enddo
           !......10^-15 fm^2/MeV..this is dsigma/dcosTdEmu



! write normalized responses on file
    open(unit=14,file=fname,status='unknown',form='formatted')
    !open(unit=15,file='xsec_dsigma_diff_Anu_1g_60_CC.out',status='unknown',form='formatted')

    do i=1,nw
       write(14,*) w(i),sum(sigma(1,:,i))
       !write(15,*) w(i),enu-w(i),sum(sigma(2,:,i))
    enddo
    close(14)
    !close(15)



    end program
    
