     subroutine nform(Q2,f_1p,f_1n,f_2p,f_2n,g_ep,G_en,g_mp,G_mn,g1,gp)

!      
!     revised January 2014, Virginia Tech
!
!     ichoice = 1  Dipole parametrization
!     ichoice = 2  Kelly's parametrization [PRC 70, 068282 (2004)]
!     ichoice = 3  BBBA parametrization [NPB Proc. Suppl., 159, 127 (2006)]
       
!
!     Q2 is given in units of fm^-2
!
     implicit none
     integer*4 :: ichoice,i,j,iaxial_choice
     real*8, parameter :: tcut = 0.0779191396,t0 = -0.7d0
     real*8 :: Q2,Q2G,G_ep,G_mp,G_en,G_mn,F_1s,F_1v,F_2s,F_2v,G_D &
    &       ,Lambdasq,hc,pm,mpi,zero,mu_p,mu_n,a_ep,a_mp,a_mn,tau &
    &       ,a_en,b_en,num_ep,den_ep,num_mp,den_mp,num_en,den_en &
    &       ,num_mn,den_mn,one,G_A,gan1,MAsq
     real*8 :: f_1p,f_1n,f_2p,f_2n,g1,gp,gs
     real*8 :: b_ep(3),b_mp(3),b_mn(3),aa_ep(4),aa_mp(4),aa_en(4) &
          &      ,aa_mn(4),bb_ep(4),bb_mp(4),bb_en(4),bb_mn(4)
     real*8 :: gep_coef(13),gmp_coef(13),gen_coef(13),gmn_coef(13)
     real*8 :: z
!     
     data a_ep , a_mp ,a_mn  / -0.24 , 0.12 , 2.33 / 
     data b_ep / 10.98 , 12.82 , 21.97 /
     data b_mp / 10.97 , 18.86 , 6.55 /
     data b_mn / 14.72 , 24.20 , 84.1 /
     data a_en , b_en / 1.70 , 3.30 /
!
     data aa_ep / 1.0 , -0.0578 , 0.0  , 0.0 /     
     data aa_mp / 1.0 ,  0.0150 , 0.0  , 0.0 /
     data aa_en / 0.0 ,  1.2500 , 1.30 , 0.0 /
     data aa_mn / 1.0 ,  1.8100 , 0.0  , 0.0 /
!
     data bb_ep / 11.10 ,  13.60 ,   33.00 ,   0.0 /
     data bb_mp / 11.10 ,  19.60 ,    7.54 ,   0.0 /
     data bb_en / 9.86  , 305.00 , -758.00 , 802.0 /
     data bb_mn / 14.10 ,  20.70 ,   68.70 ,   0.0 /
!     
     data gep_coef / 0.239163298067, -1.10985857441, 1.44438081306, 0.479569465603, &
          & -2.28689474187,  1.12663298498, 1.25061984354,-3.63102047159, 4.08221702379, &
          & 0.504097346499,  -5.08512046051,  3.96774254395,-0.981529071103 /
     data gmp_coef / 0.264142994136, -1.09530612212, 1.21855378178, 0.661136493537, -1.40567892503, &
          &  -1.35641843888, 1.44702915534, 4.2356697359, -5.33404565341, -2.91630052096,    &
          & 8.70740306757, -5.70699994375, 1.28081437589 /
     data gen_coef / 0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, 0.300445258602, &
          & -0.661888687179,-0.175639769687, 0.624691724461,-0.077684299367,-0.236003975259, 0.090401973470, &
          0.0, 0.0 /
     data gmn_coef / 0.257758326959,-1.079540642058, 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208, &
          & 2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0 /
!     
     hc = .197327d0
     pm = .938272d0    ! proton mass in GeV
     mpi = 0.139d0
     zero = 0.d0
     one = 1.d0
     mu_p = 2.79278d0   ! proton magnetic moment
     mu_n = -1.91315d0  ! neutron magnetic moment
     Q2G = Q2*hc*hc
     tau = Q2G/4./pm/pm
     gan1= 1.2694d0
!
!    Omar's value 
!    Lambdasq = 0.71
!    Rocco's value (I believe from Ciofi)
 !    Lambdasq = 0.69388d0
     Lambdasq = 0.840d0**2
     MAsq=1.0d0**2
!     Masq=1.350**2
!     Lambdasq = 0.55d0

     G_D=1.0d0/(1.0d0+Q2G/Lambdasq)**2
     G_A=-gan1/(1.0d0+Q2G/MAsq)**2

     ichoice=2
     iaxial_choice=0
!     
     if(ichoice.eq.1)then
!             
       G_ep = G_D 
       G_en = zero
       !G_en = -mu_n*Q2G*G_D/(1.+Q2G/pm**2)/(4.*pm**2)

       G_mp = mu_p*G_D
       G_mn = mu_n*G_D
     else if(ichoice.eq.2)then
       G_ep = (1. + a_ep*tau)/ &
    &         (1. + b_ep(1)*tau + b_ep(2)*tau**2 + b_ep(3)*tau**3) 
       G_en = G_D*a_en*tau/(1. + b_en*tau)
       G_mp = mu_p*(1. + a_mp*tau)/ &
    &              (1. + b_mp(1)*tau + b_mp(2)*tau**2 + b_mp(3)*tau**3)
       G_mn = mu_n*(1. + a_mn*tau)/ &
    &              (1. + b_mn(1)*tau + b_mn(2)*tau**2 + b_mn(3)*tau**3)
     else if(ichoice.eq.3)then   
       num_ep = zero
       den_ep = one
       num_mp = zero 
       den_mp = one
       num_en = zero 
       den_en = one
       num_mn = zero
       den_mn = one
       do i = 1,4
       j = i-1
       num_ep = num_ep + aa_ep(i)*tau**j
       den_ep = den_ep + bb_ep(i)*tau**i
       num_mp = num_mp + aa_mp(i)*tau**j
       den_mp = den_mp + bb_mp(i)*tau**i
       num_en = num_en + aa_en(i)*tau**j
       den_en = den_en + bb_en(i)*tau**i
       num_mn = num_mn + aa_mn(i)*tau**j
       den_mn = den_mn + bb_mn(i)*tau**i
       end do
       G_ep = num_ep/den_ep
       G_mp = mu_p*num_mp/den_mp
       G_en = num_en/den_en
       G_mn = mu_n*num_mn/den_mn
    else if(ichoice.eq.4) then
       G_ep=0.0d0
       G_en=0.0d0
       G_mp=0.0d0
       G_mn=0.0d0
       do i=1,13
          G_ep=G_ep+gep_coef(i)*z**(i-1)
          G_mp=G_mp+gmp_coef(i)*z**(i-1)*mu_p
          G_en=G_en+ gen_coef(i)*z**(i-1)
          G_mn=G_mn+gmn_coef(i)*z**(i-1)*mu_n
       enddo
     end if
        F_1p= 1.0d0/(1.0d0+tau)*(G_ep+tau*G_mp)
        F_1n= 1.0d0/(1.0d0+tau)*(G_en+tau*G_mn)
        F_2p=1.0d0/(1.0d0+tau)*(G_mp-G_ep)
        F_2n=1.0d0/(1.0d0+tau)*(G_mn-G_en)
        if(iaxial_choice.eq.0)then
          g1= G_A
        else
          call z_exp_axial(-Q2G,g1)
          g1 = -g1
          !write(6,*) 'FA = ', g1
        endif
        gp=2*g1*(pm**2)/(mpi**2 + Q2G)
        !gp = 0.0d0
!        gs=-0.15d0/(1.0d0+Q2G/MAsq)**2
        gs=-0.08d0/(1.0d0+Q2G/MAsq)**2
!
     return
     end   

     subroutine z_exp_axial(q2,FA)
          implicit none
          integer*4 :: i,kMax,kp4,kp3,kp2,kp1,kp0
          real*8 :: q2, FA, FA0
          real*8 :: z, fixed_a(0:2), a_err(0:2),t0_optimal,tcut_choice
          real*8 :: z0,zkp4,zkp3,zkp2,zkp1,denom

          ! Deuterium values
          !data kmax / 4 / 
          !fixed_a = (/0.759d0, -2.3d0, 0.6d0, 3.8d0, -2.3d0, -2.2d0, 0.9d0, 1.6d0, -0.82d0 /) !values from Aaron Meyer
          !a_err = (/.077d0 , .13d0, 1.0d0 , 2.5d0, 2.7d0, 3.1d0+ .1*3.1d0, 2.8d0, 2.4d0, 0.81d0 + .05*0.81d0/)!values from Aaron Meyer
          !data t0_optimal, tcut_choice / -0.28d0, 0.18d0 / !values from Aaron Meyer

          !Lattice values below 
          data kmax / 2 / !Lattice QCD values
          fixed_a = (/0.725d0, -1.63d0, 0.17d0/)
          a_err = (/.005d0, .03d0, .13d0/)
          data t0_optimal, tcut_choice /-0.5d0, 0.18d0 /


          call z_calc(t0_optimal,tcut_choice,q2,z)
          FA=0.0d0 
          !do i=0,kMax+4 ! For D2
          do i = 0, kMax ! For Lattice
               FA = FA + (fixed_a(i))*(z**i)
          enddo

          return
     end subroutine z_exp_axial

     subroutine z_calc(t0,tcut,qsquared,z_val)
          implicit none
          real*8, intent(in) :: t0,tcut,qsquared
          real*8 :: z_val

          z_val = (sqrt(tcut - qsquared) - sqrt(tcut - t0))/(sqrt(tcut - qsquared) + sqrt(tcut - t0))
          return
     end subroutine z_calc
