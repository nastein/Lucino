!This program constructs the nuclear current operator $j^{\mu}$ (see Ab inition caluclations eq. 17)
!$j^{\mu}_{EM} = \bar{u}(p^{'})[F_{1}\gamma^{\mu} + i\sigma^{\mu\nu}q_{\nu}F_{2}/2m_{N}]u(p)

!Useful tools that we will use to construct dirac matrices $\gamma^{\mu}$ and $\sigma^{\mu\nu}$
module dirac_matrices
    implicit none
    integer*4, private, save :: i_fl
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
    real*8, private, parameter :: pi=acos(-1.0d0)  
    real*8, private, parameter :: mpi = 134.0d0
    real*8, private, save :: mqe, qval

    complex*16, private, save :: sig(3,2,2),id(2,2),id4(4,4),up(2),down(2) !pauli matrices, identity, and 2d basis vectors
    complex*16, private, save :: up1(2,4),upp1(2,4), &
            &   ubarp1(2,4),ubarpp1(2,4) !spinors $u(p)/u(p')/\bar{u}(p)/\bar{u}(p')$
    complex*16, private, save :: gamma_mu(4,4,5),g_munu(4,4),sigma_munu(4,4,4,4)!$\gamma$ matrices, $g^{\mu\nu},\sigma^{\mu\nu}$ 
    complex*16, private, save :: q_sl(4,4)
    real*8, private, save ::  p1(4),pp1(4),q(4),qt(4)!p,p',q, initial, final and four momentum transfer 
    complex*16, private, save :: J_1(4,4,4)!nuclear current operator 
    real*8, private,save :: xmn
contains

!Now we actually construct the matrices 
subroutine dirac_matrices_in(xmn_in)
    implicit none
    integer*4 :: i,j
    real*8 :: xmd_in,xmn_in,xmpi_in
    xmn=xmn_in
    sig(:,:,:)=czero
    id(:,:)=czero
    id(1,1)=cone;id(2,2)=cone!2x2 identity
    !construct the pauli matrices here
    sig(1,1,2)=cone;sig(1,2,1)=cone
    sig(2,1,2)=-ci;sig(2,2,1)=ci
    sig(3,1,1)=cone;sig(3,2,2)=-cone
    !Set up gamma matrices 
    gamma_mu=czero
    gamma_mu(1:2,1:2,1)=id;gamma_mu(3:4,3:4,1)=-id
    id4=czero    
    id4(1:2,1:2)=id;id4(3:4,3:4)=id
    do i=2,4
      gamma_mu(1:2,3:4,i)=sig(i-1,:,:)
      gamma_mu(3:4,1:2,i)=-sig(i-1,:,:)
    enddo
    gamma_mu(1:2,3:4,5)=id
    gamma_mu(3:4,1:2,5)=id

    !Set up $g^{\mu\nu}
    g_munu=czero
    g_munu(1,1)=cone;g_munu(2,2)=-cone;g_munu(3,3)=-cone;g_munu(4,4)=-cone
    do i=1,4
       do j=1,4
          !$\sigma^{\mu\nu} = \frac{i}{2}[\gamma^{\mu},\gamma^{\nu}]
          sigma_munu(:,:,i,j)=ci*0.5d0*(matmul(gamma_mu(:,:,i),gamma_mu(:,:,j)) &
               &     -matmul(gamma_mu(:,:,j),gamma_mu(:,:,i)))
       enddo
    enddo
    !2d basis vectors
    up(1)=cone;up(2)=czero
    down(1)=czero;down(2)=cone

  end subroutine  

!Now we will define the spinors for the initial and final nucleons
subroutine define_spinors()
    implicit none
    integer*4 :: i
    complex*16 :: sigp1(2,2),sigp2(2,2),sigpp1(2,2),sigpp2(2,2)
    real*8 :: cp1,cp2,cpp1,cpp2
    sigp1=czero
    sigpp1=czero
    !.....initialize quadrispinors
    up1=czero
    upp1=czero
!.......initialize normalization factors
    !Normalization factors: 
    !$cp_{1} = \sqrt{\frac{(p(1) + xmn)}{2p(1)}}
    !$cp_{1}^{'} = \sqrt{\frac{(p^{'}(1) + xmn)}{2p^{'}(1)}}
    cp1=sqrt((p1(1)+xmn)/(2.0d0*p1(1)))
    cpp1=sqrt((pp1(1)+xmn)/(2.0d0*pp1(1)))
!.....define sigma*p
    do i=1,3
      sigp1=sigp1+sig(i,:,:)*p1(i+1)
      sigpp1=sigpp1+sig(i,:,:)*pp1(i+1)
    enddo
!.....build quadri-spinors    
    up1(1,1:2)=up(:)
    up1(1,3:4)=matmul(sigp1(:,:),up(:))/(p1(1)+xmn)
    up1(2,1:2)=down(:)
    up1(2,3:4)=matmul(sigp1(:,:),down(:))/(p1(1)+xmn)
    up1(:,:)=cp1*up1(:,:)
!
    upp1(1,1:2)=up(:)
    upp1(1,3:4)=matmul(sigpp1(:,:),up(:))/(pp1(1)+xmn)
    upp1(2,1:2)=down(:)
    upp1(2,3:4)=matmul(sigpp1(:,:),down(:))/(pp1(1)+xmn)
    upp1(:,:)=cpp1*upp1(:,:)
!
    ubarp1(1,1:2)=up(:)
    ubarp1(1,3:4)=-matmul(up(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(2,1:2)=down(:)
    ubarp1(2,3:4)=-matmul(down(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(:,:)=cp1*ubarp1(:,:)
!
    ubarpp1(1,1:2)=up(:)
    ubarpp1(1,3:4)=-matmul(up(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(2,1:2)=down(:)
    ubarpp1(2,3:4)=-matmul(down(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(:,:)=cpp1*ubarpp1(:,:)
    return
end subroutine

subroutine current_init(p1_in,pp1_in,q_in,qt_in)
    implicit none
    real*8 ::  p1_in(4),pp1_in(4),q_in(4),qt_in(4)
    p1=p1_in
    pp1=pp1_in
    q=q_in
    qt=qt_in
    return
end subroutine

!Now we will actually construct the nuclear current operator 
subroutine det_Ja(f1v,f2v,ffa,ffp)
  implicit none
    integer*4 :: mu,nu
  real*8 :: f1v,f2v,ffa,ffp 
  complex*16 :: J_1_V(4,4,4),J_1_A(4,4,4) 
  J_1_V=czero
  J_1_A=czero

  do mu=1,4
     J_1_V(:,:,mu)=czero
     J_1_A(:,:,mu)=czero
     do nu=1,4
        J_1_V(:,:,mu)=J_1_V(:,:,mu)+ci*f2v*sigma_munu(:,:,mu,nu)&
             &     *g_munu(nu,nu)*qt(nu)/2.0d0/xmn
     enddo
     J_1_V(:,:,mu)=J_1_V(:,:,mu)+f1v*gamma_mu(:,:,mu)

     !! Now add axial pieces of current 
     J_1_A(:,:,mu)=J_1_A(:,:,mu)+ffa*matmul(gamma_mu(:,:,mu),gamma_mu(:,:,5))

     J_1_A(:,:,mu)=J_1_A(:,:,mu)+ffp*gamma_mu(:,:,5)*qt(mu)/xmn

  enddo

  !! Applying current conservation first to the vector current before adding axial current
  !J_1_V(:,:,4) = (q(1)/q(4))*J_1_V(:,:,1)

  J_1 = J_1_V + J_1_A
  
end subroutine det_Ja

!Now we compute the response functions as in eq. 12
!$Rcc = R00
!$Rcl = -1/2 (R03 + R30)
!$Rll = R33
!$Rt = R11 + R22
!$Rtp = -i/2 (R12 - R21)
subroutine det_res1b(resp)
   implicit none
   integer*4 :: i1,f1,i,j
   complex*16 :: J_mu(2,2,4),J_mu_dag(2,2,4)
   complex*16 :: res(4,4)
   real*8 :: resp(5)
   

   do i1=1,2
      do f1=1,2
         do i=1,4
            J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
            J_mu_dag(f1,i1,i)=conjg(J_mu(f1,i1,i))
         enddo
      enddo
   enddo

   
   res=czero
   do i1=1,2
      do f1=1,2
         do i=1,4
            do j=1,4
              res(i,j)=res(i,j)+J_mu_dag(f1,i1,i)*J_mu(f1,i1,j)
            enddo
         enddo
      enddo
   enddo
   !write(6,*)"res(,",i,",",j,")="res(i,j)
   
   resp(1)=res(1,1)
   resp(2)=-0.5d0*(res(1,4) + res(4,1))
   resp(3)=res(4,4)
   resp(4)=res(2,2)+res(3,3)
   resp(5)=-ci*0.5d0*(res(2,3) - res(3,2))
 
  return
end subroutine det_res1b
end module dirac_matrices





