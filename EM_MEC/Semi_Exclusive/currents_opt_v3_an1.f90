module dirac_matrices
    implicit none
    integer*4, private, save :: i_fl, pair_isospin
    integer*4, private, save :: np_del
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
    real*8, private, parameter :: pi=acos(-1.0d0)    
    real*8, private, parameter :: fgnd=5.0d0,fpind=0.54d0
    real*8, private, parameter :: fstar=2.13d0, xmrho=775.8d0,ga=1.26d0,fpinn2=1.0094d0! 2.14/2.13 from JUAN, !=0.08*4.0d0*pi ARTURO
    real*8, private, save :: cv3,ca5,gep
    real*8, private, allocatable :: pdel(:),pot_del(:)
    real*8, private, parameter :: lpi=1300.0d0,lpind=1150.0d0
    real*8, private, save :: mqe, qval
    complex*16, private, save :: sig(3,2,2),id(2,2),id4(4,4)
    complex*16, save :: up(2),down(2)
    complex*16, private, save :: up1(2,4),up2(2,4),upp1(2,4),upp2(2,4), &
            &   ubarp1(2,4),ubarp2(2,4),ubarpp1(2,4),ubarpp2(2,4)
    complex*16, private, save :: uk1(2,4),ukp1(2,4), &
            &   ubark1(2,4),ubarkp1(2,4)
    complex*16, private, save :: iso(2,2)
    complex*16, private, save :: gamma_mu(4,4,5),g_munu(4,4)
    complex*16, private, save :: p1_sl(4,4),p2_sl(4,4),pp1_sl(4,4),pp2_sl(4,4), &
         &   k1_sl(4,4),k2_sl(4,4),q_sl(4,4), &
         &   Pi_k1(4,4),Pi_k2(4,4)
    real*8, private, save ::  p1_(4),p2_(4),pp1_(4),pp2_(4),q(4),k1(4),k2(4),l(4),lp(4)
    real*8, private, save ::  p1(4),p2(4),pp1(4),pp2(4)
    complex*16, private, save :: J_a_mu(4,4,4),J_b_mu(4,4,4),J_c_mu(4,4,4),J_d_mu(4,4,4)
    complex*16, private, save :: J_pif(4,4,4),J_sea1(4,4,4),J_sea2(4,4,4),J_pl1(4,4,4),J_pl2(4,4,4)    
    real*8, private,save :: xmd,xmn,xmpi,w, xmlept1, xmlept2
contains

subroutine dirac_matrices_in(xmd_in,xmn_in,xmpi_in, xmlept1_in, xmlept2_in)
    implicit none
    integer*4 :: i
    real*8 :: xmd_in,xmn_in,xmpi_in, xmlept1_in, xmlept2_in
    xmd=xmd_in
    xmn=xmn_in
    xmpi=xmpi_in
    xmlept1 = xmlept1_in
    xmlept2 = xmlept2_in
    sig(:,:,:)=czero
    id(:,:)=czero
    id(1,1)=cone;id(2,2)=cone
    sig(1,1,2)=cone;sig(1,2,1)=cone
    sig(2,1,2)=-ci;sig(2,2,1)=ci
    sig(3,1,1)=cone;sig(3,2,2)=-cone
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
    g_munu=czero
    g_munu(1,1)=cone;g_munu(2,2)=-cone;g_munu(3,3)=-cone;g_munu(4,4)=-cone
    up(1)=cone;up(2)=czero
    down(1)=czero;down(2)=cone

    iso(1,:)=up 
    iso(2,:)=down

end subroutine 

subroutine define_spinors()
    implicit none
    integer*4 :: i
    complex*16 :: sigp1(2,2),sigp2(2,2),sigpp1(2,2),sigpp2(2,2)
    real*8 :: cp1,cp2,cpp1,cpp2
    sigp1=czero
    sigp2=czero
    sigpp1=czero
    sigpp2=czero
    !.....initialize quadrispinors
    up1=czero
    up2=czero
    upp1=czero
    upp2=czero
!.......initialize normalization factors
    cp1=sqrt((p1(1)+xmn)/(2.0d0*p1(1)))
    cp2=sqrt((p2(1)+xmn)/(2.0d0*p2(1)))
    !Deleted these 1/sqrt(E) factors to go into the phase space integral
    cpp1=sqrt((pp1(1)+xmn)/(2.0d0*pp1(1))) 
    cpp2=sqrt((pp2(1)+xmn)/(2.0d0*pp2(1)))
!.....define sigma*p
    do i=1,3
      sigp1=sigp1+sig(i,:,:)*p1(i+1)
      sigp2=sigp2+sig(i,:,:)*p2(i+1)
      sigpp1=sigpp1+sig(i,:,:)*pp1(i+1)
      sigpp2=sigpp2+sig(i,:,:)*pp2(i+1)
    enddo
!.....build quadri-spinors    
    up1(1,1:2)=up(:)
    up1(1,3:4)=matmul(sigp1(:,:),up(:))/(p1(1)+xmn)
    up1(2,1:2)=down(:)
    up1(2,3:4)=matmul(sigp1(:,:),down(:))/(p1(1)+xmn)
    up1(:,:)=cp1*up1(:,:)
!
    up2(1,1:2)=up(:)
    up2(1,3:4)=matmul(sigp2(:,:),up(:))/(p2(1)+xmn)
    up2(2,1:2)=down(:)
    up2(2,3:4)=matmul(sigp2(:,:),down(:))/(p2(1)+xmn)
    up2(:,:)=cp2*up2(:,:)
!
    upp1(1,1:2)=up(:)
    upp1(1,3:4)=matmul(sigpp1(:,:),up(:))/(pp1(1)+xmn)
    upp1(2,1:2)=down(:)
    upp1(2,3:4)=matmul(sigpp1(:,:),down(:))/(pp1(1)+xmn)
    upp1(:,:)=cpp1*upp1(:,:)
!
    upp2(1,1:2)=up(:)
    upp2(1,3:4)=matmul(sigpp2(:,:),up(:))/(pp2(1)+xmn)
    upp2(2,1:2)=down(:)
    upp2(2,3:4)=matmul(sigpp2(:,:),down(:))/(pp2(1)+xmn)
    upp2(:,:)=cpp2*upp2(:,:)
!
    ubarp1(1,1:2)=up(:)
    ubarp1(1,3:4)=-matmul(up(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(2,1:2)=down(:)
    ubarp1(2,3:4)=-matmul(down(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(:,:)=cp1*ubarp1(:,:)
!
    ubarp2(1,1:2)=up(:)
    ubarp2(1,3:4)=-matmul(up(:),sigp2(:,:))/(p2(1)+xmn)
    ubarp2(2,1:2)=down(:)
    ubarp2(2,3:4)=-matmul(down(:),sigp2(:,:))/(p2(1)+xmn)
    ubarp2(:,:)=cp2*ubarp2(:,:)
!
    ubarpp1(1,1:2)=up(:)
    ubarpp1(1,3:4)=-matmul(up(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(2,1:2)=down(:)
    ubarpp1(2,3:4)=-matmul(down(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(:,:)=cpp1*ubarpp1(:,:)
!
    ubarpp2(1,1:2)=up(:)
    ubarpp2(1,3:4)=-matmul(up(:),sigpp2(:,:))/(pp2(1)+xmn)
    ubarpp2(2,1:2)=down(:)
    ubarpp2(2,3:4)=-matmul(down(:),sigpp2(:,:))/(pp2(1)+xmn)
    ubarpp2(:,:)=cpp2*ubarpp2(:,:)

    return
end subroutine

subroutine define_lept_spinors()
    implicit none
    integer*4 :: i
    complex*16 :: sigk1(2,2),sigkp1(2,2)
    real*8 :: ck1,ckp1
    sigk1=czero
    sigkp1=czero

    uk1=czero
    ukp1=czero

    ck1=sqrt((l(1)+xmlept1)/(2.0d0*l(1)))
    ckp1=sqrt((lp(1)+xmlept2)/(2.0d0*lp(1)))

    do i=1,3
      sigk1=sigk1+sig(i,:,:)*l(i+1)
      sigkp1=sigkp1+sig(i,:,:)*lp(i+1)
    enddo

    uk1(1,1:2)=up(:)
    uk1(1,3:4)=matmul(sigk1(:,:),up(:))/(l(1)+xmlept1)
    uk1(2,1:2)=down(:)
    uk1(2,3:4)=matmul(sigk1(:,:),down(:))/(l(1)+xmlept1)
    uk1(:,:)=ck1*uk1(:,:)

    ukp1(1,1:2)=up(:)
    ukp1(1,3:4)=matmul(sigkp1(:,:),up(:))/(lp(1)+xmlept2)
    ukp1(2,1:2)=down(:)
    ukp1(2,3:4)=matmul(sigkp1(:,:),down(:))/(lp(1)+xmlept2)
    ukp1(:,:)=ckp1*ukp1(:,:)

    ubark1(1,1:2)=up(:)
    ubark1(1,3:4)=-matmul(up(:),sigk1(:,:))/(l(1)+xmlept1)
    ubark1(2,1:2)=down(:)
    ubark1(2,3:4)=-matmul(down(:),sigk1(:,:))/(l(1)+xmlept1)
    ubark1(:,:)=ck1*ubark1(:,:)

    ubarkp1(1,1:2)=up(:)
    ubarkp1(1,3:4)=-matmul(up(:),sigkp1(:,:))/(lp(1)+xmlept2)
    ubarkp1(2,1:2)=down(:)
    ubarkp1(2,3:4)=-matmul(down(:),sigkp1(:,:))/(lp(1)+xmlept2)
    ubarkp1(:,:)=ckp1*ubarkp1(:,:)

    return
end subroutine

subroutine current_init(lepi_in,lepf_in,p1_in,p2_in,pp1_in,pp2_in,&
    &   q_in,w_in,gep_in,cv3_in,ca5_in,np_del_in,pdel_in,pot_del_in)
    implicit none
    integer*4 :: i,i_fl_in,iso_in,np_del_in
    real*8 :: p1_in(4),p2_in(4),pp1_in(4),pp2_in(4),q_in(4),k1_in(4),k2_in(4),w_in
    complex*16 :: t1_in(2),t2_in(2)
    real*8 :: lepi_in(4),lepf_in(4),gep_in,pdel_in(np_del_in),pot_del_in(np_del_in)
    real*8 :: cv3_in, ca5_in
    !Keep permanent copies of the momenta
    l=lepi_in
    lp=lepf_in
    p1_=p1_in
    p2_=p2_in
    pp1_=pp1_in
    pp2_=pp2_in
    q=q_in
    w=w_in
    gep=gep_in
    cv3=cv3_in
    ca5=ca5_in
    np_del=np_del_in

    q_sl=czero
    do i=1,4
        q_sl=q_sl+g_munu(i,i)*gamma_mu(:,:,i)*q(i)  
    enddo

    if (.not. allocated(pot_del)) then
        allocate(pot_del(np_del),pdel(np_del))
        pot_del = pot_del_in
        pdel = pdel_in
    endif

    return
end subroutine

subroutine had_current_init(p1_in,p2_in,pp1_in,pp2_in)
    implicit none
    integer*4 :: i
    real*8 :: p1_in(4),p2_in(4),pp1_in(4),pp2_in(4)
    !Defines hadronic pieces of the current in terms of the momentum passed in
    p1=p1_in
    p2=p2_in
    pp1=pp1_in
    pp2=pp2_in
    k1=pp1-p1
    k2=pp2-p2

    p1_sl=czero
    p2_sl=czero
    pp1_sl=czero
    pp2_sl=czero
    k1_sl=czero
    k2_sl=czero
       
    do i=1,4
       p1_sl=p1_sl+g_munu(i,i)*gamma_mu(:,:,i)*p1(i)  
       p2_sl=p2_sl+g_munu(i,i)*gamma_mu(:,:,i)*p2(i)
       pp1_sl=pp1_sl+g_munu(i,i)*gamma_mu(:,:,i)*pp1(i)  
       pp2_sl=pp2_sl+g_munu(i,i)*gamma_mu(:,:,i)*pp2(i)
       k1_sl=k1_sl+g_munu(i,i)*gamma_mu(:,:,i)*k1(i)  
       k2_sl=k2_sl+g_munu(i,i)*gamma_mu(:,:,i)*k2(i)
    enddo
    
    Pi_k1(:,:)=matmul(gamma_mu(:,:,5),k1_sl(:,:))/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)
    Pi_k2(:,:)=matmul(gamma_mu(:,:,5),k2_sl(:,:))/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)

end subroutine

subroutine det_Jpi()
   implicit none
   integer*4 :: mu
   real*8 :: fpik1,fpik2,frho1,frho2,fact
   fpik1=(lpi**2-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2))
   fpik2=(lpi**2-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2))
   frho1=0.0d0!1.0d0/(1.0d0-(k1(1)**2-sum(k1(2:4)**2))/xmrho**2)
   frho2=0.0d0!1.0d0/(1.0d0-(k2(1)**2-sum(k2(2:4)**2))/xmrho**2)
   !...this factor is needed to fulfill current conservation, see A3 Dekker
   fact=(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)*(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & *(1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & - 1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2)) &
        & - 1.0d0/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2)))
   do mu=1,4
      J_pif(:,:,mu)=gep*(k1(mu)-k2(mu))*Pi_k1(:,:)*fact
      J_sea1(:,:,mu)=-gep*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))-frho1/ga*gamma_mu(:,:,mu)!/fpik2**2
      J_sea2(:,:,mu)=gep*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))+frho2/ga*gamma_mu(:,:,mu)!/fpik1**2
      J_pl1(:,:,mu)=frho1/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
      J_pl2(:,:,mu)=-frho2/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
   enddo
  J_pif=J_pif*fpik1*fpik2*fpinn2/xmpi**2 
  J_sea1=J_sea1*fpik2**2*fpinn2/xmpi**2
  J_sea2=J_sea2*fpik1**2*fpinn2/xmpi**2
  J_pl1=J_pl1*fpinn2/xmpi**2*fpik2**2
  J_pl2=J_pl2*fpinn2/xmpi**2*fpik1**2
  return
end subroutine  


subroutine det_JaJb_JcJd()
    use mathtool
    implicit none
    integer*4 :: i,j,mu
    real*8 :: pa(4),pb(4),pc(4),pd(4),width,fpik1,fpik2,fpindk2,fpindk1
    real*8 :: ga,gb,gc,gd,pa2,pb2,pc2,pd2
    real*8 :: pot_pa,pot_pb,pot_pc,pot_pd,e_gs,e_bg
    complex*16 :: pa_sl(4,4),pb_sl(4,4),pc_sl(4,4),pd_sl(4,4)
    complex*16 :: xmd_a,xmd_b,xmd_c,xmd_d
    complex*16 :: j_a_1(4,4,4),j_a_2(4,4,4,4),RSa(4,4,4,4),RSb(4,4,4,4),j_b_1(4,4,4,4),j_b_2(4,4,4)
    complex*16 :: j_c_1(4,4,4),j_c_2(4,4,4,4),RSc(4,4,4,4),RSd(4,4,4,4),j_d_1(4,4,4,4),j_d_2(4,4,4)
    complex*16 :: J_a(4,4,4),J_b(4,4,4),J_c(4,4,4),J_d(4,4,4)
  
    !pa(1)=(e_gs-e_bg)*0.5d0+xmn+q(1)
    pa(:)=p1(:)+q(:)
    pb(:)=pp1(:)-q(:)
    pc(:)=p2(:)+q(:)
    pd(:)=pp2(:)-q(:)

    pa2=sum(pa(2:4)**2)
    pb2=sum(pb(2:4)**2)
    pc2=sum(pc(2:4)**2)
    pd2=sum(pd(2:4)**2)
    
    pot_pa=0.0d0
    pot_pb=0.0d0
    pot_pc=0.0d0
    pot_pd=0.0d0
    if(sqrt(pa2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pa2),pot_pa,1)
    if(sqrt(pb2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pb2),pot_pb,1)
    if(sqrt(pc2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pc2),pot_pc,1)
    if(sqrt(pd2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pd2),pot_pd,1)
    
    pa_sl=czero
    pb_sl=czero
    pc_sl=czero
    pd_sl=czero
    call delta_se(pa(1)**2-sum(pa(2:4)**2),ga,pot_pa)
    xmd_a=xmd-0.5d0*ci*ga
    call delta_se(pb(1)**2-sum(pb(2:4)**2),gb,pot_pb)
    xmd_b=xmd-0.5d0*ci*gb
    call delta_se(pc(1)**2-sum(pc(2:4)**2),gc,pot_pc)
    xmd_c=xmd-0.5d0*ci*gc
    call delta_se(pd(1)**2-sum(pd(2:4)**2),gd,pot_pd)
    xmd_d=xmd-0.5d0*ci*gd
    fpik1=(lpi**2-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2))
    fpik2=(lpi**2-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2))
    fpindk1=lpind**2/(lpind**2-k1(1)**2+sum(k1(2:4)**2))
    fpindk2=lpind**2/(lpind**2-k2(1)**2+sum(k2(2:4)**2))
    do i=1,4
       pa_sl=pa_sl+g_munu(i,i)*gamma_mu(:,:,i)*pa(i)
       pb_sl=pb_sl+g_munu(i,i)*gamma_mu(:,:,i)*pb(i)
       pc_sl=pc_sl+g_munu(i,i)*gamma_mu(:,:,i)*pc(i)
       pd_sl=pd_sl+g_munu(i,i)*gamma_mu(:,:,i)*pd(i)
    enddo

    ! costruisco i primi due termini della corrente a 2 corpi corrispondenti ai diagrammi a,b,c e d questo e' un passaggio intermedio,
    ! l'espressione finale di tali correnti e' data da j_a_mu, j_b_mu, j_c_mu, j_d_mu
    do i=1,4
      j_a_1(:,:,i)=k2(i)*id4(:,:)
      j_b_2(:,:,i)=k2(i)*id4(:,:)
      j_c_1(:,:,i)=k1(i)*id4(:,:)
      j_d_2(:,:,i)=k1(i)*id4(:,:)
      !!!...I AM USING THE FULL 
      do j=1,4
         RSa(:,:,i,j)=matmul(pa_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pa(i)*pa(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pa(j)-gamma_mu(:,:,j)*pa(i))/3.0d0/xmd) &
    !&    *(1.0d0/(pa(1)**2-sum(pa(2:4)**2)-xmd_a**2))
       &     *(pa(1)**2-sum(pa(2:4)**2)-xmd**2)/((pa(1)**2-sum(pa(2:4)**2)-xmd**2)**2+xmd**2*ga**2)

         RSb(:,:,i,j)=matmul(pb_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pb(i)*pb(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pb(j)-gamma_mu(:,:,j)*pb(i))/3.0d0/xmd) &
    !&    *(1.0d0/(pb(1)**2-sum(pb(2:4)**2)-xmd_b**2))
       &     *(pb(1)**2-sum(pb(2:4)**2)-xmd**2)/((pb(1)**2-sum(pb(2:4)**2)-xmd**2)**2+xmd**2*gb**2)
         RSc(:,:,i,j)=matmul(pc_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pc(i)*pc(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pc(j)-gamma_mu(:,:,j)*pc(i))/3.0d0/xmd) &
    !&    *(1.0d0/(pc(1)**2-sum(pc(2:4)**2)-xmd_c**2))
       &     *(pc(1)**2-sum(pc(2:4)**2)-xmd**2)/((pc(1)**2-sum(pc(2:4)**2)-xmd**2)**2+xmd**2*gc**2)
         RSd(:,:,i,j)=matmul(pd_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pd(i)*pd(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pd(j)-gamma_mu(:,:,j)*pd(i))/3.0d0/xmd) &
    !&    *(1.0d0/(pd(1)**2-sum(pd(2:4)**2)-xmd_d**2))
       &     *(pd(1)**2-sum(pd(2:4)**2)-xmd**2)/((pd(1)**2-sum(pd(2:4)**2)-xmd**2)**2+xmd**2*gd**2)
         J_a_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
         J_b_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
         J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
         J_d_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
    !         J_a_2(:,:,i,j)=0.5d0*cv3*matmul(q(j)*gamma_mu(:,:,i)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),gamma_mu(:,:,i))), &
    !    &                   gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
    !         J_b_1(:,:,i,j)=0.5d0*cv3*matmul(gamma_mu(:,:,5),q(i)*gamma_mu(:,:,j)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),& 
    !    &                   gamma_mu(:,:,i))))+ca5*xmn*g_munu(j,i)*id4(:,:)
    !         J_c_2(:,:,i,j)=0.5d0*cv3*matmul(q(j)*gamma_mu(:,:,i)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),gamma_mu(:,:,i))), &
    !    &                   gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
    !         J_d_1(:,:,i,j)=0.5d0*cv3*matmul(gamma_mu(:,:,5),q(i)*gamma_mu(:,:,j)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),& 
    !    &                   gamma_mu(:,:,i))))+ca5*xmn*g_munu(j,i)*id4(:,:)


      enddo
    enddo
    ! costruisco Jmua, Jmub
   do mu=1,4
      J_a(:,:,mu)=czero
      J_b(:,:,mu)=czero
      J_c(:,:,mu)=czero
      J_d(:,:,mu)=czero
      do i=1,4
         do j=1,4
            J_a(:,:,mu)=J_a(:,:,mu)+matmul(J_a_1(:,:,i)*g_munu(i,i),matmul(RSa(:,:,i,j),g_munu(j,j)*J_a_2(:,:,j,mu)))
            J_b(:,:,mu)=J_b(:,:,mu)+matmul(J_b_1(:,:,mu,i)*g_munu(i,i),matmul(RSb(:,:,i,j),g_munu(j,j)*J_b_2(:,:,j))) 
            J_c(:,:,mu)=J_c(:,:,mu)+matmul(J_c_1(:,:,i)*g_munu(i,i),matmul(RSc(:,:,i,j),g_munu(j,j)*J_c_2(:,:,j,mu)))
            J_d(:,:,mu)=J_d(:,:,mu)+matmul(J_d_1(:,:,mu,i)*g_munu(i,i),matmul(RSd(:,:,i,j),g_munu(j,j)*J_d_2(:,:,j))) 
         enddo
      enddo
   enddo

    J_a_mu=J_a*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
    J_b_mu=J_b*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
    J_c_mu=J_c*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
    J_d_mu=J_d*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
 
end subroutine

subroutine JDeltaFixed(jtot)
    implicit none
    integer*4 :: i1,i2,f1,f2,i,j,ti1,ti2,tf1,tf2
    complex*16 :: j_1(2,2),j_2(2,2)
    complex*16 :: ja_sub(2,2,4),jb_sub(2,2,4),jc_sub(2,2,4),jd_sub(2,2,4)
    complex*16 :: ja(2,2,2,2,4), jb(2,2,2,2,4)
    complex*16 :: jc(2,2,2,2,4), jd(2,2,2,2,4), jtot(2,2,2,2,2,2,2,2,4)
    complex*16 :: iso_a(2,2,2,2), iso_b(2,2,2,2), iso_c(2,2,2,2), iso_d(2,2,2,2)

    iso_a = czero
    iso_b = czero
    iso_d = czero
    iso_c = czero

    ja = czero
    jb = czero
    jc = czero
    jd = czero
    ja_sub=czero
    jb_sub=czero
    jc_sub=czero
    jd_sub=czero
    jtot=czero

    !Fill spinors for nucleons with given momenta (specified outside this function)
    call define_spinors()
    !Compute pion and delta current matrices
    call det_Jpi()
    call det_JaJb_JcJd()

    do i1=1,2
      do f1=1,2
        j_2(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k2(:,:),up2(i1,:)))
        j_1(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k1(:,:),up1(i1,:)))
        do i=1,4
            ja_sub(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
            jb_sub(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_b_mu(:,:,i),up1(i1,:)))
            jc_sub(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_c_mu(:,:,i),up2(i1,:)))
            jd_sub(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_d_mu(:,:,i),up2(i1,:)))
        enddo
      enddo
    enddo

    do i1=1,2
        do i2=1,2
            do f1=1,2
                do f2=1,2
                    do i=1,4
                        ja(f2,f1,i2,i1,i)=j_2(f2,i2)*ja_sub(f1,i1,i)
                        jb(f2,f1,i2,i1,i)=j_2(f2,i2)*jb_sub(f1,i1,i)
                        jc(f2,f1,i2,i1,i)=j_1(f1,i1)*jc_sub(f2,i2,i)
                        jd(f2,f1,i2,i1,i)=j_1(f1,i1)*jd_sub(f2,i2,i)
                    enddo
                enddo
            enddo
        enddo
    enddo

    do ti1=1,2
        do ti2=1,2
            do tf1=1,2
                do tf2=1,2
                    iso_a(tf2,tf1,ti2,ti1)=IDeltaA(iso(ti1,:),iso(ti2,:),iso(tf1,:),iso(tf2,:))
                    jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) = ja(:,:,:,:,:)*iso_a(tf2,tf1,ti2,ti1)
                    iso_b(tf2,tf1,ti2,ti1)=IDeltaB(iso(ti1,:),iso(ti2,:),iso(tf1,:),iso(tf2,:))
                    jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) = jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) + jb(:,:,:,:,:)*iso_b(tf2,tf1,ti2,ti1)
                    iso_c(tf2,tf1,ti2,ti1)=IDeltaC(iso(ti1,:),iso(ti2,:),iso(tf1,:),iso(tf2,:))
                    jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) = jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) + jc(:,:,:,:,:)*iso_c(tf2,tf1,ti2,ti1)
                    iso_d(tf2,tf1,ti2,ti1)=IDeltaD(iso(ti1,:),iso(ti2,:),iso(tf1,:),iso(tf2,:))
                    jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) = jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) + jd(:,:,:,:,:)*iso_d(tf2,tf1,ti2,ti1)
                enddo
            enddo
        enddo
    enddo

    return

end subroutine JDeltaFixed

subroutine JDelta(janti)
    implicit none
    integer*4 :: i1,i2,f1,f2,i,j,ti1,ti2,tf1,tf2
    complex*16 :: j1212(2,2,2,2,2,2,2,2,4), j1221(2,2,2,2,2,2,2,2,4)
    complex*16 :: j2121(2,2,2,2,2,2,2,2,4), j2112(2,2,2,2,2,2,2,2,4)
    complex*16 :: janti(2,2,2,2,2,2,2,2,4)

    janti=czero
    j1212=czero
    j1221=czero
    j2112=czero
    j2121=czero

    call had_current_init(p1_,p2_,pp1_,pp2_)
    call JDeltaFixed(j1212)

    !call had_current_init(p2_,p1_,pp1_,pp2_)
    !call JDeltaFixed(j2112)

    call had_current_init(p1_,p2_,pp2_,pp1_)
    call JDeltaFixed(j1221)

    !call had_current_init(p2_,p1_,pp2_,pp1_)
    !call JDeltaFixed(j2121)

    do ti1=1,2
        do ti2=1,2
            do tf1=1,2
                do tf2=1,2
                    do i1=1,2
                        do i2=1,2
                            do f1=1,2
                                do f2=1,2
                                    janti(f2,f1,i2,i1,tf2,tf1,ti2,ti1,:) = j1212(f2,f1,i2,i1,tf2,tf1,ti2,ti1,:) & 
                                    &  - j1221(f1,f2,i2,i1,tf1,tf2,ti2,ti1,:) & 
                                    &  - j2112(f2,f1,i1,i2,tf2,tf1,ti1,ti2,:) & 
                                    &  + j2121(f1,f2,i1,i2,tf1,tf2,ti1,ti2,:) 
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine JDelta

subroutine JPiFixed(jtot)
   implicit none
    integer*4 :: i1,i2,f1,f2,i,j,ti1,ti2,tf1,tf2
    complex*16 :: j_1(2,2),j_2(2,2)
    complex*16 :: js1_sub(2,2,4),js2_sub(2,2,4)
    complex*16 :: js1(2,2,2,2,4), js2(2,2,2,2,4)
    complex*16 :: jf(2,2,2,2,4), jtot(2,2,2,2,2,2,2,2,4)
    complex*16 :: iso_a(2,2,2,2)


    !Fill spinors for nucleons with given momenta (specified outside this function)
    call define_spinors()
    !Compute pion current matrices
    call det_Jpi()

    iso_a=czero

    js1 = czero
    js2= czero
    jf = czero
    js1_sub=czero
    js2_sub=czero


   do i1=1,2
      do f1=1,2
         J_2(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k2(:,:),up2(i1,:)))
         J_1(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k1(:,:),up1(i1,:)))
         do i=1,4
            js1_sub(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_sea1(:,:,i),up1(i1,:)))
            js2_sub(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_sea2(:,:,i),up2(i1,:)))
            do i2=1,2
                do f2=1,2
                    jf(f2,f1,i2,i1,i)=sum(ubarpp1(f1,:)*matmul(J_pif(:,:,i),up1(i1,:))) &
                        & *sum(ubarpp2(f2,:)*matmul(Pi_k2(:,:),up2(i2,:)))
                enddo
            enddo
        enddo
      enddo
   enddo   

    do i1=1,2
        do i2=1,2
            do f1=1,2
                do f2=1,2
                    iso_a(f2,f1,i2,i1)=-Ivz(iso(i1,:),iso(i2,:),iso(f1,:),iso(f2,:))
                    do i=1,4
                        js1(f2,f1,i2,i1,i)=j_2(f2,i2)*js1_sub(f1,i1,i)
                        js2(f2,f1,i2,i1,i)=j_1(f1,i1)*js2_sub(f2,i2,i)   
                    enddo
                enddo
            enddo
        enddo
    enddo

    do ti1=1,2
        do ti2=1,2
            do tf1=1,2
                do tf2=1,2
                    jtot(:,:,:,:,tf2,tf1,ti2,ti1,:) = iso_a(tf2,tf1,ti2,ti1)&
                        & * (js1(:,:,:,:,:) + js2(:,:,:,:,:) + jf(:,:,:,:,:))
                enddo
            enddo
        enddo
    enddo


   return
end subroutine JPiFixed

subroutine JPi(janti)
    implicit none
    integer*4 :: i1,i2,f1,f2,i,j,ti1,ti2,tf1,tf2
    complex*16 :: j1212(2,2,2,2,2,2,2,2,4), j1221(2,2,2,2,2,2,2,2,4)
    complex*16 :: j2121(2,2,2,2,2,2,2,2,4), j2112(2,2,2,2,2,2,2,2,4)
    complex*16 :: janti(2,2,2,2,2,2,2,2,4)

    janti=czero
    j1212=czero
    j1221=czero
    j2112=czero
    j2121=czero

    call had_current_init(p1_,p2_,pp1_,pp2_)
    call JPiFixed(j1212)

    !call had_current_init(p2_,p1_,pp1_,pp2_)
    !call JPiFixed(j2112)

    call had_current_init(p1_,p2_,pp2_,pp1_)
    call JPiFixed(j1221)

    !call had_current_init(p2_,p1_,pp2_,pp1_)
    !call JPiFixed(j2121)


    do ti1=1,2
        do ti2=1,2
            do tf1=1,2
                do tf2=1,2
                    do i1=1,2
                        do i2=1,2
                            do f1=1,2
                                do f2=1,2
                                    janti(f2,f1,i2,i1,tf2,tf1,ti2,ti1,:) = j1212(f2,f1,i2,i1,tf2,tf1,ti2,ti1,:) & 
                                    &  - j1221(f1,f2,i2,i1,tf1,tf2,ti2,ti1,:) & 
                                    &  - j2112(f2,f1,i1,i2,tf2,tf1,ti1,ti2,:) & 
                                    &  + j2121(f1,f2,i1,i2,tf1,tf2,ti1,ti2,:) 
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine JPi

 
subroutine lept_tens(lept)
   implicit none
   integer*4 :: i1,f1,i,j
   complex*16 :: J_mu(2,2,4),J_mu_dag(2,2,4)
   complex*16 :: lept(4,4)

   do i1=1,2
      do f1=1,2
         do i=1,4
            J_mu(f1,i1,i)=sum(ubarkp1(f1,:)*matmul(gamma_mu(:,:,i),uk1(i1,:)))/sqrt(2.0d0)
            J_mu_dag(f1,i1,i)=conjg(J_mu(f1,i1,i))
         enddo
      enddo
   enddo
   
   lept=czero
   do i1=1,2
      do f1=1,2
         do i=1,4
            do j=1,4
               lept(i,j)=lept(i,j)+J_mu_dag(f1,i1,i)*J_mu(f1,i1,j)
            enddo   
         enddo
      enddo
   enddo

   lept = lept*4.0d0
 
  return
end subroutine lept_tens

subroutine contract(hadr,lept,sig)
   implicit none
   integer*4 :: i,j
   complex*16 ::lept(4,4),hadr(4,4),sig

   sig=0.0d0
    do i=1,4
        do j=1,4
           sig = sig + g_munu(i,i)*lept(i,j)*hadr(i,j)*g_munu(j,j)
        enddo 
    enddo
    return
end subroutine contract


subroutine delta_se(pd2,width,pot)
   implicit none
   real*8 :: pd2,width,kpi,ekpi,eknuc,r2a,rfa,pot
   width=0.0d0
   !pot=0.0d0!-40.0d0
   if (pd2.ge.(xmpi+xmn)**2)then
      kpi=dsqrt(1.0d0/4.0d0/pd2*(pd2-(xmn+xmpi)**2)*(pd2-(xmn-xmpi)**2))
      ekpi=dsqrt(kpi**2+xmpi**2)
      eknuc=dsqrt(kpi**2+xmn**2)
      r2a= (eknuc-ekpi)**2 -4.0d0*kpi**2
      rfa=sqrt(0.95d0*xmn**2/(0.95d0*xmn**2-r2a))
!.....definizione Gamma
!      width=(4.0d0*fpind)**2/12.0d0/pi/xmpi**2*kpi**3/sqrt(pd2)*(xmn+eknuc)  &
!   &  *rfa**2!*(lpind**2/(lpind**2-xmpi**2))**2
!      width=120.0d0

      width=0.38/(3.0d0*xmpi**2)*kpi**3/sqrt(pd2)*(xmn+eknuc)  &
   &  *rfa**2!-pot*2.0d0!*(lpind**2/(lpind**2-xmpi**2))**2


      return
   endif
!.....Eq.3.26 Phys.Rev. C 49 2650

   return
end subroutine

function Ivz(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: Ivz

    Ivz = ci*(me(1,it1,itp1)*me(2,it2,itp2) - me(2,it1,itp1)*me(1,it2,itp2))

    return
end function Ivz


function IDeltaA(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaA, c

    c = me(3,it2,itp2)*iden(it1,itp1)

    IDeltaA = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.)
    !IDeltaA = me(1,it1,itp1)*me(1,it2,itp2)

    return
end function IDeltaA

function IDeltaADag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaADag, c

    c = me(3,it2,itp2)*iden(it1,itp1)

    IDeltaADag = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.)
    !IDeltaADag = me(1,it1,itp1)*me(1,it2,itp2)


    return
end function IDeltaADag

function IDeltaB(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaB, c

    c = me(3,it2,itp2)*iden(it1,itp1) 

    IDeltaB = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaB

function IDeltaBDag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaBDag, c

    c = me(3,it2,itp2)*iden(it1,itp1)

    IDeltaBDag = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaBDag

function IDeltaC(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaC, c

    c = me(3,it1,itp1)*iden(it2,itp2) 

    IDeltaC = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaC

function IDeltaCDag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaCDag, c

    c = me(3,it1,itp1)*iden(it2,itp2)  

    IDeltaCDag = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaCDag

function IDeltaD(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaD, c

    c = me(3,it1,itp1)*iden(it2,itp2)  

    IDeltaD = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaD

function IDeltaDDag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaDDag, c

    c = me(3,it1,itp1)*iden(it2,itp2)  

    IDeltaDDag = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaDDag

function me(i,it,itp)
    implicit none
    integer*4 :: i
    complex*16 :: me, it(2),itp(2), matrix(2)

    me = sum(itp(:)*matmul(sig(i,:,:),it))
    return
end function me


function iden(it,itp)
    implicit none 
    complex*16:: iden, it(2),itp(2)

    iden = sum(itp(:)*matmul(id,it))
    return
end function iden

end module

