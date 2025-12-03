module isospin_op
    implicit none
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
    complex*16, private, save :: sig(3,2,2),id(2,2)
    complex*16, private, save :: up(2),down(2)
    logical, private, save :: CC
contains

subroutine set_up_isospin_ops(CC_in)
    use mympi
    implicit none
    logical :: CC_in

    sig(:,:,:)=czero
    id(:,:)=czero
    id(1,1)=cone;id(2,2)=cone
    sig(1,1,2)=cone;sig(1,2,1)=cone
    sig(2,1,2)=-ci;sig(2,2,1)=ci
    sig(3,1,1)=cone;sig(3,2,2)=-cone
    up(1)=cone;up(2)=czero
    down(1)=czero;down(2)=cone
    if(CC_in.eqv..false.) then
        CC = .false.
        if (myrank().eq.0) then
            write(6,*)'Setting up isospin operators for EM Current'
        endif
    else
        CC = .true.
        if (myrank().eq.0) then
            write(6,*)'Setting up isospin operators for Charged Current'
        endif
    endif


end subroutine set_up_isospin_ops

function Iv(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: Iv

    if(CC.eqv..true.) then
        Iv = Ivplus(it1,it2,itp1,itp2)
    else
        Iv = Ivz(it1,it2,itp1,itp2)
    endif
end function Iv 

function IDeltaA(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaA

    if(CC.eqv..true.) then
        IDeltaA = IDeltaA_EW(it1,it2,itp1,itp2)
    else
        IDeltaA = IDeltaA_EM(it1,it2,itp1,itp2)
    endif
end function IDeltaA

function IDeltaB(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaB

    if(CC.eqv..true.) then
        IDeltaB = IDeltaB_EW(it1,it2,itp1,itp2)
    else
        IDeltaB = IDeltaB_EM(it1,it2,itp1,itp2)
    endif
end function IDeltaB

function IDeltaC(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaC

    if(CC.eqv..true.) then
        IDeltaC = IDeltaC_EW(it1,it2,itp1,itp2)
    else
        IDeltaC = IDeltaC_EM(it1,it2,itp1,itp2)
    endif
end function IDeltaC

function IDeltaD(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaD

    if(CC.eqv..true.) then
        IDeltaD = IDeltaD_EW(it1,it2,itp1,itp2)
    else
        IDeltaD = IDeltaD_EM(it1,it2,itp1,itp2)
    endif
end function IDeltaD

function Ivz(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: Ivz

    Ivz = ci*(me(1,it1,itp1)*me(2,it2,itp2) - me(2,it1,itp1)*me(1,it2,itp2))

    return
end function Ivz

function IDeltaA_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaA_EM, c

    c = me(3,it2,itp2)*iden(it1,itp1)

    IDeltaA_EM = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaA_EM

function IDeltaADag_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaADag_EM, c

    c = me(3,it2,itp2)*iden(it1,itp1)

    IDeltaADag_EM = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaADag_EM

function IDeltaB_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaB_EM, c

    c = me(3,it2,itp2)*iden(it1,itp1) 

    IDeltaB_EM = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaB_EM

function IDeltaBDag_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaBDag_EM, c

    c = me(3,it2,itp2)*iden(it1,itp1)

    IDeltaBDag_EM = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaBDag_EM

function IDeltaC_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaC_EM, c

    c = me(3,it1,itp1)*iden(it2,itp2) 

    IDeltaC_EM = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaC_EM

function IDeltaCDag_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaCDag_EM, c

    c = me(3,it1,itp1)*iden(it2,itp2)  

    IDeltaCDag_EM = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaCDag_EM

function IDeltaD_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaD_EM, c

    c = me(3,it1,itp1)*iden(it2,itp2)  

    IDeltaD_EM = (2.*c/3.) - (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaD_EM

function IDeltaDDag_EM(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaDDag_EM, c

    c = me(3,it1,itp1)*iden(it2,itp2)  

    IDeltaDDag_EM = (2.*c/3.) + (Ivz(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaDDag_EM

function Ivminus(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: Ivminus

    Ivminus = ci*( ( me(2,it1,itp1)*me(3,it2,itp2) - me(3,it1,itp1)*me(2,it2,itp2) ) & 
    &    - ci*( me(3,it1,itp1)*me(1,it2,itp2) - me(1,it1,itp1)*me(3,it2,itp2) ) )

    return
end function Ivminus

function Ivplus(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: Ivplus

    Ivplus = ci*( ( me(2,it1,itp1)*me(3,it2,itp2) - me(3,it1,itp1)*me(2,it2,itp2) ) & 
    &    + ci*( me(3,it1,itp1)*me(1,it2,itp2) - me(1,it1,itp1)*me(3,it2,itp2) ) )

    return
end function Ivplus


function IDeltaA_EW(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaA_EW, c

    c = (me(1,it2,itp2) + ci*me(2,it2,itp2))*iden(it1,itp1)

    IDeltaA_EW = (2.*c/3.) - (Ivplus(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaA_EW

function IDeltaB_EW(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaB_EW, c

    c = (me(1,it2,itp2) + ci*me(2,it2,itp2))*iden(it1,itp1)

    IDeltaB_EW = (2.*c/3.) + (Ivplus(it1,it2,itp1,itp2)/3.)  

    return
end function IDeltaB_EW

function IDeltaC_EW(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaC_EW, c

    c = (me(1,it1,itp1) + ci*me(2,it1,itp1))*iden(it2,itp2)

    IDeltaC_EW = (2.*c/3.) + (Ivplus(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaC_EW

function IDeltaD_EW(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaD_EW, c

    c = (me(1,it1,itp1) + ci*me(2,it1,itp1))*iden(it2,itp2)

    IDeltaD_EW = (2.*c/3.) - (Ivplus(it1,it2,itp1,itp2)/3.)  

    return
end function IDeltaD_EW

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
