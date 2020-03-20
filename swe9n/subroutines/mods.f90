module basicVars
implicit none

  integer,parameter::C_K1=4,C_K2=8    
  integer::UTMZone=44
  integer(kind=C_K1)::i,i1,i2,j,j1,j2,k,k1,k2,l,l1,l2
  integer(kind=C_K1),parameter::tf=10,ifl(10)=(/ (i,i=11,20)  /)
  integer(kind=C_K1)::iel,n1,n2,n3,n4
  integer(kind=C_K1)::na(9)
  integer(kind=C_K1)::tmpi1,tmpi2,tmpi3,tmpi4,tmpi5
  integer(kind=C_K1)::sysClk(10)    
  
  
  real(kind=C_K2)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5 
  real(kind=C_K2)::sysRate
  real(kind=C_K2)::rbfAc,rbfDc


  character(len=100)::probname,text,resumeFile
  character(len=100)::tmps(4) 


  logical::ex


  !!----------------------Constants-----------------------!!
  !! Neighbour nodes and elements limit
  integer(kind=C_K1),parameter::maxNePoi=60
  !! Gravity
  real(kind=C_K2),parameter::grav=9.81d0          
  !! PI
  real(kind=C_K2),parameter::pi=atan(1d0)*4d0
  real(kind=C_K2),parameter::deg2rad=pi/180d0
  !! Water density
  real(kind=C_K2),parameter::rhoW=1000d0
  !! Bottom Friction coefficient
  real(kind=C_K2),parameter::botCd=0.015d0
  !! Smagorinsky coefficient
  real(kind=C_K2),parameter::smCsSq=0.2d0**2
  !! RPIM size
  integer(kind=C_K1),parameter::rbfOr=6,rbfMaxNe=120
  !!--------------------End Constants---------------------!!


  !! File definitions
  ! tf=10                      !Terminal Output File
  ! ifl(1)=11                  !Input file
  ! ifl(2)=12                  !Mesh file
  ! ifl(3)=13                  !Output file
  ! ifl(4)=14                  !Cyclone track file
  ! ifl(5)=15                  !Probe Output File


end module basicVars


module jacobianModule
use basicVars
implicit none
  
  type jacbType
    real(kind=C_K2)::D(9),D11(9),D12(9),D21(9),D22(9)    
  end type jacbType

end module jacobianModule


module rpimModule
use basicVars
implicit none  

  type, public :: rpimType
    integer(kind=C_K1) :: Rid, iv, rbfN
    integer(kind=C_K1) :: jv(rbfMaxNe)
    real(kind=C_K2) :: Rpx,Rpy
    real(kind=C_K2) :: phi(rbfMaxNe)
  contains
    procedure :: init => rpimInit    
    procedure :: calPhi
    !procedure :: calcSmooth

  end type rpimType

contains

  subroutine rpimInit(this,Rid,Rpx,Rpy,iv,jv)
  use basicVars
  implicit none

    class(rpimType),intent(inout)::this
    integer(kind=C_K1),intent(in)::Rid,iv,jv(iv)
    real(kind=C_K2),intent(in)::Rpx,Rpy

    this%Rid=Rid
    this%Rpx=Rpx
    this%Rpy=Rpy
    this%iv=iv
    this%rbfN=iv+rbfOr    
    this%jv=0
    this%jv(1:iv)=jv
    this%phi=0d0
  end subroutine rpimInit


  subroutine calPhi(this,cx,cy)
  use basicVars
  implicit none

    class(rpimType),intent(inout)::this    
    integer(kind=C_K1)::errVar1
    real(kind=C_K2),intent(in)::cx(this%iv),cy(this%iv)
    real(kind=C_K2)::G(this%rbfN,this%rbfN)
    real(kind=C_K2)::invG(this%rbfN,this%rbfN)
    real(kind=C_K2)::eye(this%rbfN,this%rbfN)
    real(kind=C_K2)::Rm(this%rbfN)
    real(kind=C_K2)::xi,xj,xk,yi,yj,yk,rji,rkj
    real(kind=C_K2)::poly(rbfOr)

    xi=this%Rpx
    yi=this%Rpy

    G=0d0
    do j=1,this%iv
      xj=cx(j)
      yj=cy(j)
      rji=dsqrt((xj-xi)**2 + (yj-yi)**2)/rbfDc
      !rji=((xj-xi)**2 + (yj-yi)**2)/rbfDc**2

      do k=1,this%iv
        xk=cx(k)
        yk=cy(k)
        rkj=dsqrt((xk-xj)**2 + (yk-yj)**2)/rbfDc
        !rkj=((xk-xj)**2 + (yk-yj)**2)/rbfDc**2
        !G(k,j)=exp(-rbfAc*rkj)
        G(k,j)=(1d0-rkj)**6*(3d0+18d0*rkj+35d0*rkj**2)
      enddo 

      poly(1)=1
      poly(2)=xj
      poly(3)=yj
      poly(4)=xj*xj
      poly(5)=xj*yj
      poly(6)=yj*yj

      G(this%iv+1:this%rbfN,j)=poly
      G(j,this%iv+1:this%rbfN)=poly

      Rm(j)=(1d0-rji)**6*(3d0+18d0*rji+35d0*rji**2)
    enddo

    poly(1)=1
    poly(2)=xi
    poly(3)=yi
    poly(4)=xi*xi
    poly(5)=xi*yi
    poly(6)=yi*yi
    Rm(this%iv+1:this%rbfN)=poly


    !! Matrix inversion by Gauss Elimination
    eye=0d0
    forall(i=1:this%rbfN) eye(i,i)=1d0

    errVar1=0
    ! Forward Elimination
    do i=1,this%rbfN
      if(G(i,i).eq.0d0)then
        write(tf,*) '[ERR] Gauss Elimination diagonal zero loop1'
        errVar1=1
      endif
      do j=i+1,this%rbfN
        tmpr1=G(j,i)/G(i,i)
        G(j,:)=G(j,:)-(tmpr1*G(i,:))
        eye(j,:)=eye(j,:)-tmpr1*eye(i,:)            
      enddo
    enddo

    ! Calculating x for each b
    do i=this%rbfN,1,-1
      if(G(i,i).eq.0d0)then
        write(tf,*) '[ERR] Gauss Elimination diagonal zero loop2'
        errVar1=1
      endif
      do j=i-1,1,-1
        tmpr1=G(j,i)/G(i,i)
        G(j,:)=G(j,:)-(tmpr1*G(i,:))
        eye(j,:)=eye(j,:)-tmpr1*eye(i,:)
      enddo
    enddo
    do k=1,this%rbfN        
      do i=1,this%rbfN
          invG(i,k)=eye(i,k)/G(i,i)
      enddo
    enddo

    ! Phi calculation
    this%phi=0d0
    do i=1,this%iv
      tmpr1=0d0
      do j=1,this%rbfN
        tmpr1=tmpr1+Rm(j)*invG(j,i)
      enddo
      this%phi(i)=tmpr1
    enddo

    if(errVar1.eq.1)then
      write(tf,*)'[ERR] Failed RPIM phi calculation'
      write(tf,*)'[---] Node, iv',this%Rid,this%iv
      write(tf,*)'[---] px,py',this%Rpx,this%Rpy
      write(tf,*)'[---] jv',this%jv(1:this%iv)
      stop
    endif

  end subroutine calPhi

end module rpimModule