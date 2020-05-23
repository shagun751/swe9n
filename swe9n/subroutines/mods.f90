!!--------------------------- basicVars ---------------------------!!
module basicVars
implicit none

  integer,parameter::C_K1=4,C_K2=8    
  integer::UTMZone=44
  integer(kind=C_K1)::ifli
  integer(kind=C_K1),parameter::tf=10,ifl(10)=(/ (ifli, ifli=11,20)  /)
  integer(kind=C_K1)::sysClk(10)    
  
  real(kind=C_K2)::sysRate


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
  !!--------------------End Constants---------------------!!


  !! File definitions
  ! tf=10                      !Terminal Output File
  ! ifl(1)=11                  !Input file
  ! ifl(2)=12                  !Mesh file
  ! ifl(3)=13                  !Output file
  ! ifl(4)=14                  !Cyclone track file
  ! ifl(5)=15                  !Probe Output File


end module basicVars
!!------------------------- End basicVars -------------------------!!



!!------------------------ jacobianModule -------------------------!!
module jacobianModule
use basicVars
implicit none
  
  type jacbType
    real(kind=C_K2)::D(9),D11(9),D12(9),D21(9),D22(9)    
  end type jacbType

end module jacobianModule
!!---------------------- End jacobianModule -----------------------!!
