subroutine shapeFnc(shF,shFE,shFN,shW)
use basicVars
implicit none
  
  integer(kind=C_K1)::conv(9)
  real(kind=C_K2),intent(out)::shF(9,9),shFE(9,9),shFN(9,9)
  real(kind=C_K2),intent(out)::shW(9)
  real(kind=C_K2)::corE(9),corN(9),cnsNI

  integer(kind=C_K1):: i, j

  !! Note:
  ! For the matrices shF, shFE, shFN, the definition is
  ! shF(i,j) = value of shapeFnc of node i at node j
  ! shFE(i,j) = value of shFE of node i at node j
  ! shFN(i,j) = value of shFN of node i at node j

  !! Node numbering within an element
  ! 4----7----3
  ! |    |    |
  ! 8----9----6
  ! |    |    |
  ! 1----5----2

  conv=(/ 1, 5, 2, 6, 3, 7, 4, 8, 9 /)

  corE=(/ -1d0, 1d0, 1d0, -1d0, 0d0, 1d0, 0d0, -1d0, 0d0 /)
  corN=(/ -1d0, -1d0, 1d0, 1d0, -1d0, 0d0, 1d0, 0d0, 0d0 /)

  cnsNI=1d0/9d0
  shW=(/ 1d0, 1d0, 1d0, 1d0, 4d0, 4d0, 4d0, 4d0, 16d0 /)
  shW=shW*cnsNI

  shF=0d0
  shFE=0d0
  shFN=0d0

  do i=1,9
    shF(i,i)=1d0
    do j=1,9
      if(i.le.4)then
        shFE(i,j)=0.25d0*(2d0*corE(j) + corE(i)) &
          *(corN(j)**2 + corN(j)*corN(i))
        shFN(i,j)=0.25d0*(corE(j)**2 + corE(j)*corE(i)) &
          *(2d0*corN(j) + corN(i))
      elseif((i.eq.5).or.(i.eq.7))then
        shFE(i,j)=-corE(j)*(corN(j)**2 + corN(j)*corN(i))
        shFN(i,j)=0.5d0*(1d0-corE(j)**2) &
          *(2d0*corN(j) + corN(i))
      elseif((i.eq.6).or.(i.eq.8))then
        shFE(i,j)=0.5d0*(2d0*corE(j) + corE(i)) &
          *(1d0 - corN(j)**2)
        shFN(i,j)=-(corE(j)**2 + corE(j)*corE(i))*corN(j)
      else        
        shFE(i,j)=-2d0*corE(j)*(1d0 - corN(j)**2)
        shFN(i,j)=-2d0*(1d0-corE(j)**2)*corN(j)
      endif
    enddo
  enddo

  ! !! Debug comments
  ! write(tf,*)'[DBG] shF'
  ! do i=1,9
  !   write(tf,'(9F15.6)')shF(i,:)
  !   ! i2=conv(i)
  !   ! write(tf,'(9F15.6)')shF(i2,conv)
  ! enddo
  ! write(tf,*)'[DBG] shFE'
  ! do i=1,9
  !   write(tf,'(9F15.6)')shFE(i,:)
  !   ! i2=conv(i)
  !   ! write(tf,'(9F15.6)')shFE(i2,conv)
  ! enddo
  ! write(tf,*)'[DBG] shFN'
  ! do i=1,9
  !   write(tf,'(9F15.6)')shFN(i,:)
  !   ! i2=conv(i)
  !   ! write(tf,'(9F15.6)')shFN(i2,conv)
  ! enddo  
  ! write(tf,*)'[DBG] shW'
  ! write(tf,*)shW

end subroutine shapeFnc


subroutine jacobianInv(npt,nele,conn,coorx,coory,&
  shF,shFE,shFN,shW,jacb,eleArea)
use basicVars
use jacobianModule
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,conn(nele,9)

  real(kind=C_K2),intent(in)::coorx(npt),coory(npt)
  real(kind=C_K2),intent(in)::shF(9,9),shFE(9,9),shFN(9,9)
  real(kind=C_K2),intent(in)::shW(9)
  real(kind=C_K2)::ecx(9),ecy(9)
  type(jacbType),intent(out)::jacb(nele)
  real(kind=C_K2),intent(out)::eleArea(nele)

  integer(kind=C_K1):: iel, na(9), i, j
  real(kind=C_K2):: tmpr1, tmpr2, tmpr3, tmpr4, tmpr5

  !! Note:
  ! jacb(i)%D(j) = value of D inside element i at integration point j

  do iel=1,nele
    na=conn(iel,:)    
    ecx=coorx(na)
    ecy=coory(na)

    do i=1,9 !Looping over integration points
      tmpr1=0d0
      tmpr2=0d0
      tmpr3=0d0
      tmpr4=0d0
      do j=1,9 !Looping over element nodes
        tmpr1=tmpr1+shFE(j,i)*ecx(j)
        tmpr2=tmpr2+shFE(j,i)*ecy(j)
        tmpr3=tmpr3+shFN(j,i)*ecx(j)
        tmpr4=tmpr4+shFN(j,i)*ecy(j)
      enddo
      tmpr5=dabs(tmpr1*tmpr4 - tmpr2*tmpr3)
      if(tmpr5.lt.1d-10)then
        write(tf,*)"[ERR] Jacobian determinant < 1e-10, ele",iel,tmpr5
        stop
      endif
      jacb(iel)%D11(i)=tmpr4/tmpr5
      jacb(iel)%D12(i)=-tmpr2/tmpr5
      jacb(iel)%D21(i)=-tmpr3/tmpr5
      jacb(iel)%D22(i)=tmpr1/tmpr5            
      jacb(iel)%D(i)=tmpr5
    enddo
  enddo

  eleArea=0d0
  do iel=1,nele
    tmpr1=0d0
    do i=1,9
      tmpr1=tmpr1+(shW(i)*jacb(iel)%D(i))
    enddo
    eleArea(iel)=tmpr1
  enddo  

  ! !! Debug comments
  ! write(tf,*)"[DBG] Jacobian determinant"
  ! do iel=1,nele
  !   write(tf,'(I6,9F15.6)')iel,jacb(iel)%D(:)
  !   tmpr1=0d0
  !   do i=1,9
  !     tmpr1=tmpr1+jacb(iel)%D(i)*shW(i)      
  !   enddo
  !   write(tf,*)iel,tmpr1
  ! enddo  


end subroutine jacobianInv