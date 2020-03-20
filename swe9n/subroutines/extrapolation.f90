subroutine midpoiFindNe(npt,nele,nmidpoi,wetpoi,midpoi,&
  coorx,coory,midObj)
use basicVars
use rpimModule
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,nmidpoi
  integer(kind=C_K1),intent(in)::wetpoi(npt),midpoi(npt)  
  integer(kind=C_K1)::tmpia(npt)
  real(kind=C_K2),intent(in)::coorx(npt),coory(npt)
  real(kind=C_K2)::xi,yi
  type(rpimType),intent(inout)::midObj(npt)
  type(rpimType)::midTmp1

  tmpr5=rbfDc**2

  do l2=1,nmidpoi
    i=midpoi(l2)
    xi=coorx(i)
    yi=coory(i)

    i2=0
    do j=1,npt
      if(wetpoi(j).eq.0)cycle
      if(i.eq.j)cycle 
      tmpr1=coorx(j)-xi
      tmpr2=coory(j)-yi
      tmpr3=tmpr1**2 + tmpr2**2      
      if(tmpr3.le.tmpr5)then
        i2=i2+1
        tmpia(i2)=j
      endif
    enddo

    if(i2.gt.rbfMaxNe)then
      write(tf,*)'[ERR] Increase rbfMaxNe for node ',i
      write(tf,*)'[---] Actual vs limit ',i2,rbfMaxNe
      stop
    endif

    call midObj(i)%init(i,xi,yi,i2,tmpia(1:i2))
    call midObj(i)%calPhi(coorx(tmpia(1:i2)),coory(tmpia(1:i2)))
  enddo

  ! !! Debug comments
  ! do l2=1,nmidpoi
  !   i=midpoi(l2)

  !   midTmp1=midObj(i)
  !   write(tf,*)'---Midpoi---',i,midTmp1%iv
  !   write(tf,*)midTmp1%Rpx,midTmp1%Rpy
  !   do j=1,midTmp1%iv
  !     k=midTmp1%jv(j)
  !     write(tf,*)k,coorx(k),coory(k)
  !   enddo
  ! enddo

end subroutine midpoiFindNe