!!-------------------------findRadLinkList-------------------------!!
  subroutine findRadLinkList(ip,npt,nnzt,ivq,jvq,corx,cory,coef,&
    rad,dr,cx,cy)  
  use basicVars
  implicit none

    !! Taken from bsnq
    !! Commit : a8396ef
    !! Date   : 2020-05-22

    integer(kind=C_K1),intent(in)::ip,npt,nnzt
    integer(kind=C_K1),intent(in)::ivq(npt+1),jvq(nnzt)
    real(kind=C_K2),intent(in)::corx(npt),cory(npt),coef    
    real(kind=C_K2),intent(out)::rad,dr,cx,cy
    integer(kind=C_K1):: i2, j    

    ! Radius
    rad=0d0
    cx=corx(ip)
    cy=cory(ip)
    !k1=(ip-1)*ivq(0)
    do j=ivq(ip), ivq(ip+1)-1
      i2=jvq(j)      
      dr=(corx(i2)-cx)**2 + (cory(i2)-cy)**2
      if(rad.lt.dr) rad=dr      
    enddo
    if(rad.lt.1e-10)then
      write(tf,'(" [ERR] Check radius calculation at node",I10)')ip
      stop
    endif    
    rad=sqrt(rad)*coef

  end subroutine findRadLinkList
!!-----------------------End findRadLinkList-----------------------!!



!!-------------------------findNeiLinkList-------------------------!!
  subroutine findNeiLinkList(ip,rad,npt,nnzt,ivq,jvq,corx,cory,&
    nnMax,nn,neid,newrk,nedr)
  use basicVars
  implicit none

    !! Taken from bsnq
    !! Commit : a8396ef
    !! Date   : 2020-05-22

    !! find nearest neighbours in radius 'rad' using FEM link list

    integer(kind=C_K1),intent(in)::ip,npt,nnzt,nnMax
    integer(kind=C_K1),intent(in)::ivq(npt+1),jvq(nnzt)
    real(kind=C_K2),intent(in)::rad,corx(npt),cory(npt)

    integer(kind=C_K1),intent(out)::nn,neid(nnMax),newrk(nnMax)
    real(kind=C_K2),intent(out)::nedr(nnMax)

    integer(kind=C_K1)::i, i2, i3, j, j2, j3, k1, k2, k3
    real(kind=C_K2)::px,py,rad2,rat2

    nn=0

    px=corx(ip)
    py=cory(ip)
    rad2=rad**2

    k2=0
    k3=1
    newrk(1)=ip
    nedr(1)=0d0

    301 continue
    k1=k2+1
    k2=k3
    do j=k1,k2
      i2=newrk(j)
      !i3=(i2-1)*ivq(0)      
      do j2 = ivq(i2), (ivq(i2+1)-1)
        if(count(newrk(1:k3).eq.jvq(j2)).eq.0)then
          j3=jvq(j2)
          rat2 = ((corx(j3)-px)**2 + (cory(j3)-py)**2)/rad2
          k3=k3+1
          newrk(k3)=j3
          nedr(k3)=rat2
        endif
      enddo
    enddo

    if(maxval(nedr(k2+1:k3)).lt.1d0) goto 301 !! looping till 2x radius

    do i2=1,k3
      if(nedr(i2).gt.1d0)cycle
      nn=nn+1
      neid(nn)=newrk(i2)
    enddo

    ! write(*,'(I10)')nn
    ! write(*,'(2F20.4)')corx(ip),cory(ip)

    ! do i2=1,nn
    !   j3=neid(i2)
    !   rat2 = sqrt(((corx(j3)-px)**2 + (cory(j3)-py)**2)/rad2)
    !   write(*,'(I10,3F20.4)')j3,rat2,corx(j3),cory(j3)
    ! enddo


  end subroutine findNeiLinkList
!!-----------------------End findNeiLinkList-----------------------!!



!!-----------------------------sweMFree----------------------------!!
  subroutine sweMFree(npt, nnzt, ivf, jvf, corx, cory, pObj)
  use basicVars
  use meshFreeMod
  implicit none

    integer(kind=C_K1),intent(in)::npt, nnzt, ivf(npt+1), jvf(nnzt)
    real(kind=C_K2),intent(in)::corx(npt), cory(npt)
    type(mfPoiTyp),intent(inout)::pObj(npt)
    
    integer(kind=C_K1)::nn,err,ip
    integer(kind=C_K1),allocatable::neid(:),newrk(:)  
    real(kind=C_K2)::cx,cy,rad,coef,dr
    real(kind=C_K2),allocatable::nedr(:),phi(:),phiDx(:),phiDy(:)

    call system_clock(sysClk(6))
    allocate(neid(npt), newrk(npt), nedr(npt))
    allocate(phi(npt),phiDx(npt),phiDy(npt))

    ! Coeff of multipliciation to max rad in linkList
    ! Kept at 0.6 coz of the 9 noded element, 
    ! which is like 4 attached polys
    ! But this coef will not work for the 9th node (cell centre)
    ! For that all neighs are in 1 element
    ! Therefore coef=0.6 will return no neighs except itself
    ! Hence making coef = 2*coef if neigh = 1
    coef=0.65d0 
    do ip=1, npt

      call findRadLinkList(ip, npt, nnzt, ivf, jvf, corx, cory, &
        coef, rad, dr, cx, cy)

      call findNeiLinkList(ip, rad, npt, nnzt, ivf, &
        jvf, corx, cory, npt, nn, neid, newrk, nedr)    

      if(nn.le.3)then
        rad = rad * 2d0
        call findNeiLinkList(ip, rad, npt, nnzt, ivf, &
          jvf, corx, cory, npt, nn, neid, newrk, nedr)            
      endif      

      call mls2DDx(cx, cy, nn, rad, corx(neid(1:nn)), &
        cory(neid(1:nn)), phi(1:nn), phiDx(1:nn), phiDy(1:nn), err)

      if(err.ne.0)then
        write(tf,'(" [ERR] No MFree at node ", I10)')ip
        write(tf,'(" [---] Cx, Cy ",2F15.6)')cx,cy
      endif

      call pObj(ip)%setPoi(nn, nn, ip, cx, cy, rad, neid(1:nn), &
        phi(1:nn), phiDx(1:nn), phiDy(1:nn))

    enddo  

    deallocate(neid,newrk,nedr,phi,phiDx,phiDy)
    call system_clock(sysClk(7))
    write(tf,*)    
    write(tf,*)"[MSG] Done sweMFree"
    write(tf,'(" [SPD] ",F15.4)')1d0*(sysClk(7)-sysClk(6))/sysRate
    write(tf,*)    

  end subroutine sweMFree
!!---------------------------End sweMFree--------------------------!!