subroutine middleNodes(npl,npq,npt,nele,coorx,coory,dep,conn)
use basicVars
implicit none
  
  integer(kind=C_K1),intent(in)::npl,nele
  integer(kind=C_K1),intent(inout)::npq,npt,conn(nele,9)  
  integer(kind=C_K1),allocatable::mat1(:),mat2a(:,:),mat2b(:,:)
  integer(kind=C_K1)::mat3(4,4)
  integer(kind=C_K1)::ln1,ln2,qp
  integer(kind=C_K1):: tmpi1, iel, i, j, i2, k
  real(kind=C_K2):: tmpr1, tmpr2, tmpr3

  real(kind=C_K2),intent(inout)::coorx(npt),coory(npt),dep(npt)

  allocate(mat1(npl),mat2a(npl,maxNePoi),mat2b(npl,maxNePoi))

  tmpi1=npt
  npq=0
  npt=npl

  mat1=0
  mat2a=0
  mat2b=0
  mat3=0
  mat3(1,:)=(/ 0, 5, 0, 8 /)
  mat3(3,:)=(/ 0, 6, 0, 7 /)

  !! Edge nodes
  do iel=1,nele
    do i=1,3,2
      do j=2,4,2
        ln1=conn(iel,i)
        ln2=conn(iel,j)

        do i2=1,mat1(ln1)
          if(mat2a(ln1,i2).eq.ln2) then
            qp=mat2b(ln1,i2)
            goto 11
          endif
        enddo

        npt=npt+1
        if(npt.gt.tmpi1)then
          write(tf,*)"[ERR] npt > limit in edge nodes"
          stop
        endif
        qp=npt
        coorx(qp)=0.5d0*(coorx(ln1)+coorx(ln2))
        coory(qp)=0.5d0*(coory(ln1)+coory(ln2))
        dep(qp)=0.5d0*(dep(ln1)+dep(ln2))
        
        mat1(ln1)=mat1(ln1)+1
        mat2a(ln1,mat1(ln1))=ln2
        mat2b(ln1,mat1(ln1))=qp
        mat1(ln2)=mat1(ln2)+1
        mat2a(ln2,mat1(ln2))=ln1
        mat2b(ln2,mat1(ln2))=qp

        11 conn(iel,mat3(i,j))=qp

      enddo
    enddo
  enddo

  !! Centre node
  do iel=1,nele

    npt=npt+1
    if(npt.gt.tmpi1)then
      write(tf,*)"[ERR] npt > limit in centre nodes"
      stop
    endif
    conn(iel,9)=npt

    tmpr1=0d0
    tmpr2=0d0
    tmpr3=0d0
    do i=1,4
      k=conn(iel,i)
      tmpr1=tmpr1+coorx(k)
      tmpr2=tmpr2+coory(k)
      tmpr3=tmpr3+dep(k)
    enddo

    coorx(npt)=tmpr1*0.25d0
    coory(npt)=tmpr2*0.25d0
    dep(npt)=tmpr3*0.25d0
  enddo

  npq=npt-npl

  ! !! Debug comments
  ! write(tf,*)"Elements"
  ! do i=1,nele
  !   write(tf,'(I6," : ",9I6)')i,conn(i,:)
  ! enddo  
  ! write(tf,*)"Element node coords"
  ! do i=1,nele
  !   write(tf,'("Ele : ",I6)')i
  !   do j=1,9
  !     k=conn(i,j)
  !     write(tf,'("    : ",I6,3F15.6)')k,coorx(k),coory(k),dep(k)
  !   enddo
  ! enddo  

  deallocate(mat1,mat2a,mat2b)

  write(tf,*)"[MSG] Done generation of middle nodes"


end subroutine middleNodes



subroutine bndMiddleNodes(nele,nbnd,conn,mabnd)
use basicVars
implicit none
    
  integer(kind=C_K1),intent(in)::nele,nbnd,conn(nele,9)
  integer(kind=C_K1),intent(inout)::mabnd(nbnd,6)
  integer(kind=C_K1)::mat1(4),mat2(4,4)
  integer(kind=C_K1):: i2, n1, n2, iel, i, k2, k  

  mat1=(/ 2,3,4,1 /)    

  do i2=1,nbnd
    n1=mabnd(i2,1)
    n2=mabnd(i2,2)
    iel=mabnd(i2,3)    

    do i=1,4
      if(conn(iel,i).eq.n1) k=i
      if(conn(iel,i).eq.n2) k2=i
    enddo

    if(k2.ne.mat1(k))then
      write(tf,*)'[ERR] Chk bnd nodes not anticlk for ele',iel
      write(tf,*)'[---] conn',conn(iel,:)
      write(tf,*)'[---] mabnd',mabnd(i2,:)
      stop
    endif

    mabnd(i2,5)=conn(iel,k+4)
    mabnd(i2,6)=k
  enddo

  ! !! Debug comments
  ! write(tf,*)"[DBG] Elements"
  ! do i=1,nele
  !   write(tf,'(I6," : ",9I6)')i,conn(i,:)
  ! enddo  
  ! write(tf,*)"[DBG] Boundaries"
  ! do i=1,nbnd    
  !   write(tf,'(I6," : ",6I6)')i,mabnd(i,:)
  ! enddo
  write(tf,*)"[MSG] Done bnd middle nodes"  

end subroutine bndMiddleNodes



subroutine nodeConn(npl,npt,nele,conn,poi2poi,npoisur)
use basicVars
implicit none
  
  integer(kind=C_K1),intent(in)::npl,npt,nele,conn(nele,9)
  integer(kind=C_K1),intent(out)::poi2poi(npt,maxNePoi)
  integer(kind=C_K1),intent(out)::npoisur(npt,3)
  integer(kind=C_K1):: iel, na(9), n1, n2, i2, j, l  

  poi2poi=0
  npoisur=0

  do iel=1,nele
    na=conn(iel,:)
    do j=1,9
      n1=na(j)
      do l=1,9
        n2=na(l)
        if(n1.ne.n2) then
          do i2=1,maxNePoi
            if(poi2poi(n1,i2).eq.0) then
              npoisur(n1,1)=npoisur(n1,1)+1
              poi2poi(n1,i2)=n2
              exit    
            elseif(poi2poi(n1,i2).eq.n2) then
              exit
            elseif(i2.eq.maxNePoi) then
              write(tf,*) "[ERR] Increase maxNePoi at node ",n1
              stop
            endif      
          enddo
        endif
      enddo            
    enddo
  enddo

end subroutine nodeConn



subroutine bndNormal(npt,nbnd,nbndpoi,mabnd,coorx,coory,&
  bnd11p,bnd12p,bnd14p,bndLen,bndpNm)
use basicVars
implicit none
    
  integer(kind=C_K1),intent(in)::npt,nbnd,nbndpoi
  integer(kind=C_K1),intent(in)::mabnd(nbnd,6)
  integer(kind=C_K1),intent(out)::bnd11p(0:nbndpoi)
  integer(kind=C_K1),intent(out)::bnd12p(0:nbndpoi)
  integer(kind=C_K1),intent(out)::bnd14p(0:nbndpoi)  
  integer(kind=C_K1),allocatable::tmpia(:,:)
  real(kind=C_K2),intent(in)::coorx(npt),coory(npt)
  real(kind=C_K2),intent(out)::bndpNm(npt,2),bndLen(nbnd)

  integer(kind=C_K1):: k, i, na(9), l, j, j2, l2, l4
  integer(kind=C_K1)::bndPref(5), pref1, pref2
  real(kind=C_K2):: tmpr1, tmpr2

  allocate(tmpia(nbndpoi,2))

  bndPref = (/ 14, 11, 12, 13, 0 /)

  bnd11p=0
  bnd12p=0
  bnd14p=0
  bndpNm=0d0
  tmpia=0

  k=0
  do i=1,nbnd
    na(1)=mabnd(i,1)
    na(2)=mabnd(i,2)
    na(3)=mabnd(i,5)
    l=mabnd(i,4)  
    bndLen(i)=dsqrt( (coorx(na(2))-coorx(na(1)))**2 &
      + (coory(na(2))-coory(na(1)))**2 )  

    tmpr1=(coory(na(2))-coory(na(1)))*bndLen(i)
    tmpr2=(coorx(na(1))-coorx(na(2)))*bndLen(i)
    do j=1,3      
      bndpNm(na(j),1)=bndpNm(na(j),1)+tmpr1
      bndpNm(na(j),2)=bndpNm(na(j),2)+tmpr2

      do j2=1,k
        if(tmpia(j2,1).eq.na(j)) goto 111
      enddo
      k=k+1
      if(k.gt.nbndpoi)then
        write(tf,*)"[ERR] Check bnd not closed in loop"        
        write(tf,*)"[---] nbnd, k, nbndpoi",nbnd,k,nbndpoi
        stop
      endif
      tmpia(k,1)=na(j)
      tmpia(k,2)=l
      j2=k
      111 continue      

      do pref1=1,5
        if(tmpia(j2,2).eq.bndPref(pref1)) exit
      enddo

      do pref2=1,5
        if(l.eq.bndPref(pref2)) exit
      enddo

      tmpia(j2,2) = bndPref( min(pref1, pref2) )      

      if(tmpia(j2,2).eq.0)then
        write(tf,*)"[ERR] Bnd type 0 at node",j2
        write(tf,*)"[ERR] ",l,tmpia(j2,2)
        stop
      endif

    enddo
  enddo

  if(k.ne.nbndpoi)then
    write(tf,*)"[ERR] Check bnd not closed"    
    write(tf,*)"[---] nbnd, k, nbndpoi",nbnd,k,nbndpoi
    stop
  endif

  call normaliseNorm(npt,nbndpoi,tmpia(:,1),bndpNm)

  l=0
  l2=0
  l4=0
  do i=1,nbndpoi
    if(tmpia(i,2).eq.11)then
      l=l+1
      bnd11p(l)=tmpia(i,1)
    
    elseif(tmpia(i,2).eq.12)then
      l2=l2+1
      bnd12p(l2)=tmpia(i,1)
    
    elseif(tmpia(i,2).eq.14)then
      l4=l4+1
      bnd14p(l4)=tmpia(i,1)
    
    else
      write(tf,*)"[ERR] Unknown bndNode type",tmpia(i,:)
      stop
    endif
  enddo
  bnd11p(0)=l
  bnd12p(0)=l2
  bnd14p(0)=l4
  
  ! !!Debug comments
  ! write(tf,*)'[DBG] bnd11p',bnd11p(0)
  ! do i=1,bnd11p(0)
  !   k=bnd11p(i)
  !   tmpr1=bndpNm(k,1)
  !   tmpr2=bndpNm(k,2)
  !   tmpr3=dsqrt(tmpr1**2 +tmpr2**2)
  !   write(tf,'(I10,3F15.6)')k,bndpNm(k,:),tmpr3
  ! enddo
  ! write(tf,*)'[DBG] bnd12p',bnd12p(0)
  ! do i=1,bnd12p(0)
  !   k=bnd12p(i)
  !   tmpr1=bndpNm(k,1)
  !   tmpr2=bndpNm(k,2)
  !   tmpr3=dsqrt(tmpr1**2 +tmpr2**2)
  !   write(tf,'(I10,3F15.6)')k,bndpNm(k,:),tmpr3
  ! enddo
  ! stop

  write(tf,*)'[MSG] BndNode sorting and normal done'
end subroutine bndNormal



subroutine normaliseNorm(npt,nbndp,bndp,bndpNm)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt,nbndp,bndp(nbndp)
  real(kind=C_K2),intent(inout)::bndpNm(npt,2)

  integer(kind=C_K1):: i, k
  real(kind=C_K2):: tmpr1, tmpr2, tmpr3

  do i=1,nbndp
    k=bndp(i)
    tmpr1=bndpNm(k,1)
    tmpr2=bndpNm(k,2)
    tmpr3=dsqrt(tmpr1**2 +tmpr2**2)
    bndpNm(k,1)=tmpr1/tmpr3
    bndpNm(k,2)=tmpr2/tmpr3
  enddo

end subroutine normaliseNorm



subroutine calcCourant(npt,nele,conn,coorx,coory,dep,dt,cour)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,conn(nele,9)
  real(kind=C_K2),intent(in)::coorx(npt),coory(npt)
  real(kind=C_K2),intent(in)::dep(npt),dt
  real(kind=C_K2),intent(out)::cour(npt)
  real(kind=C_K2),allocatable::couA(:),couB(:)

  integer(kind=C_K1):: iel, na(9)
  real(kind=C_K2):: tmpr1

  allocate(couA(npt),couB(npt))

  couA=0d0
  couB=0d0
  do iel=1,nele
    na=conn(iel,:)
    tmpr1=(coorx(na(1))*coory(na(2))-coorx(na(2))*coory(na(1)))
    tmpr1=tmpr1+(coorx(na(2))*coory(na(3))-coorx(na(3))*coory(na(2)))
    tmpr1=tmpr1+(coorx(na(3))*coory(na(4))-coorx(na(4))*coory(na(3)))
    tmpr1=tmpr1+(coorx(na(4))*coory(na(1))-coorx(na(1))*coory(na(4)))
    tmpr1=abs(tmpr1)/2d0
    couA(na)=couA(na)+tmpr1**2
    couB(na)=couB(na)+tmpr1
  enddo  
  couA=dsqrt(couA/couB)/2d0 !! Coz quad nodes
  cour=dsqrt(grav*dep)*dt/couA

  deallocate(couA,couB)

end subroutine calcCourant