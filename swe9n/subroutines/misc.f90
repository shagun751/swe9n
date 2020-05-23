subroutine hex2dec(text,varlen,result,err)
implicit none
  integer(kind=4),intent(out)::result,err
  integer(kind=4),intent(in)::varlen
  character(len=1),intent(in)::text(varlen)
  !local
  integer(kind=4)::i,j
  character(len=16)::hexstring

  err=0
  hexstring="0123456789abcdef"
  result=0;
  do i=varlen,1,-1
    j=SCAN(hexstring,text(i))-1    
    if(j.lt.0) then
      err=1
    endif
    result=result+(j*16**(varlen-i))
  enddo
end subroutine hex2dec



subroutine mergeSort(maxNePoi,nNePoi,poi2poi)
implicit none
  integer(kind=4),intent(in)::maxNePoi,nNePoi
  integer(kind=4),intent(inout)::poi2poi(maxNePoi)
  integer(kind=4)::temp(maxNePoi),i,j,k,l,m
  integer(kind=4)::istart,imiddle,iend,width  

  temp=0
  width=1
  do while(width.lt.nNePoi)
    do i=1,nNePoi,(2*width)
      istart=i
      imiddle=min(i+width-1,nNePoi)
      iend=min(i+2*width-1,nNePoi)

      j=istart
      k=imiddle+1
      do l=istart,iend
        if((j.le.imiddle).and.(k.gt.iend .or. poi2poi(j).le.poi2poi(k))) then
          temp(l)=poi2poi(j)
          j=j+1
        else
          temp(l)=poi2poi(k)
          k=k+1
        endif
      enddo  
      !write(9,*)(temp(m),m=1,maxNePoi)
    enddo    
    poi2poi=temp
    width=2*width
  enddo

end subroutine mergeSort



subroutine printMat(npt,nnzt,ivf,jvf,gMat)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt,nnzt
  integer(kind=C_K1),intent(in)::ivf(npt+1),jvf(nnzt)
  real(kind=C_K2),intent(in)::gMat(nnzt)

  integer(kind=C_K1):: i, k, k2, j  

  do i=1,npt
    k=ivf(i)
    k2=ivf(i+1)-1
    write(tf,*)i,':'
    do j=k,k2
      write(tf,'(F15.6)',advance='no')gMat(j)
    enddo
    write(tf,*)
  enddo
end subroutine printMat



subroutine chkMat(npt,nnzt,ivf,jvf,gN)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt,nnzt
  integer(kind=C_K1),intent(in)::ivf(npt+1),jvf(nnzt)
  real(kind=C_K2),intent(in)::gN(nnzt)  

  integer(kind=C_K1):: i, k, k2, j, l  

  do i=1,npt
    k=ivf(i)
    k2=ivf(i+1)-1    
    l=0
    do j=k,k2
      if(abs(gN(j)).gt.1d5)then              
        l=1
        exit
      endif
    enddo
    if(l.eq.1)then
      write(tf,*)"[ERR] Node ",i
      write(tf,*)"[---] ",gN(k:k2)
    endif
  enddo

end subroutine chkMat



subroutine eleToGlob(npt,nele,nnzt,conn,ivf,jvf,elejvf9x9)
use basicVars
implicit none
  
  integer(kind=C_K1),intent(in)::npt,nele,nnzt,conn(nele,9)
  integer(kind=C_K1),intent(in)::ivf(npt+1),jvf(nnzt)
  integer(kind=C_K1),intent(out)::elejvf9x9(nele,81)

  integer(kind=C_K1):: iel, na(9), i, i2, k, k2, j, j2
  integer(kind=C_K1):: l, l2

  do iel=1,nele
    na=conn(iel,:)

    !! 9x9
    do i=1,9
      i2=na(i)
      k=ivf(i2)
      k2=ivf(i2+1)-1
      do j=1,9
        j2=na(j)
        l2=(i-1)*9+j
        do l=k,k2
          if(jvf(l).eq.j2)goto 101
        enddo
        write(tf,*)'[ERR] rh1: Check negh for node',i2
        write(tf,*)'[---] ',jvf(k:k2)
        101 continue
        elejvf9x9(iel,l2)=l
      enddo
    enddo
  enddo
  
end subroutine