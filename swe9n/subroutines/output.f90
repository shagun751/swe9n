subroutine out4NXML(probname,npoinl,npoint,nelem,code,ts,&
  conn,coorx,coory,u,v,eta,dep,wetpoi,pc,windTx,windTy,pObj)
use meshFreeMod
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem
  integer(kind=4),intent(in)::conn(nelem,9),code,ts
  integer(kind=4),intent(in)::wetpoi(npoint)
  integer(kind=4)::i,j,k

  real(kind=8),intent(in)::coorx(npoint),coory(npoint)
  real(kind=8),intent(in)::u(npoint),v(npoint),eta(npoint)
  real(kind=8),intent(in)::dep(npoint)  
  real(kind=8),intent(in)::pc(npoint),windTx(npoint),windTy(npoint)
  real(kind=8)::sDx,sDy

  type(mfPoiTyp),intent(in)::pObj(npoint)

  character(len=256),intent(in)::probname
  character(len=256)::text


  write(text,'(I10)')ts
  text=adjustl(text)
  text="Output/"//probname(1:len_trim(probname))//"_"//text(1:len_trim(text))//".vtu"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
  write(code,'(T3,a)')'<UnstructuredGrid>'
  write(code,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',npoinl,'" NumberOfCells="',nelem,'">'


  write(code,'(T5,a)')'<PointData Scalars="eta" Vectors="vel">'
  
  write(code,'(T7,a)')'<DataArray type="Float64" Name="eta" format="ascii">'
  write(code,*)eta(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="depth" format="ascii">'
  write(code,*)-dep(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="Pres" format="ascii">'
  write(code,*)pc(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="wetpoi" format="ascii">'
  write(code,'(I4)')wetpoi(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  ! write(code,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
  ! write(code,*)porH-depth(1:npoinl)
  ! write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="vel" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,'(2F20.6,F4.1)')u(i),v(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="windSh" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,'(2F20.6,F4.1)')windTx(i),windTy(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'

  ! Test derivative
  ! write(code,'(T7,a)')'<DataArray type="Float64" Name="etaD" NumberOfComponents="3" format="ascii">'  
  ! do i=1,npoinl
  !   sDx=0d0
  !   sDy=0d0
  !   do j = 1, pObj(i)%nn
  !     sDx = sDx + pObj(i)%phiDx(j) * eta( pObj(i)%neid(j) )
  !     sDy = sDy + pObj(i)%phiDy(j) * eta( pObj(i)%neid(j) )
  !   enddo
  !   write(code,'(2F20.6,F4.1)')sDx,sDy,0
  ! enddo
  ! write(code,'(T7,a)')'</DataArray>'
  
  write(code,'(T5,a)')'</PointData>'  

  write(code,'(T5,a)')'<CellData>'
  write(code,'(T5,a)')'</CellData>'

  write(code,'(T5,a)')'<Points>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,'(2F20.6,F4.1)')coorx(i),coory(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Points>'

  write(code,'(T5,a)')'<Cells>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="connectivity" format="ascii">'  
  do i=1,nelem
    write(code,*)conn(i,1:4)-1
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="offsets" format="ascii">'  
  do i=1,nelem
    write(code,*)4*i
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
  do i=1,nelem
    write(code,'(I4)')9
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Cells>'

  write(code,'(T5,a)')'</Piece>'
  write(code,'(T3,a)')'</UnstructuredGrid>'
  write(code,'(a)')'</VTKFile>'
  close(code)

end subroutine out4NXML


subroutine out8NXML(probname,npoinl,npoint,nelem,code,ts,&
  conn,coorx,coory,u,v,eta,dep)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem
  integer(kind=4),intent(in)::conn(nelem,9),code,ts
  integer(kind=4)::i,j,k,endInd

  real(kind=8),intent(in)::coorx(npoint),coory(npoint)
  real(kind=8),intent(in)::u(npoint),v(npoint),eta(npoint)
  real(kind=8),intent(in)::dep(npoint)  

  character(len=256),intent(in)::probname
  character(len=256)::text


  write(text,'(I10)')ts
  text=adjustl(text)
  text="Output/"//probname(1:len_trim(probname))//"_"//text(1:len_trim(text))//".vtu"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
  write(code,'(T3,a)')'<UnstructuredGrid>'
  write(code,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',endInd,'" NumberOfCells="',nelem,'">'

  endInd=npoint-nelem

  write(code,'(T5,a)')'<PointData Scalars="eta" Vectors="vel">'
  
  write(code,'(T7,a)')'<DataArray type="Float64" Name="eta" format="ascii">'
  write(code,*)eta(1:endInd)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="depth" format="ascii">'
  write(code,*)-dep(1:endInd)
  write(code,'(T7,a)')'</DataArray>'

  ! write(code,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
  ! write(code,*)porH-depth(1:npoinl)
  ! write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="vel" NumberOfComponents="3" format="ascii">'  
  do i=1,endInd
    write(code,*)u(i),v(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  
  write(code,'(T5,a)')'</PointData>'  

  write(code,'(T5,a)')'<CellData>'
  write(code,'(T5,a)')'</CellData>'

  write(code,'(T5,a)')'<Points>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
  do i=1,endInd
    write(code,*)coorx(i),coory(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Points>'

  write(code,'(T5,a)')'<Cells>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="connectivity" format="ascii">'  
  do i=1,nelem
    write(code,*)conn(i,1:8)-1
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="offsets" format="ascii">'  
  do i=1,nelem
    write(code,*)8*i
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
  do i=1,nelem
    write(code,*)23
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Cells>'

  write(code,'(T5,a)')'</Piece>'
  write(code,'(T3,a)')'</UnstructuredGrid>'
  write(code,'(a)')'</VTKFile>'
  close(code)

end subroutine out8NXML

subroutine out4NXMLDiff(probname,npoinl,npoint,nelem,code,ts,&
  conn,coorx,coory,depNew,depOld,depChange,cour)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem
  integer(kind=4),intent(in)::conn(nelem,9),code,ts  
  integer(kind=4)::i,j,k

  real(kind=8),intent(in)::coorx(npoint),coory(npoint)
  real(kind=8),intent(in)::depNew(npoint)
  real(kind=8),intent(in)::depOld(npoint)  
  real(kind=8),intent(in)::depChange(npoint),cour(npoint)

  character(len=256),intent(in)::probname
  character(len=256)::text


  write(text,'(I10)')ts
  text=adjustl(text)
  text="Output/"//probname(1:len_trim(probname))//"_BathySmooth"//".vtu"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
  write(code,'(T3,a)')'<UnstructuredGrid>'
  write(code,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',npoinl,'" NumberOfCells="',nelem,'">'


  write(code,'(T5,a)')'<PointData Scalars="depNew" Vectors="vel">'
  
  write(code,'(T7,a)')'<DataArray type="Float64" Name="depNew" format="ascii">'
  write(code,*)-depNew(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="depOld" format="ascii">'
  write(code,*)-depOld(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="depPercChange" format="ascii">'
  write(code,*)depChange(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="Courant" format="ascii">'
  write(code,*)cour(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  ! write(code,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
  ! write(code,*)porH-depth(1:npoinl)
  ! write(code,'(T7,a)')'</DataArray>'  
  
  write(code,'(T5,a)')'</PointData>'  

  write(code,'(T5,a)')'<CellData>'
  write(code,'(T5,a)')'</CellData>'

  write(code,'(T5,a)')'<Points>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,*)coorx(i),coory(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Points>'

  write(code,'(T5,a)')'<Cells>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="connectivity" format="ascii">'  
  do i=1,nelem
    write(code,*)conn(i,1:4)-1
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="offsets" format="ascii">'  
  do i=1,nelem
    write(code,*)4*i
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
  do i=1,nelem
    write(code,*)9
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Cells>'

  write(code,'(T5,a)')'</Piece>'
  write(code,'(T3,a)')'</UnstructuredGrid>'
  write(code,'(a)')'</VTKFile>'
  close(code)

end subroutine out4NXMLDiff