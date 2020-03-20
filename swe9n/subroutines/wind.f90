!********************************************************
! SUBROUTINE TO DETERMINE WIND STRESSES FOR A MOVING CYCLONE
!***********************************************************

! SUBROUTINE WIND(X,Y,NODE,TX,TY,THE,IT,DT,X0,Y0,N1,V,ALL,WM,R0)
!   use basicVars
!   DIMENSION X(1),Y(1),TX(1),TY(1)

!   OPEN(22,FILE='center.DAT',STATUS='UNKNOWN')    
!   PHI=THE*pi/180.
!   RHO=1.20E-03
!   AKA=1.8E-03
!   R1=2.2*R0
!   ALP=2.5*R0
!   BET=0.135*R0
!   ALL=ALL*pi/180.

!   ! Shagun modify
!   X0=X0+V*(IT-N1)*2.*DT*COS(PHI)
!   Y0=Y0+V*(IT-N1)*2.*DT*SIN(PHI)
!   ! X0=X0+V*(IT-N1)*DT*COS(PHI)
!   ! Y0=Y0+V*(IT-N1)*DT*SIN(PHI)
     
!   DO I=1,NODE
!     TX(I)=0.
!     TY(I)=0.
!     R=SQRT((X(I)-X0)**2+(Y(I)-Y0)**2)
!     IF(R.EQ.0.)GOTO 10
!     IF(R.LE.R0)AL=ALL
!     IF(R.GT.R0)AL=ALL

!     AR=R/R0
!     BR=(R0-R)/ALP  
!     CR=EXP((R0-R1)/ALP)*EXP((R1-R)/BET)

!     IF(R.LE.R0)W=WM*AR*SQRT(AR)
!     IF(R.GT.R0.AND.R.LE.R1)W=WM*EXP(BR)
!     IF(R.GT.R1)W=WM*CR
!     CSL=COS(AL)
!     SNL=SIN(AL)
!     SN=(Y(I)-Y0)/R
!     CS=(X(I)-X0)/R
!     TX(I)=-RHO*AKA*W*W*(SN*CSL+CS*SNL)
!     TY(I)=RHO*AKA*W*W*(CS*CSL-SN*SNL)
!   ENDDO
!   N1=IT
  
! END SUBROUTINE WIND

subroutine windNew2(npt,coorx,coory,windLon,windLat,windR0, &
  coriF,windWm,windA,windB,pr,windTx,windTy)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt
  real(kind=c_K2),intent(in)::coorx(npt),coory(npt)
  real(kind=c_K2),intent(in)::windLon,windLat,windR0,windWm
  real(kind=c_K2),intent(in)::windA,windB,coriF(npt)
  real(kind=c_K2),intent(out)::windTx(npt),windTy(npt),pr(npt)
  real(kind=c_K2)::windX0,windY0
  real(kind=c_K2)::rho,cd,rhoA,pDrop,al,pn,pc,dr,w
  real(kind=c_K2)::csl,snl,sn,cs,pmax,rMax

  rhoA=1.1d0
  rho=rhoA !! I think this is rho_air/rho_water
  al=20d0*pi/180d0
  pn=1.0135*10**5
  
  pDrop=rhoA*dexp(1d0)*windWm**2/windB
  pc=pn-pDrop
  rMax=windA**(1d0/windB)

  tmpr1=windLon
  tmpr2=windLat
  call utmGeo(tmpr1,tmpr2,windX0,windY0,1,UTMZone)

  !write(*,*)windA,windB,pDrop,windWm

  windTx=0d0
  windTy=0d0
  pr=0d0  

  tmpr1=windA*windB*pDrop
  do i=1,npt
    dr=dsqrt((coorx(i)-windX0)**2 + (coory(i)-windY0)**2)
    if(dr.eq.0d0)then 
      pr(i)=pc
      cycle
    endif     

    tmpr2=dr**windB
    w=tmpr1*dexp(-windA/tmpr2)/(rhoA*tmpr2)
    w=w+(dr**2)*coriF(i)**2/4d0
    w=dsqrt(w)-(dr*coriF(i)/2d0)

    csl=dcos(al)
    snl=dsin(al)
    sn=(coory(i)-windY0)/dr
    cs=(coorx(i)-windX0)/dr

    if(w.lt.50.33)then
      cd=0.0001d0*(-0.016d0*w*w+0.967d0*w+8.058d0)
    else
      cd=4.1047d0/w/w
    endif    

    windTx(i)=-rho*cd*w*w*(sn*csl+cs*snl)
    windTy(i)=rho*cd*w*w*(cs*csl-sn*snl)
    if(dr.gt.10d0*rMax)then
      windTx(i)=0d0
      windTy(i)=0d0
    endif
    
    pr(i)=pc+pDrop*dexp(-windA/tmpr2)

  enddo

end subroutine windNew2

subroutine windNew2b(npt,coorx,coory,windLon,windLat,windR0, &
  coriF,windWm,windA,windB,pr,windTx,windTy)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt
  real(kind=c_K2),intent(in)::coorx(npt),coory(npt)
  real(kind=c_K2),intent(in)::windLon,windLat,windR0,windWm
  real(kind=c_K2),intent(in)::windA,windB,coriF(npt)
  real(kind=c_K2),intent(out)::windTx(npt),windTy(npt),pr(npt)
  real(kind=c_K2)::windX0,windY0
  real(kind=c_K2)::rho,cd,rhoA,pDrop,al,pn,pc,dr,w
  real(kind=c_K2)::csl,snl,sn,cs,pmax,rMax

  rhoA=1.225d0
  rho=rhoA !! I think this is rho_air/rho_water
  al=20d0*pi/180d0
  pn=1.0135*10**5
  
  pDrop=rhoA*dexp(1d0)*windWm**2/windB
  pc=pn-pDrop
  rMax=windA**(1d0/windB)

  tmpr1=windLon
  tmpr2=windLat
  call utmGeo(tmpr1,tmpr2,windX0,windY0,1,UTMZone)

  !write(*,*)windA,windB,pDrop,windWm

  windTx=0d0
  windTy=0d0
  pr=0d0  

  tmpr1=windA*windB*pDrop
  do i=1,npt
    dr=dsqrt((coorx(i)-windX0)**2 + (coory(i)-windY0)**2)
    if(dr.eq.0d0)then 
      pr(i)=pc
      cycle
    endif     

    tmpr2=dr**windB
    w=tmpr1*dexp(-windA/tmpr2)/(rhoA*tmpr2)
    w=w+(dr**2)*coriF(i)**2/4d0
    w=dsqrt(w)-(dr*coriF(i)/2d0)

    csl=dcos(al)
    snl=dsin(al)
    sn=(coory(i)-windY0)/dr
    cs=(coorx(i)-windX0)/dr

    cd=0.001d0*max(0.067d0*w+0.75d0,3d0)

    windTx(i)=-rho*cd*w*w*(sn*csl+cs*snl)
    windTy(i)=rho*cd*w*w*(cs*csl-sn*snl)
    if(dr.gt.10d0*rMax)then
      windTx(i)=0d0
      windTy(i)=0d0
    endif
    
    pr(i)=pc+pDrop*dexp(-windA/tmpr2)

  enddo

end subroutine windNew2b

subroutine windNew3(npt,lon,lat,coorx,coory,windLon,windLat,windR0,&
  coriF,windWm,windA,windB,pr,windTx,windTy)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt
  integer(kind=C_K1)::is,secInt(6,2)
  real(kind=c_K2),intent(in)::coorx(npt),coory(npt)
  real(kind=c_K2),intent(in)::lon(npt),lat(npt)
  real(kind=c_K2),intent(in)::windLon(3),windLat(3)
  real(kind=c_K2),intent(in)::windR0(3),windWm(3)  
  real(kind=c_K2),intent(in)::windA,windB,coriF(npt)
  real(kind=c_K2),intent(out)::windTx(npt),windTy(npt),pr(npt)
  real(kind=c_K2)::windX0,windY0
  real(kind=c_K2)::rho,cd,rhoA,pDrop,pn,pc,dr,w
  real(kind=c_K2)::csl,snl,sn,cs,pmax,rMax
  real(kind=c_K2)::dir1,dir2,dirStorm,dirNode
  real(kind=c_K2),parameter::cdLimit=0.0035d0 !!Garratt Limit
  real(kind=c_K2)::secLim(6,2),secWei(3),secCd(3)

  secLim(1,:)=(/ 0d0, 40d0  /)
  secLim(2,:)=(/ 40d0, 130d0  /)
  secLim(3,:)=(/ 130d0, 170d0  /)
  secLim(4,:)=(/ 170d0, 220d0  /)
  secLim(5,:)=(/ 220d0, 260d0  /)
  secLim(6,:)=(/ 260d0, 0d0  /)

  secInt=0
  secInt(1,:)=(/ 3,1 /)
  secInt(3,:)=(/ 1,2 /)
  secInt(5,:)=(/ 2,3 /)

  rhoA=1.225d0
  rho=rhoA !! I think this is rho_air/rho_water  
  pn=1.0135*10**5
  
  pDrop=rhoA*dexp(1d0)*windWm(1)**2/windB
  pc=pn-pDrop
  rMax=windA**(1d0/windB)
  !write(*,*)windA,windB,pDrop,windWm

  tmpr1=windLon(1)
  tmpr2=windLat(1)
  call utmGeo(tmpr1,tmpr2,windX0,windY0,1,UTMZone)

  dir1=atan2((windLat(2)-windLat(3)),(windLon(2)-windLon(3)))
  dir2=atan2((windLat(1)-windLat(2)),(windLon(1)-windLon(2)))
  dirStorm=atan((sin(dir1)+sin(dir2)),(cos(dir1)+cos(dir2)))
  dirStorm=dirStorm/deg2rad
  dirStorm=90d0-dirStorm
  if(dirStorm.lt.0d0) dirStorm=dirStorm+360d0  

  windTx=0d0
  windTy=0d0
  pr=0d0  

  tmpr1=windA*windB*pDrop
  do i=1,npt
    dr=dsqrt((coorx(i)-windX0)**2 + (coory(i)-windY0)**2)
    if(dr.eq.0d0)then 
      pr(i)=pc
      cycle
    endif     

    tmpr2=dr**windB
    w=tmpr1*dexp(-windA/tmpr2)/(rhoA*tmpr2)
    w=w+(dr**2)*coriF(i)**2/4d0
    w=dsqrt(w)-(dr*coriF(i)/2d0)

    dirNode=atan2((lat(i)-windLat(1)),(lon(i)-windLon(1)))
    dirNode=dirNode/deg2rad
    dirNode=90d0-dirNode
    if(dirNode.lt.0d0) dirNode=dirNode+360d0

    !! Compute weight for each of the three sectors
    secWei=0d0
    do is=1,6
      dir1=mod(dirStorm+secLim(is,1),360d0)
      dir2=mod(dirStorm+secLim(is,2),360d0)

      if(dir1.gt.dir2)then
        if((dir1.le.dirNode).and.(dirNode.lt.360d0))then
          dir2=dir2+360d0
        elseif((0d0.le.dirNode).and.(dirNode.le.dir2))then
          dir1=dir1-360d0
        endif
      endif

      if((dir1.le.dirNode).and.(dirNode.le.dir2))then
        if(mod(is,2).eq.0)then
          secWei=0d0
          secWei(is/2)=1d0
        else
          secWei=0d0
          secWei(secInt(is,1))=1d0 - 1d0/(dir2-dir1)*(dirNode-dir1)
          secWei(secInt(is,2))=0d0 + 1d0/(dir2-dir1)*(dirNode-dir1)
        endif
      endif      
    enddo

    !! Garratt formula is used for all sectors till certain wind speed
    secCd=max(0.001d0*(0.75d0 + 0.067d0*w),cdLimit)

    !!Sector 1
    if(secCd(1).gt.0.002d0)then
      if(w.le.35d0)then
        secCd(1)=0.002d0
      elseif(w.le.45d0)then
        secCd(1)=0.0020d0 + (0.003d0-0.002d0)/(45d0-35d0)*(w-35d0)
      else
        secCd(1)=0.003d0
      endif
    endif

    !!Sector 2
    if(secCd(2).gt.0.002d0)then
      if(w.le.35d0)then
         secCd(2)=0.002d0
      elseif(w.le.45d0)then
        secCd(2)=0.0020d0 + (0.001d0-0.002d0)/(45d0-35d0)*(w-35d0)
      else
        secCd(2)=0.001d0
      endif
    endif

    !!Sector 3
    if(secCd(3).gt.0.0018d0)then
      if(w.le.25d0)then
        secCd(3)=0.0018d0
      elseif(w.le.30d0)then
        secCd(3)=0.0018d0 + (0.0045d0-0.0018d0)/(30d0-25d0)*(w-25d0)
      elseif(w.le.45d0)then
        secCd(3)=0.0045d0 + (0.0010d0-0.0045d0)/(45d0-30d0)*(w-30d0)
      else
        secCd(3)=0.001d0
      endif
    endif

    cd=secCd(1)*secWei(1) + secCd(2)*secWei(2) + secCd(3)*secWei(3)

    sn=(coory(i)-windY0)/dr
    cs=(coorx(i)-windX0)/dr

    windTx(i)=-rho*cd*w*w*sn
    windTy(i)=rho*cd*w*w*cs
    if(dr.gt.10d0*rMax)then
      windTx(i)=0d0
      windTy(i)=0d0
    endif
    
    pr(i)=pc+pDrop*dexp(-windA/tmpr2)

  enddo

end subroutine windNew3

subroutine windNew4(windDragForm,dt,npt,lon,lat,coorx,coory,&
  windLon,windLat,windR0,coriF,windWm,windA,windB,pr,&
  windTx,windTy)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt,windDragForm
  integer(kind=C_K1)::is,secInt(6,2)
  real(kind=c_K2),intent(in)::coorx(npt),coory(npt),dt
  real(kind=c_K2),intent(in)::lon(npt),lat(npt)
  real(kind=c_K2),intent(in)::windLon(3),windLat(3)
  real(kind=c_K2),intent(in)::windR0(3),windWm(3)  
  real(kind=c_K2),intent(in)::windA,windB,coriF(npt)
  real(kind=c_K2),intent(out)::windTx(npt),windTy(npt),pr(npt)
  real(kind=c_K2)::windX0,windY0,windX0t1,windY0t1
  real(kind=c_K2)::stormSpX,stormSpY
  real(kind=c_K2)::rho,cd,rhoA,pDrop,pn,pc,dr,w,wx,wy
  real(kind=c_K2)::csl,snl,sn,cs,pmax,rMax
  real(kind=c_K2)::dir1,dir2,dirStorm,dirNode
  real(kind=c_K2),parameter::cdLimit=0.0035d0 !!Garratt Limit
  real(kind=c_K2)::secLim(6,2),secWei(3),secCd(3)

  secLim(1,:)=(/ 0d0, 40d0  /)
  secLim(2,:)=(/ 40d0, 130d0  /)
  secLim(3,:)=(/ 130d0, 170d0  /)
  secLim(4,:)=(/ 170d0, 220d0  /)
  secLim(5,:)=(/ 220d0, 260d0  /)
  secLim(6,:)=(/ 260d0, 0d0  /)

  secInt=0
  secInt(1,:)=(/ 3,1 /)
  secInt(3,:)=(/ 1,2 /)
  secInt(5,:)=(/ 2,3 /)

  rhoA=1.225d0
  rho=rhoA !! I think this is rho_air/rho_water  
  pn=1.0135*10**5
  
  pDrop=rhoA*dexp(1d0)*windWm(1)**2/windB
  pc=pn-pDrop
  rMax=windA**(1d0/windB)
  !write(*,*)windA,windB,pDrop,windWm

  tmpr1=windLon(1)
  tmpr2=windLat(1)
  call utmGeo(tmpr1,tmpr2,windX0,windY0,1,UTMZone)
  tmpr1=windLon(2)
  tmpr2=windLat(2)
  call utmGeo(tmpr1,tmpr2,windX0t1,windY0t1,1,UTMZone)
  stormSpX=(windX0-windX0t1)/dt
  stormSpY=(windY0-windY0t1)/dt

  dir1=atan2((windLat(2)-windLat(3)),(windLon(2)-windLon(3)))
  dir2=atan2((windLat(1)-windLat(2)),(windLon(1)-windLon(2)))
  dirStorm=atan((sin(dir1)+sin(dir2)),(cos(dir1)+cos(dir2)))
  dirStorm=dirStorm/deg2rad
  dirStorm=90d0-dirStorm
  if(dirStorm.lt.0d0) dirStorm=dirStorm+360d0  

  windTx=0d0
  windTy=0d0
  pr=0d0  

  tmpr1=windA*windB*pDrop
  do i=1,npt
    dr=dsqrt((coorx(i)-windX0)**2 + (coory(i)-windY0)**2)
    if(dr.eq.0d0)then 
      pr(i)=pc
      cycle
    endif     

    tmpr2=dr**windB
    w=tmpr1*dexp(-windA/tmpr2)/(rhoA*tmpr2)
    w=w+(dr**2)*coriF(i)**2/4d0
    w=dsqrt(w)-(dr*coriF(i)/2d0)

    sn=(coory(i)-windY0)/dr
    cs=(coorx(i)-windX0)/dr
    wx=-w*sn+stormSpX
    wy=w*cs+stormSpY
    w=dsqrt(wx**2 + wy**2)

    if(windDragForm.eq.1)then
      !! Garratt Formulation
      cd=max(0.001d0*(0.75d0 + 0.067d0*w),cdLimit)

    elseif(windDragForm.eq.2)then
      !! Powell Formulation
      dirNode=atan2((lat(i)-windLat(1)),(lon(i)-windLon(1)))
      dirNode=dirNode/deg2rad
      dirNode=90d0-dirNode
      if(dirNode.lt.0d0) dirNode=dirNode+360d0

      !! Compute weight for each of the three sectors
      secWei=0d0
      do is=1,6
        dir1=mod(dirStorm+secLim(is,1),360d0)
        dir2=mod(dirStorm+secLim(is,2),360d0)

        if(dir1.gt.dir2)then
          if((dir1.le.dirNode).and.(dirNode.lt.360d0))then
            dir2=dir2+360d0
          elseif((0d0.le.dirNode).and.(dirNode.le.dir2))then
            dir1=dir1-360d0
          endif
        endif

        if((dir1.le.dirNode).and.(dirNode.le.dir2))then
          if(mod(is,2).eq.0)then
            secWei=0d0
            secWei(is/2)=1d0
          else
            secWei=0d0
            secWei(secInt(is,1))=1d0 - 1d0/(dir2-dir1)*(dirNode-dir1)
            secWei(secInt(is,2))=0d0 + 1d0/(dir2-dir1)*(dirNode-dir1)
          endif
        endif      
      enddo

      !! Garratt formula is used for all sectors till certain wind speed
      secCd=max(0.001d0*(0.75d0 + 0.067d0*w),cdLimit)

      !!Sector 1
      if(secCd(1).gt.0.002d0)then
        if(w.le.35d0)then
          secCd(1)=0.002d0
        elseif(w.le.45d0)then
          secCd(1)=0.0020d0 + (0.003d0-0.002d0)/(45d0-35d0)*(w-35d0)
        else
          secCd(1)=0.003d0
        endif
      endif

      !!Sector 2
      if(secCd(2).gt.0.002d0)then
        if(w.le.35d0)then
           secCd(2)=0.002d0
        elseif(w.le.45d0)then
          secCd(2)=0.0020d0 + (0.001d0-0.002d0)/(45d0-35d0)*(w-35d0)
        else
          secCd(2)=0.001d0
        endif
      endif

      !!Sector 3
      if(secCd(3).gt.0.0018d0)then
        if(w.le.25d0)then
          secCd(3)=0.0018d0
        elseif(w.le.30d0)then
          secCd(3)=0.0018d0 + (0.0045d0-0.0018d0)/(30d0-25d0)*(w-25d0)
        elseif(w.le.45d0)then
          secCd(3)=0.0045d0 + (0.0010d0-0.0045d0)/(45d0-30d0)*(w-30d0)
        else
          secCd(3)=0.001d0
        endif
      endif

      cd=secCd(1)*secWei(1) + secCd(2)*secWei(2) + secCd(3)*secWei(3)
    
    else
      !! No matching Formulation
      write(tf,'(" [ERR] No matching windDrag formulation")')
      stop
    endif

    windTx(i)=rho*cd*w*wx
    windTy(i)=rho*cd*w*wy
    if(dr.gt.10d0*rMax)then
      windTx(i)=0d0
      windTy(i)=0d0
    endif
    
    pr(i)=pc+pDrop*dexp(-windA/tmpr2)

  enddo

end subroutine windNew4

subroutine windNew4b(windDragForm,RedFacW,RedFacP,dt,npt,lon,lat,&
  coorx,coory,windLon,windLat,windR0,coriF,windWm,&
  windA,windB,pr,windTx,windTy)
use basicVars
implicit none

  integer(kind=C_K1),intent(in)::npt,windDragForm
  integer(kind=C_K1)::is,secInt(6,2)
  real(kind=c_K2),intent(in)::coorx(npt),coory(npt),dt
  real(kind=c_K2),intent(in)::lon(npt),lat(npt)
  real(kind=c_K2),intent(in)::windLon(3),windLat(3)
  real(kind=c_K2),intent(in)::windR0(3),windWm(3)  
  real(kind=c_K2),intent(in)::windA,windB,coriF(npt)
  real(kind=c_K2),intent(out)::windTx(npt),windTy(npt),pr(npt)
  real(kind=c_K2)::windX0,windY0,windX0t1,windY0t1
  real(kind=c_K2)::stormSpX,stormSpY
  real(kind=c_K2)::rho,cd,rhoA,pDrop,pn,pc,dr,w,wx,wy
  real(kind=c_K2)::csl,snl,sn,cs,pmax,rMax
  real(kind=c_K2)::dir1,dir2,dirStorm,dirNode
  real(kind=c_K2),parameter::cdLimit=0.0035d0 !!Garratt Limit
  real(kind=c_K2)::secLim(6,2),secWei(3),secCd(3)

  !! Holland formula gives geostropic winds. 
  !! The Drag coeff in for U_10 which is 10m above sea level
  !! As per Delf3D Manual a PReduce is multiplied to geostropic wind vel
  !! Pressure drop is divided by  PReduce**2
  real(kind=c_K2),intent(in)::RedFacW,RedFacP

  !! Check Delf3D Manual Section 5.2. 
  !! Wind are rotated by 20deg to account for frictional effects \_**_/
  !! This makes it look like spiralling towards the centre
  real(kind=c_K2),parameter::rotWind=20*deg2rad  
  csl=dcos(rotWind)
  snl=dsin(rotWind)

  secLim(1,:)=(/ 0d0, 40d0  /)
  secLim(2,:)=(/ 40d0, 130d0  /)
  secLim(3,:)=(/ 130d0, 170d0  /)
  secLim(4,:)=(/ 170d0, 220d0  /)
  secLim(5,:)=(/ 220d0, 260d0  /)
  secLim(6,:)=(/ 260d0, 0d0  /)

  secInt=0
  secInt(1,:)=(/ 3,1 /)
  secInt(3,:)=(/ 1,2 /)
  secInt(5,:)=(/ 2,3 /)

  rhoA=1.225d0
  rho=rhoA !! I think this is rho_air/rho_water  
  pn=1.0135*10**5
  
  pDrop=rhoA*dexp(1d0)*windWm(1)**2/windB/(RedFacP**2)
  pc=pn-pDrop
  rMax=windA**(1d0/windB)
  write(tf,'(" [CYC] ",A15,2F15.6)')"pDrop rMax :",pDrop,rMax

  tmpr1=windLon(1)
  tmpr2=windLat(1)
  call utmGeo(tmpr1,tmpr2,windX0,windY0,1,UTMZone)
  tmpr1=windLon(2)
  tmpr2=windLat(2)
  call utmGeo(tmpr1,tmpr2,windX0t1,windY0t1,1,UTMZone)
  stormSpX=(windX0-windX0t1)/dt
  stormSpY=(windY0-windY0t1)/dt

  dir1=atan2((windLat(2)-windLat(3)),(windLon(2)-windLon(3)))
  dir2=atan2((windLat(1)-windLat(2)),(windLon(1)-windLon(2)))
  dirStorm=atan((sin(dir1)+sin(dir2)),(cos(dir1)+cos(dir2)))
  dirStorm=dirStorm/deg2rad
  dirStorm=90d0-dirStorm
  if(dirStorm.lt.0d0) dirStorm=dirStorm+360d0  

  windTx=0d0
  windTy=0d0
  pr=0d0  

  tmpr1=windA*windB*pDrop

  !$OMP PARALLEL DEFAULT(shared) &
  !$OMP   PRIVATE(i,dr,tmpr2,w,sn,cs,wx,wy,cd,dirNode,&    
  !$OMP     secWei,is,dir1,dir2,secCd)
  !$OMP DO SCHEDULE(dynamic,100)
  do i=1,npt
    dr=dsqrt((coorx(i)-windX0)**2 + (coory(i)-windY0)**2)
    if(dr.eq.0d0)then 
      pr(i)=pc
      cycle
    endif     

    tmpr2=dr**windB
    w=tmpr1*dexp(-windA/tmpr2)/(rhoA*tmpr2)
    w=w+(dr**2)*coriF(i)**2/4d0
    w=(dsqrt(w)-(dr*coriF(i)/2d0))*RedFacW

    sn=(coory(i)-windY0)/dr
    cs=(coorx(i)-windX0)/dr
    wx=-w*sn+stormSpX
    wy=w*cs+stormSpY
    w=dsqrt(wx**2 + wy**2)

    if(windDragForm.eq.1)then
      !! Garratt Formulation
      cd=max(0.001d0*(0.75d0 + 0.067d0*w),cdLimit)

    elseif(windDragForm.eq.2)then
      !! Powell Formulation
      dirNode=atan2((lat(i)-windLat(1)),(lon(i)-windLon(1)))
      dirNode=dirNode/deg2rad
      dirNode=90d0-dirNode
      if(dirNode.lt.0d0) dirNode=dirNode+360d0

      !! Compute weight for each of the three sectors
      secWei=0d0
      do is=1,6
        dir1=mod(dirStorm+secLim(is,1),360d0)
        dir2=mod(dirStorm+secLim(is,2),360d0)

        if(dir1.gt.dir2)then
          if((dir1.le.dirNode).and.(dirNode.lt.360d0))then
            dir2=dir2+360d0
          elseif((0d0.le.dirNode).and.(dirNode.le.dir2))then
            dir1=dir1-360d0
          endif
        endif

        if((dir1.le.dirNode).and.(dirNode.le.dir2))then
          if(mod(is,2).eq.0)then
            secWei=0d0
            secWei(is/2)=1d0
          else
            secWei=0d0
            secWei(secInt(is,1))=1d0 - 1d0/(dir2-dir1)*(dirNode-dir1)
            secWei(secInt(is,2))=0d0 + 1d0/(dir2-dir1)*(dirNode-dir1)
          endif
        endif      
      enddo

      !! Garratt formula is used for all sectors till certain wind speed
      secCd=max(0.001d0*(0.75d0 + 0.067d0*w),cdLimit)

      !!Sector 1
      if(secCd(1).gt.0.002d0)then
        if(w.le.35d0)then
          secCd(1)=0.002d0
        elseif(w.le.45d0)then
          secCd(1)=0.0020d0 + (0.003d0-0.002d0)/(45d0-35d0)*(w-35d0)
        else
          secCd(1)=0.003d0
        endif
      endif

      !!Sector 2
      if(secCd(2).gt.0.002d0)then
        if(w.le.35d0)then
           secCd(2)=0.002d0
        elseif(w.le.45d0)then
          secCd(2)=0.0020d0 + (0.001d0-0.002d0)/(45d0-35d0)*(w-35d0)
        else
          secCd(2)=0.001d0
        endif
      endif

      !!Sector 3
      if(secCd(3).gt.0.0018d0)then
        if(w.le.25d0)then
          secCd(3)=0.0018d0
        elseif(w.le.30d0)then
          secCd(3)=0.0018d0 + (0.0045d0-0.0018d0)/(30d0-25d0)*(w-25d0)
        elseif(w.le.45d0)then
          secCd(3)=0.0045d0 + (0.0010d0-0.0045d0)/(45d0-30d0)*(w-30d0)
        else
          secCd(3)=0.001d0
        endif
      endif

      cd=secCd(1)*secWei(1) + secCd(2)*secWei(2) + secCd(3)*secWei(3)
    
    else
      !! No matching Formulation
      write(tf,'(" [ERR] No matching windDrag formulation")')
      stop
    endif

    windTx(i)=rho*cd*w*(wx*csl-wy*snl)
    windTy(i)=rho*cd*w*(wy*csl+wx*snl)
    if(dr.gt.10d0*rMax)then
      windTx(i)=0d0
      windTy(i)=0d0
    endif
    
    pr(i)=pc+pDrop*dexp(-windA/tmpr2)

  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

end subroutine windNew4b


subroutine cycTrack(cycN,cycP,cycInf,rtime,windLon,windLat,windR0,windWm)
use basicVars
implicit none
  
  integer(kind=C_K1),intent(in)::cycN
  integer(kind=C_K1),intent(inout)::cycP
  real(kind=C_K2),intent(in)::cycInf(cycN,5),rtime
  real(kind=C_K2),intent(out)::windLon,windLat,windR0,windWm

  if(cycInf(cycP+1,1).lt.rtime) cycP=cycP+1
  i=cycP
  j=cycP+1
  tmpr1=cycInf(j,1)-cycInf(i,1)
  tmpr2=rtime-cycInf(i,1)
  tmpr1=tmpr2/tmpr1
  
  k=2
  windLon=cycInf(i,k)+(cycInf(j,k)-cycInf(i,k))*tmpr1
  k=3
  windLat=cycInf(i,k)+(cycInf(j,k)-cycInf(i,k))*tmpr1
  k=4
  windR0=cycInf(i,k)+(cycInf(j,k)-cycInf(i,k))*tmpr1
  k=5
  windWm=cycInf(i,k)+(cycInf(j,k)-cycInf(i,k))*tmpr1
  
end subroutine cycTrack