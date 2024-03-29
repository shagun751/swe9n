!! Taken from bsnq
!! Commit : a8396ef
!! Date   : 2020-05-22

!!-------------------------airyWaveModule--------------------------!!
module airyWaveModule
use basicVars
implicit none

  type, public :: airyType    
    real(kind=C_K2)::T,d,H,L,w
    real(kind=C_K2)::x0,y0
    real(kind=C_K2)::k,kx,ky
    real(kind=C_K2)::thDeg,thRad,csth,snth      
  contains
    procedure ::  getEta
    procedure ::  getPQ
  end type airyType

  interface airyType
     procedure :: waveLenCalc
  end interface airyType

contains

!!-------------------------waveLenCalc-------------------------!!
  type(airyType) function waveLenCalc(inT,inD,inH,inX0,inY0,inThDeg)
  implicit none

    !! WaveLen using dispersion relation from Airy wave-theory

    integer(kind=C_K1)::iterMax,i
    real(kind=C_K2)::l0,newl,oldl,x,errLim
    real(kind=C_K2),intent(in)::inT,inD,inH,inX0,inY0,inThDeg


    waveLenCalc%T=inT
    waveLenCalc%d=inD
    waveLenCalc%H=inH
    waveLenCalc%x0=inX0
    waveLenCalc%y0=inY0
    waveLenCalc%thDeg=inThDeg
    waveLenCalc%w=2d0*pi/waveLenCalc%T
    waveLenCalc%thRad=waveLenCalc%thDeg*pi/180d0
    waveLenCalc%csth=dcos(waveLenCalc%thRad)
    waveLenCalc%snth=dsin(waveLenCalc%thRad)

    iterMax=50000
    errLim=1d-6
    l0 = (grav/2d0/pi)*(waveLenCalc%T)**2
    oldl = l0  
    do i = 1,iterMax
      newl = l0*dtanh(2d0*pi*(waveLenCalc%d)/oldl)
      !if(mod(i,100).eq.0) write(*,*)i,oldl,newl    
      x = abs(newl-oldl)
      if (x.le.errLim) then
        waveLenCalc%L = newl
        exit
      else
        oldl = newl
      end if
    end do

    if(i.ge.iterMax) then
      write(*,*)
      write(*,*)"[ERR] waveCalculator Error waveL",waveLenCalc%L
      write(*,*)
      waveLenCalc%L=-999
      waveLenCalc%k=0d0
      return
    endif

    waveLenCalc%k=2d0*pi/waveLenCalc%L
    waveLenCalc%kx=waveLenCalc%k*waveLenCalc%csth
    waveLenCalc%ky=waveLenCalc%k*waveLenCalc%snth

  end function waveLenCalc
!!-----------------------End waveLenCalc-----------------------!!



!!---------------------------getEta----------------------------!!
  subroutine getEta(b,rTime,x,y,eta)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::eta
    real(kind=C_K2)::dx,dy

    dx=(x-b%x0)
    dy=(y-b%y0)
    eta=b%H/2d0 * dsin(b%kx*dx + b%ky*dy - b%w*rTime)

  end subroutine getEta
!!-------------------------End getEta--------------------------!!



!!----------------------------getPQ----------------------------!!
  subroutine getPQ(b,rTime,x,y,p,q)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::p,q
    real(kind=C_K2)::dx,dy

    dx=(x-b%x0)
    dy=(y-b%y0)
    p=b%H/2d0 * b%w / b%k * b%csth &
      * dsin(b%kx*dx + b%ky*dy - b%w*rTime)
    q=b%H/2d0 * b%w / b%k * b%snth &
      * dsin(b%kx*dx + b%ky*dy - b%w*rTime)

  end subroutine getPQ
!!--------------------------End getPQ--------------------------!!



!!-------------------------getEtadxdy--------------------------!!
  subroutine getEtadxdy(b,rTime,x,y,detadx,detady)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::detadx,detady
    real(kind=C_K2)::dx,dy

    dx=(x-b%x0)
    dy=(y-b%y0)
    detadx=b%H/2d0 * b%kx * dcos(b%kx*dx + b%ky*dy - b%w*rTime)
    detady=b%H/2d0 * b%ky * dcos(b%kx*dx + b%ky*dy - b%w*rTime)

  end subroutine getEtadxdy
!!-----------------------End getEtadxdy------------------------!!

end module airyWaveModule
!!-----------------------End airyWaveModule------------------------!!


!!-------------------------waveFileModule--------------------------!!
module waveFileModule
use basicVars
implicit none

  type, public :: wvFileType        
    character(len=256)::fileName
    integer(kind=C_K1)::numP,posI
    real(kind=C_K2),allocatable::data(:,:)
  contains
    procedure ::  chkWaveFileInput
    procedure ::  initWaveFile
    procedure ::  initAiryFile
    procedure ::  getEta
    procedure ::  getPQ
  end type wvFileType  

contains

!!------------------------initWaveFile-------------------------!!
  subroutine initWaveFile(b)
  implicit none

    class(wvFileType),intent(inout)::b
    integer(kind=C_K1)::mf,i
    real(kind=C_K2)::tmpra(4)
    logical::ex

    inquire(file=trim(b%fileName),exist=ex)
    if(ex) then
      open(newunit=mf,file=trim(b%fileName))
    else
      write(9,*)"[ERR] Missing wave input file"
      stop
    endif

    b%numP=0
    do while(.true.)
      read(mf,*,end=11,err=12)tmpra
      b%numP=b%numP+1
      cycle
      11 exit
      12 write(9,*)"[ERR] Check wave file format"
      stop
    enddo
    write(9,'(" [INF] Number of wave file point = ",i10)')b%numP
    close(mf)

    allocate(b%data(b%numP,4))
    open(newunit=mf,file=trim(b%fileName))    
    do i=1,b%numP
      read(mf,*,end=21,err=21)b%data(i,1:4)      
    enddo
    write(9,*)"[INF] Done wave file read"

    b%posI=2

    return
    21 write(9,*)"[ERR] Check wave file format"
    stop
  end subroutine initWaveFile
!!----------------------End initWaveFile-----------------------!!



!!------------------------initAiryFile-------------------------!!
  subroutine initAiryFile(b,dt,inTotT,inT,inD,inH,inAng)
  use airyWaveModule
  implicit none

    class(wvFileType),intent(inout)::b
    integer(kind=C_K1)::i
    real(kind=C_K2),intent(in)::inTotT,inT,inD,inH,dt,inAng
    real(kind=C_K2)::t,eta,p,q,totT
    type(airyType)::wv    

    ! airyType(T,d,H,X0,Y0,thDeg)
    wv=airyType(inT,inD,inH,0d0,0d0,inAng)

    write(9,'(" [INF] ",3A15)')'T','L','d'
    write(9,'(" [---] ",3F15.6)')wv%T,wv%L,wv%d
    write(9,'(" [INF] ",A15)')'kh'
    write(9,'(" [---] ",F15.6)')wv%k*wv%d
    write(9,'(" [INF] At 2.25 WavePeriods")')
    write(9,'(" [---] ",3A15)')'Eta','P','Q'
    call wv%getEta(2.25d0*wv%T,wv%x0,wv%y0,eta)
    call wv%getPQ(2.25d0*wv%T,wv%x0,wv%y0,p,q)
    write(9,'(" [---] ",3F15.6)')eta,p,q

    totT=1.2d0*inTotT
    b%numP=floor(totT/dt)+2

    allocate(b%data(b%numP,4))
    do i=0,b%numP-1
      t=i*dt
      call wv%getEta(t,0d0,0d0,eta)
      call wv%getPQ(t,0d0,0d0,p,q)
      b%data(i+1,1)=t
      b%data(i+1,2)=eta
      b%data(i+1,3)=p
      b%data(i+1,4)=q
      !write(201,'(4F15.6)')b%data(i+1,1:4)
    enddo    

    b%posI=2

  end subroutine initAiryFile
!!----------------------End initAiryFile-----------------------!!




!!---------------------------getEta----------------------------!!
  subroutine getEta(b,rTime,eta)
  implicit none

    class(wvFileType),intent(in)::b

    integer(kind=C_K1)::k,posI

    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::eta

    posI=2
    
    if((rTime.gt.b%data(b%numP,1)).or.&
      (rTime.lt.b%data(1,1)))then
      write(9,*)"[ERR] Wave query exceeds supplied time range"
      write(9,'( "[---] ",F15.6)')rTime
      write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
      stop
    endif

    do while(b%data(posI+1,1).le.rTime)
      posI=posI+1
    enddo

    ! do while(b%data(posI-1,1).gt.rTime)
    !   if(posI.eq.2)exit
    !   posI=posI-1
    ! enddo

    k=2
    eta=b%data(posI,k) &
      + (b%data(posI+1,k)-b%data(posI-1,k)) &
      / (b%data(posI+1,1)-b%data(posI-1,1)) &
      * (rTime-b%data(posI,1))

    !write(222,'(3F20.6)')rTime,b%data(posI,1),rTime-b%data(posI,1)

  end subroutine getEta
!!-------------------------End getEta--------------------------!!



!!----------------------------getPQ----------------------------!!
  subroutine getPQ(b,rTime,p,q)
  implicit none

    class(wvFileType),intent(in)::b

    integer(kind=C_K1)::k,posI

    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::p,q

    posI=2
    
    if((rTime.gt.b%data(b%numP,1)).or.&
      (rTime.lt.b%data(1,1)))then
      write(9,*)"[ERR] Wave query exceeds supplied time range"
      write(9,'( "[---] ",F15.6)')rTime
      write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
      stop
    endif

    do while(b%data(posI+1,1).le.rTime)
      posI=posI+1
    enddo

    ! do while(b%data(posI-1,1).gt.rTime)
    !   if(posI.eq.2)exit
    !   posI=posI-1
    ! enddo

    k=3
    p=b%data(posI,k) &
      + (b%data(posI+1,k)-b%data(posI-1,k)) &
      / (b%data(posI+1,1)-b%data(posI-1,1)) &
      * (rTime-b%data(posI,1))

    k=4
    q=b%data(posI,k) &
      + (b%data(posI+1,k)-b%data(posI-1,k)) &
      / (b%data(posI+1,1)-b%data(posI-1,1)) &
      * (rTime-b%data(posI,1))

  end subroutine getPQ
!!--------------------------End getPQ--------------------------!!



!!----------------------chkWaveFileInput-----------------------!!
  subroutine chkWaveFileInput(b,rT0,rT1,dt)
  implicit none

    class(wvFileType),intent(in)::b

    integer(kind=C_K1)::i,j,mf,nT

    real(kind=C_K2),intent(in)::rT0,rT1,dt
    real(kind=C_K2)::eta,p,q,t

    open(newunit=mf,file='chkWaveFileInput.dat')

    nT=floor((rT1-rT0)/dt)+1
    do i=0,nT
      t=rT0+i*dt
      call b%getEta(t,eta)
      call b%getPQ(t,p,q)
      write(mf,'(4F20.8)')t,eta,p,q
    enddo
    close(mf)


  end subroutine chkWaveFileInput
!!--------------------End chkWaveFileInput---------------------!!


end module waveFileModule
!!-----------------------End waveFileModule------------------------!!



!!-----------------------------inletBC-----------------------------!!
subroutine inletBC(npt, nbndpoi, bnd11p, rTime, dep, bndpNm, &
  pObj, wvIn, eta, p, q)
use meshFreeMod
use airyWaveModule
implicit none

  integer(kind=C_K1),intent(in)::npt, nbndpoi
  integer(kind=C_K1),intent(in)::bnd11p(0:nbndpoi)
  real(kind=C_K2),intent(in)::bndpNm(npt,2), dep(npt), rTime
  type(mfPoiTyp),intent(in)::pObj(npt)  
  type(airyType),intent(in)::wvIn
  real(kind=C_K2),intent(inout)::eta(npt), p(npt), q(npt)

  integer(kind=C_K1)::i, j, k, neid
  real(kind=C_K2)::wvEta, nx, ny, c, etaDn, pn

  do i = 1, bnd11p(0)
    k = bnd11p(i)
    nx = bndpNm(k,1)
    ny = bndpNm(k,2)
    call wvIn%getEta(rTime, 0d0, 0d0, wvEta)
    eta(k) = wvEta
    c = dsqrt(grav*dep(k)) 

    ! etaDn = 0d0
    ! do j = 1, pObj(k)%nn
    !   neid = pObj(k)%neid(j)
    !   etaDn = etaDn + pObj(k)%phiDx(j)*eta(neid)*nx &
    !     + pObj(k)%phiDy(j)*eta(neid)*ny
    ! enddo

    pn = c*wvEta !+ c*etaDn*pObj(k)%rad

    p(k) = -pn*nx
    q(k) = -pn*ny
  enddo  

end subroutine inletBC
!!---------------------------End inletBC---------------------------!!



!!-----------------------------openBC------------------------------!!
subroutine openBC(npt, nbndpoi, bnd14p, rTime, dep, bndpNm, &
  pObj, eta, p, q)
use meshFreeMod
use airyWaveModule
implicit none

  integer(kind=C_K1),intent(in)::npt, nbndpoi
  integer(kind=C_K1),intent(in)::bnd14p(0:nbndpoi)
  real(kind=C_K2),intent(in)::bndpNm(npt,2), dep(npt), rTime
  type(mfPoiTyp),intent(in)::pObj(npt)    
  real(kind=C_K2),intent(inout)::eta(npt), p(npt), q(npt)

  integer(kind=C_K1)::i, j, k, neid
  real(kind=C_K2)::wvEta, nx, ny, c, etaDx, etaDy, etaNei, dr, pn  
  real(kind=C_K2)::tmpr3, tmpr4

  ! integer(kind=C_K1)::test(3,2),i2,j2
  ! test(1,:)=(/ 2254, 6507 /)
  ! test(2,:)=(/ 2255, 6509 /)
  ! test(3,:)=(/ 6508, 8513 /)

  do i = 1, bnd14p(0)
    k = bnd14p(i)
    nx = bndpNm(k,1)
    ny = bndpNm(k,2)    
    c = dsqrt(grav*dep(k)) 

    etaDx = 0d0
    etaDy = 0d0
    do j = 1, pObj(k)%nn
      neid = pObj(k)%neid(j)
      etaDx = etaDx + pObj(k)%phiDx(j)*eta(neid)
      etaDy = etaDy + pObj(k)%phiDy(j)*eta(neid)
    enddo
    dr = pObj(k)%rad
    etaNei = eta(k) - etaDx*(dr*nx) - etaDy*(dr*ny)

    pn = c*etaNei    

    ! if(i.eq.9)then
    !   write(8,'(4F15.6)')rTime, pObj(k)%rad/1000d0, eta(k), etaNei
    ! endif

    tmpr3 = p(k)*ny*ny - q(k)*nx*ny + pn*nx
    tmpr4 = -p(k)*nx*ny + q(k)*nx*nx + pn*ny    
    p(k) = tmpr3
    q(k) = tmpr4

    ! do i2=1,3
    !   if(test(i2,1).eq.k)then
    !     write(8,'(2I10,2F20.10)')k,test(i2,1),eta(test(i2,2)),etaNei
    !   endif
    ! enddo    

  enddo    
  ! write(8,*)
  ! write(8,*)

end subroutine openBC
!!---------------------------End openBC----------------------------!!



!!-----------------------------openBC2-----------------------------!!
subroutine openBC2(npt, nbndpoi, bnd14p, dt, dep, bndpNm, &
  pObj, etat1, eta, p, q)
use meshFreeMod
use airyWaveModule
implicit none

  integer(kind=C_K1),intent(in)::npt, nbndpoi
  integer(kind=C_K1),intent(in)::bnd14p(0:nbndpoi)
  real(kind=C_K2),intent(in)::bndpNm(npt,2), dep(npt), dt
  real(kind=C_K2),intent(in)::etat1(npt)
  type(mfPoiTyp),intent(in)::pObj(npt)
  real(kind=C_K2),intent(inout)::eta(npt), p(npt), q(npt)

  integer(kind=C_K1)::i, j, k, neid
  real(kind=C_K2)::wvEta, nx, ny, c, etaDx, etaDy, etaEst, dr, pn  
  real(kind=C_K2)::tmpr3, tmpr4

  ! integer(kind=C_K1)::test(3,2),i2,j2
  ! test(1,:)=(/ 2254, 6507 /)
  ! test(2,:)=(/ 2255, 6509 /)
  ! test(3,:)=(/ 6508, 8513 /)

  do i = 1, bnd14p(0)
    k = bnd14p(i)
    nx = bndpNm(k,1)
    ny = bndpNm(k,2)    
    c = dsqrt(grav*dep(k)) 

    etaDx = 0d0
    etaDy = 0d0
    do j = 1, pObj(k)%nn
      neid = pObj(k)%neid(j)
      etaDx = etaDx + pObj(k)%phiDx(j)*etat1(neid)
      etaDy = etaDy + pObj(k)%phiDy(j)*etat1(neid)
    enddo
    dr = c*dt
    etaEst = etat1(k) - etaDx*(dr*nx) - etaDy*(dr*ny)
    pn = c*etaEst    

    ! if(i.eq.9)then
    !   write(8,'(4F15.6)')rTime, pObj(k)%rad/1000d0, eta(k), etaEst
    ! endif

    eta(k) = etaEst
    tmpr3 = p(k)*ny*ny - q(k)*nx*ny + pn*nx
    tmpr4 = -p(k)*nx*ny + q(k)*nx*nx + pn*ny    
    p(k) = tmpr3
    q(k) = tmpr4

    ! do i2=1,3
    !   if(test(i2,1).eq.k)then
    !     write(8,'(2I10,2F20.10)')k,test(i2,1),eta(test(i2,2)),etaEst
    !   endif
    ! enddo    

  enddo    
  ! write(8,*)
  ! write(8,*)

end subroutine openBC2
!!---------------------------End openBC2---------------------------!!
