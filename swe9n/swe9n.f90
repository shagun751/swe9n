!!---------------------- Version 1.x.x -----------------------!!
!!  -> Nine-noded element
!!  -> Boundary - SemiDirect + Penalty + Gauss Seidel
!!    -> 11 - Inlet - No Absorbing
!!    -> 12 - NoSlip Wall 
!!    -> 13 - Slip Wall - Rectangular wall only
!!    -> 14 - Outlet  - Not coded
!!    -> 15 - Sponge - BOUSS2D approach generalised input
!!  -> Solver - Normalised
!!  -> Generalised input
!!  -> GMRES with guess value
!!  -> Drichlet BndCond for eta, p, q
!!  -> Paralution CSR
!!	-> XML output
!!------------------------------------------------------------!!


include 'subroutines/mods.f90'
include 'subroutines/modsMFree.f90'
include 'subroutines/eqnMatrices.f90'
include 'subroutines/femMatrices.f90'
include 'subroutines/geometry.f90'
include 'subroutines/misc.f90'
include 'subroutines/output.f90'
include 'subroutines/utmGeo.f90'
include 'subroutines/wind.f90'
include 'subroutines/sweMeshFree.f90'
include 'subroutines/modsInletBC.f90'

program sweNineNoded
use basicVars
use jacobianModule
!use rpimModule
use omp_lib
use meshFreeMod
use airyWaveModule
implicit none

!!-----------------------Declarations-------------------------!!
  integer(kind=C_K1)::i, j, k, l, tmpi1, tmpi2
  integer(kind=C_K1)::i2, j2, k2
  real(kind=C_K2)::tmpr1, tmpr2, tmpr3, tmpr4 
  character(len=256)::probname, text
  logical:: ex

  integer(kind=C_K1)::npl,npq,npt,nele,nbnd,nbndtyp,nnzt
  integer(kind=C_K1)::nbndpoi
  integer(kind=C_K1)::itime,ntime,fileOut
  integer(kind=C_K1)::cycN,cycP,ompNThread,depIter
  integer(kind=C_K1)::windDragForm
  integer(kind=C_K1),allocatable::conn(:,:),mabnd(:,:)
  integer(kind=C_K1),allocatable::poi2poi(:,:)
  integer(kind=C_K1),allocatable::bnd11p(:),bnd12p(:),bnd14p(:)
  integer(kind=C_K1),allocatable::npoisur(:,:)
  integer(kind=C_K1),allocatable::ivf(:),jvf(:)  
  integer(kind=C_K1),allocatable::wetpoi(:),midpoi(:)
  integer(kind=C_K1),allocatable::elejvf9x9(:,:),probe(:)

  type(jacbType),allocatable::jacb(:)
  real(kind=C_K2)::shF(9,9),shFE(9,9),shFN(9,9),shW(9)
  real(kind=C_K2)::dt,rtime,minDe,etaTWei(3)
  real(kind=C_K2)::windLon(4),windLat(4)
  real(kind=C_K2)::windR0(4),windWm(4),windA,windB,tau0
  real(kind=C_K2)::RedFacW,RedFacP
  real(kind=C_K2),allocatable::tmpra(:,:),bndpNm(:,:),bndLen(:)
  real(kind=C_K2),allocatable::coorx(:),coory(:),dep(:)
  real(kind=C_K2),allocatable::lat(:),lon(:),tau(:)
  real(kind=C_K2),allocatable::depTmp1(:),depTmp2(:)
  real(kind=C_K2),allocatable::eta(:),etat1(:),etat2(:)
  real(kind=C_K2),allocatable::etat1sq(:),etat0sq(:)
  real(kind=C_K2),allocatable::etaImp(:)
  real(kind=C_K2),allocatable::ht1(:),ht0(:)
  real(kind=C_K2),allocatable::ut1(:),vt1(:),uMagt1(:)
  real(kind=C_K2),allocatable::ut0(:),vt0(:),uMagt0(:)
  real(kind=C_K2),allocatable::p(:),pt1(:),pt2(:)
  real(kind=C_K2),allocatable::q(:),qt1(:),qt2(:)
  real(kind=C_K2),allocatable::jxt1(:),jyt1(:)
  real(kind=C_K2),allocatable::jxTilt1(:),jyTilt1(:) !! For Open BC
  real(kind=C_K2),allocatable::gM0(:),gK11(:),gK21(:),gK31(:)
  real(kind=C_K2),allocatable::gK41(:),gK51(:)  
  real(kind=C_K2),allocatable::gK61(:),gK71(:),gK81(:)
  real(kind=C_K2),allocatable::gD11(:),gE12(:),gD13(:)
  real(kind=C_K2),allocatable::gD14(:),gD15(:),gD16(:)
  real(kind=C_K2),allocatable::gD21(:),gE22(:),gE23(:)
  real(kind=C_K2),allocatable::gD24(:),gD25(:),gE26(:)
  real(kind=C_K2),allocatable::gD34(:),gD35(:)
  real(kind=C_K2),allocatable::gE81(:),gD88(:)
  real(kind=C_K2),allocatable::gE41(:),gE42(:),gE43(:)
  real(kind=C_K2),allocatable::gD44(:),gD54(:),gEBot(:)  

  real(kind=C_K2),allocatable::gR1(:),gR2(:),gR3(:),gR4(:)
  real(kind=C_K2),allocatable::gR5(:),gR6(:),gR7(:)
  real(kind=C_K2),allocatable::spngC(:),coriF(:)
  real(kind=C_K2),allocatable::windTx(:),windTy(:)
  real(kind=C_K2),allocatable::cycInf(:,:),pr(:),eleArea(:)
  real(kind=C_K2),allocatable::gTxx(:),gTxy(:),gTyx(:),gTyy(:)
  real(kind=C_K2),allocatable::bnTx(:),bnTy(:), bnObc(:)
  real(kind=C_K2),allocatable::probeLoc(:,:),cour(:)
  real(kind=C_K2)::lR1,lR2,lR3,lR4,lR5,lR6,lR7

  type(mfPoiTyp),allocatable::pObj(:)
  type(airyType)::wvIn

  logical::meshLonLat
!!---------------------End Declarations-----------------------!!



  call getarg(1,probname)  
  if(len_trim(probname).lt.1) then
    write(*,'(A)',advance='no')"Enter Problem Name:"
    read(*,*)probname
  endif
  write(*,*)"Problem Name: "//probname(1:len_trim(probname))    

  !! Terminal output file name
  text=probname(1:len_trim(probname))//'.rout'
  open(tf,file=text(1:len_trim(text)))
  write(tf,*)"Problem Name: "//probname(1:len_trim(probname))  

  call system_clock(sysClk(1),k)
  sysRate=real(k,8)

!!------------------------Input File--------------------------!!
  !! Input file name
  text=probname(1:len_trim(probname))//'.inp'
  inquire(file=text(1:len_trim(text)),exist=ex)
  if(ex) then
    open(ifl(1),file=text(1:len_trim(text)))    
  else
    write(*,*)"[ERR] Missing input file"
    stop
  endif

  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)dt
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)tmpr1
  ntime=floor(tmpr1/dt,4)  
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)tmpr1
  fileOut=int(tmpr1/dt,4)  
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)meshLonLat
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)tau0
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)etaTWei
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)minDe
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)depIter
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)ompNThread  
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)i
  read(ifl(1),*,end=21,err=21)text
  if(i.gt.0)then
    allocate(probe(0:i),probeLoc(i,2))
    probe(0)=i
    do j=1,probe(0)
      read(ifl(1),*,end=21,err=21)probeLoc(j,:)
    enddo
  else
    allocate(probe(0:0),probeLoc(0:0,2))
    probe(0)=0
  endif  

  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)text
  read(ifl(1),*,end=21,err=21)tmpr1, tmpr2, tmpr3
  wvIn=airyType(tmpr1, tmpr3, tmpr2, 0d0, 0d0, 0d0)
  write(tf,'(" [INF] ",3A15)')'T','L','d'
  write(tf,'(" [---] ",3F15.6)')wvIn%T,wvIn%L,wvIn%d
  write(tf,'(" [INF] ",A15)')'d/L'
  write(tf,'(" [---] ",F15.6)')wvIn%d/wvIn%L
  goto 22

  21 write(tf,*) "[ERR] Check input file format"
  stop
  22 write(tf,*) "[MSG] Input read done"
  close(ifl(1))

  !! Cyclone track information  
  text=probname(1:len_trim(probname))//'.ctk'
  inquire(file=text(1:len_trim(text)),exist=ex)
  if(ex) then
    open(ifl(4),file=text(1:len_trim(text)))    
  else
    write(*,*)"[ERR] Missing cyclone track file"
    stop
  endif

  read(ifl(4),*,end=23,err=23)text  
  read(ifl(4),*,end=23,err=23)windA,windB
  read(ifl(4),*,end=23,err=23)text  
  read(ifl(4),*,end=23,err=23)RedFacW,RedFacP
  read(ifl(4),*,end=23,err=23)text  
  read(ifl(4),*,end=23,err=23)windDragForm
  read(ifl(4),*,end=23,err=23)text  
  read(ifl(4),*,end=23,err=23)cycN
  allocate(cycInf(cycN,5))
  write(tf,'(" [MSG] Ensure that ",I10," records are available")')cycN
  read(ifl(4),*,end=23,err=23)text  
  do i=1,cycN
    read(ifl(4),*,end=23,err=23)cycInf(i,1:5)    
  enddo
  goto 24

  23 write(tf,*) "[ERR] Check cyclone track file format"
  stop
  24 write(tf,*) "[MSG] Cyclone track read done"
  !! Debug comments      
  do i=1,cycN
    write(tf,'(5F15.4)')cycInf(i,1:5)
  enddo
  close(ifl(4))
  tmpr1=windA**(1/windB)
  write(tf,'(" [INF] Cyclone constants : A B R(m)(at max Vel)")')
  write(tf,'(3F15.6)')windA,windB,tmpr1
  if(cycInf(1,1).ne.0d0)then
    write(tf,*)"[Err] Cyclone track should start at time 0"
    stop
  endif
  tmpr1=ntime*dt
  if(cycInf(cycN,1).lt.tmpr1)then
    write(tf,*)"[Err] Simulation is running for ",tmpr1," secs"
    write(tf,*)"[Err] Cyclone track available for ",&
      cycInf(cycN,1)," secs"
    stop
  endif  
  write(tf,'(" [INF] Cyclone WindDrag Formulation")')
  if(windDragForm.eq.1)then
    write(tf,'(" [---] Garratt Symmetric")')
  elseif(windDragForm.eq.2)then
    write(tf,'(" [---] Powell Asymmetric")')
  else
    write(tf,'(" [ERR] No matching formulation")')
    stop
  endif
  write(tf,*)
  cycP=1  
!!----------------------End Input File------------------------!!

!!-------------------------Mesh File--------------------------!!
  !! Mesh file name
  text=probname(1:len_trim(probname))//'.msh'
  inquire(file=text(1:len_trim(text)),exist=ex)
  if(ex) then
    open(ifl(2),file=text(1:len_trim(text)))    
  else
    write(*,*)"[ERR] Missing mesh file"
    stop
  endif

  !!--------------------------Mesh Type---------------------!!  
  read(ifl(2),*,end=11,err=11)text
  read(ifl(2),*,end=11,err=11)nele,npl
  read(ifl(2),*,end=11,err=11)text
  read(ifl(2),*,end=11,err=11)nbnd,nbndtyp
  write(tf,'(a7,a)')"[INF] ","Elemets, Linear Nodes"
  write(tf,'(a7,2I15)')"[---] ",nele,npl
  write(tf,'(a7,a)')"[INF] ","Bnd BndTyp"
  write(tf,'(a7,2I15)')"[---] ",nbnd,nbndtyp

  !! [Note]: Euler Poincare formula to find number of midnodes
  !! Ele+Nodes-Edges=chi
  !! Giving extra npl to store the quad nodes initially
  npq=2*(nele+npl)
  npt=npl+npq

  !! Nodes read & UTM conversion
  allocate(coorx(npt),coory(npt))
  read(ifl(2),*,end=11,err=11)text  
  do i=1,npl
    read(ifl(2),*,end=11,err=11)tmpr1,tmpr2    
    call utmGeo(tmpr1,tmpr2,tmpr3,tmpr4,1,UTMZone)
    coorx(i)=tmpr3
    coory(i)=tmpr4
    if(.not.meshLonLat)then
      coorx(i)=tmpr1
      coory(i)=tmpr2
    endif
  enddo  
  write(tf,*)"[MSG] Done nodes read"

  !! Elements read
  allocate(conn(nele,9))
  conn=0
  read(ifl(2),*,end=11,err=11)text
  do i=1,nele
    ![Imp Edit] : if conn is node anticlk in msh file    
    read(ifl(2),*,end=11,err=11)(conn(i,j),j=1,4)
    ![Imp Edit] : if conn is node clk in msh file
    !read(mafi(1),*,end=11,err=11)tempstr(1),tempstr(3),tempstr(2)    
  enddo
  write(tf,*)"[MSG] Done elements read"

  !! Boundary read
  !! mabnd(:,1:6)=(/ n1 n2 ele type quadNode sideNum /)
  allocate(mabnd(nbnd,6))
  mabnd=0
  k=0;
  do l=1,nbndtyp
    read(ifl(2),*,end=11,err=11)text
    read(ifl(2),*,end=11,err=11)tmpi1,tmpi2
    do i=k+1,k+tmpi2
      read(ifl(2),*,end=11,err=11)mabnd(i,1:3)      
      mabnd(i,4)=tmpi1
    enddo
    k=k+tmpi2
  enddo
  if(k.ne.nbnd) goto 14
  write(tf,*)"[MSG] Done boundaries read" 

  !! Depth read
  allocate(dep(npt))
  dep=-999
  read(ifl(2),*,end=11,err=11)text
  do i=1,npl
    read(ifl(2),*,end=11,err=11)dep(i)
  enddo
  write(tf,*)"[MSG] Done depth read"

  ! !! Debug comments
  ! write(tf,*)"[DBG] Nodes"
  ! do i=1,npl
  !   write(tf,'(I10,3F20.6)')i,coorx(i),coory(i),dep(i)
  ! enddo
  ! write(tf,*)"[DBG] Elements"
  ! do i=1,nele
  !   write(tf,*)i,conn(i,1:4)
  ! enddo
  ! write(tf,*)"[DBG] Boundaries"
  ! do i=1,nbnd
  !   write(tf,*)i,mabnd(i,1:4)
  ! enddo

  goto 12
  !!-----------------------End Mesh Type 4------------------!!
   
  11 write(tf,*) "[ERR] Check mesh file format"
  stop
  13 write(tf,*) "[ERR] hex2dec error"
  stop
  14 write(tf,*) "[ERR] Number of boundaries mismatch"
  stop
  12 close(ifl(2))
  write(tf,*) "[MSG] Successful mesh file read"
!!-----------------------End Mesh File------------------------!!


!!--------------------Generate Middle Points------------------!!
  call middleNodes(npl,npq,npt,nele,coorx,coory,dep,conn)
  call bndMiddleNodes(nele,nbnd,conn,mabnd)

  allocate(tmpra(npt,3))
  tmpra(:,1)=coorx(1:npt)
  tmpra(:,2)=coory(1:npt)
  tmpra(:,3)=dep(1:npt)
  deallocate(coorx,coory,dep)
  allocate(coorx(npt),coory(npt),dep(npt))
  coorx=tmpra(:,1)
  coory=tmpra(:,2)
  dep=tmpra(:,3)
  deallocate(tmpra)

  !! Conver UTM to LonLat for lin+quad and Coriolis Coeff
  allocate(lon(npt),lat(npt),coriF(npt))
  do i=1,npt
    call utmGeo(lon(i),lat(i),coorx(i),coory(i),2,UTMZone)
    coriF(i)=2d0*7.29212d-5*dsin(deg2rad*lat(i))    
  enddo
  if(.not.meshLonLat)then
    lon=coorx
    lat=coory
    coriF=0d0    
  endif

  nbndpoi=2*nbnd !!Lin + quad bnd nodes
  allocate(bnd11p(0:nbndpoi),bnd12p(0:nbndpoi),bnd14p(0:nbndpoi))
  allocate(bndLen(nbnd),bndpNm(npt,2))
  call bndNormal(npt,nbnd,nbndpoi,mabnd,coorx,coory,&
    bnd11p,bnd12p,bnd14p,bndLen,bndpNm)  

  ! !! Debug comments
  ! write(tf,*)"[DBG] Nodes all"
  ! do i=1,npt
  !   write(tf,'(I10,3F20.6)')i,coorx(i),coory(i),dep(i)
  ! enddo  

  write(tf,*)"[INF] LinNode, QuadNode, TotNode"
  write(tf,*)"[---] ",npl,npq,npt
  write(tf,*)"[INF] NumEle, NumBnd, NumBntTyp"
  write(tf,*)"[---] ",nele,nbnd,nbndtyp
  write(tf,*)"[INF] Time step(s), NumTStep(n), FileOut(n)"
  write(tf,'(" [---] ",F15.6,I10,I10)')dt,ntime,fileOut
  write(tf,*)"[INF] Eta Time Weight k+1, k, k-1"
  write(tf,'(" [---] ",3F15.6)')etaTWei

!!------------------End Generate Middle Points----------------!!

!!-----------------------Node Connectivity--------------------!!
  allocate(poi2poi(npt,maxNePoi),npoisur(npt,3))
  call nodeConn(npl,npt,nele,conn,poi2poi,npoisur)
  do i=1,npt
    call mergeSort(maxNePoi,npoisur(i,1),poi2poi(i,:))
  enddo

  !! CSR matrices
  allocate(ivf(npt+1))
  nnzt=0
  ivf(1)=1
  do i=1,npt
    nnzt=nnzt+npoisur(i,1)+1
    ivf(i+1)=ivf(i)+npoisur(i,1)+1
  enddo

  i2=0
  allocate(jvf(nnzt))
  do i=1,npt
    k=i2+1
    k2=i2+npoisur(i,1)
    jvf(k:k2)=poi2poi(i,1:npoisur(i,1))
    jvf(k2+1)=i
    i2=k2+1
  enddo
  if((ivf(npt+1).ne.nnzt+1).or.(i2.ne.nnzt))then
    write(tf,*)"[ERR] Check CSR matrices"
    write(tf,*)"[---] nnzt, i2, ivf(npt+1)"
    write(tf,*)"[---]",nnzt,i2,ivf(npt+1)
  endif

  ! !! Debug comments
  ! write(tf,*)'[DBG] Node connectivity'
  ! do i=1,npt
  !   write(tf,*)"Node",i,":",npoisur(i,1)
  !   write(tf,*)"----",poi2poi(i,1:npoisur(i,1))
  ! enddo
  ! write(tf,*)"[DBG] ivf"
  ! write(tf,*)ivf
  ! write(tf,*)"[DBG] jvf"
  ! do i=1,npt
  !   k=ivf(i)
  !   k2=ivf(i+1)-1
  !   write(tf,*)'Node',i,':',jvf(k:k2)
  ! enddo

  write(tf,*)'[INF] Number of nnzt = ',nnzt
  write(tf,*)'[MSG] Done node connectivity'  
  deallocate(poi2poi,npoisur)
!!---------------------End Node Connectivity------------------!!

!!------------------------Allocations-------------------------!!  
  allocate(jacb(nele),elejvf9x9(nele,81))
  allocate(tau(npt),eleArea(nele))
  allocate(eta(npt),etat1(npt),etat2(npt))
  allocate(etat1sq(npt),etat0sq(npt),ht1(npt),ht0(npt))
  allocate(etaImp(npt))
  allocate(ut1(npt),vt1(npt),uMagt1(npt))
  allocate(ut0(npt),vt0(npt),uMagt0(npt))
  allocate(p(npt),pt1(npt),pt2(npt))
  allocate(q(npt),qt1(npt),qt2(npt))
  allocate(jxt1(npt),jyt1(npt))
  allocate(jxTilt1(npt), jyTilt1(npt))
  allocate(gM0(npt),gK11(npt),gK21(npt),gK31(npt))  
  allocate(gK41(npt),gK51(npt))
  allocate(gK61(npt),gK71(npt),gK81(npt))
  allocate(gD11(nnzt),gE12(npt),gD13(nnzt))
  allocate(gD14(nnzt),gD15(nnzt),gD16(nnzt))
  allocate(gD21(nnzt),gE22(npt),gE23(npt))
  allocate(gD24(nnzt),gD25(nnzt),gE26(npt))
  allocate(gD34(nnzt),gD35(nnzt),gEBot(npt))  
  allocate(gE41(npt),gE42(npt),gE43(npt))
  allocate(gD44(nnzt),gD54(nnzt),gE81(npt),gD88(nnzt))  
  allocate(gR1(npt),gR2(npt),gR3(npt),gR4(npt))  
  allocate(gR5(npt),gR6(npt),gR7(npt))
  allocate(gTxx(npt),gTxy(npt),gTyx(npt),gTyy(npt))
  allocate(bnTx(npt),bnTy(npt), bnObc(npt))
  allocate(wetpoi(npt),midpoi(npt))  
  allocate(spngC(npt),cour(npt))
  allocate(windTx(npt),windTy(npt),pr(npt))    
  !allocate(midObj(npt))
  allocate(pObj(npt))
!!----------------------End Allocations-----------------------!!

!!----------------------Initialisation------------------------!!
  call omp_set_num_threads(ompNThread)

  rtime=0d0
  eta=0d0  
  p=0d0
  q=0d0
  etat1=0d0
  pt1=0d0
  qt1=0d0
  wetpoi=1  
  spngC=1  

  !! Changing the water depth to minDe if depth less than minDe
  k=0
  do i=1,npt
    if(dep(i).lt.minDe) then
      dep(i)=minDe
      k=k+1
    endif
  enddo
  write(tf,'(" [INF] Changed depth for ",I10," nodes")')k  

  !! Smoothing water depth
  write(tf,*)
  write(tf,'(" [MSG] Starting Depth Smoothing")')
  allocate(depTmp1(npt),depTmp2(npt))
  depTmp2=dep
  do j=1,depIter
    depTmp1=0d0
    do i=1,npt
      k=ivf(i)
      k2=ivf(i+1)-2
      depTmp1(i)=sum(dep(jvf(k:k2)))/(k2-k+1)
    enddo  
    tmpr1=sqrt(sum((depTmp1-dep)**2)/npt)
    dep=depTmp1
    write(tf,'(" [---]",I10,F15.6)')j,tmpr1
  enddo
  depTmp1=(dep-depTmp2)/depTmp2*100  !!(New-Old)/Old

  !! Curant number calculation
  call calcCourant(npt,nele,conn,coorx,coory,dep,dt,cour)
  call out4NXMLDiff(probname,npl,npt,nele,ifl(3),0,&
    conn,lon,lat,dep,depTmp2,depTmp1,cour)
  deallocate(depTmp1,depTmp2)
  write(tf,'(" [MSG] End Depth Smoothing")')
  write(tf,*)
  write(tf,'(" [INF] Maximum Courant number",F15.6)')&
    maxval(cour)
  write(tf,*)

  !! GWCE tau value
  tau=0d0
  do i=1,npt
    if(dep(i).gt.200d0)then 
      tau(i)=0.005d0
    elseif(dep(i).lt.1d0)then
      tau(i)=1d0
    else
      tau(i)=1d0/dep(i)
    endif
  enddo
  tau=tau0*tau  

  !! Probe Locations    
  do i=1,probe(0)
    tmpr1=probeLoc(i,1)
    tmpr2=probeLoc(i,2)
    tmpr4=1d10
    k=0
    do j=1,npt
      tmpr3=(lon(j)-tmpr1)**2
      tmpr3=tmpr3+(lat(j)-tmpr2)**2
      if(tmpr3.lt.tmpr4)then
        tmpr4=tmpr3
        k=j
      endif
    enddo
    probe(i)=k
  enddo  
  write(tf,'(" [INF] Probes Number = ",I10)')probe(0)
  ! write(tf,'(" [INF] ",A30,4A,30A)')"Required Lon Lat"," || ",&
  !   "Reported Lon Lat Dep"
  write(tf,'(" [INF] ",30A)')"Required Lon Lat"
  write(tf,'(" [---] ",30A)')"Reported Lon Lat Dep"
  do i=1,probe(0)
    k=probe(i)
    ! write(tf,'(" [INF] ",2F15.6," || ",3F15.6)')probeLoc(i,:),&
    !  lon(k),lat(k),dep(k)
    write(tf,'(" [INF] ",2F15.6)')probeLoc(i,:)
    write(tf,'(" [---] ",3F15.6)')lon(k),lat(k),dep(k)
  enddo
  write(tf,*)

  ! !! Sponge layer definition
  ! tmpr1=0
  ! tmpr2=5d0
  ! tmpr3=tmpr2-tmpr1
  ! do i=1,npt
  !   if((coorx(i).gt.tmpr1).and.(coorx(i).lt.tmpr2))then
  !     tmpr4=(coorx(i)-tmpr1)/tmpr3
  !     tmpr4=3d0*tmpr4**2 - 2d0*tmpr4**3
  !     spngC(i)=tmpr4
  !   endif
  ! enddo

  ! !! Gauss hump initial condition
  ! eta=0.045d0*dexp(-2d0*((coorx-3d0)**2 + (coory-3d0)**2))
  
  ! !! Solitary wave fnc2
  ! tmpr1=1d0
  ! tmpr2=0.01d0
  ! tmpr3=dsqrt(grav*(tmpr1+tmpr2))
  ! tmpr4=dsqrt(3*tmpr2/(4*(tmpr1**3)))
  ! do i=1,npt    
  !   if((coorx(i).ge.10d0).and.(coorx(i).le.70d0)) then
  !     tmpr5=tmpr2/(dcosh(tmpr4*(coorx(i)-(40d0)))**2)
  !     u(i)=tmpr3*tmpr5/tmpr1
  !     eta(i)=tmpr5
  !   endif
  ! enddo  
!!--------------------End Initialisation----------------------!!

!!---------------------Constant Matrices----------------------!!
  call shapeFnc(shF,shFE,shFN,shW)
  call jacobianInv(npt,nele,conn,coorx,coory,&
    shF,shFE,shFN,shW,jacb,eleArea)

  !! Storing the position of element matrix in global matrix
  call eleToGlob(npt,nele,nnzt,conn,ivf,jvf,elejvf9x9)  

  call left(npt,nele,conn,jacb,shW,gM0)  

  !! Continuity Equation
  gK11=gM0+(dt/2d0)*(gM0*tau)
  gE12=-gM0+(dt/2d0)*(gM0*tau)

  !! Jx Equation
  gK21=gM0
  gE22=(gM0*tau)
  gE23=(gM0*coriF)
  gE26=(gM0/rhoW)
  gEBot=-(gM0*botCd)

  !! Jy Equation
  gK31=gM0

  !! P Equation
  gK41=gM0
  gE41=gM0
  gE42=(-dt)*(gM0*tau)
  gE43=(dt)*gM0

  !! Q Equation
  gK51=gM0  

  !! Corrector Equation
  gK61=gM0
  gK71=gM0
  gK81=gK11
  gE81=2d0*gM0

  call GWCErh1(npt,nele,nnzt,conn,jacb,shF,shFE,shFN,&
    shW,tau,dep,elejvf9x9,dt,gD11,gD13,gD14,gD15,gD16,&
    gD24,gD34,gD44,gD54,gD88)

  ! call GWCErh2(npt,nele,nnzt,conn,ivf,jvf,jacb,shF,shFE,shFN,&
  !   shW,dep,elejvf9x9,ut1,vt1,gD21)

  ! write(tf,*)"[MSG] Checking gN11"
  ! call chkMat(npt,nnzt,ivf,jvf,gN11)
  ! write(tf,*)"[MSG] Checking gN12"
  ! call chkMat(npt,nnzt,ivf,jvf,gN12)
  !write(tf,*)"[MSG] Checking gN21"
  !call chkMat(npt,nnzt,ivf,jvf,gN21)
  !write(tf,*)"[MSG] Checking gN31"
  !call chkMat(npt,nnzt,ivf,jvf,gN31)

  call sweMFree(npt, nnzt, ivf, jvf, coorx, coory, pObj)  

  do i=1,bnd12p(0)
    k=bnd12p(i)
    write(tf,'(I10,F15.6)')k,dep(k)
  enddo

  !! Dry region detection
  do i=1,npt
    if(dep(i).lt.minDe) wetpoi(i)=0
  enddo  
!!-------------------End Constant Matrices--------------------!!
  
!!-----------------------Time Stepping------------------------!!  
  ! call out4NXML(probname,npl,npt,nele,ifl(3),0,&
  !   conn,lon,lat,p,q,eta,dep,wetpoi,pr,windTx,windTy)
  call out4NXML(probname,npl,npt,nele,ifl(3),0,&
    conn,lon,lat,p,q,eta,dep,wetpoi,pr,windTx,windTy,pObj)

  !! Probe file
  text='Output/WaveProbe_'//probname(1:len_trim(probname))//'.dat'
  open(ifl(5),file=text(1:len_trim(text)))
  write(ifl(5),'(5A12)')'Time(s)','Probei','ProbeiEta',&
    'ProbeiP','ProbeiQ'

  !$acc enter data copyin(conn, jacb, shF, shFE, shFN, shW)
  !$acc enter data copyin(eleArea, elejvf9x9)

  !$acc enter data copyin(dep, eta, etat1, etat2)
  !$acc enter data copyin(etat1sq, etat0sq, etaImp, etaTWei)
  !$acc enter data copyin(p, pt1, pt2, q, qt1, qt2)
  !$acc enter data copyin(ht1, ut1, vt1, uMagt1)
  !$acc enter data copyin(ht0, ut0, vt0, uMagt0)
  !$acc enter data create(gR1, gR2, gR3, gR4, gR5, gR6, gR7)
  !$acc enter data copyin(jxt1, jyt1, jxTilt1, jyTilt1)
  !$acc enter data copyin(pObj)
  !$acc enter data copyin(wvIn)

  !$acc enter data copyin(bnd11p, bnd12p, bnd14p, bndpNm)
  !$acc enter data copyin(gE12, gE41, gE42, gE43, gE81)
  !$acc enter data copyin(gE22, gE23, gE26, gEBot)    
  !$acc enter data copyin(gD11, gD13, gD14, gD15, gD16)
  !$acc enter data copyin(gD24, gD34, gD44, gD54, gD88)
  !$acc enter data copyin(ivf, jvf)
  !$acc enter data copyin(gK11, gK21, gK31, gK41, gK51)
  !$acc enter data copyin(gK61, gK71, gK81)

  !$acc enter data copyin(windTx, windTy, pr) 
  !$acc enter data copyin(gTxx, gTxy, gTyx, gTyy)
  !$acc enter data copyin(gD21, gD25, gD35)  
  !$acc enter data copyin(bnTx, bnTy, bnObc)
  do itime=1,ntime
    call system_clock(sysClk(2))
    rTime=rTime+dt
    write(tf,'(" Time : ",I10," : ",F15.6)')itime,rTime      

    !$acc parallel loop gang vector default(present) private(i)
    do i = 1, npt
      etat2(i)=etat1(i)
      pt2(i) = pt1(i)
      qt2(i) = qt1(i)
      etat1(i) = eta(i)
      pt1(i) = p(i)
      qt1(i) = q(i)
      ht1(i) = etat1(i) + dep(i)
      ut1(i) = pt1(i) / ht1(i)
      vt1(i) = qt1(i) / ht1(i)
      etat1sq(i) = etat1(i) * etat1(i)
      uMagt1(i) = dsqrt( ut1(i)**2 + vt1(i)**2 )

      gR1(i) = 0d0
      gR2(i) = 0d0
      gR3(i) = 0d0
      gR4(i) = 0d0
      gR5(i) = 0d0    
    enddo

    !$acc update self(ht1, ut1, vt1)

    !! Cyclone Track
    if(itime.lt.3)then
      call cycTrack(cycN,cycP,cycInf,(rtime),&
        windLon(1),windLat(1),windR0(1),windWm(1))
      call cycTrack(cycN,cycP,cycInf,(rtime-dt),&
        windLon(2),windLat(2),windR0(2),windWm(2))

      do i=3,4
        windLon(i)=2d0*windLon(i-1)-windLon(i-2)
        windLat(i)=2d0*windLat(i-1)-windLat(i-2)
        windR0(i)=2d0*windR0(i-1)-windR0(i-2)
        windWm(i)=2d0*windWm(i-1)-windWm(i-2)
      enddo
    else
      do i=4,2,-1
        windLon(i)=windLon(i-1)
        windLat(i)=windLat(i-1)
        windR0(i)=windR0(i-1)
        windWm(i)=windWm(i-1)
      enddo
      call cycTrack(cycN,cycP,cycInf,(rtime),&
        windLon(1),windLat(1),windR0(1),windWm(1))      
    endif
    write(tf,'(" [CYC] ",A15,2F15.6)')"Position :",windLon(1),windLat(1)
    write(tf,'(" [CYC] ",A15,F15.6)')"Max Wind :",windWm(1)
    
    call system_clock(sysClk(6))        
    call GWCErh2(npt,nele,nnzt,conn,jacb,shF,shFE,shFN,&
      shW,eleArea,elejvf9x9,ht1,ut1,vt1,gD21,gD25,gD35,&
      gTxx,gTxy,gTyx,gTyy)                
    ! call GWCErh2ACC(npt,nele,nnzt,conn,jacb,shF,shFE,shFN,&
    !   shW,eleArea,elejvf9x9,ht1,ut1,vt1,gD21,gD25,gD35,&
    !   gTxx,gTxy,gTyx,gTyy)                
    call system_clock(sysClk(7))

    !! Cyclone wind            
    call windNew4b(windDragForm,RedFacW,RedFacP,dt,npt,lon,lat,&
      coorx,coory,windLon(2:4),windLat(2:4),windR0(2:4),coriF,&
      windWm(2:4),windA,windB,pr,windTx,windTy)    

    ! gR1=0d0
    ! gR2=0d0
    ! gR3=0d0
    ! gR4=0d0
    ! gR5=0d0    

    !$acc update device(gD21, gD25, gD35)
    !$acc update device(windTx, windTy, pr)
    !$acc update device(gTxx, gTxy, gTyx, gTyy)

    !! Solving for Jx and Jy
    !! [Note] : Here the terms are not multiplied by dt    
    !$acc parallel loop gang vector default(present) &
    !$acc   private(i, j, j2, k, k2, lR2, lR3)
    do i = 1, npt

      gR2(i) = gR2(i) + (gE22(i)*pt1(i)) + (gE23(i)*qt1(i)) &
        + (gE26(i)*windTx(i)) + (gEBot(i)*uMagt1(i)*ut1(i)) &
        + gTxx(i) + gTxy(i) + 0d0*bnTx(i)
      gR3(i) = gR3(i) + (gE22(i)*qt1(i)) - (gE23(i)*pt1(i)) &
        + (gE26(i)*windTy(i)) + (gEBot(i)*uMagt1(i)*vt1(i)) &
        + gTyx(i) + gTyy(i) + 0d0*bnTy(i)    

      k=ivf(i)
      k2=ivf(i+1)-1
      lR2=0d0
      lR3=0d0
      do j=k,k2
        j2=jvf(j)        
        lR2=lR2+(gD21(j)*pt1(j2))+(gD24(j)*etat1sq(j2)) &
          +(gD25(j)*pr(j2))        
        lR3=lR3+(gD21(j)*qt1(j2))+(gD34(j)*etat1sq(j2)) &
          +(gD35(j)*pr(j2))
      enddo      
      gR2(i)=gR2(i)+lR2
      gR3(i)=gR3(i)+lR3      
    enddo
    
    !$acc parallel loop default(present) private(i)
    do i = 1, npt
      jxt1(i) = gR2(i) / gK21(i)
      jyt1(i) = gR3(i) / gK31(i)
    enddo
    

    call obcCalcJxTil(npt, dep, jxt1, jyt1, etat1, pObj, &
      jxTilt1, jyTilt1)

    !$acc update self(jxTilt1, jyTilt1)
    !!$acc update self(ht1, ut1, vt1)

    call bndInt(npt,nele,nbnd,nnzt,conn,mabnd,jacb,shF,&
      shFE,shFN,shW,eleArea,elejvf9x9,bndLen,bndpNm, &
      ht1, ut1, vt1, jxTilt1, jyTilt1, bnTx, bnTy, bnObc)

    !$acc update device(bnTx, bnTy, bnObc)
    
    ! write(8,'(I10, F20.10)')itime, rTime    
    ! do k = 1, bnd14p(0)
    !   i = bnd14p(k)      
    !   write(8,'(I10, 3F20.10)')i, jxTilt1(i), jyTilt1(i), bnObc(i)
    ! enddo
    ! write(8,*)

    !! Solving for Eta 
    !! Predictor for P and Q
    !$acc parallel loop gang vector default(present) &
    !$acc   private(i, j, j2, k, k2, lR1, lR4, lR5)
    do i=1,npt

      gR1(i) = gR1(i) + (gE12(i)*etat2(i)) + (dt*dt*bnObc(i))
      gR4(i) = gR4(i) + (gE41(i)*pt1(i)) + (gE42(i)*pt1(i)) &
        + (gE43(i)*jxt1(i))
      gR5(i) = gR5(i) + (gE41(i)*qt1(i)) + (gE42(i)*qt1(i)) &
        + (gE43(i)*jyt1(i))         

      k=ivf(i)
      k2=ivf(i+1)-1
      lR1=0d0
      lR4=0d0
      lR5=0d0
      do j=k,k2
        j2=jvf(j)        
        lR1=lR1+(gD11(j)*etat1(j2))+(gD13(j)*pt1(j2)) &
          +(gD14(j)*qt1(j2))+(gD15(j)*jxt1(j2)) &
          +(gD16(j)*jyt1(j2))        
        lR4=lR4+(gD44(j)*etat1(j2))        
        lR5=lR5+(gD54(j)*etat1(j2))
      enddo      
      gR1(i)=gR1(i)+lR1      
      gR4(i)=gR4(i)+lR4
      gR5(i)=gR5(i)+lR5      
    enddo

    !$acc parallel loop default(present) private(i)
    do i = 1, npt
      eta(i) = gR1(i) / gK11(i)
      p(i) = gR4(i) / gK41(i)
      q(i) = gR5(i) / gK51(i)
    enddo

    !! Forcing wall BC    
    !$acc parallel loop default(present) &
    !$acc   private(i, k, tmpr1, tmpr2, tmpr3, tmpr4)
    do i=1,bnd12p(0)
      !if(wetpoi(i).eq.0)cycle
      k=bnd12p(i)
      tmpr1=bndpNm(k,1)
      tmpr2=bndpNm(k,2)
      tmpr3=p(k)*tmpr2*tmpr2 - q(k)*tmpr1*tmpr2
      tmpr4=-p(k)*tmpr1*tmpr2 + q(k)*tmpr1*tmpr1
      p(k)=tmpr3
      q(k)=tmpr4      
    enddo    

    !$acc update self(eta, etat1, p, q)

    !! Forcing Inlet and Open BC
    call inletBC(npt, nbndpoi, bnd11p, rTime, dep, bndpNm, &
      pObj, wvIn, eta, p, q)    
    call openBC2(npt, nbndpoi, bnd14p, dt, dep, bndpNm, &
      pObj, etat1, eta, p, q)    

    !$acc update device(eta, etat1, p, q)

    !! Corrector Steps
    !$acc parallel loop gang vector default(present) private(i)
    do i = 1, npt
      etat0sq(i) = eta(i)*eta(i)
      ht0(i) = dep(i) + eta(i)
      ut0(i) = p(i) / ht0(i)
      vt0(i) = q(i) / ht0(i)
      uMagt0(i) = dsqrt( ut0(i)**2 + vt0(i)**2)
      etaImp(i) = etaTWei(1)*eta(i) + etaTWei(2)*etat1(i) &
        + etaTWei(3)*etat2(i)
      gR1(i) = 0d0
      gR6(i) = 0d0
      gR7(i) = 0d0    
    enddo

    !$acc update self(ht0, ut0, vt0)

    call GWCErh2(npt,nele,nnzt,conn,jacb,shF,shFE,shFN,&
      shW,eleArea,elejvf9x9,ht0,ut0,vt0,gD21,gD25,gD35,&
      gTxx,gTxy,gTyx,gTyy)     
    ! call GWCErh2ACC(npt,nele,nnzt,conn,jacb,shF,shFE,shFN,&
    !   shW,eleArea,elejvf9x9,ht0,ut0,vt0,gD21,gD25,gD35,&
    !   gTxx,gTxy,gTyx,gTyy)     

    !!$acc update self(ht0, ut0, vt0)
    call bndInt(npt,nele,nbnd,nnzt,conn,mabnd,jacb,shF,&
      shFE,shFN,shW,eleArea,elejvf9x9,bndLen,bndpNm, &
      ht0, ut0, vt0, jxTilt1, jyTilt1, bnTx, bnTy, bnObc)    

    !! Cyclone wind        
    call windNew4b(windDragForm,RedFacW,RedFacP,dt,npt,lon,lat,&
      coorx,coory,windLon(1:3),windLat(1:3),windR0(1:3),coriF,&
      windWm(1:3),windA,windB,pr,windTx,windTy)     

    !$acc update device(bnTx, bnTy, bnObc)
    !$acc update device(gD21, gD25, gD35)
    !$acc update device(windTx, windTy, pr)
    !$acc update device(gTxx, gTxy, gTyx, gTyy)

    !! Corrector for P and Q  
    !! [Note] : Do not forget to multiply by dt here 
    !!          when taking terms from Jx Jy eqn        
    !$acc parallel loop gang vector default(present) &
    !$acc   private(i, j, j2, k, k2, lR1, lR6, lR7)
    do i=1,npt

      gR1(i) = gR1(i) + (gE81(i)*etat1(i)) + (gE12(i)*etat2(i)) &
        + (dt*dt*bnObc(i))
      gR6(i) = gR6(i) + (gE41(i)*pt1(i)) + (dt*gE23(i)*q(i)) &
       + (dt*gE26(i)*windTx(i)) + (dt*gEBot(i)*uMagt0(i)*ut0(i)) &
       + dt*( gTxx(i) + gTxy(i) + 0d0*bnTx(i) )
      gR7(i) = gR7(i) + (gE41(i)*qt1(i)) - (dt*gE23(i)*p(i)) &
        + (dt*gE26(i)*windTy(i)) + (dt*gEBot(i)*uMagt0(i)*vt0(i)) &
        + dt*( gTyx(i) + gTyy(i) + 0d0*bnTy(i) )

      k=ivf(i)
      k2=ivf(i+1)-1      
      lR1=0d0
      lR6=0d0
      lR7=0d0
      do j=k,k2
        j2=jvf(j)   
        lR1=lR1+(gD13(j)*pt1(j2)) &
          +(gD14(j)*qt1(j2))+(gD15(j)*jxt1(j2)) &
          +(gD16(j)*jyt1(j2))+(gD88(j)*etaImp(j2))
        lR6=lR6+dt*((gD21(j)*p(j2))+(gD24(j)*etat0sq(j2)) &
          +(gD25(j)*pr(j2)))+(gD44(j)*eta(j2))
        lR7=lR7+dt*((gD21(j)*q(j2))+(gD34(j)*etat0sq(j2)) &
          +(gD35(j)*pr(j2)))+(gD54(j)*eta(j2))
      enddo    
      gR1(i)=gR1(i)+lR1        
      gR6(i)=gR6(i)+lR6
      gR7(i)=gR7(i)+lR7      
    enddo    

    !$acc parallel loop default(present) private(i)
    do i = 1, npt
      eta(i) = gR1(i)/gK81(i)
      p(i) = 0.5d0*( p(i) + gR6(i)/gK61(i) )
      q(i) = 0.5d0*( q(i) + gR7(i)/gK71(i) )    
    enddo

    !! Forcing wall BC    
    !$acc parallel loop default(present) &
    !$acc   private(i, k, tmpr1, tmpr2, tmpr3, tmpr4)
    do i=1,bnd12p(0)
      !if(wetpoi(i).eq.0)cycle
      k=bnd12p(i)
      tmpr1=bndpNm(k,1)
      tmpr2=bndpNm(k,2)
      tmpr3=p(k)*tmpr2*tmpr2 - q(k)*tmpr1*tmpr2
      tmpr4=-p(k)*tmpr1*tmpr2 + q(k)*tmpr1*tmpr1
      p(k)=tmpr3
      q(k)=tmpr4      
    enddo    

    !$acc update self(eta, etat1, p, q)

    !! Forcing Inlet and Open BC
    call inletBC(npt, nbndpoi, bnd11p, rTime, dep, bndpNm, &
      pObj, wvIn, eta, p, q)
    call openBC2(npt, nbndpoi, bnd14p, dt, dep, bndpNm, &
      pObj, etat1, eta, p, q)    

        
    !! Output
    if(mod(itime,fileOut).eq.0) then                      
      call out4NXML(probname,npl,npt,nele,ifl(3),itime,&
        conn,lon,lat,p,q,eta,dep,wetpoi,pr,windTx,windTy,pObj)
    endif

    !! Probe write
    write(ifl(5),'(F15.6)',advance='no')rTime
    do i=1,probe(0)
      k=probe(i)
      write(ifl(5),'(I10,3F15.6)',advance='no')k,eta(k),p(k),q(k)
    enddo
    write(ifl(5),*)

    !$acc update device(eta, etat1, p, q)

    call system_clock(sysClk(5))
    tmpr1=(sysClk(5)-sysClk(2))/sysRate
    tmpr2=(sysClk(7)-sysClk(6))/sysRate
    write(tf,'(" [SPD] ",A15,3F15.6)')"TimStp Time :",&
      tmpr1,2d0*tmpr2,2d0*tmpr2/tmpr1*100d0
    write(tf,*)

  enddo

  !$acc exit data delete(conn, jacb, shF, shFE, shFN, shW)
  !$acc exit data delete(eleArea, elejvf9x9)

  !$acc exit data delete(dep, eta, etat1, etat2)
  !$acc exit data delete(etat1sq, etat0sq, etaImp, etaTWei)
  !$acc exit data delete(p, pt1, pt2, q, qt1, qt2)
  !$acc exit data delete(ht1, ut1, vt1, uMagt1)
  !$acc exit data delete(ht0, ut0, vt0, uMagt0)
  !$acc exit data delete(gR1, gR2, gR3, gR4, gR5, gR6, gR7)
  !$acc exit data delete(jxt1, jyt1, jxTilt1, jyTilt1)
  !$acc exit data delete(pObj)
  !$acc exit data delete(wvIn)

  !$acc exit data delete(bnd11p, bnd12p, bnd14p, bndpNm)
  !$acc exit data delete(gE12, gE41, gE42, gE43, gE81)
  !$acc exit data delete(gE22, gE23, gE26, gEBot)
  !$acc exit data delete(gD11, gD13, gD14, gD15, gD16)
  !$acc exit data delete(gD24, gD34, gD44, gD54, gD88)
  !$acc exit data delete(ivf, jvf)
  !$acc exit data delete(gK11, gK21, gK31, gK41, gK51)
  !$acc exit data delete(gK61, gK71, gK81)

  !$acc exit data delete(windTx, windTy, pr) 
  !$acc exit data delete(gTxx, gTxy, gTyx, gTyy)
  !$acc exit data delete(gD21, gD25, gD35)  
  !$acc exit data delete(bnTx, bnTy, bnObc)

  call system_clock(sysClk(3))  
  write(tf,'(" [SPD] ",A15,F15.6)')"Total Duration :",&
    (sysClk(3)-sysClk(1))/sysRate
!!---------------------End Time Stepping----------------------!!

  close(ifl(5))
  write(tf,*)'[MSG] sweNineNoded execution completed'  
  close(tf)  
end program sweNineNoded
!!-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x!!
