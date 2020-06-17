subroutine left(npt,nele,conn,jacb,shW,gM1)
use basicVars
use jacobianModule
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,conn(nele,9)
  type(jacbType),intent(in)::jacb(nele)
  real(kind=C_K2),intent(in)::shW(9)
  real(kind=C_K2),intent(out)::gM1(npt)

  integer(kind=C_K1)::iel, i, i2, na(9)

  gM1=0d0
  do iel=1,nele
    na=conn(iel,:)

    do i=1,9
      i2=na(i)
      gM1(i2)=gM1(i2)+(shW(i)*jacb(iel)%D(i))
      ! gM1(i2)=gM1(i2)+(jacb(iel)%D(i))/9d0
    enddo
  enddo

  do i=1,npt
    if(gM1(i).lt.1e-5)then
      write(tf,*)'[ERR] Mass too small for node',i    
    endif
  enddo

  ! !! Debug comments
  ! write(tf,*)'[DBG] Mass array'
  ! do i=1,npt
  !   write(tf,*)i,gM1(i)
  ! enddo
end subroutine left



subroutine calcMat(shW,D,sc,Ai,Aj,lN)
use basicVars
implicit none
  
  real(kind=C_K2),intent(in)::shW(9),D(9),sc(9)
  real(kind=C_K2),intent(in)::Ai(9,9),Aj(9,9)
  real(kind=C_K2),intent(out)::lN(9,9)

  integer(kind=C_K1)::i, i2, j
  real(kind=C_K2)::tmpr1, tmpr2

  lN=0d0

  do i2=1,9 !Looping over integration pts
    tmpr1=shW(i2)*D(i2)*sc(i2)
    do i=1,9
      do j=1,9
        tmpr2=tmpr1*Ai(i,i2)*Aj(j,i2)
        lN(i,j)=lN(i,j)+tmpr2
      enddo
    enddo
  enddo

end subroutine calcMat



subroutine calcMatT2(shW,D,sc,Ai,Aj,lN)
use basicVars
implicit none
  
  real(kind=C_K2),intent(in)::shW(9),D(9),sc(9)
  real(kind=C_K2),intent(in)::Ai(9,9),Aj(9,9)
  real(kind=C_K2),intent(out)::lN(81)

  integer(kind=C_K1)::i, i2, j, l2
  real(kind=C_K2)::tmpr1, tmpr2

  lN=0d0

  do i2=1,9 !Looping over integration pts
    tmpr1=shW(i2)*D(i2)*sc(i2)
    do i=1,9
      do j=1,9
        l2=(i-1)*9+j
        tmpr2=tmpr1*Ai(i,i2)*Aj(j,i2)
        lN(l2)=lN(l2)+tmpr2
      enddo
    enddo
  enddo

end subroutine calcMatT2


subroutine rh1(npt,nele,nnzt,conn,ivf,jvf,jacb,shF,shFE,shFN,&
  shW,dep,u,v,elejvf9x9,gN11,gN12,gN22)
use basicVars
use jacobianModule
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,nnzt,conn(nele,9)
  integer(kind=C_K1),intent(in)::ivf(npt+1),jvf(nnzt)
  integer(kind=C_K1),intent(in)::elejvf9x9(nele,81)
  type(jacbType),intent(in)::jacb(nele)
  type(jacbType)::lJacb
  real(kind=C_K2),intent(in)::shF(9,9),shFE(9,9),shFN(9,9)
  real(kind=C_K2),intent(in)::shW(9),dep(npt)
  real(kind=C_K2),intent(in)::u(npt),v(npt)
  real(kind=C_K2),intent(out)::gN11(nnzt),gN12(nnzt),gN22(nnzt)  
  real(kind=C_K2)::lN11(81),lN12(81),lN(81),lN22(81)
  real(kind=C_K2)::lDep(9),lU(9),lV(9)  
  integer(kind=C_K1)::thidx(81)  

  integer(kind=C_K1)::i, i2, j, l, l2, na(9), iel  

  gN11=0d0
  gN12=0d0
  gN22=0d0  

  !$OMP PARALLEL DEFAULT(shared) &
  !$OMP   PRIVATE(i,i2,j,l,l2,iel,na,lJacb,&
  !$OMP     lN11,lN12,lN22,lN,lDep,lU,lV,thidx)
  ! !$OMP   PARALLEL DEFAULT(private) SHARED(conn,jacb,dep,u,v,ivf,jvf,gN11,gN12,gN22)
  !$OMP DO SCHEDULE(dynamic,100) !REDUCTION(+:gIEl)

  do iel=1,nele    
    na=conn(iel,:)
    lJacb=jacb(iel)
    lDep=dep(na)
    lU=u(na)
    lV=v(na)
    thidx=elejvf9x9(iel,:)
    
    lN11=0d0
    call calcMatT2(shW,lJacb%D11,lDep,shFE,shF,lN)
    lN11=lN11+lN
    call calcMatT2(shW,lJacb%D12,lDep,shFN,shF,lN)
    lN11=lN11+lN

    lN12=0d0
    call calcMatT2(shW,lJacb%D21,lDep,shFE,shF,lN)
    lN12=lN12+lN
    call calcMatT2(shW,lJacb%D22,lDep,shFN,shF,lN)
    lN12=lN12+lN

    lN22=0d0
    call calcMatT2(shW,lJacb%D11,lU,shF,shFE,lN)
    lN22=lN22+lN
    call calcMatT2(shW,lJacb%D12,lU,shF,shFN,lN)
    lN22=lN22+lN
    call calcMatT2(shW,lJacb%D21,lV,shF,shFE,lN)
    lN22=lN22+lN
    call calcMatT2(shW,lJacb%D22,lV,shF,shFN,lN)
    lN22=lN22+lN
    

    !$OMP CRITICAL
    gN11(thidx)=gN11(thidx)+lN11
    gN12(thidx)=gN12(thidx)+lN12
    gN22(thidx)=gN22(thidx)+lN22
    !$OMP END CRITICAL
    
  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL  

  ! !! Debug comments
  ! write(tf,*)'[DBG] gN11'
  ! call printMat(npt,nnzt,ivf,jvf,gN11)
  ! write(tf,*)'[DBG] gN12'
  ! call printMat(npt,nnzt,ivf,jvf,gN12)

end subroutine rh1


subroutine rh2(npt,nele,nnzt,conn,ivf,jvf,jacb,shF,shFE,shFN,&
  shW,gN21,gN31,elejvf9x9)
use basicVars
use jacobianModule
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,nnzt,conn(nele,9)
  integer(kind=C_K1),intent(in)::ivf(npt+1),jvf(nnzt)
  integer(kind=C_K1),intent(out)::elejvf9x9(nele,81)
  type(jacbType),intent(in)::jacb(nele)
  type(jacbType)::lJacb
  real(kind=C_K2),intent(in)::shF(9,9),shFE(9,9),shFN(9,9)
  real(kind=C_K2),intent(in)::shW(9)
  real(kind=C_K2),intent(out)::gN21(nnzt),gN31(nnzt)  
  real(kind=C_K2)::lN21(9,9),lN31(9,9),lN(9,9),sc(9)

  integer(kind=C_K1)::iel, na(9), i, i2, j, j2, k, k2, l, l2

  gN21=0d0
  gN31=0d0

  sc=grav

  do iel=1,nele
    na=conn(iel,:)
    lJacb=jacb(iel)

    lN21=0d0
    call calcMat(shW,lJacb%D11,sc,shF,shFE,lN)
    lN21=lN21+lN
    call calcMat(shW,lJacb%D12,sc,shF,shFN,lN)
    lN21=lN21+lN

    lN31=0d0
    call calcMat(shW,lJacb%D21,sc,shF,shFE,lN)
    lN31=lN31+lN
    call calcMat(shW,lJacb%D22,sc,shF,shFN,lN)
    lN31=lN31+lN

    !! 9x9
    do i=1,9
      i2=na(i)
      k=ivf(i2)
      k2=ivf(i2+1)-1      
      do j=1,9
        j2=na(j)
        l2=(i-1)*9+j
        do l=k,k2
          ! if(jvf(l).eq.j2) goto 101
          if(jvf(l).eq.j2) exit
        enddo
        ! write(tf,*)'[ERR] rh1: Check negh for node',i2
        ! write(tf,*)'[---] ',jvf(k:k2)
        101 continue
        gN21(l)=gN21(l)+lN21(i,j)
        gN31(l)=gN31(l)+lN31(i,j)
        elejvf9x9(iel,l2)=l
      enddo
    enddo

  enddo

  ! !! Debug comments
  ! write(tf,*)'[DBG] gN21'
  ! call printMat(npt,nnzt,ivf,jvf,gN21)
  ! write(tf,*)'[DBG] gN31'
  ! call printMat(npt,nnzt,ivf,jvf,gN31)

end subroutine rh2

subroutine GWCErh1(npt,nele,nnzt,conn,jacb,shF,shFE,&
  shFN,shW,tau,dep,elejvf9x9,dt,gD11,gD13,gD14,gD15,gD16,&
  gD24,gD34,gD44,gD54,gD88)
use basicVars
use jacobianModule
implicit none
  
  integer(kind=C_K1),intent(in)::npt,nele,nnzt,conn(nele,9)  
  integer(kind=C_K1),intent(in)::elejvf9x9(nele,81)
  integer(kind=C_K1)::locjvf9x9(81)
  type(jacbType),intent(in)::jacb(nele)
  type(jacbType)::lJacb
  real(kind=C_K2),intent(in)::shF(9,9),shFE(9,9),shFN(9,9)
  real(kind=C_K2),intent(in)::shW(9),dep(npt),tau(npt),dt
  real(kind=C_K2),intent(out)::gD11(nnzt),gD13(nnzt)
  real(kind=C_K2),intent(out)::gD14(nnzt),gD15(nnzt),gD16(nnzt)
  real(kind=C_K2),intent(out)::gD24(nnzt),gD34(nnzt)
  real(kind=C_K2),intent(out)::gD44(nnzt),gD54(nnzt)
  real(kind=C_K2),intent(out)::gD88(nnzt)
  real(kind=C_K2)::lN(81),lN00(81),lN12(81),lN13(81),lN14(81)
  real(kind=C_K2)::lN15(81),lN16(81),lN17(81),lN24(81)
  real(kind=C_K2)::lN34(81),lN42(81),lN52(81)
  real(kind=C_K2)::shFX(9,9),shFY(9,9)
  real(kind=C_K2)::lDep(9),lSc(9),lTau(9)
  real(kind=C_K2)::lTauDx(9),lTauDy(9)

  integer(kind=C_K1)::iel, na(9), i, j, k  

  gD11=0d0
  gD13=0d0
  gD14=0d0
  gD15=0d0
  gD16=0d0  
  gD24=0d0  
  gD34=0d0  
  gD44=0d0
  gD54=0d0  
  gD88=0d0

  do iel=1,nele
    na=conn(iel,:)
    lJacb=jacb(iel)
    lDep=dep(na)
    lTau=tau(na)
    locjvf9x9=elejvf9x9(iel,:)

    lTauDx=0d0
    lTauDy=0d0
    do k=1,9
      shFX(:,k)=shFE(:,k)*lJacb%D11(k)+shFN(:,k)*lJacb%D12(k)
      shFY(:,k)=shFE(:,k)*lJacb%D21(k)+shFN(:,k)*lJacb%D22(k)

      do j=1,9
        lTauDx(k)=lTauDx(k)+(shFX(j,k)*lTau(j))
        lTauDy(k)=lTauDy(k)+(shFY(j,k)*lTau(j))
      enddo
    enddo    
    
    lSc=1d0
    call calcMatT2(shW,lJacb%D,lSc,shF,shF,lN00)    

    lSc=grav*lDep
    call calcMatT2(shW,lJacb%D,lSc,shFX,shFX,lN12)    

    lSc=grav*lDep
    call calcMatT2(shW,lJacb%D,lSc,shFY,shFY,lN13)    

    lSc=lTauDx
    call calcMatT2(shW,lJacb%D,lSc,shF,shF,lN14)        

    lSc=lTauDy
    call calcMatT2(shW,lJacb%D,lSc,shF,shF,lN15)  

    lSc=1d0
    call calcMatT2(shW,lJacb%D,lSc,shFX,shF,lN16)    

    lSc=1d0
    call calcMatT2(shW,lJacb%D,lSc,shFY,shF,lN17)    

    lSc=grav/2d0
    call calcMatT2(shW,lJacb%D,lSc,shF,shFX,lN24)        
    
    lSc=grav/2d0
    call calcMatT2(shW,lJacb%D,lSc,shF,shFY,lN34)            

    lSc=grav*lDep
    call calcMatT2(shW,lJacb%D,lSc,shF,shFX,lN42)    

    lSc=grav*lDep
    call calcMatT2(shW,lJacb%D,lSc,shF,shFY,lN52)    

    !! 9x9
    gD11(locjvf9x9)=gD11(locjvf9x9)+(2d0*lN00) &
      -(dt*dt*lN12)-(dt*dt*lN13)
    gD13(locjvf9x9)=gD13(locjvf9x9)+(dt*dt*lN14)
    gD14(locjvf9x9)=gD14(locjvf9x9)+(dt*dt*lN15)
    gD15(locjvf9x9)=gD15(locjvf9x9)+(dt*dt*lN16)
    gD16(locjvf9x9)=gD16(locjvf9x9)+(dt*dt*lN17)
    gD24(locjvf9x9)=gD24(locjvf9x9)-(lN24)    
    gD34(locjvf9x9)=gD34(locjvf9x9)-(lN34)    
    gD44(locjvf9x9)=gD44(locjvf9x9)-(dt*lN42)
    gD54(locjvf9x9)=gD54(locjvf9x9)-(dt*lN52)    
    gD88(locjvf9x9)=gD88(locjvf9x9)-(dt*dt*lN12)-(dt*dt*lN13)

  enddo

end subroutine GWCErh1

subroutine GWCErh2(npt,nele,nnzt,conn,jacb,shF,shFE,&
  shFN,shW,eleArea,elejvf9x9,ht1,ut1,vt1,gD21,gD25,gD35,&
  gTxx,gTxy,gTyx,gTyy)
use basicVars
use jacobianModule
implicit none
  
  integer(kind=C_K1),intent(in)::npt,nele,nnzt,conn(nele,9)
  integer(kind=C_K1),intent(in)::elejvf9x9(nele,81)
  integer(kind=C_K1)::locjvf9x9(81)
  type(jacbType),intent(in)::jacb(nele)
  type(jacbType)::lJacb
  real(kind=C_K2),intent(in)::shF(9,9),shFE(9,9),shFN(9,9)
  real(kind=C_K2),intent(in)::shW(9),eleArea(nele)
  real(kind=C_K2),intent(in)::ut1(npt),vt1(npt),ht1(npt)
  real(kind=C_K2),intent(out)::gD21(nnzt)
  real(kind=C_K2),intent(out)::gD25(nnzt),gD35(nnzt)
  real(kind=C_K2),intent(out)::gTxx(npt),gTxy(npt)
  real(kind=C_K2),intent(out)::gTyy(npt),gTyx(npt)
  real(kind=C_K2)::lN(81),lN21(81),lN22(81)
  real(kind=C_K2)::lN26(81),lN36(81)
  real(kind=C_K2)::shFX(9,9),shFY(9,9)
  real(kind=C_K2)::lU(9),lV(9),lH(9)
  real(kind=C_K2)::lUDx(9),lVDy(9),lSc(9)
  real(kind=C_K2)::lUDy,lVDx,lnu,lAr
  real(kind=C_K2)::lTxx(9),lTxy(9),lTyy(9),lTyx(9)

  integer(kind=C_K1)::iel, k, j, na(9)  

  gD21=0d0
  gD25=0d0
  gD35=0d0
  gTxx=0d0
  gTxy=0d0
  gTyx=0d0
  gTyy=0d0

  !$OMP PARALLEL DEFAULT(shared) &
  !$OMP   PRIVATE(iel,k,j,na,lJacb,lU,lV,lH,lSc,locjvf9x9,&
  !$OMP     lUDx,lVDy,shFX,shFY,lN,lN21,lN22,lN26,lN36,&
  !$OMP     lUDy,lVDx,lnu,lAr,lTxx,lTxy,lTyx,lTyy)
  !$OMP DO SCHEDULE(dynamic,100)

  do iel=1,nele
    na=conn(iel,:)
    lJacb=jacb(iel)    
    lAr=eleArea(iel)
    lU=ut1(na)
    lV=vt1(na)
    lH=ht1(na)
    locjvf9x9=elejvf9x9(iel,:)
    
    lUDx=0d0
    lVDy=0d0
    lTxx=0d0
    lTxy=0d0
    lTyx=0d0
    lTyy=0d0
    do k=1,9
      shFX(:,k)=shFE(:,k)*lJacb%D11(k)+shFN(:,k)*lJacb%D12(k)
      shFY(:,k)=shFE(:,k)*lJacb%D21(k)+shFN(:,k)*lJacb%D22(k)

      lUDy=0d0
      lVDx=0d0
      do j=1,9
        lUDx(k)=lUDx(k)+(shFX(j,k)*lU(j))
        lVDy(k)=lVDy(k)+(shFY(j,k)*lV(j))
        lUDy=lUDy+(shFY(j,k)*lU(j))
        lVDx=lVDx+(shFX(j,k)*lV(j))
      enddo

      lnu=smCsSq*lAr &
        *dsqrt(0.5d0*(lUDy+lVDx)**2 + lUDx(k)**2 + lVDy(k)**2)
      lTxx=lTxx+(shW(k)*shFX(:,k)*lH(k)*2d0*lnu*lUDx(k)*lJacb%D(k))
      lTxy=lTxy+(shW(k)*shFY(:,k)*lH(k)*lnu*(lUDy+lVDx)*lJacb%D(k))
      lTyx=lTyx+(shW(k)*shFX(:,k)*lH(k)*lnu*(lUDy+lVDx)*lJacb%D(k))
      lTyy=lTyy+(shW(k)*shFY(:,k)*lH(k)*2d0*lnu*lVDy(k)*lJacb%D(k))
            
    enddo    

    lN21=0d0
    call calcMatT2(shW,lJacb%D,lU,shF,shFX,lN)    
    lN21=lN21+lN
    call calcMatT2(shW,lJacb%D,lUDx,shF,shF,lN)    
    lN21=lN21+lN

    lN22=0d0
    call calcMatT2(shW,lJacb%D,lV,shF,shFY,lN)    
    lN22=lN22+lN
    call calcMatT2(shW,lJacb%D,lVDy,shF,shF,lN)    
    lN22=lN22+lN

    lSc=lH/rhoW
    call calcMatT2(shW,lJacb%D,lSc,shF,shFX,lN26)
    call calcMatT2(shW,lJacb%D,lSc,shF,shFY,lN36)

    !! 9x9
    !$OMP CRITICAL
    gD21(locjvf9x9)=gD21(locjvf9x9)-(lN21+lN22)
    gD25(locjvf9x9)=gD25(locjvf9x9)-(lN26)
    gD35(locjvf9x9)=gD35(locjvf9x9)-(lN36)
    gTxx(na)=gTxx(na)-lTxx
    gTxy(na)=gTxy(na)-lTxy
    gTyx(na)=gTyx(na)-lTyx
    gTyy(na)=gTyy(na)-lTyy
    !$OMP END CRITICAL

  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

end subroutine GWCErh2



!!-----------------------------bndInt------------------------------!!
subroutine bndInt(npt,nele,nbnd,nnzt,conn,mabnd,jacb,shF,&
  shFE,shFN,shW,eleArea,elejvf9x9,bndLen,bndpNm, &
  ht1, ut1, vt1, jxTil, jyTil, bnTx, bnTy, bnObc)
use basicVars
use jacobianModule
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,nbnd,nnzt
  integer(kind=C_K1),intent(in)::conn(nele,9),mabnd(nbnd,6)
  integer(kind=C_K1),intent(in)::elejvf9x9(nele,81)
  integer(kind=C_K1)::bndIntP(4,3),lBndIntP(3),lBndTyp
  type(jacbType),intent(in)::jacb(nele)
  type(jacbType)::lJacb
  real(kind=C_K2),intent(in)::shF(9,9),shFE(9,9),shFN(9,9)
  real(kind=C_K2),intent(in)::shW(9),eleArea(nele)
  real(kind=C_K2),intent(in)::bndpNm(npt,2),bndLen(nbnd)
  real(kind=C_K2),intent(in)::ht1(npt),ut1(npt),vt1(npt)
  real(kind=C_K2),intent(in)::jxTil(npt), jyTil(npt)
  real(kind=C_K2),intent(out)::bnTx(npt), bnTy(npt), bnObc(npt)
  real(kind=C_K2)::bndIntW(3),lIntW(3)
  real(kind=C_K2)::shFX(9),shFY(9)
  !real(kind=C_K2)::pnx(9),pny(9)
  real(kind=C_K2)::pnx,pny
  real(kind=C_K2)::lU(9),lV(9),lH(9),lAr
  real(kind=C_K2)::lnu,lUDx,lUDy,lVDx,lVDy
  real(kind=C_K2)::lbnTx(9),lbnTy(9),lbnObc(9)

  integer(kind=C_K1)::iel, na(9), l, k, k2, j

  bnTx=0d0
  bnTy=0d0
  bnObc=0d0

  bndIntW=(/  1d0/6d0,  4d0/6d0,  1d0/6d0  /)
  bndIntP(1,:)=(/  1,  5,  2 /)
  bndIntP(2,:)=(/  2,  6,  3 /)
  bndIntP(3,:)=(/  3,  7,  4 /)
  bndIntP(4,:)=(/  4,  8,  1 /)


  !$OMP PARALLEL DEFAULT(shared) &
  !$OMP   PRIVATE(l,k,k2,j,iel,na,lJacb,lIntW,lBndIntP,&
  !$OMP     lU,lV,lH,pnx,pny,lUDx,lUDy,lVDx,lVDy,&
  !$OMP     shFX,shFY,lnu,lAr,lbnTx,lbnTy,lbnObc,lBndTyp)
  !$OMP DO SCHEDULE(dynamic,100)

  do l=1,nbnd
    iel=mabnd(l,3)
    na=conn(iel,:)
    lJacb=jacb(iel)
    lAr=eleArea(iel)    
    lIntW=bndIntW*bndLen(l)
    lBndIntP=bndIntP(mabnd(l,6),:)
    lBndTyp=mabnd(l,4)

    lU=ut1(na)
    lV=vt1(na)
    lH=ht1(na)
    pnx=bndpNm(na(lBndIntP(2)),1)
    pny=bndpNm(na(lBndIntP(2)),2)
    ! pnx=bndpNm(na,1)
    ! pny=bndpNm(na,2)

    lbnTx=0d0
    lbnTy=0d0
    lbnObc=0d0
    do k2=1,3
      k=lBndIntP(k2)
      shFX=shFE(:,k)*lJacb%D11(k)+shFN(:,k)*lJacb%D12(k)
      shFY=shFE(:,k)*lJacb%D21(k)+shFN(:,k)*lJacb%D22(k)

      lUDx=sum(shFX*lU)
      lUDy=sum(shFY*lU)
      lVDx=sum(shFX*lV)
      lVDy=sum(shFY*lV)      

      lnu=smCsSq*lAr &
        *dsqrt(0.5d0*(lUDy+lVDx)**2 + lUDx**2 + lVDy**2)

      lbnTx=lbnTx &
        +(lIntW(k2)*shF(:,k)*lH(k)*2d0*lnu*lUDx*pnx) &
        +(lIntW(k2)*shF(:,k)*lH(k)*lnu*(lUDy+lVDx)*pny)
      lbnTy=lbnTy &
        +(lIntW(k2)*shF(:,k)*lH(k)*lnu*(lUDy+lVDx)*pnx) &
        +(lIntW(k2)*shF(:,k)*lH(k)*2d0*lnu*lVDy*pny)             

      ! lbnTx=lbnTx &
      !   +(lIntW(k2)*shF(:,k)*lH(k)*2d0*lnu*lUDx*pnx(k)) &
      !   +(lIntW(k2)*shF(:,k)*lH(k)*lnu*(lUDy+lVDx)*pny(k))
      ! lbnTy=lbnTy &
      !   +(lIntW(k2)*shF(:,k)*lH(k)*lnu*(lUDy+lVDx)*pnx(k)) &
      !   +(lIntW(k2)*shF(:,k)*lH(k)*2d0*lnu*lVDy*pny(k))     

      if( (lBndTyp.eq.12) .or. (lBndTyp.eq.13) )then
        lbnObc = 0d0
      else
        j=na(k)
        lbnObc=lbnObc &
          +( lIntW(k2)*shF(:,k)*( jxTil(j)*pnx + jyTil(j)*pny ) )                
      endif
      
    enddo 

    !$OMP CRITICAL
    bnTx(na)=bnTx(na)+lbnTx
    bnTy(na)=bnTy(na)+lbnTy
    bnObc(na) = bnObc(na) - lbnObc
    !$OMP END CRITICAL

  enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL  

  ! iel=7989
  ! write(111,*)
  ! do i=1,9
  !   write(111,'(I10,2F15.6)')i,bnTx(conn(iel,i)),bnTy(conn(iel,i))
  ! enddo
  ! iel=7990
  ! write(112,*)
  ! do i=1,9
  !   write(112,'(I10,2F15.6)')i,bnTx(conn(iel,i)),bnTy(conn(iel,i))
  ! enddo

end subroutine bndInt
!!---------------------------End bndInt----------------------------!!



!!---------------------------obcCalcJxTil--------------------------!!
subroutine obcCalcJxTil(npt, dep, jx, jy, eta, pObj, jxTil, jyTil)
use basicVars
use meshFreeMod
implicit none

  integer(kind=C_K1),intent(in)::npt
  real(kind=C_K2),intent(in)::dep(npt), jx(npt), jy(npt), eta(npt)
  type(mfPoiTyp),intent(in)::pObj(npt)    
  real(kind=C_K2),intent(out)::jxTil(npt), jyTil(npt)

  integer(kind=C_K1)::i, j, k, neid
  real(kind=C_K2)::etaDx, etaDy

  ! jxTil = 0d0
  ! jyTil = 0d0

  !$acc parallel loop default(present) &
  !$acc   private(k, j, neid, etaDx, etaDy)
  do k = 1, npt    

    jxTil(k) = 0d0
    jyTil(k) = 0d0
    etaDx = 0d0
    etaDy = 0d0
    do j = 1, pObj(k)%nn
      neid = pObj(k)%neid(j)
      etaDx = etaDx + pObj(k)%phiDx(j)*eta(neid)
      etaDy = etaDy + pObj(k)%phiDy(j)*eta(neid)      
    enddo
    jxTil(k) = jx(k) - grav*dep(k)*etaDx
    jyTil(k) = jy(k) - grav*dep(k)*etaDy    
  enddo


end subroutine obcCalcJxTil
!!-------------------------End obcCalcJxTil------------------------!!