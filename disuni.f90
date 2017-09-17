
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine disuni(jt,yn,zn,ycn,zcn,phin,npsl,npamx, &
                 escr,kget,estr,amplim,npc,iiget,ygb,zgb)

! TO DO: write explicit what are the the entry and exit variables 

 use variables
 implicit none
 
 
! will first try to get the code rewritten and then come to initialize
! everything correctly     

 integer :: npamx,ng
 real *8, dimension(npamx) ::  yn,zn,ycn,zcn, phin, ygb, zgb
 real *8 :: escr, estr, amplim
 integer :: jt,kget,iiget,npsl,npc
 real *8 :: tin, tfi 

! Local variables
 integer :: i,ifid,iind,ivp,kfi,kin,kt,npa ,nng,npt,ii,j,kk,kkc,nn
 real *8, dimension(npamx) ::xx,zz,pp,xs2,zs2,ps2,tp,tvi,zcb,ycb,ygf,zgf,zgfc,ygfc,phisl
 real *8 :: amma,am1,a2y,a2z,am,ama,amb,amc,amme,ammp,amms,amt,amy,amz,an
 real *8 :: anint, ax,az,det,dis,escr2,rny,rnz,sc,sca,ss,tcc,tcco,ttfi,ttin
 real *8 :: ttl, yp1, ypn,ay, phisl2,si,ti,yin,zin,ttlo,w1y,w1z,wky,wkz,y1,ygf2,ygfc2
 real *8 :: dtma, xo,yo,po,ttt, xpo,z1,zgf2,zgfc2,yiin,ziin,zk,yk
 real *8 :: zo,zpo


 kin = 0
 kfi = 0
 tcc  = 0.d0
 tcco = 0.d0

  write(*,*) 'running disuni----' 
! calculate the normal of the last body panel, normal_body

  amy =  yn(npc+1) - yn(npc) 
  amz =  zn(npc+1) - zn(npc)
  am  = sqrt(amy**2 + amz**2)
  rny =  amz/am
  rnz =  -amy/am 
        
! calculate the normal of the 1st free surface panel, normal_FS

  a2y = yn(npc+2)-yn(npc+1)
  a2z = zn(npc+2)-zn(npc+1) 
  amc = (yn(npc+2)-yn(npc+1))*rny+(zn(npc+2)-zn(npc+1))*rnz 
  amt = sqrt((yn(npc+2)-yn(npc+1))**2 +(zn(npc+2)-zn(npc+1))**2)
        
! calculate the angle between normal_body and normal_FS

  sca = (-amy*a2y-amz*a2z)/(amt*am)
  anint = acos(sca)*180.d0/pi

! WARNING !!!!! Angle at intersection is enforced

  anint=90.d0

! choose a new panel length, current this is hard coded to amt (first FS panel)

  if(anint.lt.90.d0) then
    am1 = amc
  else
    am1 = amt
  endif 

  am1 = max(am1,amplim)
  write(*,*) 'jt, am1,amplim-------------',jt,am1,amplim
     
! controllo l'angolo massimo tra due pannelli:
! se e' maggiore di un valore limite riduco l'ampiezza dei pannelli
! di conseguenza

! Warning: iiget is still to be defined in the main code  
 ! if(jt.le.iiget)then
    amma = 0.d0
    do i = npc+1,npc+npsl/2
      ss = (yn(i+1)-yn(i))*(yn(i+2)-yn(i+1)) + &
           (zn(i+1)-zn(i))*(zn(i+2)-zn(i+1))
      ammp = sqrt( (yn(i+1)-yn(i))**2 + (zn(i+1)-zn(i))**2 )
      amms =sqrt((yn(i+2)-yn(i+1))**2+(zn(i+2)-zn(i+1))**2 )
      sc = ss/(ammp*amms)
      if (abs(sc).gt.1.d0) sc = sc/abs(sc)
      an = abs(acos(sc)*180.d0/pi)
      amma = max(amma,an)
    enddo
    if (amma.gt.10.d0) am1 = am1*10.d0/amma
  !endif

! Begin spline interp , here 3 maps are going to be build

  npa   = npsl+2  ! npa is the number of points in the spline
  tp(1) = 0.d0
  xx(1) = yn(npc+1)
  zz(1) = zn(npc+1)
  ama   = sqrt( (yn(npc+2)-yn(npc+1))**2 + (zn(npc+2)-zn(npc+1))**2 )
  amb   = sqrt( (yn(npc+3)-yn(npc+2))**2 + (zn(npc+3)-zn(npc+2))**2 )

! first order expansion of the FS potential on the first FS vertice,
! I am assuming phin is a vector of panels centroids potentials, 
! indexed and following the 
! same conventions we are using for ycn and zcn

 pp(1) = phin(npc+1)+(phin(npc+1)-phin(npc+2))/(ama+amb)*ama

! Prepare the data to build the spline representation of the FS (xx,zz,pp) 
 
  do i = 1,npsl
    if (i.eq.1) then
      ax = (ycn(npc+i)-yn(npc+i))
      az = (zcn(npc+i)-zn(npc+i))
    else
      ax = (ycn(npc+i)-ycn(npc+i-1))
      az = (zcn(npc+i)-zcn(npc+i-1))
    endif
      dis     = sqrt(ax*ax+az*az)
      tp(i+1) = tp(i) + dis
      xx(i+1) = ycn(npc+i)
      zz(i+1) = zcn(npc+i)
      pp(i+1) = phin(npc+i)
  enddo
  yn(npc+npsl+1) = estr

  ax        = (yn(npc+npsl+1)-ycn(npc+npsl))
  az        = (zn(npc+npsl+1)-zcn(npc+npsl))
  dis       = sqrt(ax*ax+az*az)
  tp(npa)   = tp(npa-1) + dis
  xx(npa)   = yn(npc+npsl+1)
  zz(npa)   = zn(npc+npsl+1)
  pp(npa)   = phin(npc+npsl) ! sarebbe meglio estrapolazione lineare

! Build 3 spline maps tp -> xx, tp-> zz, tp-> phi on the FS 
! Use natural splines at the end points , is this justifiable ?

  yp1 = 1.d+31
  ypn = 1.d+31 
  call spline(xx,zz,xs2,zs2,tp,yp1,ypn,npa,npamx)
  call spline1(pp,ps2,tp,yp1,ypn,npa,npamx) 

! At this point the body vertices in the jet resgion  are available
! and stored at the variables ygb,zgb, ycb and zcb (passed on entry)
! Given them find the corresponding FS nodes/centroids  on the jet by line intersection

  kk = 1
  nn = ng
  ttlo = 0.d0
  ii = npc+1 ! Warning changed the index here

! for the panels on the body of the  jet region   
  do j = 1,ng 
     yk  = ygb(ng+1-j)
     zk  = zgb(ng+1-j)
     wky = -(zgb(ng+2-j)-zgb(ng+1-j))
     wkz =   ygb(ng+2-j)-ygb(ng+1-j)
     do while (ti.le.0.d0.and.ti.gt.1.d0.and.ii.le.(npc+npsl))   
        y1 = xx(ii)
        z1 = zz(ii)
        w1y = xx(ii+1)-xx(ii)
        w1z = zz(ii+1)-zz(ii)
! find the intersection between normal body lines and free surface panels
        call intret(yin,zin,si,ti,yk,zk,wky,wkz,y1,z1,w1y,w1z)
        ii = ii + 1
     enddo
! once the panel is found the centroids will be found by interpolation
     kk = ii
     ay = yin-xx(kk)
     az = zin - zz(kk)
     ttl = tp(kk)+sqrt(ay**2+az**2)
! small js are linear interpolated , for larger ones splint is used      
      if(j.le.4)then
         ygf2 = yin
         zgf2 = zin
         yk = ycb(ng+1-j)
         zk = zcb(ng+1-j)
         ii = npc+1
         do while (ti.lt.0.d0.and.ti.gt.1.d0.and.ii.le.(npc+npsl))
            y1 = xx(ii)
            z1 = zz(ii)
            w1y = xx(ii+1)-xx(ii)
            w1z = zz(ii+1)-zz(ii)
            call intret(yiin,ziin,si,ti,yk,zk,wky,wkz,y1,z1,w1y,w1z)
            ii = ii+1
         enddo
         kkc  = ii
         ygfc2   = yiin
         zgfc2   = ziin
         phisl2  = pp(ii)*(1.d0-ti) + ti*pp(ii+1)
         ygf(j+1)= ygf2
         zgf(j+1)= zgf2
         ygfc(j) = yiin
         zgfc(j) = ziin
         phisl(j)= phisl2
      else 
         call splint(ttl,xpo,zpo,xx,zz,xs2,zs2,tp,npa+1,npamx)
         ygf(j+1) = xpo
         zgf(j+1) = zpo
         tcc = ttlo+(ttl-ttlo)*0.5d0
         call splint3(tcc,xo,zo,po,xx,zz,pp,xs2,zs2,ps2,tp,npa+1,npamx)
         ygfc(j)  = xo
         zgfc(j)  = zo
         phisl(j) = po
      endif
      ttlo = ttl
  enddo

!load the local variables in the corresponding arrays
! Need to check/debug the index, specifically the translations by npc

  tcc =0.d0
  do j = 1,ng
     yn(npc+j+1) = ygf(j+1)
     zn(npc+j+1) = zgf(j+1)
     ycn(npc+j)  = ygfc(j)
     zcn(npc+j)  = zgfc(j)
     phin(npc+j) = phisl(j)
     if(j.eq.1)then
       ay = ycn(j+npc)-yn(1+npc)
       az = zcn(j+npc)-zn(1+npc)
     else
       ay = ycn(j+npc)-ycn(j+npc-1)
       az = zcn(j+npc)-zcn(j+npc-1)
     endif
     tcc = tcc+sqrt(ay**2+az**2)
  enddo
  
! end of new bit for now

  amme = (tp(2)+tp(3))/2.d0
  if (am1.gt.amme) then
    if (amme.lt.amplim) then
      det = min(am1,amplim)
    else
      det = amme
    end if
  else
    det = am1
  end if

  if(jt.gt.iiget)then
    escr2 = 1.d0
    det = am1
  endif

! Debugging 
!   if (jt.eq.90) then
!     do i = 1,npa
!       write(72,*) tp(i),pp(i)
!     enddo
!     stop
!   endif 
   




  ivp = 0
  kt = 1
  if (ng.eq.0) ttl = 0.d0
  
  do while(ttl.lt.tp(npa)) 
! vedo quanto manca ed in base a questo fisso il det

    dtma = tp(npa) - ttl 
!choose det based on heuristics    
    if (dtma.gt.2.d0*escr*det) then
      det  = det*escr
    else if (dtma.gt.escr*det) then
      det  = dtma/2.d0
    else
      det  = dtma
    end if 
! find the abscissa positions of the centroids and vertices
! Warning tcc is not zero anymore , need to check this!
    ivp = ivp+1
    tcc = ttl + det/2.d0
    ttt = ttl + det

! - iterpol. spline centroidi e potenziale 
    call splint3(tcc,xo,zo,po,xx,zz,pp,xs2,zs2,ps2,tp,npa,npamx)

    ycn(ivp+npc) = xo
    zcn(ivp+npc) = zo
    phin(ivp+npc) = po

! interp vertices position

    call splint(ttt,xpo,zpo,xx,zz,xs2,zs2,tp,npa,npamx)

    yn(ivp+npc+1) = xpo
    zn(ivp+npc+1) = zpo
    
    ttl = ttt
  enddo ! end of the dowhile block

  npsl= ivp
  npt = npsl + npc 

  return
  end      
