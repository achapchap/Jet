
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine initial( &

!  entry variables
 krest,kvfall,pro0,ampp,pfraz,escr,estr,kffb,npamx,nfil,nbody,ybody,zbody,vfall0, &
 epsgg, &
!  variables on exit
 vfall,proat,npc,npsl,npf,kget,mget,mgeti,npt,yv,zv,yce,zce,amp,amplim, &
 tmy,tmz,rny,rnz,phi,dphi,kphi,cdold,ddold,t,jt,llf,phid,dphid,tbody, &
 ybodys2,zbodys2,ybodyn,zbodyn,kord,ksep,kmed,ksup,ng,ampli,epsg,ngo1)

 use variables
 implicit none

 integer      :: kvfall,krest,npamx,npc,npsl,npf,kffb,mget,mgeti, &
   kget,npt,k2dtax,nbody,jt,llf,ng,ngo1

 integer, dimension(npamx)   :: kphi

 real*8       :: pro0,ampp,pfraz,escr,estr,amplim,t,vfall0,vfall,proat, &
                 ampli, epsg, epsgg                  
 
 real*8, dimension(npamx)    :: yv,zv,yce,zce,tmy,tmz,rny,rnz,amp, &
                                phi,dphi,ybody,zbody,ybodyn,zbodyn,tbody,ybodys2,zbodys2,phid,dphid

! Dold Filter variables

 integer      :: nfil
 real*8, dimension(nfil)           :: ddold
 real*8, dimension(1:nfil,0:nfil)  :: cdold

! Local variables

 integer      :: i,kint
 logical      :: iint
 real*8       :: yp1,ypn,dy,dz,dl,tai,am,amy,amz,dt,dyi,dzi,rm1,rq1, &
                 rm2,rq2,tlo,yin,zin,ylo,zlo,ampl,thetai,thetaf,thel, &
                 den,den2

! Flow Separation Variables
integer, dimension(npamx)      :: kord, ksep
integer :: kmed,ksup

! First executable statement

! initialization of the Dold filter variables 

 call filiniz(nfil,cdold,ddold)

! initialization of some generally used variables

 amplim = ampp/pfraz  ! va inizializzata anche in restart

 if (krest.eq.0) then

!   initialization of time integration variables at the beginning

   t      = 0.d0 ! time
   jt     = 0    ! step of the time integration
   llf    = 0    ! output configuration number

!   initialization of the spline representation for the body contour
!   starting from the geometry in input
!   NOTE: be sure that the dataset of points in geo.in is large
!   enough to assure an accurate reconstruction

   tbody(1) = 0.d0
   do i=2,nbody
     dy = ybody(i) - ybody(i-1)
     dz = zbody(i) - zbody(i-1)
     dl = sqrt(dy**2 + dz**2)
     tbody(i) = tbody(i-1) + dl
   enddo

   yp1 = 1.e31
   ypn = 1.e31
   call spline(ybody,zbody,ybodys2,zbodys2,tbody,yp1,ypn,nbody,npamx)

!   move the body to the actual position (zbodyn(1)=pro0). Note that
!   the second derivatives are independent of the vertical translation

   proat=pro0
   do i =1,nbody
     zbodyn(i) = zbody(i)+proat
     ybodyn(i) = ybody(i)
   end do

!   define the wetted portion of the body by intersecting the still
!   liquid surface with the body contour. By using the notation
!   z=m*y+q, the equation of the still liquid surface is m=0, q=pro0

   rm1 = 0.d0
   rq1 = 0.d0 
   i  = 0
   iint = .false.    ! becomes 1 when intersection is found
   do while (.not.iint.and.i.lt.nbody)
     i = i+1
     if (abs(ybodyn(i)-ybodyn(i+1)).gt.1.d-8) then
       rm2 = (zbodyn(i)-zbodyn(i+1))/(ybodyn(i)-ybodyn(i+1))
       rq2 = zbodyn(i)-rm2*ybodyn(i)
       if (abs(rm1-rm2).gt.1.d-8) then
         yin = (rq2-rq1)/(rm1-rm2)
         zin = rm1*yin+rq1
         iint = .true.
       else
         stop 'Intersection not found: Parallel lines '
       end if
     else   ! vertical line
       yin = ybodyn(i)
       zin = rm1*yin+rq1
       iint = .true.
     end if
   enddo
  
   if (.not.iint) stop 'Intersection not found within nbody '
  
   kint = i   ! is the last body node before the intersection
   dyi  = yin-ybodyn(kint)
   dzi  = zin-zbodyn(kint)
  
!   curvilinear abscissa at the intersection point (needed for
!   spline discretization)

   tai  = tbody(kint) + sqrt( dyi**2 + dzi**2 )

!   discretization of the wetted portion of the body contour
!   at the initial time step a minimum of 6 panels is used

   npc = 6
   dt  = tai/float(npc)
   do i = 1,npc
     tlo = float(i-1)/float(npc)*tai
     call splint(tlo,ylo,zlo,ybodyn,zbodyn,ybodys2,zbodys2,tbody,nbody,npamx)
     yv(i)  = ylo
     zv(i)  = zlo 
   enddo        
   yv(npc+1) = yin
   zv(npc+1) = zin 

!   discretization of the free surface

   npsl  = log(1.d0+(estr-yin)/ampp*(escr-1.d0))/log(escr) + 1
   ampl = (estr-yin)*(1.d0-escr)/(1.d0-escr**npsl)

   do i = npc+2,npc+npsl+1
     yv(i) = yv(i-1) + ampl
     zv(i) = 0.d0
     ampl  = ampl*escr
   enddo

! just to use on the far-field boundary the same amplitude of the last 
! free surface panel

   ampl = ampl/escr

!   initialization of the far field boundary: it is located at
!   r=estr for \theta = -90,0 (the symmetry condition about y=0 is 
!   exploited
!   NOTE: at the moment the 2D dipole solution is employed, and it
!   is not valid for the axisymmetric flow

   if (kffb.eq.1) then  
     if (k2dtax.eq.1) then
       write(*,*) 'Far field behaviour not implemented for &
       axisymmetric flow'
       stop
     end if 

!    the panel size on the far field boundary is uniform and it is
!    equal to that of the last free surface panel 
          
     thetai = atan2(zv(npc+npsl+1),yv(npc+npsl+1))
     thetaf = -pi/2.d0
     npf    =  estr*(thetai-thetaf)/ampl

     do i = npc+npsl+2,npc+npsl+npf+1
       thel  = thetai+(thetaf-thetai)*real(i-npc-npsl-1)/real(npf)
       yv(i) = estr*cos(thel)
       zv(i) = estr*sin(thel)
     enddo
   end if 

   npt  = npc + npsl + npf

! At the beginning of the simulation there is no jet, nor separated free
! surface

   kget = 0
   mget = 0
   mgeti = 0

   if (npt.gt.npamx) stop 'Maximum panel number exceeded '

   if (kvfall.eq.0) vfall = vfall0

!  initialization of panel data and boundary conditions

   do i = 1,npc
     yce(i)  = (yv(i+1)+yv(i))/2.d0
     zce(i)  = (zv(i+1)+zv(i))/2.d0
     amy     =  yv(i+1) - yv(i)
     amz     =  zv(i+1) - zv(i)
     am      = sqrt(amy*amy+amz*amz)
     amp(i)  = am
     tmy(i)  = amy/am     
     tmz(i)  = amz/am     
     rny(i)  =  tmz(i)
     rnz(i)  = -tmy(i)
     dphi(i) = -vfall*rnz(i)
     kphi(i) = 0
   enddo
   do i = npc+1,npc+npsl
     yce(i)  = (yv(i+1)+yv(i))/2.d0
     zce(i)  = (zv(i+1)+zv(i))/2.d0
     amy     =  yv(i+1) - yv(i)
     amz     =  zv(i+1) - zv(i)
     am      = sqrt(amy*amy+amz*amz)
     amp(i)  = am
     tmy(i)  = amy/am     
     tmz(i)  = amz/am     
     rny(i)  =  tmz(i)
     rnz(i)  = -tmy(i)
     phi(i)  = 0.d0
     kphi(i) = 1
   enddo
  
   if (kffb.eq.1) then

     do i = npc+npsl+1,npc+npsl+npf
       yce(i)  = (yv(i+1)+yv(i))/2.d0
       zce(i)  = (zv(i+1)+zv(i))/2.d0
       amy     =  yv(i+1) - yv(i)
       amz     =  zv(i+1) - zv(i)
       am      = sqrt(amy*amy+amz*amz)
       amp(i)  = am
       tmy(i)  = amy/am     
       tmz(i)  = amz/am     
       rny(i)  =  tmz(i)
       rnz(i)  = -tmy(i)
       phid(i) = zce(i)/(yce(i)**2+zce(i)**2)
       den     = (yce(i)**2+zce(i)**2)
       den2    = den**2
       dphid(i) = -2*zce(i)*yce(i)/den2*rny(i) +   &
                  (1.d0/den-2*zce(i)**2/den2)*rnz(i)
       kphi(i) = 2
     enddo
   end if
 else    ! Restart option
   stop ' Restart option not included yet '
 end if

 do i = 1,npc+1
   write(104,*) yv(i),zv(i)
 enddo
   
 write(104,*) 
 write(104,*) 

 do i = npc+1,npc+npsl+1
   write(104,*) yv(i),zv(i)
 enddo

 write(104,*) 
 write(104,*) 

 do i = npc+npsl+1,npt+1
   write(104,*) yv(i),zv(i)
 enddo

 close(104)
!   ---- here below some quantities (mainly related to the jet and its
!   modelling) which were in the old version and have not been initialized
!   yet.
!   To be kept until the full funcionalities are restored

 kmed = 0 
 print*, kmed
 ksup = 0 
 ngo1=nbody


 do i=1,npamx
   ksep(i)=0
   kord(i)=2
 enddo

!Jet these are probably being initialized elsewhere, (TO DO check and fix it) 
 ng   = 0
 kget = 0
 ampli = 0.d0
 epsg = epsgg*ampli 


return
end
