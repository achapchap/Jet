program main
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ++                                                 ++
!     +  Fully-nonlinear BEM solver for the water entry   +
!     +    of a 2D/Axisymmetric, QUASI-arbitrary shaped   +
!     +                      body                         +
!     +   It is suitable for both, pure water entry or    +
!     +                 2D+t simulations                  +
!     ++                                                 ++
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

! definition of variables and arrays

 use variables
 implicit none

! npamx = max number of panels
! nfil  = max order of the Dold filter

 integer        :: npamx,nfil,iint
 parameter (npamx=5000, nfil=7)

! scon, sve, spre, spot: first two letters of the output files for
! configurations, velocity, pressure and potential

 character(len=2)   :: scon,svel,spre,spot

! yv,zv; yn,zn: coordinates of the panels vertices at the two 
!               Runge-Kutta levels
! yce,zce; ycn,zcn: centroids of the panels at the two 
!               Runge-Kutta levels
! tmy,tmz: components of the unit tangent vector
! rny,rnz: components of the unit normal vector (ORIENTED INWARDS) 
! amp:     panel length
! phi:     velocity potential
! dphi:    normal derivative of the velocity potential
! dpht:    tangential derivative of the velocity potential
! vym1,vzm1;vym2,vzm2: velocity components of the panel centroids 
!                at the two Runge-Kutta levels
! ybody,zbody;ybodyn,zbodyn: coordinate of the body geometry (2D contour in the
!                2D+t case) at the original and actual positions
! phid: dipole solution (unit amplitude) for 2D asymptotic behaviour
! dphid: dipole normal derivative  (unit amplitude) for 2D asymptotic behaviour
! tbody,ybodys2,zbodys2: variables for spline reconstruction of the body
!               contour
! ycb,zcb,ygb,zgb : centroids/vertices of the jet region on the body
! tgb, tcb : respective abscissas of ygb,zgb and ycb,zcb
! tc : is the centroid abscissa used by ridis6 (twice) and splver2.f in the old code
!      its meaning is still to be clarified. in particular whats the diff between tc and tcb ? 
! tn : seems to be related to tc but its a vertex abcissa 


 real*8, dimension(npamx) :: yv,zv,yn,zn,yce,zce,ycn,zcn, &
         tmy,tmz,rny,rnz,amp,phi,phin,dphi,dpht,vym1,vym2,vzm1, &
         vzm2,ybody,zbody,ybodyn,zbodyn,tbody,ybodys2,zbodys2,phid, dphid,&
         ycb,zcb,ygb,zgb,tgb,tcb,tc,tn 
! vfall0: vertical prescribed velocity (it can be either constant or variable
!        depending on kvfall)
! vfall: actual entry velocity (derived by vfrall0 and kvfall)
! pro0: initial submergence of the apex of the body
! ampp: initial FS panel size near the body contour
! pfraz: ratio ampp/smallest panel size
! ancut: angle of cut of the jet
! escr: growth factor of panel size (bulk region)
! estr: domain size (maximum x-coordinate on the right)
! tend: simulation time
! frdt: CFL limit (generally 0.25) 
! amplim: minimum panel amplitude for discretization
! dt: time step
! proat: actual submergence of the apex of the body
! gfrac : fraction of the jet kept when jet_angle<cut_angle (gfrac <= 1)
! rmg : number of panels in the jet base, given as input TO BE CONFIRMED
! espgg : inferior limit of free surface panels close to the body (amplitude fraction)
! eskg : stretching factor of the jet panels
! ampli: is the length of body where the jet model is acting upon (TO BE CONFIRMED) 
! frtend : maximum allowable time step as a fraction of the total simulation time
! flux : real variable to check the mass conservation inthe system
! coeffdi: time dependent coefficient of the dipole solutio
! eskkk:  stretching pannelli dal vertice (TO BE READ FROM INPUT!)     
! frint: fraction of the transition region (usually 0.1, TO BE READ FROM INPUT) 
! frint: matching region (0=no matching, frint<0.9))
! ramii  = multiple of amii, griglia bordi punto separazione 
! ramiii  = multiple of amii, griglia vertice corpo 
! TODO: Add them to the input file and check wha is the variable that is doing the role of amii here !






 real*8            ::  vfall0,vfall,pro0,ampp,pfraz,ancut,escr,estr, &
                       tend,frdt,amplim,t,dt,proat,ang,di,yc, gfrac,rmg, &
                       epsgg, eskg, ampli,frtend,flux,coeffdi, epsg, eskkk, &
                       frint, ramii, ramiii 

! The following arrays are used in the Jet modelling part
! TODO: describe each array and its content
 real*8, dimension(-npamx:npamx) :: xigs,zegs,xigb,zegb,xigf,zegf



! npc: number of panels on the body contour 
! npsl: number of panel on the free surface
! npf:  number of panel along the far field boundary
! nbody: number of points of the body geometry (2D contour in case of
!      2D+t)
! krest: restart option (0: start from undisturbed solution, 1:
!       start from a previous solution given as a "restart" file)
! kffb: far field asymp solution (0=no, 1=yes), not yet implemented for k2dtax=1
! k2dtax: 2D+t option (0: vertical 2D water entry, 1: vertical
!      Axisymmetric solution; 2: 2D+t solution)
! kvfall: type of entry velocity (0: constant entry velocity, 1:
!    variable entry velocity). In the 2D+t case, 0 can be used for
!    steady planing and constant inclination of the keel 
!    (e.g. wedge shaped), otherwise 1 is needed
! ksta: number of time steps between two successive outputs
! kget: kget=1 means jet modelling activated
! mget: number of elements in the modelled portion of the jet
! mgeti: number of elements in the intermediate region
! ift: time steps for the action of the Dold filter (typically 4)
! iford: order of the Dold filter (typically 3)
! llf: index of the output files
! jt : index time step
! jjget: minimum number of iterations needed before start the jet modelling
! ng    : number of panels in the jet model ,it is init here but defined and set at  shallo.f90 
! nngo : number of panels in the jet model when it is activated (ng*kget)
! ngo1 : nbody + number of separated points (as incremented in splver2) of the body spline representation 
! nsep : 0 or 1 is a separation flag
! kord: boolean 1/2 related to flow separation not really sure what is its role
! ksep: boolean 0/1 identify which panels are separated 
! nn1old: last value of nn1 saved recursively by ridis6   
! nnold:  "    "        nn       " " 
 integer         :: npc,npsl,npf,mget,mgeti,nbody,npt,krest,k2dtax, &
                    kvfall,ksta,kffb,kget,ift,iford,llf,jt,iiget
 integer         :: jjget,ng,nngo,ngo1,nsep, nn1old, nnold

! kphi: index 0 for Neumann BC, 1 for Dirichlet, 2 for far field
!     boundary, 3 panel belong to the modelled regionfor specific use (e.g. jet modelling)

 integer, dimension(npamx)     :: kphi

! Dold filter variables

 real*8, dimension(nfil)            :: ddold
 real*8, dimension(1:nfil,0:nfil)   :: cdold


! Variables related to the flow separation
integer, dimension(npamx)      :: kord, ksep
integer :: kmed,ksup






! Local Variables

 integer                       :: i,ip,mm, npc_entry,npc_exit, nng

 real*8                        :: deph2

 real*8, dimension(npamx)      :: depn1,depn2 

  
 real*8 ::am,amy,amz

! First executable statement

! Input of data for time integration and body shape

 call input (krest,k2dtax,vfall0,kvfall,pro0,ampp,pfraz,ancut,  &
    escr,estr,kffb,tend,frdt,ksta,scon,svel,spot,spre,ift,iford, &
    npamx,nbody,ybody,zbody,jjget,gfrac,rmg,epsgg,eskg,frtend, eskkk,frint,&
    ramii, ramiii)
 
! Initialization of the domain boundary and boundary conditions

 call initial (krest,kvfall,pro0,ampp,pfraz,escr,estr,kffb,npamx, &
    nfil,nbody,ybody,zbody,vfall0,epsgg,vfall,proat,npc,npsl,npf,kget,mget,&
    mgeti,npt,yv,zv,yce,zce,amp,amplim,tmy,tmz,rny,rnz,phi,dphi,kphi,&
    cdold,ddold,t,jt,llf,phid,dphid,tbody,ybodys2,zbodys2,ybodyn,zbodyn,&
    kord,ksep,kmed, ksup,ng,ampli,epsg,ngo1)
 
 ! initialize nng
 nng = int(kget*ng)


 if (krest.eq.0) then

! Solution of the BVP
       
 !call dipole(yv,zv,npc,npsl,npf,kffb,kphi,estr,phid,dphid,npamx) 
 call solver(k2dtax,kffb,yv,zv,amp,kphi,phid,dphid,npc,npsl,npf, &
                 kget,mget,mgeti,npamx,phi,dphi,coeffdi)
write(*,*) 'coeffdi------=',coeffdi 
!stop   
 
! Compute the tangential velocity

   call caveta(yce,zce,amp,phi,npc,npsl,npf,kget,mget,mgeti,npamx,dpht)

! Compute the velocity components along the body boundary and free
! surface

   do i = 1,npc+npsl
     vym1(i) = dpht(i)*tmy(i)+dphi(i)*rny(i)
     vzm1(i) = dpht(i)*tmz(i)+dphi(i)*rnz(i)
   enddo

! ----------------- File 0
  call output(yv,zv,yce,zce,phi,dphi,vym1,vzm1,npc,npsl,npf,scon,svel, &
              spot,spre,llf,npamx,t)

 else

   write(*,*) ' Restart option not implemented yet '
   stop

 end if

 do while (t.le.tend)

   jt = jt+1

! Assign the body velocity: 
! NOTE!!: here appropriate call to a routine for the computation of
! the vertical velocity should be included in the more general case.
! The parameter kvfall will be used to set the time history
! At this stage, the velocity is just set to the initial and constant
! value. 

   vfall = vfall0   ! to be changed in case of 2D+t approximation on in
                   ! case of dynamic impact

! Call to assign the time step for the next iteration
!
! In the first iterations vfall is increased smoothly towards vfall0 
! NOTE: formally the same criterion should be used at t=0 in initial.
! Here we use this condition only here as this was done in the old code
! 2D+t. This is not a problem for the solution as before the jet
! development it is not reliable anyhow.

   if(jt.le.30) vfall = vfall0*float(jt)/float(30)   
  
   call check(yv,zv,yce,zce,amp,vym1,vzm1,frdt,frtend,tend,t,jt,npc,npsl, &
              npf,npamx,kget,vfall,zbodyn,nbody,dt)


   t  = t+dt

   write(*,*) 'Simulation Summary------- t, dt, jt ', t,dt,jt

! First Runge-Kutta step: the following part needs to be updated to
! include the free surface portion detached from the body
!
! Displacement of the centroids of free surface panels and update the 
! velocity potential on them, predictor RK step
   do ip = npc+1,npc+npsl
     deph2      = vym1(ip)**2 + vzm1(ip)**2
     depn1(ip)  = deph2/2.d0   ! here the gravity term should be included
     ycn(ip)    = yce(ip) + vym1(ip)*dt
     zcn(ip)    = zce(ip) + vzm1(ip)*dt
     phin(ip)   = phi(ip) + depn1(ip)*dt
   enddo

! Body moved down by vfall*dt. 

   proat = proat-vfall*dt

! define the new position of the body geometry. NOTE: for the spline
! reconstruction, the second derivatives will remain constant

   do i =1,nbody
     zbodyn(i) = zbody(i)+proat
     ybodyn(i) = ybody(i)
   end do

! reconstruct the position of the panel vertices starting from a spline
! interpolation of the panel centroids for the FS. 
! The intersection with the body contour and the far field points are 
! determined as well.

 call splver(ycn,zcn,npc,npsl,proat,estr,tbody,ybodyn,zbodyn,ybodys2,zbodys2,nbody, &
        kffb,iint,ampli,npamx,yn,zn)


! discretization of the wetted body contour. NOTE: in the intermediate
! RK step (first argument is 0) the number of panels is NOT varied 
 
! call ridis(0,proat,escr,npc,npsl,iint,npamx,ybodyn,zbodyn,nbody,ybodys2,zbodys2,& 
!                 yn,zn,ycn,zcn,phin,tbody)


 call ridis6(0,ng,proat,kget,ygb,zgb,&
             escr,npc,npt,yn,zn,ycn,zcn,ampli,&
             ybodyn,zbodyn,ybodys2,zbodys2,tbody,nbody,tgb,iint,ngo1,&
             nsep,ksep,phin,ycb,zcb,tcb,jt,nngo,&
             npsl,di,ang,tc,kord,frint,tn,&
             nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)

! reinitialization of the tangents and normals to the body contour, of
! the panel amplitude and of the boundary conditions on the body
  
 call nortan(vfall,npamx,yn,zn,kget,ng,ygb,zgb, &
             amp,tmy,tmz,rny,rnz,npc,npsl,dphi)


! solution of the boundary value problem based on the new panel
! distribution (yn,zn) and velocity potential (phin)

! Here we need to break into 2 cases kget =0/1


 if (kget == 0)   then

   call dipole(yn,zn,npc,npsl,npf,kffb,kphi,estr,phid,dphid,npamx,amp) 

   call solver( k2dtax,kffb,yn,zn,amp,kphi,phid,dphid,npc,npsl,npf, &
                kget,mget,mgeti,npamx,phi,dphi,coeffdi)
 
   write(*,*) 'coeffdi------=',coeffdi
 
 else  ! (TO BE completed)
   
   call get(ng,npc,npamx,npsl,jt,ycn,zcn,yn,zn, &
                          xigs,zegs,xigb,zegb,xigf,zegf)

   !call solv22(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl, &
   !            ycsl,zcsl,dphi,phisl,dphtb,dphtbsl,phb, &
   !            xigs,zegs,xigb,zegb,xigf,zegf, &
   !            ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl, &
   !            ph,dphn,dphnb,a,b,c,d,e,phi,phib,dphisl,dphibsl,rl, &
   !            mb,mf,mt,m,n,nt,ntt,xi,ze,xis,zes,jt,ksep,kse,kord,kor) 
 
   
   !call calsol  
 endif


write(*,*) 'coeffdi------=',coeffdi 
!stop
! compute tangential velocity components - dpht, normal velocities dphi
! are already returned from the BVP solution


    
 call caveta(ycn,zcn,amp,phin,npc,npsl,npf,kget,mget,mgeti,npamx,dpht)

! Compute the velocity components along the body boundary and free
! surface for the 2ond RK step, v = v_t + v_n

   do i = 1,npc+nng+npsl
     vym2(i) = dpht(i)*tmy(i)+dphi(i)*rny(i)
     vzm2(i) = dpht(i)*tmz(i)+dphi(i)*rnz(i)
   enddo

! DEBUG OUTPUT 

! ----------------- File 1

 call output(yn,zn,ycn,zcn,phin,dphi,vym2,vzm2,npc,npsl,npf,scon,svel, &
              spot,spre,llf,npamx,t)

!   if(jt==1) stop
 
! Displacement of the free surface centroids and potential for the 2ond 
! RK step.
! NOTE: this part will have to be updated to deal with the flow
! separated part, as well as with the modelled part of the jet

! using the old indexation the FS nodes are in the position npc+nng+1

   do ip = npc+nng+1,npc+npsl
     deph2      = vym2(ip)**2 + vzm2(ip)**2
     depn2(ip)  = deph2/2.d0   ! here the gravity term should be included
     ycn(ip)    = yce(ip) + 0.5*(vym1(ip)+vym2(ip))*dt
     zcn(ip)    = zce(ip) + 0.5*(vzm1(ip)+vzm2(ip))*dt
     phin(ip)   = phi(ip) + 0.5*(depn1(ip)+depn2(ip))*dt
   enddo

! after the FS centroids are moved, the panel vertices need to be updated,
! so does the free surface - body intersection, this call returns (yn,zn)
! updated.   

! ----------------- File 2
call output(yn,zn,ycn,zcn,phin,dphi,vym2,vzm2,npc,npsl,npf,scon,svel, &
              spot,spre,llf,npamx,t)

! TO DO: splver will need to be updated to tackle the flow separation
call splver(ycn,zcn,npc,npsl,proat,estr,tbody,ybodyn,zbodyn,ybodys2,zbodys2,nbody, &
        kffb,iint,ampli,npamx,yn,zn)

! DEBUG OUTPUT: NOTE, here the velocities have not been updated yet

!if (jt==jjget+2) stop

! the panel vertices are modified by splver, call nortan to set 
! the BCs to the new panels (yn,zn)  

call nortan(vfall,npamx,yn,zn,kget,ng,ygb,zgb, &                
            amp,tmy,tmz,rny,rnz,npc,npsl,dphi)

!call nortan(vfall,npamx,yn,zn,amp,tmy,tmz,rny,rnz,npc,npsl,dphi)

! identification of the jet region along the body (return ng,ygb,zgb,ygcb,zcgb)
!Warning will need to change this to include separation, it will be a output of splver2
 

 ngo1=nbody
 call shallo(yn,zn,ycn,zcn,kget,jt,npsl,npc,npamx, &
             ng,ampli,rmg,jjget,phin,dphi,epsgg,eskg,gfrac, &
            ybodyn,zbodyn,ybodys2,zbodys2,tbody,tgb,nbody,ngo1,nsep, &
             ycb,zcb,ygb,zgb,tcb,nngo)
 

 nng = ng*kget  ! Update nng in case Jet model is activated very important 
 print*, 'ng================', ng,nng,kget
  
 



! ----------------- File 3
 call output(yn,zn,ycn,zcn,phin,dphi,vym2,vzm2,npc,npsl,npf,scon,svel, &
              spot,spre,llf,npamx,t)

 
! the free surface discretization is updated to guarantee adequate
! accuracy 
 call disuni(jt,yn,zn,ycn,zcn,phin,npsl,npamx,escr,kget,estr, &
                amplim,npc,iiget,ygb,zgb)
   
! DEBUG OUTPUT: NOTE, here the velocities have not been updated yet

! WARNING !!! 
!  if we swicth 0-.1 ridis is cutting all the panels on the body, need 
! to investigate 

 npc_entry =npc

 !call ridis(1,proat,escr,npc,npsl,iint,npamx,ybodyn,zbodyn,nbody,ybodys2,zbodys2,& 
 !             yn,zn,ycn,zcn,phin,tbody)

 call ridis6(0,ng,proat,kget,ygb,zgb,&
             escr,npc,npt,yn,zn,ycn,zcn,ampli,&
             ybodyn,zbodyn,ybodys2,zbodys2,tbody,nbody,tgb,iint,ngo1,&
             nsep,ksep,phin,ycb,zcb,tcb,jt,nngo,&
             npsl,di,ang,tc,kord,frint,tn,&
             nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)

 npc_exit = npc


  do ip = 1, npc + nng
     yv(ip) =  yn(ip)
     zv(ip) =  zn(ip) 
     yce(ip) = ycn(ip)
     zce(ip) = zcn(ip)
     kphi(ip) = 0
  enddo

   !yv(npc+1) = yn(npc+1)
   !zv(npc+1) = zn(npc+1)
 


! ----------------- File 4 : note that here dphi might be wrong if rdis changed 
!the body panels, this is fixed in what follows though
call output(yn,zn,ycn,zcn,phin,dphi,vym2,vzm2,npc,npsl,npf,scon,svel, &
              spot,spre,llf,npamx,t)

! copy the the panels coordinates (yn,zn,ycn,zcn,phin) -> (yv, zv,yce,zce,phi)
! and set Neuman BCs

write(*,*) 'npc, npc_exit,npc_entry, npsl ------', npc, npc_exit,npc_entry,npsl

!if (npc_exit-npc_entry>5) stop 
!if (jt == 17) stop     


! set Dirichlet BCs on the free surface

   do ip=npc+nng+1,npc+nng+npsl
     yv(ip) = yn(ip)
     zv(ip) = zn(ip)
     yce(ip) = ycn(ip)
     zce(ip) = zcn(ip)
     kphi(ip) = 1
     phi(ip) = phin(ip) 
   enddo
   yv(npc+nng+npsl+1) = yn(npc+nng+npsl+1)
   zv(npc+nng+npsl+1) = zn(npc+nng+npsl+1)

! Free surface filtering
   
!   call output(yv,zv,yce,zce,phi,dphi,vym1,vzm1,npc,npsl,scon,svel, &
!              spot,spre,llf,npamx)

   mm = mod(jt,ift) ! filter every 4 time steps

! with separation this call changes as well see main.f lines 635-647
   if(mm.eq.0) then
     call doldfil1(yv,npc+nng+1,npc+nng+npsl+1,npamx,iford,nfil,cdold,ddold)
     call doldfil1(zv,npc+nng+1,npc+nng+npsl+1,npamx,iford,nfil,cdold,ddold)
     call doldfil1(phi,npc+nng+1,npc+nng+npsl,npamx,iford,nfil,cdold,ddold)
   endif

! the panel vertices are modified by doldfill1, call nortan to set 
! the BCs to the new panels (yv,zv)  
!   call nortan(vfall,npamx,yv,zv,amp,tmy,tmz,rny,rnz,npc,npsl,dphi)

   call nortan(vfall,npamx,yv,zv,kget,ng,ygb,zgb, &                
               amp,tmy,tmz,rny,rnz,npc,npsl,dphi)
! solution of the boundary value problem based on the panel
! distribution (yv,zv) and velocity potential (phi)


! ----------------- File 5 
  call output(yv,zv,yce,zce,phi,dphi,vym1,vzm1,npc,npsl,npf,scon,svel, &
             spot,spre,llf,npamx,t)

!  write(*,*) 'llf---------------',llf
!  if (jt == 12)  stop 

! call recalc_dphtsl

  if (kget == 0)   then

     call dipole(yv,zv,npc,npsl,npf,kffb,kphi,estr,phid,dphid,npamx,amp)

     call solver( k2dtax,kffb,yv,zv,amp,kphi,phid,dphid,npc,npsl,npf, &
                kget,mget,mgeti,npamx,phi,dphi,coeffdi)

     write(*,*) 'coeffdi------=',coeffdi

   else  ! (TO BE completed)
   call  get(ng,npc,npamx,npsl,jt,ycn,zcn,yn,zn, &
            xigs,zegs,xigb,zegb,xigf,zegf)

   !call solv22(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl, &
   !            ycsl,zcsl,dphi,phisl,dphtb,dphtbsl,phb, &
   !            xigs,zegs,xigb,zegb,xigf,zegf, &
   !            ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl, &
   !            ph,dphn,dphnb,a,b,c,d,e,phi,phib,dphisl,dphibsl,rl, &
   !            mb,mf,mt,m,n,nt,ntt,xi,ze,xis,zes,jt,ksep,kse,kord,kor) 


   !call calsol  
   endif
  

! check flow conservation
  flux = 0.0    
  do i = 1,npc+nng+npsl 
  flux = flux + amp(i)*dphi(i)
  enddo

  write(*,*) 'jt,flux---------------' , jt, flux 

! compute tangential velocity components - dpht, normal velocities dphi
! are already returned from the BVP solution
    
  call caveta(yce,zce,amp,phi,npc,npsl,npf,kget,mget,mgeti,npamx,dpht)

! compose velocities in the y,z direction respectively for the next cycle

   do i = 1,npc+nng+npsl
     vym1(i) = dpht(i)*tmy(i)+dphi(i)*rny(i)
     vzm1(i) = dpht(i)*tmz(i)+dphi(i)*rnz(i)
   enddo
 
! -------------File 6
  call output(yv,zv,yce,zce,phi,dphi,vym1,vzm1,npc,npsl,npf,scon,svel, &
             spot,spre,llf,npamx,t)



! in the development phase the solution at the intermediate step is
! printed out and the procedure stops 

 enddo  ! end of the time step cycle

stop
end
