
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

 subroutine check(   &
!     entry variables
 yv,zv,yc,zc,amp,vy,vz,frdt,frtend,tend,t,jt,npc,npsl,npf,npamx,kget,  &
 vfall,zbodyn,nbody,    &
!     variables on entry and exit      (to be completed when the jet
!                                model will be included)
 dt)

 implicit none

 integer                 :: npc,npsl,npf,npamx,nbody,kget,jt,i

 real*8                  :: frdt,t,vfall,dt, frtend, tend
 real*8                  :: cm1,cq1,vvm,velco,frdtt,dtsu,dtk,dta,dista,dtc

 real*8, dimension (npamx)  :: yv,zv,yc,zc,vy,vz,amp,zbodyn

! Then, the time step will be adjusted so that it will
! within (0.90:1.10) times the time step at the previous iteration

! First executable statement
 write(*,*) 'Check-----'
 dtk = dt
 dta = 1.d+10
 frdtt = frdt

 write(*,*) 'JT -----', jt
 if(jt.le.10) frdtt=frdt*jt/10  ! The reduced frdt is used during the
!                                 first 10 iterations
 write(*,*) 'FRDTT ', frdtt

! The routine computes the time step to be used at the next iteration.
! Basically, it starts from the CFL limit, so that the maximum
! displacement of each panel in less than "frdt" times the panel length.

 do i = npc+1,npc+npsl
   vvm = sqrt(vy(i)*vy(i) + vz(i)*vz(i))
   dtc = frdtt*amp(i)/vvm
   dta = min(dtc,dta)

 enddo
 write(*,*) 'amp, vy,vz-----',amp(npc+1), vy(npc+1), vz(npc+1), frdtt 
dtc = frtend*tend
write(*,*) 'dtc dta -------------',dtc,dta
dta = min(dta,dtc)

! Hence, additional constraints are introduced. The first one concerns
! a limit to the time step to avoid that a panel centroid penetrates the
! body contour. 

 cm1 = (zv(npc+1)-zv(npc))/(yv(npc+1)-yv(npc))
 cq1 = zv(npc)-cm1*yv(npc) 

! NOTE: maybe a further check of the following condition is needed when
! the jet modelling is introduced

 dista = cq1+cm1*yc(npc+1)-zc(npc+1)    
 velco = vz(npc+1)-cm1*vy(npc+1)

! the maximum time step (dtsu) is set to 20% of the ratio between
! distance and normal velocity component

 if (velco.gt.0.d0.and.kget.eq.0) then   ! the control operated before
                                         ! the jet modelling is activated 
   dtsu = 0.2*dista/velco
   if (dtsu.lt.dta) then    

! the time step to be used is 80% of dtsu (a further safety margin)

     dta = 0.8*dtsu
   end if
 end if

 dt = dta

! in the old version of the code, there was the additional constraint
! when the jet separates. At present, that constraint is not included.

! the code stops when the time step is too small
 write(*,*) 'vfall----' , vfall 
 if (dt.lt.1.d-15) stop 'Code stops: too small time step '
 
! before closing the routine the vertical position of the body geometry is
! updated by using the actual velocity vfall

 do i=1,nbody
   zbodyn(i)=zbodyn(i)-vfall*dt
 enddo

 return
 end
