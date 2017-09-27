
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine nortan( &
! on entry variables                      
                        vfall,npamx,yn,zn,kget,ng,ygb,zgb, &
! on entry and exit variables, these get overwritten                 
                         amp,tmy,tmz,rny,rnz,npc,npsl,dphi)
 
 implicit none

 integer :: npc, npsl, npamx, kget,ng
 real*8, dimension(npamx) :: yn, zn, amp, tmy, tmz, rny, rnz, dphi,ygb,zgb
 real*8 :: vfall
     
! local variables: 
 real *8 :: amy, amz, am
 integer :: ip,nng

 nng = kget*ng

 write(*,*) '--> nortan......'

! Body nodes update normals and body  boundary condition

 do ip      = 1,npc
   amy      = yn(ip+1) - yn(ip)
   amz      = zn(ip+1) - zn(ip)
   am       = sqrt(amy*amy + amz*amz)
   amp(ip)  =  am
   tmy(ip)  =  amy/am
   tmz(ip)  =  amz/am
   rny(ip)  =  tmz(ip)
   rnz(ip)  = -tmy(ip)
   dphi(ip) = -vfall*rnz(ip)
 enddo

! Jet nodes update only the amplitudes on body side 
! and assume no point is separated 
  do ip = 1,nng
    amy = ygb(ip+1)-ygb(ip)
    amz = zgb(ip+1)-zgb(ip)
    am  = sqrt(amy*amy + amz*amz)
    amp(npc+ip) = am
    tmy(npc+ip) =  amy/am
    tmz(npc+ip) =  amz/am
    rny(npc+ip) =  tmz(npc+ip)
    rnz(npc+ip) = -tmy(npc+ip)
    dphi(npc+ip)=  -vfall*rnz(npc+ip)
  enddo


! Same for free surface nodes, i.e update normals, FS BC is integrated in time

 do ip     = npc+1,npc+npsl
   amy     = yn(ip+1) - yn(ip)
   amz     = zn(ip+1) - zn(ip)
   am      = sqrt(amy*amy + amz*amz)
   amp(ip) = am
   tmy(ip) =  amy/am
   tmz(ip) =  amz/am
   rny(ip) =  tmz(ip)
   rnz(ip) = -tmy(ip)
 enddo



! Below the old code  updated the  vectors for the second RK interation 
! For the moment we are doing this in main.f90  
!        if(k2.eq.1)then
!          do ip=1,npc+nng
!            yv(ip) = yn(ip)
!            zv(ip) = zn(ip)
!            yce(ip) = ycn(ip)
!            zce(ip) = zcn(ip)
!            phi(ip) = phin(ip)
!          enddo
!          yv(npc+nng+1) = yn(npc+nng+1)
!          zv(npc+nng+1) = zn(npc+nng+1)
!          yv(npc+1) = yn(npc+1)
!          zv(npc+1) = zn(npc+1)
!          do ip=1,npsl
!            ysl(ip) = ynsl(ip)
!            zsl(ip) = znsl(ip)
!            ycsl(ip) = ycnsl(ip)
!            zcsl(ip) = zcnsl(ip)
!            phisl(ip) = phinsl(ip)
!          enddo
!          ysl(npsl+1) = ynsl(npsl+1)
!          zsl(npsl+1) = znsl(npsl+1)
!        endif

 return        
 end
