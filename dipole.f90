subroutine dipole(yv,zv,npc,npsl,npf,kphi,phid,dphid,npamx)
! calculate the far field vertical dipole potential and normal derivative
! and also enforce the boundary condition flag on kphi
  implicit none
  integer :: npamx
  real*8, dimension (npamx)  :: yv,zv,phid,dphid
  integer, dimension (npamx)    :: kphi
  integer :: npc,npsl,npf

! locals 
  integer :: i
  real*8 :: yc,zc,dy,dz,amp,tmy,tmz,rny,rnz,den

 do i = npc+npsl+1,npc+npsl+npf
   yc  = (yv(i+1)+yv(i))/2.d0
   zc  = (zv(i+1)+zv(i))/2.d0
   dy     =  yv(i+1) - yv(i)
   dz     =  zv(i+1) - zv(i)
   amp      = sqrt(dy*dy+dz*dz)
   tmy     = dy/amp
   tmz     = dz/amp
   rny     =  tmz
   rnz     = -tmy
   phid(i)  = yc/(yc**2+zc**2)
   ! dphid = dphid/dx *nx + dphid/dy*ny
   den = ((yc**2+zc**2)**2)
   dphid(i) = ((zc**2-yc**2)/den)*rny+ ((-2*zc*yc)/den)*rnz 
   kphi(i) = 2
 enddo

 return
end
