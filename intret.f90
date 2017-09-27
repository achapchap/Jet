
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine intret(&
!      entry variables
x1,y1,w1x,w1y,x2,y2,w2x,w2y, &
!      variables on exit
xi,yi,s,t)

! finds the intersection point (xi,yi) between two straight lines 
! given in the form (x1,y1) + s(w1x,w1y) and (x2,y2) + t(w2x,w2y) 
! w1 and w2 are not necessarily of unitary modulus

  implicit none

  real*8    :: xi,yi,s,t,x1,y1,w1x,w1y,x2,y2,w2x,w2y

! local variables

  real*8    :: sca,scan,eps,a1,a2

! First executable statement

  eps = 1.d-16

  a1 = w1x**2+w1y**2 
  a2 = w2x**2+w2y**2 
  if (min(a1,a2).lt.eps) then
    write(*,*) 'Error intret: too small segments '
    stop
  endif  

! scala product of one segment by the normal to the other

  sca = w1y*w2x - w1x*w2y
  scan= sca/(sqrt( (w1x**2+w1y**2)*(w2x**2+w2y**2)))
  if (abs(scan).le.eps) then
    write(*,*) 'Error intret: Intersection between parallel lines '
    stop
  endif

! intersection point can be identified in 

  s = ((w2y*x1 - w2x*y1) - (w2y*x2 - w2x*y2))/sca

  t = ((w1y*x1 - w1x*y1) - (w1y*x2 - w1x*y2))/sca

  xi = x1 + s*w1x
  yi = y1 + s*w1y

  return
end
