
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
subroutine poli3( &
! variable on entry
a,b,c,x1,x2,x3)

! ---------------------------------------------------------------     
!
! -- compute the zeros (real ones) for a third order equation 
!
!     x**3 + a*(x**2) + b*x + c = 0
!     
!     y = x + a/3 
!
!     Reference : G.Gorni, L'equazione di terzo grado,
!                 notes fo the course on Mathematics by G. Gorni, Univ.
!                 Udine available  online 
!-----------------------------------------------------------------

  use variables
  implicit none   

  real*8         :: a, b, c, x1, x2, x3, delta, p, q , theta, &
                    u1, v1
! --- solution
      p = -(a**2)/3.d0 + b
      q =  2.d0*(a**3)/27.d0 - a*b/3.d0 + c
      write(*,*) ' p q ',p,q

      delta = (q**2)/4.d0 + (p**3)/27.d0
      
      if(delta.ge.0.d0)then
        u1  = (-q/2.d0 + sqrt(delta))
        if(u1.lt.0d0)then
          u1 = -(-u1)**(1.d0/3.d0)
        else
          u1 = u1**(1.d0/3.d0)
        endif
        v1  = (-q/2.d0 - sqrt(delta)) 
        if(v1.lt.0d0)then
          v1 = -(-v1)**(1.d0/3.d0)
        else
          v1 = v1**(1.d0/3.d0)
        endif
        x1  = -a/3.d0 + u1 + v1
        x2  = 1.d31
        x3  = 1.d31
      else
        if(q.le.0.d0)then
          theta = atan(-2.d0*sqrt(-delta)/q)
        else
          theta = pi + atan(-2.d0*sqrt(-delta)/q)
        endif
        x1 = -a/3.d0 + 2.d0*sqrt(-p/3.d0)*cos(theta/3.d0)       
        x2 = -a/3.d0 + 2.d0*sqrt(-p/3.d0)*cos((theta+2.d0*pi)/3.d0)
        x3 = -a/3.d0 + 2.d0*sqrt(-p/3.d0)*cos((theta+4.d0*pi)/3.d0)
      endif

return
end


      
