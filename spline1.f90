
!     ======================================================

subroutine spline1(xs,xs2,tp,yp1,ypn,n,nmax)

 implicit none

 integer          :: n,nmax

 real*8, dimension(nmax)    :: xs,xs2,tp

! Local variables

 integer                      :: iv,i,k
 real*8                       :: yp1,ypn,p,sig,qn,un
 real*8, dimension(nmax)      :: u

! First executable statement

! Interpolation of xs

 if (yp1.gt.1.d+30) then
   xs2(1) = 0.d0
   u(1)  = 0.d0
 else
   xs2(1) = -0.5d0
   u(1)  = (3.d0/(tp(2)-tp(1)))*((xs(2)-xs(1))/(tp(2)-tp(1)) - yp1)
 end if
 do i = 2,n-1
   sig     = (tp(i)-tp(i-1))/(tp(i+1)-tp(i-1))
   p       = sig*xs2(i-1)+2.d0
   xs2(i)   = (sig-1.d0)/p
   u(i)    = (6.d0*((xs(i+1)-xs(i))/(tp(i+1)-tp(i))-(xs(i)-xs(i-1)) &
          /(tp(i)-tp(i-1)))/(tp(i+1)-tp(i-1)) - sig*u(i-1))/p
 enddo
 if (ypn.gt.1.d+30) then
   qn = 0.d0
   un = 0.d0
 else
   qn = 0.5d0
   un = (3.d0/(tp(n) - tp(n-1)))*(ypn - (xs(n)-xs(n-1))/(tp(n)-tp(n-1)))
 end if
 xs2(n) = (un - qn*u(n-1))/(qn*xs2(n-1)+1.d0)
 do k = n-1,1,-1
   xs2(k) = xs2(k)*xs2(k+1)+u(k)
 enddo

return 
end


