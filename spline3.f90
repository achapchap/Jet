
c     ======================================================

subroutine spline3(xs,ys,zs,xs2,ys2,zs2,tp,yp1,ypn,n,npmax)

 implicit none
 real *8, dimension(npmax):: xs,ys,zs
 real *8, dimension(npmax):: xs2,ys2,zs2,tp
 real *8, dimension(npmax+1):: x,y,y2,u
      


!     - interpolazione componente xs

      do iv = 1,n
        x(iv) = tp(iv)
        y(iv) = xs(iv)
      enddo
      if (yp1.gt.1.d+30) then
        y2(1) = 0.d0
        u(1)  = 0.d0
      else
        y2(1) = -0.5d0
        u(1)  = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1)) - yp1)
      end if
      do 11 i = 2,n-1
        sig     = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p       = sig*y2(i-1)+2.d0
        y2(i)   = (sig-1.d0)/p
        u(i)    = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &            /(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig*u(i-1))/p
   11 continue
      if (ypn.gt.1.d+30) then
        qn = 0.d0
        un = 0.d0
      else
        qn = 0.5d0
        un = (3.d0/(x(n) - x(n-1)))*(ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
      end if
      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k = n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
   12 continue
      do iv = 1,n
        xs2(iv) = y2(iv)
      enddo

!     - interpolazione componente ys

      do iv = 1,n
        x(iv) = tp(iv)
        y(iv) = ys(iv)
      enddo
      if (yp1.gt.1.d+30) then
        y2(1) = 0.d0
        u(1)  = 0.d0
      else
        y2(1) = -0.5d0
        u(1)  = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1)) - yp1)
      end if
      do 21 i = 2,n-1
        sig     = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p       = sig*y2(i-1)+2.d0
        y2(i)   = (sig-1.d0)/p
        u(i)    = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &            /(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig*u(i-1))/p
   21 continue
      if (ypn.gt.1.d+30) then
        qn = 0.d0
        un = 0.d0
      else
        qn = 0.5d0
        un = (3.d0/(x(n) - x(n-1)))*(ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
      end if
      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 22 k = n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
   22 continue
      do iv = 1,n
        ys2(iv) = y2(iv)
      enddo

!     - interpolazione componente zs - potential

      do iv = 1,n
        x(iv) = tp(iv)
        y(iv) = zs(iv)
      enddo
      if (yp1.gt.1.d+30) then
        y2(1) = 0.d0
        u(1)  = 0.d0
      else
        y2(1) = -0.5d0
        u(1)  = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1)) - yp1)
      end if
      do 31 i = 2,n-1
        sig     = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p       = sig*y2(i-1)+2.d0
        y2(i)   = (sig-1.d0)/p
        u(i)    = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &            /(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig*u(i-1))/p
   31 continue
      if (ypn.gt.1.d+30) then
        qn = 0.d0
        un = 0.d0
      else
        qn = 0.5d0
        un = (3.d0/(x(n) - x(n-1)))*(ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
      end if
      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 32 k = n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
   32 continue
      do iv = 1,n
        zs2(iv) = y2(iv)
      enddo

      return
      end
