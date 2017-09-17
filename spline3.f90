

subroutine spline3(xs,ys,zs,xs2,ys2,zs2,tp,yp1,ypn,n,npmax)

 implicit none
 
 integer :: n,npmax
 real *8, dimension(npmax):: xs,ys,zs
 real *8, dimension(npmax):: xs2,ys2,zs2,tp
 real *8 :: yp1,ypn 
 
! interpolate xs, ys free surface points
 call spline(xs,ys,xs2,zs2,tp,yp1,ypn,n,npmax)

! interpolate the potential only
 call spline1(zs,zs2,tp,yp1,ypn,n,npmax)

return
end
