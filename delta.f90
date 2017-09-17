
 function delta(f,c,d,iford)

 implicit none
 integer*4 :: i, iford
 real*8 :: f(-iford:iford)
 real*8 :: delta
 real*8 :: c(1:7,0:7) 
 real*8 ::  d(1:7)

 delta = f(0)*c(iford,0)
 do i=1,iford
   delta = delta + (f(i)+f(-i))*c(iford,i)
 enddo
 delta = delta/d(iford)

return
end
