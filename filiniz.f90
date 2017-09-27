
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine filiniz (nfil,c,d)

! Initialization of the coefficients for the Dold filter

 implicit none

 integer         :: nfil
 real*8, dimension(1:nfil,0:nfil)      :: c
 real*8, dimension(nfil)               :: d

 real*8 a(1:7,0:7)
 data a(1,0),a(1,1),a(1,2),a(1,3),a(1,4),a(1,5),a(1,6),a(1,7) &
    /2.d0  ,-1.d0 ,0.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0/
 data a(2,0),a(2,1),a(2,2),a(2,3),a(2,4),a(2,5),a(2,6),a(2,7) &
    /6.d0  ,-4.d0 ,1.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0  ,0.d0/
 data a(3,0),a(3,1),a(3,2),a(3,3),a(3,4),a(3,5),a(3,6),a(3,7) &
    /20.d0 ,-15.d0,6.d0  ,-1.d0 ,0.d0  ,0.d0  ,0.d0  ,0.d0/
 data a(4,0),a(4,1),a(4,2),a(4,3),a(4,4),a(4,5),a(4,6),a(4,7) &
    /70.d0 ,-56.d0,28.d0 ,-8.d0 ,1.d0  ,0.d0  ,0.d0  ,0.d0/
 data a(5,0),a(5,1) ,a(5,2),a(5,3),a(5,4),a(5,5),a(5,6),a(5,7) &
    /252.d0,-210.d0,120.d0,-45.d0,10.d0 ,-1.d0 ,0.d0  ,0.d0/
 data a(6,0),a(6,1) ,a(6,2),a(6,3) ,a(6,4),a(6,5),a(6,6),a(6,7) &
    /924.d0,-792.d0,495.d0,-220.d0,66.d0 ,-12.d0,1.d0  ,0.d0/
 data a(7,0) ,a(7,1)  ,a(7,2) ,a(7,3)  ,a(7,4),a(7,5),a(7,6),a(7,7) &
    /3432.d0,-3003.d0,2002.d0,-1001.d0,364.d0,-91.d0,14.d0 ,-1.d0/

 integer*4, dimension(nfil)         :: e

 integer*4            :: i,j

 if (nfil.gt.7) stop 'dimensions exceeded in filiniz '
 
 e = (/ (2*i, i=1,nfil) /)

 do i=1,nfil
   d(i) = 2.d0**e(i)
   do j=0,nfil
     c(i,j) = a(i,j)
   enddo
 enddo

return
end
