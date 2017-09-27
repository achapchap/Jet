
 subroutine doldfil1 (v,ni,nf,npamx,iford,nfil,c,d)

!filter the values of v

 implicit none 

 integer :: ni,nf,ntt,npamx,iford, nfil
 real*8 v(npamx+1),vaux(npamx+1)
 real*8    delta
 real*8 :: c(1:nfil,0:nfil)
 real*8 :: d(1:nfil) 

 integer*4 i,j

! First executable statement

!Make a copy vaux of v
 vaux = v

 do i=ni+iford,nf-iford
   v(i) = v(i) - delta (vaux(i-iford),c,d,iford)
 enddo


return
end
