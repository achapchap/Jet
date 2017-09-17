
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine output(yv,zv,yc,zc,phi,dphi,vy,vz,npc,npsl,npf,scon,svel, &
           spot,spre,llf,npmax,t)

 implicit none

 integer                  :: npc,npsl,npf,llf,npmax
 real*8                   :: t
 real*8, dimension(npmax) :: yv,zv,yc,zc,phi,dphi,vy,vz
 character(len=2)         :: scon,svel,spot,spre

! locals
 integer :: i
 
! receive 8 vectors and print them on the text file 
! 'scon'XXXXXX contains the vertices
! 'spot'XXXXXX contains the centroids, potential and dphi/dn
! 'svel'XXXXXX contains the centroids and velocity components 

! character(len=8) :: fmt ! format descriptor for saving the results
! fmt = '(I5.5)' ! integer with 5 zeros at the left 

 character(len=5) :: x1

! First Executable statement

 write(x1,'(I5.5)') llf

! output of panel vertices -cd

 open(5,file=scon//trim(x1),status='unknown')
   write(5,*) '# Solution at time ',t
   do i = 1,npc+1
     write(5,*) yv(i),zv(i)
   enddo
   write(5,*)
   write(5,*)
   do i = npc+1,npc+npsl+1
     write(5,*) yv(i),zv(i)
   enddo
   write(5,*)
   write(5,*)
   do i = npc+npsl+1,npc+npsl+npf+1
     write(5,*) yv(i),zv(i)
   enddo
 close(5)

! output of panel centroids and potential - dd

 open(5,file=spot//trim(x1),status='unknown')
   write(5,*) '# Solution at time ',t
   do i = 1,npc
     write(5,4) yc(i),zc(i),phi(i),dphi(i)
   enddo
   write(5,*)
   write(5,*)
   do i = npc+1,npc+npsl
     write(5,4) yc(i),zc(i),phi(i),dphi(i)
   enddo
 close(5)

! output of panel centroids and velocity  - vd

 open(5,file=svel//trim(x1),status='unknown')
   write(5,*) '# Solution at time ',t
   do i = 1,npc
     write(5,4) yc(i),zc(i),vy(i),vz(i)
   enddo
   write(5,*)
   write(5,*)
   do i = npc+1,npc+npsl
     write(5,4) yc(i),zc(i),vy(i),vz(i)
   enddo
 close(5)

4 format(1x,4(e16.8,1x))

 llf = llf+1
 return
endsubroutine

