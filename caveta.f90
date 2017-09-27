
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine caveta( &
!      entry variables
 yc,zc,amp,phi,npc,npsl,npf,kget,mget,mgeti,npamx,  &

!      variables on exit

 dpht)

 implicit none

 integer               :: npc,npsl,npf,nps,kget,&
                          mget,mgeti,npamx

 real*8, dimension (npamx)  :: yc,zc,amp,phi,dpht

! Local variables

 integer               :: i
 real*8                :: dff, dfb

! First executable statement

! The routine computes the tangential velocity based on a second order 
! finite difference scheme. Forward and backward schemes are used at the
! extremes

! Body side

 dpht(1) = 2.d0*(phi(2)-phi(1))/(amp(2)+amp(1))
 
 do i = 2,npc-1
   dff = 2.d0*(phi(i+1)-phi(i))/(amp(i+1)+amp(i))
   dfb = 2.d0*(phi(i)-phi(i-1))/(amp(i-1)+amp(i))
   dpht(i) = (dff+dfb)/2.d0
 enddo

 dpht(npc) = 2.d0*(phi(npc)-phi(npc-1))/(amp(npc)+amp(npc-1))

! free surface side 
! NOTE: to be modify to treat the modelled part of the jet 

 dpht(npc+1) = 2.d0*(phi(npc+2)-phi(npc+1))/ &
                         (amp(npc+2)+amp(npc+1))
 
 do i = npc+2,npc+npsl-1
   dff = 2.d0*(phi(i+1)-phi(i))/(amp(i+1)+amp(i))
   dfb = 2.d0*(phi(i)-phi(i-1))/(amp(i-1)+amp(i))
   dpht(i) = (dff+dfb)/2.d0
 enddo

 dpht(npc+npsl) = 2.d0*(phi(npc+npsl)-phi(npc+npsl-1))/ & 
                            (amp(npc+npsl)+amp(npc+npsl-1))

 return
 end
