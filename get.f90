!
subroutine get(nget,npc,npamx,npsl,jt,ycn,zcn,yn,zn, &
                          xigs,zegs,xigb,zegb,xigf,zegf)

  implicit none
  integer :: jt,nget,npc,npamx,npsl,nng
  real*8, dimension(npamx) :: ycn,zcn,yn,zn
  real*8, dimension(-npamx:npamx) :: xigs,zegs,xigb,zegb,xigf,zegf

! local stuff
  integer :: i, nb  
  real*8, dimension(npamx) :: ycnsl,zcnsl

  write(*,*) '--> running get..........'

! Parser from the old to the new code , will get rid of it once 
! functionalities are restored 

 nng  = nget
! just parse,  but now the FS vectors  are defined only locally
  do i=1,npsl
    ycnsl(i) = ycn(npc+nng+i)
    zcnsl(i) = zcn(npc+nng+i)
  enddo
!hack ends



  nb = 2
! - calcolo xig  xigs  zegs zeb zef (per solv2)
  do i= -nb,nget
    xigb(i) = ycn(npc-nget+i)
    zegb(i) = zcn(npc-nget+i)

! New FS code
    xigf(i) = ycnsl(nget-i+1)
    zegf(i) = zcnsl(nget-i+1)
    !xigf(i) = ycn(npc+1+i) 
    !zegf(i) = zcn(npc+1+i)

  enddo
  write(95,*) '# ', jt , nget
  do i = -nb+1,nget
    xigs(i) = 0.25d0*(xigb(i)+xigb(i-1)+xigf(i)+xigf(i-1))
    zegs(i) = 0.25d0*(zegb(i)+zegb(i-1)+zegf(i)+zegf(i-1)) 
    write(95,'(i6,10d15.6)')  &
            i,xigs(i),zegs(i), &
            xigb(i),zegb(i),xigf(i),zegf(i), &
            xigb(i-1),zegb(i-1),xigf(i-1),zegf(i-1)
  enddo  
  write(95,*)
  write(95,*)
  return
  end
