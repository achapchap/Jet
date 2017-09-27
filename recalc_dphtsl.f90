subroutine  recalc_dphtsl(ycn,zcn,phi,npamx,npc,npsl,nng,&
             dphtsl)


  integer :: npc,npamx,npsl,nng
  real*8, dimension(npamx) :: ycn,zcn,phi,dphtsl

! locals 
  real*8, dimension(npamx) :: ycsl,zcsl
  real*8 :: dff,dfb,tb,tf    


! small hack 

! just parse but now the FS vectors  are defined only locally
  do i=1,npsl
    ycsl(i) = ycn(npc+nng+i)
    zcsl(i) = zcn(npc+nng+i)
  enddo
!hack ends

  do i=1,npsl
    if(i.eq.1)then
      tf = sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
      dff= (phisl(i+1)-phisl(i))/tf
      dphtsl(i) = dff
    elseif(i.eq.npsl)then
      tb = sqrt( (ycsl(i-1)-ycsl(i))**2 + (zcsl(i-1)-zcsl(i))**2 )
      dfb= (phisl(i)-phisl(i-1))/tb
      dphtsl(i) = dfb
    else
      tf = sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
      tb = sqrt( (ycsl(i-1)-ycsl(i))**2 + (zcsl(i-1)-zcsl(i))**2 )
      dff= (phisl(i+1)-phisl(i))/tf
      dfb= (phisl(i)-phisl(i-1))/tb
      dphtsl(i) = 0.5d0*(dff+dfb)
    endif
  enddo

return
end


