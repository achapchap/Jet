
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine ridis(&
! entry variables
                 kk,proat,escr,npc,npsl,iint,npamx,ybodyn,zbodyn,nbody,ybodys2,zbodys2, &
! entry/exit these get overwritten     
                 yn,zn,ycn,zcn,phin,tbody)

! Description of the routine:      
! This routine discretizes the body mesh based on the amplitude of 
! the first free surface panel. 

implicit none

integer  :: k,npc,npsl,npamx,iint,nbody
real*8   :: proat,escr
real*8, dimension(npamx) :: yn,zn,ycn,zcn,phin,tbody,tn,ybodyn,zbodyn,ybodys2,zbodys2
   
! local variables:

real*8  :: yy,zz,amii,ay,az,aa,amdi,tt,r,dr,yf,zfi,eskk,zf,y,z
integer :: nn,ninc,i,kk,ip,no

! First executable statement

 eskk = escr
 ninc = 0

 yf = 0.d0
 zf = proat

!!!!!!     if (ng.eq.0) then
   yy   = yn(npc+1)
   zz   = zn(npc+1) 
   amii = sqrt((yn(npc+2)-yn(npc+1))**2+ (zn(npc+2)-zn(npc+1))**2 )
   ay   = yy - ybodyn(iint)  ! this is submerged length in y direction
   az   = zz - zbodyn(iint)  ! same for vertical diretcion
   aa   = sqrt(ay**2+az**2)
   amdi = tbody(iint) + aa
   tt   = amdi 
!!     else
!!   write(*,*) ' Not implemented yet '
!!   stop
!! end if
   
! compute the number of panels to be used along the wetted body contour 
! by assuming the last body panel (i.e. the one arriving at the 
! intersection with the free surface) to be of the same length of the first 
! free surface panel (i.e. amii) 
! nn is the number of panels that are needed to discretize the wetted
! part of the body contour using amii as initial panel amplitude and eskk
! as growing factor

  nn = int( log( 1.d0+(eskk-1.d0)*amdi/amii)/log(eskk) )

!!!!!! no = npc-ng*(1-kget)

  no   = npc 
  ninc = nn-no

! WARNING: TO BE modified when the jet will be included

  if (ninc.lt.0) then
    ninc = 0 
    nn   = no
  end if

  if(kk.eq.0) then
    ninc = 0
    nn   = npc
  else
 ! displace the Free surface panels along the arrays in order to write the 
 ! new coordinates of the body panels we are removing/adding   
    if (ninc.lt.0) then
      write(*,*) 'Cut ',ninc,' panels from the body '

      do i = npc+1,npc+npsl
        yn(i+ninc) = yn(i)
        zn(i+ninc) = zn(i)
        ycn(i+ninc) = ycn(i)
        zcn(i+ninc) = zcn(i)
        phin(i+ninc) = phin(i)
      enddo
      yn(npc+npsl+ninc+1) = yn(npc+npsl+1)
      zn(npc+npsl+ninc+1) = zn(npc+npsl+1)

    else if (ninc.gt.0)then
      write(*,*) 'Add ',ninc,' panels to the body ' 

      yn(npc+npsl+ninc+1) = yn(npc+npsl+1)
      zn(npc+npsl+ninc+1) = zn(npc+npsl+1)
      do i = npc+npsl,npc+1,-1
        yn(i+ninc) = yn(i)
        zn(i+ninc) = zn(i)
        ycn(i+ninc) = ycn(i)
        zcn(i+ninc) = zcn(i)
        phin(i+ninc) = phin(i)
      enddo

    endif

  end if

! compute the exact value of the amplitude of the first panel
  
  amii = amdi*(1.d0-eskk)/(1.d0-eskk**nn)
!  write(*,*) 'amii, amdi, eskk,nn', amii, amdi, eskk,nn   
!  write(*,*) '              ip,   r, tbody  ,       yn(ip),      zn(ip),     y ,       z'
! discretization of the body
  y = 0.
  z = 0.
  r   = amdi    ! curvilinear abscissa at the FS-body intersection 
  dr  = amii
! write new vertices on the body 
  do ip = nn,1,-1
    r = r-dr  
    if(ip.eq.1) r = 0.d0
    call splint(r,y,z,ybodyn,zbodyn,ybodys2,zbodys2,tbody,nbody,npamx)
 !   write(*,*)   ip,r,tbody(ip),yn(ip),zn(ip),y,z


    tn(ip)   = r
    yn(ip)   = y 
    zn(ip)   = z
    dr       = dr*eskk
  enddo 
  tn(nn+1) = amdi

  npc = nn 

! write new centroids on the body

  tt = 0.d0
  do i = 1,npc
    tt = 0.5d0*(tn(i)+tn(i+1)) 
    call splint(tt,y,z,ybodyn,zbodyn,ybodys2,zbodys2,tbody,nbody,npamx)
    ycn(i)  = y 
    zcn(i)  = z
  enddo

return
end
