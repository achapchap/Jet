
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine solver( &
!      entry variables
 k2dtax,kffb,yv,zv,amp,kphi,phid,dphid,npc,npsl,npf, &
 kget,mget,mgeti,npamx,  &
!      variables on entry and exit  (to be completed when the jet model
!                                    will be included)
 phi,dphi,coeffdi)

 use variables 
 implicit none

 integer               :: k2dtax,kffb,npc,npsl,npf,nps,npt,kget,&
                          mget,mgeti,npamx

 integer, dimension(npamx)  :: kphi

 real*8, dimension (npamx)  :: yv,zv,amp,phi,dphi, phid, dphid

 real*8  :: coeffdi
! Local variables

 integer         :: i, j, nunk, npe
 real*8          :: fint1, fins1,fint2, fins2, yy, zz, ys, yp, zs, zp, &
                    yss, yps, zss, zps, idum, a_diag

 real*8, dimension(:,:), allocatable :: aa
 real*8, dimension(:), allocatable  :: bb
 integer, dimension(:), allocatable  :: indx
 character(len=1)      :: trans

! First executable statement

 write(*,*) ' Far field asymptotics not implemented yet!! '
 write(*,*) 'npf----------',npf
! count of the number of unknowns (+1*kffb accounts for the FF boundary
! condition, if used, +5*mget is the number of unknowns of the modelled
! part of the jet, +10 is for some additional unknowns)


! TO DO : check the size of nunk and see if it makes sense for npf =  0  and npf !=0
 nunk = npc+npsl+npf+1*kffb+5*mget+1

! allocate the arrays 

 allocate(aa(nunk,nunk))
 allocate(bb(nunk))
 allocate(indx(nunk))

! coefficient matrix: boundary integral equation part

  aa = 0.d0
  bb = 0.d0

! unknowns are ordered as follows:
!  - velocity potential phi on the body side, including the modelled
!    part, if any (npc)
!  - normal derivative dphi on the separated portion of the jet,
!    including the modelled part, if any (npse)
!  - normal derivative dphi on the free surface on the right of the jet
!    tip, including the modelled part, if any (npsl)
!  - unknown (either phi or dphi) for panels at the far field
!    boundary (nff)
!  - dipole constant for the ff behaviour (if used)
!  - coefficients of the hybrid jet model

! At present it is assumed that there is no modelled part for the jet

!TO DO : implement a simple subroutine that calculates the dipole potential for the far field
! i.e it will initialize the vectors dphid and phid for the dipole far field approx 

  do i = 1,npc+npsl+npf
    yy = 0.5*(yv(i)+yv(i+1))
    zz = 0.5*(zv(i)+zv(i+1)) 
    do j = 1,npc+npsl+npf
      yp = yv(j)
      ys = yv(j+1)
      zp = zv(j)
      zs = zv(j+1)

! enforce the symmetry with respect to the x=0 axis

      yps= -yv(j+1)
      yss= -yv(j)
      zps= zv(j+1)
      zss= zv(j)

! influence coefficient of panel + image panel 
!      write(*,*) 'calculating influence-----',i,j,yy,zz,yp,zp,ys,zs,amp(j)
      call finte(yy,zz,yp,zp,ys,zs,amp(j),fint1,fint2)
      call finte(yy,zz,yps,zps,yss,zss,amp(j),fins1,fins2)

! built of the matrix and know term

      if (kphi(j).eq.0) then   ! normal velocity is known on j
        aa(i,j) = - (fint2+fins2)
        bb(i)   = bb(i) + (fint1+fins1)*dphi(j)
      else if (kphi(j).eq.1) then ! velocity potential is known on j
        aa(i,j) = - (fint1+fins1)
        bb(i)   = bb(i) + (fint2+fins2)*phi(j)
! FF Dipole term
      else if (kphi(j)==2) then
        aa(i,j) = -(fint1+fins1)
        aa(i,npc+npsl+npf+1) = aa(i,npc+npsl+npf+1)-(fint2+fins2)*phid(j)
      end if
    enddo

! self-induced contribution

    if (kphi(i).eq.0) then
      aa(i,i) = aa(i,i) + pi 
    else if (kphi(i).eq.1) then
      bb(i) = bb(i) - pi*phi(i)
    else if (kphi(i)==2) then
       aa(i,npc+npsl+npf+1) = aa(i,npc+npsl+npf+1)+pi*phid(i)
    end if
  enddo 

! last line of the far field approx
  if (kffb==1) then
    do j =1, npc+npsl+npf+1
      aa(npc+npsl+npf+1,j) = 0.d0
    enddo

    a_diag  = 0.0d0
    do j =npc+npsl+1, npc+npsl+npf
      aa(npc+npsl+npf+1,j) = -amp(j)
      a_diag = a_diag + dphid(j)*amp(j)
    enddo
    aa(npc+npsl+npf+1,npc+npsl+npf+1) = a_diag
    bb(npc+npsl+npf+1) = 0.0d0
    npe = npc+npsl+npf+1
  else
    npe = npc+npsl
  endif
! solution of the Boundary value problem


  call DGETRF(npe,npe,aa,nunk,indx,idum)
  trans='N'
  call DGETRS(trans,npe,1,aa,nunk,indx,bb,nunk,idum)

! copy on the solution on phi and dphi

  do i = 1,npc+npsl+npf
    if (kphi(i).eq.0) then
      phi(i) = bb(i)
    else if (kphi(i).eq.1) then
      dphi(i) = bb(i)
    else if (kphi(i) == 2) then
! f_n on the farfield panels
      dphi(i) = bb(i)
    end if
  enddo

! Dipole coefficient
  if (kffb==1) then 
    coeffdi = bb(npc+npsl+npf+1)
  else
    coeffdi = 0.0 
  endif  
 
 return
 end
