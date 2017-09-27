
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine input (krest,k2dtax,vfall0,kvfall,pro0,ampp,pfraz,ancut, &
 escr,estr,kffb,tend,frdt,ksta,scon,svel,spot,spre,ift,iford, &
 npamx,nbody,ybody,zbody,jjget,gfrac,rmg,epsgg,eskg,frtend, eskkk,frint,&
 ramii, ramiii)

 implicit none

! npamx = max number of panels on entry

 integer         :: npamx

! ybody,zbody; coordinate of the body geometry (2D contour in the
!           2D+t case) at the original and actual positions

 real*8, dimension(npamx)   :: ybody,zbody

! vfall0: vertical velocity (it can be either constant of variable
!         depending on kvfall)
! pro0: initial submergence of the apex of the body
! ampp: initial FS panel size near the body contour
! pfraz: ratio ampp/smallest panel size
! ancut: angle of cut of the jet
! escr: growth factor of panel size (bulk region)
! estr: domain size (maximum x-coordinate on the right)
! kffb: asymptotic behaviour at far field(0=no, 1=yes)
! tend: simulation time
! frdt: CFL limit (generally 0.25)

 real*8          :: vfall0,pro0,ampp,pfraz,ancut,escr,estr,tend,frdt,frtend, & 
                    eskkk, frint, ramii, ramiii

! nbody: number of points of the body geometry (2D contour in case of
!      2D+t)
! krest: restart option (0: start from undisturbed solution, 1:
!      start from a previous solution given as a "restart" file)
! k2dtax: 2D+t/ax option (0:2D water entry, 1: AX water entry, 2: 2D+t)
! kvfall: type of entry velocity (0: constant entry velocity, 1:
!      variable entry velocity). In the 2D+t case, 0 can be used for
!      steady planing and constant inclination of the keel 
!      (e.g. wedge shaped), otherwise 1 is needed
! ksta: number of time steps between two successive outputs
! ift: time steps for the action of the Dold filter (typically 4)
! iford: order of the Dold filter (typically 3)
! jjget: minimum number of iterations needed before start the jet modelling
! gfrac : fraction of the jet kept when jet_angle<cut_angle (gfrac <= 1)
! rmg : number of panels in the jet base ????, given as input
! espgg :  inferior limit of free surface panels close to the body (amplitude fraction)
! eskg :  stretching factor of the jet panels

 integer  :: nbody,jjget
 integer  :: krest,k2dtax,kvfall,ksta,ift,iford,kffb
 real*8 :: gfrac, rmg, epsgg,eskg

! scon, sve, spre, spot: first two letters of the output files for
! configurations, velocity, pressure and potential

 character(len=2)  :: scon,svel,spre,spot

! local variables

 integer :: ig

! First executive statement

 write(*,*) '--> read input data for discretization and time &
 integration '

 open(8,file='slam.inp',status='UNKNOWN')
   read(8,*) krest
   read(8,*) k2dtax
   read(8,*) vfall0
   read(8,*) kvfall
   read(8,*) pro0 
   read(8,*) ampp 
   read(8,*) pfraz
   read(8,*) ancut
   read(8,*) escr
   read(8,*) estr
   read(8,*) kffb
   read(8,*) tend 
   read(8,*) frdt 
   read(8,*) ksta
   read(8,'(a)') scon
   read(8,'(a)') svel
   read(8,'(a)') spot
   read(8,'(a)') spre
   read(8,*) ift 
   read(8,*) iford
   read(8,*) jjget
   read(8,*) gfrac
   read(8,*) rmg
   read(8,*) epsgg
   read(8,*) eskg
   read(8,*) frtend
!   read(8,*) eskkk
!   read(8,*) frint
     
 close(8)

 write(*,*) '--> read body geometry '

 open(unit=7,file='geo.in',status='unknown')
   read(7,*) nbody
   do ig=1,nbody
     read(7,*) ybody(ig), zbody(ig)
   enddo
 close(7)

! To be fixed and read from file  
 eskkk = 1.1d0  
 frint = 0.25d0 
 ramii = 1.0d0
 ramii = 1.0d0
return
end
