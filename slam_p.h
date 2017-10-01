!     --  file da includere nel programma slam.f  --

!     - contiene la definizione delle variabili e i parametri per il
!       dimensionamento dei vettori

      implicit real*8 (a-h,o-z)
      character*2 scon,svel,spre,spot
      character*5 sfor

!     parameter (npamx=3000, ntmx = 400000 )
      INTEGER, PARAMETER :: npamx = 3000,  ntmx = 400000 
      common/costanti/pi
