
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine splver( &
!      entry variables
                  xc,zc,npc,npsl,proat,estr,tg,xg,zg,xgs2,zgs2,ngo,& 
                  kffb,iint,ampli,npamx,&
!      variable on exit
                  xp,zp)

! The routine determines the position of the panel vertices of the FS starting
! from the coordinates of the centroids and using a spline interpolation.

 implicit none

 integer      :: iint,ngo,npc,npamx,npsl,kffb
 real*8       :: ampli,proat,estr

 real*8, dimension (npamx)  :: xc,zc,xp,zp,xgs2,zgs2,xg,zg,tg

! Local variables

 integer      :: i,im,ng0
 real*8       :: xi,a,b,c,ax,az,d,dd, xp2,zp2
 real*8       :: delta,den,dis,dt,eps0,epsp,f1,f2,fc,ff,fs1,fs2
 real*8       :: r1,r2,t0,t1,tt1,tt2,tt3,tts,ttt,w1x,w1z,w2x,w2z,x1,x2
 real*8       :: xd1,xd2,xip,xt1,xt2,yp1,ypn,z1,z2,zd1,zd2,zi,zip
 real*8       :: zt1,zt2,coeffdi,c1,c2,cd
 real*8, dimension (npamx)  :: xx,zz,tp,xs2,zs2

 write(*,*) '--> splver2........'

! First executable statement

! The routine computes the vertices of the free surface panels starting
! from a spline interpolation of the centroids

  epsp = 1.d-10

! initialization of the coordinate along the free surface (for spline)
! the abscissa is 0 at the centroid of the first free surface panel

  tp(1) = 0.d0 
  xx(1) = xc(npc+1)
  zz(1) = zc(npc+1)

  do i = 1,npsl-1
    ax      = (xc(npc+i+1)-xc(npc+i))
    az      = (zc(npc+i+1)-zc(npc+i))
    dis     = sqrt(ax*ax+az*az)
    tp(i+1) = tp(i) + dis
    xx(i+1) = xc(npc+i+1)
    zz(i+1) = zc(npc+i+1)
  enddo

! computation of the spline coefficients

  yp1 = 1.d+31
  ypn = 1.d+31
  call spline(xx,zz,xs2,zs2,tp,yp1,ypn,npsl,npamx)

! The vertex at the junction between two successive panel is located 
! at the midpoint of the corresponding centroids. Note: panel j is
! bounded by vertex j and j+1

  do i = 2,npsl
    tts = (tp(i)+tp(i-1))/2.d0
    call splint(tts,xp(npc+i),zp(npc+i),xx,zz,xs2,zs2,tp,npsl,npamx)
  enddo


! Location of the Free-surface/body intersection
! NOTE: this is only necessary BEFORE separation. At the moment
! REATTACHMENT is not accounted for

! The intersection is found by intersecting the line which passes 
! through the first two CENTROIDS of the free surface with the body segments.
! The free surface line is written as (x1,y1) + s(w1x,w1y). The segment
! between the two centroids will have s negative then. The
! intersection should be before (x1,y1), and thus s should be positive
! on exit of intret

  x1  = xc(npc+1)
  z1  = zc(npc+1)
  w1x = xc(npc+1)-xc(npc+2)
  w1z = zc(npc+1)-zc(npc+2)
 
! As the body is represented by a spline, the spline segment where the
! intersection is expected is identified by the line-line intersection.
! Hence, the spline is used to locate the intersection point accurately
! In order to reduce the number of operations, the line-line
! intersection is not done on the entire body segments but only those
! starting from iint-10 to ngo. iint is index of the body segment where the 
! intersection has been found last time

  ng0   = max(iint-10,1)
  iint  =0
  i     = ng0-1

!! NOTE: TO BE MODIFIED in case of separation this 

  do while (iint.eq.0.and.i.le.ngo-1)
    i   = i+1
    x2  = xg(i)
    z2  = zg(i)
    w2x = xg(i+1) - xg(i)
    w2z = zg(i+1) - zg(i)
    call intret(x1,z1,w1x,w1z,x2,z2,w2x,w2z,xi,zi,r1,r2)

! the procedure proceeds until r2 is within (0,1). Hence the
! intersection will be with the segment "iint" given in the global
! body representation (xg,zg,ngo)
     
    if (r2.ge.0.d0.and.r2.le.1.d0)then   
      iint = i

! -- look for intersection between straightline and body spline

      dt   = tg(i+1)-tg(i)
      eps0 = 0.01d0
      xt1  = xg(i)
      zt1  = zg(i)
      xt2  = xg(i+1)
      zt2  = zg(i+1)

      xd1  = xgs2(i)
      zd1  = zgs2(i)
      xd2  = xgs2(i+1)
      zd2  = zgs2(i+1)

      f1   = w1x*zt1 - w1z*xt1
      f2   = w1x*zt2 - w1z*xt2
      fs1  = w1x*zd1 - w1z*xd1
      fs2  = w1x*zd2 - w1z*xd2
      fc   = w1x*z1  - w1z*x1
      den  = fs1 - fs2 
      dd   = abs(1.d0/6.d0*(dt**2)*den) 

! check if the third order coefficient x**3 is < eps

      if(dd.le.epsp)then

        a = (dt**2)*fs2/2.d0
        b = f1 - f2 - (dt**2)*fs1/6.d0 - (dt**2)*fs2/3.d0 
        c = f2 - fc
        if(abs(a).le.epsp)then
          write(*,*) 'aaa! eps piccolo'
          if(abs(b).le.epsp)then
             write(*,*) 'ATTENZIONE!!! x=-c/b, b molto piccolo !!!'
          endif 
          tt1 = -c/b
          tt2 = 1.d31
          tt3 = 1.d31
        else
          delta = b**2 - 4.d0*a*c
          if(delta.lt.0.d0)then
            tt1 = 1.d31
            tt2 = 1.d31
            tt3 = 1.d31
          else
            tt1 = (-b + sqrt(delta))/(2.d0*a)
            tt2 = (-b - sqrt(delta))/(2.d0*a)
            tt3 = 1.d31
          endif
        endif
      else
        a = 3.d0*fs2/den
        b = ( 6.d0/(dt**2)*(f1-f2) - (fs1+2.d0*fs2) )/(fs1-fs2) 
        c = 6.d0/(dt**2)*(f2-fc)/(fs1-fs2)
        call poli3(a,b,c,tt1,tt2,tt3)
      endif

      im = 0
      t0 = -eps0
      t1 = 1.d0+eps0
      if(tt1.ge.t0.and.tt1.le.t1) then
        ttt = tt1
        im  = im + 1
      endif
      if(tt2.ge.t0.and.tt2.le.t1) then
        im  = im + 1
        ttt = tt2
      endif
      if(tt3.ge.t0.and.tt3.le.t1) then
        im  = im + 1
        ttt = tt3
      endif
      if(im.gt.1) write(*,*) 'splver, double intersection '
      if(im.eq.0) write(*,*) 'splver, intersezione non trovata'
        c  = (ttt**3 - ttt)*(dt**2)/6.d0  
        d  = (-(ttt**3) + 3.d0*(ttt**2) - 2.d0*ttt)*(dt**2)/6.d0
        xip = ttt*xt1 + (1.d0-ttt)*xt2 + c*xd1 + d*xd2
        zip = ttt*zt1 + (1.d0-ttt)*zt2 + c*zd1 + d*zd2
    endif
  enddo
  
  if (iint.eq.0.and.i.eq.ngo) then
    write(*,*) ' Intersection not found in Splver, please check '
    stop
  end if

! assign the coordinates to the body-free surface intersection

  xi = xip
  zi = zip
  xp(npc+1) = xi
  zp(npc+1) = zi

! locate the last free surface vertex 
! use the asymptotic behaviour (C/x^2) to assign it. The constant C is
! taken as the average of that computed at the last 3 centroids

  c1 = zc(npc+npsl-10)*(xc(npc+npsl-10)**2)
  c2 = zc(npc+npsl-5)*(xc(npc+npsl-5)**2)
  cd = (c1+c2)/2.d0

  zc(npc+npsl) = cd/(xc(npc+npsl)**2)
  xp(npc+npsl+1) = estr
  zp(npc+npsl+1) = cd/(xp(npc+npsl+1)**2)

end subroutine
