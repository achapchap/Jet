
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine solv22(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl, &
                 ycsl,zcsl,dphi,phisl,dphtb,dphtbsl,phb, &
                 xigs,zegs,xigb,zegb,xigf,zegf, &
                 ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl, &
                 ph,dphn,dphnb,a,b,c,d,e,phi,phib,dphisl,dphibsl,rl, &    
                 mb,mf,mt,m,n,nt,ntt,xi,ze,xis,zes,jt,ksep,kse,kord,kor)



!      include"slam_p.h"

      implicit real*8 (a-h,o-z)
      parameter (npamx=3000, ntmx = 400000 )
      common/costanti/pi

      dimension dphn(npamx),dphi(npamx),ph(npamx),dphtb(npamx)
      dimension dphtsl(npamx)
      dimension phisl(npamx),phb(npamx),dphnb(npamx),dphtbsl(npamx)  
      dimension dphisl(npamx),dphibsl(npamx),phi(npamx),phib(npamx)
      dimension yc(npamx),zc(npamx),yce(npamx),zce(npamx)
      dimension yb(npamx),zb(npamx),yf(npamx),zf(npamx)
      dimension yv(npamx),zv(npamx),ysl(npamx),zsl(npamx)
      dimension ycsl(npamx),zcsl(npamx),am(npamx)
      dimension xig(-npamx:npamx),xigs(-npamx:npamx)
      dimension zeb(-npamx:npamx),zegs(-npamx:npamx)
      dimension zef(-npamx:npamx)
      dimension xi(npamx),xis(npamx),ze(npamx),zes(npamx)
      dimension aa(npamx,npamx),bb(npamx),appo(npamx)
      dimension ff(npamx),rl(npamx)
      dimension a(npamx),b(npamx),c(npamx),d(npamx),e(npamx)
      dimension indx(npamx)
      dimension xigb(-npamx:npamx),zegb(-npamx:npamx)
      dimension xigf(-npamx:npamx),zegf(-npamx:npamx)
      dimension ty(npamx) ,tz(npamx),ry(npamx),rz(npamx) 
      dimension tmy(npamx),tmz(npamx),rny(npamx),rnz(npamx) 
      dimension tmysl(npamx),tmzsl(npamx),rnysl(npamx),rnzsl(npamx)
      dimension ksep(npamx),kse(npamx),kor(npamx),kord(npamx) 
      parameter (naux=3*npamx)
      character*1 trans

! variables description (only the more obscure for the moment):
! phi,dphi: velocity potential and normal derivative on the body
! panels
! phisl,dphisl: velocity potential and normal derivative on the
! free-surface panels
! phinsl: free-surface velocity potential at the second R-K level
! phin: velocity potential in the body + jet  region at the second R-K level
! dpht,dphtsl: tangential derivative of the velocity potential on
! body and free-surface, respectively

! dphtbsl,dphn,dphnb,phb,ph,phib: these arrays are used in an overlapping
! region between the bulk of the fluid and the modelled part of
! the jet in order to smooth the transition between the two regions
! They have to be further understood but it seems they are
! only used inside solv22 and solv22p (They get modified here)

! hp,xj,ze,xis,zes: these variables seems still related to the jet
! modelling but not all of them are used at present. Better
! keep them in the debugging phase

! dpttsl,dpb,dptb,dpntbsl: these variables are used in solv22p and
! seems to be related to the solution of the Laplace equation
! for the time derivative of the velocity potential but their
! use isn't clear

! vxi,dpntt: vxi is the tangential velocity in the jet and dpntt is
! the second derivative of the velocity potential. They seem
! to be used for the pressure

! kord, kor: are indices used in the modelled part of the jet and/or
! in the transition region, but it is not clear how they work

! rl: matching factor in the matching region
! a1,b1,c1,d1,e1: coefficients of the FEM model in the modelled
! region of the jet

! a2,b2,c2,d2,e2: as before but for the time derivative of the
! velocity potential (pressure solution)

!frint   = matching region (0=no matching, frint<0.9))

! Important vectors on exit , where the solution is saved are"
! ph, phb, dphnn, dphnb and a,b,c,d,e 



  m  = int(frint*ng)
  n  = ng - m    ! this defines the part where the pure FEM model is used

  mb = npc+ng-(m+n)   ! There is an apparent difference in the numbering
            ! of the panels on the body and on the free surface 
            ! side. Indeed, npc only accounts for the panels on the 
            ! body side belonging to the "bulk" of the fluid whereas the 
            ! panels in the modelled portion of the jet are not
            ! included in npc. 
            ! This explains why the equation mb = npc+ng-(m+n) has been
            ! used in spite it always returns mb=npc.

  mf = npsl-(m+n)       
  mt = mb + mf
  nt = mt + 2*m
  ntt= nt + 2*n
  write(*,*) npc,ng,npsl
  write(*,*) mb,m,n
  write(*,*) mf,m,n
  write(*,*) mt,nt,ntt

! - Preamble!!!

  do i=1,mb+m+n
    if(i.le.mb)then     ! this should be on the body surface but it is
                        ! not clear why npc should inclde the jet part 
                        ! whereas npsl should not
      dphn(i) = dphi(i)
      dphnb(i)= dphi(i)
      ph(i)   = phi(i)
      phb(i)  = phi(i)
      kse(i)  = ksep(i)
      kor(i)  = kord(i) ! where does kord comes from is still not understood !
      yb(i)   = yv(i)
      zb(i)   = zv(i)
      yf(i)   = yv(i+1)
      zf(i)   = zv(i+1)
      yc(i)   = 0.5d0*(yf(i)+yb(i))
      zc(i)   = 0.5d0*(zf(i)+zb(i))
      ry(i)   = rny(i)
      rz(i)   = rnz(i)
      ty(i)   = tmy(i)
      tz(i)   = tmz(i)
!      write(71,'(i6,10d15.7)') i,yc(i),zc(i), &
!                yb(i),zb(i),yf(i),zf(i), &
!                ry(i),rz(i),ty(i),tz(i)

    elseif(i.ge.mb+1.and.i.le.mb+m)then    ! considering the above
                                           ! definition of mb, this should 
                                           ! involve panel on the free surface 
                                           ! side

                ! the first time this part is entered, i=mb+1
                ! and thus (mt+i-mb) is equal to (mb+mf+mb+1-mb)=(mb+mf+1)
                ! the last time this part is enterd, i=mb+m and thus
                ! (mt+i-mb) is equal to (mb+mf+mb+m-mb)=(mb+mf+m)

      dphn(mt+i-mb) = dphi(i)
      dphnb(mt+i-mb)= dphi(i)
      ph(mt+i-mb)   = phi(i)
      phb(mt+i-mb)  = phi(i)
      kse(mt+i-mb)  = ksep(i)
      kor(mt+i-mb)  = kord(i) ! same question as above ...
      yb(mt+i-mb)  = yv(i)
      zb(mt+i-mb)  = zv(i)
      yf(mt+i-mb)  = yv(i+1)
      zf(mt+i-mb)  = zv(i+1)
      yc(mt+i-mb)  = 0.5d0*(yf(mt+i-mb)+yb(mt+i-mb))
      zc(mt+i-mb)  = 0.5d0*(zf(mt+i-mb)+zb(mt+i-mb))
      ry(mt+i-mb)= rny(i)
      rz(mt+i-mb)= rnz(i)
      ty(mt+i-mb)= tmy(i)
      tz(mt+i-mb)= tmz(i)
!      write(72,'(i6,10d15.7)') mt+i-mb,yc(mt+i-mb),zc(mt+i-mb), &
!                              yb(mt+i-mb),zb(mt+i-mb),yf(mt+i-mb),zf(mt+i-mb), &
!                              dphn(mt+i-mb),dphnb(mt+i-mb) &
!                              ry(mt+i-mb),rz(mt+i-mb),ty(mt+i-mb),tz(mt+i-mb)
                       
    else
                   ! the fisrt time this part is enterd, i=mb+m+1, and
                   ! then (nt+i-(mb+m)) is equal to
                   ! (mb+mf+2*m+mb+m+1-mb-m) which is (mf+mb+2*m+1)
                   ! the last time this part is enteerd, i=mb+m+n
                   ! then (nt+i-(mb+m)) is equal to (mb+mf+2*m+mb+m+n-mb-m) 
                   ! which is (mf+mb+2*m+n)

      dphn(nt+i-(mb+m)) = dphi(i)
      dphnb(nt+i-(mb+m))= dphi(i)
      ph(nt+i-(mb+m))   = phi(i)
      phb(nt+i-(mb+m))  = phi(i)
      kse(nt+i-(mb+m))  = ksep(i)
      kor(nt+i-(mb+m))  = kord(i)
      yb(nt+i-(mb+m))  = yv(i)
      zb(nt+i-(mb+m))  = zv(i)
      yf(nt+i-(mb+m))  = yv(i+1)
      zf(nt+i-(mb+m))  = zv(i+1)
      yc(nt+i-(mb+m))  = 0.5d0*(yf(nt+i-(mb+m))+yb(nt+i-(mb+m)))
      zc(nt+i-(mb+m))  = 0.5d0*(zf(nt+i-(mb+m))+zb(nt+i-(mb+m)))
      ry(nt+i-(mb+m))  = rny(i)
      rz(nt+i-(mb+m))  = rnz(i)
      ty(nt+i-(mb+m))  = tmy(i)
      tz(nt+i-(mb+m))  = tmz(i)
!      write(73,'(i6,10d15.7)') &
!                nt+i-(mb+m),yc(nt+i-(mb+m)),zc(nt+i-(mb+m)), &
!                yb(mt+i-(mb+m)),zb(mt+i-(mb+m)),yf(mt+i-(mb+m)),zf(mt+i-(mb+m)), &
!                dphn(nt+i-(mb+m)),dphnb(nt+i-(mb+m))  &
!                ry(nt+i-(mb+m)),rz(nt+i-(mb+m)),ty(nt+i-(mb+m)),tz(nt+i-(mb+m))
    endif
  enddo

!      write(71,*)
!      write(71,*)
!      write(72,*)
!      write(72,*)
!      write(73,*)
!      write(73,*)

! TO DO: reindex the loop  (to mb+n+m to mf+m+n+mb+n+m ? )

  do i=1,mf+m+n
    if(i.le.n)then      ! this spans over the FEM part of the jet

                        ! the first time this part is entered, i=1 and
                        ! thus (nt+n+i) is equal to (mb+mf+2*m+n+1) which 
                        ! is just the element next top the last one in the 
                        ! previous loop

      dphtbsl(i)   = dphtsl(i)
      dphtb(nt+n+i)= dphtbsl(i)
      ph(nt+n+i)   = phisl(i)
      phb(nt+n+i)  = phisl(i)
      yb(nt+n+i)  = ysl(i)
      zb(nt+n+i)  = zsl(i)
      yf(nt+n+i)  = ysl(i+1)
      zf(nt+n+i)  = zsl(i+1)
      yc(nt+n+i)  = 0.5d0*(yf(nt+n+i)+yb(nt+n+i))
      zc(nt+n+i)  = 0.5d0*(zf(nt+n+i)+zb(nt+n+i))
      ry(nt+n+i)  = rnysl(i)
      rz(nt+n+i)  = rnzsl(i)
      ty(nt+n+i)  = tmysl(i)
      tz(nt+n+i)  = tmzsl(i)
!      write(87,'(i6,10d15.7)') nt+n+i,yc(nt+n+i),zc(nt+n+i), &
!                          yb(nt+n+i),zb(nt+n+i),yf(nt+n+i),zf(nt+n+i), &
!                          dphtb(nt+n+i),ph(nt+n+i),phb(nt+n+i) &
!                          ry(nt+n+i),rz(nt+n+i),ty(nt+n+i),tz(nt+n+i) &
    elseif(i.ge.n+1.and.i.le.n+m)then

                        ! the first time this part is enterd, i=n+1 and
                        ! thus (mt+m+(i-n))=(mb+mf+m+(n+1-n))=(mb+mf+m+1)
                        ! the last time i=n+m and thus(mt+m+(i-n)) is
                        ! equal to (mb+mf+m+(n+m-n)) = (mb+mf+2*m)

      dphtbsl(i)       = dphtsl(i)
      dphtb(mt+m+(i-n))= dphtbsl(i)
!     dpht(mt+m+(i-n)) = dphtsl(i)
      ph(mt+m+(i-n))   = phisl(i)
      phb(mt+m+(i-n))   = phisl(i)
!     yc(mt+m+(i-n))   = ycsl(i)
!     zc(mt+m+(i-n))   = zcsl(i)
      yb(mt+m+(i-n))   = ysl(i)
      zb(mt+m+(i-n))   = zsl(i)
      yf(mt+m+(i-n))   = ysl(i+1)
      zf(mt+m+(i-n))   = zsl(i+1)
      yc(mt+m+(i-n))   = 0.5d0*(yf(mt+m+(i-n))+yb(mt+m+(i-n)))
      zc(mt+m+(i-n))   = 0.5d0*(zf(mt+m+(i-n))+zb(mt+m+(i-n)))
      ry(mt+m+(i-n))   = rnysl(i)
      rz(mt+m+(i-n))   = rnzsl(i)
      ty(mt+m+(i-n))   = tmysl(i)
      tz(mt+m+(i-n))   = tmzsl(i)
!     write(82,'(i6,10d15.7)') mt+m+i-n,yc(mt+m+i-n),zc(mt+m+i-n), &
!                 yb(mt+m+i-n),zb(mt+m+i-n),yf(mt+m+i-n),zf(mt+m+i-n), &
!                 dphtb(mt+m+i-n), ph(mt+m+i-n), phb(mt+m+i-n), &
!                 ry(mt+m+i-n),rz(mt+m+i-n),ty(mt+m+i-n),tz(mt+m+i-n)
    else

              ! the first time this part is entered, i=n+m+1 and thus
              ! (mb+i-(n+m)) = (mb+n+m+1-n-m) = mb+1

      dphtbsl(i)       = dphtsl(i)
      dphtb(mb+i-(n+m))= dphtbsl(i)
      ph(mb+i-(n+m))   = phisl(i)
      phb(mb+i-(n+m))   = phisl(i)

      yb(mb+i-(n+m))   = ysl(i)
      zb(mb+i-(n+m))   = zsl(i)
      yf(mb+i-(n+m))   = ysl(i+1)
      zf(mb+i-(n+m))   = zsl(i+1)

      yc(mb+i-(n+m))   = 0.5d0*(yf(mb+i-(n+m))+yb(mb+i-(n+m)))
      zc(mb+i-(n+m))   = 0.5d0*(zf(mb+i-(n+m))+zb(mb+i-(n+m)))
      ry(mb+i-(n+m))   = rnysl(i)
      rz(mb+i-(n+m))   = rnzsl(i)
      ty(mb+i-(n+m))   = tmysl(i)
      tz(mb+i-(n+m))   = tmzsl(i)
!      write(89,'(i6,10d15.7)') &
!             mb+i-(n+m),yc(mb+i-(n+m)),zc(mb+i-(m+n)), &
!             yb(mb+i-(n+m)),zb(mb+i-(m+n)),yf(mb+i-(n+m)),zf(mb+i-(m+n)), &
!             dphtb(mb+i-(m+n)), ph(mb+i-(m+n)), phb(mb+i-(m+n)), &
!            ry(mb+i-(n+m)),rz(mb+i-(n+m)),ty(mb+i-(n+m)),tz(mb+i-(n+m))
    endif
  enddo

! write(89,*)
! write(89,*)
! write(82,*)
! write(82,*)
  write(87,*)
  write(87,*)









! - recalculate the amplitudes by using the new indexation of the boundary
!   panels

  do i=1,ntt
    am(i)=sqrt( (yf(i)-yb(i))**2 + (zf(i)-zb(i))**2 )
  enddo






! - Jet modelling
  write(75,*) '# jt ',jt
  write(76,*) '# jt ',jt

  i0 = ng-(m+n)      ! why ? isnt ng very close to m+n ? 

                     ! INDEED i0 should be always 0 by definition of n
                     ! as (n=ng-m)

  do i=i0-1,ng     ! so, i always starts from -1

    if(i.le.i0)then   ! this is entered when i=-1 and i=0
                      ! When i=-1, we have ii=mb-1 and jj = mb+2
                      ! When i=0, ii=mb and jj=mb+1
                      ! something need to be clarified yet!!!
      ii = mb+i-i0
      jj = mb-(i-i0)+1
    elseif(i.ge.i0+1.and.i.le.i0+m)then   ! Here i ranges from 1 to m
                            ! and thus ii=mt+i and jj=nt-i+1
      ii = mt+i-i0     ! SEE part 3 of the sketch (whatsapp)
      jj = nt-(i-i0)+1  ! SEE part 4 of the sketch (whatsapp)

    else      ! Here i ranges from m+1 to m+n (=ng!) 
              ! and then ii = nt-m+i and jj=ntt+m+1-i
      ii = nt+(i-i0-m)  ! SEE part 5 of the sketch (whatsapp)
      jj = ntt-(i-i0-m)+1   ! SEE part 6 of the sketch (whatsapp)
    endif






! here below the coordinates \xi and \zeta are assigned on
! the body (xi(ii),ze(ii)) and on the free surface (xi(jj),ze(jj)) 
! sides of the modelled part of the jet
! Similarly, the coordinates of the center of the control volume
! (xis(ii),zis(ii)) are also assigned

! According to the notation used in the "get" routine, xigb(0)
! corresponds to xcn(npc), xigb(1) to xcn(npc+1)  and so forth. 
! The last statement needs to be better verified

    xi(ii)  = xigb(i)
    xis(ii) = xigs(i)
    ze(ii)  = zegb(i)
    zes(ii) = zegs(i)
!   hp(jj)  = dh(i)

    xi(jj)  = xigf(i)
    xis(jj) = xigs(i)
    ze(jj)  = zegf(i)
    zes(jj) = zegs(i)
    write(75,'(2i6,4d15.5)') ii,i,xi(ii),xis(ii), &
                            ze(ii),zes(ii)
    write(76,'(2i6,4d15.5)') jj,i,xi(jj),xis(jj), &
                            ze(jj),zes(jj)
  enddo
  write(75,*)
  write(75,*)
  write(75,*)
  write(76,*)
  write(76,*)
! - fattore di matching rl non va con n=0 m>1







!!!! ----------------- Things are more or less clear up to this point


! Here below the coefficient l_i used in equations 14 and 15 of the JEM
! paper are initialized. The coefficient varies linearly and it is 1 at
! the matching with the full BEM region and 0 at the end of the transition
! region. The variation is assigned based on the distance rather than on
! the panel number

  do i=1,nt
    rl(i)=0.d0
  enddo

  do i=nt+1,ntt   ! rl is 1 on portions 5 and 6 of the sketch
    rl(i)=1.d0
  enddo


! here rl is assigned on portions 3 and 4. As a first step rl is
! assigned as the distance from the beginning of the modelled part of the
! jet until the end of the transition region. 
! Then, the total length of the jet amt is computed and rl is divided
! by amt to assure that rl will range from 0 to 1 in the transition region

  do i=1,m  
    ii  = mt+i
    if(i.eq.1)then
      dam = 0.5d0*(am(ii)+am(mb))   ! dam is the distance between two
            ! successives panel centroids or vertices of the elements
    else
      dam = 0.5d0*(am(ii)+am(ii-1))
    endif
    rl(ii)   = rl(ii-1) + dam
  enddo

! When n=0 then rl reaches 1 at the tip of the jet (which is the last
! element) and thus dam needs to be evaluated differently
! otherwise the same procedure as above applies

  if(n.eq.0)then
    dam = am(mt+m)
  else
    dam = 0.5d0*( am(mt+m) + am(nt+1) )
  endif

! Here below rl is scaled by the length of the transition region amt

  amt = rl(ii) + dam
  write(67,*) ' # jt ',jt
  do i=1,m
    rl(mt+i)   = rl(mt+i)/amt

! Once the value on the body side is computed, the same value is
! assigned on the corresponding panel on the free surface side
    rl(nt-i+1) = rl(mt+i)
    write(67,'(3i6,d15.6)') i,mt+i,nt-i+1,rl(mt+i)
  enddo
  write(67,*)
  write(67,*)





















! This is the begning of the important bit, 
! where the coupled BEM-FEM influence matrix aa is actually built
! Part 1 - initialize  the matrix with zero entries ------------------------------------------------------------------------
! note that an auxiliary bb vector, ff, is also initialized
! also note the number of equations (npe) is defined here
! remember the index definitions in the begnning,
! so there are mb+mf+2*m equation for the bulk (BEM + transition region) 
! plus 5*(n+m) unknow coefficients for the jet part 

! An additional equation would be necessary in case the FF condition is
! exploited. Such additional equation is not considered at the moment in
! this verion of the code


  npe = nt+5*(n+m)
  do i=1,npe
    ff(i) = 0.d0
    bb(i) = 0.d0
    do j=1,npe
      aa(i,j) = 0.d0
    enddo
  enddo

! Part 2 - construct the BEM influence matrix, lines 1 to nt (which is
! part 1 to 4 of the sketch)
! here i is the row index whereas k is the column (or unknown ) index

! As to the unknowns: on the wetted body surface phi is unknown and phi is
! also considered as unknown on part 3 (i.e. the transition region where 
! k is between mb+1 and mb+m) for all panels which are not
! separated from the solid contour. For the separated part the velocity
! potential is assigned and the normal derivative is unknown (this depends
! on kse)
! On the free surface side (both in zone 2 and 4) phi is assigned and
! dphi is unknown

  do i = 1,nt   ! define the rows of the linear system on zone 1 to 4
    yy = yc(i)    ! collocation point located and the panel midpoint of
                  ! the BEM discretization. The coordinates yc and zc
                  ! have been assigned at the beginning of the routine
    zz = zc(i)

    do k = 1,ntt  ! span over all panels from zone 1 to 6. 
           ! NOTE that in the modelled part of the jet a contribution 
           ! from the FEM discretization should also be included (which
           ! is partial in the transition region). The way in which that
           ! contribution is included is explained in the following

      yp = yb(k)
      ys = yf(k)
      zp = zb(k)
      zs = zf(k)
      yps= -yf(k)
      yss= -yb(k)
      zps= zf(k)
      zss= zb(k)

! - calculate the influence coefficient of panel k including the 
!   symmetry (by image method about y=0)

     call finte(yy,zz,yp,zp,ys,zs,am(k),fint1,fint2,0.d0)
     call finte(yy,zz,yps,zps,yss,zss,am(k),fins1,fins2,0.d0)
     gik  =  (fint1+fins1)/(2.d0*pi)    ! this is the integral of green
                 ! funtion G
     dgik = -(fint2+fins2)/(2.d0*pi)    ! this is the integral of normal
                 ! derivative of the green \partial G/\partial n

! Zone 1 of the sketch (panels on bulk region on the body side)
     if(k.le.mb)then


! to make the code clear, I propose to rewrite (assuming kse(k) is boolean) 
! Perfect!! 

        aa(i,k) = aa(i,k) + dgik*(1-kse(k)) - gik*kse(k)
        ff(i)   = ff(i) + dphn(k)*gik*(1-kse(k)) - ph(k)*dgik*kse(k)

! Zone 2 (panels on the bulk region on the free surface side)

     elseif(k.ge.mb+1.and.k.le.mt)then
! free surface element
        aa(i,k) = aa(i,k) - gik
        ff(i)   = ff(i) - ph(k)*dgik

! Zone 3 (transition region - body side) 

     elseif(k.ge.mt+1.and.k.le.mt+m)then
        rlk  = rl(k)
        xik  = xi(k)
        xisk = xis(k)
        zek  = ze(k)
        zesk = zes(k)





! If panel k is not separated the dphi/dn is assigned whereas the
! velocity potential is unknown





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  IMPORTANT NOTE !!!!!

!! IT WOULD BE IMPORTANT TO CHECK THE NOTATION USED IN THE FINITE
!! ELEMENTS. IN FIGURE 2 OF THE JEM PAPER, THE ELEMENT i IS LIMITED BY
!! VERTICES i-1 AND i and it is not clear if a different assumption is
!! used here

        if(kse(k).ne.1)then

! If the velocity potential is unknown, we use equation 14 and than take
! the first term (i.e. the one multiplied by 1-l_i) and use that as the
! influence coefficient of the velocity potential of panel k (see
! definition of the influence coefficient a(i,k)

! effect of potential of panel k on panel i

            aa(i,k)                = aa(i,k)+ (1.d0-rlk)*dgik !/1000.d0 

! But then we have to add the term phi_i^J (multiplied by l_i). For
! phi_i^J we use equation 9 collocated at the centroid of panel k (which
! is located at xi(k)).
! This introduce the new unknowns (a_i, b_i, c_i, d_i, e_i). 
! When substituting equation 9 into the Green expansion, the coefficient
! A_i will be multiplied by dgik (which will be the corresponding
! coefficient). Similarly, b_i will have a coefficient dgik*(xik-xisk),
! the latter term coming directly from eq. 9 of the JEM paper, and so
! forth. Because of equation 14, all the coefficient will have to be
! multiplied by rlk (see below)

            aik = dgik
            bik = dgik*(xik-xisk) 
            cik = dgik*(zek-zesk) 
            dik = 0.5d0*dgik*( (xik-xisk)**2-(zek-zesk)**2 )  
            eik = dgik*(xik-xisk)*(zek-zesk)  


! The unknown a_i is positioned after the unknowns phi/dphi of the
! elemets in the transition region, and thus at k+2m 
! As the a_i are (m+n) unknowns, the unknowns b_i will be shifted by
! (n+m), c_i by 2*(n+m) and so forth

            aa(i,k+2*m)          = aa(i,k+2*m)         + rlk*aik
            aa(i,(n+m)+k+2*m)    = aa(i,(n+m)+k+2*m)   + rlk*bik
            aa(i,2*(n+m)+k+2*m)  = aa(i,2*(n+m)+k+2*m) + rlk*cik
            aa(i,3*(n+m)+k+2*m)  = aa(i,3*(n+m)+k+2*m) + rlk*dik
            aa(i,4*(n+m)+k+2*m)  = aa(i,4*(n+m)+k+2*m) + rlk*eik
            ff(i)  = ff(i) + dphn(k)*gik

        else     ! panel is detached from the body surface 
                 ! in this case the velocity potential is assigned
                 ! whereas the normal derivative is unknown

            ! Hence, the influence coefficient of the normal
            ! derivative of the velocity potential k on the unknown 
            ! variable of panel i will be

            aa(i,k) = aa(i,k)-(1.d0-rlk)*gik

            ! for the FEM representation, the normal derivative of
            ! equation 9 has to be computed. In this case A_i disappears 
            ! and the other terms are computed as (d/dy*ny+d/dz*nz)

            ryk  = ry(k)
            rzk  = rz(k)
            bik  = ryk*gik
            cik  = rzk*gik
            dik  = gik*( ryk*(xik-xisk) - rzk*(zek-zesk) )  
            eik  = gik*( ryk*(zek-zesk) + rzk*(xik-xisk) )   
            aa(i,(n+m)+k+2*m)    = aa(i,(n+m)+k+2*m)   - rlk*bik
            aa(i,2*(n+m)+k+2*m)  = aa(i,2*(n+m)+k+2*m) - rlk*cik
            aa(i,3*(n+m)+k+2*m)  = aa(i,3*(n+m)+k+2*m) - rlk*dik
            aa(i,4*(n+m)+k+2*m)  = aa(i,4*(n+m)+k+2*m) - rlk*eik
            ff(i)   = ff(i) - ph(k)*dgik 
        endif

! Zone 4
     elseif(k.ge.mt+m+1.and.k.le.nt)then

         ! In Zone 4 (i.e. transition region on the free surface side) 
         ! the influence coefficients should be computed as for Zone 3
         ! in case of separated panels

       kk   = mt+(nt-k)+1
       rlk  = rl(k)

!!!!!!!!!!!!!!!!!!!!   IMPORTANT
!!! NOTE: here xis(kk) is used whereas xis(k) is used in Zone 3. This is
!   simply because the midpoint of the Element is based on the indexation of
!   the jet panels taken on the body side. So, kk is the index of the
!   panel belonging to the same Element on the body side (there could be a
!   +/- 1 variation because of the orientation of the numbering)
!   It would be good to have a sketch for this

       xik  = xi(k)
       xiskk= xis(kk)
       zek  = ze(k)
       zeskk= zes(kk)
       ryk  = ry(k)
       rzk  = rz(k)
       bik  = ryk*gik
       cik  = rzk*gik
       dik  = gik*( ryk*(xik-xiskk) - rzk*(zek-zeskk) )
       eik  = gik*( ryk*(zek-zeskk) + rzk*(xik-xiskk) )
       aa(i,(n+m)+kk+2*m)    = aa(i,(n+m)+kk+2*m)   - rlk*bik
       aa(i,2*(n+m)+kk+2*m)  = aa(i,2*(n+m)+kk+2*m) - rlk*cik
       aa(i,3*(n+m)+kk+2*m)  = aa(i,3*(n+m)+kk+2*m) - rlk*dik
       aa(i,4*(n+m)+kk+2*m)  = aa(i,4*(n+m)+kk+2*m) - rlk*eik
       aa(i,k) = aa(i,k)-(1.d0-rlk)*gik
       ff(i)   = ff(i) - ph(k)*dgik

! Zone 5:
! In zone 5 the influence coefficients are the same as in Zone 3 but for
! the fact that now rlk is 1 and thus there is no direct effect of the
! potential (or normal derivative) of panel k on panel i and there are
! only the other effects related to the coefficients of the expansion 9


!! The indexation used in the FE needs to be better understood.
!! The coefficient a_i (as well as all other b_i, c_i ...) appear when
!! the FE panel on the body and on the FS side are involved. It is
! important to better understant the relation between the index ii on the
! body side and jj on the FS side. It should be very simple but a sketch
! would be of help

     elseif(k.ge.nt+1.and.k.le.nt+n)then

! - Zone 5

       xik  = xi(k)
       xisk = xis(k)
       zek  = ze(k)
       zesk = zes(k)
       if(kse(k).ne.1)then
         aik = dgik
         bik = dgik*(xik-xisk)
         cik = dgik*(zek-zesk)
         dik = 0.5d0*dgik*( (xik-xisk)**2-(zek-zesk)**2 )
         eik = dgik*(xik-xisk)*(zek-zesk)
         aa(i,k+m)          = aa(i,k+m)         + aik
         aa(i,k+(n+m)+m)    = aa(i,k+(n+m)+m)   + bik
         aa(i,k+2*(n+m)+m)  = aa(i,k+2*(n+m)+m) + cik
         aa(i,k+3*(n+m)+m)  = aa(i,k+3*(n+m)+m) + dik
         aa(i,k+4*(n+m)+m)  = aa(i,k+4*(n+m)+m) + eik
         ff(i)  = ff(i) + dphn(k)*gik

       else
         ryk  = ry(k)
         rzk  = rz(k)
         bik  = ryk*gik
         cik  = rzk*gik
         dik  = gik*( ryk*(xik-xisk) - rzk*(zek-zesk) )
         eik  = gik*( ryk*(zek-zesk) + rzk*(xik-xisk) )
         aa(i,k+(n+m)+m)  = aa(i,k+(n+m)+m)   - bik
         aa(i,k+2*(n+m)+m)= aa(i,k+2*(n+m)+m) - cik
         aa(i,k+3*(n+m)+m)= aa(i,k+3*(n+m)+m) - dik
         aa(i,k+4*(n+m)+m)= aa(i,k+4*(n+m)+m) - eik
         ff(i)   = ff(i) - ph(k)*dgik
       endif

! Zone 6
     elseif(k.ge.nt+n+1)then

! here again (as before for zone 4) the relation between kk and k has to
! better understood. But things are more or less clear anyhow

       kk   = nt+(ntt-k)+1
       xik  = xi(k)
       xiskk= xis(kk)
       zek  = ze(k)
       zeskk= zes(kk)
       ryk  = ry(k)
       rzk  = rz(k)
       bik  = ryk*gik
       cik  = rzk*gik
       dik  = gik*( ryk*(xik-xiskk) - rzk*(zek-zeskk) )
       eik  = gik*( ryk*(zek-zeskk) + rzk*(xik-xiskk) )
       aa(i,kk+(n+m)+m)  = aa(i,kk+(n+m)+m)   - bik
       aa(i,kk+2*(n+m)+m)= aa(i,kk+2*(n+m)+m) - cik
       aa(i,kk+3*(n+m)+m)= aa(i,kk+3*(n+m)+m) - dik
       aa(i,kk+4*(n+m)+m)= aa(i,kk+4*(n+m)+m) - eik
       ff(i)   = ff(i) - ph(k)*dgik
     endif
    enddo ! end of the loop in k


!-----------------   Arrived here


















! - self induced calculations

    if (i.le.mb) then
      if(kse(i).ne.1)then
        aa(i,i) = aa(i,i) + 0.5d0  !/1000.d0
      else
!   write(98,*) 'ecc5!',i ,jt,ph(i)
        ff(i)   = ff(i) - 0.5d0*ph(i)
      endif

    elseif(i.ge.mb+1.and.i.le.mt)then
          ff(i)   = ff(i) - 0.5d0*ph(i)
    elseif(i.ge.mt+1.and.i.le.mt+m)then

    if(kse(i).ne.1)then
      rli = rl(i)
      xii = xi(i)
      xisi= xis(i)
      zei = ze(i)
      zesi= zes(i)
      aa(i,i)             = aa(i,i) + (1.d0-rli)*0.5d0  !/1000.d0
      aa(i,nt+i-mt)       = aa(i,nt+i-mt)  +  rli*0.5d0
      aa(i,nt+(m+n)+i-mt) = aa(i,nt+(m+n)+i-mt) + rli*0.5d0*(xii-xisi)
      aa(i,nt+2*(m+n)+i-mt)= aa(i,nt+2*(m+n)+i-mt) + rli*0.5d0*(zei-zesi)
      aa(i,nt+3*(m+n)+i-mt)= aa(i,nt+3*(m+n)+i-mt)+rli*0.25d0*((xii-xisi)**2-(zei-zesi)**2)
      aa(i,nt+4*(m+n)+i-mt)=aa(i,nt+4*(m+n)+i-mt) + rli*0.5d0*(xii-xisi)*(zei-zesi)
    else
!     write(98,*) 'ecc6!',i ,jt,ph(i)
      ff(i)   = ff(i) - 0.5d0*ph(i)
    endif
    elseif(i.ge.mt+m+1.and.i.le.nt)then
      ff(i)   = ff(i) - 0.5d0*ph(i)
    endif
  enddo  ! finishes the loop BEM in i, BEM influenced matrix is built



!Part 3: The Jet FEM model , lines from nt+1,nt+5(m+n)
! Attention: i here is not the line index itself, but the control volume index (V_i)
!i1: is the first panel on the body side 
!i : second panel on the body side which is the same index number as the control volume (V_i) 
!j : index of the second panel on the free surface side
!j2 : index of the first panel on the free surface side
!i2: (TO UNDERSTAND BETTER), seems it is used only for debugging purposes in the dxib



!First step: find the the corresponding vertex indexes
  do 100, i = mt+1, nt+n

! Region 1: Transition region in the body from mt+1 to mt+m
    if(i.ge.mt+1.and.i.le.mt+m)then
! If I got this right, it finds the corresponding index  of ith FEM element in the transition
! region of the body
      if(i.eq.mt+1)then
        i1 = mb
        i2 = mb-1
        j  = mb+1
        j1 = nt
      elseif(i.eq.mt+2)then
        i1 = i-1
        i2 = mb
        j  = nt-(i-mt)+2
        j1 = j-1
      else
        i1 = i-1
        i2 = i-2
        j  = nt-(i-mt)+2
        j1 = j-1
      endif
      ii  = nt+i-mt
      iii = ii
      iij = ii
      mn  = m+n
      k1  = m

! Region 2: Jet modelled part on the body from mt+2m+1 to mt+2m+n
    elseif(i.ge.nt+1.and.i.le.nt+n)then

      if(i.eq.nt+1)then
        if(m.eq.0)then
          i1 = mb
          i2 = mb-1
          j  = mb+1
          j1 = ntt
        else
          i1 = mt+m
          i2 = mt+m-1
          j  = mt+m+1
          j1 = ntt
        endif

      elseif(i.eq.nt+2)then
        if(m.eq.0)then
          i1 = i-1
          i2 = mb
          j  = ntt-(i-nt)+2
          j1 = j-1
        else
          i1 = i-1
          i2 = mt+m
          j  = ntt-(i-nt)+2
          j1 = j-1
        endif
      else
        i1 = i-1
        i2 = i-2
        j  = ntt-(i-nt)+2
        j1 = j-1
      endif
      ii  = nt+i-nt
      iii = ii+5*m
      iij = ii+m
      mn  = m+n
      k1  = n
    else
      goto 100
    endif

! Second Step: Indexes are found can load the temp variables
! and calculate the influence matrix coefficients of the 
! full modelled part

    rli  = rl(i)
    rli1 = rl(i1)
    drl  = rli1-rli
    dxif = (xi(i)-xi(i1))     !*1000.d0
    dxib = (xi(i1)-xi(i2))    !*1000.d0

! xii: x coordinate of the ith point on the body
    xii  = xi(i)
! xii1: x coordinate of the (i-1)th point on the body
    xii1 = xi(i1)
! xisi: x coordinate of the ith control volume center point 
    xisi = xis(i)
! xisi1: x  coordinate of the (i-1)th  control volume point, needed to enforce 
! matching between adjacent elements
    xisi1= xis(i1)    
    zei  = ze(i)
    zei1 = ze(i1) 
! zei: z coordinate of the ith point on the body
    zesi = zes(i)
! zesi1: x  coordinate of the (i-1)th  control volume point, needed to enforce 
! matching between adjacent elements
    zesi1= zes(i1)
! same as above but for the free surface
    xij  = xi(j)
    xij1 = xi(j1)
    zej  = ze(j)
    zej1 = ze(j1)
!        sqh  = sqrt(1.d0+hp(i1)**2)
!        hpi1 = hp(i1)
    ryj  = ry(j)
    rzj  = rz(j)
    ryi1 = ry(i1)
    rzi1 = rz(i1)
    ryi  = ry(i)
    rzi  = rz(i)

! - phi on the free surface  (j)
    aa(iii+k1,iij)       = aa(iii+k1,iij)       + 1.d0
    aa(iii+k1,iij+mn)    = aa(iii+k1,iij+mn)    + (xij-xisi)
    aa(iii+k1,iij+2*mn)  = aa(iii+k1,iij+2*mn)  + (zej-zesi)
    aa(iii+k1,iij+3*mn)  = aa(iii+k1,iij+3*mn)  + 0.5d0*((xij-xisi)**2-(zej-zesi)**2)
    aa(iii+k1,iij+4*mn)  = aa(iii+k1,iij+4*mn)  + (xij-xisi)*(zej-zesi)
    ff(iii+k1)  = ff(iii+k1) + ph(j)

! - phi on the free surface (j1)
    aa(iii+4*k1,iij)       = aa(iii+4*k1,iij)       + 1.d0
    aa(iii+4*k1,iij+mn)    = aa(iii+4*k1,iij+mn)    + (xij1-xisi)
    aa(iii+4*k1,iij+2*mn)  = aa(iii+4*k1,iij+2*mn)  + (zej1-zesi)
    aa(iii+4*k1,iij+3*mn)  = aa(iii+4*k1,iij+3*mn)  +0.5d0*((xij1-xisi)**2-(zej1-zesi)**2)
    aa(iii+4*k1,iij+4*mn)  = aa(iii+4*k1,iij+4*mn)  +(xij1-xisi)*(zej1-zesi)
    ff(iii+4*k1)  = ff(iii+4*k1) + ph(j1)

! TO DO: find out what kor(i) actually is !

    if(kor(i).eq.2)then

! - continuita vn sulla SL (j)
      aa(iii,j)         = aa(iii,j)         + drl
      aa(iii,iij+mn)    = aa(iii,iij+mn)    + rli*ryj
      aa(iii,iij-1+mn)  = aa(iii,iij-1+mn)  - rli1*ryj
      aa(iii,iij+2*mn)  = aa(iii,iij+2*mn)  + rli*rzj
      aa(iii,iij-1+2*mn)= aa(iii,iij-1+2*mn)- rli1*rzj
      aa(iii,iij+3*mn)  = aa(iii,iij+3*mn) +rli*( ryj*(xij-xisi) - rzj*(zej-zesi) )
      aa(iii,iij-1+3*mn)= aa(iii,iij-1+3*mn)-rli1*( ryj*(xij-xisi1) - rzj*(zej-zesi1) )
      aa(iii,iij+4*mn)  = aa(iii,iij+4*mn) +rli*( ryj*(zej-zesi) + rzj*(xij-xisi) )
      aa(iii,iij-1+4*mn)= aa(iii,iij-1+4*mn)- rli1*( ryj*(zej-zesi1)+ rzj*(xij-xisi1) )
      ff(iii)       = ff(iii)      + 0.0d0
! - vn sul corpo (i1) o phi sulla SL (se separa)
      if(kse(i1).ne.1)then
        aa(iii+2*k1,iij+mn)   = aa(iii+2*k1,iij+mn)    + ryi1 
        aa(iii+2*k1,iij+2*mn) = aa(iii+2*k1,iij+2*mn)  + rzi1 
        aa(iii+2*k1,iij+3*mn) = aa(iii+2*k1,iij+3*mn)  + ryi1*(xii1-xisi) - rzi1*(zei1-zesi)
        aa(iii+2*k1,iij+4*mn) = aa(iii+2*k1,iij+4*mn)  + ryi1*(zei1-zesi) + rzi1*(xii1-xisi)
        ff(iii+2*k1)  = ff(iii+2*k1) + dphn(i1)
      else
!        write(98,*) 'eccolo10!',i1,nt+n,jt
!        write(98,*)  xii1,zei1,ph(i1)
        aa(iii+2*k1,iij)       = aa(iii+2*k1,iij)       + 1.d0 
        aa(iii+2*k1,iij+mn)    = aa(iii+2*k1,iij+mn)    + (xii1-xisi)
        aa(iii+2*k1,iij+2*mn)  = aa(iii+2*k1,iij+2*mn)  + (zei1-zesi)
        aa(iii+2*k1,iij+3*mn)  = aa(iii+2*k1,iij+3*mn)  + 0.5d0*((xii1-xisi)**2-(zei1-zesi)**2)
        aa(iii+2*k1,iij+4*mn)  = aa(iii+2*k1,iij+4*mn)  + (xii1-xisi)*(zei1-zesi)
        ff(iii+2*k1)  = ff(iii+2*k1) + ph(i1)       
      endif
! - vn sul corpo (i)o phi sulla SL (se separa)
      if(kse(i).ne.1)then
        aa(iii+3*k1,iij+mn)   = aa(iii+3*k1,iij+mn)    + ryi 
        aa(iii+3*k1,iij+2*mn) = aa(iii+3*k1,iij+2*mn)  + rzi 
        aa(iii+3*k1,iij+3*mn) = aa(iii+3*k1,iij+3*mn)  + ryi*(xii-xisi) - rzi*(zei-zesi)
        aa(iii+3*k1,iij+4*mn) = aa(iii+3*k1,iij+4*mn)  + ryi*(zei-zesi) + rzi*(xii-xisi)
        ff(iii+3*k1)  = ff(iii+3*k1) + dphn(i)
      else
!        write(98,*) 'eccolo11!',i,nt+n,jt
!        write(98,*) xii,zei,ph(i)
        aa(iii+3*k1,iij)       = aa(iii+3*k1,iij)       + 1.d0 
        aa(iii+3*k1,iij+mn)    = aa(iii+3*k1,iij+mn)    + (xii-xisi)
        aa(iii+3*k1,iij+2*mn)  = aa(iii+3*k1,iij+2*mn)  + (zei-zesi)
        aa(iii+3*k1,iij+3*mn)  = aa(iii+3*k1,iij+3*mn)  + 0.5d0*((xii-xisi)**2-(zei-zesi)**2)
        aa(iii+3*k1,iij+4*mn)  = aa(iii+3*k1,iij+4*mn)  + (xii-xisi)*(zei-zesi)
        ff(iii+3*k1)  = ff(iii+3*k1) + ph(i)       
      endif

    else

! - vn sul corpo (i1) o continuita di phi sul corpo (se separa)
      if(kse(i1).ne.1)then
        aa(iii,iij+mn)   = aa(iii,iij+mn)    + ryi1 
        aa(iii,iij+2*mn) = aa(iii,iij+2*mn)  + rzi1 
        aa(iii,iij+3*mn) = aa(iii,iij+3*mn)  + ryi1*(xii1-xisi) - rzi1*(zei1-zesi)
        aa(iii,iij+4*mn) = aa(iii,iij+4*mn)  + ryi1*(zei1-zesi) + rzi1*(xii1-xisi)
        ff(iii)  = ff(iii) + dphn(i1)
      else
        aa(iii,iij)        = aa(iii,iij)        + rli  
        aa(iii,iij-1)      = aa(iii,iij-1)      - rli1 
        aa(iii,iij+mn)     = aa(iii,iij+mn)     + rli*(xii1-xisi) 
        aa(iii,iij-1+mn)   = aa(iii,iij-1+mn)   - rli1*(xii1-xisi1) 
        aa(iii,iij+2*mn)   = aa(iii,iij+2*mn)   + rli*(zei1-zesi) 
        aa(iii,iij-1+2*mn) = aa(iii,iij-1+2*mn) - rli1*(zei1-zesi1) 
        aa(iii,iij+3*mn)   = aa(iii,iij+3*mn)   + 0.5d0*rli*((xii1-xisi)**2-(zei1-zesi)**2)
        aa(iii,iij-1+3*mn) = aa(iii,iij-1+3*mn) - 0.5d0*rli1*((xii1-xisi1)**2-(zei1-zesi1)**2)
        aa(iii,iij+4*mn)   = aa(iii,iij+4*mn)   + rli*(xii1-xisi)*(zei1-zesi)
        aa(iii,iij-1+4*mn) = aa(iii,iij-1+4*mn) - rli1*(xii1-xisi1)*(zei1-zesi1)
        ff(iii)  = ff(iii) - drl*ph(i1) 
      endif
! - nullita dei coefficienti d ed e
      aa(iii+2*k1,iij+3*mn) = aa(iii+2*k1,iij+3*mn)  + 1.d0
      aa(iii+3*k1,iij+4*mn) = aa(iii+3*k1,iij+4*mn)  + 1.d0
      ff(iii+2*k1)  = ff(iii+2*k1) + 0.d0
      ff(iii+3*k1)  = ff(iii+3*k1) + 0.d0
    endif
  100  continue ! end the do that started up there

  do i=1,npe
    bb(i) = ff(i)
  enddo

!     -------------   SOLUZIONE SISTEMA LINEARE CON LAPACK -------------

  call DGETRF(npe,npe,aa,npamx,indx,idum)
  trans='N'
!      rmax = 0.d0
!      rmin = 1.d31
!      do i=1,npe
!        rmax = max(rmax,abs(aa(i,i)))
!        rmin = min(rmin,abs(aa(i,i)))
!      enddo
!      write(*,*) 'rr ',rmax,rmin
  call DGETRS(trans,npe,1,aa,npamx,indx,bb,npamx,idum)

!     -------------   FINE SOLUZIONE SISTEMA LINEARE       -------------

!    risistemo arrays dopo soluzione sistema lineare
  do i=1,mb
    if(kse(i).ne.1)then
      ph(i) = bb(i)
      phb(i) = bb(i)
    else
      dphn(i) = bb(i)
      dphnb(i) = bb(i)
    endif
!        write(76,'(i6,4d15.7)') i,yc(i),zc(i),ph(i),dphn(i)
  enddo

!        write(76,*)
!Free surface solution dphn
  do i=mb+1,mb+mf
    dphn(i) = bb(i)
    dphnb(i) = bb(i)
!        write(77,'(i6,4d15.7)') i,yc(i),zc(i),dphn(i),ph(i)
  enddo

! Body transition region 
  do i=mt+1,mt+m
    if(kse(i).ne.1)then
      phb(i) = bb(i)
    else
      dphn(i) = bb(i)
      dphnb(i) = bb(i)
    endif
!        write(76,'(i6,4d15.7)') i,yc(i),zc(i),phb(i),dphn(i)
  enddo
!        write(77,*)
  do i=mt+m+1,nt
    dphnb(i) = bb(i)
!        write(77,'(i6,4d15.7)') i,yc(i),zc(i),dphnb(i),ph(i)
  enddo


!      write(76,*)
!      write(76,*)
!      write(77,*)
!      write(77,*)
!      do i =nt+1,nt+(m+n)
!        ii = i-nt
!        a(ii) = bb(i)
!        b(ii) = bb(i+m+n)
!        c(ii) = bb(i+2*(m+n))
!        d(ii) = bb(i+3*(m+n))
!        e(ii) = bb(i+4*(m+n))
!        write(80,'(i6,5d25.17)') ii,a(ii),b(ii),c(ii),d(ii),e(ii)
!      enddo
!      write(88,*) '# jt ',jt
  do i = mt+1,mt+m
    j  = nt-(i-mt)+1
    ii = i-mt
    a(i) = bb(nt+ii)
    b(i) = bb(nt+(m+n)+ii)
    c(i) = bb(nt+2*(m+n)+ii)
    d(i) = bb(nt+3*(m+n)+ii)
    e(i) = bb(nt+4*(m+n)+ii)
    a(j) = a(i)
    b(j) = b(i)
    c(j) = c(i)
    d(j) = d(i)
    e(j) = e(i)
    write(88,'(i6,5d15.7)') i,a(i),b(i),c(i),d(i),e(i)
  enddo
! write(88,*)
  do i = nt+1,nt+n
    j  = ntt-(i-nt)+1
    ii = i-nt+m
    a(i) = bb(nt+ii)
    b(i) = bb(nt+(m+n)+ii)
    c(i) = bb(nt+2*(m+n)+ii)
    d(i) = bb(nt+3*(m+n)+ii)
    e(i) = bb(nt+4*(m+n)+ii)
    a(j) = a(i)
    b(j) = b(i)
    c(j) = c(i)
    d(j) = d(i)
    e(j) = e(i)
    write(88,'(i6,5d15.7)') i,a(i),b(i),c(i),d(i),e(i)
  enddo
  write(88,*)
  write(88,*)
      
! riposiziono negli array d uscita
!      do i = 1,mb+m
!        if(i.le.mb)then
!          phi(i) = ph(i)
!          phib(i) = ph(i)
!        elseif(i.ge.mb+1.and.i.le.mb+m)then
!          phib(i) = phb(mt+i-mb)
!        endif
!      enddo
!      do i = n+1,mf+m+n
!        if(i.ge.n+1.and.i.le.n+m)then
!          dphibsl(i) = dphnb(mt+m+(i-n))
!        else
!          dphibsl(i)  = dphn(mb+i-(n+m))
!          dphisl(i)  = dphn(mb+i-(n+m))
!        endif
!      enddo

!  scrivo un po di controlli
!      do i=1,m
!        if(i.eq.1)then
!         ii1 = mb
!         jj  = mt+1
!        else
!         ii1 = mt+i-1
!         jj  = nt-i+2
!        endif
!        ii = mt+i
!        dxi  = xi(ii1)-xis(ii)
!        dze  = ze(ii1)-zes(ii)
!        phig = a(ii)+b(ii)*dxi+c(ii)*dze+0.5d0*d(ii)*(dxi**2-dze**2)+
!     #         e(ii)*dxi*dze
!        write(92,'(4d15.7)') yc(ii1),zc(ii1),phb(ii1),phig
!        dxi  = xi(jj)-xis(ii)
!        dze  = ze(jj)-zes(ii)
!        phig = a(ii)+b(ii)*dxi+c(ii)*dze+0.5d0*d(ii)*(dxi**2-dze**2)+
!     #         e(ii)*dxi*dze
!        write(93,'(4d15.7)') yc(jj ),zc(jj ),ph(jj ),phig
!       enddo
!       do i=1,m+n
!        if(i.eq.1.and.m.eq.0)then
!         ii1 = mb
!         jj  = mb+1
!         jj1 = ntt
!         ii  = nt+i
!         i1  = i-1
!        elseif(i.eq.1.and.m.ne.0)then
!         ii1 = mb
!         jj  = mb+1
!         jj1 = nt
!         ii  = mt+i
!         i1  = i-1
!        elseif(i.ge.2.and.i.le.m)then
!         ii  = mt+i
!         ii1 = ii-1
!         jj  = nt-i+2
!         jj1 = jj-1
!         i1  = i-1
!        else
!         ii  = nt+i-m
!         ii1 = ii-1 
!         jj  = ntt-(i-m)+2
!         jj1 = jj-1 
!         i1  = i-1
!        endif
!        sqh=sqrt(1.d0+hp(ii)**2)
!        rlii = rl(ii)
!        a1 =  a(ii)
!        b1 =  b(ii)
!        c1 =  c(ii)
!        d1 =  d(ii)
!        e1 =  e(ii)
!        dxi = xi(ii)-xis(ii)
!        dze = ze(ii)-zes(ii)
!        phig = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!        vxig  = b1+d1*dxi+e1*dze
!        vzeg  = c1-d1*dze+e1*dxi
!        dxi = xi(jj1)-xis(ii)
!        dze = ze(jj1)-zes(ii)
!        phif = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!        vxif  = b1+d1*dxi+e1*dze
!        vzef  = c1-d1*dze+e1*dxi
!        vnf   = (vxif*hp(ii)-vzef)/sqh
!        vtf   = (-vxif-hp(ii)*vzef)/sqh
!        a1 =  a(ii+1)
!        b1 =  b(ii+1)
!        c1 =  c(ii+1)
!        d1 =  d(ii+1)
!        e1 =  e(ii+1)
!        dxi = xi(ii)-xis(ii+1)
!        dze = ze(ii)-zes(ii+1)
!        phigg = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!        vxigg = b1+d1*dxi+e1*dze
!        vzegg = c1-d1*dze+e1*dxi
!        dxi = xi(jj1)-xis(ii+1)
!        dze = ze(jj1)-zes(ii+1)
!        phiff = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!        vxiff = b1+d1*dxi+e1*dze
!        vzeff = c1-d1*dze+e1*dxi
!        vnff  = (vxiff*hp(ii)-vzeff)/sqh
!        vtff  = (-vxiff-hp(ii)*vzeff)/sqh
!        phii = (1.d0-rlii)*phb(ii)+rlii*phig
!        dphii= (1.d0-rlii)*dphnb(jj1)+rlii*vnff
!        write(62,'(6d15.6)')   yc(ii),zc(ii), phb(ii),phig,phigg,phii
!        write(64,'(8d15.7)') 
!     #      yc(ii),zc(ii), dphtb(ii),vxig,vxigg,dphn(ii),vzeg,vzegg
!        write(63,'(5d15.6)') yc(jj1),zc(jj1),ph(jj1),phif,phiff
!        write(65,'(9d15.7)') 
!     #      yc(jj1),zc(jj1), dphtb(jj1),vtf,vtff,dphnb(jj1),vnf,vnff,
!     #      dphii
!     #      yc(jj1),zc(jj1), dphtb(jj1),vxif,vxiff,dphn(jj1),vzef,vzeff
!      enddo
!      write(62,*)
!      write(62,*)
!      write(63,*)
!      write(63,*)
!      write(64,*)
!      write(64,*)
!      write(65,*)
!      write(65,*)
!            
!
  return
  end
