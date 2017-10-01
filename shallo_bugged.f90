
subroutine  shallo(ynsl,znsl,ygb,zgb,kget,jt,npsl, &
                    ng,ampli,rmg,iiget,jjget, &
                    ycnsl,zcnsl,phinsl,dphisl,epsg,epsgg,eskg,gfrac, &
                    ygn,zgn,ygs2,zgs2,tg,tgb,nbody,ngo1,nsep, &
                     ycb,zcb,tcb,nngo,jind,jfid,jend)



!jind,jfid,jend -> shifted indexes of the control strip 
!To be completed, but:
! on exit:  ng,nngo, ygb,zgb,tcb
       
 implicit none
 integer, parameter :: npamx = 3000
! common/costanti/pi 
 integer :: nsep,ng,ngo1, npc,jjget,jt,kget    
 integer :: npsl, nbody, nngo, jind,jfid,jend
 real*8, dimension(npamx) :: ygb,zgb,ycb,zcb,tg,tgb,tcb
 real*8, dimension(npamx) :: phin,dphin,ygn,zgn,ygs2,zgs2 
 real*8 :: epsg,epsgg,eskg,gfrac, ampli
  
! local variables:
 integer :: kk, i,kbi,j,ii,kb,kbb,ki,kkin,max_iters,iiget,nng
 real*8 :: y1,z1,y2,z2,w2y,w2z,ww,w1y,w1z,yi,zi,ti,si,epss,yo,yp1,ypn
 real*8 :: dy,dyo,ampg,amplio,ay,az,dr,rc,rmg,tb,tbb,tbi,wby,wbz,y,r 
 real*8 :: yb,ybb,ybi,yy,z,zb,zbb,zbi,zz
 real*8, dimension(npamx) :: ycnsl,zcnsl,ynsl,znsl,phinsl, dphisl

 
 write(*,*) '-->running  shallo......,ng===' , ng

! Parse the data

 nng  = ng*kget

! eskg is the growing factor for the panel size in the jet region


! the routine calculates the limiting distance to be used to remove the 
! panels which are too close.
! The limiting distance is based on ampli (which is computed in the
! following). The limiting distance is increased when the jet is
! separated.

! Define iiget as a function of jjget 
! iiget = jjget-4
 
!!              ---- > TO BE UNDERSTOOD: what is epsgg??
!!              ---- > TO BE UNDERSTOOD: why epsg is 3 times in case of flow
!!                                separation

!espgg : inferior limit of free surface panels close to the body (amplitude fraction, of the 10**-3)

 if(nsep.eq.0)then    
   epsg = epsgg*ampli
 else
   epsg = 3.d0*epsgg*ampli
 endif

 epss = 1.d-10    ! TO BE UNDERSTOOD!! it seems something to avoid
                  ! problems with machine precision

! epss : is used in a if statement down below in order to make the segment a little bigger  


 nngo =  kget*ng    ! nngo indicates the number of panels in the jet
                    ! when the jet model is activated (kget=1)

 ! if the flow is not separated, kk is initialize to be used in the
 ! following.

 if(nsep.lt.1)then  
   kk=0  

! ------------> DA CAPIRE: it is not clear what kbi represents and where
!               it is used

   kbi  = 0 
    
   do i=1,ng     ! sweep up to ng (which is assignet in the following by
                 ! a previous call to the routine splver2 ????)

! centroids in the free surface within the jet zone
! WARNING: differently from the previous version of the code (here npc
! does not include the panels in the modelled part of the jet)
! we want to access the first free surface panels


     y1  = ycnsl(i) ! ycnsl(i) in old code 
     z1  = zcnsl(i)
     
!! WARNING: the part below has to be rearranged yet. Even though it
!! works, it is TO BE UNDERSTOOD what NGO1 indicates. It seems NGO1
!! is essentially NBODY before flow separation but when there is flow
!! separation it is not clear how the spline of the body side is modified.

!! It seems that NGO1 is used for an EXTENDED representation of the body
!! spline when flow separation occur. However that needs TO BE UNDERSTOOD

     do j=1,ngo1-1          ! chiarire la questione ngo1
       y2   = ygn(j)     ! segmento del corpo base 
       z2   = zgn(j)
       w2y  = ygn(j+1)-ygn(j)   ! DY of the body panel
       w2z  = zgn(j+1)-zgn(j)
       ww   = sqrt(w2y**2+w2z**2)
       w1y  = -w2z/ww     ! y-component of normal to body panel
                          ! (basic representation)
       w1z  =  w2y/ww

! localizzo l'intersersezione tra la retta ortogonale al corpo
!  passante per y1,z1 e il segmento w2y,w2z. In pratica identifico
! il segmento del corpo si trova in corrispondenza del punto y1,z1
! (questo avviene solo se ti è tra 0 e 1). NOTA: questa cosa potrebbe
! dare problemi nel caso di corpi convessi/concavi con una
! discretizzazione di base non adeguata

       call intret(y1,z1,w1y,w1z,y2,z2,w2y,w2z,yi,zi,si,ti)
       if(ti.ge.0.d0.and.ti.le.1.d0)then
 ! se il punto è all'interno del pannello individuato e la sua distanza
 ! normale "si" è minore di epsg, assegno a kk l'indice di quel centroide
          if(abs(si).lt.epsg) kk=i
       endif
     enddo
   enddo

! se succede che un centroide viene a trovarsi ad una distanza dal corpo
! minore di epsg, tolgo tutti i pannelli di superficie libera anvanti 
! (e li sostituisco con uno solo... se ho capito bene)



!DEbugging
!   if(jt==100) then
!     write(*,*) 'si, epsg,kk=', si,epsg,kk    
!     stop
!   endif 

   if(kk.gt.0)then
     write(*,*) ' AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa '
     write(*,*) 'shallo 0, butto ',kk,' pannelli from the free surface'
     !stop
     do i=kk+1,npsl
        ycnsl(i-kk) = ycnsl(i)
        zcnsl(i-kk) = zcnsl(i)
        phinsl(i-kk)= phinsl(i)
        dphisl(i-kk)= dphisl(i)
        ynsl(i-kk)  = ynsl(i)
        znsl(i-kk)  = znsl(i)
!       ycn(npc+ng+i-kk) = ycn(npc+ng+i)
!       zcn(npc+ng+i-kk) = zcn(npc+ng+i)
!       phin(npc+ng+i-kk)= phin(npc+ng+i)
!       dphin(npc+ng+i-kk)= dphin(npc+ng+i) 
!       yn(npc+ng+i-kk)  = yn(npc+ng+i) 
!       zn(npc+ng+i-kk)  = zn(npc+ng+i)
     enddo
!     yn(npc+ng+npsl+1-kk)  = yn(npc+ng+npsl+1) 
!     zn(npc+ng+npsl+1-kk)  = zn(npc+ng+npsl+1)
      ynsl(npsl+1-kk)  = ynsl(npsl+1)
      znsl(npsl+1-kk)  = znsl(npsl+1)


! --  shift indici striscia di controllo
      if(jind.lt.kk+1.and.jend.eq.0)then
         write(*,*)'attenzione, la striscia se ne va!!!'
         write(*,*) 'jend shallo = 1'
         jend=1
       elseif(jend.eq.0)then
         jind=jind-kk
         jfid=jfid-kk
       endif
! end of control slip procedure 
   endif
   
   npsl = npsl - kk

   


! NOTA !!!
! alla fine della ridistribuzione, il vertice yn(npc+1) sarà il primo
! vertice del pannello kk che è quello con centroide vicino al corpo
! At this point panels have been deleted and the new intersection point 
! needs to be found. 


!   y1 = yn(npc+ng+1)  
!   z1 = zn(npc+ng+1)
    y1 = ynsl(1)
    z1 = znsl(1)

   kkin = 9999
   j = 1
   max_iters = ngo1-1

! scorro i nodi base della geometria del corpo
   do j=1,ngo1-1
     y2   = ygn(j)
     z2   = zgn(j)
! definisco il segmento del corpo
     w2y  = ygn(j+1)-ygn(j)
     w2z  = zgn(j+1)-zgn(j)
     ww   = sqrt(w2y**2+w2z**2)
! normale al corpo passante per il punto dove si taglia 
! dove e li si vanno a  mettere un certo numero minimo (nb) pannelli 
! sul truncation
     w1y  = -w2z/ww
     w1z  =  w2y/ww
! locate the intersection with the body
    
     call intret(y1,z1,w1y,w1z,y2,z2,w2y,w2z,yi,zi,si,ti)
   
     if(ti.ge.(0.d0-epss).and.ti.le.(1.d0+epss))then
       kb = j
       yb = yi
       zb = zi
       ay = yi-ygn(kb)
       az = zi-zgn(kb)
       tb = tg(kb)+sqrt(ay**2+az**2)
       kkin = 0
       exit
     endif
   enddo

!!!   si aggiunge epss per evitare problemi con la precisione e fare il
!segmento leggermente più ampio
 
 else   ! entro qui se c'è separazione e.g nsep = 1 
           ! non capisco bene in cosa si differenzia dal precedente
   write(*,*) 'epsg nsep 1',epsg,epsgg,ampli,nsep
   kk   = 0
   kbi  = 0
   do i=1,ng
! we are using the old code convention here , the first panel on the free surface 
! is (or should be) npc+ng+1
    ! y1  = ycn(npc+ng+i)
    ! z1  = zcn(npc+ng+i)
       y1  = ycnsl(i)
       z1  = zcnsl(i) 
       do j=1,ngo1-1
       y2   = ygn(j)
       z2   = zgn(j)
       w2y  = ygn(j+1)-ygn(j)
       w2z  = zgn(j+1)-zgn(j)
       ww   = sqrt(w2y**2+w2z**2)
       w1y  = -w2z/ww
       w1z  =  w2y/ww
       call intret(y1,z1,w1y,w1z,y2,z2,w2y,w2z,yi,zi,si,ti)
       if(ti.ge.0.d0.and.ti.le.1.d0)then
         if(abs(si).lt.epsg) kk=i
       endif
     enddo
   enddo

 ! Fino a qui ha fatto la stessa identica cosa del caso senza
 ! separazione. Probabilmente l'estensione nbody---> ngo1 include la parte
 ! separata (ma ancora non mi è ben chiaro. Dovrei farmi stampare un
 ! esempio in un caso concreto... magari lo faccio


   if(kk.gt.0)then
     write(*,*) 'shallo 1, butto ',kk,' pannelli'

 ! Come sopra rimuovo i pannelli nella parte del tip oltre ynsl(kk)
 ! quel punto diventerà ynsl(1)
!     do i=kk+1,npsl
!       ycn(npc+i-kk) = ycn(npc+i)
!       zcn(npc+i-kk) = zcn(npc+i)
!       phin(npc+i-kk)= phin(npc+i)
!       dphin(npc+i-kk)= dphin(npc+i) 
!       yn(npc+i-kk)  = yn(npc+i) 
!       zn(npc+i-kk)  = zn(npc+i)
!     enddo
!     yn(npc+npsl+1-kk)  = yn(npc+npsl+1) 
!     zn(npc+npsl+1-kk)  = zn(npc+npsl+1)
!   endif

!   npsl = npsl - kk
!   ngo1 = ngo1 - kk


     do i=kk+1,npsl
       ycnsl(i-kk) = ycnsl(i)
       zcnsl(i-kk) = zcnsl(i)
       phinsl(i-kk)= phinsl(i)
       dphisl(i-kk)= dphisl(i)
       ynsl(i-kk)  = ynsl(i)
       znsl(i-kk)  = znsl(i)
     enddo
       ynsl(npsl+1-kk)  = ynsl(npsl+1)
       znsl(npsl+1-kk)  = znsl(npsl+1)
   endif
   npsl = npsl - kk
   ngo1 = ngo1 - kk



!!!! Vedi che in questo caso toglie i kk pannelli anche a ngo1. Però è
!! strano... come puoi essere sicuro che non scendi sotto ngo (che è un
!! dato di ingresso???)
! - punto di intersezione corpo SL

!   y1 = yn(npc+1)
!   z1 = zn(npc+1)

   y1 = ynsl(1)
   z1 = znsl(1)
   kkin = 9999
   j = ngo1-1
   y2   = ygn(j)
   z2   = zgn(j)
   w2y  = ygn(j+1)-ygn(j)
   w2z  = zgn(j+1)-zgn(j)
   ww   = sqrt(w2y**2+w2z**2)
   w1y  = -w2z/ww
   w1z  =  w2y/ww

!   come prima cerco l'interzesione della normale per il punto ynsl(1) 

   call intret(y1,z1,w1y,w1z,y2,z2,w2y,w2z,yi,zi,si,ti)

! ridefinisco la spline tenendo conto del ridotto numero di nodi di
! controllo

   kb = ngo1-1 
   yb = yi 
   zb = zi
   ay = yi-ygn(kb)
   az = zi-zgn(kb)
   tb = tg(kb)+sqrt(ay**2+az**2)
   tg(ngo1)  = tb
   ygn(ngo1) = yb
   zgn(ngo1) = zb
   kkin = 0
   call spline(ygn,zgn,ygs2,zgs2,tg,yp1,ypn,ngo1,npamx)
! For Debug purposes:
   write(*,*) 'punti usat i',j,nbody,ti
   write(*,*) ygn(j),zgn(j)
   write(*,*) ygn(j+1),zgn(j+1)
   write(*,*) ynsl(1),znsl(1)
   write(*,*) yb,zb
!end debug
 endif ! end the if of nsep==0 or 1 !

! A questo punto ha tagliato la parte che poteva dare origine a delle
! instabilità. Vediamo come assegna ng e kget



! -------------- > WARNING: there are several paramters in the part above
! which need to be better understood

! - Identification of the MODELLED PORTION OF THE JET: 
!   don't get inside for jt < iiget  (iiget is about 50-60)
 
 if(jt.gt.iiget)then 

! Note: in this version a different indexation is employed. Here npc
! indicates the number of panel on the BODY side (regardless if they are
! on the bulk or in the modelled part of the jet and independently if they
! are separated from the body or not)

   dyo = 1.d0

   do i=2,npsl

! compute dy of the free surface panel
!     dy = yn(npc+i+1)-yn(npc+i)
      dy = ynsl(i+1)-ynsl(i)


! ng is set as the panel after the minimum in y of the free surface profile

     if (dy.gt.0.d0.and.dyo.le.0.d0) then 
        ng=i+1   ! NOTE: this is intentionally different from the previous
               ! version of the code (it was i+1)
               ! see SKETCH 1
        write(*,*) 'shallo ---- ng,ngo1,nbody = ', ng,ngo1,nbody  
     endif     

     dyo = dy
   enddo 

! kget is set to 1 only when ng > 0 and jt large enough (larger than
! jjget)
! Honestly, it is not clear why iiget and jjget are needed. This is
! something deserving further understanding
   print*, 'shallo ----- jt,jjget', jt, jjget
   if (ng.gt.0.and.jt.gt.jjget) then 
       kget = 1
       print*,'shallo changes kget----------->', kget
   endif 

! - Discretization (or rediscretization if the jet was already
!   initialized) of the MODELLED PART of the jet

   if (ng.gt.0) then

! y1 and z1 are the coordinate of the vertex (on the free surface side) 
! just after the point of inversion (or y-minimum )
!     y1  = yn(npc+ng+1)
!     z1  = zn(npc+ng+1)

     y1  = ynsl(ng+1)
     z1  = znsl(ng+1)


! identify the panel on the body which is at the root of y1,z1. The
! curvilinear abscissa of the point (tbb) is approximated as the sum of
! the curvilinear abscissa at the previous vertex of the body panel plus
! the distance from the vertex to the point ybb,zbb

     do j=1,ngo1-1
       y2   = ygn(j)
       z2   = zgn(j)
       w2y  = ygn(j+1)-ygn(j)
       w2z  = zgn(j+1)-zgn(j)
       w1y  = -w2z
       w1z  =  w2y
       print*, 'shallo calls intret'
       call intret(y1,z1,w1y,w1z,y2,z2,w2y,w2z,yi,zi,si,ti)
       if(ti.ge.0.d0.and.ti.le.1.d0) then
         kbb = j    ! index of body panel (WARNING: not so if ngo1>ngo)
         ybb = yi
         zbb = zi
         ay  = yi-ygn(kbb)
         az  = zi-zgn(kbb)
         tbb = tg(kbb)+sqrt(ay**2+az**2)
         write(*,*) 'kbb,tbb,tbody(kbb),tb,kb, ti',kbb,tbb,tg(kbb),tb,kb,ti
         exit
       endif
     enddo
     
     ampg = gfrac*(tb-tbb)
     tbi  = tb-ampg
     call splint(tbi,ybi,zbi,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
     ygb(1) = ybi
     zgb(1) = zbi
     tgb(1) = tbi
     write(*,*) 'ampg,tb,tbb',ampg,tbb,tb

! the jet model is used in the body region with abscissa between tbi and
! tb. Before jet modelling, tb is the tip of the jet or the first vertex
! on the free surface)
! After the jet model is activated, it is the inversection of the normal
! to the body (or to the separated free surface)
! passing trough the truncated part  with the body
! (or with the separareted free surface) 

! from now on, the modelled part of the jet will be identified by
! ygb(i), zgb(i) and tgb(i) along the body surface


! - intersezione con la SL (per determinare ampli e quindi ng)
     do j=1,ngo1
       if(tg(j).gt.tbi)then
         kbi = j-1
         exit
       endif
     enddo

! kbi is the index of the body panel around the beginning opf the
! jet modelling region

     if(kbi.eq.0) then
       write(*,*) 'intersezione non trovata SJALLO'
       write(*,*) 'kbi,',kbi,tbi,ampg,tb
       stop
     endif
     
     wby = -(zgn(kbi+1)-zgn(kbi))
     wbz =   ygn(kbi+1)-ygn(kbi)

! look for the point of the free surface where to start the jet model




     do i=1,ng
!       y1  = yn(npc+i)
!       z1  = zn(npc+i)
!       w1y = yn(npc+i+1)-yn(npc+i)
!       w1z = zn(npc+i+1)-zn(npc+i)
       y1  = ynsl(i)
       z1  = znsl(i)
       w1y = ynsl(i+1)-ynsl(i)
       w1z = znsl(i+1)-znsl(i)
       call intret(y1,z1,w1y,w1z,ybi,zbi,wby,wbz,yi,zi,si,ti)

       if(si.ge.0.d0.and.si.le.1.d0) exit
     enddo

 ! ki is the index of the free surface panel where the jet model starts

     ki = i
     write(*,*) 'jet models starts at -------ki===',ki 
! yi,zi are the coordinate of the intersection of the normal to the
! body passing at ybi,zbi with the free surface 
! ay and az are the component of the segment connecting the body and
! free surface at the root of the jet

     ay = ybi-yi
     az = zbi-zi
     write(*,*) 'ng,ampli-------1! ',ng,ampli

! store in amplio the value of ampli at the previous iteration
     amplio = ampli
! ampli is the size of the amplitude of the root panel divided by the
! number of opanels that will be used at the root (rmg)

     ampli = sqrt(ay**2+az**2)/rmg
     write(*,*) 'ng,ampli-------2! ',ng,ampli

! to avoid an excessive variation, ampli is taken as an average betwwn
! the value at the previous iteration and that at present one

     ampli = 0.5*ampli + 0.5*amplio

! the parameter ampli is used to discretize the modelled part of the
! jet
! NOTE: in this case the panel size grows towards the tip






!- calcolo ng
     if (eskg.ne.1.d0) then
       ng = int( log( 1.d0+(eskg-1.d0)*ampg/ampli)/log(eskg) )
       write(*,*) 'eskg,ng,ampg----' ,eskg,ng,ampg
! compute the exact value of ampli for discretization with ng panels
       if(ng.le.1) ng = 2
       ampli = ampg*(1.d0-eskg)/(1.d0-eskg**ng)
     else
! uniform panel size
       ng   = int(ampg/ampli)
       if(ng.le.1) ng = 2
       ampg = float(ng)*ampli 
     endif
   
     write(*,*) 'ng,ampli-------3',ng,ampli
! - array pti corpo nel getto, vertici e centroidi
! Initialize the discretization in the modelled part of the jet

! dr is the initial panel amplitude and r is the initial abscissa along
! the body (or separated jet portion)

     dr=ampli
     r =tbi
     write(29,*) '# jt ',jt,nsep
     write(28,*) '# jt ',jt,nsep
     write(29,'(i4,3d15.7)') 1,ygb(1),zgb(1),tgb(1)
    
     do i = 1,ng
       r = r+dr
       rc = r-0.5d0*dr
! as in disuni, r is the coordinate of the next vertes, rc is the
! coordinate of the next centroid

! locate the vertex
       call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx)
! locate the centroid
       call splint(rc,yy,zz,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx)

! ygb is the coordinate on the body (or on the separated jet) of the
! vertex
! ycb is the coordinate of the centroi
       ygb(i+1) = y
       zgb(i+1) = z
       tgb(i+1) = r
       ycb(i)   = yy
       zcb(i)   = zz
       tcb(i)   = rc
       dr  = dr*eskg
! write jet vertices to 29  
     write(29,'(i4,3d15.7)') i+1,ygb(i+1),zgb(i+1),tgb(i+1)
! write jet centroids to 28
     write(28,'(i4,3d15.7)') i,ycb(i),zcb(i),tcb(i)
     enddo

! in case of flow separation, the last vertex of the modelled jet
! correspond to the tip of the jet whereas the last centroid corresponds
! to the first-last vertex

     if(nsep.eq.1)then
       ygb(ng+1) = ygn(ngo1)
       zgb(ng+1) = zgn(ngo1)
       tgb(ng+1) = tg(ngo1) 
       ycb(ng) = ygn(ngo1-1)
       zcb(ng) = zgn(ngo1-1)
       tcb(ng) = tg(ngo1-1)
     endif 
     write(28,*)
     write(29,*)
     write(29,*)
     write(28,*)

   endif
 endif

return
end 
