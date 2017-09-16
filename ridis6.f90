
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ridis6(k,ng,proat,kget,ygb,zgb,&
                       escr,npc,npt,yn,zn,ycn,zcn,ampli,&
                       ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,&
                       nsep,ksep,phin,ycb,zcb,tcb,jt,tysl,nngo,&
                       phinsl,npsl,di,ang,tc,kord,frint,tn,&
                       nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)

! Description of the variables:
! k =0/1 : recalculate the jet vertices / dont recalaculate
! ng : number of elements in the modelled parte of the jet
! kget 0/1 : jet model not activated /  activated
! ynsl, znsl : coordinates of the  vertices of the free surface
! ygb,zgb : Tobe understood
! escr : growth factor of the discretization 
! npc : number of panels on  the body not considering the ones on the jet
! npt : total number of panels 
! yn,zn : coordinates of the vertices of the body (To be confirmed !)
! ycn,zcn : coordinates of the centroids of the body (To be confirmed !)  
! ampli: is the length of body where the jet model is acting upon (TO BE CONFIRMED) 
! ygn,zgn,ygs2,zgs2,tg : coordinates, 2ond derivatives  and absicissa of the body spline representation of the body 
! ngo: number of points in the spline representation of the body that are not separated (exact when there is no sep)
! ngo1 : ngo + number of separated points (as incremented in splver2) of the body spline representation
! nsep =  0/1 in case of no flow separation/flow separation
! nsepo: dont know what it is , but its not being used 
! ksep =  array with 0 for non separated points / 1 for separated
! phin : potential vector on the body 
! ycnsl,zcnsl : coordinates of the free surface vertices
! ycb,zcb,tcb: centroids  and abscissa of the jet modelled part of the body 
! both are being modified by the current function. 
! jt : interation number of the main loop
! tysl: not being used in the present version 
!ne: fixed at 2 at line 79
!phinsl : potential vector the free surface panels
!npsl: number of panels on the free surface including the ng elements on the free surface , see skecht
! On exit (TO BE CONFIRMED)
!di: length of the vertice 
!ang: angle of the vertex in radians  
!tc: centroid abscissa , TO BE UNDERSTOOD FURTHER  
!frint: fraction of the transition region 
!tn: vertex abscissa , TO BE UNNDERSTOOD FURTHER 
!nnold: last value of nn saved recursively by the this function
!nn1old: last value of nn1 saved recursively by the this function
!ramii: factor used with amii to define the size of the first panel to be meshed by distri a1
!ramiii: "     "  a2
!eskkk: another growth factor of the discretization  
!nngo: running the code its found that its  an integer with the same value as ng; its value is set in a call to prefin.f, 
!in principle it stores the old  value of nng and is used to check if the number of elements in the jet 
!has changed or not (nng is just kget*ng)
!kmed and ksup : are initialized with zero at intial.f and changed to 1 only in ridis5.f (retaled to the separated jet)
!nsep, nsepo: boolean flags related separation 0/1 (nsepo was used in ridis4 but it is not being used here)
!kord: boolean 1/2 related to flow separation not really sure what is its role 
 

!
!      include"slam_p.h"
      implicit real*8 (a-h,o-z)
      parameter (npamx=3000, ntmx = 400000 )
      common/costanti/pi

      dimension ygb(npamx),zgb(npamx),ynsl(npamx),znsl(npamx)
      dimension yn(npamx),zn(npamx),ycn(npamx),zcn(npamx)
      dimension ygn(npamx),zgn(npamx),ygs2(npamx),zgs2(npamx)
      dimension tg(npamx),tgb(npamx),ksep(npamx) 
      dimension phin(npamx)
      dimension tn(npamx),tpo(npamx),phpo(npamx),phpo2(npamx)
      dimension ycnsl(npamx),zcnsl(npamx)
      dimension ycb(npamx),zcb(npamx),tcb(npamx)
      dimension phinsl(npamx),tc(npamx),kord(npamx)
      dimension rft(npamx)
!
      write(*,*) '--> ridis6...........'

! Small hack split FS centroids/panels in another vector 
! old code convention (TO Be corrected once the code is validated)
  nng  = ng*kget
! just parse but now the FS vectors  are defined only locally
  do i=1,npsl  
    ycnsl(i) = ycn(npc+nng+i)
    zcnsl(i) = zcn(npc+nng+i)
    ynsl(i) = yn(npc+nng+i) 
    znsl(i) = zn(npc+nng+i)
  enddo 
! last vertice
  ynsl(npsl+1) = yn(npc+nng+npsl+1) 
  znsl(npsl+1) = zn(npc+nng+npsl+1)
! hack ends

      ne   = 2 
      kk   = k
      eskk = escr
      ninc = 0
      nng  = ng*kget
      epss = 1.d-10
!
      yf = 0.d0
      zf = proat
      if(ng.eq.0)then
        yy   = ynsl(1)
        zz   = znsl(1) 
! in case ng =0 , amii is the amplitude of the first FS paneel
        amii = sqrt((ynsl(2)-ynsl(1))**2+ (znsl(2)-znsl(1))**2 )
        ay   = yy - ygn(iint) 
        az   = zz - zgn(iint)
        aa   = sqrt(ay**2+az**2)
        amdi = tg(iint) + aa
        tt   = amdi 
      write(*,*) '--> ammi amdi 1 .',amii,amdi
      else
        if(kk.eq.0)then
! - con kk=0 ricalcolo vertici getto, con kk=1 no
!
! -  ricalcolo centroidi getto
        write(25,*) ' # jt ng',jt,ng


! On notation:
!ngo  : number of points on the body spline representation  which are not separated 
!ngo1 : ngo + number of separated points  


! Step 1:
! Idea: Find the intersection between a line passing on FS vertice (perpendicular to a body panel) and the 
! The coordinates of this point  are yb,zb  and the corresponding body length is up to the sep. point is ta
! These are saved in ycb(jj), zcb(jj) and tcb(jj)
      do ii=1,ng
          y1  = ycnsl(ii)
          z1  = zcnsl(ii)
          jj  = ng-ii+1 ! runs from ng to 1 (index of the vector ycb,zcb,tcb)
          jjj  = npc+jj ! runs from npc+ng to npc+1 (from the apex of the jet to the bulk on the body side) 
          if(ksep(jjj).ne.1)then
!          else
      do i=1,ngo1-1       ! loop in the body points (separation points are included)
             w1y = zgn(i+1) - zgn(i)
             w1z = -(ygn(i+1) - ygn(i))
             w2y = ygn(i+1) - ygn(i)
             w2z = zgn(i+1) - zgn(i)
             y2  = ygn(i)
             z2  = zgn(i)
             call intret(yi,zi,r1,r2,y1,z1,w1y,w1z,y2,z2,w2y,w2z) ! intersect 2 lines  (y1,z1) + lambda(w1y,w1z) (FS point, normal to the body) and y2,z2 +c* (w2y,w2z) (body point, tangent)
             if(r2.ge.(0.d0-epss).and.r2.le.(1.d0+epss))then
               kint = i ! save the body index where the lines intersect
               ay = yi - ygn(kint) ! calculate the amplitude increment
               az = zi - zgn(kint)
               aa = sqrt(ay**2+az**2)
               ta = tg(kint) + aa ! update body absicissa before the up to the sep point
               call splint(ta,yb,zb,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
               ! given the abscissa find the corresponding centroids and save them in ycb,zcb, tcb (TO UNDERSTAND WHERE THIS VECTORS ARE USED)
               ! In principle , this is definning new coordinates on the body  part where the flow is not separated yet
               ! print this jj, ycb(jj) to check
               ycb(jj) = yb 
               zcb(jj) = zb 
               tcb(jj) = ta
               write(25,'(i4,3d15.7,i4)') jj,yb,zb,ta,0
               goto 999 
             endif
          enddo
! - DA INSERIRE LA DEFINIZIONE DEI PUNTI GETTO NON RIGRIGLIATI (ksep=1)
          else !(ksep==1) 


! in the case the flow is separated just copy the arrays, no need for rearangement 
!TO DO : check the indexation bearing in mind  the new convention we are following        
! So ycb receives body centroid coordinates of the separated points
             ycb(jj) = ycn(jjj) 
             zcb(jj) = zcn(jjj) 
             tcb(jj) = tc(jjj) 
!cc               write(25,'(i4,3d15.7,i4)') jj,ycb(jj),zcb(jj),tcb(jj),1
          endif ! end ksep if
  999     continue
        enddo
         write(25,*)
         write(25,*)
! - calcolo vertici getto
!cc        write(26,*) ' # jt ng',jt,ng
!        write(26,'(i4,3d15.7)') 1,ygb(1),zgb(1),tgb(1)

! Redistribute the  jet vertices by interpolation of the abscissa mid point on the body
! note that tcb was corrected/updated in the precedent loop
      do i=1,ng-1
          tb = 0.5d0*(tcb(i)+tcb(i+1))
          call splint(tb,yb,zb,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          ygb(i+1) = yb 
          zgb(i+1) = zb 
          tgb(i+1) = tb 
!cc          write(26,'(i4,3d15.7)') i+1,yb,zb,tb
      enddo
        if(nsep.eq.0) then
! enforce ygb in case of no separation 
          ygb(ng+1) = ynsl(1) 
          zgb(ng+1) = znsl(1)
        else
! In case of separation the last vertex of the jet is the ...... 
          ygb(ng+1) = ygn(ngo1) 
          zgb(ng+1) = zgn(ngo1)
        endif
! update body abscissa once the last point of the jet (ng+1) has been found 
        dy = ygb(ng+1)-ygb(ng) 
        dz = zgb(ng+1)-zgb(ng)
        dd = sqrt(dy**2+dz**2)
        tgb(ng+1) = tgb(ng)+dd 
!cc        write(26,'(i4,3d15.7)') ng+1,ygb(ng+1),zgb(ng+1),tgb(ng+1)
!cc        write(26,*)
!cc        write(26,*)
!
        endif !ends if(kk==0) 
! - se kk=1 ygb e zgb  (e ampli) e tgb gia' calcolati in shallo
! dont recalculate the jet vertices in this case

        yy = ycb(1)
        zz = zcb(1)
        tt = tcb(1)
        amii = ampli
        amdi = tcb(1)-0.5d0*amii ! abscissa - half of the amplitude (ampli) calculated in shallo.f90
      endif !(end if kk=0/1)
      nn1 = 0

! In the block that follows it seems to be using the amplitudes to
! redistribute the nodes on body and jet part , this looks like another check
! to be clarified: What is amdi, ramii, amii, ramiii and ksup ?
! amii - in case the flow is not sep. is the amplitude of the first FS panel
! amdi - a function of tcb and ampli , calculated in shallo.f90
! ramii is sitill to be understood
 
      if(amdi-20.d0*ramii*amii.gt.tg(ngo).or.ksup.eq.1)then
! this check needs to be understood,  tg -> global body abscissa , what is tg(ngo) ?   

! - non voglio che cominci con kk=0
        if(ksup.eq.0.and.kk.eq.0) goto 98  ! kk == 0 , jet vertices are being recalculated
        ksup = 1
        write(*,*) 'POPOPO' 
        t1 = tg(ngo)
        t0 = 0.5d0*ramiii*amii
        tl = t1-t0
        a1 = ramiii*amii/tl
        a2 =  ramii*amii/tl
        et1 = eskkk
        et2 = eskkk
! This is the first time nn appears, it seems that it is initialized in
! line 343 when kmed and kup ==  0  
        n1 = nn/2
        n2 = nn-n1
! Distribute centroids   
        if(kk.eq.0) then
          kkk = 1
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn,kkk)
          if(nn.ne.nnold)then
            nn  =  nnold
            n2  =  nn - n1
            if(n2.lt.0)then 
              write(*,*) 'ATTENTO A DISTRI, ridis'
            endif 
            call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn,0)
          endif
        else
          kkk = 1
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn,kkk)
        endif 
        tn(1) = 0.d0
        yn(1) = ygn(1) 
        zn(1) = zgn(1) 
        write(55,*) '# jt ' ,jt,nn,kkk 
          write(55,'(i4,3d15.7)') 1,yn(1),zn(1), 0.

! After redistributing the centroids ,  we know nn and now we loop to find the vertex
! Note: The meaning of nn is yet to be clarified
      do i = 1,nn
          ip     = i+1
          dtn    = 0.5d0*(rft(i)+rft(i+1))*tl  
          tn(ip) = t0 + dtn
          call splint(tn(ip),y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
! Interpolate and rewrite the vertices       
          yn(ip) = y 
          zn(ip) = z 
          write(55,'(i4,3d15.7)') ip,y,z,rft(i)
        enddo
        write(55,*)
        r1 = tcb(1)
        r0 = tg(ngo)
        rl = r1-r0
        a1 = ramii*amii/rl
        a2 = 1.d0*amii/rl
        et1 = eskkk
        et2 = eskkk
        n1 = nn/2
        n2 = nn-n1
        if(kk.eq.0) then
          kkk = 1
! redistribute points again, whats the diff between this bit and the one above ?    
! Seems the whole point of the code that follows is to find nn1 rather nn now
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn1,kkk)
          if(nn1.ne.nn1old)then
            nn1  =  nn1old
            n2   =  nn1 - n1
            if(n2.lt.0)then 
              write(*,*) 'ATTENTO A DISTRI 2, ridis'
            endif 
            call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn1,0)
          endif
        else
          kkk = 1
          call distri(a1,a2,rt1,rt2,et1,et2,n1,n2,rft,nn1,kkk)
        endif 
        write(55,*) ' # ciao' ,a1,a2,nn1,kkk
! new vertex interpolation now from 1 to nn1
      do i = 1,nn1
          ip     = nn+1+i ! goes from nn+1 to nn+1+nn1, so this looks to be other region
          dtn    = 0.5d0*(rft(i)+rft(i+1))*rl  
          tn(ip) = r0 + dtn
! find the vertices at the nn+1 up to nn+nn1+1        
          call splint(tn(ip),y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          yn(ip) = y 
          zn(ip) = z 
          write(55,'(i4,3d15.7)') ip,y,z,rft(i)
        enddo
        write(55,*) 
        write(55,*) 
! number of panels to be incremented (there were originally npc panels ) 
        ninc   = nn + nn1 + 1  + ng*(1-kget) - npc 
        nn1old = nn1 
        nnold  = nn 
!
      elseif(amdi-4.d0*amii.gt.tg(ngo).or.kmed.eq.1)then
!  - non voglio che cominci con kk=0
  98    continue  ! if ksup == 0 and kk == 0 we end up here:
        if(kmed.eq.0.and.kk.eq.0) goto 99 ! yet another case to be investigated 
        write(*,*) 'PIPIPI'
        kmed  = 1
        amdi  = tcb(1) - amii
        amdi1 = amdi-tg(ngo)
! Solve the Geom Progression to find the number of nodes to be distributed in this case
! nn1 meshes a region of size tcb(1) - amii-tg(ngo) i.e (amdi1)  with amii panels x growth factor
        nn1   = int( log( 1.d0+(eskk-1.d0)*amdi1/amii)/log(eskk) )
        if(kk.eq.0.and.kmed.eq.1) nn1 = nn1old
        if(nn1.ne.0) then
           amii1 = amdi1*(1.d0-eskk)/(1.d0-eskk**nn1)
          write(*,*) 'PIPIPI nn1 = ',nn1
        else
          write(*,*) 'PIPIPI nn1 = 0'
        endif
        nn1old = nn1
        if(nn1.ne.0)  nn1  = nn1+1
!
        amii0  = amii1*eskk**(nn1-1) 
        amdi0  = tg(ngo)-0.5d0*amii0
! nn meshes a region of size tg(ngo) - 0.5*amii1*eskk**(nn1-1), using eskk as a growth factor
        nn     = int( log( 1.d0+(eskk-1.d0)*amdi0/amii0)/log(eskk) )
        if(kk.eq.0) nn = nnold
        amii0  = amdi0*(1.d0-eskk)/(1.d0-eskk**nn)

! These are number of panels that are going to added
! TO DO : understand the diff between nn and nn1 ? seems one could be on the jet whereas the other could be on the body/ transtion region
        ninc   = nn + nn1 + ng*(1-kget) - npc
        nnold  = nn

!Set the last vertex and work backwars to the body up to the nn+1 vertex 
        dr           = amii 
        r            = amdi + 0.5d0*dr 
        tn(nn+nn1+1) = r
        call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
! last vertex onnce the centriods have been redistributed
        yn(nn+nn1+1) = y 
        zn(nn+nn1+1) = z
        ygb(1)       = y
        zgb(1)       = z
        tgb(1)       = r
        r = r - dr
        dr = amii1 
        write(*,*) 'RID ', amdi1,amii1,nn1
        write(*,*) 'RID ', amdi0,amii0,nn  
        write(28,'(2i4,3d15.7)') jt,nn+nn1+1,r, y,z
      do ip = nn+nn1,nn+1,-1  ! region nn1
          tn(ip)   = r 
          call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          yn(ip)   = y 
          zn(ip)   = z
          ksep(ip) = 1
          if(ip.eq.nn+1) ksep(ip)=2
          r        = r - dr
          dr       = dr*eskk
          write(28,'(2i4,3d15.7)') jt,ip,r,y,z
        enddo
          write(28,'(2i4,3d15.7)') 
        r             = amdi0 
        dr            = amii0 
! now finish up and work it up to the body origin
      do ip = nn,1,-1
          r        = r-dr
          if(ip.eq.1) r = 0.d0
          tn(ip)   = r 
          call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
          yn(ip)   = y 
          zn(ip)   = z
          ksep(ip) = 0
          dr       = dr*eskk
          write(28,'(2i4,3d15.7)') jt,ip,r,y,z
        enddo
          write(28,*) 
          write(28,*) 

      else
 
 99     write(*,*) 'PUPUPU' ! in case kmed.eq.0.and.kk.eq.0 we end up here.

! - !! eskk =! 1.d0  !!!

! definition of nn , **important** , number of elements to m mesh amdi using amplitude amii (TO BE CHECKED)
      nn = int( log( 1.d0+(eskk-1.d0)*amdi/amii)/log(eskk) )
      no = npc - ng*(1-kget)
! - calcolo differenza
! panels to be added/cut 
      ninc = nn+ng*(1-kget)-npc
! - non cambio se kk=0
      if(kk.eq.0.and.ng.gt.0) then
        ninc = 0
        nn   = no
      endif
! before the jet modelling kicks in 
! no panel is ever cut 
      if(ninc.lt.0.and.ng.eq.0)then
        nn   = no
        ninc = 0
      elseif(ninc.lt.0)then
        write(*,*) 'taglio   ',ninc,' pannelli dal corpo'
      elseif(ninc.gt.0)then
        write(*,*) 'aggiungo ',ninc,' pannelli sul corpo'
      else
        write(*,*) 'npc costante, ninc= ',ninc,ng,kget
        write(*,*) yy,zz
      endif
      nnold = nn
! - ricalcolo amii
! amii is a factor related to the panel amplitude  and 
! is probably (TO check) being calculated
! recursively by this subroutine
      amii = amdi*(1.d0-eskk)/(1.d0-eskk**nn)
      write(*,*) 'amii amdi2 ',amii,amdi,nn

! - griglio bulk (in this case there is no jet portion to be meshed)
      r        = amdi 
      tn(nn+1) = r
      call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
      yn(nn+1) = y 
      zn(nn+1) = z
      ygb(1)   = y
      zgb(1)   = z
      tgb(1)   = r
      dr       = amii
!cc      write(27,'(2i4,3d15.7)') jt,nn+1,r, y,z
! nn indeed seems to be related to number of panels in the bulk region of the body 
! so that nn1 should relate the number of points added or cut in the jet modelled part (To be checked!)
      do ip = nn,1,-1
        r        = r-dr
        if(ip.eq.1) r = 0.d0
        tn(ip)   = r 
        call splint(r,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        yn(ip)   = y 
        zn(ip)   = z
        dr       = dr*eskk
!cc        write(27,'(2i4,3d15.7)') jt,ip,r,y,z
      enddo
!
      endif 

! At this point it seems the body is remeshed 
!Therefore npc can be updated 

!cc        write(27,*)
      npco = npc
      npc  = nn  + nn1 +  ng*(1-kget)
      npt  = npt + ninc  
! - aggiungo getto a yn zn, se ng>0 e kget = 0
      if(ng.gt.0)then
      do i=2,ng+1
          yn(nn+nn1+i) = ygb(i)  
          zn(nn+nn1+i) = zgb(i)
          tn(nn+nn1+i) = tgb(i)
!cc        write(27,'(2i4,3d15.7)') jt,npc+i,tgb(i), ygb(i),zgb(i)
        enddo
      endif 



!cc      write(27,*) 
!cc      write(27,*)
!
      if(kk.eq.0)then
!
! - aggiorno solo centroidi del corpo, se kk=0
      tt = 0.d0
      do i = 1,npc+nng
       if(ksep(i).ne.1)then
        tt = 0.5d0*(tn(i)+tn(i+1)) 
        call splint(tt,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        ycn(i)  = y 
        zcn(i)  = z
!cc        write(25,'(2i4,2d15.7)') jt,i,y,z 
       endif
      enddo
!cc      write(25,*) 
!cc      write(25,*)
!
      else
!
! --  ascissa curvilinea
      tpo(npco+nngo+1) = tn(npc+nng+1)
!      npo    = npco+nngo+1
! First time nngo appears here, it is comming as an input to the subroutine as well
! To Do : check what it means/is ?  
      npo    = npco+nngo
      dy     = ycn(1)-yn(1)
      dz     = zcn(1)-zn(1)
      tt     = sqrt(dy**2+dz**2) 
      tpo(1) = tt
      phpo(1)= phin(1) 
      write(*,*) 'ppppp0 ',npco,nngo  
!cc      write(24,*) '# jt',jt,nsep

! Build a potential representation as function 
! of curviliniar abiscissa , i.e a table : tpo x phpo
      do i=2,npco+nngo
        dy = ycn(i)-ycn(i-1)
        dz = zcn(i)-zcn(i-1)
        dd     = sqrt(dy**2+dz**2)
        tt     = tt+dd
        tpo(i) = tt
        phpo(i)= phin(i)
!cc        write(24,'(i4,4d15.7)') i,ycn(i),zcn(i),tpo(i),phpo(i)
      enddo
!cc      write(24,*)
!cc      write(24,*)
! -- estrapolo potenziale:
! Very close expansion to what we have done at ridis.f90 really
      write(*,*) 'ppppp0 '  
      dtt  = tpo(npco+nngo)-tpo(npco+nngo-1)
      write(*,*) dtt,tpo(npco+nngo+1),tpo(npco+nngo)
      phia = phin(npco+nngo)
      phib = phin(npco+nngo-1)
      write(*,*) 'ppppp1 '  
      write(*,*) phia,phib
      delphi = (phia-phib)/dtt 
      write(*,*) 'ppppp2 ' ,ng 
      write(*,*) phia,delphi,tpo(npco+nngo+1),tpo(npco+nngo)
      phie = phia+delphi*(tpo(npco+nngo+1)-tpo(npco+nngo))
      phpo(npco+nngo+1) = phie
      write(*,*) npco+nngo+1,tpo(npco+nngo+1),phie
!c INTERPOLAZ SPLINE POTENZIALKE>>
      yp1 = 1.d31
      ypn = 1.d31
      call spline1(phpo,phpo2,tpo,yp1,ypn,npo,npamx)
! -- aggiorno divisione corpo SL, e tutti i centroidi ,se kk=1
      tt  = 0.d0
!      tysl = 1d31
!      if(ng.eq.0)then
!c      ycn(npc+nng)  = ycn(npco+nngo)
!c      zcn(npc+nng)  = zcn(npco+nngo)
!c      phin(npc+nng) = phin(npco+nngo)
!      else
!      ycn(npc+nng)  = ycb(nng)
!      zcn(npc+nng)  = zcb(nng)
!      endif
!      tt = tcb(nng)
!      call splint1(tt,ph,phpo,phpo2,tpo,npo,npamx)
!      phin(npc+nng) = ph
!cc      write(23,*) '# jt',jt,nsep

!This parte is related to jet and introduces the same nomenclature as in solv22 
! m  -> transition region
! n   -> fully modelled part
 
      m  = int(frint*ng)
      n  = ng - m
! how do we know that npc+nng is the  separation point really ? 
! what is nng by the way  ? (To do!)
! mysterious kord array is finally found ...      
      if(nsep.eq.1) ksep(npc+nng) = 1

      tto = 0.d0
! body loop  up to the separation point 
      do i = 1,npc+nng
        if(i.le.npc)then
        tt = 0.5d0*(tn(i)+tn(i+1))
        else
        tt = tcb(i-npc)
        endif 
        call splint(tt,y,z,ygn,zgn,ygs2,zgs2,tg,ngo1,npamx,0)
        if(nsep.eq.1)then
        epst = (tt-tto)/3.d0
        if(tt.ge.tg(ngo)-epst)then
          if(ksep(i-1).eq.0)then
          !special treatment for the for the first separated point         
            ksep(i) = 2
            kord(i) = 2
            ichin   = i ! this variable is new , check from where does it come from !
!          elseif(y.lt.tysl)then
!            ksep(i) = 2
          elseif(i-npc.le.m+0.25d0*n)then
          ! for separated points in the region m+0.25n (trasition plus the first quarter of the jet modelled part)
            ksep(i) = 1
            kord(i) = 2
          else
          ! other separation points 
            ksep(i) = 1
            kord(i) = 1
          endif 
        else
!no separation case 
          ksep(i) = 0
          kord(i) = 2
        endif
        endif
!        endif
! Assign centroids returned by splint
        ycn(i)  = y 
        zcn(i)  = z
        tc(i)   = tt
        tto     = tt
        if(i.gt.1)then
! Interpolate the potential        
        call splint1(tt,ph,phpo,phpo2,tpo,npo,npamx)
        phin(i) = ph 
        endif
!cc        write(23,'(i4,5d15.7,i4)') i,tt,tg(ngo),y,z,ph,ksep(i)
      enddo
!cc      write(23,'(i4,5d15.7,i4)') npc+nng,tpo(npo),tg(ngo),ycn(npc+nng),
!cc     #      zcn(npc+nng),phin(npc+nng),ksep(npc+nng)
!cc      write(23,*)
!cc      write(23,*)


! This is Important as it tackles flow separation
! - controllo sulla separazione
      if(nsep.eq.0)then
! My feeling is that ygn(ngo) is coordinate of the body knuckle 
! if that is the case the separation check boils down to an if statment on the y 
! direction, i.e if the flow passed the knuckle then it is classified as separated 
        if(ycn(npc+nng-ne-1).gt.ygn(ngo))then
           write(*,*) 'RIDIS SEPARA .. ',nsep,npc+nng
           write(*,*) ne,ygn(ngo) ! TO DO : what  ngo and  nng  ?  ne  = 2 (line 30) and does not change throughtout the code
           write(*,*) ycn(npc+nng-ne-1),zcn(npc+nng-ne-1) 
           write(*,*) ycn(npc+nng-ne),zcn(npc+nng-ne) 
           write(*,*) ycn(npc+nng-ne+1),zcn(npc+nng-ne+1) 
           write(*,*) ycn(npc+nng),zcn(npc+nng) 
           nsep = 1
           write(*,*) 'The code does not treat separation yet!'
           stop 
          
!           ng   = ng-ne
!           nng  = nng-ne
!           yn(npc+nng+1) = yn(npc+nng+1+ne)
!           zn(npc+nng+1) = zn(npc+nng+1+ne)
!           ksep(npc+ng)   = 2
!           ksep(npc+ng-1) = 2
!           npsl = npsl - ne
!           do i = 1,npsl
!             ycnsl(i)  = ycnsl(i+ne)
!             zcnsl(i)  = zcnsl(i+ne)
!             ynsl(i+1) = ynsl(i+1+ne)
!             znsl(i+1) = znsl(i+1+ne)
!             phinsl(i) = phinsl(i+ne)
!           enddo

! Classify the points in the vicinity of and the sepration point
      do i=npc+ng-ne,npc+ng
              ksep(i) = 1
              kord(i) = 1
            enddo
! not sure what point is it, but looks special , need to figure out what ne stands for really.
            ksep(npc+ng-ne-1) = 2
            kord(npc+ng-ne-1) = 2
!            kord(npc+ng-ne-1) = 1
! - lunghezza e angolo del vertice 
           ay = ynsl(1)-ycn(npc+ng)  ! in the old indexation npc+ng is the tip of the jet (pointed tip)    
           az = znsl(1)-zcn(npc+ng)  ! vertex - center
           aa = sqrt(ay**2+az**2)      
           by = ycnsl(1)-ycn(npc+ng) ! center - center     
           bz = zcnsl(1)-zcn(npc+ng) 
           bb = sqrt(by**2+bz**2)      
           sc = (ay*by + az*bz)/(aa*bb)   ! inverse of the cosine (dot product definition)
           th = acos(sc)          ! th is angle in radians
           di = aa                 
           ang= th
           write(*,*) 'e OORA ' ,di,ang*180./pi
        endif
      endif
!
      endif     
!
      return
      end
