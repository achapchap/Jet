      program main
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ++                                                 ++
!     +  programma di calcolo del flusso generato dalla   +
!     +    caduta di un wedge avente un deadrise angle    +
!     +               pari ad alfa (gradi)                +
!     ++                                                 ++
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

!     - File inclusi: 

!       -- "slam_p.h" in cui sono riportati la definizione 
!                     delle variabili e i parametri per il 
!                     dimensionamento,

      include"slam_p.h"
      include"slam_f.h"
 
!     Key assumptions about variables and declaration are as follows:
!     1 - implicit real*8 (a-h,o-z) is used in slam_p.hi, so be aware       
!     2 - any variables beginning with the letters I, J, K, L, M, N are 
!     assumed to be INTEGER by the compiler     

      dimension yv(npamx),zv(npamx),yn(npamx),zn(npamx)
      dimension ysl(npamx),zsl(npamx),ynsl(npamx),znsl(npamx)
      dimension phi(npamx),phisl(npamx),phinsl(npamx)
      dimension ygb(npamx),zgb(npamx)
      dimension dddphit(-npamx:npamx),ddphit(-npamx:npamx)
      dimension xi(-npamx:npamx),xin(-npamx:npamx)
      dimension h(-npamx:npamx),dh(-npamx:npamx),hh(-npamx:npamx)
      dimension dphi(npamx),dphisl(npamx)
      dimension dpht(npamx),dphtsl(npamx)
      dimension yce(npamx),zce(npamx),ycn(npamx),zcn(npamx)
      dimension ycsl(npamx),zcsl(npamx),ycnsl(npamx),zcnsl(npamx)
      dimension tmy(npamx),tmz(npamx),rny(npamx),rnz(npamx)
      dimension tmysl(npamx),tmzsl(npamx),rnysl(npamx),rnzsl(npamx)
      dimension amp(npamx),ampsl(npamx)
      dimension vym1(npamx),vzm1(npamx),vym2(npamx),vzm2(npamx)
      dimension vymsl1(npamx),vzmsl1(npamx),vymsl2(npamx),vzmsl2(npamx)
      dimension depn1(npamx),depn2(npamx),zcr(npamx)
      dimension pre(npamx)
      dimension dpt(npamx),dpnt(npamx),dptsl(npamx),dpntsl(npamx)
      dimension dpht2(npamx)
      dimension pres(npamx),yco(npamx),zco(npamx)
      dimension vyo(npamx),vzo(npamx),phio(npamx) 
      dimension dphtbsl(npamx),phb(npamx),ph(npamx),dphn(npamx)
      dimension dphnb(npamx) 
      dimension phib(npamx),dphibsl(npamx),rl(npamx) 
      dimension a1(npamx),b1(npamx),c1(npamx),d1(npamx),e1(npamx)
      dimension xig(-npamx:npamx),xigs(-npamx:npamx)
      dimension zeb(-npamx:npamx),zegs(-npamx:npamx),zef(-npamx:npamx)
      dimension hp(npamx),xj(npamx),ze(npamx),xis(npamx),zes(npamx)
      dimension dpttsl(npamx),dpb(npamx),dptb(npamx),dpntbsl(npamx)
      dimension a2(npamx),b2(npamx),c2(npamx),d2(npamx),e2(npamx)
      dimension pre2(npamx),dpt2(npamx)
      dimension vxi(npamx),dpntt(npamx)
      dimension xigb(-npamx:npamx),zegb(-npamx:npamx)
      dimension xigf(-npamx:npamx),zegf(-npamx:npamx)
      dimension ry(npamx),rz(npamx),ty(npamx),tz(npamx)
      dimension ygn(npamx),zgn(npamx),ygs2(npamx),zgs2(npamx)
      dimension tg(npamx),tgb(npamx),ksep(npamx),kse(npamx) 
      dimension phin(npamx),phg(npamx),phgs2(npamx)
      dimension depn1s(npamx),depn2s(npamx)
      dimension tcb(npamx),ycb(npamx),zcb(npamx),tc(npamx)
      dimension kord(npamx),kor(npamx),tn(npamx)
!
!
      call input(vfall0,ancut,iiget,jjget,frint,rmg,epsgg,eskg,&
     &       gfrac,pro0 ,frdt,ampp,pfraz,escr,estr,tend,ksta,scon,svel,&
     &       spot,spre,idis,ift,iford,frtend,ramii,ramiii,eskkk)
      call initial(t,amplim,vfall,llf,spo0,pro0,npc,npt,npsl,&
     &           kget,ng,yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,amp,ampsl,&
     &           rny,rnz,rnysl,rnzsl,tmy,tmz,tmysl,tmzsl,ampp,vfall0,&
     &           dphi,phisl,pfraz,estr,escr,epsg,epsgg,ampli,&
     &           ygn,zgn,ygs2,zgs2,tg,ngo,ngo1,nsep,nsepo,ksep,kord,&
     &         nnold,nn1old,frin,frfi,jind,tin,jfid,tfi,jend,kmed,ksup)

         write(*,*) 'VFALL0, VFALL ',vfall0,vfall
!
! - soluzione, calcolo velocita
      call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,&
     &                  phi,dphi,phisl,dphisl,jt)
!
      call calvel(npsl,npc,ng,kget,mb,mf,mt,m,n,nt,ntt,&
     &            phi,phib,phb,dpht,dphi,dphisl,dphnb,dphtsl,phisl,&
     &            rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,&
     &            amp,vym1,vzm1,rny,rnz,tmy,tmz,&
     &            ampsl,vymsl1,vzmsl1,rnysl,rnzsl,tmysl,tmzsl,vxi,&
     &            ry,rz,ty,tz,kse)


     
!
! --------- CALCOLO LUNGHEZZA STRISCIA DI CONTROLLO
          rrl=0.d0
      do i=jind+1,jfid-1
            dl= sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
            rrl= rrl + dl
          enddo
          dlin= (1.d0-tin)*sqrt( (ycsl(jind+1)-ycsl(jind))**2 +&
     &                    (zcsl(jind+1)-zcsl(jind))**2 )
          dlfi= tfi*sqrt( (ycsl(jfid+1)-ycsl(jfid))**2 +&
     &                    (zcsl(jfid+1)-zcsl(jfid))**2 )
          rrl= rrl + dlin + dlfi
          yin = ycsl(jind)+tin*(ycsl(jind+1)-ycsl(jind))
          zin = zcsl(jind)+tin*(zcsl(jind+1)-zcsl(jind))
          yfi = ycsl(jfid)+tfi*(ycsl(jfid+1)-ycsl(jfid))
          zfi = zcsl(jfid)+tfi*(zcsl(jfid+1)-zcsl(jfid))
          write(77,'(i4,2d15.6,2i4,4d15.6)') 0 ,0.d0,rrl,jind,jfid+1,&
     &               yin,zin,yfi,zfi                             
           
!-----------
!
      call stampa(vfall,ng,kget,npc,npsl,jt,t,dt,frdt,llf,&
     &                  scon,svel,spot,spre,phi,dphi,phisl,dphisl, &
     &                 yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,&
     &                 vym1,vzm1,vymsl1,vzmsl1,&
     &                 dpht,dphtsl,dpt,dpnt,pre,pre2,pres,vxi,&
     &                 jend,jind,jfid,yin,zin,yfi,zfi)



!
!----------------------------------------------------------------------
!- inizio integrazione temporale
      proat = zv(1)
          write(*,*) ' PROAT prima ',proat
      jt    = 0
      do while (t.lt.tend)
        jt=jt+1
        vfall = vfall0
        if(jt.lt.30) vfall=vfall0*float(jt)/float(30) !E LA PRIMA ITERAZIONE ?
!- calcolo passo di integrazione, e volume
        write(*,*) 'prima check, nsep',nsep
        
        write(*,*) 'After VFALL0, VFALL ',vfall0,vfall
        call check(yv,zv,npc,npsl,dt,t,jt,ampsl,frdt,vymsl1,vzmsl1,&
     &             ycsl,zcsl,kget,frtend,tend,zgn,ngo,vfall,ngo1,nsep)
        write(*,*) 'dt jt ',dt,jt,vfall

!- spostamento (predictor) dei centroidi 
        write(*,*) 'predictor',npsl
!cc        write(15,*) '#jt ',jt,ng
!cc        write(16,*) '#jt ',jt,ng
!      if(jt.eq.254)then
!        ksep(npc+ng)   =1
!        ksep(npc+ng-1) =1
!        ksep(npc+ng-2) =1
!        ksep(npc+ng-3) =1
!        ksep(npc+ng-4) =1
!      endif

          write(*,*) 'After2 VFALL0, VFALL ',vfall0,vfall
      do ip = 1,npc+ng
!cc          write(15,'(2d15.7,i4)') vym1(ip),vzm1(ip),ksep(ip)
          if(nsep.eq.1)then
          if(ksep(ip).ne.0)then
          deph2      = vym1(ip)**2+ vzm1(ip)**2
          depn1s(ip) = deph2/2.d0
          ycn(ip)    = yce(ip) + vym1(ip)*dt/2.d0
          zcn(ip)    = zce(ip) + vzm1(ip)*dt/2.d0
          phin(ip)   = phi(ip) + depn1s(ip)*dt/2.d0
          else 
          phin(ip)   = phi(ip) 
          endif 
!cc          write(16,'(2d15.7,i4)') phi(ip),phin(ip),ksep(ip)
          if(ksep(ip).eq.2)then
            ay  = ycn(ip)  - ygn(ngo)
            ay1 = ycn(ip-1)- ygn(ngo)
            if(ay1.gt.0.d0)then
              ksep(ip) = 1
            else
              ycn(ip)  = yce(ip)
              zcn(ip)  = zce(ip)
              phin(ip) = phi(ip) 
            endif 
          endif 
          else
              ycn(ip)  = yce(ip)
              zcn(ip)  = zce(ip)
              phin(ip) = phi(ip)
          endif 
        enddo
!cc        write(15,*) 
!cc        write(16,*) 
      do ip = 1,npsl
!cc          write(15,'(2d15.7,i4)') vymsl1(ip),vzmsl1(ip),ksep(ip)
          deph2        = vymsl1(ip)**2 + vzmsl1(ip)**2
          depn1(ip)    = deph2/2.d0
          ycnsl(ip)    = ycsl(ip) + vymsl1(ip)*dt
          zcnsl(ip)    = zcsl(ip) + vzmsl1(ip)*dt
          phinsl(ip)   = phisl(ip)+ depn1(ip)*dt
!cc          write(16,'(2d15.7)') phisl(ip),phinsl(ip)
        enddo
!cc        write(15,*) 
!cc        write(15,*) 
!cc        write(16,*) 
!cc        write(16,*) 
!- sposto vertice corpo
        proat = zv(1) - vfall*dt
          write(*,*) ' PROAT, VFALL, DT ',proat,vfall,dt
!
!cckkk        write(30,*) '# jt ',jt,ng
!cckkk        do i=1,npc+ng
!cckkk          write(30,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(30,*)
!cckkk        write(30,*)
!cckkk        write(31,*) '# jt ',jt
!cckkk        do i=1,npc+ng+1
!cckkk          write(31,*) yn(i),zn(i),ksep(i)
!cckkk        enddo
!cckkk        write(31,*)
!cckkk        write(31,*)
!cckkk        write(32,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(32,'(4d15.7)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(32,*)
!cckkk        write(32,*)
!cckkk        write(33,*) '# jt ',jt
!cckkk        do i=1,npsl+1   
!cckkk          write(33,*) ynsl(i),znsl(i)
!cckkk        enddo
!cckkk        write(33,*)
!cckkk        write(33,*)
!cckkk        write(34,*) '# jt ',jt
!cckkk        do i=1,ngo1     
!cckkk          write(34,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
!cckkk        enddo
!cckkk        write(34,*)
!cckkk        write(34,*)
!- ricostruisco la configurazione di tentativo dei vertici
        call splver2(ycnsl,zcnsl,ynsl,znsl,npsl,proat,&
     &               kget,estr,ng,ygn,zgn,ygs2,zgs2,tg,ngo,iint,jt,&
     &               ampli,ngo1,nsep,nsepo,ycn,zcn,npc,phg,phgs2,&
     &               phin,ksep,di,ang,tc)

!
!cckkk        write(40,*) '# jt ',jt,ng
!cckkk        do i=1,npc+ng
!cckkk          write(40,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(40,*)
!cckkk        write(40,*)
!cckkk        write(41,*) '# jt ',jt
!cckkk        do i=1,npc+ng+1
!cckkk          write(41,*) yn(i),zn(i),ksep(i)
!cckkk        enddo
!cckkk        write(41,*)
!cckkk        write(41,*)
!cckkk        write(42,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(42,'(4d15.7)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(42,*)
!cckkk        write(42,*)
!cckkk        write(43,*) '# jt ',jt,ng
!cckkk        do i=1,npsl+1   
!cckkk          write(43,*) ynsl(i),znsl(i)
!cckkk        enddo
!cckkk        write(43,*)
!cckkk        write(43,*)
!cckkk        write(44,*) '# jt ',jt
!cckkk        do i=1,ngo1     
!cckkk          write(44,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
!cckkk        enddo
!cckkk        write(44,*)
!cckkk        write(44,*)
!- rigriglio corpo
           write(*,*) ' AAAA - NG input ',ng
!        call ridis5(0,ng,proat,kget,ynsl,znsl,ygb,zgb,
!     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
!     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
!     #                  nsep,nsepo,ksep,phg,phgs2,phin,
!     #                  ycnsl,zcnsl,ycb,zcb,tcb,jt,tysl,nngo,
!     #                  ne,phinsl,npsl,di,ang,tc,kord,frint,tn,
!     #                  nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)


      call ridis6(0,ng,proat,kget,ygb,zgb,ynsl,znsl,ycnsl,zcnsl,phinsl,&
                       escr,npc,npt,yn,zn,ycn,zcn,ampli,&
                       ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,&
                       nsep,ksep,phin,ycb,zcb,tcb,jt,nngo,&
                       npsl,di,ang,tc,kord,frint,tn,&
                       nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)




!- normali tangenti  ampiezze e boundary conditions
        call nortan(vfall,yn,zn,amp,tmy,tmz,rny,rnz,npc,&
     &            ynsl,znsl,ampsl,tmysl,tmzsl,rnysl,rnzsl,npsl,&
     &            dphi,kget,ng,ycn,zcn,yce,zce,yv,zv,&
     &            ycnsl,zcnsl,ycsl,zcsl,ysl,zsl,phinsl,phisl,&
     &            phin,phi,ygb,zgb,0,ksep,jt)
!- trattamneto getto
        if(kget.eq.1)then
            call  get(ng,npc,ycnsl,zcnsl,ycn,zcn,yn,zn,  &
     &               xigs,zegs,xigb,zegb,xigf,zegf) 
!
        endif
!
!cckkk        write(50,*) '# jt ',jt
!cckkk        do i=1,npc+ng
!cckkk          write(50,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(50,*)
!cckkk        write(50,*)
!cckkk        write(51,*) '# jt ',jt
!cckkk        do i=1,npc+ng+1
!cckkk          write(51,*) yn(i),zn(i),ksep(i)
!cckkk        enddo
!cckkk        write(51,*)
!cckkk        write(51,*)
!cckkk        write(52,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(52,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(52,*)
!cckkk        write(52,*)
!cckkk        write(53,*) '# jt ',jt
!cckkk        do i=1,npsl+1   
!cckkk          write(53,*) ynsl(i),znsl(i)
!cckkk        enddo
!cckkk        write(53,*)
!cckkk        write(53,*)
!cckkk        write(54,*) '# jt ',jt
!cckkk        do i=1,ngo1     
!cckkk          write(54,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
!cckkk        enddo
!cckkk        write(54,*)
!cckkk        write(54,*)
        
!- risolvo il problema nella configurazione intermedia
        print*,'t_old ',t
        t   = t + dt
        nng = kget*ng
        print*,'t_new ',t
        if(kget.eq.0)then
          call solv1(npc,nng,npsl,yn,zn,ygb,zgb,ynsl,znsl,amp,ampsl,&
     &                  phin,dphi,phinsl,dphisl,jt)
        elseif(kget.eq.1)then
!          call solv2(frint,ng,npc,npsl,ycn,zcn,yn,zn,ynsl,znsl,
!     #            ycnsl,zcnsl,dphi,phinsl,dphtsl,dphtbsl,phb,
!     #           xigs,zegs,xigb,zegb,xigf,zegf,
!     #           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
!     #          ph,dphn,dphnb,a1,b1,c1,d1,e1,phi,phib,dphisl,dphibsl,rl,
!     #          mb,mf,mt,m,n,nt,ntt,xj,ze,xis,zes,jt)
          call solv22(frint,ng,npc,npsl,ycn,zcn,yn,zn,ynsl,znsl,&
     &            ycnsl,zcnsl,dphi,phinsl,dphtsl,dphtbsl,phb,&
     &           xigs,zegs,xigb,zegb,xigf,zegf,&
     &           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,&
     &        ph,dphn,dphnb,a1,b1,c1,d1,e1,phin,phib,dphisl,dphibsl,rl,&
     &          mb,mf,mt,m,n,nt,ntt,xj,ze,xis,zes,jt,ksep,kse,kord,kor)
          call calsol(mb,mf,mt,m,n,nt,ntt,phin,phib,phb,dphisl,&
     &           dphnb,rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,ry,rz,kse,dphi)
        endif
!- calcolo velocita
          call calvel(npsl,npc,ng,kget,mb,mf,mt,m,n,nt,ntt,&
     &            phin,phib,phb,dpht,dphi,dphisl,dphnb,dphtsl,phinsl,&
     &              rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,&
     &              amp,vym2,vzm2,rny,rnz,tmy,tmz,&
     &              ampsl,vymsl2,vzmsl2,rnysl,rnzsl,tmysl,tmzsl,vxi,&
     &              ry,rz,ty,tz,kse)
!- spostamento definitivo dei centroidi
        write(*,*) 'corrector '
        tysl = 1.d31
      do ip = 1,npc+ng
          if(ksep(ip).eq.1)then
          deph2      = vym2(ip)**2 + vzm2(ip)**2
          depn2s(ip) = deph2/2.d0
          ycn(ip)    = yce(ip)+(vym1(ip)+vym2(ip))*dt/2.d0
          zcn(ip)    = zce(ip)+(vzm1(ip)+vzm2(ip))*dt/2.d0
          phin(ip)   = phi(ip)+(depn1s(ip)+depn2s(ip))*dt/2.d0
          tysl       = min(tysl,ycn(ip))
          endif 
        enddo
      do ip = 1,npsl
          deph2       = vymsl2(ip)**2 + vzmsl2(ip)**2
          depn2(ip)   = deph2/2.d0
          ycnsl(ip)   = ycsl(ip)+(vymsl1(ip)+vymsl2(ip))*dt/2.d0
          zcnsl(ip)   = zcsl(ip)+(vzmsl1(ip)+vzmsl2(ip))*dt/2.d0
          phinsl(ip)  = phisl(ip)+(depn1(ip)+depn2(ip))*dt/2.d0
        enddo



!- ricostruisco la configurazione della superficie libera
        call splver2(ycnsl,zcnsl,ynsl,znsl,npsl,proat,&
     &               kget,estr,ng,ygn,zgn,ygs2,zgs2,tg,ngo,iint,jt,&
     &               ampli,ngo1,nsep,nsepo,ycn,zcn,npc,phg,phgs2,&
     &               phin,ksep,di,ang,tc)
!- calcolo angolo intersezion, event. separo (e/o rigriglio) il getto
!cckkk        write(60,*) '# jt ',jt,ng,ngo1
!cckkk        do i=1,npc+ng
!cckkk          write(60,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(60,*)
!cckkk        write(60,*)
!cckkk        write(61,*) '# jt ',jt
!cckkk        do i=1,npc+ng+1
!cckkk          write(61,*) yn(i),zn(i),ksep(i)
!cckkk        enddo
!cckkk        write(61,*)
!cckkk        write(61,*)
!cckkk        write(62,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(62,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(62,*)
!cckkk        write(62,*)
!cckkk        write(63,*) '# jt ',jt
!cckkk        do i=1,npsl+1   
!cckkk          write(63,*) ynsl(i),znsl(i)
!cckkk        enddo
!cckkk        write(63,*)
!cckkk        write(63,*)
!cckkk        write(64,*) '# jt ',jt
!cckkk        do i=1,ngo1     
!cckkk          write(64,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
!cckkk        enddo
!cckkk        write(64,*)
!cckkk        write(64,*)

        write(*,*) 'jind jifid 0',jind,jfid
        call shallo(ynsl,znsl,ygb,zgb,kget,jt,npsl,&
     &               ng,ampli,rmg,iiget,jjget,&
     &               ycnsl,zcnsl,phinsl,dphisl,epsg,epsgg,eskg,gfrac,&
     &               ygn,zgn,ygs2,zgs2,tg,tgb,ngo,ngo1,nsep,&
     &               ycb,zcb,tcb,nngo,jind,jfid,jend)
       


!cckkk        write(70,*) '# jt ',jt,ng
!cckkk        do i=1,npc+ng
!cckkk          write(70,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(70,*)
!cckkk        write(70,*)
!cckkk        write(71,*) '# jt ',jt
!cckkk        do i=1,npc+ng+1
!cckkk          write(71,*) yn(i),zn(i),ksep(i)
!cckkk        enddo
!cckkk        write(71,*)
!cckkk        write(71,*)
!cckkk        write(72,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(72,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(72,*)
!cckkk        write(72,*)
!cckkk        write(73,*) '# jt ',jt
!cckkk        do i=1,npsl+1   
!cckkk          write(73,*) ynsl(i),znsl(i)
!cckkk        enddo
!cckkk        write(73,*)
!cckkk        write(73,*)
!cckkk        write(74,*) '# jt ',jt
!cckkk        do i=1,ngo1     
!cckkk          write(74,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
!cckkk        enddo
!cckkk        write(74,*)
!cckkk        write(74,*)
!- redistribuzione centroidi e vertici SL
        if (mmm.eq.0.or.jt.le.jjget) then
        call disun2(jt,yn,zn,ng,ynsl,znsl,ycnsl,zcnsl,phinsl,npsl,&
     &    ygb,zgb,escr,kget,estr,amplim,ampli,npt,npc,iiget,&
     &    ycb,zcb,jind,jfid,tin,tfi,jend)
!- ridefinizione (e ridiscretizzazione) dei pannelli sul corpo
        write(*,*) 'jind jifid 2',jind,jfid
        nng=ng*kget
!cckkk        write(80,*) '# jt ',jt,ng
!cckkk        do i=1,npc+nngo
!cckkk          write(80,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(80,*)
!cckkk        write(80,*)
!cckkk        write(81,*) '# jt ',jt
!cckkk        do i=1,npc+nngo+1
!cckkk          write(81,*) yn(i),zn(i),ksep(i)
!cckkk        enddo
!cckkk        write(81,*)
!cckkk        write(81,*)
!cckkk        write(82,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(82,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(82,*)
!cckkk        write(82,*)
!cckkk        write(83,*) '# jt ',jt
!cckkk        do i=1,npsl+1   
!cckkk          write(83,*) ynsl(i),znsl(i)
!cckkk        enddo
!cckkk        write(83,*)
!cckkk        write(83,*)
!cckkk        write(84,*) '# jt ',jt
!cckkk        do i=1,ngo1     
!cckkk          write(84,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
!cckkk        enddo
!cckkk        write(84,*)
!cckkk        write(84,*)
!        if(jt.eq.iiget+1) stop
        nng=ng*kget
        mmm=mod(jt,idis)
           write(*,*) ' BBBB - NG input ',ng
!        call ridis5(1,ng,proat,kget,ynsl,znsl,ygb,zgb,
!     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
!     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
!     #                  nsep,nsepo,ksep,phg,phgs2,phin,
!     #                  ycnsl,zcnsl,ycb,zcb,tcb,jt,tysl,nngo,
!     #                  ne,phinsl,npsl,di,ang,tc,kord,frint,tn,
!     #                  nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)

      call ridis6(1,ng,proat,kget,ygb,zgb,ynsl,znsl,ycnsl,zcnsl,phinsl,&
                       escr,npc,npt,yn,zn,ycn,zcn,ampli,&
                       ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,&
                       nsep,ksep,phin,ycb,zcb,tcb,jt,nngo,&
                       npsl,di,ang,tc,kord,frint,tn,&
                       nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)

         





        else
!- ridefinizione (e ridiscretizzazione) dei pannelli sul corpo
           write(*,*) ' CCCC - NG input ',ng
!        call ridis5(0,ng,proat,kget,ynsl,znsl,ygb,zgb,
!     #                  escr,npc,npt,yn,zn,ycn,zcn,ampli,
!     #                  ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,
!     #                  nsep,nsepo,ksep,phg,phgs2,phin,
!     #                  ycnsl,zcnsl,ycb,zcb,tcb,jt,tysl,nngo,
!     #                  ne,phinsl,npsl,di,ang,tc,kord,frint,tn,
!     #                  nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup)
 

      call ridis6(0,ng,proat,kget,ygb,zgb,ynsl,znsl,ycnsl,zcnsl,phinsl,&
                       escr,npc,npt,yn,zn,ycn,zcn,ampli,&
                       ygn,zgn,ygs2,zgs2,tg,ngo,tgb,iint,ngo1,&
                       nsep,ksep,phin,ycb,zcb,tcb,jt,nngo,&
                       npsl,di,ang,tc,kord,frint,tn,&
                       nnold,nn1old,ramii,ramiii,eskkk,kmed,ksup) 

   





        end if
        nng=ng*kget
!- filtro la superficie libera
        mm = mod(jt,ift)
        if (mm.eq.0) then
          if(nsep.eq.1)then
            call doldfil1(yn,npc+2,npc+nng-2,npc+nng+1)
            call doldfil1(zn,npc+2,npc+nng-2,npc+nng+1)
          endif
          call doldfil1(ynsl,1+iford,npsl-10,npsl+1)
          call doldfil1(znsl,1+iford,npsl-10,npsl+1)
          call doldfil1(phinsl,1+iford,npsl-10,npsl)
        end if
          call nortan(vfall,yn,zn,amp,tmy,tmz,rny,rnz,npc,&
     &            ynsl,znsl,ampsl,tmysl,tmzsl,rnysl,rnzsl,npsl,&
     &            dphi,kget,ng,ycn,zcn,yce,zce,yv,zv,&
     &            ycnsl,zcnsl,ycsl,zcsl,ysl,zsl,phinsl,phisl,&
     &            phin,phi,ygb,zgb,1,ksep,jt)
!cckkk        write(90,*) '# jt ',jt,ng
!cckkk        do i=1,npc+ng
!cckkk          write(90,'(4d16.8,i4)') ycn(i),zcn(i),phin(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(90,*)
!cckkk        write(90,*)
!cckkk        write(91,*) '# jt ',jt
!cckkk        do i=1,npc+ng+1
!cckkk          write(91,*) yn(i),zn(i),ksep(i)
!cckkk        enddo
!cckkk        write(91,*)
!cckkk        write(91,*)
!cckkk        write(92,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(92,'(4d16.8)') ycnsl(i),zcnsl(i),phinsl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(92,*)
!cckkk        write(92,*)
!cckkk        write(93,*) '# jt ',jt
!cckkk        do i=1,npsl+1   
!cckkk          write(93,*) ynsl(i),znsl(i)
!cckkk        enddo
!cckkk        write(93,*)
!cckkk        write(93,*)
!cckkk        write(94,*) '# jt ',jt
!cckkk        do i=1,ngo1     
!cckkk          write(94,'(3d15.6,i4)') ygn(i),zgn(i),tg(i),i
!cckkk        enddo
!cckkk        write(94,*)
!cckkk        write(94,*)
!- trattamento getto 
        if(kget.eq.1)then
            call  get(ng,npc,ycsl,zcsl,ycn,zcn,yn,zn,  &
     &               xigs,zegs,xigb,zegb,xigf,zegf) 
        endif
        nng = kget*ng
! - chiamo il solutore:
! - prima ricalcolo dphtsl (non so se serve, dipende dalle BCs
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


!        write(*,*) 'Before Solver------JT',jt
!        write(*,*) 'Body-----,npc', npc
!        do ip =1,npc   
!        write(*,*) yce(ip), zce(ip), phi(ip), dphi(ip) 
!        enddo
!        write(*,*) 'FS--------------'
!        do ip = 1,10 
!        write(*,*) ycnsl(ip), zcnsl(ip), phisl(ip), dphisl(ip) 
!        enddo

        if(jt.le.jjget)then
          call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,&
     &                  phi,dphi,phisl,dphisl,jt)
        else
!          call solv2(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,
!     #           ycsl,zcsl,dphi,phisl,dphtsl,dphtbsl,phb,
!     #           xigs,zegs,xigb,zegb,xigf,zegf,
!     #           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
!     #           ph,dphn,dphnb,a1,b1,c1,d1,e1,
!     #           phi,phib,dphisl,dphibsl,rl,mb,mf,mt,m,n,nt,ntt,
!     #           xj,ze,xis,zes,jt)
          call solv22(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,&
     &           ycsl,zcsl,dphi,phisl,dphtsl,dphtbsl,phb,&
     &           xigs,zegs,xigb,zegb,xigf,zegf,&
     &           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,&
     &           ph,dphn,dphnb,a1,b1,c1,d1,e1,&
     &           phi,phib,dphisl,dphibsl,rl,mb,mf,mt,m,n,nt,ntt,&
     &           xj,ze,xis,zes,jt,ksep,kse,kord,kor)
          call calsol(mb,mf,mt,m,n,nt,ntt,phi,phib,phb,dphisl,&
     &                  dphnb,rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,&
     &                 ry,rz,kse,dphi)
        endif
!- velocita   
        call calvel(npsl,npc,ng,kget,mb,mf,mt,m,n,nt,ntt,&
     &                 phi,phib,phb,dpht,dphi,dphisl,dphnb,dphtsl,phisl,&
     &                 rl,a1,b1,c1,d1,e1,xj,ze,xis,zes,&
     &                 amp,vym1,vzm1,rny,rnz,tmy,tmz,&
     &                 ampsl,vymsl1,vzmsl1,rnysl,rnzsl,tmysl,tmzsl,vxi,&
     &                 ry,rz,ty,tz,kse)


       
        write(*,*) 'JT---------',jt
      do ip = 1,10 
        write(*,*) ycnsl(ip), zcnsl(ip), phisl(ip), dphisl(ip) 
        enddo
       
!
!cckkk        write(10,*) '# jt ',jt,ng
!cckkk        do i=1,npc+ng
!cckkk          write(10,'(4d16.8,i4)') yce(i),zce(i),phi(i),dphi(i),ksep(i)
!cckkk        enddo
!cckkk        write(10,*)
!cckkk        write(10,*)
!cckkk        write(11,*) '# jt ',jt
!cckkk        do i=1,npc+ng+1
!cckkk          write(11,*) yv(i),zv(i),ksep(i)
!cckkk        enddo
!cckkk        write(11,*)
!cckkk        write(11,*)
!cckkk        write(12,*) '# jt ',jt,ng
!cckkk        do i=1,npsl  
!cckkk          write(12,'(4d16.8)') ycsl(i),zcsl(i),phisl(i),dphisl(i)
!cckkk        enddo
!cckkk        write(12,*)
!cckkk        write(12,*)
!cckkk        write(13,*) '# jt ',jt
!cckkk        do i=1,npsl+1   
!cckkk          write(13,*) ysl(i),zsl(i)
!cckkk        enddo
!cckkk        write(13,*)
!cckkk        write(13,*)
!
! -------- LUNGHEZZA STRISCIA DI CONTROLLO -----
!
        if(jend.eq.0)then
          rrl=0.d0
      do i=jind+1,jfid-1
            dl= sqrt( (ycsl(i+1)-ycsl(i))**2 + (zcsl(i+1)-zcsl(i))**2 )
            rrl= rrl + dl
          enddo
          dlin= (1.d0-tin)*sqrt( (ycsl(jind+1)-ycsl(jind))**2 +&
     &                    (zcsl(jind+1)-zcsl(jind))**2 )
          dlfi= tfi*sqrt( (ycsl(jfid+1)-ycsl(jfid))**2 +&
     &                    (zcsl(jfid+1)-zcsl(jfid))**2 )
          rrl= rrl + dlin + dlfi
          yin = ycsl(jind)+tin*(ycsl(jind+1)-ycsl(jind))
          zin = zcsl(jind)+tin*(zcsl(jind+1)-zcsl(jind))
          yfi = ycsl(jfid)+tfi*(ycsl(jfid+1)-ycsl(jfid))
          zfi = zcsl(jfid)+tfi*(zcsl(jfid+1)-zcsl(jfid))
          write(77,'(i4,d15.6,d25.15,2i4,4d15.6)') jt,t,rrl,jind,jfid+1,&
     &               yin,zin,yfi,zfi                             
        endif
! ---------
!
!
!
! - PRESSIONE ------------------------------------------------------ 
! -- calcolo derivata tangenziale della velocita tangenziale
      nng   = ng*kget
      dphtu    = -vfall*san
      dphtd    = (phi(2)-phi(1))/(0.5d0*(amp(1)+amp(2)))
      dpht2(1) = ( dphtd - dphtu )/amp(1)
      dpnt(1)  = dpht2(1)*dphi(1)
      dpntt(1)  = dpht2(1)*dphi(1)
      i = 1 
      do i = 2,npc+nng-1
        amu = 0.5d0*(amp(i-1)+amp(i))
        amd = 0.5d0*(amp(i+1)+amp(i))
        dpht2u    = (dpht(i)-dpht(i-1))/amu
        dpht2d    = (dpht(i+1)-dpht(i))/amd
        dpht2(i)  = 0.5d0*(dpht2u+dpht2d)
        vxiu      = (vxi(i)-vxi(i-1))/amu
        vxid      = (vxi(i+1)-vxi(i))/amd
        dpht22    = 0.5d0*(vxiu+vxid)
        dpnt(i)   = dpht2(i)*dphi(i)
        if(i.ge.npc+1)then 
         dpntt(i)   = dpht22*dphi(i)
         else
         dpntt(i)   = dpnt(i)
        endif  
      enddo     
      dpht2(npc+nng)=(dpht(npc+nng)-dpht(npc+nng-1))/&
     &               (0.5d0*(amp(npc+nng)+amp(npc+nng-1)))
      dpnt(npc+nng) = dpht2(npc+nng)*dphi(npc+nng) 
      dpntt(npc+nng)= dpnt(npc+nng) 
! -- calcolo dphi/dt sulla SL
      do i = 1,npsl
        dptsl(i) = -(dphisl(i)**2 + dphtsl(i)**2)/2.d0
      enddo
      do i = 1,npc+nng 
        dpt(i)  = -(dphi(i)**2 + dpht(i)**2)/2.d0
        dpt2(i) = dpt(i)
      enddo
!
      call bcpre(ng,kget,vfall,npc,npsl,tmz,amp,phi,dphi,dpht,&
     &                 vxi,dphisl,dphtsl,yce,zce,dpnt,dpntt,dptsl,jt,tn,&
     &                 ygn,zgn,ygs2,zgs2,tg,ngo1,tc)
!
      if(kget.eq.0)then
        call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,&
     &                  dpt2,dpnt,dptsl,dpntsl,jt)
      else
!          call solv2(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,
!     #           ycsl,zcsl,dphi,phisl,dphtsl,dphtbsl,phb,
!     #           xigs,zegs,xigb,zegb,xigf,zegf,
!     #           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
!     #           ph,dphn,dphnb,a1,b1,c1,d1,e1,
!     #           phi,phib,dphisl,dphibsl,rl,mb,mf,mt,m,n,nt,ntt,
!     #           xj,ze,xis,zes,jt)
!        call solv2p(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,
!     #              ycsl,zcsl,dpnt,dptsl,dpttsl,dpttsl,dpb,
!     #              xigs,zegs,xigb,zegb,xigf,zegf,
!     #              ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,
!     #              ph,dphn,dphnb,a2,b2,c2,d2,e2,dpt,dptb,
!     #              dpntsl,dpntbsl,rl,mb,mf,mt,m,n,nt,ntt,xj,ze,xis,
!     #              zes,jt)
        call solv22p(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,&
     &              ycsl,zcsl,dpnt,dptsl,dpttsl,dpttsl,dpb,&
     &              xigs,zegs,xigb,zegb,xigf,zegf,&
     &              ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,&
     &              ph,dphn,dphnb,a2,b2,c2,d2,e2,dpt,dptb,&
     &              dpntsl,dpntbsl,rl,mb,mf,mt,m,n,nt,ntt,xj,ze,xis,&
     &              zes,jt,ksep,kse,kord,kor)
!
        call calsolp(mb,mf,mt,m,n,nt,ntt,dpt,dptb,dpb,dpntsl,dphnb,&
     &             rl,a2,b2,c2,d2,e2,xj,ze,xis,zes,ry,rz,kse,dpnt)
        call solv1(npc,nng,npsl,yv,zv,ygb,zgb,ysl,zsl,amp,ampsl,&
     &                  dpt2,dpnt,dptsl,dpntsl,jt)
      endif
!
! - calcolo con Dphi/Dt
      call prefin(npc,nng,npco,nngo,vfall,dt,dpht,dphi,tmy,tmz,&
     &                  rny,rnz,vyo,vzo,yce,zce,yco,zco,phi,phio,&
     &                  vym1,vzm1,pres)
! -- calcolo pressione e forza
      ff  = 0.d0
      ff2 = 0.d0
      do i=1,npc+nng
        p1 = -dpt(i)
        p11 = -dpt2(i)
        p2 = -(dphi(i)**2+dpht(i)**2)/2.d0
        pre(i) = p1+p2
        pre2(i) = p11+p2
        ff = ff + pre(i)*amp(i)*rnz(i)
        ff2 = ff2 + pres(i)*amp(i)*rnz(i)
!        write(41,'(i6,6d15.6)') i,yce(i),zce(i),p1,-p2,pre(i),pre2(i)
        write(41,'(i6,6f14.6)') i,yce(i),zce(i),p1,-p2,pre(i),pre2(i)
        if(i.eq.npc) write(41,*)
      enddo
      write(41,*)
      write(41,*)
      ff  = -2.d0*ff
      ff2 = -2.d0*ff2
      write(99,'(i10,3d15.7)') jt,t,ff,ff2
! FINE PRESSIONE -------------------------------------------------------
! - stampo
      if(mod(jt,ksta).eq.0.or.jt.le.200)then

      call stampa(vfall,ng,kget,npc,npsl,jt,t,dt,frdt,llf,&
     &                  scon,svel,spot,spre,phi,dphi,phisl,dphisl, &
     &                 yv,zv,yce,zce,ysl,zsl,ycsl,zcsl,&
     &                 vym1,vzm1,vymsl1,vzmsl1,&
     &                 dpht,dphtsl,dpt,dpnt,pre,pre2,pres,vxi,&
     &                 jend,jind,jfid,yin,zin,yfi,zfi)

         write(*,*) 'JT-----------------------------' , jt
         write(181,*) '# jt,ngo,ngo1,ng,nng',jt,ngo,ngo1,ng,nng

         if ((ngo1-ngo).eq.10) then

      do l = 1,ngo
             write(181,*) ygn(l),zgn(l)
           enddo
             write(181,*)
             write(181,*)
      do l = ngo,ngo1
             write(181,*) ygn(l),zgn(l)
           enddo
             close(181)
                stop 



         end if

      endif
    
!
      nng = ng*kget
!      if(jt.eq.jjget+1)stop
! ------------------------------------  CONTROLLI
!c - B+
!      write(31,*) '# ',jt,zv(1)
!      do i=1,npc+1
!        write(31,'(2d15.7)') yv(i),zv(i)
!      enddo
!      write(31,*)  
!      write(31,*)  
!      write(32,*) '# ',jt
!      do i=1,npc
!        write(32,'(8d15.7)') yce(i),
!     #       zce(i),phi(i),dphi(i),
!     #            dpht(i),vym2(i),vzm2(i) ,amp(i)  
!      enddo
!      write(32,*)
!      do i =1,nng
!        write(32,'(8d15.7)') 0.5d0*(ygb(i)+ygb(i+1)),
!     #        0.5d0*(zgb(i)+zgb(i+1)),phi(npc+i),dphi(npc+i),
!     #            dpht(npc+i),vym2(npc+i),vzm2(npc+i) ,amp(npc+i)  
!      enddo
!      write(32,*)  
!      write(32,*)  
!      write(33,*) '# ',jt
!      do i=1,npsl+1
!        write(33,'(4d15.7)') ysl(i),zsl(i)
!      enddo
!      write(33,*)  
!      write(33,*)  
!      write(34,*) '# ',jt
!      do i=1,npsl
!        write(34,'(8d15.7)') 0.5d0*(ysl(i)+ysl(i+1)),
!     #        0.5d0*(zsl(i)+zsl(i+1)),phisl(i),dphisl(i),
!     #            dphtsl(i),vymsl2(i),vzmsl2(i) ,ampsl(i)  
!      enddo
!      write(34,*)  
!      write(34,*)
!-----------------------------------------------  FINE CONTROLLI
      enddo
!
      stop
      end
