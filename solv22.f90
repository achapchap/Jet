
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine solv22(frint,ng,npc,npsl,yce,zce,yv,zv,ysl,zsl,&
     &           ycsl,zcsl,dphi,phisl,dphtsl,dphtbsl,phb,&
     &           xigs,zegs,xigb,zegb,xigf,zegf,&
     &           ty,tz,ry,rz,tmy,tmz,rny,rnz,tmysl,tmzsl,rnysl,rnzsl,&
     &           ph,dphn,dphnb,a,b,c,d,e,phi,phib,dphisl,dphibsl,rl,    &
     &           mb,mf,mt,m,n,nt,ntt,xi,ze,xis,zes,jt,ksep,kse,kord,kor)
!
      include"slam_p.h"
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
!
!      dphn dphi yc zc yce zce yb zb yf zf yv zv ysl zsl ycsl zcsl
!      ph dphtb phisl
!      dh xig xigs zeb zebs zef zefs hp xi xis ze zes 
!      am aa ff bb rl
!      phb dphnb a b c d e 
      write(*,*) '--> solv2.............'
!
! - riorganizzo array punti yf yb zy zb yc zc am e phi dphn dpht
!                           dphnb = deriv normale bulk
!                           dphtb = deriv tangenz bulk
!                           phb   = phi nel bulk
      m  = int(frint*ng)
      n  = ng - m
      mb = npc+ng-(m+n)
      mf = npsl-(m+n)
      mt = mb + mf
      nt = mt + 2*m
      ntt= nt + 2*n
      write(*,*) npc,ng,npsl
      write(*,*) mb,m,n
      write(*,*) mf,m,n
      write(*,*) mt,nt,ntt
!
      do i=1,mb+m+n
        if(i.le.mb)then
          dphn(i) = dphi(i)
          dphnb(i)= dphi(i)
          ph(i)   = phi(i)
          phb(i)  = phi(i)
          kse(i)  = ksep(i)
          kor(i)  = kord(i)
!          yc(i)   = yce(i)
!          zc(i)   = zce(i)
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
!          write(71,'(i6,10d15.7)') i,yc(i),zc(i),
!     #                yb(i),zb(i),yf(i),zf(i), 
!     #                ry(i),rz(i),ty(i),tz(i) 
        elseif(i.ge.mb+1.and.i.le.mb+m)then
          dphn(mt+i-mb) = dphi(i)
          dphnb(mt+i-mb)= dphi(i)
          ph(mt+i-mb)   = phi(i)
          phb(mt+i-mb)  = phi(i)
          kse(mt+i-mb)  = ksep(i)
          kor(mt+i-mb)  = kord(i)
!          yc(mt+i-mb ) = yce(i)
!          zc(mt+i-mb)  = zce(i)
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
!          write(72,'(i6,10d15.7)') mt+i-mb,yc(mt+i-mb),zc(mt+i-mb),
!     #                yb(mt+i-mb),zb(mt+i-mb),yf(mt+i-mb),zf(mt+i-mb),
!     #                dphn(mt+i-mb),dphnb(mt+i-mb)  
!     #                ry(mt+i-mb),rz(mt+i-mb),ty(mt+i-mb),tz(mt+i-mb) 
                       
        else
          dphn(nt+i-(mb+m)) = dphi(i) 
          dphnb(nt+i-(mb+m))= dphi(i) 
          ph(nt+i-(mb+m))   = phi(i) 
          phb(nt+i-(mb+m))  = phi(i) 
          kse(nt+i-(mb+m))  = ksep(i) 
          kor(nt+i-(mb+m))  = kord(i) 
!          yc(nt+i-(mb+m))  = yce(i)
!          zc(nt+i-(mb+m))  = zce(i)
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
!          write(73,'(i6,10d15.7)')
!     #    nt+i-(mb+m),yc(nt+i-(mb+m)),zc(nt+i-(mb+m)), 
!     #  yb(mt+i-(mb+m)),zb(mt+i-(mb+m)),yf(mt+i-(mb+m)),zf(mt+i-(mb+m)),
!     #           dphn(nt+i-(mb+m)),dphnb(nt+i-(mb+m))   
!     #  ry(nt+i-(mb+m)),rz(nt+i-(mb+m)),ty(nt+i-(mb+m)),tz(nt+i-(mb+m)) 
        endif  
      enddo    
!      write(71,*)
!      write(71,*)
!      write(72,*)
!      write(72,*)
!      write(73,*)
!      write(73,*)
      do i=1,mf+m+n
        if(i.le.n)then
          dphtbsl(i)   = dphtsl(i)
          dphtb(nt+n+i)= dphtbsl(i)
!          dpht(nt+n+i)= dphtsl(i)
          ph(nt+n+i)   = phisl(i)
          phb(nt+n+i)  = phisl(i)
!          yc(nt+n+i)  = ycsl(i)
!          zc(nt+n+i)  = zcsl(i)
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
          write(87,'(i6,10d15.7)') nt+n+i,yc(nt+n+i),zc(nt+n+i),&
     &                yb(nt+n+i),zb(nt+n+i),yf(nt+n+i),zf(nt+n+i),    &
!     #                dphtb(nt+n+i),ph(nt+n+i),phb(nt+n+i) 
     &             ry(nt+n+i),rz(nt+n+i),ty(nt+n+i),tz(nt+n+i)    
        elseif(i.ge.n+1.and.i.le.n+m)then
          dphtbsl(i)       = dphtsl(i)
          dphtb(mt+m+(i-n))= dphtbsl(i)
!          dpht(mt+m+(i-n)) = dphtsl(i)
          ph(mt+m+(i-n))   = phisl(i)
          phb(mt+m+(i-n))   = phisl(i)
!          yc(mt+m+(i-n))   = ycsl(i)
!          zc(mt+m+(i-n))   = zcsl(i)
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
!         write(82,'(i6,10d15.7)') mt+m+i-n,yc(mt+m+i-n),zc(mt+m+i-n),  
!     #             yb(mt+m+i-n),zb(mt+m+i-n),yf(mt+m+i-n),zf(mt+m+i-n), 
!     #            dphtb(mt+m+i-n), ph(mt+m+i-n), phb(mt+m+i-n)     
!     #          ry(mt+m+i-n),rz(mt+m+i-n),ty(mt+m+i-n),tz(mt+m+i-n)    
        else
          dphtbsl(i)       = dphtsl(i)
          dphtb(mb+i-(n+m))= dphtbsl(i)
!          dpht(mb+i-(n+m)) = dphtsl(i)
          ph(mb+i-(n+m))   = phisl(i)
          phb(mb+i-(n+m))   = phisl(i)
!          yc(mb+i-(n+m))   = ycsl(i)
!          zc(mb+i-(n+m))   = zcsl(i)
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
!          write(89,'(i6,10d15.7)')
!     #        mb+i-(n+m),yc(mb+i-(n+m)),zc(mb+i-(m+n)),  
!     #      yb(mb+i-(n+m)),zb(mb+i-(m+n)),yf(mb+i-(n+m)),zf(mb+i-(m+n)),
!     #       dphtb(mb+i-(m+n)), ph(mb+i-(m+n)), phb(mb+i-(m+n))     
!     #      ry(mb+i-(n+m)),rz(mb+i-(n+m)),ty(mb+i-(n+m)),tz(mb+i-(n+m))
        endif
      enddo
!      write(89,*)
!      write(89,*)
!      write(82,*)
!      write(82,*)
      write(87,*)
      write(87,*)
! - ampiezze
      do i=1,ntt
        am(i)=sqrt( (yf(i)-yb(i))**2 + (zf(i)-zb(i))**2 ) 
      enddo
! - getto
      write(75,*) '# jt ',jt
      write(76,*) '# jt ',jt
      i0 = ng-(m+n)
      do i=i0-1,ng
        if(i.le.i0)then
          ii = mb+i-i0
          jj = mb-(i-i0)+1 
        elseif(i.ge.i0+1.and.i.le.i0+m)then
          ii = mt+i-i0
          jj = nt-(i-i0)+1
        else
          ii = nt+(i-i0-m)
          jj = ntt-(i-i0-m)+1
        endif
!        hp(ii) = dh(i)
        xi(ii)  = xigb(i)
        xis(ii) = xigs(i)
        ze(ii)  = zegb(i)
        zes(ii) = zegs(i)
!        hp(jj)  = dh(i)
        xi(jj)  = xigf(i)
        xis(jj) = xigs(i)
        ze(jj)  = zegf(i)
        zes(jj) = zegs(i)
        write(75,'(2i6,4d15.5)') ii,i,xi(ii),xis(ii),&
     &              ze(ii),zes(ii)
        write(76,'(2i6,4d15.5)') jj,i,xi(jj),xis(jj),&
     &              ze(jj),zes(jj)
      enddo
      write(75,*)
      write(75,*)
      write(75,*) 
      write(76,*)
      write(76,*)
! - fattore di matching rl non va con n=0 m>1
      do i=1,nt
        rl(i)=0.d0
      enddo
      do i=nt+1,ntt
        rl(i)=1.d0
      enddo
      do i=1,m
        ii  = mt+i
        if(i.eq.1)then
          dam = 0.5d0*(am(ii)+am(mb))
        else
          dam = 0.5d0*(am(ii)+am(ii-1))
        endif
        rl(ii)   = rl(ii-1) + dam 
      enddo
      if(n.eq.0)then 
        dam = am(mt+m)
      else
        dam = 0.5d0*( am(mt+m) + am(nt+1) )
      endif
      amt = rl(ii) + dam
      write(67,*) ' # jt ',jt
      do i=1,m
        rl(mt+i)   = rl(mt+i)/amt 
        rl(nt-i+1) = rl(mt+i)
        write(67,'(3i6,d15.6)') i,mt+i,nt-i+1,rl(mt+i)
      enddo
      write(67,*)
      write(67,*)
!      
!
!
! - costituisco la matrice:
!
!1 - inizializzo matrice
      npe = nt+5*(n+m)
      do i=1,npe
        ff(i) = 0.d0
        bb(i) = 0.d0
      do j=1,npe
        aa(i,j) = 0.d0
      enddo
      enddo
!2 - parte BEM, righe =1,nt
! - qui i e indice di riga  NB: ntt = NT negli appunti
      do i = 1,nt
        yy = yc(i)
        zz = zc(i)
      do k = 1,ntt
          yp = yb(k)
          ys = yf(k)
          zp = zb(k)
          zs = zf(k)
          yps= -yf(k)
          yss= -yb(k)
          zps= zf(k)
          zss= zb(k)
! - coeff influenza, pannello + simmetrico rispetto a y=0
          call finte(yy,zz,yp,zp,ys,zs,am(k),fint1,fint2,0.d0)
          call finte(yy,zz,yps,zps,yss,zss,am(k),fins1,fins2,0.d0)
          gik  =  (fint1+fins1)/(2.d0*pi)
          dgik = -(fint2+fins2)/(2.d0*pi)
          if(k.le.mb)then
            if(kse(k).ne.1)then
            aa(i,k) = aa(i,k) + dgik    !/1000.d0
            ff(i)   = ff(i) + dphn(k)*gik
            else
!            write(98,*) 'ecc1!', k,nt+n,jt
            aa(i,k) = aa(i,k) - gik
            ff(i)   = ff(i) - ph(k)*dgik
            endif 
          elseif(k.ge.mb+1.and.k.le.mt)then
            aa(i,k) = aa(i,k) - gik
            ff(i)   = ff(i) - ph(k)*dgik 
          elseif(k.ge.mt+1.and.k.le.mt+m)then
            rlk  = rl(k)
!            hpk  = hp(k)
!            sqh  = sqrt(1.d0+hpk**2)
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
            aa(i,k)                = aa(i,k)+ (1.d0-rlk)*dgik !/1000.d0 
            aa(i,nt+k-mt)          = aa(i,nt+k-mt)         + rlk*aik
            aa(i,nt+(n+m)+k-mt)    = aa(i,nt+(n+m)+k-mt)   + rlk*bik
            aa(i,nt+2*(n+m)+k-mt)  = aa(i,nt+2*(n+m)+k-mt) + rlk*cik
            aa(i,nt+3*(n+m)+k-mt)  = aa(i,nt+3*(n+m)+k-mt) + rlk*dik
            aa(i,nt+4*(n+m)+k-mt)  = aa(i,nt+4*(n+m)+k-mt) + rlk*eik
            ff(i)  = ff(i) + dphn(k)*gik
            else
!            write(98,*) 'eccolo0!', k,nt+n,jt
            ryk  = ry(k)
            rzk  = rz(k)
            bik  = ryk*gik
            cik  = rzk*gik
            dik  = gik*( ryk*(xik-xisk) - rzk*(zek-zesk) )  
            eik  = gik*( ryk*(zek-zesk) + rzk*(xik-xisk) )   
            aa(i,nt+(n+m)+k-mt)    = aa(i,nt+(n+m)+k-mt)   - rlk*bik
            aa(i,nt+2*(n+m)+k-mt)  = aa(i,nt+2*(n+m)+k-mt) - rlk*cik
            aa(i,nt+3*(n+m)+k-mt)  = aa(i,nt+3*(n+m)+k-mt) - rlk*dik
            aa(i,nt+4*(n+m)+k-mt)  = aa(i,nt+4*(n+m)+k-mt) - rlk*eik
            aa(i,k) = aa(i,k)-(1.d0-rlk)*gik
            ff(i)   = ff(i) - ph(k)*dgik 
            endif
          elseif(k.ge.mt+m+1.and.k.le.nt)then
            kk   = mt+(nt-k)+1
            rlk  = rl(k)
!            hpk  = hp(k)
!            sqh  = sqrt(1.d0+hpk**2)
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
!           bik  = + hpk/sqh*gik
!           cik  = - 1.d0/sqh*gik 
!            dik  = + gik/sqh*( hpk*(xik-xiskk) + (zek-zeskk) )
!            eik  = + gik/sqh*( hpk*(zek-zeskk) - (xik-xiskk) )
            aa(i,nt+(n+m)+kk-mt)    = aa(i,nt+(n+m)+kk-mt)   - rlk*bik
            aa(i,nt+2*(n+m)+kk-mt)  = aa(i,nt+2*(n+m)+kk-mt) - rlk*cik
            aa(i,nt+3*(n+m)+kk-mt)  = aa(i,nt+3*(n+m)+kk-mt) - rlk*dik
            aa(i,nt+4*(n+m)+kk-mt)  = aa(i,nt+4*(n+m)+kk-mt) - rlk*eik
            aa(i,k) = aa(i,k)-(1.d0-rlk)*gik
            ff(i)   = ff(i) - ph(k)*dgik 
          elseif(k.ge.nt+1.and.k.le.nt+n)then
!            hpk  = hp(k)
!            sqh  = sqrt(1.d0+hpk**2)
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
            aa(i,nt+m+k-nt)          = aa(i,nt+m+k-nt)         + aik
            aa(i,nt+(n+m)+m+k-nt)    = aa(i,nt+(n+m)+m+k-nt)   + bik
            aa(i,nt+2*(n+m)+m+k-nt)  = aa(i,nt+2*(n+m)+m+k-nt) + cik
            aa(i,nt+3*(n+m)+m+k-nt)  = aa(i,nt+3*(n+m)+m+k-nt) + dik
            aa(i,nt+4*(n+m)+m+k-nt)  = aa(i,nt+4*(n+m)+m+k-nt) + eik
            ff(i)  = ff(i) + dphn(k)*gik
            else
!            write(98,*) 'eccolo1!', k,nt+n,jt,ph(k)
            ryk  = ry(k)
            rzk  = rz(k)
!            write(98,*) 'ph rn',ph(k),ryk,rzk 
            bik  = ryk*gik
            cik  = rzk*gik
            dik  = gik*( ryk*(xik-xisk) - rzk*(zek-zesk) )  
            eik  = gik*( ryk*(zek-zesk) + rzk*(xik-xisk) )   
            aa(i,nt+(n+m)+m+k-nt)  = aa(i,nt+(n+m)+m+k-nt)   - bik
            aa(i,nt+2*(n+m)+m+k-nt)= aa(i,nt+2*(n+m)+m+k-nt) - cik
            aa(i,nt+3*(n+m)+m+k-nt)= aa(i,nt+3*(n+m)+m+k-nt) - dik
            aa(i,nt+4*(n+m)+m+k-nt)= aa(i,nt+4*(n+m)+m+k-nt) - eik
            ff(i)   = ff(i) - ph(k)*dgik 
            endif
          elseif(k.ge.nt+n+1)then
            kk   = nt+(ntt-k)+1
!            hpk  = hp(k)
!            sqh  = sqrt(1.d0+hpk**2)
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
!            bik  = + hpk/sqh*gik
!            cik  = - 1.d0/sqh*gik 
!            dik  = + gik/sqh*( hpk*(xik-xiskk) + (zek-zeskk) )
!            eik  = + gik/sqh*( hpk*(zek-zeskk) - (xik-xiskk) )
            aa(i,nt+(n+m)+m+kk-nt)  = aa(i,nt+(n+m)+m+kk-nt)   - bik
            aa(i,nt+2*(n+m)+m+kk-nt)= aa(i,nt+2*(n+m)+m+kk-nt) - cik
            aa(i,nt+3*(n+m)+m+kk-nt)= aa(i,nt+3*(n+m)+m+kk-nt) - dik
            aa(i,nt+4*(n+m)+m+kk-nt)= aa(i,nt+4*(n+m)+m+kk-nt) - eik
            ff(i)   = ff(i) - ph(k)*dgik 
          endif
        enddo
! - termini autoindotti
        if (i.le.mb) then
          if(kse(i).ne.1)then 
          aa(i,i) = aa(i,i) + 0.5d0  !/1000.d0 
          else
!            write(98,*) 'ecc5!',i ,jt,ph(i)
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
          aa(i,nt+(m+n)+i-mt) = aa(i,nt+(m+n)+i-mt) + &
     &                                             rli*0.5d0*(xii-xisi) 
          aa(i,nt+2*(m+n)+i-mt)= aa(i,nt+2*(m+n)+i-mt) + &
     &                                            rli*0.5d0*(zei-zesi) 
          aa(i,nt+3*(m+n)+i-mt)= aa(i,nt+3*(m+n)+i-mt) + &
     &                          rli*0.25d0*((xii-xisi)**2-(zei-zesi)**2)
          aa(i,nt+4*(m+n)+i-mt)=aa(i,nt+4*(m+n)+i-mt) + &
     &                          rli*0.5d0*(xii-xisi)*(zei-zesi)
          else
!            write(98,*) 'ecc6!',i ,jt,ph(i)
          ff(i)   = ff(i) - 0.5d0*ph(i)
          endif
        elseif(i.ge.mt+m+1.and.i.le.nt)then
          ff(i)   = ff(i) - 0.5d0*ph(i)
        endif 
      enddo
!
!3 - parte getto ,  righe =nt+1,nt+5(m+n)
! - NB : qui i non e indice di riga, bensi indice di pannello (cella) del getto
!
      do 100, i = mt+1, nt+n ! Was loop 100
!
        if(i.ge.mt+1.and.i.le.mt+m)then 
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
!
        elseif(i.ge.nt+1.and.i.le.nt+n)then
!
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
!
        else 
!
          goto 100
        endif
!
        rli  = rl(i)
        rli1 = rl(i1)
        drl  = rli1-rli
        dxif = (xi(i)-xi(i1))     !*1000.d0
        dxib = (xi(i1)-xi(i2))    !*1000.d0
        xii  = xi(i)
        xii1 = xi(i1)
        xisi = xis(i)
        xisi1= xis(i1) 
        zei  = ze(i)
        zei1 = ze(i1)
        zesi = zes(i)
        zesi1= zes(i1) 
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
!        write(99,'(''drl '',3i6,7d15.6)')i,i1,i2,rli,
!     #                      ryj,rzj,ryi1,rzi1,ryi,rzi
!        write(98,'(''drl '',3i6,5d15.6)')i,i1,i2,drl,zesi,zesi1,zei,zei1
!        write(97,'(''drl '',3i6,5d15.6)')i,i1,i2,drl,xisi,xisi1,xij,xij1
!        write(96,'(''drl '',3i6,5d15.6)')i,i1,i2,drl,zesi,zesi1,zej,zej1
! - continuita phi sul corpo (i1)
!        aa(iii,i1)         = aa(iii,i1)         + drl
!        aa(iii,iij)        = aa(iii,iij)        + rli  
!        aa(iii,iij-1)      = aa(iii,iij-1)      - rli1 
!        aa(iii,iij+mn)     = aa(iii,iij+mn)     + rli*(xii1-xisi) 
!        aa(iii,iij-1+mn)   = aa(iii,iij-1+mn)   - rli1*(xii1-xisi1) 
!        aa(iii,iij+2*mn)   = aa(iii,iij+2*mn)   + rli*(zei1-zesi) 
!        aa(iii,iij-1+2*mn) = aa(iii,iij-1+2*mn) - rli1*(zei1-zesi1) 
!        aa(iii,iij+3*mn)   = aa(iii,iij+3*mn)   + 
!     #                       0.5d0*rli*((xii1-xisi)**2-(zei1-zesi)**2)
!        aa(iii,iij-1+3*mn) = aa(iii,iij-1+3*mn) - 
!     #                      0.5d0*rli1*((xii1-xisi1)**2-(zei1-zesi1)**2)
!        aa(iii,iij+4*mn)   = aa(iii,iij+4*mn)   +
!     #                      rli*(xii1-xisi)*(zei1-zesi)
!        aa(iii,iij-1+4*mn) = aa(iii,iij-1+4*mn) - 
!     #                      rli1*(xii1-xisi1)*(zei1-zesi1)
!- continuita vtau sul corpo (i1)  NON FUNZIONA
!        if(i.eq.nt+1)then
!          aa(iii,i2)         = aa(iii,i2)         + 0.5d0/dxib
!          aa(iii,i1)         = aa(iii,i1)         +
!     #                                  0.5d0/dxif-0.5d0/dxib
!          aa(iii,iij)       = aa(iii,iij)       - 0.5d0/dxif 
!          aa(iii,iij+mn)    = aa(iii,iij+mn)    - 0.5d0/dxif*(xii1-xisi)
!          aa(iii,iij+2*mn)  = aa(iii,iij+2*mn)  - 0.5d0/dxif*(zei1-zesi)
!          aa(iii,iij+3*mn)  = aa(iii,iij+3*mn)  -
!     #                  0.5d0/dxif*0.5d0*((xii1-xisi)**2-(zei1-zesi)**2)
!          aa(iii,iij+4*mn)  = aa(iii,iij+4*mn)  -
!     #                      0.5d0/dxif*(xii1-xisi)*(zei1-zesi)
!          aa(iii,iij+mn)      = aa(iii,iij+mn)      + rli 
!          aa(iii,iij+3*mn)    = aa(iii,iij+3*mn)    + rli*(xii1-xisi) 
!          aa(iii,iij+4*mn)    = aa(iii,iij+4*mn)    + rli*(zei1-zesi) 
!          aa(iii,i2)         = aa(iii,i2)         + 0.5d0/dxib
!          aa(iii,i1)         = aa(iii,i1)         - 0.5d0/dxib 
!        elseif(i.ge.mt+1.and.i.le.mt+m)then
!          aa(iii,i)          = aa(iii,i)          - 0.5d0/dxif
!        else 
!          aa(iii,iij+mn)      = aa(iii,iij+mn)      + rli 
!          aa(iii,iij+3*mn)    = aa(iii,iij+3*mn)    + rli*(xii1-xisi) 
!          aa(iii,iij+4*mn)    = aa(iii,iij+4*mn)    + rli*(zei1-zesi) 
!          aa(iii,iij-1+mn)    = aa(iii,iij-1+mn)    - rli1 
!          aa(iii,iij-1+3*mn)  = aa(iii,iij-1+3*mn)  - rli1*(xii1-xisi1) 
!          aa(iii,iij-1+4*mn)  = aa(iii,iij-1+4*mn)  - rli1*(zei1-zesi1) 
!        endif
! - continuita phi sulla SL (j)
!        aa(iii+k1,iij)       = aa(iii+k1,iij)        + rli  
!        aa(iii+k1,iij-1)     = aa(iii+k1,iij-1)      - rli1 
!        aa(iii+k1,iij+mn)    =aa(iii+k1,iij+mn)    +rli*(xij-xisi) 
!        aa(iii+k1,iij-1+mn)  =aa(iii+k1,iij-1+mn)  -rli1*(xij-xisi1)
!        aa(iii+k1,iij+2*mn)  =aa(iii+k1,iij+2*mn)  +rli*(zej-zesi) 
!        aa(iii+k1,iij-1+2*mn)=aa(iii+k1,iij-1+2*mn)-rli1*(zej-zesi1)
!        aa(iii+k1,iij+3*mn)  =aa(iii+k1,iij+3*mn)  + 
!     #                       0.5d0*rli*((xij-xisi)**2-(zej-zesi)**2)
!        aa(iii+k1,iij-1+3*mn) = aa(iii+k1,iij-1+3*mn) - 
!     #                      0.5d0*rli1*((xij-xisi1)**2-(zej-zesi1)**2)
!        aa(iii+k1,iij+4*mn)   = aa(iii+k1,iij+4*mn)   +
!     #                           rli*(xij-xisi)*(zej-zesi)
!        aa(iii+k1,iij-1+4*mn) = aa(iii+k1,iij-1+4*mn) - 
!     #                           rli1*(xij-xisi1)*(zej-zesi1)
!        ff(iii+ k1 )  = ff(iii+k1)   + 0.d0 
!        ff(iii+ k1 )  = ff(iii+k1)   - drl*ph(j)
! 
! - phi sulla SL  (j)
        aa(iii+k1,iij)       = aa(iii+k1,iij)       + 1.d0 
        aa(iii+k1,iij+mn)    = aa(iii+k1,iij+mn)    + (xij-xisi)
        aa(iii+k1,iij+2*mn)  = aa(iii+k1,iij+2*mn)  + (zej-zesi)
        aa(iii+k1,iij+3*mn)  = aa(iii+k1,iij+3*mn)  +&
     &                             0.5d0*((xij-xisi)**2-(zej-zesi)**2)
        aa(iii+k1,iij+4*mn)  = aa(iii+k1,iij+4*mn)  +&
     &                             (xij-xisi)*(zej-zesi)
        ff(iii+k1)  = ff(iii+k1) + ph(j)       
! - phi sulla SL (j1)
        aa(iii+4*k1,iij)       = aa(iii+4*k1,iij)       + 1.d0 
        aa(iii+4*k1,iij+mn)    = aa(iii+4*k1,iij+mn)    + (xij1-xisi)
        aa(iii+4*k1,iij+2*mn)  = aa(iii+4*k1,iij+2*mn)  + (zej1-zesi)
        aa(iii+4*k1,iij+3*mn)  = aa(iii+4*k1,iij+3*mn)  +&
     &                             0.5d0*((xij1-xisi)**2-(zej1-zesi)**2)
        aa(iii+4*k1,iij+4*mn)  = aa(iii+4*k1,iij+4*mn)  +&
     &                             (xij1-xisi)*(zej1-zesi)
        ff(iii+4*k1)  = ff(iii+4*k1) + ph(j1)       
!
        if(kor(i).eq.2)then       
!
! - continuita vn sulla SL (j)
        aa(iii,j)         = aa(iii,j)         + drl           
        aa(iii,iij+mn)    = aa(iii,iij+mn)    + rli*ryj 
        aa(iii,iij-1+mn)  = aa(iii,iij-1+mn)  - rli1*ryj 
        aa(iii,iij+2*mn)  = aa(iii,iij+2*mn)  + rli*rzj 
        aa(iii,iij-1+2*mn)= aa(iii,iij-1+2*mn)- rli1*rzj 
        aa(iii,iij+3*mn)  = aa(iii,iij+3*mn)  +&
     &                         rli*( ryj*(xij-xisi) - rzj*(zej-zesi) ) 
        aa(iii,iij-1+3*mn)= aa(iii,iij-1+3*mn)- &
     &                        rli1*( ryj*(xij-xisi1) - rzj*(zej-zesi1) )
        aa(iii,iij+4*mn)  = aa(iii,iij+4*mn)  +&
     &                         rli*( ryj*(zej-zesi) + rzj*(xij-xisi) ) 
        aa(iii,iij-1+4*mn)= aa(iii,iij-1+4*mn)- &
     &                        rli1*( ryj*(zej-zesi1)+ rzj*(xij-xisi1) )
        ff(iii)       = ff(iii)      + 0.0d0
! - vn sul corpo (i1) o phi sulla SL (se separa)
        if(kse(i1).ne.1)then
        aa(iii+2*k1,iij+mn)   = aa(iii+2*k1,iij+mn)    + ryi1 
        aa(iii+2*k1,iij+2*mn) = aa(iii+2*k1,iij+2*mn)  + rzi1 
        aa(iii+2*k1,iij+3*mn) = aa(iii+2*k1,iij+3*mn)  +&
     &                            ryi1*(xii1-xisi) - rzi1*(zei1-zesi)
        aa(iii+2*k1,iij+4*mn) = aa(iii+2*k1,iij+4*mn)  +&
     &                            ryi1*(zei1-zesi) + rzi1*(xii1-xisi)
        ff(iii+2*k1)  = ff(iii+2*k1) + dphn(i1)
        else
!        write(98,*) 'eccolo10!',i1,nt+n,jt
!        write(98,*)  xii1,zei1,ph(i1)
        aa(iii+2*k1,iij)       = aa(iii+2*k1,iij)       + 1.d0 
        aa(iii+2*k1,iij+mn)    = aa(iii+2*k1,iij+mn)    + (xii1-xisi)
        aa(iii+2*k1,iij+2*mn)  = aa(iii+2*k1,iij+2*mn)  + (zei1-zesi)
        aa(iii+2*k1,iij+3*mn)  = aa(iii+2*k1,iij+3*mn)  +&
     &                             0.5d0*((xii1-xisi)**2-(zei1-zesi)**2)
        aa(iii+2*k1,iij+4*mn)  = aa(iii+2*k1,iij+4*mn)  +&
     &                             (xii1-xisi)*(zei1-zesi)
        ff(iii+2*k1)  = ff(iii+2*k1) + ph(i1)       
        endif
! - vn sul corpo (i)o phi sulla SL (se separa) 
        if(kse(i).ne.1)then
        aa(iii+3*k1,iij+mn)   = aa(iii+3*k1,iij+mn)    + ryi 
        aa(iii+3*k1,iij+2*mn) = aa(iii+3*k1,iij+2*mn)  + rzi 
        aa(iii+3*k1,iij+3*mn) = aa(iii+3*k1,iij+3*mn)  +&
     &                            ryi*(xii-xisi) - rzi*(zei-zesi)
        aa(iii+3*k1,iij+4*mn) = aa(iii+3*k1,iij+4*mn)  +&
     &                            ryi*(zei-zesi) + rzi*(xii-xisi)
        ff(iii+3*k1)  = ff(iii+3*k1) + dphn(i)
        else
!        write(98,*) 'eccolo11!',i,nt+n,jt
!        write(98,*) xii,zei,ph(i)
        aa(iii+3*k1,iij)       = aa(iii+3*k1,iij)       + 1.d0 
        aa(iii+3*k1,iij+mn)    = aa(iii+3*k1,iij+mn)    + (xii-xisi)
        aa(iii+3*k1,iij+2*mn)  = aa(iii+3*k1,iij+2*mn)  + (zei-zesi)
        aa(iii+3*k1,iij+3*mn)  = aa(iii+3*k1,iij+3*mn)  +&
     &                             0.5d0*((xii-xisi)**2-(zei-zesi)**2)
        aa(iii+3*k1,iij+4*mn)  = aa(iii+3*k1,iij+4*mn)  +&
     &                             (xii-xisi)*(zei-zesi)
        ff(iii+3*k1)  = ff(iii+3*k1) + ph(i)       
        endif
!
        else
!
! - vn sul corpo (i1) o continuita di phi sul corpo (se separa) 
        if(kse(i1).ne.1)then
        aa(iii,iij+mn)   = aa(iii,iij+mn)    + ryi1 
        aa(iii,iij+2*mn) = aa(iii,iij+2*mn)  + rzi1 
        aa(iii,iij+3*mn) = aa(iii,iij+3*mn)  +&
     &                            ryi1*(xii1-xisi) - rzi1*(zei1-zesi)
        aa(iii,iij+4*mn) = aa(iii,iij+4*mn)  +&
     &                            ryi1*(zei1-zesi) + rzi1*(xii1-xisi)
        ff(iii)  = ff(iii) + dphn(i1)
        else
        aa(iii,iij)        = aa(iii,iij)        + rli  
        aa(iii,iij-1)      = aa(iii,iij-1)      - rli1 
        aa(iii,iij+mn)     = aa(iii,iij+mn)     + rli*(xii1-xisi) 
        aa(iii,iij-1+mn)   = aa(iii,iij-1+mn)   - rli1*(xii1-xisi1) 
        aa(iii,iij+2*mn)   = aa(iii,iij+2*mn)   + rli*(zei1-zesi) 
        aa(iii,iij-1+2*mn) = aa(iii,iij-1+2*mn) - rli1*(zei1-zesi1) 
        aa(iii,iij+3*mn)   = aa(iii,iij+3*mn)   + &
     &                       0.5d0*rli*((xii1-xisi)**2-(zei1-zesi)**2)
        aa(iii,iij-1+3*mn) = aa(iii,iij-1+3*mn) - &
     &                      0.5d0*rli1*((xii1-xisi1)**2-(zei1-zesi1)**2)
        aa(iii,iij+4*mn)   = aa(iii,iij+4*mn)   +&
     &                      rli*(xii1-xisi)*(zei1-zesi)
        aa(iii,iij-1+4*mn) = aa(iii,iij-1+4*mn) - &
     &                      rli1*(xii1-xisi1)*(zei1-zesi1)
        ff(iii)  = ff(iii) - drl*ph(i1) 
        endif
! - nullita dei coefficienti d ed e
        aa(iii+2*k1,iij+3*mn) = aa(iii+2*k1,iij+3*mn)  + 1.d0 
        aa(iii+3*k1,iij+4*mn) = aa(iii+3*k1,iij+4*mn)  + 1.d0 
        ff(iii+2*k1)  = ff(iii+2*k1) + 0.d0 
        ff(iii+3*k1)  = ff(iii+3*k1) + 0.d0 
!
        endif
!
!
 100  continue
!      end do ! Was loop 100
!    
      do i=1,npe
        bb(i) = ff(i)
      enddo
!
!      if(jt.eq.252)then
!      do i=1,npe
!      sss = 0.d0
!      do j=1,npe
!        sss = max(sss,abs(aa(i,j)))
!        write(90,'(2i6,d15.7)') i,j,aa(i,j)
!      enddo
!      write(90,*)
!      write(90,*)
!      write(91,'(i6,2d15.6)') i,yc(i),bb(i)
!      enddo
!      write(91,*)
!      write(91,*)
!      endif
!
! ----------------- SOLUZIONE CON ROUTINE GAUSS ------------------------
!      do i=1,npe
!      aa(i,npe+1)=bb(i)
!      end do
!      call gauss(aa,npamx,npe,npe+1,1.d-12)
!      do i=1,npe
!      bb(i)=aa(i,npe+1)
!      end do 



!
!     -------------   SOLUZIONE SISTEMA LINEARE CON ESSL   -------------

!      call dgefcd(aa,npamx,npe,indx,1,rcond,det,appo,npamx)
!      call dges(aa,npamx,npe,indx,bb,0) 

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
!
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
      do i=mb+1,mb+mf
        dphn(i) = bb(i)
        dphnb(i) = bb(i)
!        write(77,'(i6,4d15.7)') i,yc(i),zc(i),dphn(i),ph(i)
      enddo
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
!cc      write(88,*) '# jt ',jt
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
!cc        write(88,'(i6,5d15.7)') i,a(i),b(i),c(i),d(i),e(i)
      enddo
!cc      write(88,*)
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
!cc        write(88,'(i6,5d15.7)') i,a(i),b(i),c(i),d(i),e(i)
      enddo
!cc      write(88,*)
!cc      write(88,*) 
      
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

!cc* -  scrivo un po di controlli
!cc      do i=1,m
!cc        if(i.eq.1)then
!cc         ii1 = mb
!cc         jj  = mt+1
!cc        else
!cc         ii1 = mt+i-1
!cc         jj  = nt-i+2
!cc        endif
!cc        ii = mt+i
!cc        dxi  = xi(ii1)-xis(ii)
!cc        dze  = ze(ii1)-zes(ii)
!cc        phig = a(ii)+b(ii)*dxi+c(ii)*dze+0.5d0*d(ii)*(dxi**2-dze**2)+
!cc     #         e(ii)*dxi*dze 
!ccc        write(92,'(4d15.7)') yc(ii1),zc(ii1),phb(ii1),phig 
!cc        dxi  = xi(jj)-xis(ii)
!cc        dze  = ze(jj)-zes(ii)
!cc        phig = a(ii)+b(ii)*dxi+c(ii)*dze+0.5d0*d(ii)*(dxi**2-dze**2)+
!cc     #         e(ii)*dxi*dze 
!ccc        write(93,'(4d15.7)') yc(jj ),zc(jj ),ph(jj ),phig 
!cc       enddo         
!cc       do i=1,m+n
!cc        if(i.eq.1.and.m.eq.0)then
!cc         ii1 = mb
!cc         jj  = mb+1
!cc         jj1 = ntt
!cc         ii  = nt+i
!cc         i1  = i-1
!cc        elseif(i.eq.1.and.m.ne.0)then
!cc         ii1 = mb
!cc         jj  = mb+1
!cc         jj1 = nt 
!cc         ii  = mt+i
!cc         i1  = i-1
!cc        elseif(i.ge.2.and.i.le.m)then
!cc         ii  = mt+i
!cc         ii1 = ii-1 
!cc         jj  = nt-i+2
!cc         jj1 = jj-1 
!cc         i1  = i-1
!cc        else
!cc         ii  = nt+i-m
!cc         ii1 = ii-1 
!cc         jj  = ntt-(i-m)+2
!cc         jj1 = jj-1 
!cc         i1  = i-1
!cc        endif
!cc        sqh=sqrt(1.d0+hp(ii)**2)
!cc        rlii = rl(ii)
!cc        a1 =  a(ii)
!cc        b1 =  b(ii)
!cc        c1 =  c(ii)
!cc        d1 =  d(ii)
!cc        e1 =  e(ii)
!cc        dxi = xi(ii)-xis(ii)
!cc        dze = ze(ii)-zes(ii)
!cc        phig = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!cc        vxig  = b1+d1*dxi+e1*dze
!cc        vzeg  = c1-d1*dze+e1*dxi
!cc        dxi = xi(jj1)-xis(ii)
!cc        dze = ze(jj1)-zes(ii)
!cc        phif = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!cc        vxif  = b1+d1*dxi+e1*dze
!cc        vzef  = c1-d1*dze+e1*dxi
!cc        vnf   = (vxif*hp(ii)-vzef)/sqh
!cc        vtf   = (-vxif-hp(ii)*vzef)/sqh
!cc        a1 =  a(ii+1)
!cc        b1 =  b(ii+1)
!cc        c1 =  c(ii+1)
!cc        d1 =  d(ii+1)
!cc        e1 =  e(ii+1)
!cc        dxi = xi(ii)-xis(ii+1)
!cc        dze = ze(ii)-zes(ii+1)
!cc        phigg = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!cc        vxigg = b1+d1*dxi+e1*dze
!cc        vzegg = c1-d1*dze+e1*dxi
!cc        dxi = xi(jj1)-xis(ii+1)
!cc        dze = ze(jj1)-zes(ii+1)
!cc        phiff = a1+b1*dxi+c1*dze+0.5d0*d1*(dxi**2-dze**2)+e1*dxi*dze
!cc        vxiff = b1+d1*dxi+e1*dze
!cc        vzeff = c1-d1*dze+e1*dxi
!cc        vnff  = (vxiff*hp(ii)-vzeff)/sqh
!cc        vtff  = (-vxiff-hp(ii)*vzeff)/sqh
!cc        phii = (1.d0-rlii)*phb(ii)+rlii*phig
!cc        dphii= (1.d0-rlii)*dphnb(jj1)+rlii*vnff
!cc        write(62,'(6d15.6)')   yc(ii),zc(ii), phb(ii),phig,phigg,phii
!cc        write(64,'(8d15.7)') 
!cc     #      yc(ii),zc(ii), dphtb(ii),vxig,vxigg,dphn(ii),vzeg,vzegg
!cc        write(63,'(5d15.6)') yc(jj1),zc(jj1),ph(jj1),phif,phiff
!cc        write(65,'(9d15.7)') 
!cc     #      yc(jj1),zc(jj1), dphtb(jj1),vtf,vtff,dphnb(jj1),vnf,vnff,
!cc     #      dphii
!ccc     #      yc(jj1),zc(jj1), dphtb(jj1),vxif,vxiff,dphn(jj1),vzef,vzeff
!cc      enddo
!cc      write(62,*)
!cc      write(62,*)
!cc      write(63,*)
!cc      write(63,*)
!cc      write(64,*)
!cc      write(64,*)
!cc      write(65,*)
!cc      write(65,*)
!cc            
!cc*
      return
      end
