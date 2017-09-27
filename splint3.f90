
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine splint3(tts,xpo,ypo,zpo,xs,ys,zs,xs2,ys2,zs2,tp,n,nmax)

implicit none

 integer n,nmax
 real*8 tts,xpo,ypo,zpo
 real*8 xs(nmax),ys(nmax),zs(nmax),xs2(nmax),ys2(nmax),zs2(nmax),tp(nmax)

! Local Variables

 integer klo,khi,k
 real*8 h,a,b,a1,b1,c1,uh

! First executable statement

 klo=1
 khi=n
 do while (khi-klo.gt.1)
   k=(khi+klo)/2
   if(tp(k).gt.tts)then
     khi=k
   else
     klo=k
   endif
 enddo 

 h=tp(khi)-tp(klo)
 if (h.eq.0.d0) stop 'bad xa input in splint'
 uh=1.d0/h
 a=(tp(khi)-tts)*uh
 b=(tts-tp(klo))*uh
 a1=(a**3-a)
 b1=(b**3-b)
 c1=(h**2)/6.d0

! coordinates of the interpolated point

 xpo=a*xs(klo)+b*xs(khi)+(a1*xs2(klo)+b1*xs2(khi))*c1

 ypo=a*ys(klo)+b*ys(khi)+(a1*ys2(klo)+b1*ys2(khi))*c1

 zpo=a*zs(klo)+b*zs(khi)+(a1*zs2(klo)+b1*zs2(khi))*c1

return
end
