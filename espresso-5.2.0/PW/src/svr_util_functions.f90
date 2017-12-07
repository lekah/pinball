function sumnoises(nn)
   !
   USE kinds,     ONLY : DP
   !
   implicit none
   integer, intent(in) :: nn
   !
   ! returns the sum of n independent gaussian noises squared
   ! (i.e. equivalent to summing the square of the return values of nn calls to gasdev) 
   !
   real(DP) sumnoises
   real*8, external :: gamdev,gasdev
   !
   if(nn==0) then
     sumnoises=0.0
   elseif(nn==1) then
     sumnoises=gasdev()**2
   elseif(modulo(nn,2)==0) then
     sumnoises=2.0*gamdev(nn/2)
   else
     sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
   end if
   !
endfunction sumnoises
!
function gamdev (ia)
   !
   ! gamma-distributed random number, implemented as described in numerical recipes
   !  
   USE kinds,     ONLY : DP
   !
   implicit none
   integer, intent(in) :: ia
   real(DP) gamdev
   integer j
   real*8 am,e,s,v1,v2,x,y
   real*8, external :: ran1
   !
   if(ia.lt.1) pause 'bad argument in gamdev'
     !
   if(ia.lt.6)then
     !
     x=1.0
     !
     do 11 j=1,ia
      !
      x=x*ran1()
      !
  11  continue
      !
      x=-log(x)
   else
      !  
  22  v1=2.*ran1()-1.0
      v2=2.*ran1()-1.0
      !
      if(v1**2+v2**2.gt.1.) goto 22
      !
      y=v2/v1
      am=ia-1
      s=sqrt(2.*am+1.)
      x=s*y+am
      !
      if(x.le.0.) goto 22
      !
      e=(1.+y**2)*exp(am*log(x/am)-s*y)
      ! 
      if(ran1().gt.e) goto 22
   endif
   ! 
   gamdev=x
   !
endfunction gamdev
   !
   FUNCTION gasdev()
   !
   ! gaussian-distributed random number, implemented as described in numerical recipes
   ! 
   USE kinds,     ONLY : DP
   implicit none
    REAL(DP) gasdev
   integer, save :: iset = 0
   real*8, save :: gset
   real*8, external :: ran1
   real*8 fac,rsq,v1,v2
   if(iset==0) then
 33   v1=2.*ran1()-1.0d0
      v2=2.*ran1()-1.0d0
      rsq=v1**2+v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 33
      fac=sqrt(-2.*log(rsq)/rsq)
      gset=v1*fac
      gasdev=v2*fac
      iset=1
  else
      gasdev=gset
      iset=0
  end if
  END FUNCTION gasdev
  ! 
  FUNCTION ran1()
  !
  ! random number generator
  !
   USE kinds,     ONLY : DP
  INTEGER IA,IM,IQ,IR,NTAB,NDIV
  REAL ran1,AM,EPS,RNMX
  PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
            NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER j,k,iv(NTAB),iy
  SAVE iv,iy
  DATA iv /NTAB*0/, iy /0/
  INTEGER, SAVE :: idum=0 !! ATTENTION: THE SEED IS HARDCODED
  if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do 11 j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
   endif
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
   if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
  END FUNCTION ran1
