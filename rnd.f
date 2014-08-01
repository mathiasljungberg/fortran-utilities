
c randomfunktion
      real(8) function rnd(seed)
      implicit none
      integer a,m,q,p,n,ndiv,j,k,seed
      real(8) rm,rmax
c     m=2**31-1 and m=a*q+p
      parameter (a=16807, m=2147483647, rm=1.0d0/m)
      parameter (q=127773, p=2836, n=32, ndiv=1+(m-1)/n)
      parameter (rmax=1.0d0-1.2d-7)
      integer r(n),r0,r1
      save r,r0,r1
      logical first
      data r/n*0/,first/.true./
c     initialize table of random numbers
      if (first) then
         first=.false.
         r1=abs(seed)
         do j=n+8,1,-1
            k=r1/q
            r1=a*(r1-k*q)-p*k
            if (r1.lt.0.) r1=r1+m
            if (j.le.n) r(j)=r1
         enddo
         r0=r(1)
      endif
c     beginning when not initializing
c     compute r1=mod(a*r1,m) without overflows
      k=r1/q
      r1=a*(r1-k*q)-p*k
      if (r1.lt.0) r1=r1+m
      j=1+r0/ndiv
      r0=r(j)
      r(j)=r1
      rnd=min(rm*r0,rmax)
      end
 
     
