C     Functions from numerical recepies adjusted to real(8)
C     Adjusted to factorials up to 1000!
C     Mathias Ljungberg, 2006-10-01

      FUNCTION factln(n)
      INTEGER n
      REAL(8) factln
CU    USES gammln
      REAL(8) a(1000),gammln
      SAVE a
      DATA a/1000*-1./
      if (n.lt.0) then
         write(6,*) 'negative factorial in factln'
         stop
      endif
      if (n.le.999) then
        if (a(n+1).lt.0.) a(n+1)=gammln(n+1.d0)
        factln=a(n+1)
      else
        factln=gammln(n+1.d0)
      endif
      return
      END


      FUNCTION factrl(n)
      INTEGER n
      REAL(8) factrl
CU    USES gammln
      INTEGER j,ntop
      REAL(8) a(1000),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      if (n.lt.0) then
        write(6,*) 'negative factorial in factrl'
        stop
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.999) then
        do 11 j=ntop+1,n
          a(j+1)=j*a(j)
11      continue
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(n+1.d0))
      endif
      return
      END


      
      FUNCTION gammln(xx)
      REAL(8) gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END


      FUNCTION gamma(xx)
      REAL(8) gamma,gammln,xx
C     Mathias Ljungberg, 2006-09-23       
      gamma = exp(gammln(xx))
      RETURN
      END


      
