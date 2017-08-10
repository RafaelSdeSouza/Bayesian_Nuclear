CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       PROGRAM numRates
c      numerical calculation of t(d,n) reaction rates;
c      rates are calculated from MCMC resonance samples and
c      numerical S-factor calculated with external code; 
c      if itest=1, rate samples are
c      written to output for select temperatures for plotting
c
c      INPUT
c      numRates.in:  input with info on nuclei and integration
c                    parameters       
c      numRates.dat: resonance samples from MCMC
c 
c      OUTPUT  
c      numRates.out:     output of 16, 50, 84 percentiles for the rates 
c                        on grid of 60 temperatures
c      numRates_smp.out: rate samples at given T;
c                        line 1: selected T values at which rates are 
c                                calculated
c                        line 2,3: logmu and logsigma at those T
c                        lines 4...: rate samples at those T
c
c      !!! Always make sure that the JAGS output is formatted
c          properly, i.e., it can be read by numRates !!!
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     ---------------------------------------------------
c     DECLARATIONS
c     ---------------------------------------------------
      integer i,j,jj,nt,nsamp,nq16,nq50,nq84
      PARAMETER (nt=60)

      real*8 alpha(100000),beta(100000),gamma(100000)
      real*8 xxx(100000)
      real*8 t(nt),nasv(100000,nt)
      real*8 sumxxx,logmu,sumyyy,logsigma,logmu2,logsigma2
      
      character*16 title
      character*132 dummy
            
      integer nbad,nok
      integer kmax,kount
      real*8 y0,y,xi,xf,eps,h1,hmin,hmax,hsmal,hbig
      real*8 dxsav,xp(1000),yp(1000),x,ss,dydx
      real*8 avog,pi,muc2,c,factor1
      real*8 ttest(6),hlogmu(6),hlogsigma(6)
      
      integer z0,z1,itest,ihelp(6)
      real*8 M0,M1,mue,kt

      common/intd/kmax,kount,dxsav,xp,yp

      common/cca/mue,kt,alpha,beta,gamma
      common/cci/z0,z1,jj  
      
      data ttest/0.001d0,0.016d0,0.05d0,0.1d0,1.0d0,10.0d0/  ! T in 10**9 K 
c     ---------------------------------------------------
c     CONSTANTS
c     ---------------------------------------------------
      avog=6.02214129d23        ! Avogadro's constant
      pi=3.141592653589793d0
      muc2=931.494061d0         ! m_u*c^2 in MeV
      c=299792458d2             ! speed of light in cm/s
c     ---------------------------------------------------
c     INPUT
c     ---------------------------------------------------
      open(unit=9,file='numRates.in',status='unknown')
      read(9,'(a16)') title
      read(9,*) z0,z1           ! charge projectile, target
      read(9,*) M0              ! mass projectile in u
      read(9,*) M1              ! mass target in u
      read(9,*) nsamp           ! number of samples from MCMC output
      read(9,*) itest           ! output rate samples (0=no/1=yes) 
      read(9,*) y0              ! y-initial
      read(9,*) xi              ! E_min (MeV)
      read(9,*) xf              ! E_max (MeV)
      read(9,*) eps             ! accuracy
      read(9,*) h1              ! initial step size (MeV)
      read(9,*) hmin            ! minimum step size (MeV)
      read(9,*) hmax            ! maximum step size (MeV)
      read(9,*) dxsav           ! dx for intermediate output
      close(9)

      write(6,*) 'Always check if JAGS output has right format!'
      write(6,*) 'Print first sample and compare with input...'
      write(6,*)
      
c     read MCMC samples
      open(unit=8,file='numRates.dat',status='unknown')
      
c      read(8,'(a100)') dummy
      do 432 j=1,nsamp          ! j labels samples
         read(8,*,end=2222) alpha(j),beta(j),gamma(j)      
 432  continue
 2222 continue
c      write(6,*) j-1,nsamp
      if(j-1.ne.nsamp)then
         write(6,*) 'you are trying to read too many samples'
         stop
      endif
c!!!!!!!!!!!!!!!
c      write first sample to screen to check if input was ok
      write(6,*) ' first sample: alpha, beta, gamma'
      write(6,*) alpha(1)
      write(6,*) beta(1)
      write(6,*) gamma(1)
c!!!!!!!!!!!!!!!
      close (8)

 9999 format(9x,f12.6,2x,f11.7,2x,1pe13.6)      
c     ---------------------------------------------------
c     GRID OF TEMPERATURES
c     ---------------------------------------------------
      t(1)=0.001d0
      t(2)=0.002d0
      t(3)=0.003d0
      t(4)=0.004d0
      t(5)=0.005d0
      t(6)=0.006d0
      t(7)=0.007d0
      t(8)=0.008d0
      t(9)=0.009d0
      t(10)=0.010d0
      t(11)=0.011d0
      t(12)=0.012d0
      t(13)=0.013d0
      t(14)=0.014d0
      t(15)=0.015d0
      t(16)=0.016d0
      t(17)=0.018d0
      t(18)=0.020d0
      t(19)=0.025d0
      t(20)=0.030d0
      t(21)=0.040d0
      t(22)=0.050d0
      t(23)=0.060d0
      t(24)=0.070d0
      t(25)=0.080d0
      t(26)=0.09d0
      t(27)=0.10d0
      t(28)=0.11d0
      t(29)=0.12d0
      t(30)=0.13d0
      t(31)=0.14d0
      t(32)=0.15d0
      t(33)=0.16d0
      t(34)=0.18d0
      t(35)=0.20d0
      t(36)=0.25d0
      t(37)=0.30d0
      t(38)=0.35d0
      t(39)=0.40d0
      t(40)=0.45d0
      t(41)=0.50d0
      t(42)=0.60d0
      t(43)=0.70d0
      t(44)=0.80d0
      t(45)=0.90d0
      t(46)=1.0d0
      t(47)=1.25d0
      t(48)=1.5d0
      t(49)=1.75d0
      t(50)=2.0d0
      t(51)=2.5d0
      t(52)=3.0d0
      t(53)=3.5d0
      t(54)=4.0d0
      t(55)=5.0d0
      t(56)=6.0d0
      t(57)=7.0d0
      t(58)=8.0d0
      t(59)=9.0d0
      t(60)=10.0d0
c     ---------------------------------------------------
c     CALCULATION OF MC RATES
c     ---------------------------------------------------
      mue=M0*M1/(M0+M1)         ! reduced mass
c     last factor is for conversion barn --> cm^2 
      factor1=dsqrt(8.0d0/(pi*mue*muc2))*c*avog*1.0d-24

      do 300 jj=1,nsamp
         do 200 i=1,nt             ! i labels temperatures         
            kt=0.086173324*t(i)    ! kT in MeV
            kmax=1000
            hsmal=h1
            hbig=h1
            y=y0

c           integrate rates numerically 
            call rungekutta(y,xi,xf,eps,h1,hmin,hmax,hsmal,
     &         hbig,nok,nbad)
c           integral returned as y
            nasv(jj,i)=factor1*y*(kt**(-1.5d0))
 200     continue
         write(6,*) jj
 300  continue 
 
c     ---------------------------------------------------
c     OUTPUT MC RATES AND LOGNORMAL PARAMETERS
c     ---------------------------------------------------
      open(unit=10,file='numRates.out',status='unknown')
c      open(unit=15,file='fu.out',status='unknown')

      write(10,'(a16)') title
      write(10,*) '       charge projectile=',z0
      write(10,*) '           charge target=',z1
      write(10,*) ' nuclear mass projectile=',M0
      write(10,*) '     nuclear mass target=',M1
      write(10,*) '       number of samples=',nsamp
      write(10,*)
      write(10,8334) 

c     project 2D into 1D array of samples xxx
c     i: temperatures; j: samples
      do 650 i=1,nt             ! loop over temperatures
         do 651 j=1,nsamp
            xxx(j)=nasv(j,i)
 651     continue
c        sort rate samples into ascending order
         call sort(nsamp,xxx)

         nq50=int(0.50d0*nsamp)
         nq16=int(0.16d0*nsamp)
         nq84=int(0.84d0*nsamp)

c        calculate lognormal parameters mu and sigma directly from
c        samples [mu and sigma of Gaussian for ln(x)]
         sumxxx=0.0d0
         do 751 j=1,nsamp
            sumxxx=sumxxx+dlog(xxx(j))
 751     continue
         logmu=sumxxx/nsamp

         sumyyy=0.0d0
         do 752 j=1,nsamp
            sumyyy=sumyyy+(dlog(xxx(j))-logmu)**2.0d0
 752     continue
         logsigma=dsqrt(sumyyy/nsamp)

c         fu=dexp(logsigma)
c         write(10,'(f5.3)') fu
         
c        for tests, output lognormal parameters directly from rates;
c        see Longland et al. 2010 paper
ccc         logmu2=dlog(xxx(nq50))
ccc         logsigma2=dlog(dsqrt(xxx(nq84)/xxx(nq16)))

         write(10,8333) t(i),xxx(nq16),xxx(nq50),xxx(nq84),
     &       logmu,logsigma,dexp(logsigma)
c         write(15,'(1f6.3,3x,f5.3)') t(i),fu
     
           if(itest.eq.1)then
              do 735 k=1,6
                 if(t(i).eq.ttest(k))then
                    hlogmu(k)=logmu
                    hlogsigma(k)=logsigma
                    ihelp(k)=i
c                    write(6,*) ihelp(k)
                 endif
 735          continue 
           endif
 
 650  continue      ! next temperature
      close(10)
c      close(15)

c     ---------------------------------------------------
c     OUTPUT RATE SAMPLES AT SELECT TEMPERATURES (itest=1)
c     ---------------------------------------------------
      if(itest.eq.1)then
      
      open(unit=13,file='numRates_samp.out',status='unknown')

      write(13,8335) (ttest(k),k=1,6)
      write(13,8336) (hlogmu(k),k=1,6)
      write(13,8336) (hlogsigma(k),k=1,6)

c     ihelp: integer labeling the temperature (i=1...60)
      do 835 j=1,nsamp
         write(13,8337) (nasv(j,ihelp(k)),k=1,6)      
 835  continue
      
      close(13) 
           
      endif
      
 8333 format(1f6.3,3(3x,1pe9.3),3x,1pe10.3,3x,1pe10.3,3x,1pe10.3) 

 8334 format(1x,'T9',8x,'low',9x,'med',9x,'high',8x,'logmu',7x,
     &   'logsigma',7x,'f.u.')
 8335 format(6(1f6.3,3x))
 8336 format(6(1pe12.5,3x))
 8337 format(6(1pe9.3,3x))
c     ---------------------------------------------------
c     OUTPUT FIRST S-FACTOR SAMPLE TO SCREEN
c     ---------------------------------------------------
c      write(6,*) 
c      write(6,*) ' first rate samples'
c      do 31 i=1,nt
c         write(6,8888) t(i),nasv(1,i) 
c 31   continue
c 8888 format(1f6.3,3x,1pe9.3)
c     ---------------------------------------------------
      STOP   
      END

c     ******************************************************
c     SUBROUTINES AND FUNCTIONS
c     ******************************************************
      subroutine integrateme(x,y,dydx)
c     ---------------------------------------------------
c     function to be integrated: dydx 
c     ******************************************************
c     full integral in the following short test should be  
C     equal to 1.0 for the Gaussian
c      real*8 x,y,dydx
c      real*8 c,w,pie
c      c=100.0
c      w=0.001
c      pie=4.0d+00*datan(1.0d+00)
c      dydx=1.0/w/dsqrt(2.0*pie)*dexp(-(x-c)**2/(2.0*w**2))
c      return
c      end
c     ******************************************************

      integer z0,z1,jj
      real*8 x,y,dydx,sfac
      real*8 twopieta,mue,kt
      real*8 alpha(100000),beta(100000),gamma(100000)

      common/cca/mue,kt,alpha,beta,gamma
      common/cci/z0,z1,jj      

c     cm energy: x  
      twopieta=0.98951013d0*z0*z1*dsqrt(mue/x)    ! 2*pi*eta

      call tdn(alpha(jj),beta(jj),gamma(jj),x,sfac)
c      write(6,*) sfac
      dydx=dexp(-twopieta)*dexp(-x/kt)*sfac

      return
      end
c     ******************************************************
      subroutine rungekutta(ystart,x1,x2,eps,h1,hmin,hmax,hsmal,
     &         hbig,nok,nbad)
c     ******************************************************
      implicit none 
      integer nbad,nok,maxstp
      real*8 eps,h1,hmin,hmax,hsmal,hbig,x1,x2,ystart,tiny
      parameter (maxstp=1000000,tiny=1.0d-60)
      integer kmax,kount,nstp
      real*8 dxsav,h,hdid,hnext,x,xsav,dydx
      real*8 xp(1000),y,yp(1000),yscal

      common/intd/kmax,kount,dxsav,xp,yp
        
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      y=ystart
      if (kmax.gt.0) xsav=x-2.*dxsav
c
c     maxstp is the maximum number of steps, go at most that far
      do 16 nstp=1,maxstp
         call integrateme(x,y,dydx)
         yscal=dabs(y)+dabs(h*dydx)+tiny
         if (kmax.gt.0) then
            if (dabs(x-xsav).gt.dabs(dxsav)) then
               if (kount.lt.kmax-1)then
                  kount=kount+1
                  xp(kount)=x
                  yp(kount)=y
                  xsav=x
               endif
            endif
         endif
         if ((x+h-x2)*(x+h-x1).gt.0.0) h=x2-x
         call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext)
         if (hdid.eq.h) then
            nok=nok+1
         else
            nbad=nbad+1
         endif
         if (hdid.lt.hsmal) hsmal=hdid
         if (hdid.gt.hbig)  hbig=hdid
c
c        kick out of the loop if integration completes before
c        maxsteps encountered
         if ((x-x2)*(x2-x1).ge.0.0) then
            ystart=y
            if (kmax.ne.0) then
               kount=kount+1
               xp(kount)=x
               yp(kount)=y
            endif
            return
         endif
         if (dabs(hnext).lt.hmin) then
            write(6,*) 'stepsize smaller than minimum'
            stop
         endif
         h=hnext 
         if (h.gt.hmax) h=hmax
 16   enddo
  9   format(2x,f12.3,2x,d12.5)
      write(6,*) 'too many steps'
      return
      end
c     ******************************************************
      subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext)
c     ---------------------------------------------------
c     rkqs does the adaptive stepsize control
c     straight from numerical recipes
c     ******************************************************
      implicit none
      real*8 eps,hdid,hnext,htry,x,dydx,y,yscal
      real*8 errmax,h,htemp,xnew,yerr,ytemp,safety,pgrow,pshrnk,errcon
      parameter (safety=.9,pgrow=-.2,pshrnk=-.25,errcon=1.89e-4)
      h=htry
c 
c     take a baby step and test the water
 1    call rkck(y,dydx,x,h,ytemp,yerr)
c
c     set integration scale, based on result
      errmax=0.
      errmax=max(errmax,dabs(yerr/yscal))
      errmax=errmax/eps   
c
c     evaluate integration stepsize to see if its kosher
      if (errmax.gt.1.) then
         htemp=safety*h*(errmax**pshrnk)
         h=sign(max(dabs(htemp),0.1*dabs(h)),h)
         xnew=x+h
         if (xnew.eq.x)then
            write(6,*) 'stepsize underflow in rkqs'
            stop
         endif
         goto 1
      else
         if (errmax.gt.errcon) then
           hnext=safety*h*(errmax**pgrow)
         else
           hnext=5.*h
         endif
         hdid=h
         x=x+h
         y=ytemp   
         return
      endif
      end
c     ******************************************************
      subroutine rkck(y,dydx,x,h,yout,yerr)
c     ---------------------------------------------------
c     The rkck routine takes 1 Runge-Kutta step
c     see Abramowitz and Stegun Sec. 25.5 (p. 896)
c     ******************************************************
      implicit none
      real*8 h,x,dydx,y,yerr,yout
      real*8 ak2,ak3,ak4,ak5,ak6,ytemp   
      real*8 a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,
     &   b52,b53,b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,
     &   dc1,dc3,dc4,dc5,dc6
      parameter (a2=.2,a3=.3,a4=.6,a5=1.,a6=.875,b21=.2,b31=3./40.,
     &    b32=9./40.,b41=0.3,b42=-.9,b43=1.2,b51=-11./54.,b52=2.5,
     &    b53=-70./27.,b54=35./27.,b61=1631./55296.,b62=175./512.,
     &    b63=575./13824.,b64=44275./110592.,b65=253./4096.,
     &    c1=37./378.,c3=250./621.,c4=125./594.,c6=512./1771.,
     &    dc1=c1-2825./27648.,dc3=c3-18575./48384,
     &    dc4=c4-13525./55296.,dc5=-277./14336.,dc6=c6-0.25)
      ytemp=y+b21*h*dydx
      call integrateme(x+a2*h,ytemp,ak2) 
      ytemp=y+h*(b31*dydx+b32*ak2)
      call integrateme(x+a3*h,ytemp,ak3)
      ytemp=y+h*(b41*dydx+b42*ak2+b43*ak3)
      call integrateme(x+a4*h,ytemp,ak4)
      ytemp=y+h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
      call integrateme(x+a5*h,ytemp,ak5)
      ytemp=y+h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)
      call integrateme(x+a6*h,ytemp,ak6)
      yout=y+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6)
      yerr=h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)
      return
      end
c     ******************************************************
      subroutine sort(n,arr)
c     ******************************************************
      integer n,m,nstack
      real*8 arr(n)
      parameter (m=7,nstack=50)
c        Sorts an array arr(1:n) into ascending numerical order using the
c        Quicksort algorithm. n is input; arr is replaced on output by its
c        sorted rearrangement.
c        Parameters: M is the size of subarrays sorted by straight insertion
c        and NSTACK is the required auxiliary storage.        
      integer i,ir,j,jstack,k,l,istack(nstack)
      real a,temp
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.m)then
         do 12 j=l+1,ir
	        a=arr(j)
            do 11 i=j-1,1,-1
	           if(arr(i).le.a)goto 2
	           arr(i+1)=arr(i)
 11	        continue
	        i=l-1
  2	        arr(i+1)=a
 12	     continue

	     if(jstack.eq.0)return
	     ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
          k=(l+ir)/2
	      temp=arr(k)
	      arr(k)=arr(l+1)
	      arr(l+1)=temp
	      if(arr(l).gt.arr(ir))then
	         temp=arr(l)
	         arr(l)=arr(ir)
	         arr(ir)=temp
	      endif
	      if(arr(l+1).gt.arr(ir))then
	         temp=arr(l+1)
	         arr(l+1)=arr(ir)
	         arr(ir)=temp
	      endif
	      if(arr(l).gt.arr(l+1))then
	         temp=arr(l)
	         arr(l)=arr(l+1)
	         arr(l+1)=temp
	      endif
	      i=l+1
	      j=ir
	      a=arr(l+1)
 3	      continue
             i=i+1
	      if(arr(i).lt.a)goto 3
 4        continue
             j=j-1
	      if(arr(j).gt.a)goto 4
	      if(j.lt.i)goto 5
	      temp=arr(i)
	      arr(i)=arr(j)
	      arr(j)=temp
	      goto 3
 5        arr(l+1)=arr(j)
          arr(j)=a
	      jstack=jstack+2
	      if(jstack.gt.nstack)then
             write(6,*) 'nstack too small in sort2'
             stop
          endif
	      if(ir-i+1.ge.j-l)then
	         istack(jstack)=ir
	         istack(jstack-1)=i
	         ir=j-1
	      else
	         istack(jstack)=j-1
	         istack(jstack-1)=l
	         l=i
	      endif
      endif
      goto 1
      end
c     ******************************************************      
c     ******************************************************      
      subroutine tdn(e1,gi,gf,ecm,sfac)
*-------------------------------------------------------------------*
* Calcule sig(E) et S(E) pour tdn comme dans Barker PRC (1997) 2646 *
* se superpose exactement au S(E) de DAACV04                        *
*-------------------------------------------------------------------*
      implicit real*8(a-h,o-z)
      parameter (cs_pi=0.31415926535897932384626433832795028d+01)
      parameter (cs_2pi=2.d0*cs_pi)
      data hm/20.9/
      data ai,af,li,lf/6,5,0,2/ ! Bar97
! donnees H3(d,n)
      data a1i,a2i,z1i,z2i/3,2,1,1/,a1f,a2f,z1f,z2f/4,1,2,0/
      data q/17.59/ ! Bar97

      pi=acos(-1.0d0)
      rmui=a1i*a2i/(a1i+a2i)
      rmuf=a1f*a2f/(a1f+a2f)
      hmi=hm/rmui
      hmf=hm/rmuf
      eta0i=1.44*z1i*z2i/(2*hmi)
      eta0f=1.44*z1f*z2f/(2*hmf)

      xki=sqrt(e1/hmi)
      xkf=sqrt((e1+q)/hmf)
      call coul(li,xki*ai,eta0i/xki,pli,bi)
      call coul(lf,xkf*af,eta0f/xkf,plf,bf)
      
      xki=sqrt(ecm/hmi)
      xkf=sqrt((ecm+q)/hmf)
      etai=eta0i/xki
      etaf=eta0f/xkf
      call coul(li,xki*ai,etai,pli,si)
      call coul(lf,xkf*af,etaf,plf,sf)
      gitot=2*gi*pli
      gftot=2*gf*plf
      del=-gi*(si-bi)-gf*(sf-bf)
      sdp2=gitot*gftot/((e1+del-ecm)**2+((gitot+gftot)/2)**2)
      sig=2*pi*sdp2/(3*xki**2)
      
      sigma = sig / 100
      
      sfac = sigma * exp (cs_2pi*etai) * ecm  

      end
c     ******************************************************      
      subroutine coul(l,ak,eta,pl,sl)
c     ******************************************************      
      implicit real*8(a-h,o-z)
      dimension fc(100),dfc(100),gc(100),dgc(100)
c     ak=rho, eta=eta, l=l,...
      call coufra(ak,eta,l,l,fc,dfc,gc,dgc)

      pl=ak/(fc(l+1)**2+gc(l+1)**2)
      sl=pl*(fc(l+1)*dfc(l+1)+gc(l+1)*dgc(l+1))
      return
      end
 
c     ******************************************************      
      SUBROUTINE COUFRA(RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP)
C*** FONCTIONS COULOMBIENNES CALCULEES EN R = RHO PAR LA METHODE DES FRA
C*** CONTINUES DE STEED. MINL ET MAXL CORRESPONDENT AUX VRAIES VALEURS D
C*** VOIR BARNETT, FENG, STEED ET GOLDFARB, CPC 1974 *******************
      IMPLICIT double precision (A-H,O-Z) 
      double precision K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(10),FCP(10),GC(10),GCP(10)
      DATA ACCUR,STEP/1.D-7,100.0D0/
      PACE = STEP 
      ACC = ACCUR 
      R = RHO 
      KTR = 1 
      LMAX = MAXL 
      LMIN1 = MINL+1  
      XLL1 = MINL*LMIN1 
      ETA2 = ETA**2 
      TURN = ETA+SQRT(ETA2+XLL1)
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.D-6) KTR = -1  
      KTRP = KTR
      GO TO 2 
    1 R = TURN
      TF = F
      TFP = FP
      LMAX = MINL 
      KTRP = 1
    2 ETAR = ETA*R
   21 RHO2=R*R
      PL = LMAX+1 
      PMX = PL+0.5D0  
C** FRACTION CONTINUE POUR FP(MAXL)/F(MAXL) ; XL=F ; XLPRIME=FP ********
      FP = ETA/PL+PL/R
      DK = ETAR+ETAR  
      DEL = 0
      D = 0
      F = 1 
      K = (PL*PL-PL+ETAR)*(PL+PL-1) 
      IF(PL*PL+PL+ETAR.NE.0.) GO TO 3 
      R = R*1.0000001D0 
      GO TO 2 
    3 H = (PL*PL+ETA2)*(1-PL*PL)*RHO2 
      K = K+DK+PL*PL*6
      D = 1/(D*H+K) 
      DEL = DEL*(D*K-1) 
      IF(PL.LT.PMX) DEL = -R*(PL*PL+ETA2)*(PL+1)*D/PL 
      PL = PL+1 
      FP = FP+DEL 
      IF(D.LT.0) F = -F
      IF(PL.GT.20000.0D0) GO TO 11
      IF(ABS(DEL/FP).GE.ACC) GO TO 3
      FP = F*FP 
      IF(LMAX.EQ.MINL) GO TO 5  
      FC(LMAX+1) = F  
      FCP(LMAX+1) = FP
C*** RECURRENCE ARRIERE POUR F ET FP ; GC,GCP UTILISES POUR STOCKAGE ***
      L = LMAX
      DO 4 LP=LMIN1,LMAX
      PL = L
      GC(L+1) = ETA/PL+PL/R 
      GCP(L+1) = SQRT(ETA2+PL*PL)/PL
      FC(L) =(GC(L+1)*FC(L+1)+FCP(L+1))/GCP(L+1)
      FCP(L) = GC(L+1)*FC(L)-GCP(L+1)*FC(L+1) 
    4 L = L-1 
      F = FC(LMIN1) 
      FP = FCP(LMIN1) 
    5 IF(KTRP.EQ.-1) GO TO 1
C*** MEME CALCUL POUR R = TURN SI RHO.LT.TURN 
C*** P + I.Q CALCULE EN MINL , EQUATION (32)
      P = 0
      Q = R-ETA 
      PL = 0 
      AR = -(ETA2+XLL1) 
      AI = ETA
      BR = Q+Q
      BI = 2
      WI = ETA+ETA
      DR = BR/(BR*BR+BI*BI) 
      DI = -BI/(BR*BR+BI*BI)
      DP = -(AR*DI+AI*DR) 
      DQ = AR*DR-AI*DI
    6 P = P+DP
      Q = Q+DQ
      PL = PL+2 
      AR = AR+PL
      AI = AI+WI
      BI = BI+2 
      D = AR*DR-AI*DI+BR
      DI = AI*DR+AR*DI+BI 
      T = 1/(D*D+DI*DI) 
      DR = T*D
      DI = -T*DI
      H = BR*DR-BI*DI-1
      K = BI*DR+BR*DI 
      T = DP*H-DQ*K 
      DQ = DP*K+DQ*H  
      DP = T
      IF(PL.GT.46000.0D0) GO TO 11
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6
      P = P/R 
      Q = Q/R 
C*** CALCUL DE FP,G,GP, ET NORMALISATION DE F EN L = MINL **************
      G = (FP-P*F)/Q  
      GP = P*G-Q*F
      W = 1/SQRT(FP*G-F*GP) 
      G = W*G 
      GP = W*GP 
      IF(KTR.EQ.1) GO TO 8
      F = TF
      FP = TFP
      LMAX = MAXL 
C*** CALCUL DE G(MINL) ET GP(MINL) PAR INTEGRATION RUNGE-KUTTA A PARTIR 
C***         VOIR FOX ET MAYERS(1968) PG 202
      IF(RHO.LT.0.2D0*TURN) PACE = 999.0D0
      R3=1.0D0/3.0D0  
      H = (RHO-TURN)/(PACE+1) 
      H2 = H/2
      I2 = INT(PACE+0.001D0)
      ETAH = ETA*H
      H2LL = H2*XLL1  
      S = (ETAH+H2LL/R)/R-H2
    7 RH2 = R+H2
      T = (ETAH+H2LL/RH2)/RH2-H2
      K1 = H2*GP
      M1 = S*G
      K2 = H2*(GP+M1) 
      M2 = T*(G+K1) 
      K3 = H*(GP+M2)  
      M3 = T*(G+K2) 
      M3 = M3+M3
      K4 = H2*(GP+M3) 
      RH = R+H
      S = (ETAH+H2LL/RH)/RH-H2  
      M4 = S*(G+K3) 
      G = G+(K1+K2+K2+K3+K4)*R3 
      GP = GP+(M1+M2+M2+M3+M4)*R3 
      R = RH
      I2 = I2-1 
      IF(ABS(GP).GT.1.D35) GO TO 11
      IF(I2.GE.0) GO TO 7 
      W = 1/(FP*G-F*GP) 
C*** RECURRENCE AVANT A PARTIR DE GC(MINL) ET GCP(MINL) 
C*** RENORMALISATION DE FC ET FCP POUR CHAQUE VALEUR DE L **************
    8 GC(LMIN1) = G 
      GCP(LMIN1) = GP 
      IF(LMAX.EQ.MINL) GO TO 10 
      DO 9 L=LMIN1,LMAX 
      T = GC(L+1) 
      GC(L+1) = (GC(L)*GC(L+1)-GCP(L))/GCP(L+1) 
      GCP(L+1) = GC(L)*GCP(L+1)-GC(L+1)*T 
      FC(L+1) = W*FC(L+1) 
    9 FCP(L+1) = W*FCP(L+1) 
      FC(LMIN1) = W*FC(LMIN1) 
      FCP(LMIN1) = W*FCP(LMIN1) 
      RETURN
   10 FC(LMIN1) = W*F 
      FCP(LMIN1) = W*FP 
      RETURN
   11 W = 0
      G = 0
      GP = 0 
      GO TO 8 
      END 

    
      
      
      
