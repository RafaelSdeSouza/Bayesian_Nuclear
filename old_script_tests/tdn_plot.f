      program tdn_plot
*-------------------------------------------------------------------*
* Calcule sig(E) et S(E) pour tdn comme dans Barker PRC (1997) 2646 *
* se superpose exactement au S(E) de DAACV04                        *
*-------------------------------------------------------------------*
      integer ret
      ret=0
      call tdn_pSub(ret)
      end

c     -----------------------------------------------------------
      subroutine tdn_pSub(ret)

      implicit real*8(a-h,o-z)
      integer ret
      parameter (cs_pi=0.31415926535897932384626433832795028d+01)
      parameter (cs_2pi=2.d0*cs_pi)
      data hm/20.9/
!     "hm" is Pierre's constant (hbar.c)^2/(2u) i.e. (197.3)^2/(2*931) = 20.9
      data ai,af,li,lf/6,5,0,2/ ! Bar97

      data a1i,a2i,z1i,z2i/3,2,1,1/
      data a1f,a2f,z1f,z2f/4,1,2,0/
      data q/17.59/ ! Bar97

      open(unit=11,file='tdn.out',status='unknown')	
      open(unit=12,file='tdn.in',status='unknown')	

      read(12,*) e1
      read(12,*) gi     ! reduced width for deuteron
      read(12,*) gf     ! reduced width for neutron

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
      
!      write(6,2000)bi,bf
 2000 format('bi,bf=',1p,2e14.4)

      do  ie=1,1000
      
      ecm = ie * 0.001
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

      write(11,2001)ecm,sfac
 2001 format(1p,4e14.4)
	      
      end do ! ie
      close(11)
      close(12)

      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine coul(l,ak,eta,pl,sl)
      implicit real*8(a-h,o-z)
      dimension fc(100),dfc(100),gc(100),dgc(100)
      call coufra(ak,eta,l,l,fc,dfc,gc,dgc)
      pl=ak/(fc(l+1)**2+gc(l+1)**2)
      sl=pl*(fc(l+1)*dfc(l+1)+gc(l+1)*dgc(l+1))
      return
      end

*COUFRA 
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

