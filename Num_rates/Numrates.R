# Numrates


#    CONSTANTS
#    ---------------------------------------------------
 avog = 6.02214129e23        # Avogadro's constant
 pi = 3.141592653589793
 muc2 = 931.494061         # m_u*c^2 in MeV
 c = 299792458            # speed of light in cm/s
 
cp <- 2                    # charge projectile
ct <- 1                    # charge target
mpu <-   3.01493216        # mass projectile in u
mtu <- 2.01355332           # mass target in u
ns <- 20                    #number of samples from MCMC output
os <- 0                   # output rate samples (0=no/1=yes) for plotting pdf's 
numRates.in <- c(cp,ct,mpu,mtu,ns,os)
mcmcdat <- read.table("numRates.dat",head=T) 
Tgrid <- exp(seq(log(1e-3),log(10),length.out =  100)) 

# 0.0                 ! y_initial [y0]
# 0.001             ! E_min (MeV) [xi] ***
 #1.2                 ! E_max (MeV) [xf] ***
 #1.0d-10             ! accuracy [eps]
 #1.0d-10             ! initial step size (MeV) [h1] ***
# 0.0                 ! minimum step size (MeV) [hmin]
#1.0d-03             ! maximum step size (MeV) [hmax]
 #10.0                ! dx for intermediate output [dxsav]


mue <- M0*M1/(M0+M1)        # reduced mass
#    last factor is for conversion barn --> cm^2 
factor1 = sqrt(8.0/(pi*mue*muc2))*c*avog*1.0e-24

do 300 jj=1, nsamp
do 200 i=1, nt            # i labels temperatures         

kt=0.086173324*t(i)   # kT in MeV
kmax=1000
hsmal=h1
hbig=h1
y=y0

#          integrate rates numerically 
rungekutta(y,xi,xf,eps,h1,hmin,hmax,hsmal,hbig,nok,nbad)
#          integral returned as y
nasv(jj,i)=factor1*y*(kt**(-1.5))

