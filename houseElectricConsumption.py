#
# houseElectricConsumption.py
# python version of the matlab code by Jan with some slight modifications to 
# introduce interhouse variability rather than randomness in the time series of 
# each house
#
# generates power consumption profile roughly
#  matching the data reported in 
#
#  S.K. Firth, K.J. Lomas, A.J. Wright and R. Wall: "Identifying trends 
#  in the use of domestic appliances from household electricity
#  consumption measurements", Energy and Buildings vol. 40(5), pp.
#  926-936, 2008 .
#
#  The basic profile is matched to the total power consumption values 
#  recorded at one house over a np.single day at five-minute intervals, 
#  provided as an example in the paper. 
#  
#  Usage:
#
#  houseElectricConsumption(), by itself, returns the basic load profile
#  measured in Watts, with a sampling period of 5 min.
#
#  houseElectricConsumption(Pvar) adds normal distibuted random 
#  zero-mean noise with variance Pvar (units of Watt) to the basic load 
#  profile Pnom. 
#  The noise is truncated such that the consumption is never below 
#  0.5*Pnom[t] at any sample point t in order to model consumption
#  from continuously operating appliances and standby appliances such
#  as computers, electric clocks, TVs etc. Be aware that specifying a 
#  large variance generally tends to increase the average power 
#  consumption throughout the day.
#
#  houseElectricConsumption(Pvar, PP, Ton, Toff) adds the effect of 
#  air conditioning, large refridgerators and similar to the base load 
#  profile (so-called "cold appliances" in Firth et. al.). 
#  PP is the peak power amplitude (Watt), while Ton is the on-time 
#  (samples) and Toff is the off-time (samples) of the equipment. 
#  Set Pvar = 0 to recover the load profile without noise.
#
#  houseElectricConsumption(Pvar, PP, Ton, Toff, Ts) returns a noisy 
#  load profile sampled at the sample rate Ts (given in minutes). 
#  Set Pvar = 0 to recover the load profile without noise, sampled at 
#  period Ts.
#  Set PP = 0 to recover the load profile (with or without noise) without
#  "cold appliances", sampled at time Ts.
#
#  For each of the variations above, 
#  [P, T] = houseElectricConsumption( ... ) returns a time vector T 
#  (measured in minutes from 00:00) along with the power consumption 
#  profile P.
#

def houseElectricConsumption():

    pvar = 0; 
    pp = 0;
    ton = 6;
    toff = 8;
    ts = 5;
    # Check arguments
    #if nargin > 0, pvar = Pvar; end
    #if nargin > 3,  
        #pp = PP; 
        #ton = Ton; 
        #toff = Toff;
    #end
    #if nargin > 4
        #ts = Ts;
    #end
    import numpy as np
    import scipy.signal as scisig 
    ##
    np.random.seed(4)
    multFac = 0.005
    dayFac = 0.001
    aftFac =0.002 
    eveningFac = 0.002 
    # Base load profile
    #pNight = 500*(1+multFac*np.random.beta(1,3))       # Night-time base consumption
    #pMorning = 2000*(1+multFac*np.random.beta(2,5))    # Morning high load
    #pAfternoon = 1000*(1+aftFac*np.random.beta(.5,.5)) ;  # Afternoon base consumption
    #pCooking = 2000*(1+dayFac*np.random.beta(5,3))    # "Bump" around 5 PM
    #pEvening = 1500*(1+eveningFac*np.random.beta(.5,.5));    # Evening base load
    pNight = 500*(1+multFac*np.random.randn())       # Night-time base consumption
    pMorning = 2000*(1+multFac*np.random.randn())    # Morning high load
    pAfternoon = 1000*(1+aftFac*np.random.randn()) # Afternoon base consumption
    pCooking = 2000*(1+dayFac*np.random.randn()) # "Bump" around 5 PM
    pEvening = 1500*(1+eveningFac*np.random.randn())# Evening base load
    # Time vectors
    tNight = np.arange(1,6*60,ts)
    tMorning =np.arange((6*60+ts),(9*60),ts)
    tMidday = np.arange((9*60+ts),(15*60),ts)
    tEarlyEvening = np.arange((15*60+ts),18*60,ts)
    tEvening = np.arange((18*60+ts),(22*60),ts)
    tLateEvening = np.arange((22*60+ts),(24*60),ts)
    morningProfile = pNight + (pAfternoon-pNight)*np.arange(1,len(tMorning))/len(tMorning) + (pMorning-pAfternoon)*np.sin(np.pi*np.arange(1,len(tMorning))/len(tMorning));
    eveningProfile = pAfternoon +  (pEvening-pAfternoon)*np.arange(1,len(tEarlyEvening))/len(tEarlyEvening) + (pCooking-pEvening)*np.sin(np.pi*np.arange(1,len(tEarlyEvening))/len(tEarlyEvening));
    lateEveningProfile = pEvening + (pNight-pEvening)*np.arange(1,len(tLateEvening))/len(tLateEvening);
    t = np.hstack([tNight, tMorning, tMidday, tEarlyEvening, tEvening, tLateEvening]).T
    pvec = np.vstack([pNight*np.ones((len(tNight),1)),np.reshape(morningProfile.T,(len(morningProfile),1)),pAfternoon*np.ones((len(tMidday),1)),np.reshape(eveningProfile,(len(eveningProfile),1)),pEvening*np.ones((len(tEvening),1)),np.reshape(lateEveningProfile,(len(lateEveningProfile),1))])
    # Filter the profile so that it looks nice and smooth
    # (the filter coefficients were found by hand tuning)
    #print p
    #return p
    #p = 
    powervec = scisig.lfilter([0.2],[1,0.150000], pvec-pNight) + pNight
    plen = len(powervec)
    #import pdb; pdb.set_trace()
    # now we compute the average power consumption each hour
    listOfPlevels = np.split(powervec[0:(plen/24)*24],24)
    avgPower = map(lambda x: np.mean(x),listOfPlevels)
    #indptr =  np.linspace(0,plen,24).astype(np.uint)
    #avgp = map(lambda n1,n2: np.mean(pvec[indptr[n1]:indptr[min([n2,len(pvec)-1])]]),indptr)
            #scisig.convolve([1 -0.7],[1 -0.5]),pvec-pNight) + pNight
    #filter(0.2,conv([1 -0.7],[1 -0.5]),p-pNight) + pNight;
    #print p
    #print avgPower
    return avgPower 



if __name__ == "__main__":
    houseElectricConsumption()

