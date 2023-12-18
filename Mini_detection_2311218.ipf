Function Mini_detection_230901()

// Setting up variables

variable i,j,k,l,m
variable startblankpoint, endblankpoint,tau,rise, amp, peak, bas, avgamp,EPSCtime,V_levelx, fiftypercent, fibervolleyamp, fibervolleytime, peaktopeak
wave W_coef

make/o/n=100000 minitimes=NAN,miniamps=NAN, wavecounter=NAN, minionset_X=NAN, minionset_Y=NAN			//holds all times and amps. wavecounter for raster plot

variable vtrainstart = 1000 // Start of stimulus within a trace in ms
variable baselinestart = 0.8999
variable baselineend = 0.9999
variable peakwindowstart = (Vtrainstart/1000+0.001)
variable peakwindowend = (Vtrainstart/1000+0.01)


string Num_of_traces = wavelist("0.1*",";","") // Raw traces always were named "0.1 Hz stim_X"

variable width = 0.015
variable length = itemsinlist(Num_of_traces)
variable level_DIF = 0
variable level = 0


make/o/n=100000 minitimes=NAN,miniamps=NAN, wavecounter=NAN		//holds all times and amps. wavecounter for raster plot
make/o/n=(length) Noisewave=NAN
make/o/n=(length) Noisewave_DIF=NAN
make/N=10000/O Wavecounter
make/N=10000/O Timings



// Generate average EPSC waveform

for (i=0;i<length;i+=1)

string currentwave_string = stringfromlist(i,Num_of_traces)

if (i==0)
duplicate/O $currentwave_string, Avg_wave 
Avg_wave = 0
duplicate/O/R=[0,x2pnt(Avg_wave,width)] Avg_wave, Avg_sEPSC 
endif
wave current_wave = $currentwave_string
Avg_wave += current_wave

endfor

Avg_wave /= length


// Find stats for average EPSC


wavestats/q/R=(Vtrainstart/1000+.002,Vtrainstart/1000+.007) Avg_wave         					//find the average EPSC peak
avgamp = V_min
peak = V_minloc
        
FindLevel/q/R=(Vtrainstart/1000 +.001, Vtrainstart/1000+.007)/edge=2 Avg_wave, 0.2*avgamp 		//find the 20% rise to align AR and also to calculate onset
EPSCtime=V_levelx

FindLevel/q/R=(Vtrainstart/1000 +.001, Vtrainstart/1000+.007)/edge=2 Avg_wave, 0.8*avgamp 		//find the 80% to calculate rise
rise=(V_levelx-EPSCtime)*1000

FindLevel/q/R=(Vtrainstart/1000 +.003, Vtrainstart/1000+.02)/edge=1 Avg_wave, 0.5*avgamp 	//find the 50% decay time to report kinetics
fiftypercent=(V_levelX-peak)*1000

FindLevel/q/R=(Vtrainstart/1000 +.003, Vtrainstart/1000+.02)/edge=1 Avg_wave, 0.9*avgamp 	//find the 10% decay time to subtract the EPSC
endblankpoint=X2pnt(Avg_wave,V_levelX)
startblankpoint=X2pnt(Avg_wave,Vtrainstart/1000 - .003)

Differentiate Avg_wave/D=avgwave_DIF
Differentiate Avg_wave/D=avgwave_DIF2

wavestats/q/R=(Vtrainstart/1000 +.001, Vtrainstart/1000+.002) avgwave_DIF2
fibervolleyamp=-Avg_wave(V_maxloc)

wavestats/q/R=(V_maxloc, V_maxloc+.002) Avg_wave
fibervolleyamp+=V_max
fibervolleytime=V_maxloc       																		//find time from the peak of FV to peak of EPSC
peaktopeak=(peak-fibervolleytime)*1000 


wavestats/q/R=(Vtrainstart/1000-0.010,Vtrainstart/1000-0.001) Avg_wave 
avgamp -= (V_avg)

avgamp *= -1


// Mini detetion 


for(i=0;i<length;i+=1)

duplicate/o $stringfromlist(i,Num_of_traces),temp_wave,fit_temp_wave

temp_wave *= 1E12
fit_temp_wave *= 1E12


wavestats/Q/R=(Vtrainstart/1000-0.1,Vtrainstart/1000) /Q temp_wave
 
variable base = V_avg


wavestats/Q/R=(Vtrainstart/1000+0.002,Vtrainstart/1000+0.01) /Q temp_wave 

variable outputX = V_minloc
variable outputY = V_min
variable decayY = (outputY-base)*90/100+base
	
findlevel/Q/R=(OutputX,OutputX+0.02) temp_wave,decayY

variable decayX = V_LevelX


K0 = 0;K1 = (-avgamp/2);K2 = .002; K3 = (-avgamp/2); K4 = .04
CurveFit/Q/L=(V_npnts) dblexp_XOffset temp_wave(decayX,) /D=fit_temp_wave	//fit the decay phase of theEPSC

if (V_flag==1)
CurveFit/Q/L=(V_npnts) exp_XOffset temp_wave(decayX,) /D=fit_temp_wave	//fit the decay phase of theEPSC
V_flag=0
endif

wavestats/q/R=(decayX,) fit_temp_wave	//find the final value of th eEPSC fit (baseline current)

for (l=0;l<=Vtrainstart/1000;l+=deltax(temp_wave))
fit_temp_wave [x2pnt(temp_wave,l)] = V_avg
endfor

temp_wave-=fit_temp_wave


duplicate/O/R=((peak-0.1),(peak+1.5)) temp_wave, temp_wave_corr
setscale/P x,0,(deltax(temp_wave_corr)),"s",temp_wave_corr

variable filterfreq = 1000 * deltax(temp_wave_corr)

//remove residual 60 Hz noise
duplicate/o temp_wave_corr, temp_wave_base
jf_FastBPFilter(temp_wave_base, 59, 61,100)
temp_wave_corr -= temp_wave_base

//remove low frequency noise
duplicate/o temp_wave_corr, temp_wave_base
jf_fastfilter(temp_wave_base, 4, 80)
temp_wave_corr -= temp_wave_base


//low-pass filter at 0.8 kHz
jf_fastfilter(temp_wave_corr, 800,3) 


Differentiate temp_wave_corr/D=temp_wave_DIF


duplicate/o temp_wave, threshold
threshold=-level_dif


wavestats/Q/R=(0.5,1.5) temp_wave_DIF
level_DIF = V_rms * 2
Noisewave_DIF [i] = V_rms * 2 

wavestats/Q/R=(1.1,1.5) temp_wave_corr
level = V_rms * 2
Noisewave [i] = V_rms * 2

// set minimal amplitude thresholds

if( level < 35) // 35 pA for climbing fibers, 20 pA for parallel fibers
level = 35
endif

if( level_DIF < 35000) // 35000 pA/s for climbing fibers, 20000 pA/s for parallel fibers
level_DIF = 35000
endif

//Find all the events that meet threshold
FindLevels/q/D=times/edge=2 temp_wave_DIF,-level_DIF
// store the event times and assign them a counter (trial number)

variable Levels_found = V_levelsfound

for(j=0;j<V_levelsfound;j+=1)

if (j<(V_levelsfound-1))
if (times [j] + 0.001 > times [j+1])
continue
endif
endif

WaveStats/q/R=(times[j],times[j]+0.001) temp_wave_corr
minitimes[k] = V_minloc
times[j] = V_minloc 
miniamps[k]=-V_min - temp_wave_corr(times[j]-0.001)



minionset_Y [k] = (V_min - temp_wave_corr(times[j]-0.001)) * 0.2 + temp_wave_corr(times[j]-0.001)
findlevel/Q/R=((minitimes[k]),(minitimes[k]-0.001)) temp_wave_corr,minionset_Y [k]
minionset_X[k] = V_LevelX




if(miniamps[k]<level)//level)
miniamps[k] = NaN
minitimes[k] = NaN
times[j] = NaN
continue
endif

if (Times[j] < 0.11)
	if (Times[j] > (0.09))
	miniamps[k] = NaN
	minitimes[k] = NaN
	Times[j] = NaN
	continue
	else
	wavecounter[k]=i
	k+=1
	endif
else
wavecounter[k]=i
k+=1
endif




for(m=0;m<=(width/deltax(Avg_sEPSC)+1);m+=1)
Avg_sEPSC[m] += temp_wave_corr[x2pnt(temp_wave_corr,V_LevelX)-0.001/deltax(temp_wave_corr)+m]
endfor

endfor

endfor		//end of the for loop for trials


redimension/N=(k) wavecounter,minitimes,miniamps,minionset_X,minionset_Y

// Find parameters for bookkeeping 

Avg_sEPSC /= (k)

variable Avg_sEPSC_base = Avg_sEPSC[0]

Avg_sEPSC -= Avg_sEPSC_base

WaveStats/q Avg_sEPSC
setscale/P x,(-V_minloc),deltax(Avg_sEPSC),"s",Avg_sEPSC

variable sEPSC_amp = V_min


FindLevel/q/R=(0,-0.002) Avg_sEPSC, 0.2*sEPSC_amp 		//find the 20% rise to align AR and also to calculate onset
variable sEPSCtime=V_levelx

FindLevel/q/R=(0,-0.002) Avg_sEPSC, 0.8*sEPSC_amp 		//find the 80% to calculate rise
variable/g sEPSC_rise=(V_levelx-sEPSCtime)*1000



FindLevel/q/R=(0,0.05) Avg_sEPSC, 0.8*sEPSC_amp 		//find the 20% rise to align AR and also to calculate onset
variable sEPSCfitstart=V_levelx

CurveFit/Q  exp_XOffset Avg_sEPSC(sEPSCfitstart,) 

variable/g sEPSC_tau = W_coef[2] * 1000


// Find parameters for bookkeeping 

Make/o/n=1 minihist			//histogram for the 100 ms after stimulation

variable binsize = 0.001 // 1 ms bins
variable Duration = (numpnts(temp_wave_corr)-1) * deltax(temp_wave_corr)

Histogram/B={0,.001,1600} minitimes,minihist

duplicate/o minihist,minihistx
minihistx=x

wavestats/q minihist


minihist/=avgamp/1000				//normalize to the size of EPSC (proxy for number of synapses stimulated)
minihist/=(length)		//normalize to the number of trials


Make/o/n=1 miniAMPhist		//histogram for the quantal sizes
Histogram/B={0,1,500} miniamps, miniAMPhist


print sEPSC_rise
print sEPSC_tau


End
