import matplotlib.pyplot as plt
import numpy as np 

input1 = "eff_layers/sector13_layer0_fixed.txt"
input2 = "ped_layers/sector13_layer0_threshold_vs_channel.txt"

handle1 = open(input1)
handle2 = open(input2) 

channels_one = []
channels_two = []
thresholds = [] 
effs = [] 

lines1 = handle1.readlines()
lines2 = handle2.readlines()

for line in lines1:
	splitLine = line.split()
	channel_1 = splitLine[0]
	channels_one.append(channel_1)
	eff = splitLine[1]
	effs.append(eff)

for otherLine in lines2:
	otherSplit = otherLine.split()
	channel_2 = otherSplit[0] 
	channels_two.append(channel_2) 
	thresh = otherSplit[1]
	thresholds.append(thresh) 


# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].scatter(channels_one, effs)
axarr[0].set_title("Sector -3, Layer 1")
axarr[1].scatter(channels_two, thresholds)

axarr[0].axis([-50, 1, 0.0, 1.0]) 
#axarr[1].axis([-1, 50, 0, 32]) 

major_ticks = np.arange(-50, 0, 5)                                              
minor_ticks = np.arange(-50, 0, 1)
axarr[0].set_xticks(major_ticks)                                                       
axarr[0].set_xticks(minor_ticks, minor=True)

axarr[1].set_xlabel("Phi Channel")
axarr[1].set_ylabel("Threshold")

axarr[0].set_ylabel("Phi Efficiency")

plt.savefig("sector13_layer0.png")



                                            

                                           
