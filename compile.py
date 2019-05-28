#compile code for project 6
import os

os.system("rm project6")

blocksize = [16, 32, 64]
numTrails = [16384, 32768, 65536, 131072, 262144, 524288]

print("NumTrials\tBlockSize\tPerformance\tNumHits\tFrequency")

for b in blocksize:
	for n in numTrails:
		cmd = "touch project6.cu"
		os.system(cmd)
		cmd="make CXXFLAGS=-DBLOCKSIZE=%d CXXFLAGS+=-DNUMTRIALS=%d > make_output.txt" % (b, n)
		os.system(cmd)
		cmd="./project6 > stdout.txt"
		os.system(cmd)

