import sys
import numpy as np
import gzip


input_bed=sys.argv[1]
inbed = []
with gzip.open(input_bed, "rt", encoding="utf-8") as bed:
	for line in bed:
		arr = line.rstrip().split("\t")
		inbed.append(list([arr[0],int(arr[1]),int(arr[2])])) #trans to three cols(chr,site,methylation)

input_mb=sys.argv[2]
inmb = []
with open(input_mb, "r") as bed:
        for line in bed:
                arr = line.rstrip().split("\t")
                inmb.append(list([arr[0],int(arr[1]),int(arr[2])])) #trans to three cols(chr,site,methylation)


fsrlist=[]
fsrlist.append(list(["chr","start","end","short_ratio","long_ratio"]))

k=0
for i in range(len(inmb)):
	n=0
	s=0
	l=0
	nnn=0
	for j in range(k,len(inbed)):
		#print(list([inbed[j],inmb[i]]))
		#print(inbed[j][0]==inmb[i][0])
		#print(inbed[j][1]>inmb[i][1])
		#print(inbed[j][1]<=inmb[i][2])
		#a=(inbed[j][0]==inmb[i][0]) & (inbed[j][1]>inmb[i][1]) & (inbed[j][1]<=inmb[i][2])
		#if((inbed[j][0]!=inmb[i][0])):
		#	next
		if((inbed[j][0]==inmb[i][0]) & (inbed[j][1]>inmb[i][1]) & (inbed[j][1]>inmb[i][2])):
			break
		if((inbed[j][0]==inmb[i][0]) & (inbed[j][1]>=inmb[i][1]) & (inbed[j][1]<=inmb[i][2])):
			n=n+1
			nnn=nnn+1
			#print(n)
			#print(inbed[j][2]<=150);print(inbed[j][2]>=151)
			if((inbed[j][2]>=90) & (inbed[j][2]<=150)):		#100-150bp
				s=s+1
				#print("s"+str(s))
			elif((inbed[j][2]>=151) & (inbed[j][2]<=220)):	#151-220bp
				l=l+1
				#print("l"+str(l))
			if(j+1>=len(inbed)):
				shor=s/n
				lon=l/n
				fsrlist.append(list([inmb[i][0],inmb[i][1],inmb[i][2],format(shor, '.3f'),format(lon, '.3f')]))
				k=k+n
				break
			elif((inbed[j+1][0]!=inmb[i][0]) | (inbed[j+1][1]>inmb[i][2])):
				shor=s/n
				lon=l/n
				fsrlist.append(list([inmb[i][0],inmb[i][1],inmb[i][2],format(shor, '.3f'),format(lon, '.3f')]))
				k=k+n
				#print(list([j,i,inbed[j],inmb[i],nnn,shor,lon,int(n)]))
				break
			
		


result=fsrlist
f = open(sys.argv[3], 'w')
for i in range(len(result)):
    for j in range(len(result[i])):
        f.write(str(result[i][j]))  # 
        if (j < len(result[i]) - 1):
            f.write('\t')
    f.write('\n')  #
f.close()


