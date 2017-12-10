import Functions as F
from PIL import Image
import random
import os
import time


######### VARIOGRAM-BASED DIRECT SAMPLING WITH PATCHES ###############

class Pixel:
    def __init__(self,SimTrack,SimQual,NeighbsUsed,Zipcode,
    IsConditional, IsPerimeter, SimNumber):
        self.SimTrack = SimTrack
        self.SimQual = SimQual
        self.NeighbsUsed = NeighbsUsed
        self.Zipcode = Zipcode
        self.IsConditional = IsConditional
        self.IsPerimeter = IsPerimeter
        self.SimNumber = SimNumber


# USER-SET PARAMETERS - CHANGE THESE PARAMETERS FOR DIFFERENT RESULTS
ImgName = "Channels125x125"
SG_Name = "SG"
WithoutDistWeight = True
All_ReqNeighbs = [16, 16, 16]; # Required number of matching neighbors in TI required to copy pixel
PatchSize = [1, 3, 9];
threshold = 0.0; # threshold of 0 means ALL neighbors must match exactly
delta = 0.5
PatchLevels = [0.0, 0.04, 0.2, 1.0]
SimOrder_nConditional = 2
MaxDist = 40
MaxArea = 9
Cols = 30
Rows = 30

# LOAD AN IMAGE, RETRIEVE IMPORTANT PARAMETERS
TI = Image.open(ImgName+".png")
PixTI = TI.load()
ColsTI = TI.size[0]
RowsTI = TI.size[1]
TI_Tracker = [[0]*RowsTI for i in range(ColsTI)]
#TI_Pixels = TI.load()

# DEFINE IMPORTANT STUFF
random.seed(a=None)
nConditional = 1

# Retrieve variogram-information for weighting
Directions = [[1, 0], [2, 1], [1, 1], [1, 2], [0, 1], [-1, 2], [-1, 1], [2, -1]]
Directions, Angles, Splines = F.Variogram_Splines(ImgName, Directions, MaxDist)

# COUNTERS AND LOGICALS
Tic = time.clock()
starttime = time.clock()
i = nConditional
n_PixelsSimmed = nConditional
NTotPixels = Rows*Cols
Begin_Patching = PatchLevels[1]*NTotPixels
Time_counter=1

# CALCULATES PATCH-LOOKS AND BEST CONDITIONALS
Patches, All_BestNeighbs = F.FindPatchNeighbs(PatchSize, MaxDist, Splines, Angles, delta)
i_Patch_Level = 0
while PatchLevels[i_Patch_Level+1]<n_PixelsSimmed/float(NTotPixels):
    i_Patch_Level += 1
Patch = Patches[i_Patch_Level];
BestNeighbs = All_BestNeighbs[i_Patch_Level]
ReqNeighbs = All_ReqNeighbs[i_Patch_Level]

# MAKE SIMGRID (SG) AND A PIXELTRACKER
Pixels = [[Pixel(0, 1, [0]*ReqNeighbs, 0, 0, 0, 0) for i in range(Rows)] for  i in range(Cols)]
SimOrder = [0]*(Rows*Cols)
SG = Image.new( 'RGB', (Cols,Rows), "white")
PixSG = SG.load()
with open(SG_Name + ".txt","r") as f:
    n_PixelsSimmed = 0
    for line in f:
        Line = line.split()
        x = int(Line[0])
        y = int(Line[1])
        PixSG[x, y] = (int(Line[2]), int(Line[3]), int(Line[4]))
        Pixels[x][y].SimTrack = 1
        Pixels[x][y].IsConditional = 1
        Pixels[x][y].SimNumber = n_PixelsSimmed+1
        n_PixelsSimmed += 1

# LOAD SIMORDER
SimOrder = [0]*(Cols*Rows)
with open("SimOrder.txt", "r") as f:
    c = 0
    for line in f:
        Line = line.split()
        SimOrder[c] = [int(Line[0]), int(Line[1])]
        c += 1

# Begin simulating
i = n_PixelsSimmed
while n_PixelsSimmed < NTotPixels:

    # CHOOSE PIXEL TO SIMULATE
    PixToSim = SimOrder[i];

    # UPDATE PATCH-FORM
    if PatchLevels[i_Patch_Level+1]<n_PixelsSimmed/float(NTotPixels):
        i_Patch_Level += 1
        Patch = Patches[i_Patch_Level];
        BestNeighbs = All_BestNeighbs[i_Patch_Level]
        ReqNeighbs = All_ReqNeighbs[i_Patch_Level]

    # CONTINUE TO NEXT PIXEL IS ALREADY SIMULATED
    if (Pixels[PixToSim[0]][PixToSim[1]].SimTrack == 1 and i+1 != len(SimOrder)):
        i+=1
        continue

    # CHECKS IF ENCLOSED PATCH
    Area, Checked, Perimeter, Asym_Patch = F.FindUnsimmedArea(Pixels, [], PixToSim, 0, MaxArea, [], [], PixToSim, Cols, Rows)

    # FINDS N NEAREST PIXELS IN SG
    if Area <= MaxArea:
        LagVecs = Perimeter
        Z = [0]*len(LagVecs)
        Weights = [1.0]*len(LagVecs)
        a, b, c, d = 0, 0, 0, 0
        for ac in range(len(LagVecs)):
            Z[ac] = PixSG[PixToSim[0]+LagVecs[ac][0], PixToSim[1]+LagVecs[ac][1]]
            if a > LagVecs[ac][1]: a = LagVecs[ac][1]
            if b < LagVecs[ac][0]: b = LagVecs[ac][0]
            if c < LagVecs[ac][1]: c = LagVecs[ac][1]
            if d > LagVecs[ac][0]: d = LagVecs[ac][0]
        abcd = [a, b, c, d]
    else:
        (SG, LagVecs, abcd, Z, QuadCount, Weights) = F.FindNearestPixels(SG, BestNeighbs, PixToSim, ReqNeighbs, Pixels, n_PixelsSimmed, delta, WithoutDistWeight)

    # FINDS A NODE IN TI WITH SAME NEAREST PIXELS AS IN SG
    SimmedPixel, Checker, Node, SimQual = F.SimGridNode(TI, LagVecs, abcd,
    Z, threshold, Weights)

    if SimmedPixel == (255, 255, 255): SimmedPixel = (255, 0, 0)

    # INSERT IN SG
    if Checker == 1:
        if Area <= MaxArea: P = Asym_Patch
        else: P = Patch
        Lap_Start = time.clock()
        PixSG, Pixels, n_PixelsSimmed, TI_Tracker = F.InsertPatch(PixSG, Cols, Rows, PixToSim,
        PixTI, ColsTI, RowsTI, LagVecs, Pixels, Node, P,
        n_PixelsSimmed, TI_Tracker, SimQual)
    else:
        SimOrder.append(PixToSim)

    # Print progression
    if n_PixelsSimmed >= NTotPixels*(float(Time_counter)/10):
        Toc = time.clock()
        Time_counter+=1
        print ""
        print "Simulated ", n_PixelsSimmed, " of ", NTotPixels, " pixels."
        print "Laptime: ", Toc-Tic
        Tic = Toc

    # Moves to next pixel
    i += 1

SG.show()
print "Your program ran for",time.clock()-starttime, "seconds"
