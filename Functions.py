import os
import time
import random
from PIL import Image
import ROOT as R
from array import array

def Length(Vector):
    return (Vector[0]**2 + Vector[1]**2)**0.5
def Angle(Vector):
    Ang = (180/3.14159)*R.acos(Vector[0]/Length(Vector))
    if Vector[1]<0:
        return 180 - Ang
    return Ang
def FindPatchNeighbs(PatchSizes, MaxDist, Splines, Angles, delta):
    # Calculates all patches being used in simulation before actually starting
    BestNeighbs_Original = HighestCorrelation(MaxDist-5**0.5, Splines, Angles, delta)

    Patches = []; BestNeighbs = [];
    for PatchSize in PatchSizes:
        BestNeighbs_Copy = BestNeighbs_Original[:]
        Patch = BestNeighbs_Copy[0:PatchSize-1]

        for i in range(len(BestNeighbs)):
            BestNeighbs_Copy = AddPatchCorrelations(BestNeighbs_Copy, Patch, Splines, Angles, delta)

        Patch = [[v[0], v[1]] for v in Patch]

        Patches.append(Patch); BestNeighbs.append(BestNeighbs_Copy[:])

    return Patches, BestNeighbs
def Variogram_Splines(TI_Name, Directions, MaxDist):

    if not os.path.exists(TI_Name+"_Variograms_Dist"+str(MaxDist)+".txt"):
        print "Variograms not calculated. Calculation begun"
        Tic = time.clock()
        Variogram_MaxDist_Write_Data(TI_Name, Directions, MaxDist)
        print "Variogram Calculation time: ", time.clock()-Tic
    Splines = []
    Directions = []
    Angles = []
    Vario = []
    Distances = []

    with open(TI_Name+"_Variograms_Dist"+str(MaxDist)+".txt","r") as f:
        for line in f:
            Line = line.split()
            Directions.append([float(Line[0]), float(Line[1])])
            Angles.append(Angle(Directions[-1]))
            Dist = []
            Var = []
            n = 2
            while n<len(Line):
                Dist.append(float(Line[n]))
                n+=1
                Var.append(float(Line[n]))
                n+=1
            Distances.append(Dist)
            Vario.append(Var)

    for n in range(len(Vario)):
        x = array("f", Distances[n])
        y = array("f", Vario[n])
        graph = R.TGraphErrors(len(Vario[n]),x,y)
        Splines.append(R.TSpline3("spline"+str(n),graph))

    return Directions, Angles, Splines
def Variogram_MaxDist_Write_Data(TI_Name, Directions, MaxDist):
    TI = Image.open(TI_Name+".png")
    Vario, Distance = Variogram_MaxDist(TI, Directions, MaxDist)


    with open(TI_Name+"_Variograms_Dist"+str(MaxDist)+".txt", "w") as f:
        for i in range(len(Directions)):
            f.write("%s %s "%(Directions[i][0], Directions[i][1]))
            for n in range(len(Vario[i])):
                f.write("%s %s "%(Distance[i][n], Vario[i][n]))
            f.write("\n")
def Variogram_MaxDist(Image, Directions, MaxDist):
    Pixels = Image.load()
    Cols = Image.size[0]
    Rows = Image.size[1]

    x_Max = Cols-1
    y_Max = Rows-1

    MaxSteps = [int(MaxDist/Length(Dir)) for Dir in Directions]

    Count=[[0]*(MaxSteps[i]+1) for i in range(len(MaxSteps))]
    Var= [[0]*(MaxSteps[i]+1) for i in range(len(MaxSteps))]
    Tot = [[0]*(MaxSteps[i]+1) for i in range(len(MaxSteps))]
    Dists = []#[[0]*(MaxSteps[i]+1) for i in range(len(MaxSteps))]
    Variogram = [[0]*(MaxSteps[i]+1) for i in range(len(MaxSteps))]



    Ave = 0
    Variance = 0
    Total = 0.0

    for x0 in range(Cols):
        for y0 in range(Rows):
            Start = 1;
            if Pixels[x0, y0] == (0, 0, 0): Start = 0
            Ave+=Start
            Variance+=Start**2
            Total+=1.0
            for n_dir in range(len(Directions)):

                for n in range(1, len(Count[n_dir])):
                    x = x0 + Directions[n_dir][0]*n;
                    y = y0 + Directions[n_dir][1]*n

                    if x < 0 or y < 0 or x > x_Max or y > y_Max: continue
                    End = 1
                    if Pixels[x,y] == (0, 0, 0): End = 0

                    Count[n_dir][n]+=Start-End
                    Var[n_dir][n]+=(Start-End)**2
                    Tot[n_dir][n]+=1

    for n_dir in range(len(Directions)):
        Dists.append([n*Length(Directions[n_dir]) for n in range(MaxSteps[n_dir]+1)])
        for i in range(1, len(Count[n_dir])):
            if Tot[n_dir][i] == 0: break
            Count[n_dir][i] = Count[n_dir][i]/float(Tot[n_dir][i])
            Var[n_dir][i] = Var[n_dir][i]/float(Tot[n_dir][i])
            Variogram[n_dir][i] = (Var[n_dir][i]-Count[n_dir][i]**2)/2.0

    # Ave = Ave/Total; Variance = Variance/Total
    # Variance = Variance - Ave**2
    # Distance = [i*(Direction[0]**2+Direction[1]**2)**0.5 for i in range(Steps+1)]
    #
    return Variogram, Dists
def Corr(x, y, Splines, Angles, delta):
    L = Length([x, y])
    if x== 0 and y == 0:
        return 0
    Angle = (180/3.14159)*R.acos(x/L)
    if y<0:
        Angle = 180 - Angle

    i = 1
    while True:
        if i == len(Angles):
            Lower = Angles[i-1]
            Lower_Index = i-1
            Upper = 180.0
            Upper_Index = 0
            break
        if Angle > Angles[i]:
            i+=1
        else:
            Lower = Angles[i-1]
            Lower_Index = i-1
            Upper = Angles[i]
            Upper_Index = i
            break

    Lower_Diff = Angle-Lower;
    Upper_Diff = Upper-Angle;
    Sum_Diff = Lower_Diff+Upper_Diff
    Lower_Weight = Upper_Diff/Sum_Diff
    Upper_Weight = Lower_Diff/Sum_Diff

    Weight = ((abs(Splines[Lower_Index].Eval(L))**0.5)*Lower_Weight + (abs(Splines[Upper_Index].Eval(L))**0.5)*Upper_Weight)*(Length([x, y]))**delta
    if Weight == 0.0:
        return 10**3
    return 1/Weight
def HighestCorrelation(Within_Radius, Splines, Angles, delta):

    MaxR = int(Within_Radius)
    HCN = [] # Highest Correlation Neighbours


    Weight = []

    for x in range(0, MaxR+1, 1):
        y = 0

        while Length([x, y]) <= Within_Radius:
            if y!=0: HCN.append([x, y, Corr(x, y,Splines, Angles, delta)])
            if x!=0: HCN.append([-x, y, Corr(-x, y,Splines, Angles, delta)])
            if y!=0: HCN.append([-x,-y, Corr(-x, -y,Splines, Angles, delta)])
            if x!=0: HCN.append([x, -y, Corr(x, -y,Splines, Angles, delta)])
            y+=1

    #del HCN_NE[0]
    HCN = sorted(HCN, key = lambda x: x[2], reverse = True)

    return HCN
def AddPatchCorrelations(List, Patch, Splines, Angles, delta):
    #To_Add = List[0:Patchsize-1]
    #print len(To_Add)
    #print To_Add
    for V in Patch:
        for Entry in List:
            if Entry[0] == V[0] and Entry[1]==V[1]: continue
            x = Entry[0]-V[0]
            y = Entry[1]-V[1]
            Correlation = Corr(x, y, Splines, Angles, delta)
            if Correlation > Entry[2]: Entry[2] = Correlation

    return sorted(List, key = lambda x: x[2], reverse = True)
def FindUnsimmedArea(Pixels, Checked, Pixel, Area, Max_Area, Perimeter, Patch, Orig_Pixel, Cols, Rows):

    # Only enters function if unsimmed pixel
    x0 = Pixel[0]; y0 = Pixel[1]
    Patch.append([x0-Orig_Pixel[0], y0-Orig_Pixel[1]])
    Area += 1
    Checked.append([x0, y0])

    for x in range(-1,2,1):
        for y in range(-1, 2, 1):

            # Return if area too big, checking itself or outside grid
            if Area > Max_Area:
                return Area, Checked, Perimeter, Patch
            if x == 0 and y == 0: continue
            if x0+x >= Cols or y0+y >= Rows or x0+x<0 or y0+y < 0: continue

            # Only continue of pixel hasnt been checked
            if [x0+x, y0+y] not in Checked:

                # Continue of not "direct" neighbor, but perhaps add as perimeter
                if x != 0 and y != 0:
                    if Pixels[x0+x][y0+y].SimTrack == 1:
                        Perimeter.append([x0+x-Orig_Pixel[0], y0+y-Orig_Pixel[1]])
                        Checked.append([x0+x, y0+y])
                    continue

                # If already simmed, add to checked
                if Pixels[x0+x][y0+y].SimTrack == 1:
                    Checked.append([x0+x, y0+y])
                    Perimeter.append([x0+x-Orig_Pixel[0], y0+y-Orig_Pixel[1]])
                #
                else:
                    Area, Checked, Perimeter, Patch = FindUnsimmedArea(Pixels, Checked, [x0+x, y0+y], Area, Max_Area, Perimeter, Patch, Orig_Pixel, Cols, Rows)

    return Area, Checked, Perimeter, Patch
def FindNearestPixels(SG, BestNeighbs, PixToSim,
ReqNeighbs, Pixels, n_PixSimmed, delta, WithoutDistWeight):
    # Define stuff and paint pixel being simmed
    x0 = PixToSim[0]; y0 = PixToSim[1];

    x_max = SG.size[0]-x0;
    y_max = SG.size[1]-y0;

    SG_Pixels = SG.load()
    LagVecs, Z, Init_Weights = [0]*min(ReqNeighbs,n_PixSimmed), [0]*min(ReqNeighbs,n_PixSimmed), [0]*min(ReqNeighbs,n_PixSimmed);
    QuadCount = [0, 0, 0, 0]
    count, VecCount, a, b, c, d = 0, 0, 0, 0, 0, 0;

    while LagVecs[-1] == 0:
        if count == len(BestNeighbs):
            # print "Entered this shit"
            LagVecs = [Entry for Entry in LagVecs if Entry != 0]
            Init_Weights = [Entry for Entry in Init_Weights if Entry != 0]
            break

        x = BestNeighbs[count][0]; y = BestNeighbs[count][1];

        if WithoutDistWeight:
            weight = BestNeighbs[count][2]*Length([x, y])**delta
        else:
            weight = BestNeighbs[count][2]

        count+=1

        # Ensures not to check outside SG
        if x >= x_max or x < -x0 or y >= y_max or y < -y0:
            continue

        # If a "best" pixel is found, save the lag vector
        # Find a, b, c, d to define search window now
        if Pixels[x0+x][y0+y].SimTrack == 1:
        # if Pixels[x0+x][y0+y] == 1:

            # Note which quadrant it is in
            if x > 0 and y >= 0:
                QuadCount[0]+=1
            elif x<=0 and y>0:
                QuadCount[1]+=1
            elif x<0 and y<=0:
                QuadCount[2]+=1
            else:
                QuadCount[3]+=1

            LagVecs[VecCount] = [x,y]
            Init_Weights[VecCount] = weight
            Z[VecCount] = SG_Pixels[x0+x,y0+y]

            if a > y: a = y
            if b < x: b = x
            if c < y: c = y
            if d > x: d = x

            VecCount+=1
            # SG_Pixels[x0+x,y0+y] = Purple

    Weights = Init_Weights
    abcd = [a, b, c, d]

    return SG, LagVecs, abcd, Z, QuadCount, Weights
def SimGridNode(TI, LagVecs, abcd, Z, threshold, Weights):
    # If no perfect found, take best
    Checker = 1

    # CALCULATE search window
    TI_Pixels = TI.load()
    X = TI.size[0]; Y = TI.size[1]
    SW_xDim = X-(abcd[1]-abcd[3])
    SW_yDim = Y-(abcd[2]-abcd[0])

    # Ensure SW is well defined
    if SW_xDim < 1 or SW_yDim < 1:
        SimmedPixel = 0
        Node = 0
        Best = 0
        return SimmedPixel, Checker, Node, Best

    # Make list with candidates, transform list to be in SW
    SW = [[i%SW_xDim - abcd[3], i/SW_xDim - abcd[0]] for i in range(SW_xDim*SW_yDim)]

    # LINEARLY scan search window starting at random location
    # ASSUMES GRAYSCALE TI
    Start = random.randint(0,len(SW)-1)
    Best = 1.0
    count = 0
    for n in range(Start, -len(SW) + Start, -1):
        # Calculates a "distance" and keeps the best
        # Breaks if perfect match
        Distance = CalcDistBW(TI_Pixels, LagVecs, Z, SW[n], Weights)

        if Distance < Best:
            Best = Distance
            Node = [SW[n][0], SW[n][1]]
            SimmedPixel = TI_Pixels[SW[n][0], SW[n][1]]
        #  Distance
        if Distance <= threshold:

            Checker = 1
            break

    return SimmedPixel, Checker, Node, Best
def CalcDistBW(TI_Pixels, LagVecs, Z, Node, Weights):
    Length = 0; x0 = Node[0]; y0 = Node[1]

    for k in range(len(LagVecs)):
        x = LagVecs[k][0]; y = LagVecs[k][1]

        if TI_Pixels[x0+x,y0+y][0] != Z[k][0]:
            Length+=1*Weights[k]

    if len(Weights) == 0:
        Length = 0.0
    else:
        Length = Length/sum(Weights)

    return Length
def InsertPatch(Pixels_SG, XMax_SG, YMax_SG, PixelToSim, Pixels_TI, XMax_TI,
YMax_TI, LagVecs, SimTracker, TI_Node, NearestNeighbors, n_PixelsSimmed,
TI_Tracker, SimQual):
    x0_SG = PixelToSim[0]; y0_SG = PixelToSim[1]
    x0_TI = TI_Node[0]; y0_TI = TI_Node[1]

    # Insert the matching node itself first
    if Pixels_TI[x0_TI, y0_TI] == (255, 255, 255):
        Pixels_SG[x0_SG, y0_SG] = (255, 0, 0)
    else:
        Pixels_SG[x0_SG, y0_SG] = Pixels_TI[x0_TI, y0_TI]

    SimTracker[x0_SG][y0_SG].SimTrack = 1
    SimTracker[x0_SG][y0_SG].NeighbsUsed = LagVecs
    SimTracker[x0_SG][y0_SG].SimQual = SimQual
    SimTracker[x0_SG][y0_SG].SimNumber = n_PixelsSimmed


    TI_Tracker[x0_TI][y0_TI]+=1
    n_PixelsSimmed+=1

    # Insert remaining patch
    for n in NearestNeighbors:
        x_SG = x0_SG+n[0]; y_SG = y0_SG+n[1]
        x_TI = x0_TI+n[0]; y_TI = y0_TI+n[1]

        # Ensure patch is within TI and SG and no overwritings are made
        if (x_SG < 0 or y_SG < 0 or x_TI < 0 or y_TI < 0 or x_SG >= XMax_SG
        or y_SG >= YMax_SG or x_TI >= XMax_TI or y_TI >= YMax_TI):
        #or SimTracker[x_SG][y_SG]==1):
            continue

        # Does NOT overwrite earlier simulations
        if SimTracker[x_SG][y_SG].SimTrack == 0:
            if Pixels_TI[x_TI, y_TI] == (255, 255, 255):
                Pixels_SG[x_SG, y_SG] = (255, 0, 0)
            else:
                Pixels_SG[x_SG, y_SG] = Pixels_TI[x_TI, y_TI]

            n_PixelsSimmed+=1
            SimTracker[x_SG][y_SG].SimNumber = n_PixelsSimmed


        SimTracker[x_SG][y_SG].SimTrack = 1
        TI_Tracker[x_TI][y_TI]+=1

    return Pixels_SG, SimTracker, n_PixelsSimmed, TI_Tracker
