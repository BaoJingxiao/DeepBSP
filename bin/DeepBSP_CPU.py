#!/usr/bin/env python
# coding: utf-8

#Import
import sys
import torch
from torch.utils.data import Dataset,DataLoader
import netCDF4
import numpy
import copy

#Settings
ComplexListFile=sys.argv[1]
NetParametersPath=sys.argv[2]
NetParametersNames=["Net_Data_1.pt","Net_Data_2.pt","Net_Data_3.pt","Net_Data_4.pt","Net_Data_5.pt"]
ResultFile=sys.argv[3]
BenchSize=int(sys.argv[4])
BenchNumWorkers=int(sys.argv[5])
TotalThreads=int(sys.argv[6])
GridNum=24
RotTypes=24

#Define function to rotate voxel data 90,180,270 degrees counterclockwise around X axis
def Rot_X(GridNum,InData,RotIndex):
    OutData=copy.deepcopy(InData)
    for y in range(GridNum):
        for z in range(GridNum):
            if (RotIndex==1):   #90 degrees
                OutData[:,:,GridNum-1-z,y]=InData[:,:,y,z]
            elif (RotIndex==2): #180 degrees
                OutData[:,:,GridNum-1-y,GridNum-1-z]=InData[:,:,y,z]
            else:               #270 degrees
                OutData[:,:,z,GridNum-1-y]=InData[:,:,y,z]
    return(OutData)

#Define function to rotate voxel data 90,180,270 degrees counterclockwise around Y axis
def Rot_Y(GridNum,InData,RotIndex):
    OutData=copy.deepcopy(InData)
    for x in range(GridNum):
        for z in range(GridNum):
            if (RotIndex==1):   #90 degrees
                OutData[:,z,:,GridNum-1-x]=InData[:,x,:,z]
            elif (RotIndex==2): #180 degrees
                OutData[:,GridNum-1-x,:,GridNum-1-z]=InData[:,x,:,z]
            else:               #270 degrees
                OutData[:,GridNum-1-z,:,x]=InData[:,x,:,z]
    return(OutData)

#Define function to rotate voxel data 90,180,270 degrees counterclockwise around Z axis
def Rot_Z(GridNum,InData,RotIndex):
    OutData=copy.deepcopy(InData)
    for x in range(GridNum):
        for y in range(GridNum):
            if (RotIndex==1):   #90 degrees
                OutData[:,GridNum-1-y,x,:]=InData[:,x,y,:]
            elif (RotIndex==2): #180 degrees
                OutData[:,GridNum-1-x,GridNum-1-y,:]=InData[:,x,y,:]
            else:               #270 degrees
                OutData[:,y,GridNum-1-x,:]=InData[:,x,y,:]
    return(OutData)

#Define function to rotate voxel data
def RotData(GridNum,InData,RotIndex):
    if (RotIndex//4==0):
        Temp=copy.deepcopy(InData)
        if (RotIndex%4==0):
            OutData=Temp                  #Front:1,Top:2
        elif (RotIndex%4==1):
            OutData=Rot_Z(GridNum,Temp,1) #Front:1,Top:4
        elif (RotIndex%4==2):
            OutData=Rot_Z(GridNum,Temp,2) #Front:1,Top:5
        else:
            OutData=Rot_Z(GridNum,Temp,3) #Front:1,Top:3
    elif (RotIndex//4==1):
        Temp=Rot_X(GridNum,InData,1)
        if (RotIndex%4==0):
            OutData=Temp                  #Front:2,Top:6
        elif (RotIndex%4==1):
            OutData=Rot_Z(GridNum,Temp,1) #Front:2,Top:4
        elif (RotIndex%4==2):
            OutData=Rot_Z(GridNum,Temp,2) #Front:2,Top:1
        else:
            OutData=Rot_Z(GridNum,Temp,3) #Front:2,Top:3
    elif (RotIndex//4==2):
        Temp=Rot_Y(GridNum,InData,1)
        if (RotIndex%4==0):
            OutData=Temp                  #Front:3,Top:2
        elif (RotIndex%4==1):
            OutData=Rot_Z(GridNum,Temp,1) #Front:3,Top:1
        elif (RotIndex%4==2):
            OutData=Rot_Z(GridNum,Temp,2) #Front:3,Top:5
        else:
            OutData=Rot_Z(GridNum,Temp,3) #Front:3,Top:6
    elif (RotIndex//4==3):
        Temp=Rot_Y(GridNum,InData,3)
        if (RotIndex%4==0):
            OutData=Temp                  #Front:4,Top:2
        elif (RotIndex%4==1):
            OutData=Rot_Z(GridNum,Temp,1) #Front:4,Top:6
        elif (RotIndex%4==2):
            OutData=Rot_Z(GridNum,Temp,2) #Front:4,Top:5
        else:
            OutData=Rot_Z(GridNum,Temp,3) #Front:4,Top:1
    elif (RotIndex//4==4):
        Temp=Rot_X(GridNum,InData,3)
        if (RotIndex%4==0):
            OutData=Temp                  #Front:5,Top:1
        elif (RotIndex%4==1):
            OutData=Rot_Z(GridNum,Temp,1) #Front:5,Top:4
        elif (RotIndex%4==2):
            OutData=Rot_Z(GridNum,Temp,2) #Front:5,Top:6
        else:
            OutData=Rot_Z(GridNum,Temp,3) #Front:5,Top:3
    else:
        Temp=Rot_X(GridNum,InData,2)
        if (RotIndex%4==0):
            OutData=Temp                  #Front:6,Top:5
        elif (RotIndex%4==1):
            OutData=Rot_Z(GridNum,Temp,1) #Front:6,Top:4
        elif (RotIndex%4==2):
            OutData=Rot_Z(GridNum,Temp,2) #Front:6,Top:2
        else:
            OutData=Rot_Z(GridNum,Temp,3) #Front:6,Top:3
    return(OutData)

#Define GetData class
class GetVoxelData(Dataset):
    def __init__(self,ComplexListFile,GridNum,RotIndex):
        self.ComNames=numpy.loadtxt(ComplexListFile,skiprows=0,usecols=[0],dtype=str,ndmin=1)
        self.VoxelDataFiles=numpy.loadtxt(ComplexListFile,skiprows=0,usecols=[1],dtype=str,ndmin=1)
        self.GridNum=GridNum
        self.RotIndex=RotIndex
        self.Constant=2**-15
    def __getitem__(self,index):
        ComName=self.ComNames[index]
        VoxelDataFile=self.VoxelDataFiles[index]
        #Read VoxelData file data
        VoxelDataSet=netCDF4.Dataset(VoxelDataFile)
        TempData=numpy.array(VoxelDataSet.variables['Data'][:],dtype=numpy.float32)
        TempData=TempData*self.Constant+0.5
        #Rotate the data around xyz
        VoxelData=RotData(self.GridNum,TempData,self.RotIndex[index])
        #Return data
        return(VoxelData)
    def __len__(self):
        return(self.ComNames.shape[0])

#Define Squeezenet's FireNet
class FireNet(torch.nn.Module):
    def __init__(self,inplanes,squeeze_planes,expand1_planes,expand3_planes):
        super(FireNet,self).__init__()
        self.inplanes=inplanes
        self.squeeze=torch.nn.Conv3d(inplanes,squeeze_planes,kernel_size=1)
        self.squeeze_active=torch.nn.ReLU(inplace=True)
        self.expand1=torch.nn.Conv3d(squeeze_planes,expand1_planes,kernel_size=1)
        self.expand1_active=torch.nn.ReLU(inplace=True)
        self.expand3=torch.nn.Conv3d(squeeze_planes,expand3_planes,kernel_size=3,padding=1)
        self.expand3_active=torch.nn.ReLU(inplace=True)
    def forward(self,x):
        x=self.squeeze_active(self.squeeze(x))
        return(torch.cat([self.expand1_active(self.expand1(x)),self.expand3_active(self.expand3(x))],1))
    
#Define kdeep network structure
class KdeepNet(torch.nn.Module):
    def __init__(self):
        super(KdeepNet,self).__init__()
        self.Net=torch.nn.Sequential(
            torch.nn.Conv3d(in_channels=16,out_channels=96,kernel_size=7,stride=2,padding=3),
            torch.nn.ReLU(inplace=True),
            FireNet(96,16,64,64),
            FireNet(128,16,64,64),
            FireNet(128,32,128,128),
            torch.nn.MaxPool3d(kernel_size=4,stride=2),
            FireNet(256,32,128,128),
            FireNet(256,48,192,192),
            FireNet(384,48,192,192),
            FireNet(384,64,256,256),
            torch.nn.AvgPool3d(kernel_size=3,stride=2),
            torch.nn.modules.Flatten(),
            torch.nn.Dropout(0.5),
            torch.nn.Linear(4096,1)
        )
    def forward(self,x):
        x=self.Net(x)
        return(x)

#Define function to predict with the network
def Predict(Net,GridNum,RotTypes,ComplexListFile,BenchSize,BenchNumWorkers):
    #Get RotIndexList
    TempComNames=numpy.loadtxt(ComplexListFile,skiprows=0,usecols=[0],dtype=str,ndmin=1)
    TempComNum=TempComNames.shape[0]
    RotIndexList=numpy.zeros((TempComNum,RotTypes),dtype=numpy.int8)
    for n in range(TempComNum):
        RotIndexList[n,:]=range(0,RotTypes)
    #Initialize arrays to store calculated values
    TempCal=numpy.zeros((TempComNum,RotTypes),dtype=numpy.float32)
    #Calculate predictions for each RotTypes
    with torch.no_grad():
        for RotIndex in range(RotTypes):
            TempDataSet=GetVoxelData(ComplexListFile,GridNum,RotIndexList[:,RotIndex])
            TempDataLoader=DataLoader(TempDataSet,batch_size=BenchSize,shuffle=False,num_workers=BenchNumWorkers,drop_last=False)
            Temp1=0
            for DataData in TempDataLoader:
                DataNum=len(DataData)
                DataData=DataData.float()
                pred=Net(DataData)
                TempCal[Temp1:(Temp1+DataNum),RotIndex]=pred.cpu().detach().numpy().transpose()
                Temp1=Temp1+DataNum
    #Average among all RotTypes
    AvgTempCal=numpy.mean(TempCal,axis=1)
    #Return result
    return(AvgTempCal)

#Set number of threads
torch.set_num_threads(TotalThreads)

#Get number of poses
TempComNames=numpy.loadtxt(ComplexListFile,skiprows=0,usecols=[0],dtype=str,ndmin=1)
PoseNum=TempComNames.shape[0]

#Predict with each network
Flag=0
CalResult=numpy.zeros([len(NetParametersNames),PoseNum])
for NetParametersName in NetParametersNames:
    #Load the network
    MyNet=KdeepNet()
    Temp=torch.load(''.join([NetParametersPath,'/',NetParametersName]),map_location='cpu')
    MyNet.load_state_dict(Temp['model'])
    #Set network to evaluate model
    MyNet.eval()
    #Predict with current network
    CalResult[Flag,:]=Predict(MyNet,GridNum,RotTypes,ComplexListFile,BenchSize,BenchNumWorkers)
    Flag=Flag+1

#Print averaged results to result file
OutFile=open(ResultFile,"w")
for n in range(PoseNum):
    print('%s %8.4f' %(TempComNames[n],numpy.mean(CalResult[:,n])),file=OutFile)
OutFile.close()

