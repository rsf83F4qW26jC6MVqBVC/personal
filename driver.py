#!/opt/local/bin/python

import sys
from datetime import datetime
from dateutil.parser import parse
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as pdfwriter


#==========================================================================
def pause(location=None):
    if (location==None):
        input("press the <ENTER> key to continue . . .")
    else:
        input("(DEBUG:"+str(location)+") press the <ENTER> key to continue . . .")


#==========================================================================
class mesh:

    def __init__(self,NumEdgeElems,PlotInit=False):
        print("initializing mesh . . . ",end="")
        self.coords,self.elems=self.InitMesh(NumEdgeElems)
        # print("PlotInit in \"mesh\": {}".format(PlotInit))
        if (PlotInit):
            # pause("starting \"mesh\" initialization plot")
            fig=plt.figure("mesh initialization",figsize=(9,9))
            self.plot()
            plt.show()
            # pause("finishing \"mesh\" initialization plot")
        print("done")

    def InitMesh(self,NumEdgeElems):

        coords=np.zeros(shape=((NumEdgeElems+1)*(NumEdgeElems+1),2))
        elems=np.zeros(shape=(4,NumEdgeElems*NumEdgeElems),dtype=np.int64)

        h=1.0/NumEdgeElems
        iNode=-1
        iElem=-1
        for j in range(NumEdgeElems+1):
            for i in range(NumEdgeElems+1):
                # print("i={}".format(i))
                iNode=iNode+1
                coords[iNode,0]=i*h
                coords[iNode,1]=j*h
                if (i>0 and j>0):
                    iElem=iElem+1
                    elems[0,iElem]=iNode-(NumEdgeElems+1)-1
                    elems[1,iElem]=iNode-(NumEdgeElems+1)
                    elems[2,iElem]=iNode
                    elems[3,iElem]=iNode-1
                # if (16*6<iElem and iElem<16*7):
                    # print("iElem={}, ElemNodes={}".format(iElem,elems[:,iElem]))

        return coords,elems

    def plot(self):

        NumElems=self.elems.shape[1]
        # print("NumElems={}".format(NumElems))
        NumElemNodes=self.elems.shape[0]

        xplot=list()
        yplot=list()

        for iElem in range(NumElems):
        # for iElem in range(16*7-7):
            # print("i={}".format(i))
            for iNode in range(NumElemNodes):
                NodeNum=self.elems[iNode,iElem]
                # print("iElem={}, iNode={}, NodeNum={}".format(iElem,iNode,NodeNum))
                xplot.append(self.coords[NodeNum,0])
                yplot.append(self.coords[NodeNum,1])
            NodeNum=self.elems[0,iElem]
            xplot.append(self.coords[NodeNum,0])
            yplot.append(self.coords[NodeNum,1])
            # print("iElem={}, x_coords={}, y_coords={}".format(iElem,xplot,yplot))

            # plt.plot(xplot,yplot,marker='o',color="black")
            plt.plot(xplot,yplot,color="black")

            # delete elements of a list, but leave the list itself intact
            xplot.clear()
            yplot.clear()

        plt.grid(True)
        ax=plt.gca()
        ax.set_aspect('equal')


#==========================================================================
class fields:

    def __init__(self,InitField,BCField,m,PlotInit=False):
        print("initializing field . . . ",end="")
        self.fieldNP1=self.InitField(InitField,BCField,m.coords)
        if (PlotInit):
            # pause("beginning of PlotInit")
            fig=plt.figure("field initialization",figsize=(9,9))
            # print("field figure number: {}".format(fig.number))
            self.plot(m)
            # pause("after field initialization plot . . .")
            # print("fieldNP1: {}".format(self.fieldNP1))
            plt.show()
        print("done")
        print("")

    def InitField(self,InitField,BCField,coords):

        NumNodes=coords.shape[0]
        FieldNP1=np.zeros(shape=(NumNodes))
        for iNode in range(NumNodes):
            FieldNP1[iNode]=InitField
            if (coords[iNode,1]<1.0e-05):
                FieldNP1[iNode]=BCField
        # print("within InitField:")
        # print("InitField={}".format(FieldNP1))

        return FieldNP1

    def plot(self,m):

        # stupid ass arguments . . .
        NumNodes=m.coords.shape[0]
        NumEdgeNodes=int(math.sqrt(NumNodes))
        xplot=np.arange(NumEdgeNodes)/(NumEdgeNodes-1)
        yplot=xplot
        xplot1,yplot1=np.meshgrid(xplot,yplot)
        zplot=self.fieldNP1.reshape(NumEdgeNodes,NumEdgeNodes)
        # zplot=self.fieldNP1
        # print("after call to reshape:")
        # print("fieldNP1: {}".format(self.fieldNP1))
        # print("zplot: {}".format(zplot))

        plt.contourf(xplot1,yplot1,zplot)
        # plt.show()

        ax=plt.gca()
        ax.set_aspect('equal')


#==========================================================================
# fig=plt.figure(figsize=(9,9))



# xplot=list()
# yplot=list()
# legend_handles=list()
# legend_names=list()
# for i in range(40):
#     print("i={}".format(i+1))

#     x=i
#     term1=(1+math.exp(-0.3*(x-20)))**(-1/7)
#     f=1.0*(1.0-term1)

#     xplot.append(x)
#     yplot.append(f)


# ph, = plt.plot(xplot,yplot,marker='o')
# legend_handles.append(ph)


# date_string=datetime.now()
# date_string=datetime.strftime(date_string,"%Y-%m-%d")
# print("title date string: {}".format(date_string))
# plt.title('title title title, as of '+date_string)
# plt.xlabel('x label')
# plt.ylabel('y label')
# plt.grid(True)
# plt.xticks([1,30,60,90,120,150,180,210,240,270,300,330,365])
# print("legend handles: {}".format(legend_handles))
# print("legend names: {}".format(legend_names))
# plt.legend(legend_handles,legend_names)
# plt.show()

# print("")



#==========================================================================
def plotter(m,f,time,PlotMesh,PlotField):
    if (PlotMesh):
        m.plot()
    if (PlotField):
        f.plot(m)
    plt.title("time={:12.5e}".format(time))
    plt.show()


#==========================================================================
# user inputs
#==========================================================================
NumEdgeElems=13

InitialTemp=303.0
BCTemp=505.0

MaxNumSteps=10^10
EndTime=1.0
dt=0.1
dt=0.5000000000000001

MaxNonlinIters=20


#==========================================================================
# initialization
#==========================================================================
m=mesh(NumEdgeElems,PlotInit=False)
f=fields(InitialTemp,BCTemp,m,PlotInit=False)

fig=plt.figure("primary figure window",figsize=(9,9))
plotter(m,f,0.0,True,True)


#==========================================================================
# main time stepping loop
#==========================================================================
StepNum=0
Time=0.0
while True:
    StepNum=StepNum+1

    # deal with the time stepping
    Time=Time+dt

    print("Starting step {:8d} at time {:12.5e}".format(StepNum,Time))

    iNonlin=0
    while True:
        iNonlin=iNonlin+1

        # UpdateElems
        # SolveLinearSystem
        # UpdateField
        # CheckConvergnce
    
        if (iNonlin==MaxNonlinIters):
            print("ERROR: nonlinear solve not converging, quitting . . .")
            print("")
            sys.exit(0)

    print("")

    # termination criteria
    # number of steps
    if (StepNum==MaxNumSteps):
        break
    # end of time
    if (Time>=EndTime):
        break



