#!/opt/local/bin/python

import sys
from datetime import datetime
from dateutil.parser import parse
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages as pdfwriter


#==========================================================================
def pause(location=None):
    if (location==None):
        input("press the <ENTER> key to continue . . .")
    else:
        input("(DEBUG:"+str(location)+") press the <ENTER> key to continue . . .")


#==========================================================================
class ParametricSpace:

    def __init__(self,NumIntPoints):
        print("initializing parametric space . . . ",end="")
        self.NumIntPoints=NumIntPoints
        self.p,self.w=self.DefineIntPoints(NumIntPoints)
        self.N,self.dNdxi,self.dNdeta=self.EvalInterpFuncs(NumIntPoints,self.p,4)
        print("done")

    def DefineIntPoints(self,NumIntPoints):
        p=np.empty(shape=(2,NumIntPoints),dtype=np.float64)
        w=np.empty(shape=(NumIntPoints),dtype=np.float64)
        if (NumIntPoints==4):
            iip=0
            p[0,iip]=-1.0/sqrt(3.0)
            p[1,iip]=-1.0/sqrt(3.0)
            w[iip]=1.0

            iip=1
            p[0,iip]=+1.0/sqrt(3.0)
            p[1,iip]=-1.0/sqrt(3.0)
            w[iip]=1.0

            iip=2
            p[0,iip]=+1.0/sqrt(3.0)
            p[1,iip]=+1.0/sqrt(3.0)
            w[iip]=1.0

            iip=3
            p[0,iip]=-1.0/sqrt(3.0)
            p[1,iip]=+1.0/sqrt(3.0)
            w[iip]=1.0

        else:
            print("ERROR: unknown number of integration points, quitting . . .")
            sys.exit(-1)
        return p,w
            

    def EvalInterpFuncs(self,NumIntPoints,p,NumNodes):
        N=np.empty(shape=(NumNodes,NumIntPoints),dtype=np.float64)
        dNdxi=np.empty(shape=(NumNodes,NumIntPoints),dtype=np.float64)
        dNdeta=np.empty(shape=(NumNodes,NumIntPoints),dtype=np.float64)
        for iip in range(NumIntPoints):
            xi=p[0,iip]
            eta=p[1,iip]

            N[0,iip]=1.0/4.0*(1.0-xi)*(1.0-eta)
            N[1,iip]=1.0/4.0*(1.0+xi)*(1.0-eta)
            N[2,iip]=1.0/4.0*(1.0+xi)*(1.0+eta)
            N[3,iip]=1.0/4.0*(1.0-xi)*(1.0+eta)

            dNdxi[0,iip]=1.0/4.0*(-1.0)*(1.0-eta)
            dNdxi[1,iip]=1.0/4.0*(+1.0)*(1.0-eta)
            dNdxi[2,iip]=1.0/4.0*(+1.0)*(1.0+eta)
            dNdxi[3,iip]=1.0/4.0*(-1.0)*(1.0+eta)

            dNdeta[0,iip]=1.0/4.0*(1.0-xi)*(-1.0)
            dNdeta[1,iip]=1.0/4.0*(1.0+xi)*(-1.0)
            dNdeta[2,iip]=1.0/4.0*(1.0+xi)*(+1.0)
            dNdeta[3,iip]=1.0/4.0*(1.0-xi)*(+1.0)

        return N,dNdxi,dNdeta


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


    def UpdateKinematics(self,f):
    
        NumElems=self.elems.shape[1]

        for iElem in range(NumElems):
            ElemNodes=m.elems[:,iElem]
            # NodeCoords=
            # NodeField=
# dNdx: elements,IntPoints,N,coords
# dTdx: dNdx, T
            for iip in range(NumIntPoints):        
                dNdxi=ps.dNdxi[:,iip]
                dNdeta=ps.dNdeta[:,iip]

                dx_dxi=
                dy_dxi=
                dx_deta=
                dy_deta=

                for iNode in range(NumNodes):
                    dN_dxi[0]=dNdxi[iNode]
                    dN_dxi[1]=dNdeta[iNode]
            #     J=[dx_dxi   dy_dxi]
            #        dx_deta  dy_deta]
            #     dNdx=np.linslg.solve(J,dNdxi)


        self.dNdx=dNdx
        self.dfdx=dfdx


#==========================================================================
class material:

    def __init__(self,rho,cp):
        self.rho=rho
        self.cp=cp

    # def InitMesh(self,NumEdgeElems):

    #     coords=np.zeros(shape=((NumEdgeElems+1)*(NumEdgeElems+1),2))
    #     elems=np.zeros(shape=(4,NumEdgeElems*NumEdgeElems),dtype=np.int64)

    #     h=1.0/NumEdgeElems
    #     iNode=-1
    #     iElem=-1
    #     for j in range(NumEdgeElems+1):
    #         for i in range(NumEdgeElems+1):
    #             # print("i={}".format(i))
    #             iNode=iNode+1
    #             coords[iNode,0]=i*h
    #             coords[iNode,1]=j*h
    #             if (i>0 and j>0):
    #                 iElem=iElem+1
    #                 elems[0,iElem]=iNode-(NumEdgeElems+1)-1
    #                 elems[1,iElem]=iNode-(NumEdgeElems+1)
    #                 elems[2,iElem]=iNode
    #                 elems[3,iElem]=iNode-1
    #             # if (16*6<iElem and iElem<16*7):
    #                 # print("iElem={}, ElemNodes={}".format(iElem,elems[:,iElem]))

    #     return coords,elems


#==========================================================================
class fields:

    def __init__(self,InitField,BCField,m,PlotInit=False):
        print("initializing field . . . ",end="")
        self.InitField(InitField,BCField,m.coords)
        if (PlotInit):
            # pause("beginning of PlotInit")
            fig=plt.figure("field initialization",figsize=(9,9))
            # print("field figure number: {}".format(fig.number))
            self.plot(m)
            # pause("after field initialization plot . . .")
            # print("fieldNP1: {}".format(self.fieldNP1))
            plt.show()
        print("done")

    def InitField(self,InitField,BCField,coords):

        NumNodes=coords.shape[0]
        DOFNums=np.empty(shape=(NumNodes),dtype=np.int32)
        fNP1=np.empty(shape=(NumNodes),dtype=np.float64)
        fN=np.empty(shape=(NumNodes),dtype=np.float64)
        fdotNP1=np.empty(shape=(NumNodes),dtype=np.float64)
        fdotN=np.empty(shape=(NumNodes),dtype=np.float64)
        NumDOFs=0
        NumActiveDOFs=0
        NumDirichletDOFs=0
        for iNode in range(NumNodes):
            NumDOFs=NumDOFs+1
            if (coords[iNode,1]<1.0e-05):
                NumDirichletDOFs=NumDirichletDOFs+1
                DOFNums[iNode]=-NumDirichletDOFs
                fNP1[iNode]=BCField
                fN[iNode]=BCField
                fdotNP1[iNode]=0.0
                fdotN[iNode]=0.0
            else:
                NumActiveDOFs=NumActiveDOFs+1
                DOFNums[iNode]=NumActiveDOFs
                fNP1[iNode]=InitField
                fN[iNode]=InitField
                fdotNP1[iNode]=0.0
                fdotN[iNode]=0.0

        self.DOFNums=DOFNums
        self.NumActiveDOFs=NumActiveDOFs
        self.fNP1=fNP1
        self.fN=fN
        self.fdotNP1=fdotNP1
        self.fdotN=fdotN

    def plot(self,m):

        # stupid ass arguments . . .
        NumNodes=m.coords.shape[0]
        NumEdgeNodes=int(sqrt(NumNodes))
        xplot=np.arange(NumEdgeNodes)/(NumEdgeNodes-1)
        yplot=xplot
        xplot1,yplot1=np.meshgrid(xplot,yplot)
        zplot=self.fNP1.reshape(NumEdgeNodes,NumEdgeNodes)
        # print("after call to reshape:")
        # print("fNP1: {}".format(self.fNP1))
        # print("fdotNP1: {}".format(self.fdotNP1))
        # print("zplot: {}".format(zplot))

        cs=plt.contourf(xplot1,yplot1,zplot,100,cmap=cm.jet)
        cb=plt.colorbar(cs)
        cb.ax.set_ylabel('temperature [K]')

        ax=plt.gca()
        ax.set_aspect('equal')

    def update(self,df,m,ti,dt):
        NumNodes=m.coords.shape[0]
        for iNode in range(NumNodes):
            DOFNum=self.DOFNums[iNode]
            if (DOFNum>0):
                self.fNP1[iNode]=self.fNP1[iNode]+df[DOFNum-1]
        self.fdotNP1=ti.update(dt,self)


#==========================================================================
class TimeInt:

    def __init__(self,TimeIntType,alpha=0.0):
        self.type=TimeIntType
        self.InitTimeInt(TimeIntType,alpha)

    def InitTimeInt(self,TimeIntType,alpha):
        if (TimeIntType==1):
            # generalized trapezoidat
            self.alpha=alpha

    def dTdot_dT(self,dt):
        if (self.type==1):
            alpha=self.alpha
            dTdot_dT=1.0/(alpha*dt)

        return dTdot_dT

    def update(self,dt,f):
        if (self.type==1):
            alpha=self.alpha
            fdotNP1=(f.fNP1-f.fN-dt*(1.0-alpha)*f.fdotN)/(alpha*dt)

        return fdotNP1


#==========================================================================
class LinearSystem:

    def __init__(self,f):
        print("initializing linear system . . . ",end="")
        self.RHS=np.empty(shape=(f.NumActiveDOFs),dtype=np.float64)
        self.LHS=np.empty(shape=(f.NumActiveDOFs,f.NumActiveDOFs),dtype=np.float64)
        print("done")
        print("")

    def zero(self):
        self.RHS.fill(0.0)
        self.LHS.fill(0.0)

    def solve(self):
        df=np.linalg.solve(self.LHS,self.RHS)
        return df


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
def UpdateStorageTerm(ps,m,mat,f,ls):

    NumElemNodes=m.elems.shape[0]
    NumElems=m.elems.shape[1]

    ElemNodes=np.empty(shape=(NumElemNodes),dtype=np.int32)
    ElemRHS=np.empty(shape=(NumElemNodes),dtype=np.float64)
    ElemLHS=np.empty(shape=(NumElemNodes,NumElemNodes),dtype=np.float64)
    # print("DOFNums={}".format(f.DOFNums))

    for iElem in range(NumElems):
        print("updating storage term for iElem={} now . . .".format(iElem))
        PrintFlag=False
        if (iElem%7==0):
            PrintFlag=True
        ElemNodes=m.elems[:,iElem]
        NodeCoords=m.coords[ElemNodes,:]
        NodeTemps=f.fNP1[ElemNodes]
        NodeRates=f.fdotNP1[ElemNodes]
        if (PrintFlag):
            print("ElemNodes(iElem={})={}".format(iElem,ElemNodes))
        ElemRHS.fill(0.0)
        ElemLHS.fill(0.0)
        for iip in range(NumIntPoints):
            # int point field and rate
            # IntPointTemp=
            IntPointRate=np.dot(ps.N[:,iip],NodeRates)
            # print("IntPointRate={}".format(IntPointRate))
            # print("dt={}".format(dt))
            dTdot_dT=ti.dTdot_dT(dt)
            print("dTdot_dT={}".format(dTdot_dT))
            Nlocal=ps.N[:,iip]

            # RHS
            ElemRHS=ElemRHS+Nlocal*mat.rho*mat.cp*IntPointRate
            # LHS
            ElemLHS=ElemLHS+np.outer(Nlocal,Nlocal)*mat.rho*mat.cp*dTdot_dT

        if (PrintFlag):
            print("ElemRHS(iElem={})={}".format(iElem,ElemRHS))
            print("ElemLHS(iElem={})={}".format(iElem,ElemLHS))
        ElemDOFs=f.DOFNums[ElemNodes]
        if (PrintFlag):
            print("ElemDOFs(iElem={})={}".format(iElem,ElemDOFs))
        for iDOF in range(NumElemNodes):
            iDOFNum=ElemDOFs[iDOF]
            ls.RHS[iDOFNum-1]=ls.RHS[iDOFNum-1]+ElemRHS[iDOF]
            for jDOF in range(NumElemNodes):
                jDOFNum=ElemDOFs[jDOF]
                ls.LHS[iDOFNum-1,jDOFNum-1]=ls.LHS[iDOFNum-1,jDOFNum-1]+ElemLHS[iDOF,jDOF]
        if (PrintFlag):
            print("")
        # if (iElem==33):
        #     sys.exit()
        print("finished updating storage term for iElem={}".format(iElem))
        print("")


#==========================================================================
def UpdateDiffusionTerm(ps,m,mat,f,ls):

    NumElemNodes=m.elems.shape[0]
    NumElems=m.elems.shape[1]

    ElemNodes=np.empty(shape=(NumElemNodes),dtype=np.int32)
    ElemRHS=np.empty(shape=(NumElemNodes),dtype=np.float64)
    ElemLHS=np.empty(shape=(NumElemNodes,NumElemNodes),dtype=np.float64)
    # print("DOFNums={}".format(f.DOFNums))

    for iElem in range(NumElems):
        PrintFlag=False
        if (iElem%7==0):
            PrintFlag=True
        if (PrintFlag):
            print("updating diffusion term for iElem={} now . . .".format(iElem))
        ElemNodes=m.elems[:,iElem]
        NodeCoords=m.coords[ElemNodes,:]
        NodeTemps=f.fNP1[ElemNodes]
        NodeRates=f.fdotNP1[ElemNodes]
        if (PrintFlag):
            print("ElemNodes(iElem={})={}".format(iElem,ElemNodes))
        ElemRHS.fill(0.0)
        ElemLHS.fill(0.0)
        for iip in range(NumIntPoints):
            # int point field and rate
            # IntPointTemp=
            IntPointRate=np.dot(ps.N[:,iip],NodeRates)
            # print("IntPointRate={}".format(IntPointRate))
            # print("dt={}".format(dt))
            dTdot_dT=ti.dTdot_dT(dt)
            print("dTdot_dT={}".format(dTdot_dT))
            Nlocal=ps.N[:,iip]

            m.UpdateKinematics(f)

            # RHS
            ElemRHS=ElemRHS+dNlocal_dx*mat.k*dTdx
            # LHS
            ElemLHS=ElemLHS+np.outer(dNlocal_dx,dNlocal_dx)*mat.k

        if (PrintFlag):
            print("ElemRHS(iElem={})={}".format(iElem,ElemRHS))
            print("ElemLHS(iElem={})={}".format(iElem,ElemLHS))
        ElemDOFs=f.DOFNums[ElemNodes]
        if (PrintFlag):
            print("ElemDOFs(iElem={})={}".format(iElem,ElemDOFs))
        for iDOF in range(NumElemNodes):
            iDOFNum=ElemDOFs[iDOF]
            ls.RHS[iDOFNum-1]=ls.RHS[iDOFNum-1]+ElemRHS[iDOF]
            for jDOF in range(NumElemNodes):
                jDOFNum=ElemDOFs[jDOF]
                ls.LHS[iDOFNum-1,jDOFNum-1]=ls.LHS[iDOFNum-1,jDOFNum-1]+ElemLHS[iDOF,jDOF]
        if (PrintFlag):
            print("")
        # if (iElem==33):
        #     sys.exit()
        if (PrintFlag):
            print("finished updating diffusion term for iElem={}".format(iElem))
            print("")



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
NumEdgeElems=3
NumIntPoints=4

# 316L
rho=8.0e-06 # kg/mm^3
cp=500.0 # J/kg-K

InitialTemp=303.0
BCTemp=505.0

MaxNumSteps=10^10
MaxNumSteps=2
EndTime=1.0
dt=0.1
TimeIntType=1 # generalized trapezoidal
alpha=1.0

MaxNonlinIters=2



#==========================================================================
# initialization
#==========================================================================
ps=ParametricSpace(NumIntPoints)
m=mesh(NumEdgeElems,PlotInit=False)
mat=material(rho,cp)
f=fields(InitialTemp,BCTemp,m,PlotInit=False)
ti=TimeInt(TimeIntType,alpha)
ls=LinearSystem(f)

fig=plt.figure("primary figure window",figsize=(11,9))
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

        ls.zero()
        UpdateStorageTerm(ps,m,mat,f,ls)
        UpdateDiffusionTerm(ps,m,mat,f,ls)
        # UpdateSourceTerm
        dT=ls.solve()
        print("solution increment={}".format(dT))
        f.update(dT,m,ti,dt)
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



