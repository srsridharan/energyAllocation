"""
This is a program for performing the aggregator house mapping
11/25/2012
Author: s^2

v1
"""


import cplex as cplex
import numpy as np
from cplex.exceptions import CplexError
import scipy as scpy
from scipy import sparse
import scipy.linalg as scipylin
import sys
import pdb



import cPickle as cpk 


def mainFn():
    """
    read and set values required prior to calling the allocator routine
    """
    constNumAgg = np.int64(30)
    constNumHouses = np.int64(30000)
    numHoursInDay = 24
    displayStatus =0 
    hMat = np.eye(constNumHouses,numHoursInDay)
        
    import houseElectricConsumption as hEConsumption
    #import tables

    #h5file = tables.openFile("powerData.h5", mode = "w", title = "Test file")
    #group = h5file.createGroup("/", 'detector', 'Detector information')
    #table = h5file.createTable(group, 'readout', Particle, "Readout example")

    output = open('powerDataBinary30_1.pkl', 'wb')

    

    ## Now obtain the house data
    for k in xrange(constNumHouses):
        hMat[k,:] = hEConsumption.houseElectricConsumption() 
   
    # now read the total power production data
    totalPowerDataFile = open('./totalPowerFromAllAgg.txt', 'r'); d1= totalPowerDataFile.read(); totalPowerDataFile.close()
    totalPowerData= d1.split('\r') # split the lines using the location of the carraige return
    aggregatorMat = map(lambda x: np.mean(constNumHouses*hMat[20:30,x]), range(24)) # this creates an artificial aggregator power
        # array which is used for the simulation below. It computes a feasible aggregator supply as (number of houses * power reqd for one house)
    gammaMatDict = dict()
    cpk.dump([constNumAgg,constNumHouses,numHoursInDay],output)
    # now for the first time instant create hMeanVec
    np.random.seed(18)
    for stepCtr in xrange(numHoursInDay):
        hMeanVec = np.reshape(hMat[:,stepCtr],(constNumHouses,1))
        aArray =  np.reshape([1.0*aggregatorMat[stepCtr]/constNumAgg]*constNumAgg,(constNumAgg,1))
        try:
            flex_weights= np.abs(np.random.beta(2,2,(constNumHouses,1)))
            gammaMat =  aggregatorAgg(constNumAgg,constNumHouses  , hMeanVec = hMeanVec,aggMeanVec = aArray,flex_weights = flex_weights)

            print "step number: ", stepCtr
            if displayStatus:
                print gammaMat
                #pdb.set_trace()
            cpk.dump(hMeanVec,output)
            cpk.dump(flex_weights,output)
            cpk.dump(aArray,output)
            cpk.dump(gammaMat, output)
            output.flush()
        except:
            cpk.dump(None,output)
            cpk.dump(None,output)
            cpk.dump(None,output)
            output.flush()
            print 'failed'


    output.close()

def LinearTermInObjectiveFn(constNumAgg,constNumHouses,hMeanVec,aggMeanVec):
    """
    This function generates linear terms in the LP objective 
    Input:
        @constNumAgg:
        @constNumHouses
        @hMeanVec: vector of average house power demands for each house
        @aggMeanVec: vector of aggregator demands 
    Output:
        linearTerm: generated linear terms for the objective function

    """

    sparseAggMeanVec = sparse.csc_matrix(aggMeanVec)
    linTH = -2*sparseAggMeanVec.T*H
    linTH = linTH.tocsc()
    linData = linTH.data
    linObjRowInd,linColInd =  linTH.nonzero()

    linearTerm = zip(linColInd,linData)
        
    return linearTerm

def aggregatorAgg(constNumAgg,constNumHouses , hMeanVec ,aggMeanVec,flex_weights):
    """
    This is the function where the allocation is performed
    Input: 
        hMeanVec : vector of house consumption values
        aggMeanVec: vector of aggregator mean supply values
    Output:
        gammaMat: sparse matrix of aggregator allocation for houses. gammaMat(i,j) = 1 iff house j is allocated to aggregator i

    
    The problem involves an unknown gamma matrix:
        Gamma = [gamma_1, gamma_2... gamma_m
                gamma_{m+1}, ...    gamma_{2m},

                .....
                ...
                gamma{(n-1)*m+1}.......  gamma_{n*m}
                ]
    The vectorized form of this is taken to be:
        x = [gamma_1; 
            gamma_2; 
            .... 
            gamma_{m*n}
            ]

        NOTE: the vectorization is performed row wise instead of columnwise. this is done for ease of formulating the constraints on the 
        sum of the Gamma matrix over each row being equal to 1.

    The optimization problem being solved has the form:
        minimize_{x} {Sigma*x}    where x is the \mathop{vec}(\Gamma^T) (NOTE: also possible to modify this by adding 
            (upperbound*agg Demand - houseDemand*gamma)+ (houseDemand*gamma - lowerbound*aggDemand) with possible weightings along each dimension  )
        constraints:1.  The constraint is that the sum of all rows of the Gamma matrix must be  = [1,1,1...]
                    2. LowerBound*AggregatorDemand<=HouseDemandMatrix*vec(GammaMatrix) <= UpperBound*AggregatorDemand



    """

    my_prob = cplex.Cplex()
    flex_weightsTranspose = flex_weights.T
    agg_Sigma =  np.ones((constNumAgg,1))
    UB = 1.05
    LB = .95
    my_prob.objective.set_sense(my_prob.objective.sense.minimize)
    my_prob.set_problem_type(type=  my_prob.problem_type.MILP)

    my_prob.variables.add(names = [str(i) for i in range(constNumHouses*constNumAgg)])
    my_prob.variables.set_types(zip(range(constNumHouses*constNumAgg),[my_prob.variables.type.binary]*(constNumHouses*constNumAgg)))
    linearT = sparse.csc_matrix(np.reshape((agg_Sigma*flex_weightsTranspose).T,(1,len(agg_Sigma)*len(flex_weights))).squeeze())
    linData = linearT.data
    linRInd,linCInd = linearT.nonzero()
    linearTerm = zip(linCInd,linData)
    my_prob.objective.set_linear(linearTerm) 

    #now set the linear constraints 
    my_sense = ["E"]*constNumHouses+ ["R"]*constNumAgg


    rhsVal= LB*aggMeanVec
    rhsVal =  rhsVal.squeeze().tolist()

    sumConstraint = [np.float(1)]*constNumHouses+rhsVal
    range_valuesVal = [0.0]*constNumHouses + ((UB-LB)*aggMeanVec).squeeze().tolist()
    
    #generate the matrix for the equality constraints on the fact that the rows in Gamma sum to [1,1,...1]

    m1 = sparse.eye(constNumHouses,constNumHouses,dtype = np.float)
    eqMat1  = reduce(lambda y1,y2:sparse.hstack([y1,y2]), [m1]*constNumAgg)
    
    
    H =reduce(lambda x,y: sparse.bmat( [[x,None],[None,y]] ), [hMeanVec.T]*constNumAgg)
    eqMat = sparse.vstack([eqMat1,H]) 
    

    sumValues = eqMat.data
    linEqRowInd, linEqColInd = eqMat.nonzero()
    eqTuple = zip(linEqRowInd,linEqColInd,sumValues)
    my_prob.linear_constraints.add(senses = my_sense,rhs = sumConstraint, range_values = range_valuesVal)
    my_prob.linear_constraints.set_coefficients(eqTuple) 
    
    # now in order to make the count on the constraint indices start from the correct number:

    #my_prob.parameters.mip.strategy.probe.set(3)
    #my_prob.parameters.mip.cuts.mircut.set(1)
    #my_prob.parameters.mip.cuts.mcfcut.set(1)  
    #my_prob.parameters.mip.tolerances.mipgap.set(.03) 
    my_prob.parameters.mip.display.set(0) 
    my_prob.solve()
    # solution.get_status() returns an integer code
    #print "Solution status = " , my_prob.solution.get_status(), ":",
    # the following line prints the corresponding string

    gammavec = my_prob.solution.get_values()
    gammaMat = np.reshape(gammavec,(constNumAgg,constNumHouses))
    return gammaMat
       
if __name__ == "__main__":
    mainFn()

