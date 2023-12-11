import os, sys, math
from os import sep, close
import numpy as np
from operator import itemgetter
import openturns as ot

def center_reduce(data_matrix):
    """
    data_matrix : numpy array, columns are different variables
    """
    if(len(data_matrix.shape)>1):
        Nvar = data_matrix.shape[1]
        centeredDataMatrix = np.zeros((len(data_matrix),Nvar))
    else:
        Nvar = 1
        centeredDataMatrix = np.zeros((len(data_matrix)))
    #Center and Reduce
    means = np.zeros((Nvar))
    stds = np.zeros((Nvar))
    for j in range(Nvar):
        if(Nvar>1):
            means[j] = np.mean(data_matrix[:,j])
            stds[j] = np.std(data_matrix[:,j])
            centeredDataMatrix[:,j] = (data_matrix[:,j] - means[j])/stds[j]
        else:
            means[j] = np.mean(data_matrix[:])
            stds[j] = np.std(data_matrix[:])
            centeredDataMatrix[:] = (data_matrix[:] - means[j])/stds[j]

    return means, stds, centeredDataMatrix

def subjective_center_reduce(data_matrix,means,stds):
    """
    data_matrix : numpy array, columns are different variables
    """
    if(len(data_matrix.shape)>1):
        Nvar = data_matrix.shape[1]
        centeredDataMatrix = np.zeros((len(data_matrix),Nvar))
    else:
        Nvar = 1
        centeredDataMatrix = np.zeros((len(data_matrix)))
    for j in range(Nvar):
        if(Nvar>1):
            centeredDataMatrix[:,j] = (data_matrix[:,j] - means[j])/stds[j]
        else:
            centeredDataMatrix[:] = (data_matrix[:] - means[j])/stds[j]

    return centeredDataMatrix


class StatisticalSet():
    """\
    class 
    """
    def __init__(self, inputVariable=np.array([0]), outputVariable=np.array([0])):
        self.helpText = ""
        self.errorMsg = ""

        #
        self.inputVariable = inputVariable      #Numpy array: shape=[setSize,inputDimension]
        self.outputVariable = outputVariable     #Numpy array: shape=[setSize,outputDimension]
        #
        self.setSize = 0             #Integer: Number of individuals
        self.inputDimension = 0      #Integer: Number of input variables
        self.outputDimension = 0     #Integer: Number of output variables; if 2D field: number of points
               
        #Calculate adequate dimensions
        self.setSize = self.inputVariable.shape[0]
        if(len(self.inputVariable.shape)>1):
            self.inputDimension = self.inputVariable.shape[1]
        else:
            self.inputDimension = 1
        #
        if(len(self.outputVariable.shape)>1):
            self.outputDimension = self.outputVariable.shape[1]
        else:
            self.outputDimension = 1
            
    def change_array_values(self,inputVariable, outputVariable):
        self.inputVariable = inputVariable      #Numpy array: shape=[setSize,inputDimension]
        self.outputVariable = outputVariable     #Numpy array: shape=[setSize,outputDimension]
        
    def inherit_from_set_with_indexes(self, indexes):
        """
        Create new set from another one with only elements of given indexes
        
        :param indexes: array of indexes to be stored
        :type indexes: numpy array
        
        :return sub_statistical_set: the created set 
        :type sub_statistical_set: StatisticalSet()
        """
                
        if(self.inputDimension>1):
            inputVariable = self.inputVariable[indexes,:]
        else:
            inputVariable = self.inputVariable[indexes]
            
        if(self.outputDimension>1):
            outputVariable = self.outputVariable[indexes,:]
        else:
            outputVariable = self.outputVariable[indexes]
            
        sub_statistical_set = StatisticalSet(inputVariable=inputVariable, outputVariable=outputVariable)
        
        return sub_statistical_set
             
class ModelAccuracy():
    """\
    class 
    """
    def __init__(self):
        self.helpText = ""
        self.errorMsg = ""

        #Construction phase
        #self.statisticalSet = None      #Inherits from statisticalSet class 
                                        #Set on which the accuracy is calculated (training/prediction set, etc.)
        #
        self.residuals = None           #numpy array, Residuals on statistical set members, length= statisticalSet.setSize
        self.relative_residuals = None   #numpy array, relative residuals on statistical set members, length= statisticalSet.setSize
        self.squared_residuals = None    #numpy array, squared Residuals on statistical set members
        self.absolute_residuals = None   #numpy array, absolute Residuals on statistical set members
        
        #
        self.MAE = None                 #float, positive. Mean Absolute Error (the mean of absolute residuals)
        self.relative_MAE = None         #float, positive. Mean Absolute Error/Mean absolute values of the statistical set output
        self.MSE = None                 #float, positive. Mean Squared Error (the mean of squared residuals)
        self.RMSE = None                #float, positive. Root Mean Squared Error (the root of MSE)
        self.relative_RMSE = None        #float, positive. RMSE/(the mean of squares of the statistical set output values)
        
        #
        self.relativeEmpiricalError= None    #float, positive. MSE/(the variance of the statistical set output values) 
        self.determination_coefficient= None    #float R2 = 1-relativeEmpiricalError 
        
        
        
class OTChaos():
    """\
    class 
    """
    def __init__(self, fullSet=StatisticalSet(), training_size=0, \
                                            test_size=0, prediction_size=0, \
                                            number_of_sets=0):
        self.helpText = ""
        self.errorMsg = ""
        #
        
        self.fullSet = fullSet
        
        if(test_size+prediction_size+training_size>self.fullSet.setSize):
            print("WARNING: the given sizes exceeds the full size. Only the fullSet attribute is set in the OtChaos object.")
        elif(training_size==0): 
            print("WARNING: training size was not given. Only the fullSet attribute is set in the OtChaos object.")
        else:
            if(test_size+prediction_size+training_size>self.fullSet.setSize):
                print("WARNING: the full set is not covered.")
            #Set on which the model is being fitted
            training_indexes = np.arange(training_size)
            self.trainingSet = self.fullSet.inherit_from_set_with_indexes(training_indexes)
            #        
            if(test_size==0 and prediction_size==0):
                if(number_of_sets!=1):
                    print("WARNING: test_size and prediction_size were not given. Only the fullSet and trainingSet attributes are set.")
                    print("       : try again or adjust number_of_sets argument to ignore message.")
                
            elif(prediction_size==0):
                test_indexes = np.arange(test_size) + training_size
                self.testSet = self.fullSet.inherit_from_set_with_indexes(test_indexes)
                self.predictionSet = self.fullSet.inherit_from_set_with_indexes(test_indexes)
                if(number_of_sets!=2):
                    print("WARNING: prediction_size was not given. Attributes predictionSet and testSet are equal.")
                    print("       : try again or adjust number_of_sets argument to ignore message.")
                    
            elif(test_size==0):    
                prediction_indexes = np.arange(prediction_size) + training_size
                self.predictionSet = self.fullSet.inherit_from_set_with_indexes(prediction_indexes)
                self.testSet = self.fullSet.inherit_from_set_with_indexes(prediction_indexes)
                if(number_of_sets!=2):
                    print("WARNING: test_size was not given. Attributes predictionSet and testSet are equal.")
                    print("       : try again or adjust number_of_sets argument to ignore message.")
            else:
                test_indexes = np.arange(test_size) + training_size
                self.testSet = self.fullSet.inherit_from_set_with_indexes(test_indexes)
                prediction_indexes = np.arange(prediction_size) + test_size + training_size
                self.predictionSet = self.fullSet.inherit_from_set_with_indexes(prediction_indexes)
               
                    
        #
        self.trainingAccuracy = ModelAccuracy()    #Accuracy indicators are calculated on trainingSet
        self.testAccuracy = ModelAccuracy()        #Accuracy indicators are calculated on testSet
        self.predictionAccuracy = ModelAccuracy()  #Accuracy indicators are calculated on predictionSet
        
        #Calibration choices
        self.polynomialDegree= None     #integer
        self.marginals = None           #List of strings: marginals for input. Choices: Gaussian, Uniform, Kernel
        self.copula = None              #Modeling the dependencies in the input's space. Choice: independent, guassian
        #Algorithm choices
        self.truncationMethod= None     #string, method for the expansion trunction. Choices: Linear, Hyperbolic
        self.fittingStrategy = None     #string, method for expansion's coefficients calibration. Choices: OLS, LARS
        
        #output result, openturns objects
        self.metamodel = None 
        self.polynomialChaosResult = None 
        self.enumerateFunction = None
        
    def change_set_values(self,trainingSet, testSet, predictionSet):
        self.trainingSet = trainingSet      #Numpy array: shape=[setSize,inputDimension]
        self.testSet = testSet     #Numpy array: shape=[setSize,outputDimension]        
        self.predictionSet = predictionSet     #Numpy array: shape=[setSize,outputDimension]
        
    def center_reduce_inputs_with_training_set(self):
        """
        description
        
        :param self: explique
        :return bla: 
        """
        #Center and reduce trainingSet
        means, stds, self.trainingSet.inputVariable = center_reduce(self.trainingSet.inputVariable)
        
        #Center and reduce predictionSet and testSet with means and variance from trainingSet
        self.fullSet.inputVariable = subjective_center_reduce(self.fullSet.inputVariable, means, stds)
        self.predictionSet.inputVariable = subjective_center_reduce(self.predictionSet.inputVariable, means, stds)
        self.testSet.inputVariable = subjective_center_reduce(self.testSet.inputVariable, means, stds)
        
        return means, stds
        
    def center_reduce_output(self):
        means, stds, self.trainingSet.outputVariable = center_reduce(self.trainingSet.outputVariable)
        self.predictionSet.outputVariable = subjective_center_reduce(self.predictionSet.outputVariable, means, stds)
        self.testSet.outputVariable = subjective_center_reduce(self.testSet.outputVariable, means, stds)
        
        
    def construct_PCE_ot(self):
        trainingGerm = self.trainingSet.inputVariable
        trainingOutput = self.trainingSet.outputVariable
        marginals = self.marginals
        degree = self.polynomialDegree
        ##########################INPUTS##########################
        ##########################################################
        Nt = len(trainingGerm)
        if(len(trainingGerm.shape) > 1):
            Nvar = trainingGerm.shape[1]
        else:
            Nvar = 1
        
        #Define Sample
        outputSample = ot.Sample(Nt,1)
        for i in range(Nt):
            outputSample[i,0] = trainingOutput[i]
    
        #Define Collection and PDFs
        polyColl = ot.PolynomialFamilyCollection(Nvar)    
        collection = ot.DistributionCollection(Nvar)
        marginal = {}
        UncorrelatedInputSample = ot.Sample(Nt,Nvar)
    
        if(Nvar>1):
            for i in range(Nvar):
                varSample = ot.Sample(Nt,1)
                for j in range(Nt):
                    varSample[j,0] = trainingGerm[j,i]
                    UncorrelatedInputSample[j,i] = trainingGerm[j,i]
                minValue = varSample.getMin()[0]
                maxValue = varSample.getMax()[0]
                if(self.marginals[i]=="gaussian" or self.marginals[i]=="normal"):
                    marginal[i] = ot.NormalFactory().build(varSample)
                elif(self.marginals[i]=="uniform"):
                    marginal[i] = ot.Uniform(minValue-minValue/100.,maxValue+maxValue/100.)
                elif(self.marginals[i]=="kernel"):
                    marginal[i] = ot.KernelSmoothing().build(varSample)
                else:
                    print("WARNING: couldn't find distribution '"+str(self.marginals[i])+"', applied kernel smoothing instead")
                    marginal[i] = ot.KernelSmoothing().build(varSample)
                    
                collection[i] = ot.Distribution(marginal[i])
        else:
            varSample = ot.Sample(Nt,1)
            for j in range(Nt):
                varSample[j,0] = trainingGerm[j]
                UncorrelatedInputSample[j,0] = trainingGerm[j]
            minValue = varSample.getMin()[0]
            maxValue = varSample.getMax()[0]
            if(self.marginals[i]=="gaussian" or self.marginals[i]=="normal"):  
                marginal[0] = ot.NormalFactory().build(varSample)
            elif(self.marginals[i]=="uniform"):
                marginal[0] = ot.Uniform(minValue-minValue/100.,maxValue+maxValue/100.)
            elif(self.marginals[i]=="kernel"):
                marginal[0] = ot.KernelSmoothing().build(varSample)
            else:
                print("WARNING: couldn't find distribution '"+str(self.marginals[i])+"', applied kernel smoothing instead")
                marginal[0] = ot.KernelSmoothing().build(varSample)
            collection[0] = ot.Distribution(marginal[0])
        
        if(self.copula=="independent"):
            copula = ot.IndependentCopula(Nvar)
        elif(self.copula=="gaussian" or self.copula=="normal"):
            inputSample = ot.Sample(trainingGerm)
            copula = ot.NormalCopulaFactory().build(inputSample)
        else:
            print("WARNING: couldn't find copula '"+str(self.copula)+"', applied independent copula instead")
            copula = ot.IndependentCopula(Nvar)
            
        #UncorrelatedInputDistribution = ot.ComposedDistribution(collection,ot.Copula(copula))
        UncorrelatedInputDistribution = ot.ComposedDistribution(collection,copula)
    
        #Calcul des polynomes du chaos
        for v in range(0,Nvar):
            marginalv=UncorrelatedInputDistribution.getMarginal(v)            
            if(self.marginals[i]=="kernel"):
                #Works with arbitrary PDF
                basisAlgorithm = ot.AdaptiveStieltjesAlgorithm(marginalv)
                polyColl[v] = ot.StandardDistributionPolynomialFactory(basisAlgorithm)
            else:
                #Works with standard PDF: gaussian, uniform, ..
                polyColl[v] = ot.StandardDistributionPolynomialFactory(marginalv)
                
        # Definition de la numerotation des coefficients des polynomes du chaos
        enumerateFunction = ot.LinearEnumerateFunction(Nvar)
        #enumerateFunction = HyperbolicAnisotropicEnumerateFunction(Nvar,0.4)
        # Creation de la base des polynomes multivaries en fonction de la numerotation
        #                     et des bases desiree
        multivariateBasis = ot.OrthogonalProductPolynomialFactory(polyColl,enumerateFunction)
        # Number of PC terms
        P = enumerateFunction.getStrataCumulatedCardinal(degree)
        #Troncature
        adaptativeStrategy = ot.FixedStrategy(multivariateBasis,P)
        #Evaluation Strategy : LARS
        basisSequenceFactory = ot.LARS()
        fittingAlgorithm = ot.CorrectedLeaveOneOut()
        approximationAlgorithm = ot.LeastSquaresMetaModelSelectionFactory(basisSequenceFactory, fittingAlgorithm)
        
        #Methode d'approximation des coefficients de la decomposition en PC
        #weights = NumericalPoint(UncorrelatedInputSample.getSize(),1.)
        projectionStrategy = ot.LeastSquaresStrategy(UncorrelatedInputSample, outputSample, approximationAlgorithm)
        #projectionStrategy = LeastSquaresStrategy(UncorrelatedInputSample, weights, outputSample)
        #algo = FunctionalChaosAlgorithm(UncorrelatedInputSample, weights, outputSample, UncorrelatedInputDistribution, adaptativeStrategy, projectionStrategy)
        algo = ot.FunctionalChaosAlgorithm(UncorrelatedInputSample, outputSample, UncorrelatedInputDistribution, adaptativeStrategy, projectionStrategy)
        algo.run()
        polynomialChaosResult = algo.getResult()
        metamodel = polynomialChaosResult.getMetaModel()
    
        self.metamodel = metamodel
        self.polynomialChaosResult = polynomialChaosResult
        self.enumerateFunction = enumerateFunction
        
        
         
    def estimate_on_statistical_set(self,statisticalSet):
        """
        Mean Squared Error
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
        
        return metaEvaluation
    
    def MSE_on_statistical_set(self,statisticalSet):
        """
        Mean Squared Error
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
        MSE = 0. #Sum of Squared Difference
        for i in range(Np):
            MSE = MSE + (metaEvaluation[i] - predictionOutput[i])**2
        MSE = MSE/Np
        
        return MSE
    
    def MAE_on_statistical_set(self,statisticalSet):
        """
        Mean Absolute Error
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
        MAE = 0. #Sum of Squared Difference
        for i in range(Np):
            MAE = MAE + np.abs(metaEvaluation[i] - predictionOutput[i])
        MAE = MAE/Np
        
        return MAE
    
    def relative_MAE_on_statistical_set(self,statisticalSet):
        """
        Mean Absolute Error
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
        MAE = 0. #Mean Absolute Err
        SA = 0. #Sum of Absolutes
        for i in range(Np):
            MAE = MAE + np.abs(metaEvaluation[i] - predictionOutput[i])
            SA = SA + np.abs(predictionOutput[i])
        
        relative_MAE = MAE/SA
        
        return relative_MAE
    
    def RMSE_on_statistical_set(self,statisticalSet):
        """
        Root Mean Squared Error
        """
        
        MSE = self.MSE_on_statistical_set(statisticalSet)
        RMSE = np.sqrt(MSE)
        
        return RMSE
    
    def relative_RMSE_on_statistical_set(self,statisticalSet):
        """
        Relative Root Mean Squared Error
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
        RMSE = 0. #Sum of Squared Difference
        SS = 0. #Sum of Squares
        for i in range(Np):
            RMSE = RMSE + (metaEvaluation[i] - predictionOutput[i])**2
            SS = SS + (predictionOutput[i])**2
        
        relative_RMSE = np.sqrt(RMSE/SS)
        
        return relative_RMSE
    
    def relative_empirical_error_on_statistical_set(self,statisticalSet):
        """ compute relative empirical error as in Blatman 2009
        
        :param statisticalSet: set on which the model error is evaluated
        :type statisticalSet: StatisticalSet() class 
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
        SSD = 0. #Sum of Squared Difference
        variance = 0.
        for i in range(Np):
            SSD = SSD + (metaEvaluation[i] - predictionOutput[i])**2
            variance = variance + (predictionOutput[i] - np.mean(predictionOutput))**2
        SSD = SSD/Np
        variance = variance/(Np-1)
        empiricalErr = SSD/variance #relative
        return empiricalErr

    def determination_coefficient_on_statistical_set(self,statisticalSet):
        return 1 - self.relative_empirical_error_on_statistical_set(statisticalSet)


    def residuals_on_statistical_set(self,statisticalSet):
        """
        residuals
        
        :return residuals: table of residuals on statistical set
        :type residuals: numpy array
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
                
        residuals = np.zeros((Np))
        for i in range(Np):
            residuals[i] = predictionOutput[i]-metaEvaluation[i]
                
        return residuals
    
    def relative_residuals_on_statistical_set(self,statisticalSet):
        """
        residuals
        
        :return relative residuals: table of relative residuals on statistical set
        :type residuals: numpy array
        """
        metamodel = self.metamodel
        
        predictionGerm = statisticalSet.inputVariable
        predictionOutput = statisticalSet.outputVariable
        
        Np = len(predictionGerm)
        metaEvaluation = np.zeros((Np))
        for i in range(Np):
            if(len(predictionGerm.shape)>1):
                metaEvaluation[i] =  metamodel(predictionGerm)[i,0]
            else:
                metaEvaluation[i] =  metamodel([predictionGerm[i]])[0]
                
        relative_residuals = np.zeros((Np))
        for i in range(Np):
            relative_residuals[i] = (predictionOutput[i]-metaEvaluation[i])/np.abs(predictionOutput[i])
                
        return relative_residuals
    
    def absolute_residuals_on_statistical_set(self,statisticalSet):
        """
        Absolute Residuals
        
        :return absolute_residuals: table of absolute residuals on statistical set
        :type absolute_residuals: numpy array
        """
        absolute_residuals = self.residuals_on_statistical_set(statisticalSet)
        
        absolute_residuals = np.abs(absolute_residuals)
        
        return absolute_residuals

    def squared_residuals_on_statistical_set(self,statisticalSet):
        """
        Squared Residuals
        
        :return squared_residuals: table of absolute residuals on statistical set
        :type squared_residuals: numpy array
        """
        squared_residuals = self.residuals_on_statistical_set(statisticalSet)
        
        squared_residuals = squared_residuals**2
        
        return squared_residuals


    def compute_relative_empirical_error(self):
        self.trainingAccuracy.relativeEmpiricalError = self.relative_empirical_error_on_statistical_set(self.trainingSet)
        self.testAccuracy.relativeEmpiricalError = self.relative_empirical_error_on_statistical_set(self.testSet)
        self.predictionAccuracy.relativeEmpiricalError = self.relative_empirical_error_on_statistical_set(self.predictionSet)
    
    def compute_determination_coefficient(self):
        self.trainingAccuracy.determination_coefficient = self.determination_coefficient_on_statistical_set(self.trainingSet)
        self.testAccuracy.determination_coefficient = self.determination_coefficient_on_statistical_set(self.testSet)
        self.predictionAccuracy.determination_coefficient = self.determination_coefficient_on_statistical_set(self.predictionSet)
        
    def compute_MAE(self):
        self.trainingAccuracy.MAE = self.MAE_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.MAE = self.MAE_on_statistical_set(self.predictionSet)
        self.testAccuracy.MAE = self.MAE_on_statistical_set(self.testSet)
        
    def compute_relative_MAE(self):
        self.trainingAccuracy.relative_MAE = self.relative_MAE_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.relative_MAE = self.relative_MAE_on_statistical_set(self.predictionSet)
        self.testAccuracy.relative_MAE = self.relative_MAE_on_statistical_set(self.testSet)
        
    def compute_MSE(self):
        self.trainingAccuracy.MSE = self.MSE_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.MSE = self.MSE_on_statistical_set(self.predictionSet)
        self.testAccuracy.MSE = self.MSE_on_statistical_set(self.testSet)
        
    def compute_RMSE(self):
        self.trainingAccuracy.RMSE = self.RMSE_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.RMSE = self.RMSE_on_statistical_set(self.predictionSet)
        self.testAccuracy.RMSE = self.RMSE_on_statistical_set(self.testSet)
        
    def compute_relative_RMSE(self):
        self.trainingAccuracy.relative_RMSE = self.relative_RMSE_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.relative_RMSE = self.relative_RMSE_on_statistical_set(self.predictionSet)
        self.testAccuracy.relative_RMSE = self.relative_RMSE_on_statistical_set(self.testSet)
    
    
    def compute_residuals(self):
        self.trainingAccuracy.residuals = self.residuals_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.residuals = self.residuals_on_statistical_set(self.predictionSet)
        self.testAccuracy.residuals = self.residuals_on_statistical_set(self.testSet)

    def compute_relative_residuals(self):
        self.trainingAccuracy.relative_residuals = self.relative_residuals_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.relative_residuals = self.relative_residuals_on_statistical_set(self.predictionSet)
        self.testAccuracy.relative_residuals = self.relative_residuals_on_statistical_set(self.testSet)
        
    def compute_squared_residuals(self):
        self.trainingAccuracy.squared_residuals = self.squared_residuals_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.squared_residuals = self.squared_residuals_on_statistical_set(self.predictionSet)
        self.testAccuracy.squared_residuals = self.squared_residuals_on_statistical_set(self.testSet)
        
    def compute_absolute_residuals(self):
        self.trainingAccuracy.absolute_residuals = self.absolute_residuals_on_statistical_set(self.trainingSet)
        self.predictionAccuracy.absolute_residuals = self.absolute_residuals_on_statistical_set(self.predictionSet)
        self.testAccuracy.absolute_residuals = self.absolute_residuals_on_statistical_set(self.testSet)
        
        
    def compute_all_errors(self):
        self.compute_relative_empirical_error()
        #self.compute_determination_coefficient()
        self.compute_MAE()
        self.compute_relative_MAE()
        self.compute_MSE()
        self.compute_RMSE()
        self.compute_relative_RMSE()
        
    def compute_all_residuals(self):
        self.compute_residuals()
        self.compute_relative_residuals()
        self.compute_absolute_residuals()
        self.compute_squared_residuals()
        
                
    def get_error(self, criteria='Relative Empirical Error'):
        """
        Returns the value of the calculated error for current model
        
        :param criteria: the error type 
        :type criteria: string
        
        :return the error
        :type float
        """
        #Choice of criteria
        if(criteria=='Relative Empirical Error'):
            return self.trainingAccuracy.relativeEmpiricalError, \
                self.testAccuracy.relativeEmpiricalError, \
                self.predictionAccuracy.relativeEmpiricalError 
        elif(criteria=='MSE'):
            return self.trainingAccuracy.MSE, \
                self.testAccuracy.MSE, \
                self.predictionAccuracy.MSE
        elif(criteria=='RMSE'):
            return self.trainingAccuracy.RMSE, self.testAccuracy.RMSE, \
                    self.predictionAccuracy.RMSE
        elif(criteria=='Relative RMSE'): 
            return self.trainingAccuracy.relative_RMSE, \
                self.testAccuracy.relative_RMSE, \
                self.predictionAccuracy.relative_RMSE
        elif(criteria=='MAE'): 
            return self.trainingAccuracy.MAE, \
                self.testAccuracy.MAE, \
                self.predictionAccuracy.MAE
        elif(criteria=='Relative MAE'):
            return self.trainingAccuracy.relative_MAE, \
                self.testAccuracy.relative_MAE, \
                self.predictionAccuracy.relative_MAE
        elif(criteria=='Determination Coefficient'):
            return self.trainingAccuracy.determination_coefficient, \
                self.testAccuracy.determination_coefficient, \
                self.predictionAccuracy.determination_coefficient
            
            
    def find_optimum(self,maximal_degree, criteria='Relative Empirical Error', \
                     decrease_rate=0., decrease_margin=0., \
                     stop_if_error_growth=False, \
                     compute_errors="", compute_residuals=""):
        """
        iterates on polynomial degrees to find optimum
        
        :param maximal_degree: maximum for loop
        :param decrease_rate: model of degree P is better if
                 error(degree P) < optimal_error(degree P') - 
                             decrease_rate*optimal_error(degree P')
        :param decrease_margin: model of degree P is better if
                 error(degree P) < optimal_error(degree P') - decrease_margin
                              
        :param stop_if_error_growth: if  error(degree P) > optimal_error(degree P')
                                        then optimum is P-1
        
        :type maximal_degree: integer
        :type decrease_rate: float 
        :type stop_if_error_growth: boolean
        
        :return optimal_degree: optimal degree that minimizes the criteria error
        :return prediction_error_table: table of errors on the prediction set, 
                            for each degree (rows) in the following order (columns): 
                            relative empirical error, determination coefficient, 
                            MSE, RMSE, relative_RMSE, MAE, relative_MAE
        :return training_error_table: same as prediction_error_table, but 
                            calculated on the training set
                            
        :return error_name_list : list of error names 
        
        :type optimal_degree: integer
        :type error_table: two-dimensional numpy array
        """
        minimize_error = 1000.
        maximize_error = -1000
        
        prediction_error_table = np.zeros((maximal_degree+1, 7))  
        training_error_table = np.zeros((maximal_degree+1, 7))  
        error_name_list = ['Relative Empirical Error', 'Determination Coefficient', 
                            'MSE', 'RMSE', 'Relative RMSE', 'MAE', 'Relative MAE']
                
        for PC_deg in range(maximal_degree+1):
            self.polynomialDegree= PC_deg
            self.construct_PCE_ot()
            if(compute_errors=="all"):
                self.compute_all_errors()
            elif(compute_errors=="Relative Empirical Error"):
                self.compute_relative_empirical_error()
            if(compute_residuals=="all"):
                self.compute_all_residuals()
                        
            prediction_error_table[PC_deg,:] = np.array([self.predictionAccuracy.relativeEmpiricalError, \
                                             self.predictionAccuracy.determination_coefficient, \
                                             self.predictionAccuracy.MSE, self.predictionAccuracy.RMSE, \
                                             self.predictionAccuracy.relative_RMSE, self.predictionAccuracy.MAE, \
                                             self.predictionAccuracy.relative_MAE])
            
            training_error_table[PC_deg,:] = np.array([self.trainingAccuracy.relativeEmpiricalError, \
                                             self.trainingAccuracy.determination_coefficient, \
                                             self.trainingAccuracy.MSE, self.trainingAccuracy.RMSE, \
                                             self.trainingAccuracy.relative_RMSE, self.trainingAccuracy.MAE, \
                                             self.trainingAccuracy.relative_MAE])
            
            #Choice of criteria
            if(criteria=='Relative Empirical Error'):
                action='minimize'
                criteria_error = self.predictionAccuracy.relativeEmpiricalError
            if(criteria=='MSE'):
                action='minimize'
                criteria_error = self.predictionAccuracy.MSE
            if(criteria=='RMSE'):
                action='minimize'
                criteria_error = self.predictionAccuracy.RMSE
            if(criteria=='Relative RMSE'):
                action='minimize'
                criteria_error = self.predictionAccuracy.relative_RMSE 
            if(criteria=='MAE'):
                action='minimize'
                criteria_error = self.predictionAccuracy.MAE 
            if(criteria=='Relative MAE'):
                action='minimize'
                criteria_error = self.predictionAccuracy.relative_MAE
                
            elif(criteria=='Determination Coefficient'):
                action='maximize'
                criteria_error = self.predictionAccuracy.determination_coefficient
                
            #is optimum?
            if(action=='minimize'):
                if(criteria_error<minimize_error - decrease_rate*minimize_error and \
                   criteria_error<minimize_error - decrease_margin):
                    optimal_degree = PC_deg
                    minimize_error = criteria_error
                if(stop_if_error_growth):
                    if(criteria_error > minimize_error):
                        break;
                           
            elif(action=='maximize'):
                if(criteria_error>maximize_error + decrease_rate*maximize_error and \
                   criteria_error<maximize_error + decrease_margin):
                    optimal_degree = PC_deg
                    maximize_error = criteria_error
                if(stop_if_error_growth):
                    if(criteria_error < maximize_error):
                        break;
                           
            del(self.metamodel)
            del(self.polynomialChaosResult)
            del(self.enumerateFunction)
        
        return optimal_degree, training_error_table,prediction_error_table,error_name_list
    
    
    def string_model(self, variable_name_list):
        """
        returns the string function with adequate variables names
        
        :param variable_name_list: list of variable names (acronyms, etc.)
        :type variable_name_list: string
        
        :return model_string: the written PCE with adequate variables names
        :type model_string: string
        """
        
        variable_number = len(variable_name_list)
        model_string = str(self.polynomialChaosResult.getComposedMetaModel())
        
        for v in range(variable_number-1,-1,-1):
            model_string = model_string.replace('x'+str(v),variable_name_list[v])
    
        return model_string
