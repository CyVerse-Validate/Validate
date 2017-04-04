#Brandyn Deffinbaugh
#Simulate Phenotype Value from PED file



import argparse
from sklearn.grid_search import GridSearchCV
from sklearn.linear_model import Ridge
import numpy as np
from sklearn.preprocessing import Imputer
import os


#For all missing Phenotype Values please mark them as NaN

'''
Example PED file rows
1 1 0 0 3 60.57 A A B B A A A A A A B B A
1 2 0 0 3 NaN B B A A A A A A A A A A A
1 3 0 0 3 40.42 B B A A A A 0 0 A A A A B 
'''

class RidgePredict:
        def __init__(self):
                self.arguments = self.checkArgs()
                self.predicted = self.RidgeRegression(self.input,self.output)

                
        def checkArgs(self):
                parser = argparse.ArgumentParser(description = 'PED file')
                parser.add_argument("-i","--input",required = True)
                parser.add_argument("-o","--output",required = True)
                args = parser.parse_args()
                self.input = args.input
                self.output = args.output
                


        #Reads in the data file and parses the input
        def inputParse(self,file):
                pheno = []
                geno = []
                with open(file) as infile:
                        for line in infile:
                                line = line.strip('\n')
                                line = line.replace(' ','\t',6)
                                line = line.split('\t')
                                genoHold = line[6]
                                genoHold = genoHold.split(' ')
                                phenoHold = line[5]
                                geno.append(genoHold)
                                pheno.append(phenoHold)
                return(pheno,geno)


        #Finds the best alpha value out of the given data from the alphas array.
        def alphaOptimization(self,genoMiss,pheno):
                alphas = np.array([0.1,0.01,0.001,0.0001,0.5,0.05,0.2,0.002,0.0002,1,0])
                model = Ridge()
                grid = GridSearchCV(estimator = model, param_grid=dict(alpha=alphas))
                grid.fit(genoMiss,pheno)
                grid.best_score_
                return grid.best_estimator_.alpha


        '''
        Maps all geno to integer values in a range between 0 and the number of unique
        types of geno characters in the given PED file.
        Example:
        1 3 0 0 3 40.42 B B A A A A 0 0 A A A A B
        Unique characters: 3; B A 0
        Would be mapped to 0,1,2
        1 3 0 0 3 40.42 2 2 0 0 0 0 1 1 0 0 0 0 2
        Removes rows that are missing Phenotype Values and then creates a Ridge Regression
        model which then predicts the missing Phenotype Values.
        '''

        def RidgeRegression(self,filename,outputFile):
            pheno,geno = self.inputParse(filename)
            for row in geno:
                        if len(row)%2 !=0:
                                return "Rows are not even."
            maxGeno = max(geno)
            allGeno = list(set(maxGeno))
            encoder = [i for i in range(len(allGeno))]
            lengthGeno = len(geno)
            length = len(geno)
            lenInnerGeno = len(geno[0])
            genoMake = [0 for x in range(len(allGeno))]
            dictionary = dict(zip(allGeno,encoder))
            for i in range(length):
                    for x in range(lenInnerGeno):
                            geno[i][x] = dictionary[geno[i][x]]
            phenoNaN = []
            for i in range(len(pheno)):
                if pheno[i] == 'NaN':
                    phenoNaN.append(i)
            phenoNaN.reverse()
            for i in phenoNaN:
                del pheno[i]
            genoMiss = []
            for i in range(len(geno)):
                if i not in phenoNaN:
                    genoMiss.append(geno[i])
            pheno = [float(i) for i in pheno]    
            alpha = self.alphaOptimization(genoMiss,pheno)
            clf = Ridge(alpha = alpha)
            clf.fit(genoMiss,pheno)
            predicted = clf.predict(geno)
            predicted = np.transpose(predicted)
            np.savetxt(outputFile,np.transpose(predicted))

if __name__ == '__main__':
        R = RidgePredict()
        
