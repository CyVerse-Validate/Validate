import json
from string import Template

__author__ = "Michael J. Suggs // mjs3607@uncw.edu"


def fastlmm_json(data_name, input_dict):
    """Generates an Agave-compatible JSON for FaST-LMM.

    This JSON is stored in a file named "FaST-LMM_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    :param input_dict: A dictionary generated by the "parse_inputs" function
    """
    """Debug Statement"""
    print "Creating FaST-LMM JSON."

    output = None

    json_text = {
        "name": "FaST-LMM_{}".format(data_name),
        "appId": "FaST-LMM-hpc-2.07",
        "batchQueue": "normal",
        "nodeCount": "1",
        "processorsPerNode": "1",
        "memoryPerNode": "4",
        "maxRunTime": "01:00:00",
        "archive": True,
        "inputs": input_dict,
        "parameters": {
            "output": "full_results"
        },
    }

    # file_name = 'FaST-LMM_{}.json'.format(data_name)
    # with open(file_name, 'w') as f:
    #     json.dumps(json_text, f)
    # print type(json.dumps(json_text))
    # print json.dumps(json_text)
    # return json.dumps(json_text, indent=4, separators=(',', ': '))
    return json_text, output


def ridge_json(data_name, input_dict):
    """Generates an Agave-compatible JSON for RidgePredict.

    This JSON is stored in a file named "RidgePredict_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    :param input_dict: A dictionary generated by the "parse_inputs" function
    """
    """Debug Statement"""
    print "Creating Ridge JSON."

    output = "RidgePredict_{}.ped".format(data_name)

    json_text = {
        "name": "RidgePredict_{}".format(data_name),
        "appId": "RidgePredict-1.1",
        "batchQueue": "normal",
        "nodeCount": "1",
        "processorsPerNode": "1",
        "memoryPerNode": "4",
        "maxRunTime": "01:00:00",
        "archive": True,
        "inputs": {
            "inputPed": "{}".format(input_dict["inputPED"])
        },
        "parameters": {
            "outputPed": "RidgePredict_{}.ped".format(data_name)
        }
    }

    # file_name = 'RidgePredict_{}.json'.format(data_name)
    # with open(file_name, 'w') as f:
    #     json.dumps(json_text, f)
    return json.dumps(json_text, indent=4, separators=(',', ': ')), output


def bayesr_json(data_name, input_dict):
    """Generates an Agave-compatible JSON for BayesR.

    This JSON is stored in a file named "BayesR_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    :param input_dict: A dictionary generated by the "parse_inputs" function
    """
    """Debug Statement"""
    print "Creating BayesR JSON."

    output = "BayesR_{}".format(data_name)

    json_text = {
        "name": "BayesR_{}".format(data_name),
        "appId": "bayesR-2.02",
        "processorsPerNode": 16,
        "requestedTime": "48:00:00",
        "memoryPerNode": 32,
        "nodeCount": 16,
        "archive": True,
        "inputs": input_dict,
        "parameters": {
            "output": "BayesR_{}".format(data_name)
        }
    }

    # file_name = 'BayesR_{}.json'.format(data_name)
    # with open(file_name, 'w') as f:
    #     json.dumps(json_text, f)
    return json.dumps(json_text, indent=4, separators=(',', ': ')), output


def plink_json(data_name):
    """Generates an Agave-compatible JSON for PLINK.

    This JSON is stored in a file named "PLINK_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    """
    """Debug Statement"""
    print "Creating PLINK JSON."

    output = None

    json_text = {
        "name": "Plink_{}".format(data_name),
        "appId": "PLINK-hpc-1.07",
        "processorsPerNode": 16,
        "requestedTime": "01:00:00",
        "memoryPerNode": 32,
        "nodeCount": 1,
        "archive": True,
        "inputs": {
            "file.list": data_name
        },
        "parameters": {
            "out.basename": "success"
        }
    }

    # file_name = 'PLINK_{}.json'.format(data_name)
    # with open(file_name, 'w') as f:
    #     json.dumps(json_text, f)
    return json.dumps(json_text, indent=4, separators=(',', ': ')), output


#TODO make current pipeline implementation work with this input format
def qxpak_json(data_name, input_dict):
    """"""
    """Debug Statement"""
    print "Creating QxPak JSON."

    output = None

    json_text = {
        "name": "QxPak_{}".format(data_name),
        "appId": "qxpak-stampede-5.05",
        "batchQueue": "normal",
        "nodeCount": "1",
        "processorsPerNode": "1",
        "memoryPerNode": "4",
        "maxRunTime": "01:00:00",
        "archive": True,
        "inputs": input_dict,
        "parameter": {
            "output": "test"
        }
    }

    return json.dumps(json_text, indent=4, separators=(',', ': ')), output

#TODO find Gemma JSON?
def gemma_json(data_name, input_dict):
    """Note, GEMMA is SERIAL.

    Inputs:
        inputBED
        inputBIM
        inputFAM
        Phenotype
        Annotation
        MeanGenotype
        RelatednessMatrix
        RelatednessEigenvalue
        RelatednessEigenvector
        Covariate

    Parameters:
        PLINK:  True if PLINK format is being used.
        BIMBAM: True if BIMBAM format is being used.
        gk:     Type of relatedness matrix to calculate. 1 is centered,
            2 is standaredized.
        lmm:    Which frequentist test to use for univariate mixed model analysis.
            1 preforms Wald Test, 2 preforms likelihood ratio test,
            3 preforms score test, and 4 preforms all three tests.
        PhenotypeNumbers:   Numbers indicating which phenotype GEMMA should look
            at for association tests in the multivariate linear mixed model.
        bslmm:  Which type of Bayesian sparse mixed model to fit. 1 runs a
            standard BSLMM, 2 fits a ridge regression/GBLUP, and 3 fits a
            probit BSLMM.
        Type:   Type of anaysis desired. 1 computes relatedness matrix, 2
            computes eigen decomposition of relatedness matrix, 3 computes a
            [univariate/multivariate (?)] mixed model, and 4 is Bayesian spare
            linear mixed model.
        output: Name of the analysis output file.
    """

    """Debug Statement"""
    print "Creating Gemma JSON."

    output = None

    json_text = {
        "name": "GEMMA_{}".format(data_name),
        "appId": "GEMMA-0.94.4",
        "batchQueue": "normal",
        "nodeCount": "1",
        "processorsPerNode": "1",
        "memoryPerNode": "4",
        "maxRunTime": "01:00:00",
        "archive": True,
        "inputs": input_dict,
        "parameters":{
            "PLINK": "0",
            "BIMBAM":"0"
        }
    }

    return json_text, output

def puma_json(data_name, input_dict):
    output = None

    json_text = {
        "name": "puma_{}".format(data_name),
        "appId": "Puma-1.0",
        "batchQueue": "normal",
        "nodeCount": "1",
        "processorsPerNode": "1",
        "memoryPerNode": "4",
        "maxRunTime": "01:00:00",
        "archive": True,
        "inputs": input_dict,
        "parameters":{
            "regression": "LINEAR",
            "penalty":"LASSO",
            "name":data_name
        }
    }

    return json_text, output


def get_json_parameters():
    """

    Returns:


    """
    pass


def make_gwas_json(selected_gwas, dataset_name, inputs):
    """Selected GWAS is a list of binary values with each element representing
        whether a given gwas should be run or not.

    :param selected_gwas: List of GWAS applications to run the data through
    :return:
    """
    if selected_gwas == 0:
        return fastlmm_json(dataset_name, inputs)
    elif selected_gwas == 1:
        return ridge_json(dataset_name, inputs)
    elif selected_gwas == 2:
        return bayesr_json(dataset_name, inputs)
    elif selected_gwas == 3:
        return plink_json(dataset_name)
    elif selected_gwas == 4:
        return qxpak_json(dataset_name, inputs)
    elif selected_gwas == 5:
        return gemma_json(dataset_name, inputs)
    elif selected_gwas == 6:
        return puma_json(dataset_name, inputs)
    else:
        print "No such GWAS"

    # return {
    #     0: fastlmm_json(dataset_name, inputs),
    #     1: ridge_json(dataset_name, inputs),
    #     2: bayesr_json(dataset_name, inputs),
    #     3: plink_json(dataset_name),
    #     4: qxpak_json(dataset_name, inputs),
    #     5: gemma_json(dataset_name, inputs),
    #     6: puma_json(dataset_name, inputs)
    # }[selected_gwas]

    # return switcher.get(selected_gwas)


def make_winnow_json(job_name, gwas_output_folder, ote_file):
    """" Builds JSON files for running Winnow via Agave.

    Built from JSON file here: https://raw.githubusercontent.com/CyVerse-Validate/Stampede-Files/master/json/winnow-job.json
    :param job_name: The name of the GWAS program frmo which the output came.
    :param gwas_output_folder: A folder containing GWAS outputs.
    :param ote_file: OTE known-truth file for truth testing.
    :return:
    """

    json_text = {
        "name" : 'Winnow_{}'.format(job_name),
        "appId" : "Winnow-1.0.3",
        "batchQueue": "normal",
        "nodeCount": "1",
        "processorsPerNode": "1",
        "memoryPerNode": "4",
        "maxRunTime": "01:00:00",
        "archive" : "true",
        "inputs" : {
            "Class" : ote_file,
            "Folder" : gwas_output_folder
        },
        "parameters" : {
            "SNP": "SNP",
            "Score": "P",
            "Filename": "TestResults",
            "beta": "BETA",
            "threshold": "0.01",
            "seper": "whitespace",
            "kttype": "OTE"
        }
    }

    # file_name = 'winnow-{}.json'.format(job_name)
    # with open(file_name, 'w') as f:
    #     json_text = json.dumps(json_text, f)
    return json_text