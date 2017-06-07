import json

__author__ = "Michael J. Suggs // mjs3607@uncw.edu"


def fastlmm_json(data_name, input_dict):
    """Generates an Agave-compatible JSON for FaST-LMM.

    This JSON is stored in a file named "FaST-LMM_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    :param input_dict: A dictionary generated by the "parse_inputs" function
    """
    json_text = {
        "jobName": "FaST-LMM_{}".format(data_name),
        "softwareName": "FaST-LMM-hpc-2.07u1",
        "requestedTime": "02:00:00",
        "archive": True,
            "inputs": input_dict,
        "parameters": {
            "output": "full_results"
        }
    }

    file_name = 'FaST-LMM_{}.json'.format(data_name)
    with open(file_name, 'w') as f:
        json.dumps(json_text, f)
    return file_name


def ridge_json(data_name, input_dict):
    """Generates an Agave-compatible JSON for RidgePredict.

    This JSON is stored in a file named "RidgePredict_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    :param input_dict: A dictionary generated by the "parse_inputs" function
    """
    json_text = {
        "jobName": "RidgePredict_{}".format(data_name),
        "softwareName": "RidgePredict-1.1u1",
        "requestedTime": "00:30:00",
        "archive": True,
        "inputs": input_dict,
        "parameters": {
            "outputPed": "RidgePredict_{}.ped".format(data_name)
        }
    }

    file_name = 'RidgePredict_{}.json'.format(data_name)
    with open(file_name, 'w') as f:
        json.dumps(json_text, f)
    return file_name


def bayesr_json(data_name, input_dict):
    """Generates an Agave-compatible JSON for BayesR.

    This JSON is stored in a file named "BayesR_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    :param input_dict: A dictionary generated by the "parse_inputs" function
    """
    json_text = {
        "jobName": "BayesR_{}".format(data_name),
        "softwareName": "bayesR-2.00u1",
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

    file_name = 'BayesR_{}.json'.format(data_name)
    with open(file_name, 'w') as f:
        json.dumps(json_text, f)
    return file_name


def plink_json(data_name):
    """Generates an Agave-compatible JSON for PLINK.

    This JSON is stored in a file named "PLINK_<dataset_name>.json".

    :param data_name: name of the dataset (without any extensions)
    """
    json_text = {
        "jobName": "Plink_{}".format(data_name),
        "softwareName": "PLINK-hpc-1.07u2",
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

    file_name = 'PLINK_{}.json'.format(data_name)
    with open(file_name, 'w') as f:
        json.dumps(json_text, f)
    return file_name


#TODO make current pipeline implementation work with this input format
def qxpak_json(data_name, input_dict):
    json_text = {
        "name": "QxPak_{}".format(data_name),
        "appId": "qxpak-stampede-5.05u2",
        "archive": True,
        "inputs": {
            "parameter_file": "agave://data.iplantcollaborative.org/home/shared/iplantcollaborative/example_data/qxpak/input/parameterFile.par",
            "data_file": "agave://data.iplantcollaborative.org/home/shared/iplantcollaborative/example_data/qxpak/input/dataFile.dat",
            "pedigree_file": "agave://data.iplantcollaborative.org/home/shared/iplantcollaborative/example_data/qxpak/input/pedigreeFile.ped",
            "marker_file": "agave://data.iplantcollaborative.org/home/shared/iplantcollaborative/example_data/qxpak/input/markerFile.mkr"
        },
        "parameter": {
            "output": "test"
        }

    }

#TODO find Gemma JSON?
def gemma_json(data_name, input_dict):
    pass


def make_gwas_json(selected_gwas, dataset_name, inputs):
    """Selected GWAS is a list of binary values with each element representing
        whether a given gwas should be run or not.

    :param selected_gwas: List of GWAS applications to run the data through
    :return:
    """
    switcher = {
        0: fastlmm_json(dataset_name, inputs),
        1: ridge_json(dataset_name, inputs),
        2: bayesr_json(dataset_name, inputs),
        3: plink_json(dataset_name),
        4: qxpak_json(dataset_name, inputs),
        5: gemma_json(dataset_name, inputs)
    }

    return switcher.get(selected_gwas)


def make_winnow_json(job_name, gwas_output_folder, ote_file):
    """" Builds JSON files for running Winnow via Agave.

    Built from JSON file here: https://raw.githubusercontent.com/CyVerse-Validate/Stampede-Files/master/json/winnow-job.json
    :param job_name: The name of the GWAS program frmo which the output came.
    :param gwas_output_folder: A folder containing GWAS outputs.
    :param ote_file: OTE known-truth file for truth testing.
    :return:
    """

    json_text = {
        "jobName" : 'Winnow_{}'.format(job_name),
        "softwareName" : "Winnow-0.9.2u1",
        "requestedTime" : "02:00:00",
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

    file_name = 'winnow-{}.json'.format(job_name)
    with open(file_name, 'w') as f:
        json_text = json.dumps(json_text, f)
    return file_name