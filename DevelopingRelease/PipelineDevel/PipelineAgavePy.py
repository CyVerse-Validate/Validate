import agavepy.agave as a
import argparse
from getpass import getpass
import json

__author__ = "Michael J. Suggs"
__credits__ = ["Michael Suggs"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Michael Suggs"
__email__ = "mjs3607@uncw.edu"
__status__ = "Development"


class Pipeline:
    """Handles automatic job submission and data-handling for the Validate Pipeline.

    Attributes:
        data_folder : str
                    User-provided Agave-style folder location in string format
                    detailing the data to be processed by the Pipeline.
        inputs : dict
        known_truth : str
                    OTE file found within the provided data folder.
        dataset_name : str
                    Name of the dataset - this is extracted from provided data folder name.
        desired_gwas : tuple(bool)
                    Tuple of booleans collected via command-line arguments
                    detailing the GWAS to be run.
        running_jobs : dict{ 'id' : 'archivePath' }
                    Dictionary of currently-running jobs as { 'id' : 'archivePath' }
        finished_gwas : dict{ 'id' : list[RemoteFile] }
                    Dictionary of finished GWAS jobs as { 'id' : list[RemoteFile] }
        output_folders: dict{ 'id' : list[RemoteFile] }
                    Each ID is associated with a subdirectory within the
                    Validate directory on the Data Store.
    """

    def __init__(self):
        # self.user_home = Agave.filesApi. TODO get user home directory
        self.agave = None
        self.username = ""          # User's username
        self.data_folder = ""       # User defined folder from the Datastore
        self.input_format = ""      # Type of inputs the user is providing
        self.inputs = {}
        self.known_truth = ""
        self.dataset_name = ""
        self.desired_gwas = ()      # Tuple of booleans
        self.running_jobs = {}      # Running jobs dictionary with the format { 'id' : 'archivePath' }
        self.finished_gwas = {}     # Finished jobs dictionary with the format { 'id' : list[RemoteFile] }
        self.output_folders = {}    # Output folder location within the 'Validate' directory

        self.agave_initialization()

    def checkArgs(self):
        """Checks all command-line arguments provided in the command-line call.

        :return:
        """
        # TODO add simulation || prediction || gwas
        # if simulation, make simulated data and folder for data before validating
        parser = argparse.ArgumentParser()
        parser.add_argument("-u", "--username", type=str,
                            help="Username used for Agave Services.")
        parser.add_argument("-i", "--InFormat", type=str,
                            help="Input format:\n"
                                 "\tp for PED/MAP\n"
                                 "\tb for BIM/BED/FAM\n"
                                 "\tt for TPED/TFAM")
        parser.add_argument("-f", "--Folder", type=str,
                            help="Folder to be pipelined. This folder should"
                                 "contain all input data as well as the known-truth"
                                 "file for the given data set.")
        parser.add_argument("-l", "--fastlmm", type=bool,
                            help="\"True\" if FaST-LMM is to be run.")
        parser.add_argument("-r", "--ridge", type=bool,
                            help="\"True\" if Ridge is to be run.")
        parser.add_argument("-b", "--bayes", type=bool,
                            help="\"True\" if BayesR is to be run.")
        parser.add_argument("-p", "--plink", type=bool,
                            help="\"True\" if PLINK is to be run.")
        parser.add_argument("-q", "--qxpak", type=bool,
                            help="\"True\" if QxPak is to be run.")
        parser.add_argument("-g", "--gemma", type=bool,
                            help="\"True\" if Gemma is to be run.")
        # TODO get parameters for each GWAS somehow - potentially JSON?

        args = parser.parse_args()
        self.username = args.username
        self.input_format = args.InFormat
        self.data_folder = args.Folder
        self.desired_gwas = tuple([args.fastlmm, args.ridge, args.bayes, args.plink,
                                   args.qxpak, args.gemma])

    def agave_initialization(self):
        password = getpass()
        self.agave = a.Agave(api_server='https://agave.iplantc.org',
                             username=self.username, password=password,
                             verify=False)

        # Check for an existing Agave client
        # If none exists, a new one is created
        client = [cl for cl in self.agave.clients.list() if cl['name'] == 'pipelineClient']
        if not client:
            client = self.agave.clients.create(body={'clientName': 'pipelineClient'})

        # Creating CyVerse DE storage system connections
        cyverse_de = {
            "id": "de.cyverseorg.",
            "name": "Agavepy test storage system",
            "status": "UP",
            "type": "STORAGE",
            "description": "Storage system for agavepy test",
            "site": "de.cyverse.org",
            "storage": {
                "host": STORAGE_IP, # 128.196.254.62 (?)
                "port": 22,
                "protocol": "SFTP",
                "rootDir": "/",
                "homeDir": "/home/{}".format(self.username),
                "auth": {
                    "username": STORAGE_USERNAME,
                    "type": "PASSWORD",
                    "password": STORAGE_PASSWORD
                }
            }
        }
        self.agave.systems.add(body=cyverse_de);

    def parse_inputs(self):
        """Grabs the known-truth and all input files from the given directory.

        :return:
        """
        # TODO pull phenotype file too
        file_list = [f for f in self.agave.files.list(
            systemId='de.cyverse.org', filePath='/home/{}'.format(self.username))]

        # TODO error check if no input file is found
        # If the input data was declared as PED/MAP
        if self.input_format == 'p':
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = file
                elif ".ped" in file.name:
                    self.inputs['inputPED'] = file
                elif ".map" in file.name:
                    self.inputs['inputMAP'] = file

        # If the input data was declared as BIM/BED/FAM
        elif self.input_format == 'b':
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = file
                elif '.bed' in file.name:
                    self.inputs['inputBED'] = file
                elif '.bim' in file.name:
                    self.inputs['inputBIM'] = file
                elif '.fam' in file.name:
                    self.inputs['inputFAM'] = file

        # If the input data was declared as TPED/TMAP
        else:
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = file
                elif ".tped" in file.name:
                    self.inputs['inputTPED'] = file
                elif ".tfam" in file.name:
                    self.inputs['inputTFAM'] = file
