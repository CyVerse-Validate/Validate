import agavepy.agave as a
import argparse
import requests
from getpass import getpass
import time
import json
import JsonBuilder
from requests.packages.urllib3.exceptions import InsecureRequestWarning

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
        # Cleaning up console output
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

        # self.user_home = Agave.filesApi. TODO get user home directory
        self.agave = None
        self.access_token = None

        self.username = 'mjs3607'
        self.password = 'Pixi3-W0rk-F00tpath'
        self.home_dir = '/{}'.format(self.username)

        self.data_folder = ""       # User defined folder from the Datastore
        self.input_format = ""      # Type of inputs the user is providing
        self.inputs = {}
        self.known_truth = ""
        self.dataset_name = ""
        self.desired_gwas = ()      # Tuple of booleans
        self.running_jobs = {}      # Running jobs dictionary with the format { 'id' : 'archivePath' }
        self.finished_gwas = {}     # Finished jobs dictionary with the format { 'id' : list[RemoteFile] }
        self.output_folders = {}    # Output folder location within the 'Validate' directory

        print "Beginning Agave init."
        self.agave_initialization()
        # self.cyverse_test()
        self.validate()

    def agave_initialization(self):
        # password = getpass()
        self.agave = a.Agave(api_server='https://agave.iplantc.org',
                             username=self.username, password=self.password,
                             verify=False)

        # Check for an existing Agave client
        # If none exists, a new one is created
        self.client = [cl for cl in self.agave.clients.list() if cl['name'] == 'pipelineClient']
        if not self.client:
            print "Making new Agave client."
            self.client = self.agave.clients.create(
                body={'clientName': 'pipelineClient'})
        else:
            print "Pipeline client extant."
            self.agave.clients.delete(clientName='pipelineClient')
            self.client = self.agave.clients.create(body={'clientName': 'pipelineClient'})

        self.access_token = self.agave.token

    def cyverse_test(self):
        file_list = [f for f in self.agave.files.list(
            systemId='data.iplantcollaborative.org', filePath=self.home_dir)]

        print "\n\nFiles:"
        for file in file_list:
            print file

    def validate(self):
        """Main run method for the Validate Workflow.

        :return:
        """
        # TODO add option for the user to upload their own data from their local machines
        # TODO add simulation option
        self.checkArgs()
        self.parse_inputs()
        self.build_jsons()

        # TODO build gwas files from the given folder.
        # self.gwas_submission()
        # self.finished_gwas = self.poll_jobs(self.running_jobs)
        # self.make_output_folders(self.parse_archives(self.finished_gwas))
        # TODO Equalise outputs
        # TODO Make Winnow JSONs
        # TODO Submit Winnow
        # TODO Retrieve Winnow outputs

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

    def build_jsons(self):
        """Builds JSONs from the parsed input information.

        If there are no finished GWAS jobs, this method will default to building
        JSONs for Winnow from the GWAS output files.

        :return:
        """
        # if not self.finished_gwas:
        self.gwas_jsons = []
        for gwas in self.desired_gwas:
            if gwas:
                self.gwas_jsons += JsonBuilder.make_gwas_json(
                    self.desired_gwas.index(gwas), self.dataset_name, self.inputs)
                # else:
                #     Make Winnow JSONs here.
                # self.winnow_jsons = []
                # for gwas_output in self.finished_gwas:
                #     pass

    def gwas_submission(self):
        """Submits all provided JSON files via the Agave REST APIs.
        All currently running jobs are stored in the 'running_jobs' dictionary,
        which has the format of { 'job[id]' : 'Job' } with 'Job' being the
        Job.py response from the server given upon job submission.

        :return:
        """
        for json in self.gwas_jsons:
            job = self.agave.jobs.submit(json)
            self.running_jobs[job['id']] = job

    def poll_jobs(self, job_dict):
        """Polls all running jobs for status until completion.

        Running jobs are stored in the running_jobs dictionary. Once finished,
        jobs are removed and added to the finished_gwas dictionary for easy
        tracking and handling for running through Winnow.
        """
        #TODO look at notifications instead?
        #TODO keep data on the datastore - no downloading
        # Instead of polling jobs, check if a directory has been created in the archive
        sleep_time = 60
        finished_jobs = {}

        # Iterating through all JobIDs and polling until there are no more jobs
        while not job_dict.keys():
            for job_id in job_dict.keys:
                job_status = self.agave.jobs.getStatus(job_id)

                # If the curernt given job is finished, download output and
                # remove it from the running queue. The finished job and its
                # output list is added to the finished_gwas dictionary.
                # TODO Check for failed jobs / etc.
                if (job_status['status'] == "FINISHED"):
                    # TODO use archive instead of downloading
                    # agave://data.iplantcollaborative.org/<user-home>/archive/jobs/job-<jobID>
                    finished_jobs[job_id] = job_dict[job_id]
                    del job_dict[job_id]

            # Sleep before repolling Agave. Max sleep time is 1 hour.
            time.sleep(sleep_time)
            sleep_time *= 2 if sleep_time <= 1800 else sleep_time

        return finished_jobs

    def parse_archives(self, jobid_dict):
        """Parses job output archives on the Discovery Environment via Agave.

        :param jobid_dict: Dictionary with the Agave job-id as keys.
        :return:
        """
        out_extensions = ['.out.txt', '.freq', '.gv', '.hyp', '.log', '.model',
                          '.param']
        job_output_dict = {}

        # Iterates thorugh all Job-IDs and collects the output files for each
        for jobid in jobid_dict.keys:
            job_outputs = []
            remote_file_list = self.agave.files.list(
                filePath=self.home_dir + "/archive/jobs/job-{}".format(jobid))

            for file in remote_file_list:
                if any(out_extensions) in file['name']:
                    job_outputs.append(file['path'])

            job_output_dict[jobid] = job_outputs

        return job_output_dict

    def make_output_folders(self, job_outputs):
        """Makes folders of all outputs for each finished GWAS job.

        :param job_outputs: Dictionary of { 'jobid' : list(output_paths) }
        :return:
        """
        # TODO add timestamp to each validate 'run' OR user-provided name
        self.output_folders = {}
        self.agave.files
        # {"action":"mkdir","path":"$NEWDIR"}
        FilesApi.manageOnDefaultSystem(sourcefilePath='.', action='mkdir',
                                             filePath='Validate')
        FilesApi.manageOnDefaultSystem(sourcefilePath='./Validate',
                                             action='mkdir', filePath='GWAS Outputs')

        # Loops through all finished Job IDs and creates a subdirectory within
        # the 'Validate GWAS Outputs' folder defined above simply named after
        # each Job ID. All outputs are stored here for easy access.
        for jobid in job_outputs.keys:
            FilesApi.manageOnDefaultSystem(
                sourcefilePath='./Validate/GWAS Outputs', action='mkdir',
                filePath=jobid)

            # Copying all job outputs from the system archive directory to the
            # newly created subdirectory, leaving the original archive as is.
            for file in job_outputs[jobid]:
                FilesApi.manageOnDefaultSystem(
                    sourcefilePath='./Validate/GWAS Outputs/{}'.format(jobid),
                    action='cp', filePath=file)

            # TODO get Validate GWAS Outputs full path
            self.output_folders[jobid] = FilesApi.listOnDefaultSystem(
                filePath='./Validate/GWAS Outputs/').swaggerTypes['result']

if __name__ == '__main__':
    Pipeline()