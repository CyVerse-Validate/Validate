import json
import os
import agavepy.agave as a
import argparse
from uuid import getnode as get_mac

import datetime
import requests
from getpass import getpass
import time

import sys
from dateutil.tz import tzoffset

import JsonBuilder
from requests.packages.urllib3.exceptions import InsecureRequestWarning
from time import gmtime, strftime

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
        self.agave_token = None
        self.access_token = None
        self.run_date = strftime("%Y-%m-%d_%H:%M:%S", gmtime())

        self.username = ""
        self.password = ""
        self.systemid = ""
        self.iplant = 'https://agave.iplantc.org/files/v2/media/system/data.iplantcollaborative.org'
        self.home_dir = ""
        self.home_full = ""
        self.validate_dir = ""

        self.data_folder = ""       # User defined folder from the Datastore
        self.input_format = ""      # Type of inputs the user is providing
        self.inputs = {}
        self.known_truth = ""
        self.dataset_name = ""
        self.desired_gwas = ()      # Tuple of booleans
        self.running_jobs = {}      # Running jobs dictionary with the format { 'id' : 'archivePath' }
        self.finished_gwas = {}     # Finished jobs dictionary with the format { 'id' : list[RemoteFile] }
        self.output_folders = {}    # Output folder location within the 'Validate' directory
        self.clientName = ""
        self.checkArgs()

        print "Beginning Agave init."
        self.agave_initialization()
        # self.cyverse_test()
        self.validate()

    def agave_initialization(self):
        """Initializes Agave client via AgavePy for executing desired jobs on Stampede.

        :return:
        """

        if self.username is None:
            self.username = raw_input("Username:")

        # If no password is passed via arguments, prompts the user to enter it
        # in their terminal. This method securely handles the password.
        if self.password is None:
            self.password = getpass()

        # Unique client name is generated using the user's username and MAC
        # address. This is done to allow for multiple pipeline clients on
        # separate systems. Takes the format "<username>-pipelineClient-<mac>"
        self.clientName = "{}-pipelineClient-{}".format(self.username, get_mac()) \
            if (self.clientName == None) else self.clientName

        # TODO delete debug statement.
        print "{a}:\t{b}".format(a="clientName", b=self.clientName)

        # # Establishing connection with Agave using the user's allocation username and password
        # self.agave = a.Agave(api_server='https://agave.iplantc.org',
        #                      username=self.username, password=self.password,
        #                      client_name=self.clientName, verify=False)

        # if a.recover(self.clientName) == None:
        #     #except a.AgaveError:
        #     print "No pipeline client cached in the local agavepy db. Searching for existing client..."
        #     self.client = [cl for cl in self.agave.clients.list() if cl['name'] == self.clientName]
        #     if not self.client:
        #         print "Making new pipeline client: " + self.clientName
        #         self.client = self.agave.clients.create(
        #             body={'clientName': self.clientName})
        #     else:
        #         print "Existing pipline client found, but not stored locally. Recreating: " + self.clientName
        #         self.agave.clients.delete(clientName=self.clientName)
        #         self.client = self.agave.clients.create(body={'clientName': self.clientName})

        try:
            # TODO try instantiating without providing clientname - then search after.
            self.agave = a.Agave(api_server='https://agave.iplantc.org',
                                 username=self.username, password=self.password,
                                 client_name=self.clientName, verify=False)
            # a.recover(self.clientName)

        except KeyError:
            # TODO adapt client checks and make new client if needed.
            self.agave = a.Agave(api_server='https://agave.iplantc.org',
                                 username=self.username, password=self.password,
                                 verify=False)

            print "No pipeline client cached in the local agavepy db. Searching for existing client..."
            self.client = [cl for cl in self.agave.clients.list() if
                           cl['name'] == self.clientName]
            if not self.client:
                print "Making new pipeline client: " + self.clientName
                self.client = self.agave.clients.create(
                    body={'clientName': self.clientName})
            else:
                print "Existing pipline client found, but not stored locally. Recreating: " + self.clientName
                self.agave.clients.delete(clientName=self.clientName)
                self.client = self.agave.clients.create(
                    body={'clientName': self.clientName})

        # Original Validate Pipeline method.
        # Check for an existing Agave client
        # If none exists, a new one is created
        # for clnt in self.agave.clients.list():
        #     print clnt
        # print self.agave.token
        # self.client = [cl for cl in self.agave.clients.list() if cl['name'] == self.clientName]
        # if not self.client:
        #     print "Making new Agave client: " + self.clientName
        #     self.client = self.agave.clients.create(
        #         body={'clientName': self.clientName})
        # else:
        #     print "Pipeline client extant."
        #     self.agave.clients.delete(clientName='pipelineClient')
        #     self.client = self.agave.clients.create(body={'clientName': self.clientName})

        # TODO fix hardcoding for home_full, etc.
        # Grabbing token and links from the created Agave client
        self.agave_token = self.agave.token
        self.access_token = self.agave_token.token_info['access_token']
        # print self.agave_token.token_info
        self.home_dir = '/{}'.format(self.username)
        self.home_full = self.iplant + self.home_dir
        self.validate_dir = '/Validate/{}'.format(self.run_date)

    # TODO delete test method
    def cyverse_test(self):
        """Test method for checking connection to CyVerse via AgavePy.

        :return:
        """
        job_test = {}
        job_test['530158803913141785-242-ac113-0001-007'] = {
            u'status': u'PENDING', u'inputs': {u'tped': [
            u'agave://data.iplantcollaborative.org/mjs3607/PumaTest/X_test.tped'],
                                            u'tfam': [
                                                u'agave://data.iplantcollaborative.org/mjs3607/PumaTest/X_test.tfam']},
         u'localId': None, u'memoryPerNode': 4.0,
         u'archiveSystem': u'data.iplantcollaborative.org',
         u'processorsPerNode': 1,
         u'submitTime': datetime.datetime(2017, 7, 26, 14, 16, 52,
                                          tzinfo=tzoffset(None, -18000)),
         u'executionSystem': u'stampede.tacc.utexas.edu', u'startTime': None,
         u'appId': u'Puma-1.0u1', u'owner': u'mjs3607',
         u'id': u'530158803913141785-242ac113-0001-007', u'retries': 0,
         u'name': u'PumaTest',
         u'parameters': {u'penalty': u'LASSO', u'name': u'PumaTest',
                         u'regression': u'LINEAR'}, u'batchQueue': u'serial',
         u'nodeCount': 1, u'lastModified': u'2017-07-26T14:16:52.000-05:00',
         u'created': u'2017-07-26T14:16:52.000-05:00',
         u'archivePath': u'mjs3607/archive/jobs/job-530158803913141785-242ac113-0001-007',
         u'archive': True, u'outputPath': None, u'maxRunTime': u'00:05:00',
         u'_links': {u'notifications': {
             u'href': u'https://agave.iplantc.org/notifications/v2/?associatedUuid=530158803913141785-242ac113-0001-007'},
                     u'archiveSystem': {
                         u'href': u'https://agave.iplantc.org/systems/v2/data.iplantcollaborative.org'},
                     u'notification': [], u'self': {
                 u'href': u'https://agave.iplantc.org/jobs/v2/530158803913141785-242ac113-0001-007'},
                     u'metadata': {
                         u'href': u'https://agave.iplantc.org/meta/v2/data/?q=%7B%22associationIds%22%3A%22530158803913141785-242ac113-0001-007%22%7D'},
                     u'archiveData': {
                         u'href': u'https://agave.iplantc.org/jobs/v2/530158803913141785-242ac113-0001-007/outputs/listings'},
                     u'executionSystem': {
                         u'href': u'https://agave.iplantc.org/systems/v2/stampede.tacc.utexas.edu'},
                     u'owner': {
                         u'href': u'https://agave.iplantc.org/profiles/v2/mjs3607'},
                     u'history': {
                         u'href': u'https://agave.iplantc.org/jobs/v2/530158803913141785-242ac113-0001-007/history'},
                     u'app': {
                         u'href': u'https://agave.iplantc.org/apps/v2/Puma-1.0u1'},
                     u'permissions': {
                         u'href': u'https://agave.iplantc.org/jobs/v2/530158803913141785-242ac113-0001-007/pems'}},
         u'endTime': None}

        print "Access Token:    {}".format(self.agave_token.token_info['access_token'])
        print "Refresh Token:   {}".format(self.agave_token.token_info['refresh_token'])
        print "API Key:     {}".format(self.agave_token.api_key)
        print "API Secret:  {}\n".format(self.agave_token.api_secret)
        print "Home dir:    {}".format(self.home_dir)
        print "Home full:   {}".format(self.home_full)
        print "Data folder: {}".format(self.data_folder)
        print "Full path:   {}".format(self.home_full + self.data_folder)
        print "Validate:    {}".format(self.validate_dir)
        print "ValidateF:   {}".format(self.home_full + self.validate_dir)
        print "\n\n\n"

        # systems = self.agave.systems.list()
        # for system in systems:
        #     print system
        #     print "\n\n\n"

        print "Systems:"
        system_list = self.agave.systems.list()
        for system in system_list:
            print "{}\n".format(system)

        print "\n\nFiles:"
        file_list = [f for f in self.agave.files.list(
            systemId=self.systemid, filePath=self.data_folder)]
        for file in file_list:
            print "agave://{}{}".format(file['system'], file['path'])

        # print "\n\nApps:"
        # apps_list = self.agave.apps.list()
        # for app in apps_list:
        #     print "{}:\n".format(app['id'])

        # remote_file_list = self.agave.files.list(
        #     systemId=self.systemid,
        #     filePath=job_test['530158803913141785-242-ac113-0001-007']['archivePath'])
        #
        # for f in remote_file_list:
        #     print f

    def validate(self):
        """Main run method for the Validate Pipeline.

        This method handles all stages of the Pipeline as defined below:
            1. Parsing of inputs and retrieving file links from the desired
                storage system.
            2. Building GWAS JSON files for individual jobs in preparation to
                submit to an execution system via Agave. Note, the execution
                system for most apps is contained within the app definitions
                stored on Agave.
            3. Submits the previously generated GWAS JSON files and handles the
                polling of submitted jobs via "poll_jobs".
            4. When all jobs have finished, outputs are moved from the job
                archive system into the Validate folder within the user's home
                directory on their chosen storage system. Within this folder,
                individual runs of Validate are encapsulated within their own
                folders to allow for the storage/archive of multiple runs.
            5. Winnow JSONs are created and submitted for each previously run
                GWAS job, as long as known truth files exist for the given
                dataset. Currently, Winnow only accepts OTE files.

        Returns:

        """
        # TODO add option for the user to upload their own data from their local machines
        # TODO add simulation option
        self.parse_inputs()

        # TODO add with_refresh to all methods.
        if self.gwas:
            self.build_jsons()

            # TODO build gwas files from the given folder.
            self.gwas_submission()
            self.finished_gwas = self.poll_jobs(self.running_jobs)
            self.make_output_folders(self.parse_archives(self.finished_gwas))

        if self.winnow:
            self.winnow_submission(self.create_winnow_jsons(self.output_folders))
            self.finished_winnow = self.poll_jobs(self.running_jobs)
            self.make_winnow_folders(self.parse_archives(self.finished_winnow,
                                                         winnow=True))
        # TODO Equalise outputs

    def checkArgs(self):
        """Method for parsing the provided command-line arguments.

        Attributes:
            -u, --username:
                String. The user's username for their preferred login system.
            --password:
                String. The user's password for their preferred login system.
            -i, --InFormat:
                String. Format of the files the user is submitting for analysis.
                Options are "p" for PED/MAP, "b" for BIM/BED/FAM, and "t" for
                TPED/TFAM.
            -f, --Folder:
                String. Data folder for analysis. This should contain all input
                data (of the format defined by -i, although others may exist
                there as well but will not be analysed) as well as the
                known-truth file.
            -pno, --pheno:
                String. Name (with extension) of the phenotype file. Requirement
                of this argument/file is specific to different GWAS definitions.
            -lmm, --fastlmm:
                Binary value. If True, FaST-LMM will be run on the dataset.
            -rdg, --ridge:
                Binary value. If True, RidgePredict will be run on the dataset.
            -bay, --bayes:
                Binary value. If True, BayesR will be run on the dataset.
            -plk, --plink:
                Binary value. If True, PLINK will be run on the dataset.
            -qxp, --qxpak:
                Binary value. If True, QxPak will be run on the dataset.
            -gma, --gemma:
                Binary value. If True, Gemma will be run on the dataset.
            -pma, --puma:
                Binary value. If True, Puma will be run on the dataset.
            --no_gwas:
                Binary value. If True, no GWAS will be run. Instead, the
                provided input folder will be interpreted as a folder for
                Winnow and all GWAS stages will be omitted.
            --no_winnow:
                Binary value. If true, Winnow will not be run. GWAS analyses
                will be submitted and handled as defined above, but no Winnow
                jobs will be created or submitted after data is moved to the
                Validate folder.
            -sys, --system:
                String. System ID for the user's preferred filesystem. This is
                the system where all data will be retrieved from and stored to.
            -ote, --truth:
                String. Optional. Used to give an absolute path to the known
                truth file for the given dataset if it does not reside in the
                provided data folder (-f or --Folder).
            -c, --client_name:
                String. Allows the user to define an existing client for usage
                with the Validate Pipeline. If this client is not found, this
                string will instead be used to create a new client with the
                provided name.

        Returns:

        """
        # TODO add simulation || prediction || gwas
        # if simulation, make simulated data and folder for data before validating
        parser = argparse.ArgumentParser()
        parser.add_argument("-u", "--username", type=str,
                            default=os.getenv("AGAVE_USERNAME"),
                            help="Username used for Agave services.")
        parser.add_argument("--password", type=str,
                            default=os.getenv("AGAVE_PASSWORD"),
                            help="Password used for Agave services.")
        parser.add_argument("-i", "--InFormat", type=str,
                            help="Input format:\n"
                                 "\tp for PED/MAP\n"
                                 "\tb for BIM/BED/FAM\n"
                                 "\tt for TPED/TFAM")
        parser.add_argument("-f", "--Folder", type=str,
                            help="Folder to be pipelined. This folder should"
                                 "contain all input data as well as the known-truth"
                                 "file for the given data set and be a path"
                                 "relative to the \"/iplant/home\" directory.")
        parser.add_argument("-pno", "--pheno", type=str,
                            help="Name (including extension) for the covariate"
                                 "file for the given dataset.")
        parser.add_argument("-lmm", "--fastlmm", action="store_true",
                            help="Apply this flag if FaST-LMM is to be run.")
        parser.add_argument("-rdg", "--ridge", action="store_true",
                            help="Apply this flag if Ridge is to be run.")
        parser.add_argument("-bay", "--bayes", action="store_true",
                            help="Apply this flag if BayesR is to be run.")
        parser.add_argument("-plk", "--plink", action="store_true",
                            help="Apply this flag if PLINK is to be run.")
        parser.add_argument("-qxp", "--qxpak", action="store_true",
                            help="Apply this flag if QxPak is to be run.")
        parser.add_argument("-gma", "--gemma", action="store_true",
                            help="Apply this flag if Gemma is to be run.")
        parser.add_argument("-pma", "--puma", action="store_true",
                            help="Apply this flag if Puma is to be run.")
        parser.add_argument("--no_gwas", action="store_false",
                            help="Apply this flag if gwas is not to be run. In"
                                 " this scenario, the input folder will be"
                                 " taken as the Winnow folder.")
        parser.add_argument("--no_winnow", action="store_false",
                            help="Apply this flag if Winnow is not to be run. "
                                 "In this scenario, GWAS output data will be "
                                 "moved to the Validate directory and the "
                                 "process will complete.")
        parser.add_argument("-sys", "--system", type=str,
                            default="data.iplantcollaborative.org",
                            help="System ID for the desired filesystem. Defaults"
                                 "to data.iplantcollaborative.org if nonne is"
                                 "provided.")
        parser.add_argument("-ote", "--truth", type=str,
                            help="If the known-truth file to be used for Winnow"
                                 "validation is not located within the given"
                                 "input folder, it may be specificed separately"
                                 "here.")
        parser.add_argument("-c", "--client_name", type=str, default=None,
                            help="Name of existing Validate Pipeline client."
                                 "If none is supplied, checks the Agave cache."
                                 "If name is supplied but no client exists, a "
                                 "new client with the provided name is created.")

        # TODO get parameters for each GWAS somehow - potentially JSON?
        # TODO Add option for user to pass in their own JSONs instead of a folder

        args = parser.parse_args()
        self.gwas = args.no_gwas
        self.winnow = args.no_winnow
        self.username = args.username
        self.password = args.password
        self.input_format = args.InFormat
        self.data_folder = args.Folder if args.Folder.startswith("/") else "/" + args.Folder
        self.desired_gwas = tuple([args.fastlmm, args.ridge, args.bayes, args.plink,
                                   args.qxpak, args.gemma, args.puma])
        self.dataset_name = self.data_folder.split("/")[-1]
        self.systemid = args.system
        self.known_truth = "agave://{}".format(os.path.join(self.systemid, args.truth))
        self.clientName = args.client_name
        # TODO Add option for expandable apps?

    def parse_inputs(self):
        """Grabs individual links to all input/OTE files defined by the user.

        Begins by listing all files within the provided data folder (-f or
        --Folder) and, depending on the provided input format (-i, --InFormat),
        creates Agave-suitable links to all files matching this input format. If
        any ".ote" files are found, these links are also generated.

        Returns:

        """
        # TODO pull phenotype file too
        file_list = [f for f in self.agave.files.list(
            systemId=self.systemid, filePath=self.data_folder)]

        # TODO error check if no input file is found
        # If the input data was declared as PED/MAP
        # TODO Restrict to ONE set ONLY.
        # TODO Allow for explicit filename passing.
        if self.input_format == 'p':
            for file in file_list:
                if ".ote" in file.name:
                    print "Found ote file: {}".format(file.name)
                    self.known_truth = "agave://{}{}".format(file['system'], file['path'])
                elif ".ped" in file.name:
                    self.inputs['inputPED'] = "agave://{}{}".format(file['system'], file['path'])
                elif ".map" in file.name:
                    self.inputs['inputMAP'] = "agave://{}{}".format(file['system'], file['path'])

        # If the input data was declared as BIM/BED/FAM
        elif self.input_format == 'b':
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = "agave://{}{}".format(file['system'], file['path'])
                elif '.bed' in file.name:
                    self.inputs['inputBED'] = "agave://{}{}".format(file['system'], file['path'])
                elif '.bim' in file.name:
                    self.inputs['inputBIM'] = "agave://{}{}".format(file['system'], file['path'])
                elif '.fam' in file.name:
                    self.inputs['inputFAM'] = "agave://{}{}".format(file['system'], file['path'])

        # If the input data was declared as TPED/TMAP
        elif self.input_format == 't':
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = "agave://{}{}".format(file['system'], file['path'])
                elif ".tped" in file.name:
                    self.inputs['tped'] = "agave://{}{}".format(file['system'], file['path'])
                elif ".tfam" in file.name:
                    self.inputs['tfam'] = "agave://{}{}".format(file['system'], file['path'])

        # TODO catch this earlier (within arg checking)
        else:
            raise ValueError("Incorrect input format declared.")

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
                print("desired_gwas:\t{}\ngwas index:\t{}".format(self.desired_gwas, self.desired_gwas.index(gwas)))
                json, _ = JsonBuilder.make_gwas_json(
                    self.desired_gwas.index(gwas), self.dataset_name, self.inputs)
                self.gwas_jsons.append(json)
                # else:
                #     Make Winnow JSONs here.
                # self.winnow_jsons = []
                # for gwas_output in self.finished_gwas:
                #     pass

        # TODO delete print
        for JSON in self.gwas_jsons:
            print "{}".format(JSON)

    def gwas_submission(self):
        """Submits all provided JSON files via AgavePy.
        All currently running jobs are stored in the 'running_jobs' dictionary,
        which has the format of { 'job[id]' : 'job' } with 'Job' being the
        response from the server given upon job submission.

        :return:
        """
        # TODO delete print statement
        print "Beginning submission...\n"

        for JSON in self.gwas_jsons:
            print "Submitting: {}".format(json.dumps(JSON, indent=4, separators=(',', ': ')))
            job = self.agave.jobs.submit(body=JSON)
            self.running_jobs[job['id']] = job
            print "Submitted: {}".format(job['id'])

        # TODO delete print statements
        for jb in self.running_jobs.keys():
            print "{}\t\t{}".format(jb, self.running_jobs[jb])

        print "\nSubmission finished\n\n"

    def poll_jobs(self, job_dict):
        """Polls all running jobs for status until completion.

        Running jobs are stored in the running_jobs dictionary. Once finished,
        jobs are removed and added to the finished_gwas dictionary for easy
        tracking and handling for running through Winnow.
        """
        # Instead of polling jobs, check if a directory has been created in the archive
        sleep_time = 60
        finished_jobs = {}

        # Iterating through all JobIDs and polling until there are no more jobs
        while job_dict.keys():
            print "Begin polling... {} remaining.".format(len(job_dict.keys()))

            for job_id in job_dict.keys():
                job_status = self.agave.jobs.getStatus(jobId=job_id)

                # If the curernt given job is finished, download output and
                # remove it from the running queue. The finished job and its
                # output list is added to the finished_gwas dictionary.
                # TODO Check for failed jobs / etc.
                if (job_status['status'] == "FINISHED"):
                    print "{} has finished.".format(job_id)
                    finished_jobs[job_id] = job_dict[job_id]
                    del job_dict[job_id]

                elif (job_status['status'] == "FAILED"):
                    print "{} has FAILED. Removing...".format(job_id)
                    del job_dict[job_id]

                else:
                    print "{} is unfinished. Currently {}.".format(job_id, job_status['status'])

            # Sleep before repolling Agave. Max sleep time is 1 hour.
            if job_dict.keys():
                time.sleep(sleep_time)

            # sleep_time *= 1.25 if sleep_time <= 3600 else sleep_time
            print("Jobs finished:  {}".format(finished_jobs.keys()))
            print "Jobs remaining: {}\n".format(job_dict.keys())

        return finished_jobs

    def parse_archives(self, jobid_dict, winnow=False):
        """Parses job output archives on the Discovery Environment via Agave.

        :param jobid_dict: Dictionary with the Agave job-id as keys.
        :return:
        """
        # if not winnow:
        #     out_extensions = ['.out.txt', '.freq', '.gv', '.hyp', '.model',
        #                       '.param', '.R']
        # else:
        #     out_extensions = ['.txt']
        job_output_dict = {}

        # Iterates thorugh all Job-IDs and collects the output files for each
        for jobid in jobid_dict.keys():
            job_outputs = []
            remote_file_list = self.agave.files.list(
                systemId=self.systemid,
                filePath=jobid_dict[jobid]['archivePath'])

            for file in remote_file_list:
                # if any(ext in file['name'] for ext in out_extensions):
                job_outputs.append(file['path'])
                print("\tAppending output file {}".format(file['name']))

            job_output_dict[jobid] = job_outputs

        return job_output_dict

    # def parse_archives(self, jobid_dict, winnow=False):
    #     # excluded_ext = [".log", ".err", ".out"]
    #     job_output_dict = {}
    #
    #     for jobid in jobid_dict.keys():
    #         job_outputs = []
    #         remote_file_list = self.agave.files.list(
    #             systemId=self.systemid,
    #             filePath=jobid_dict[jobid]['archivePath'])
    #
    #         for file in remote_file_list:
    #             job_outputs.append(file['path'])
    #             print("\tAppending output file {}".format(file['name']))
    #
    #         job_output_dict[jobid] = job_outputs
    #
    #     return job_output_dict


    # TODO Translate
    def make_output_folders(self, job_outputs):
        """Makes folders of all outputs for each finished GWAS job.

        :param job_outputs: Dictionary of { 'jobid' : list(output_paths) }
        :return:
        """
        # TODO add timestamp to each validate 'run' OR user-provided name
        self.output_folders = {}

        newf = self.agave.files.manage(systemId=self.systemid,
                           filePath='{}'.format(self.home_dir),
                           body={'action':'mkdir',
                                       'path':'{}'.format(self.validate_dir)})
        print newf

        # Loops through all finished Job IDs and creates a subdirectory within
        # the 'Validate GWAS Outputs' folder defined above simply named after
        # each Job ID. All outputs are stored here for easy access.
        for jobid in job_outputs.keys():
            newf = self.agave.files.manage(systemId=self.systemid,
                           filePath='{}{}'.format(self.home_dir, self.validate_dir),
                           body={'action':'mkdir',
                                       'path':"GWAS/{}".format(jobid)})
            print newf

            # Copying all job outputs from the system archive directory to the
            # newly created subdirectory, leaving the original archive as is.
            for file in job_outputs[jobid]:
                print "Copying {} to {}/GWAS/{}".format(file, self.validate_dir, jobid)
                copyf = self.agave.files.manage(systemId=self.systemid,
                           filePath='{}'.format(file),
                           body={'action':'copy',
                                       'path':"{}{}/GWAS/{}".format(
                                           self.home_dir, self.validate_dir, jobid)})
                print copyf

            # TODO get Validate GWAS Outputs full path
            self.output_folders[jobid] = "agave://{}{}{}/GWAS/{}".format(
                self.systemid, self.home_dir, self.validate_dir, jobid)

            # Checking files in the Validate directory after copying
            print "Files in {}/GWAS/ after copy:".format(self.home_dir + self.validate_dir)
            validate_files = [f for f in self.agave.files.list(
                systemId=self.systemid,
                filePath="{}/GWAS/{}".format(self.home_dir + self.validate_dir, jobid))]

            for file in validate_files:
                print file

    def create_winnow_jsons(self, output_folders):
        # Files are located in /iplant/home/<user-dir>/archive/jobs/job-<jobid>
        winnow_jsons = []

        for jobid in output_folders.keys():
            winnow_js = JsonBuilder.make_winnow_json(
                jobid, output_folders[jobid], self.known_truth)

            print json.dumps(winnow_js, indent=4, separators=(',', ': '))

            winnow_jsons.append(winnow_js)

        return winnow_jsons

    def winnow_submission(self, winnow_jsons):
        #FaST-LMM uses the base filename for the input files for its output
        # e.g./ <input-file-name>.out.txt
        #
        # BayesR gives 6 outputs: .txt.frq, .txt.gv, .txt.hyp, .txt.log,
        #  .txt.model, and .txt.param
        #
        # Ridge gibt eine einzelner Ausgabedatei -- velleicht PED?
        #
        # PLINK gives plain-text, space-delimited output files.
        #
        # QxPak single output file - qxpak.out.
        #
        # Gemma ???
        # for JSON in winnow_jsons:
        #     print "Submitting: {}".format(JSON)
        #     winnow_submission = self.agave.jobs.submit(body=JSON)
        #     self.running_jobs[winnow_submission['id']] = winnow_submission
        #     print "Submitted: {}".format(JSON)

        for JSON in winnow_jsons:
            print "Submitting: {}".format(json.dumps(JSON, indent=4, separators=(',', ': ')))
            job = self.agave.jobs.submit(body=JSON)

            self.running_jobs[job['id']] = job
            print "Submitted: {}".format(job['id'])

        # TODO delete print statements
        for jb in self.running_jobs.keys():
            print "{}\t\t{}".format(jb, self.running_jobs[jb])

        print "\nSubmission finished\n\n"

    def make_winnow_folders(self, job_outputs):
        """Makes folders of all outputs for each finished GWAS job.

        :param job_outputs: Dictionary of { 'jobid' : list(output_paths) }
        :return:
        """
        # TODO add timestamp to each validate 'run' OR user-provided name

        # Loops through all finished Job IDs and creates a subdirectory within
        # the 'Validate Winnow Outputs' folder defined above simply named after
        # each Job ID. All outputs are stored here for easy access.
        for jobid in job_outputs.keys():
            newf = self.agave.files.manage(systemId=self.systemid,
                           filePath='{}{}'.format(self.home_dir, self.validate_dir),
                           body={'action':'mkdir',
                                       'path':"Winnow/{}".format(jobid)})
            print newf

            # Copying all job outputs from the system archive directory to the
            # newly created subdirectory, leaving the original archive as is.
            for file in job_outputs[jobid]:
                print "Copying {} to {}/Winnow/{}".format(file, self.validate_dir, jobid)
                copyf = self.agave.files.manage(systemId=self.systemid,
                           filePath='{}'.format(file),
                           body={'action':'copy',
                                       'path':"{}{}/Winnow/{}".format(
                                           self.home_dir, self.validate_dir, jobid)})
                print copyf

            # TODO get Validate GWAS Outputs full path
            self.output_folders[jobid] = "agave://{}{}/GWAS/{}".format(
                self.home_dir, self.validate_dir, jobid)

            # Checking files in the Validate directory after copying
            print "Files in {}/Winnow/ after copy:".format(self.home_dir + self.validate_dir)
            validate_files = [f for f in self.agave.files.list(
                systemId=self.systemid,
                filePath="{}/Winnow/{}".format(self.home_dir + self.validate_dir, jobid))]

            for file in validate_files:
                print file


if __name__ == '__main__':
    Pipeline()
