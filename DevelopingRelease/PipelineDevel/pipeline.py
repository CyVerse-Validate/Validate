import argparse, os, time
import JsonBuilder
import AgavePythonSDK as Agave

__author__ = "Michael J. Suggs::mjs3607@uncw.edu"
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
    """
    # TODO Change pipeline around to take user folder location

    def __init__(self):
        self.data_folder = ""       # User defined folder from the Datastore
        self.inputs = {}
        self.known_truth = ""
        self.dataset_name = ""
        self.desired_gwas = ()      # Tuple of booleans
        self.running_jobs = {}      # Running jobs dictionary with the format { 'id' : 'archivePath' }
        self.finished_gwas = {}     # Finished jobs dictionary with the format { 'id' : list[RemoteFile] }

        self.checkArgs()
        self.parse_inputs()
        self.build_jsons()

        # TODO build gwas files from the given folder.
        self.gwas_submission()
        self.poll_jobs()
        # TODO Equalise outputs
        # TODO Make Winnow JSONs
        # TODO Submit Winnow
        # TODO Retrieve Winnow outputs

    def checkArgs(self):
        """Checks all command-line arguments provided in the command-line call.

        :return:
        """
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", "--InFormat", type=chr(1),
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

        args = parser.parse_args()
        self.input_format = args.InFormat
        self.data_folder = args.Folder
        self.desired_gwas = tuple([args.lmm, args.ridge, args.bayes, args.plink,
                                   args.qxpak, args.gemma])

    def parse_inputs(self):
        """Grabs the known-truth and all input files from the given directory.

        :return:
        """
        file_list = Agave.FilesApi.listOnDefaultSystem(self.data_folder).swaggerTypes['result']

        # If the input data was declared as PED/MAP
        if self.input_format == 'p':
            for file in file_list:
                if ".ote" in file.swaggerTypes['name']:
                    self.known_truth = file.swaggerTypes['path']
                elif ".ped" in file.swaggerTypes['name']:
                    self.inputs['inputPED'] = file.swaggerTypes['path']
                elif ".map" in file.swaggerTypes['name']:
                    self.inputs['inputMAP'] = file.swaggerTypes['path']

        # If the input data was declared as BIM/BED/FAM
        elif self.input_format == 'b':
            for file in file_list:
                if ".ote" in file.swaggerTypes['name']:
                    self.known_truth = file.swaggerTypes['path']
                elif '.bed' in file.swaggerTypes['name']:
                    self.inputs['inputBED'] = file.swaggerTypes['path']
                elif '.bim' in file.swaggerTypes['name']:
                    self.inputs['inputBIM'] = file.swaggerTypes['path']
                elif 'fam' in file.swaggerTypes['name']:
                    self.inputs['inputFAM'] = file.swaggerTypes['path']

        # If the input data was declared as TPED/TMAP
        else:
            for file in file_list:
                if ".ote" in file.swaggerTypes['name']:
                    self.known_truth = file.swaggerTypes['path']
                elif ".tped" in file.swaggerTypes['name']:
                    self.inputs['inputTPED'] = file.swaggerTypes['path']
                elif ".tmap" in file.swaggerTypes['name']:
                    self.inputs['inputTMAP'] = file.swaggerTypes['path']

    def build_jsons(self):
        """Builds JSONs from the parsed input information.

        If there are no finished GWAS jobs, this method will default to building
        JSONs for Winnow from the GWAS output files.

        :return:
        """
        if not self.finished_gwas:
            self.gwas_jsons = []
            for gwas in self.desired_gwas:
                if gwas:
                    self.gwas_jsons += JsonBuilder.make_gwas_json(
                        self.desired_gwas.index(gwas), self.dataset_name, self.inputs)
        else:
            # Make Winnow JSONs here.
            self.winnow_jsons = []
            for gwas_output in self.finished_gwas:
                pass

    def gwas_submission(self):
        """Submits all provided JSON files via the Agave REST APIs.
        All currently running jobs are stored in the 'running_jobs' dictionary,
        which has the format of { 'job[id]' : 'Job' } with 'Job' being the
        Job.py response from the server given upon job submission.

        :return:
        """
        for json in self.gwas_json:
            job = Agave.JobsApi.submit(json).swaggerTypes['result']
            self.running_jobs[job.swaggerTypes['id']] = job.swaggerTypes['archivePath']

    def poll_jobs(self):
        """Polls all running jobs for status until completion.

        Running jobs are stored in the running_jobs dictionary. Once finished,
        jobs are removed and added to the finished_gwas dictionary for easy
        tracking and handling for running through Winnow.
        """
        #TODO look at notifications instead?
        #TODO keep data on the datastore - no downloading
        sleep_time = 60

        # Iterating through all JobIDs and polling until there are no more jobs
        while not self.running_jobs.keys():
            for job_id in self.running_jobs.keys:
                job_status = Agave.JobsApi.getStatus(job_id)

                # If the curernt given job is finished, download output and
                # remove it from the running queue. The finished job and its
                # output list is added to the finished_gwas dictionary.
                if ((job_status.swaggerTypes['result'].swaggerTypes['status'] == "FINISHED") or
                        (job_status.swaggerTypes['status'] == "FINISHED")):
                    # TODO use archive instead of downloading
                    # agave://data.iplantcollaborative.org/<user-home>/archive/jobs/job-<jobID>
                    self.finished_gwas[job_id] = self.running_jobs[job_id]
                    del self.running_jobs[job_id]

            # Sleep before repolling Agave. Max sleep time is 1 hour.
            time.sleep(sleep_time)
            sleep_time *= 2 if sleep_time <= 3600 else sleep_time

    def download_outputs(self, job_id):
        # TODO avoid - work via Agave & Discovery Environment instead
        # Make a new directory for the current job's output
        if not os.path.exists():
            os.makedirs(job_id)
        os.chdir(job_id)

        # Getting the list of outputs for a given job and downloading
        self.finished_gwas[job_id] = Agave.JobsApi.listOutputs(job_id, self.running_jobs[job_id])
        Agave.JobsApi.downloadOutput(job_id)
        del self.running_jobs[job_id]

        os.chdir("../")

    def parse_archives(self, jobid_dict):
        """Parses job output archives on the Discovery Environment via Agave.

        :param jobid_dict: Dictionary with the Agave job-id as keys.
        :return:
        """
        out_extensions = ['.out.txt', '.freq', '.gv', '.hyp', '.log', '.model',
                          '.param']
        job_output_dict = {}

        # Iterates thorugh all Job-IDs and collects the output type jobs.
        for jobid in jobid_dict.keys:
            job_outputs = []
            remote_file_list = Agave.FilesApi.listOnDefaultSystem(
                filePath="/archive/jobs/job-{}".format(jobid)).swaggerTypes['result']

            for file in remote_file_list:
                if file.swaggerTypes['format'] in out_extensions:
                    job_outputs.append(file.swaggerTypes['path'])

            job_output_dict[jobid] = job_outputs

        return job_output_dict

    def make_output_folders(self, job_outputs):
        """Makes folders of all outputs for each finished GWAS job.

        :param job_outputs: Dictionary of { 'jobid' : list(output_paths) }
        :return:
        """
        output_folders = {}
        Agave.FilesApi.manageOnDefaultSystem(sourcefilePath='.', action='mkdir',
                                             filePath='Validate GWAS Outputs')

        # Loops through all finished Job IDs and creates a subdirectory within
        # the 'Validate GWAS Outputs' folder defined above simply named after
        # each Job ID. All outputs are stored here for easy access.
        for jobid in job_outputs.keys:
            Agave.FilesApi.manageOnDefaultSystem(
                sourcefilePath='./Validate GWAS Outputs', action='mkdir',
                filePath=jobid)

            # Moving all job outputs from the system archive directory to the
            # newly created subdirectory.
            for file in job_outputs[jobid]:
                Agave.FilesApi.manageOnDefaultSystem(
                    sourcefilePath='./Validate GWAS Outputs/{}'.format(jobid),
                    action='cp', filePath=file)

    def create_winnow_jsons(self):
        # Files are located in /iplant/home/<user-dir>/archive/jobs/job-<jobid>
        for jobid, remotelist in self.finished_gwas:
            JsonBuilder.make_winnow_json(jobid, gwas_output_folder=, ote_file=self.known_truth)

    def winnow_submission(self):
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
        pass

    def demonstrate_submission(self):
        pass


if __name__ == '__main__':
    Pipeline()
