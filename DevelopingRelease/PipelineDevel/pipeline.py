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

# TODO: Move pipeline to a new directory to contain all outputs
class Pipeline:
    """Handles job submission for the Validate Pipeline."""
    # TODO Change pipeline around to take user folder location

    def __init__(self):
        self.data_folder = ""       # User defined folder from the Datastore
        self.inputs = {}
        self.known_truth = ""
        self.dataset_name = ""
        self.desired_gwas = ()      # Tuple of booleans
        self.running_jobs = {}      # Running jobs dictionary with the format { 'id' : 'outputPath' }
        self.finished_jobs = {}     # Finished jobs dictionary with the format { 'id' : list[RemoteFile] }

        self.checkArgs()

        # TODO build gwas files from the given folder.
        self.gwas_submission()
        self.poll_jobs()

    def checkArgs(self):
        """Checks all command-line arguments provided in the command-line call.

        :return:
        """
        parser = argparse.ArgumentParser()
        # parser.add_argument("-g", "--gwas", help="GWAS JSON files for Agave", required=True)
        # parser.add_argument("-t", "--OTE", help="OTE file location on the Discovery Environment", required=True)
        parser.add_argument("-i", "--InFormat", type=chr(1),
                            help="Input format:\n"
                                 "\tp for PED/MAP\n"
                                 "\tb for BIM/BED/FAM\n"
                                 "\tt for TPED/TFAM")
        parser.add_argument("-f", "--Folder", type=str,
                            help="Folder to be pipelined.")
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
        # self.gwas_json = args.gwas
        # self.ote_location = args.ote
        self.input_format = args.InFormat
        self.data_folder = args.Folder
        self.desired_gwas = tuple([args.lmm, args.ridge, args.bayes, args.plink,
                                   args.qxpak, args.gemma])

    def parse_inputs(self):
        file_list = Agave.FilesApi.listOnDefaultSystem(self.data_folder).result

        if self.input_format == 'p':
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = file.path
                elif ".ped" in file.name:
                    self.inputs['inputPED'] = file.path
                elif ".map" in file.name:
                    self.inputs['inputMAP'] = file.path
        elif self.input_format == 'b':
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = file.path
                elif '.bed' in file.name:
                    self.inputs['inputBED'] = file.path
                elif '.bim' in file.name:
                    self.inputs['inputBIM'] = file.path
                elif 'fam' in file.name:
                    self.inputs['inputFAM'] = file.path
        else:
            for file in file_list:
                if ".ote" in file.name:
                    self.known_truth = file.path
                elif ".tped" in file.name:
                    self.inputs['inputTPED'] = file.path
                elif ".tmap" in file.name:
                    self.inputs['inputTMAP'] = file.path

        # for file in file_list:
        #     if ".ote" in file.name:
        #         self.known_truth = file.path
        #     elif self.input_format == 'p':
        #         self.inputs['inputPED'] =
        #         self.inputs['inputMAP'] =
        #     elif self.input_format == 'b':
        #         self.inputs['inputBED'] =
        #         self.inputs['inputBIM'] =
        #         self.inputs['inputFAM'] =
        #     elif self.input_format == 't':
        #         self.inputs['inputTPED'] =
        #         self.inputs['inputTFAM'] =
            # elif any([x in file.name for x in [".ped", ".map"]]):
            #     self.dataset_name, tmp = file.split(".")
            #     self.inputs.append(file)
            # elif any([x in file.name for x in [".bim", ".bam", ".fam"]]):
            #     self.dataset_name, tmp = file.split(".")
            #     self.inputs.append(file)
            # elif any([x in file.name for x in [".tped", ".tfam"]]):
            #     self.dataset_name, tmp = file.split(".")
            #     self.inputs.append(file)

        # for file in file_list:
        #     if file.type.lower() == "file":
        #         name, ext = file.name.split(".")
        #         if ext.lower() == ".ote":
        #             self.known_truth = file
        #         elif ext

    def build_jsons(self):
        self.gwas_jsons = []
        for gwas in self.desired_gwas:
            if gwas:
                self.gwas_jsons += JsonBuilder.make_gwas_json(self.desired_gwas.index(gwas),
                                                              self.data_folder)

    def gwas_submission(self):
        """Submits all provided JSON files via the Agave REST APIs.

        :return:
        """
        for json in self.gwas_json:
            job = Agave.JobsApi.submit(json).result
            self.running_jobs[job.id] = job.outputPath

    def poll_jobs(self):
        """Polls all running jobs for status until completion.
        Jobs are removed once they are finished.
        """
        #TODO look at notifications instead?
        sleep_time = 60
        os.makedirs("Outputs")
        os.chdir("Outputs")

        # Iterating through all JobIDs and polling until there are no more jobs
        while not self.running_jobs.keys():
            for job_id in self.running_jobs:
                job_status = Agave.JobsApi.getStatus(job_id)

                # If the curernt given job is finished, download output and
                # remove it from the running queue. The finished job and its
                # output list is added to the finished_jobs dictionary.
                if (job_status.result.status == "FINISHED") or (job_status.status == "FINISHED"):
                    # TODO use archive instead of downloading
                    # agave://data.iplantcollaborative.org/<user-home>/archive/jobs/job-<jobID>
                    self.download_outputs(job_id)

            # Sleep before repolling Agave. Max sleep time is 1 hour.
            time.sleep(sleep_time)
            sleep_time *= 2 if sleep_time <= 3600 else sleep_time

        os.chdir(self.submission_dir)

    def download_outputs(self, job_id):
        #TODO avoid - work via Agave & Discovery Environment instead
        # Make a new directory for the current job's output
        if not os.path.exists():
            os.makedirs(job_id)
        os.chdir(job_id)

        # Getting the list of outputs for a given job and downloading
        self.finished_jobs[job_id] = Agave.JobsApi.listOutputs(job_id, self.running_jobs[job_id])
        Agave.JobsApi.downloadOutput(job_id)
        del self.running_jobs[job_id]

        os.chdir("../")


    def create_winnow_jsons(self):
        # Files are located in /iplant/home/<user-dir>/archive/jobs/job-<jobid>
        for jobid, remotelist in self.finished_jobs:
            JsonBuilder.make_winnow_json(jobid, gwas_folder=, ote_file=)

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
        # Gemma
        pass

    def demonstrate_submission(self):
        pass


if __name__ == '__main__':
    Pipeline()
