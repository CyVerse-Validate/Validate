import argparse, os, time
from AgavePythonSDK import JobsApi

__author__ = "Michael J. Suggs::mjs3607@uncw.edu"
__credits__ = ["Michael Suggs"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Michael Suggs"
__email__ = "mjs3607@uncw.edu"
__status__ = "Development"

# TODO: Move pipeline to a new directory to contain all outputs
class Pipeline:
    def __init__(self):
        # Running jobs dictionary with the format { 'id' : 'outputPath' }
        self.running_jobs = {}

        # Finished jobs dictionary with the format { 'id' : list[RemoteFile] }
        # The remote file list gives the relative path of output files.
        self.finished_jobs = {}
        self.submission_dir = os.getcwd()

        self.checkArgs()
        self.gwas_submission()
        self.poll_jobs()

    def checkArgs(self):
        """Checks all command-line arguments provided in the command-line call.

        :return:
        """
        parser = argparse.ArgumentParser()
        parser.add_argument("-g", "--gwas", help="GWAS JSON files for Agave", required=True)
        args = parser.parse_args()
        self.gwas_json = args.gwas

    def gwas_submission(self):
        """Submits all provided JSON files via the Agave REST APIs.

        :return:
        """
        for json in self.gwas_json:
            job = JobsApi.submit(json).result
            self.running_jobs[job.id] = job.outputPath

    def poll_jobs(self):
        """Polls all running jobs for status until completion.
        Jobs are removed once they are finished.
        """
        sleep_time = 60
        os.makedirs("Outputs")
        os.chdir("Outputs")

        # Iterating through all JobIDs and polling until there are no more jobs
        while not self.running_jobs.keys():
            for job_id in self.running_jobs:
                job_status = JobsApi.getStatus(job_id)

                # If the curernt given job is finished, download output and
                # remove it from the running queue. The finished job and its
                # output list is added to the finished_jobs dictionary.
                if (job_status.result.status == "FINISHED") or (job_status.status == "FINISHED"):
                    self.download_outputs(job_id)

            # Sleep before repolling Agave. Max sleep time is 1 hour.
            time.sleep(sleep_time)
            sleep_time *= 2 if sleep_time <= 3600 else sleep_time

        os.chdir(self.submission_dir)

    def download_outputs(self, job_id):
        # Make a new directory for the current job's output
        if not os.path.exists():
            os.makedirs(job_id)
        os.chdir(job_id)

        # Getting the list of outputs for a given job and downloading
        self.finished_jobs[job_id] = JobsApi.listOutputs(job_id, self.running_jobs[job_id])
        JobsApi.downloadOutput(job_id)
        del self.running_jobs[job_id]

        os.chdir("../")


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

    def create_winnow_jsons(self):
        pass

    def demonstrate_submission(self):
        pass


if __name__ == '__main__':
    Pipeline()
