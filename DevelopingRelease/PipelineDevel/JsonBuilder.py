import json

__author__ = "Michael J. Suggs // mjs3607@uncw.edu"

class JsonBuilder:
    """Summary of class here.

    Longer class information....
    Longer class information....

    Attributes:
        likes_spam: A boolean indicating if we like SPAM or not.
        eggs: An integer count of the eggs we have laid.
    """

    def __init__(self):
        pass

    def make_gwas_json(self, type):
        """Makes an Agave-API compatible JSON for running jobs.

        :param type: The type of GWAS JSON to generate.
        :return:
        """
        pass

    def fastlmm_json(self):
        pass

    def make_winnow_json(self, gwas_name, gwas_folder, ote_file):
        json_text = {
            "jobName" : gwas_name,
            "softwareName" : "Winnow-0.9.2u1",
            "requestedTime" : "02:00:00",
            "archive" : "true",
            "inputs" : {
                "Class" : ote_file,
                "Folder" : gwas_folder
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

        with open('{}.json'.format(gwas_name), 'w') as f:
            json_text = json.dumps(json_text, f)


if __name__ == '__main__':
    pass
