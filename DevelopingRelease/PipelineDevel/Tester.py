import json
import JsonBuilder
import unittest

class Tester:

    def __init__(self):
        self.dataset_name = 'dongwong'
        self.desired_gwas = (True, None, None, None, None, None)
        self.inputs = {
            "inputMAP": "https://agave.iplantc.org/files/v2/media/system/data.iplantcollaborative.org//mjs3607/fastlmm_data/dongwang.map",
            "inputPED": "https://agave.iplantc.org/files/v2/media/system/data.iplantcollaborative.org//mjs3607/fastlmm_data/dongwang.ped"
        }


    def json_test(self):
        self.gwas_jsons = []

        print "\n\n"
        for k in self.inputs.keys():
            print json.dumps(self.inputs)
        print "\n\n"

        for gwas in self.desired_gwas:
            if gwas:
                self.gwas_jsons.append(JsonBuilder.make_gwas_json(
                    self.desired_gwas.index(gwas), self.dataset_name, self.inputs))
                # else:
                #     Make Winnow JSONs here.
                # self.winnow_jsons = []
                # for gwas_output in self.finished_gwas:
                #     pass

        # TODO delete print
        for js in self.gwas_jsons:
            print js


if __name__ == '__main__':
    tester = Tester()
    tester.json_test()
