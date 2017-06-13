from agavepy.agave import Agave

__author__ = "Michael J. Suggs // mjs3607@uncw.edu"


class Pipeline:
    def __init__(self, username, password):
        agave_tenant = Agave(api_server='https://agave.iplantc.org',
                             username=username, password=password)
