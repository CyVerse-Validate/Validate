import argparse


class CommandLine:
    def __init__(self):
        self.init_graphics()
        self.args = self.check_args()

    @staticmethod
    def init_graphics():
        print "###################################################################"
        print "###                                                            ####"
        print "###                 File Merger for Validate                   ####"
        print "###                                                            ####"
        print "###################################################################"

    @staticmethod
    def usage():
        options = ['launcher', 'bayesr', 'alphasim']
        print ''
        print 'Command-line usage help menu'
        print '--help or -h to see the help menu'
        print '--verbose or -v for verbose mode'
        print '--output or -o to specify the output file prefix'
        print ', '.join(options), 'to specify the mode'
        print '\nLauncher Arguments'
        print '\t--folder or -f to specify the folder containing GWAS outputs'
        print 'BayesR Arguments'
        print '\t--bim or -b to specify the BIM file used in BayesR'
        print '\t--param or -p to specify the .param output file from BayesR'
        print 'AlphaSim Arguments'
        print '\t--snp or -s to specify the SNP information output file from AlphaSim'
        print '\nExample Usage: '
        print 'python merger --output mergedOut bayesr --bim sim.bim --param sim.param'

    def check_args(self):
        parser = argparse.ArgumentParser(description='File merger for Validate apps')
        parser.add_argument('-v', '--verbose', help='Trigger verbose mode', action='store_true', default=False)
        parser.add_argument('-o', '--output', help='Prefix for output files', type=str, default='merged_output')

        subparsers = parser.add_subparsers(help='Merge mode', dest='mode')
        self.add_options(subparsers)
        args = parser.parse_args()
        return args

    @property
    def args(self):
        return self.args

    def add_options(self, sub):
        self.add_alpha_sim_options(sub)
        self.add_bayes_r_options(sub)
        self.add_launcher_options(sub)

    @staticmethod
    def add_launcher_options(sub):
        launcher_parser = sub.add_parser('launcher')
        launcher_parser.add_argument('-f', '--folder', help='The folder containing all GWAS outputs', required=True)

    @staticmethod
    def add_bayes_r_options(sub):
        bayes_r_parser = sub.add_parser('bayesr')
        bayes_r_parser.add_argument('-b', '--bim', help='BIM file used in BayesR (.bim)', required=True)
        bayes_r_parser.add_argument('-p', '--param', help='BayesR param output file (.param)', required=True)

    @staticmethod
    def add_alpha_sim_options(sub):
        alpha_sim_parser = sub.add_parser('alphasim')
        alpha_sim_parser.add_argument('-s', '--snp', help='SNP information file from AlphaSim', required=True)
