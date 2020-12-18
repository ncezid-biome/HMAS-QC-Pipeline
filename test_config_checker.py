import unittest
from configparser import ConfigParser
import os
import config_checker as cc

class TestConfig_checker(unittest.TestCase):

    def setUp(self):
        """ set up test variables"""
        cfg_file = os.getcwd() + "\\settings.ini"  # assume the init file is in current directory

        self.config = ConfigParser()
        self.config.read(cfg_file)

        return


    def test_hasAllOptions(self):
        """Test all sections have all required options.
           but I notice config.get() methods in mpy_batch.py used fallback values, except 'prefix' in rename_param
           section.
        """
        pass

    def test_hasIntVal(self):
        """Test those options that must have int values.
        They are: processors,bdiffs,pdiffs,insert. maxambig,maxlength. pdiffs,rdiffs. nseqs
        param_dict = {'contigs_params': ['processors', 'bdiffs', 'pdiffs', 'insert'],
                      'screen_params': ['maxambig', 'maxlength'],
                      'pcr_params': ['pdiffs', 'rdiffs'],
                      'rare_seqs_param': ['nseqs']}
        'init' file is assumed in the current directory.
        """
        self.assertTrue(cc.hasIntVal(self.config,'contigs_params','insert'))
        #self.assertFalse(cc.hasIntVal(self.config, 'contigs_params', 'insert'))

        self.assertTrue(cc.hasIntVal(self.config,'rare_seqs_param','nseqs'))
        #self.assertFalse(cc.hasIntVal(self.config, 'rare_seqs_param', 'nseqs'))

        return

    def test_hasAllSections(self):
        pass

    def test_dirFileExists(self):
        """
        Test essential directories and files options appear in the init file and exist .
        They're: input_dir, output_dir, batch_file, oligos in [file_inputs]
        'init' file is assumed in the current directory.
        """
        #self.assertTrue(cc.dirFileExists(self.config,'file_inputs','input_dir'))
        self.assertFalse(cc.dirFileExists(self.config,'file_inputs','input_dir'))
        self.assertFalse(cc.dirFileExists(self.config,'file_inputs','no_such_dir'))

        return


if __name__ == '__main__':
    unittest.main()
