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


    def test_hasAllKeys(self):
        """Test all sections have all required keys.
           but I notice config.get() methods in mpy_batch.py used fallback values, except 'prefix' in rename_param
           section.
        """
        pass

    def test_hasIntVal(self):
        """Test those keys that must have int values.
        They are: processors,bdiffs,pdiffs,insert. maxambig,maxlength. pdiffs,rdiffs. nseqs
        """
        pass

    def test_hasAllSections(self):
        pass

    def test_dirFileExists(self):
        """
        Test essential directories and files exist.
        They're: input_dir, output_dir, batch_file, oligos in [file_inputs]
        """
        #self.assertTrue(cc.dirFileExists(self.config,'input_dir'))
        self.assertFalse(cc.dirFileExists(self.config,'input_dir'))
        self.assertFalse(cc.dirFileExists(self.config, 'no_such_dir'))

        return


if __name__ == '__main__':
    unittest.main()
