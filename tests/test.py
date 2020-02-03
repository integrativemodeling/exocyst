import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))


class Tests(unittest.TestCase):


    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'production_script'))
        p = subprocess.check_call(["python", "exocyst.py", "--test"])
        
    

if __name__ == '__main__':
    unittest.main()
