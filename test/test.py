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
        
    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'archiving'))
        if os.path.exists("exocyst.cif"):
            os.unlink("exocyst.cif")
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
                ["python", "exocyst_mmcif.py", "--mmcif", "--dry-run"], env=env)
        # Check output file
        self._check_mmcif_file('exocyst.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 4)
        self.assertEqual(s.citations[0].doi, '10.1002/pro.3863')
        self.assertEqual(len(s.software), 3)
        self.assertEqual(len(s.orphan_starting_models), 11)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 1 model
        self.assertEqual(sum(len(x) for x in state1), 1)
        # Check # of spheres and atoms in each model
        m = state1[0][0]
        self.assertEqual(len(m._spheres), 3382)
        self.assertEqual(len(m._atoms), 0)
        # Should be 1 ensemble
        self.assertEqual([e.num_models for e in s.ensembles], [9669])
        # Two sets of crosslinks, and an EM map restraint
        xl1, xl2, em = s.restraints
        self.assertEqual(len(xl1.experimental_cross_links), 256)
        self.assertEqual(len(xl1.cross_links), 256)
        self.assertEqual(xl1.linker.auth_name, 'DSS')
        self.assertEqual(len(xl2.experimental_cross_links), 178)
        self.assertEqual(len(xl2.cross_links), 178)
        self.assertEqual(xl2.linker.auth_name, 'EDC')
        self.assertEqual(em.dataset.location.path,
                         'data/emd_21226.map.mrc.gmm.200.txt')
        self.assertEqual(em.dataset.parents[0].location.access_code,
                         'EMD-21226')


if __name__ == '__main__':
    unittest.main()
