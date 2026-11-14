import unittest
import numpy as np
from src.retune_biological_models import stopping_mask


class StoppingSplitTests(unittest.TestCase):
    def test_same_locus_never_crosses_stopping_split(self):
        groups=np.repeat([f'1:{i}' for i in range(1000)],2)
        mask=stopping_mask(groups)
        self.assertTrue(mask.any())
        self.assertTrue((~mask).any())
        self.assertFalse(set(groups[mask]) & set(groups[~mask]))
        np.testing.assert_array_equal(mask,stopping_mask(groups))


if __name__=='__main__':
    unittest.main()
