import unittest
import prm2gmx

class CoreUnittests(unittest.TestCase):

    def test_convertBondline(self):
        out = prm2gmx.bondline("nmh   mgc    50.00   2.025")
        #print out
        self.assertEqual(out, "NM MC 2.025 20920")
    def test_vstrongBond(self):
        out = prm2gmx.bondline("cqq   ct3    270.0    1.510", precision = 7)
        self.assertEqual(out, "CO C3 1.51 112968")
    def test_noatomBond(self):
        out = prm2gmx.bondline("cq2   cqo    250.0    1.495", precision = 7)
        self.assertEqual(out, None)
    def test_checkManageBond(self):
        out = prm2gmx.managelines("nmh   mgc    50.00   2.025", "BOND")
        #print out
        self.assertEqual(out, "NM MC 2.025 20920")

    def test_convertBendline(self):
        out = prm2gmx.bendline("nmh   mgc   nmh       50.0    178.9")
        #print out
        self.assertEqual(out, "NM MC NM 178.9 209.2")

    def test_convertTorsionline(self):
        out = prm2gmx.torsionline("x    fe   nb   x      0.00      4")
        #print out
        self.assertEqual(out, "x fe nb x 180.0 0.0 4")

    def test_torsionNegativeDefault(self):
        out = prm2gmx.torsionline("x    ct1  ct2  x       -0.156       3")
        self.assertEqual(out, "x C1 C2 x 0.0 0.652 3")
    def test_torsionNegativeNondefault(self):
        out = prm2gmx.torsionline("cqq  cq2  o1c  ct3     -0.632       1   180.0")
        self.assertEqual(out, "FAIL")

    def test_convertImproperline(self):
        out = prm2gmx.improperline("x    x    crb   cab     1.100     2.0  180.0   cosine")
        #print out
        self.assertEqual(out, "x x CI CD 180.0 4.602 2")


    def test_checkManageEnd(self):
        mode = "hooray!"
        out = prm2gmx.managelines("END\n", mode)
        self.assertEqual(out, None)
        self.assertEqual(mode, "")


if __name__ == "__main__":
    unittest.main()
