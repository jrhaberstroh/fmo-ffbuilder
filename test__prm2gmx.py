import unittest
import prm2gmx

class prmUnittests(unittest.TestCase):
    def setUp(self):
        prm2gmx.niceformat=False
        prm2gmx.convertToTwo=True
        prm2gmx.convertToCharmm=False
        prm2gmx.addSuffix=False
        prm2gmx.VDW_scale = 1.0
        prm2gmx.kCal2kJ = 4.184
        prm2gmx.bond_scale = 1.0
        prm2gmx.angle_scale = 1.0
        prm2gmx.dihedral_scale = 1.0

    def test_convertBondline(self):
        out = prm2gmx.bondline("nmh   mgc    50.00   2.025", precision = 6)
        #print out
        self.assertEqual(out, "NM MC 1 0.2025 20920.")
    def test_vstrongBond(self):
        out = prm2gmx.bondline("cqq   ct3    270.0    1.510", precision = 6)
        self.assertEqual(out, "CO C3 1 0.151 112968")
    # Specification changed; test deprecated. Return expected is now simply
    # the line, but to-upper
    #def test_noatomBond(self):
    #    out = prm2gmx.bondline("cq2   cqo    250.0    1.495", precision = 6)
    #    self.assertEqual(out, None)
    def test_checkManageBond(self):
        out, trash = prm2gmx.manageprmlines("nmh   mgc    50.00   2.025", "BOND", precision=6)
        #print out
        self.assertEqual(out, "NM MC 1 0.2025 20920.")

    def test_convertBendline(self):
        out = prm2gmx.bendline("nmh   mgc   nmh       50.0    178.9")
        #print out
        self.assertEqual(out, "NM MC NM 1 178.9 209.2")

    def test_convertTorsionline(self):
        out = prm2gmx.torsionline("x    fe   nb   x      0.00      4")
        #print out
        self.assertEqual(out, "X FE NB X 9 180.0 0.0 4")

    def test_torsionNegativeDefault(self):
        out = prm2gmx.torsionline("x    ct1  ct2  x       -0.156       3")
        self.assertEqual(out, "X C1 C2 X 9 0.0 0.6527 3")
# 2015-03-02: Removed becase feature was already removed from prm2gmx.py 
#             without explanation. Forcefield must conform with this
#             convention already?
#    def test_torsionNegativeNondefault(self):
#        out = prm2gmx.torsionline("cqq  cq2  o1c  ct3     -0.632       1   180.0")
#        self.assertEqual(out, "FAIL")

    def test_convertImproperline(self):
        out = prm2gmx.improperline("x    x    crb   cab     1.100     2.0  180.0   cosine", precision=5)
        #print out
        self.assertEqual(out, "X X CI CD 4 180.0 4.602 2")

    def test_convertNonbondedline(self):
        out = prm2gmx.nonbondedline("hs      0.600   0.016 0.000  0.000   1.008")
        self.assertEqual(out, "HS 1 1.008 0.0 A 0.0534 0.0669")
    
    def test_checkRounding(self):
        out = prm2gmx.nonbondedline("hs      0.600   0.016 0.000  0.000   1.008")
        self.assertEqual(out, "HS 1 1.008 0.0 A 0.0534 0.0669")

    def test_checkManageEnd(self):
        mode = "hooray!"
        out, mode = prm2gmx.manageprmlines("END\n", mode)
        self.assertEqual(out, None)
        self.assertEqual(mode, None)



class tpgUnitTests(unittest.TestCase):
    def setUp(self):
        lines = "RESIDUE  bcl ( Total Charge =   0.0 ) \n atoms \n group \n mg     mgc     0.1340 \n group  \n cha    csb     0.0650    \n group    \n chb    cab    -0.2240                                              \n hb     ha0     0.1150\n"
        self.lines = lines.split('\n')
        self.lines_bond = " bonds                                                              \nnb   mg \nnc   mg \nnd   mg \nna   mg \nnb   c1b".split('\n')
        prm2gmx.convertToTwo = True
        prm2gmx.convertToCharmm=False
        prm2gmx.convertToBCL=False

    def test_genericChargeLine(self):
        self.lines[3]
        out = prm2gmx.generic_line(self.lines[3].split(), 2, [0], [1.], [float])
        self.assertEqual(out, "MG MC 0.134")
    def test_chargeLine(self):
        self.lines[3]
        out = prm2gmx.chargeline(self.lines[3], 1)
        self.assertEqual(out, "MG MC 0.134 1")

    def test_manageLineSimple(self):
        mode = "atoms"
        group = 2
        test_lines = [None, "MG MC 0.134 3", None, "CA CB 0.0650"]
        test_modes = ["atoms"] * 4
        test_group = [3, 3, 4, 4]
        for i,line in enumerate(self.lines[2:5]):
            newline, mode, group = prm2gmx.managetpglines(line, mode, group)
            self.assertEqual(test_lines[i], newline)
            self.assertEqual(test_modes[i], mode)
            self.assertEqual(test_group[i], group)
    def test_manageLineMode(self):
        mode = None
        group = 0
        test_lines = ['[BCL]', ' [ atoms ]', None, "MG MC 0.134 1", None, "CA CB 0.0650   2"]
        test_modes = ["RESIDUE"] + ["atoms"] * 4
        test_group = [0, 0, 1, 1, 2, 2]
        for i,line in enumerate(self.lines[0:5]):
            newline, mode, group = prm2gmx.managetpglines(line, mode, group)
            self.assertEqual(test_lines[i], newline)
            self.assertEqual(test_modes[i], mode)
            self.assertEqual(test_group[i], group)
    def test_TLA_NegativeIn(self):
        bad_in = "cma    ct3    -0.2750\n"
        out = prm2gmx.chargeline(bad_in,1, precision=6)
        self.assertEqual(out, "CMA C3 -0.275 1")
        out = prm2gmx.chargeline(bad_in,1, precision=5)
        self.assertNotEqual(out, "CMA C3 -0.275 1")
        


    def test_manageBlank(self):
        newline,mode,group= prm2gmx.managetpglines("",None,0)
        self.assertEqual(newline, None)

    def test_manageBonds(self):
        mode = None
        group = 0
        test_lines = [' [ bonds ]', 'NB MG']
        test_modes = ["bonds"]*2
        test_group = [None, None]
        for i,line in enumerate(self.lines_bond[0:1]):
            newline, mode, group = prm2gmx.managetpglines(line, mode, group)
            self.assertEqual(test_lines[i], newline)
            self.assertEqual(test_modes[i], mode)
            self.assertEqual(test_group[i], group)

    def test_manageSwitchMode(self):
        mode = "Hella Cool"
        group = 0
        test_lines = [None,' [ bonds ]']
        test_modes = [None, 'bonds']
        test_group = [None, None]
        endline =["end"]
        endline.append(self.lines_bond[0])
        for i,line in enumerate(endline[0:1]):
            newline, mode, group = prm2gmx.managetpglines(line, mode, group)
            self.assertEqual(test_lines[i], newline)
            self.assertEqual(test_modes[i], mode)
            self.assertEqual(test_group[i], group)


class toCharmmUnitTests(unittest.TestCase):
    def setUp(self):
        prm2gmx.convertToTwo = False
        prm2gmx.niceformat=False
        prm2gmx.convertToCharmm= True

    def test_atomtypeOnly(self):
        line = "hc     ha0     0.1680   "
        out = prm2gmx.chargeline(line, 1, precision=6)
        self.assertEqual(out, "HC HA0 0.168 1")



if __name__ == "__main__":
    unittest.main()
