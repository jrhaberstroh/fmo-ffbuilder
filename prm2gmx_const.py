import yaml

twoletter = { 'cab' : 'CD', 'cbb' : 'CE', 'cnb' : 'CF', 'cpb' : 'CG', 'cqb' : 'CH', 'crb' : 'CI', 'csb' : 'CJ', 'ccs' : 'CL', 'cqq' : 'CO', 'c2k' : 'CP', 'c2a' : 'C5', 'c2e' : 'CU', 'cq2' : 'CX', 'ct1' : 'C1', 'ct2' : 'C2', 'ct3' : 'C3', 'nmh' : 'NM', 'mgc' : 'MC', 'ha0' : 'HQ', 'o2c' : 'OC', 'o1c' : 'O1', 'rb'  : 'Rb', 'br'  : 'Br', 'cl'  : 'Cl', 'na'  : 'Na', 'cs'  : 'Cs', 'li'  : 'Li', 'zn'  : 'Zn'};

file_handle = open("dat/amber2charmm.yml", 'r')
amber2charmm = yaml.safe_load(file_handle)
file_handle.close()

# Conversion for atomtypes, from three-letter into two-letter. I don't know why she swallowed
#   the fly; perhaps she'll die.
# atomnames remain arbitrary numbers of letters
#
#
#     'c2a' used to convert to CS, but this was changed to C5 for fear of overlap with Cs

mass2num = { "1.008"  : 1, \
	     "6.940"  : 3, \
	     "12.010" : 6, \
	     "12.011" : 6, \
	     "14.010" : 7, \
	     "15.999" : 8, \
	     "16.000" : 8, \
	     "16.00000" : 8, \
	     "19.000" : 9, \
	     "22.990" : 11,\
	     "24.305" : 12,\
	     "30.970" : 15,\
	     '32.060' : 16,\
	     '35.450' : 17,\
	     "39.100" : 19,\
	     "40.080" : 20,\
	     "55.847" : 26,\
	     "63.550" : 29,\
	     "85.470" : 37,\
	     "126.900": 53,\
	     "131.000": 54,\
	     "132.910": 55 };

atom2mass = { "H" : "1.008",  \
              "C" : "12.010", \
              "N" : "14.010", \
              "O" : "16.000", \
              "M" : "24.305"}


rtp_header = """; Header taken from amber03 aminoacids.rtp
[ bondedtypes ]
; Column 1 : default bondtype
; Column 2 : default angletype
; Column 3 : default proper dihedraltype
; Column 4 : default improper dihedraltype
; Column 5 : This controls the generation of dihedrals from the bonding.
;            All possible dihedrals are generated automatically. A value of
;            1 here means that all these are retained. A value of
;            0 here requires generated dihedrals be removed if
;              * there are any dihedrals on the same central atoms
;                specified in the residue topology, or
;              * there are other identical generated dihedrals
;                sharing the same central atoms, or
;              * there are other generated dihedrals sharing the
;                same central bond that have fewer hydrogen atoms
; Column 6 : number of neighbors to exclude from non-bonded interactions
; Column 7 : 1 = generate 1,4 interactions between pairs of hydrogen atoms
;            0 = do not generate such
; Column 8 : 1 = remove proper dihedrals if found centered on the same
;                bond as an improper dihedral
;            0 = do not generate such



; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih
     1       1          9          4        1         3      1     0\n\n\n"""
