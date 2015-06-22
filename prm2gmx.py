import sys
import traceback
import numpy as np
import yaml
from prm2gmx_const import *


niceformat = True
convertToTwo = False
convertToCharmm = False

addSuffix = True
atomtypeSuffix = "BC"

kCal2kJ = 4.184
nm2A    = 10.
VDW_scale = 2.0
bond_scale = 2.0
angle_scale = 2.0
dihedral_scale = .5

def set_pos(string_in, pos):
    assert pos >= len(string_in)
    return string_in + " "*(pos - len(string_in))

def generic_line(line_arr, entries_per_line, order, units, datatype=None, precision=6, atomname=None):
    ''' 
    generic_line()

    Converts line_arr into the correct format.



    order determines the order of the quantitative arguments
    units determines their units
    dataype determines the conversion for output (default float)
    precision determines how many chars to print in quantitative args
    
    e.g.
    IN :: [q1, q2, q3]
    order = [1, 2, 0]
    units = [1., 2., 3.]
    OUT :: [1 * q2, 2 * q3, 3 * q1]

    atomname is an optional array of bools, of length entries_per_line. Bools with False will be converted
        to CHARMM names, and bools with True will be skipped in conversion. This protects atomnames but
        allows conversion of atomtypes.

    Returns a string of line_out.
    '''

    if datatype==None:
        datatype=[float]*len(units)
    if atomname==None:
        atomname=[False]*entries_per_line

    assert (len(line_arr) == (entries_per_line + len(order)))
    assert (len(units)    == len(order))
    assert (len(datatype) == len(order))
    assert (len(atomname) == entries_per_line)

    line_out = ""

    column_spacing = 8
    spacer_index = 1

    for i in xrange(entries_per_line) :
        assert(convertToTwo + convertToCharmm + addSuffix == 1)

#        if (twoletter.get(line_arr[i])==None and len(line_arr[i]) >= 3):
#            return None
        if convertToTwo:
            line_out += (line_arr[i].upper() if twoletter.get(line_arr[i])==None else twoletter[line_arr[i]])
        elif addSuffix:
            if atomname[i] == False and line_arr[i].upper() != "X":
                line_out += line_arr[i].upper() + atomtypeSuffix
            else:
                line_out += line_arr[i].upper()
        elif convertToCharmm:
            if atomname[i] == False:
                line_out += amber2charmm[line_arr[i].upper()] if amber2charmm.get(line_arr[i].upper())!=None else line_arr[i].upper()
            else:
                line_out += line_arr[i].upper()
        if niceformat:
            line_out = set_pos(line_out, column_spacing * spacer_index)
            spacer_index += 1
        else:
            line_out += " "

    peak_pos = column_spacing*spacer_index

    column_spacing = 13
    spacer_index = 1
    sign_char = 0

    for i in xrange(len(order)):
        offset = entries_per_line
        if datatype[i] == str:
            val = datatype[i](line_arr[offset +order[i]])
        else:
            val = datatype[i](float(line_arr[offset + order[i]]) * units[i])
            if val < 0:
                sign_char = 1
            if not ((abs(val) < float("1"+"0"*(precision-sign_char)) \
                     and abs(val) >= float("."+"0"*(precision-sign_char)+"1")\
                     ) or val == 0):
                raise RuntimeError("Quantity too large or too small (but non-zero) appeared (has more than "+str(precision)+" chars): "+str(val)+"\n\tline = "+" ".join(line_arr))
        # Check that there are no more than 5 chars to display
        val_full = str((val))
        val = val_full[0:min(precision, len(val_full))]

        if niceformat:
            line_out = set_pos(line_out, peak_pos + column_spacing*spacer_index)
            spacer_index += 1
        line_out += val
        if not niceformat:
            line_out += " "
    return line_out.strip()


def bondline(line_in, precision = 6):
    line_arr = line_in.split() + [1] # Harmonic bond FF type parameter for GROMACS
       #  type  bondlen theta
    order = [2,    1,    0]
    conv =  [1, 1./nm2A, bond_scale * kCal2kJ * nm2A**2]
    types = [int, float, float]
    return generic_line(line_arr, 2, order, conv, types, precision = precision);

def bendline(line_in, precision = 6):
    line_arr = line_in.split() + [1] # Harmonic angle FF type parameter for GROMACS
       #  type   theta    force    
    order = [2,    1,       0]
    conv =  [1.,   1.,  angle_scale * kCal2kJ]
    types = [int,  float, float]
    return generic_line(line_arr, 3, order, conv, types, precision = precision);


def append_phase(line_arr):
    # Append the default angle value if the angle is missing
    # Also, convert the force to a positive number.
    entries_per_line = 4;
    force_pos = 4;

    if float(line_arr[force_pos]) < 0:
        if len(line_arr) == (entries_per_line + 3):
            return False
    if len(line_arr) == (entries_per_line + 3) - 1:
       line_arr.append(0. if (float(line_arr[force_pos]) < 0) else 180.);
    line_arr[force_pos] = str(abs(float(line_arr[force_pos])))
    assert (len(line_arr) == (entries_per_line + 3))
    return True

def torsionline(line_in, precision = 6, improper = False):
    line_arr = line_in.split() 
    append_success = append_phase(line_arr)
    if not append_success:
        return None
    if improper:
        line_arr += [4] # Periodic improper dihedral FF type parameter for GROMACS
    else:
        line_arr += [9] # Periodic proper dihedral FF type parameter for GROMACS

           #  type phase  force    nw
    order =    [3,   2,     0,       1]
    units =    [1,   1.,    dihedral_scale*kCal2kJ, 1]
    datatype = [int, float, float, int]
    entries_per_line = 4;

    # NOTE: Removed failure on torsion lines with 
    # Append the default angle value if the angle is missing
    # Also, convert the force to a positive number.
    #if float(line_arr[force_pos]) < 0:
    #    if len(line_arr) == (entries_per_line + 3):
    #        return "FAIL"
    #if len(line_arr) == (entries_per_line + 3) - 1:
    #   line_arr.append(0. if (float(line_arr[force_pos]) < 0) else 180.);
    #line_arr[force_pos] = str(abs(float(line_arr[force_pos])))
    #assert (len(line_arr) == (entries_per_line + 3))


    return generic_line(line_arr, 4, order, units, datatype, precision = precision)

def improperline(line_in, precision = 6):
    return torsionline(" ".join(line_in.split()[0:-1]), precision = precision, improper = True)

def chargeline(line_in, group, precision = 6):
    return generic_line(line_in.split()+[str(group)], 2, [0, 1], [1., 1.], [float, int], precision=precision, atomname=[True, False])

def nonbondedline(line_in, precision = 6):
    sixth_root_onehalf = .890899
    order =    [2,    3,     4,     5,    0,     1]
    units =    [1.  , 1.,   1.,    1., VDW_scale * sixth_root_onehalf/nm2A, kCal2kJ]
    datatype = [int , float, float, str, float, float]
    dat_arr = line_in.split()
    if mass2num.get(dat_arr[-1]) == None:
        return None
    # Use the first three columns (atomname, R*i [VDW radius], ei [VDW energy]), 
    #   then append the atomic number (using mass lookup table)
    #   then append the mass
    #   then append 0
    #   then append A
    dat_arr = dat_arr[0:3] + [str(mass2num[dat_arr[-1]]), dat_arr[-1], '0.000', 'A']
    return generic_line(dat_arr, 1, order, units, datatype, precision=precision)


def manageprmlines(line, mode, precision = 7):
    line = line.strip()
    newline = None

    if (line == ""):
        pass
    elif (line[0] == "#"):
        newline = "; " + line[1:]
    elif (line == "END"):
        mode = None
    elif (mode == None):
        mode = line
        if   (mode == "BOND"):
            newline = "[bondtypes]\n; A B bond.type  spring  bondlen  (kJ, nm)" 
        elif (mode == "BENDINGS"):
            newline = "[angletypes]\n; A B C bond.type  spring   angle (kJ, nm)"
        elif (mode == "TORSION PROPER"):
            newline = "[dihedraltypes]\n; A B C D  dihed.type  phase  ampl  freq (kJ, nm)"
        elif (mode == "TORSION IMPROPER"):
            newline = "[dihedraltypes]\n; A B C D  dihed.type  phase  ampl  freq (kJ, nm)"
        elif (mode == "NONBONDED MIXRULE"):
            newline = "[ atomtypes ]\n; name at.num mass charge ptype sigma epsilon (kJ, nm)"
        else:
            raise RuntimeError("No acceptable mode active. \n\tmode = "+str(mode))
   
    elif (mode == "BOND"):
        newline = bondline(line, precision)
    elif (mode == "BENDINGS"):
        newline =  bendline(line, precision)
    elif (mode == "TORSION PROPER"):
        newline =  torsionline(line, precision)
    elif (mode == "TORSION IMPROPER"):
        newline =  improperline(line, precision)
    elif (mode == "NONBONDED MIXRULE"):
        newline =  nonbondedline(line, precision)

    else:
        raise RuntimeError("No acceptable active mode or line data.\n\tmode = "+str(mode)+"\n\tline = "+line)

    return newline, mode


def managetpglines(line, mode, group = None, precision = 7):
    if mode == "RESIDUE":
        mode = None
    done = False
    newline = None
    line = line.strip()
    if (line == ""):
        pass
    elif (line[0] == "#"):
        newline = "; " + line[1:]
    elif line == "end":
        mode = None
        group = None
    elif mode == None:
        line_arr = line.split()
        if line_arr[0] == "RESIDUE":
            newline = "[" + line_arr[1].upper() + "]"
            mode = "RESIDUE"
        elif line_arr[0] == "atoms":
            newline = ' [ atoms ]'
            mode = "atoms"
            group = 0
        elif line_arr[0] == "bonds":
            newline = ' [ bonds ]'
            mode = "bonds"
            group = None
        elif line_arr[0] == "imphd":
            newline = ' [ impropers ]'
            mode = "imphd"
            group = None
        elif line_arr[0] == "termatom":
            pass
        elif line_arr[0] == "RESIDUE_END":
            pass
        else:
            raise RuntimeError("Bad line with no active mode: "+line)
    elif mode == "atoms":
        if line == "group":
            group += 1 
        else:
            newline = chargeline(line, group, precision) 
            if newline == None:
                print "BAD PROGRAM, SAD PROGRAM"
                print line
    elif mode == "bonds":
        newline = line.upper()
    elif mode == "imphd":
        newline = line.upper()
    else:
        raise RuntimeError("No acceptable active mode or line data.\n\tmode = "+str(mode)+"\n\tline = "+line)

    return newline, mode, group


class FileType:
    PRM = 1
    TPG = 2

default_precision = 10

def ManageFFFiles(filename, filetype, fields):
    file_out = ""
    with open(filename, 'r') as f:
        mode = None
        group = None
        for l in f:
            try:
                if   filetype == FileType.PRM:
                    newline, mode = manageprmlines(l,mode, precision = default_precision)
                elif filetype == FileType.TPG:
                    newline, mode, group = managetpglines(l, mode, group, precision = default_precision)
                if newline != None and (mode in fields):
                    file_out += newline + "\n"
            except RuntimeError as error_info:
                print "Aborting on RuntimeError encounter"
                print error_info
                raise RuntimeError(file_out)
            except AssertionError:
                _,_,tb = sys.exc_info()
                traceback.print_tb(tb) # Fixed format

                tbInfo = traceback.extract_tb(tb)
                filename,line,func,text = tbInfo[-1]
                print ('An error occurred on line ' + str(line) + ' in statement ' + text)
                print ("Bad handling for line = "+l)
                raise RuntimeError(file_out)
    return file_out

prm_all_fields=["BOND","BENDINGS","TORSION PROPER","TORSION IMPROPER","NONBONDED MIXRULE"]
tpg_all_fields = ["RESIDUE","atoms","bonds", "imphd"]



def generate_ffbonded(prm_filename):
    fields = set(prm_all_fields[0:-1])
    ffbonded = ManageFFFiles(prm_filename, FileType.PRM, fields)
    return ffbonded

def generate_ffnonbonded(prm_filename):
    fields = set(["NONBONDED MIXRULE"])
    ffnonbonded = ManageFFFiles(prm_filename, FileType.PRM, fields)
    return ffnonbonded

def generate_rtp(tpg_filename):
    fields = set(tpg_all_fields)
    rtp = ManageFFFiles(tpg_filename, FileType.TPG, fields)
    return rtp_header + rtp

def generate_atm(tpg_filename):
    fields = set(["atoms"])
    rtp = ManageFFFiles(tpg_filename, FileType.TPG, fields)
    lines = [ x.split() for x in rtp.split("\n")]
    for x in lines:
        if len(x) != 4:
            lines.remove(x)
    lines = np.array(lines)
    lines = np.array(list(set(lines[:,1])))
    lines.shape = (len(lines),)
    mass_arr = np.array([atom2mass.get(x[0]) for x in lines])
    return_arr = np.column_stack((lines, mass_arr))
    return_str = [set_pos(str(x[0]), 6) + str(x[1]) for x in return_arr]
    return_str = "\n".join(return_str)
    return return_str


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-AMBER94prm',   required=True,  
            help='Forcefield paramers (.prm) file for AMBER94 forcefield')
    parser.add_argument('-AMBER94tpg',   required=True,  
            help='Forcefield topology (.tpg) file for AMBER94 forcefield')
    parser.add_argument('-GMXbonded',     required=True,
            help='Forcefield output for GROMACS')
    parser.add_argument('-GMXnonbonded',  required=True,
            help='Forcefield output for GROMACS')
    parser.add_argument('-GMXrtp',       required=True,
            help='Forcefield output for GROMACS')
    parser.add_argument('-GMXatomtypes',  required=True,
            help='Forcefield output for GROMACS')
    parser.add_argument('-suffix', default='',
            help='Suffix to append to atomtypes')
    parser.add_argument('-bond_scale', type=float,
            help='Scale for all bond energies')
    parser.add_argument('-angle_scale', type=float,
            help='Scale for all angle energies')
    parser.add_argument('-dihedral_scale', type=float,
            help='Scale for all dihedral energies')
    parser.add_argument('-VDW_scale', type = float,
            help='Scale for VDW length parameters')
    args = parser.parse_args()

    atomtypeSuffix = args.suffix
    if not args.bond_scale is None:
        bond_scale = args.bond_scale
    if not args.angle_scale is None:
        angle_scale = args.angle_scale
    if not args.dihedral_scale is None:
        dihedral_scale = args.dihedral_scale
    if not args.VDW_scale is None:
        VDW_scale  = args.VDW_scale

    if args.GMXbonded:
        ffbonded = generate_ffbonded(args.AMBER94prm)
        with open(args.GMXbonded,'w') as f:
            f.write(ffbonded)
    if args.GMXnonbonded:
        ffnonbonded = generate_ffnonbonded(args.AMBER94prm)
        with open(args.GMXnonbonded,'w') as f:
            f.write(ffnonbonded)
    rtp = generate_rtp(args.AMBER94tpg)
    with open(args.GMXrtp,'w') as f:
        f.write(rtp)
    atm = generate_atm(args.AMBER94tpg)
    with open(args.GMXatomtypes,'w') as f:
        f.write(atm)

