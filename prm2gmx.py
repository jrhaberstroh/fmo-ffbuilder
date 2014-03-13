
twoletter = { 'cab' : 'CD',\
	      'cbb' : 'CE',\
	      'cnb' : 'CF',\
	      'cpb' : 'CG',\
	      'cqb' : 'CH',\
	      'crb' : 'CI',\
	      'csb' : 'CJ',\
	      'ccs' : 'CL',\
	      'cqq' : 'CO',\
	      'c2k' : 'CP',\
	      'c2a' : 'C5',\
              'c2e' : 'CU',\
	      'cq2' : 'CX',\
	      'ct1' : 'C1',\
	      'ct2' : 'C2',\
	      'ct3' : 'C3',\
	      'nmh' : 'NM',\
	      'mgc' : 'MC',\
	      'ha0' : 'HQ',\
	      'o2c' : 'OC',\
	      'o1c' : 'O1',\
	      'rb'  : 'Rb',\
	      'br'  : 'Br',\
	      'cl'  : 'Cl',\
	      'na'  : 'Na',\
	      'cs'  : 'Cs',\
	      'li'  : 'Li',\
	      'zn'  : 'Zn'};

# Silly Notes:
#     'c2a' used to convert to CS, but this was changed to C5 for fear of overlap with Cs

mass2num = { "1.008"  : 1, \
	     "6.940"  : 3, \
	     "12.010" : 6, \
	     "12.011" : 6, \
	     "14.010" : 7, \
	     "15.999" : 8, \
	     "16.000" : 8, \
	     "19.000" : 9, \
	     "22.990" : 11,\
	     "24.305" : 12,\
	     "30.970" : 15,\
	     '32.060' : 16,\
	     '35.450' : 17,\
	     "39.100" : 19,\
	     "55.847" : 26,\
	     "85.470" : 37,\
	     "126.900": 53,\
	     "131.000": 54,\
	     "132.910": 55 };


kCal2kJ = 4.184
nm2A    = 10.

def generic_line(line_arr, entries_per_line, order, units, datatype=None, precision=6):
    ''' 
    generic_line()

    Converts line_arr into the correct format.

    order determines the order of the quantitative arguments
    units determines their units
    dataype determines the conversion for output (default float)
    precision determines how many chars to print in quantitative args

    Returns a string of line_out.
    '''
    if datatype==None:
        datatype=[float]*len(units)

    assert (len(line_arr) == (entries_per_line + len(order)))
    assert (len(units)    == len(order))
    assert (len(datatype) == len(order))

    line_out = ""

    for i in xrange(entries_per_line) :
        if (twoletter.get(line_arr[i])==None and len(line_arr[i]) >= 3):
            return None
        line_out += (line_arr[i] if twoletter.get(line_arr[i])==None else twoletter[line_arr[i]]) + " "

    for i in xrange(len(order)):
        offset = entries_per_line
        val = datatype[i](float(line_arr[offset + order[i]]) * units[i])
        # Check that there are no more than 5 chars to display
        if not ((val < float("1"+"0"*(precision-1)) and val >= float("."+"0"*(precision-2)+"1"))\
                or val == 0):
            raise RuntimeError("Quantity too large or too small (but non-zero) appeared (has more than 5 chars): "+str(val))
        val = str((val))
        val = val[0:min(precision-1, len(val))] 
        if i != len(order) -1:
            val += " "
        line_out += val

    return line_out


def bondline(line_in, precision = 6):
    return generic_line(line_in.split(), 2, [1,0], [1, kCal2kJ * nm2A**2], precision = precision);

def bendline(line_in, precision = 6):
    return generic_line(line_in.split(), 3, [1,0], [1.,kCal2kJ], precision = precision);

def torsionline(line_in, precision = 6):
    order = [2,0,1]
    units = [1., kCal2kJ, 1]
    datatype = [float, float, int]
    entries_per_line = 4;
    force_pos = entries_per_line
    line_arr = line_in.split()

    # Append the default angle value if the angle is missing
    # Also, convert the force to a positive number.
    if float(line_arr[force_pos]) < 0:
        if len(line_arr) == (entries_per_line + 3):
            return "FAIL"
    if len(line_arr) == (entries_per_line + 3) - 1:
       line_arr.append(0. if (float(line_arr[force_pos]) < 0) else 180.);
    line_arr[force_pos] = str(abs(float(line_arr[force_pos])))
    assert (len(line_arr) == (entries_per_line + 3))

    return generic_line(line_arr, 4, order, units, datatype, precision = precision)

def improperline(line_in, precision = 6):
    return torsionline(" ".join(line_in.split()[0:-1]), precision = precision)


def managelines(line, mode):
    if (line == "END\n"):
        mode = ""
        return None
    elif (mode == ""):
        mode = line
    elif (mode == "BOND"):
        return bondline(line)
    elif (mode == "BEND"):
        return bendline(line)
    elif (mode == "TORSION PROPER"):
        return torsionline(line)
    elif (mode == "TORSION IMPROPER"):
        return improperline(line)
    else:
        raise RuntimeError("No acceptable mode active. mode = "+str(mode))


def prm2gmx_mode(filename, bondtype):
    with open(filename, 'r') as f:
        mode = ""
        for l in f:
            managelines(l,mode)




#if __name__ == "__main__":


