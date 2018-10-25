 #! /usr/bin/env python

from itertools import product
from subprocess import run

def add_to_file(operator, osc_energy, order, nmax, jmax,
                lmin_i, lmax_i, lmin_f, lmax_f,
                si, sf, ji, jf, ti, tf):
    sep = ""
    jval = sep.join(["j", jmax])
    osc_val = sep.join([osc_energy, "MeV"])
    sep = "."
    fname = sep.join([operator, order, osc_val, jval, 'dat'])
    out = open(fname, "a+")
    out.write(f"# s = {si}, j = {ji}, t = {ti}, ")
    out.write(f"s' = {sf}, j' = {jf}, t' = {tf}\n")
    run(["./matelm.exe", operator, osc_energy, order, nmax, lmin_i,
         lmax_i, lmin_f, lmax_f, si, sf, ji, jf, ti, tf], stdout=out)
    out.write("\n")
    out.close()    

def generate_all_channels(operator, order, osc_energy, nmax, jmax):
    for (ji, jf) in product(range(jmax + 1), range(jmax + 1)):
        for (si, sf, ti, tf) in product(range(2), range(2), range(2), range(2)):
            lmin_i, lmax_i = ji - si, ji + si
            lmin_f, lmax_f = jf - sf, jf + sf
            if (lmin_i >= 0 and lmin_f >= 0):
                add_to_file(operator, order, osc_energy, nmax, jmax
                            lmin_i, lmax_i, lmin_f, lmax_f,
                            si, sf, ji, jf, ti, tf)
                
def main():
    run(["make", "clean"])
    run("make")
    
    operator = input("Operator? (rsq)")
    osc_energy = input("Oscillator Energy? ")
    nmax = input("Enter the maximum radial quantum number. ")

    sep = ""
    is_coupled, li = None, None
    lmin_i, lmax_i = 0, 0

    while (is_coupled != "y" and is_coupled != "n"):
        is_coupled = input("Is the intial channel coupled? (y/n) ")
        if (is_coupled == "y"):
            lmin_i = int(input("Enter the mininum value of l. "))
            lmax_i = lmin_i + 2
            li = sep.join(["l", str(lmin_i), str(lmax_i)])
        elif (is_coupled == "n"):
            lmin_i = int(input("Enter the value of l. "))
            lmax_i = lmin_i
            li = sep.join(["l", str(lmin_i)])
        else:
            print("Please enter y or n.")

    is_coupled, lf = None, None
    lmin_f, lmax_f = 0, 0

    while (is_coupled != "y" and is_coupled != "n"):
        is_coupled = input("Is the final channel coupled? (y/n) ")
        if (is_coupled == "y"):
            lmin_f = int(input("Enter the mininum value of l'. "))
            lmax_f = lmin_f + 2
            lf = sep.join(["lp", str(lmin_f), str(lmax_f)])
        elif (is_coupled == "n"):
            lmin_f = int(input("Enter the value of l'. "))
            lmax_f = lmin_f
            lf = sep.join(["lp", str(lmin_f)])
        else:
            print("Please enter y or n.")    
    
    ji = input("Enter the value of j. ")
    jf = input("Enter the value of j'. ")

    fn_ji = sep.join(["j", ji])
    fn_jf = sep.join(["jp", jf])
    osc_energy = sep.join([osc_energy, "MeV"])

    sep = "."
    cname = sep.join([operator, li, fn_ji, lf, fn_jf])

    filename = sep.join([cname, "bare", osc_energy, "dat"])
    f = open(filename, "w");
    f.write("# ISU Two Body Reduced Matrix Elements\n")
    f.write("# Operator = ")
    f.write(operator)
    f.write("\n")
    f.write("# Order = LO\n")
    f.write("# Relevant LECs: None\n")
    f.write("# Oscillator Energy: 20 MeV\n")
    f.write(f"# s = 1, j = {ji}, t = 0, s' = 1, j' = {jf}, t' = 0\n")
    f.close()
    f = open(filename, "a")
    run(["./matelm.exe", operator, "bare", osc_energy, nmax, str(lmin_i),
         str(lmax_i), str(lmin_f), str(lmax_f), ji, jf], stdout=f)
    f.close()

    filename = sep.join([cname, "lo", osc_energy, "dat"])
    f = open(filename, "w");
    f.write("# ISU Two Body Reduced Matrix Elements\n")
    f.write("# Operator = ")
    f.write(operator)
    f.write("\n")
    f.write("# Order = LO\n")
    f.write("# Relevant LECs: None\n")
    f.write("# Oscillator Energy: 20 MeV\n")
    f.write(f"# s = 1, j = {ji}, t = 0, s' = 1, j' = {jf}, t' = 0\n")
    f.close()
    f = open(filename, "a")
    run(["./matelm.exe", operator, "lo", osc_energy, nmax, str(lmin_i),
         str(lmax_i), str(lmin_f), str(lmax_f), ji, jf], stdout=f)
    f.close()

    filename = sep.join([cname, "n2lo", osc_energy, "dat"])
    f = open(filename, "w");
    f.write("# ISU Two Body Reduced Matrix Elements\n")
    f.write("# Operator = ")
    f.write(operator)
    f.write("\n")
    f.write("# Order = N2LO\n")
    f.write("# Relevant LECs: None\n")
    f.write("# Oscillator Energy: 20 MeV\n")
    f.write(f"# s = 1, j = {ji}, t = 0, s' = 1, j' = {jf}, t' = 0\n")
    f.close()
    f = open(filename, "a")
    run(["./matelm.exe", operator, "n2lo", osc_energy, nmax, str(lmin_i),
         str(lmax_i), str(lmin_f), str(lmax_f), ji, jf], stdout=f)
    f.close()

    filename = sep.join([cname, "all", osc_energy, "dat"])
    f = open(filename, "w");
    f.write("# ISU Two Body Reduced Matrix Elements\n")
    f.write("# Operator = ")
    f.write(operator)
    f.write("\n")
    f.write("# Order = LO\n")
    f.write("# Relevant LECs: None\n")
    f.write("# Oscillator Energy: 20 MeV\n")
    f.write(f"# s = 1, j = {ji}, t = 0, s' = 1, j' = {jf}, t' = 0\n")
    f.close()
    f = open(filename, "a")
    run(["./matelm.exe", operator, "all", osc_energy, nmax, str(lmin_i),
         str(lmax_i), str(lmin_f), str(lmax_f), ji, jf], stdout=f)
    f.close()

if __name__ == "__main__":
    main()
