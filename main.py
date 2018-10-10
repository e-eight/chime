#! /usr/bin/env python

from subprocess import run

def main():
    run(["make", "clean"])
    run("make")
    
    operator = input("Operator? (rsq, e2, m1) ")
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
