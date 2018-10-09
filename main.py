from subprocess import run

def main():
    f = open("rsq_lo.dat", "w");
    f.write("# ISU Two Body Reduced Matrix Elements\n")
    f.write("# Operator = r^2\n")
    f.write("# Order = LO\n")
    f.write("# Relevant LECs: None\n")
    f.write("# Oscillator Energy: 20 MeV\n")
    f.write("# s = 1, j = 1, t = 0, s' = 1, j' = 1, t' = 0\n")
    f.close()
    f = open("rsq_lo.dat", "a")
    run("./matelm.exe", stdout=f)
    f.close()

if __name__ == "__main__":
    main()
