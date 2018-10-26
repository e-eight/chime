#! /usr/bin/env python

from itertools import product
from subprocess import run

class MatrixFile:
    """
    A class that generates the files containing the matrix elements.
    """

    def __init__(self, operator, osc_energy, nmax, jmax):
        self.operator = operator
        self.osc_energy = osc_energy
        self.nmax = nmax
        self.jmax = jmax

    @classmethod
    def from_input(cls):
        operator = input('Operator: ')
        osc_energy = input('Oscillator energy: ')
        nmax = input('Maximum quanta: ')
        jmax = input('Maximum angular momentum: ')
        return cls(operator, osc_energy, nmax, jmax)

    def generate_files(self):

        def create_file():
            sep = ""
            nval = sep.join(['nmax', self.nmax])
            jval = sep.join(['jmax', self.jmax])
            osc_val = sep.join(['hw', self.osc_energy])
            sep = "."
            fname = sep.join([self.operator, order, osc_val, nval, jval, 'dat'])
            datafile = open(fname, 'w')
            datafile.write('# ISU two body reduced matrix elements\n')
            datafile.write(f'# {self.operator.upper()}\n\n')
            datafile.close()
            return fname

        def add_to_file(fname):
            datafile = open(fname, 'a')
            datafile.write(f'# si, ji, ti, sf, jf, tf, hw(MeV)\n')
            datafile.write(f'# ni li nf lf ME\n')
            datafile.write(f'{si} {ji} {ti} {sf} {jf} {tf} {self.osc_energy}\n')
            datafile.close()
            datafile = open(fname, 'a')
            run(["./matelm.exe", self.operator, order, self.osc_energy, 
                 self.nmax, str(lmin_i), str(lmax_i), str(lmin_f), str(lmax_f),
                 str(si), str(sf), str(ji), str(jf), str(ti), str(tf)],
                stdout=datafile)
            datafile.write('\n')
            datafile.close()

        orders = ["bare", "lo", "nlo", "n2lo", "all"]

        for order in orders:
            filename = create_file()
            for (ji, jf) in product(range(int(self.jmax) + 1), repeat=2):
                for (si, sf, ti, tf) in product(range(2), repeat=4):
                    lmin_i, lmax_i = ji - si, ji + si
                    lmin_f, lmax_f = jf - sf, jf + sf
                    if (lmin_i >= 0 and lmin_f >= 0):
                        add_to_file(filename)
            
def main():
    run(["make", "clean"])
    run("make")

    mf = MatrixFile.from_input()
    mf.generate_files()

if __name__ == "__main__":
    main()
