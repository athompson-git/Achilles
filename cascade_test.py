""" Run nuchic over a range of parameters. """

import os
from mpi4py import MPI
from string import Template
import ast
import numpy as np

from absl import flags, app

from nuchic.main import NuChic

FLAGS = flags.FLAGS

PARAMETERS = ['beam_energy']


def create_run(parameters, fname):
    """ Create a run card from a dictionary of parameters. """
    template = open('run_template.yml')
    src = Template(template.read())
    template.close()
    result = src.substitute(parameters)
    with open(fname, 'w') as run_card:
        for line in result:
            run_card.write(line)


def main(argv):
    """ Loop over all the parameters and generate results. """
    del argv

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    energies = None
    xsec_p = None
    xsec_n = None
    if rank == 0:
        if not os.path.exists('Results'):
            os.mkdir('Results')
        if not os.path.exists('runs'):
            os.mkdir('runs')

        emin = 250
        emax = 2500
        energies = np.linspace(emin, emax, (emax-emin)/10+1, dtype='float64')
        xsec_p = np.zeros_like(energies)
        xsec_n = np.zeros_like(energies)
        energies_split = np.array_split(energies, size)
        energies_size = [len(energies_split[i]) for i in range(len(energies_split))]
        energies_disp = np.insert(np.cumsum(energies_size), 0, 0)[0:-1]
    else:
        energies_split = None
        energies_size = None
        energies_disp = None

    energies_size = comm.bcast(energies_size, root = 0)
    energies_disp = comm.bcast(energies_disp, root = 0)
    energies_local = np.zeros(energies_size[rank])
    local_xsec_p = np.zeros(energies_size[rank])
    local_xsec_n = np.zeros(energies_size[rank])
    comm.Scatterv([energies, energies_size, energies_disp, MPI.DOUBLE], energies_local, root = 0)

    for i, energy in enumerate(energies_local):
        print(energy)
        create_run({PARAMETERS[0]: energy}, 'runs/run_{}.yml'.format(energy))
        main = NuChic('runs/run_{}.yml'.format(energy))
        xsec = main.calc_cross_section()
        local_xsec_p[i] = xsec[0]
        local_xsec_n[i] = xsec[1]

    comm.Barrier()
    comm.Gatherv(local_xsec_p, [xsec_p, energies_size, energies_disp, MPI.DOUBLE], root = 0)
    comm.Gatherv(local_xsec_n, [xsec_n, energies_size, energies_disp, MPI.DOUBLE], root = 0)

    if rank == 0:
        with open(os.path.join('Results', 'xsec.txt'), 'w') as output:
            for i, energy in enumerate(energies):
                output.write('{} {} {}\n'.format(energy, xsec_p[i], xsec_n[i]))


if __name__ == '__main__':
    app.run(main)
