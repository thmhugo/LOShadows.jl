import perceval as pcvl
import numpy as np
from perceval.utils import DensityMatrix
from perceval import Simulator
from perceval.backends import SLOSBackend
from perceval.algorithm import Sampler
from multiprocessing import Pool

import scipy

import tqdm as tq

np.set_printoptions(precision=3, suppress=True)
np.set_printoptions(threshold=np.inf)
from perceval.utils import StateVector


mzi = (
    pcvl.BS()
    // (0, pcvl.PS(phi=pcvl.Parameter("φ_a")))
    // pcvl.BS()
    // (1, pcvl.PS(phi=pcvl.Parameter("φ_b")))
)



def generate_all_fock_states(m, n):
    if n == 0:
        return [[0] * m]
    if m == 1:
        return [[n]]

    fock_states = []
    for i in range(n + 1):
        fock_states += [[i] + state for state in generate_all_fock_states(m - 1, n - i)]
    return fock_states




def get_state_density_matrix(circuit, input_state):
    try:
        input_state = list_repr_to_statevector(input_state)
    except:
        pass
    simulator = Simulator(SLOSBackend())
    simulator.set_circuit(circuit)
    input_dm = DensityMatrix.from_svd(input_state)
    return simulator.evolve_density_matrix(input_dm)


def partial_size(m, n, k):
    return sum([int(scipy.special.binom(m + i - 1, i)) for i in range(k + 1)])



def list_repr_to_statevector(lst):
    """necessary because statevector are not pickleable"""
    return sum([c * StateVector(str) for (c, str) in lst], StateVector())


def list_repr_to_basicstate(lst):
    """necessary because statevector are not pickleable"""
    return pcvl.BasicState(lst[0][1])


def direct_sum(A, B):
    A = np.array(A)
    B = np.array(B)

    rows_A, cols_A = A.shape
    rows_B, cols_B = B.shape

    result = np.zeros((rows_A + rows_B, cols_A + cols_B), dtype=A.dtype)

    result[:rows_A, :cols_A] = A

    result[rows_A:, cols_A:] = B

    return result


def direct_sum_arr(M):
    D = M[0]
    for i in range(1, len(M)):
        D = direct_sum(D, M[i])
    return D


def collect_samples(circuit, n, input_state, n_samples, simulator):
    """
    Input state given as a list representing a Fock basis state
    """
    # simulator = pcvl.Processor("SLOS")
    simulator.set_circuit(circuit)
    simulator.min_detected_photons_filter(n)
    simulator.with_input(list_repr_to_basicstate(input_state))
    sampler = Sampler(simulator, max_shots_per_call=n_samples)
    return sampler.sample_count(n_samples)


class Snapshot:
    def __init__(
        self,
        m,
        n,
        n_samples = 1,
        unknown_dm_mat=None,
        input_state=None,
        unknown_U=None,
        processor="SLOS",
        brightness=None,
        g2=None,
        indistinguishability=None,
    ):
        self.m = m
        self.n = n
        self.n_samples = 1
        self.unknown_dm_mat = unknown_dm_mat
        self.input_state = input_state
        self.unknown_U = unknown_U
        self.processor = processor
        self.brightness = (brightness,)
        self.g2 = (g2,)
        self.indistinguishability = (indistinguishability,)
        
    def collect_snapshots(self, K):
        """ """
        pool = Pool(10)

        snapshots = []

        with tq.tqdm(total=K, position=0, leave=False) as pbar:
            for s in tq.tqdm(
                pool.imap(self.get_snapshot, range(K)),
                total=K,
                position=0,
                leave=True,
            ):
                snapshots.append(s)

        return snapshots

    def get_snapshot(self, i):

        U = pcvl.Matrix.random_unitary(self.m)

        if self.unknown_dm_mat is not None:
            simulator = pcvl.Simulator(SLOSBackend())
            simulator.set_circuit(pcvl.Unitary(U))
            evolved_dm = simulator.evolve_density_matrix(
                DensityMatrix(
                    self.unknown_dm_mat, m=self.m, n_max=self.n, check_hermitian=False
                )
            )

            samples = evolved_dm.sample(self.n_samples)

        else:
            # if self.input_state is not None and self.unknown_U is not None:

            noise_model = pcvl.NoiseModel(
                brightness=self.brightness[0],
                g2=self.g2[0],
                indistinguishability=self.indistinguishability[0],
            )

            processor = (
                pcvl.Processor("SLOS", self.m, noise=noise_model)
                if self.processor == "SLOS"
                else pcvl.RemoteProcessor(self.processor)
            )

            samples = collect_samples(
                pcvl.Unitary(pcvl.Matrix(np.matmul(U, self.unknown_U))),
                self.n,
                self.input_state,
                self.n_samples,
                processor,
            )["results"]

        Udag = U.conj().T
        circuit_rand_udag = pcvl.Unitary(Udag)

        dms = [
            count * get_state_density_matrix(circuit_rand_udag, sample).mat.toarray()
            for (sample, count) in samples.items()
        ]

        rho_hat = sum(dms) * 1 / sum(samples.values())

        return rho_hat
