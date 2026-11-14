from scipy.sparse import csr_matrix
from numpy.linalg import eigvals
import matplotlib.pyplot as plt


from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
import matplotlib.pyplot as plt


def reorder(K, freedofs=None):
    """ Reorder a ngsolve sparse matrix using the reverse Cuthill-McKee algorithm """
    K = to_scipy(K, freedofs)
    perm = reverse_cuthill_mckee(K)
    K = K[perm, :][:, perm]
    return K

def to_scipy(K, freedofs=None):
    """ Convert a ngsolve sparse matrix to a scipy sparse matrix """
    rows, cols, vals = K.COO()
    K = csr_matrix((vals, (rows, cols)))
    if freedofs is not None:
        K = K[freedofs, :]
        K = K[:, freedofs]
    return K

def spectrum(bilinearform, freedofs=None):
    """ Get the eigenvalues of a ngsolve sparse matrix """
    K = bilinearform.Assemble().mat
    freedofs = bilinearform.space.FreeDofs()
    K = to_scipy(K, freedofs)
    return eigvals(K.toarray())


def spy(K, freedofs=None, *args, **kwargs):
    """ Visualize the sparsity pattern of a ngsolve sparse matrix """
    K = to_scipy(K, freedofs)
    plt.spy(K, *args, **kwargs)