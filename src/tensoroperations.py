# pylint: disable=invalid-name

"""
Module to gather all additional tensor operations not present in UFL.
This especially includes all 3D operations and operations on tensors with
rank above 2.
"""

import dolfin as df
import ufl

def dev3d(mat):
    "2d deviatoric part of actually 3d matrix"
    return (
        0.5 * (mat + ufl.transpose(mat))
        - (1/3) * ufl.tr(mat) * ufl.Identity(2)
    )

def sym3(rank3):
    """
    Returns the symmetric part of a rank-3 tensor
    Henning p231ff
    """
    i, j, k = ufl.indices(3)
    symm_ijk = 1/6 * (
        # all permutations
        + rank3[i, j, k]
        + rank3[i, k, j]
        + rank3[j, i, k]
        + rank3[j, k, i]
        + rank3[k, i, j]
        + rank3[k, j, i]
    )
    return ufl.as_tensor(symm_ijk, (i, j, k))

def dev3(rank3):
    """
    Return the deviator of a rank-3 tensor
    Henning p231ff
    Manuels lecture p51
    """
    i, j, k, l = ufl.indices(4)
    delta = df.Identity(3)

    sym_ijk = sym3(rank3)[i, j, k]
    trace_ijk = 1/5 * (
        + sym3(rank3)[i, l, l] * delta[j, k]
        + sym3(rank3)[l, j, l] * delta[i, k]
        + sym3(rank3)[l, l, k] * delta[i, j]
    )
    tracefree_ijk = sym_ijk - trace_ijk
    return ufl.as_tensor(tracefree_ijk, (i, j, k))

def gen3dTF2(rank2):
    r"""
    Returns the synthetic 3D version
    :math:`A \in \mathbb{R}^{3 \times 3}`
    of a 2D rank-2 tensor
    :math:`B \in \mathbb{R}^{2 \times 2}`.

    .. math::

        B = \begin{pmatrix}
                b_{xx} & b_{xy} \\
                b_{yx} & b_{yy}
            \end{pmatrix}

        A = \begin{pmatrix}
                b_{xx} & b_{xy} & 0                \\
                b_{yx} & b_{yy} & 0                \\
                0      & 0      & -(b_{yx}+b_{yy})
            \end{pmatrix}
    """
    return df.as_tensor([
        [rank2[0, 0], rank2[0, 1], 0],
        [rank2[1, 0], rank2[1, 1], 0],
        [0, 0, -rank2[0, 0]-rank2[1, 1]]
    ])

def grad3dOf2(rank2):
    """
    Returns the 3D version gradient of a 3D synthetic tracefree tensor,
    created from a 2D rank-2 tensor.
    """
    grad2d = df.grad(rank2)
    dim3 = df.as_tensor([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
    ])
    grad3d = df.as_tensor([
        grad2d[:, :, 0], # pylint: disable=unsubscriptable-object
        grad2d[:, :, 1], # pylint: disable=unsubscriptable-object
        dim3[:, :]
    ])
    return grad3d
