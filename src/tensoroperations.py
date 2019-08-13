# pylint: disable=invalid-name

"""
Module to gather all additional tensor operations not present in UFL.
This especially includes all 3D operations and operations on tensors with
rank above 2.

*References*

.. [STR2005] Struchtrup, H. (2005). Macroscopic transport equations for rarefied gas flows. In Macroscopic Transport Equations for Rarefied Gas Flows. Springer, Berlin, Heidelberg.

.. [TOR2018] Torrilhon, M. et al. “Kinetic Theory of Non-Equilibrium Gas Flows: Theory and Computations”. Lecture Notes. Indian Institute of Technology Madras, 2018.
"""

import dolfin as df
import ufl

def dev3d2(rank2_2d):
    r"""
    Returns the 3D-synthetic 2D deviator
    :math:`B \in \mathbb{R}^{2 \times 2}`
    of the 2D 2-tensor
    :math:`A \in \mathbb{R}^{2 \times 2}`.

    .. math::

        A = \begin{pmatrix}
                a_{xx} & a_{xy} \\
                a_{yx} & a_{yy}
            \end{pmatrix}

        B = (A)_\mathrm{dev} = \frac{1}{2} (A)_\mathrm{sym}
            - \frac{1}{3} \mathrm{tr}(A) I_{2 \times 2}
    """
    symm = rank2_2d + ufl.transpose(rank2_2d)
    return (
        1/2 * symm
        - (1/3) * ufl.tr(rank2_2d) * ufl.Identity(2)
    )

def sym3d3(rank3_3d):
    r"""
    Returns the symmetric part of a 3D rank-3 tensor [STR2005]_ [TOR2018]_.

    .. math::

        A_{(i j k)}=\frac{1}{6}\left(A_{i j k}+A_{i k j}+A_{j i k}+A_{j k i}+A_{k i j}+A_{k j i}\right)
    """
    i, j, k = ufl.indices(3)
    symm_ijk = 1/6 * (
        # all permutations
        + rank3_3d[i, j, k]
        + rank3_3d[i, k, j]
        + rank3_3d[j, i, k]
        + rank3_3d[j, k, i]
        + rank3_3d[k, i, j]
        + rank3_3d[k, j, i]
    )
    return ufl.as_tensor(symm_ijk, (i, j, k))

def dev3dOfSym3(rank3_3d):
    r"""
    Return the deviator of a symmetric rank-3 tensor [STR2005]_ [TOR2018]_.
    Henning p231ff
    Manuels lecture p51
    """
    i, j, k, l = ufl.indices(4)
    delta = df.Identity(3)

    sym_ijk = sym3d3(rank3_3d)[i, j, k]
    trace_ijk = 1/5 * (
        + sym3d3(rank3_3d)[i, l, l] * delta[j, k]
        + sym3d3(rank3_3d)[l, j, l] * delta[i, k]
        + sym3d3(rank3_3d)[l, l, k] * delta[i, j]
    )
    tracefree_ijk = sym_ijk - trace_ijk
    return ufl.as_tensor(tracefree_ijk, (i, j, k))

def gen3dTF2(rank2_2d):
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
        [rank2_2d[0, 0], rank2_2d[0, 1], 0],
        [rank2_2d[1, 0], rank2_2d[1, 1], 0],
        [0, 0, -rank2_2d[0, 0]-rank2_2d[1, 1]]
    ])

def gen3d2(rank2_2d):
    r"""
    Returns the synthetic 3D version
    :math:`A \in \mathbb{R}^{3 \times 3}`
    of a 2D rank-2 tensor
    :math:`B \in \mathbb{R}^{2 \times 2}`.
    The 3D-components are set to zero.

    .. math::

        B = \begin{pmatrix}
                b_{xx} & b_{xy} \\
                b_{yx} & b_{yy}
            \end{pmatrix}

        A = \begin{pmatrix}
                b_{xx} & b_{xy} & 0                \\
                b_{yx} & b_{yy} & 0                \\
                0      & 0      & 0
            \end{pmatrix}
    """
    return df.as_tensor([
        [rank2_2d[0, 0], rank2_2d[0, 1], 0],
        [rank2_2d[1, 0], rank2_2d[1, 1], 0],
        [0, 0, 0]
    ])

def grad3dOf2(rank2_3d):
    """
    Returns the 3D version gradient of a 3D synthetic tracefree tensor,
    created from a 2D rank-2 tensor.
    """
    grad2d = df.grad(rank2_3d)
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
