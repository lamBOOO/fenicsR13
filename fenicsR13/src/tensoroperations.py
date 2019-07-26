# pylint: disable=invalid-name

"""
Module to gather all additional tensor operations not present in UFL.
"""

import dolfin as df
import ufl


def dev3d(mat):
    "2d deviatoric part of actually 3d matrix"
    return (
        0.5 * (mat + ufl.transpose(mat))
        - (1/3) * ufl.tr(mat) * ufl.Identity(2)
    )

def innerOfTracefree2(rank2_1, rank2_2):
    """
    Return the 3D inner prodcut of two symmetric tracefree
    rank-2 tensors (two-dimensional) as used the weak form of
    westerkamp2019.
    """
    return (
        df.inner(rank2_1, rank2_2)
        # part from u_zz=-(u_xx+u_yy) contribution
        + rank2_1[0][0] * rank2_2[0][0]
        + rank2_1[0][0] * rank2_2[1][1]
        + rank2_1[1][1] * rank2_2[0][0]
        + rank2_1[1][1] * rank2_2[1][1]
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
    """
    i, j, k, l = ufl.indices(4)
    delta = df.Identity(3)

    sym_ijk = sym3(rank3)[i, j, k]
    trace_ijk = 1/5 * (
        + sym3(rank3)[i, l, l] * delta[j, k]
        + sym3(rank3)[j, l, l] * delta[i, k]
        + sym3(rank3)[k, l, l] * delta[i, j]
    )
    tracefree_ijk = sym_ijk - trace_ijk
    return ufl.as_tensor(tracefree_ijk, (i, j, k))

def innerOfDevOfGrad2AndGrad2(sigma, psi):
    r"""
    Implements the inner product of the deviator of a symmetric,
    tracefree rank-2 tensor with a symmetric, tracefree rank-2 tensor.

    .. math::

        (\nabla \underline{\underline{\sigma}})_{\mathrm{dev}} :
            \nabla \underline{\underline{\psi}}
    """
    hardcoded = False
    if hardcoded: # pylint: disable=no-else-return
        return psi[1, 1].dx(1)*((3*sigma[1, 1].dx(1))/5. - sigma[0, 1].dx(0)/5. - sigma[1, 0].dx(0)/5.) + (-psi[0, 0].dx(1) - psi[1, 1].dx(1))*(-sigma[0, 0].dx(1)/3. - (7*sigma[1, 1].dx(1))/15. - sigma[0, 1].dx(0)/15. - sigma[1, 0].dx(0)/15.) + psi[0, 0].dx(1)*(sigma[0, 0].dx(1)/3. - (2*sigma[1, 1].dx(1))/15. + (4*sigma[0, 1].dx(0))/15. + (4*sigma[1, 0].dx(0))/15.) + psi[0, 1].dx(1)*((4*sigma[0, 1].dx(1))/15. + (4*sigma[1, 0].dx(1))/15. - (2*sigma[0, 0].dx(0))/15. + sigma[1, 1].dx(0)/3.) + psi[1, 0].dx(1)*((4*sigma[0, 1].dx(1))/15. + (4*sigma[1, 0].dx(1))/15. - (2*sigma[0, 0].dx(0))/15. + sigma[1, 1].dx(0)/3.) + (-sigma[0, 1].dx(1)/5. - sigma[1, 0].dx(1)/5. + (3*sigma[0, 0].dx(0))/5.)*psi[0, 0].dx(0) + (sigma[0, 0].dx(1)/3. - (2*sigma[1, 1].dx(1))/15. + (4*sigma[0, 1].dx(0))/15. + (4*sigma[1, 0].dx(0))/15.)*psi[0, 1].dx(0) + (sigma[0, 0].dx(1)/3. - (2*sigma[1, 1].dx(1))/15. + (4*sigma[0, 1].dx(0))/15. + (4*sigma[1, 0].dx(0))/15.)*psi[1, 0].dx(0) + (-sigma[0, 1].dx(1)/15. - sigma[1, 0].dx(1)/15. - (7*sigma[0, 0].dx(0))/15. - sigma[1, 1].dx(0)/3.)*(-psi[0, 0].dx(0) - psi[1, 1].dx(0)) + ((4*sigma[0, 1].dx(1))/15. + (4*sigma[1, 0].dx(1))/15. - (2*sigma[0, 0].dx(0))/15. + sigma[1, 1].dx(0)/3.)*psi[1, 1].dx(0)
    else:
        return df.inner(
            dev3(grad3dOf2(gen3dTracefreeTensor(sigma))),
            grad3dOf2(gen3dTracefreeTensor(psi))
            # dev3(grad3dOf2(gen3dTracefreeTensor(psi))) # same
        )

def gen3dTracefreeTensor(rank2):
    r"""
    Returns the synthetic 3D version
    :math:`A \in \mathbb{R}^{3 \times 3}`
    of a 2D rank-2 tensor
    :math:`B \in \mathbb{R}^{2 \times 2}`.

    .. math::

        B = \begin{pmatrix}
                b_{xx} & b_{xy} \
                b_{yx} & b_{yy} \
            \end{pmatrix}

        A = \begin{pmatrix}
                b_{xx} & b_{xy} & 0                \
                b_{yx} & b_{yy} & 0                \
                0      & 0      & -(b_{yx}+b_{yy}) \
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
