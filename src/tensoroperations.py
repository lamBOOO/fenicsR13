# pylint: disable=invalid-name

"""
Module to gather all additional tensor operations not present in UFL.

This especially includes all 3D operations and operations on tensors with rank
above 2.

.. [STR2005] Struchtrup, H. (2005). Macroscopic transport equations for rarefied
    gas flows. In Macroscopic Transport Equations for Rarefied Gas Flows.
    Springer, Berlin, Heidelberg.

.. [TOR2018] Torrilhon, M. et al. (2018). “Kinetic Theory of Non-Equilibrium
    Gas Flows:
    Theory and Computations”. Lecture Notes. Indian Institute of Technology
    Madras.
"""

import dolfin as df
import ufl

def dev3d2(rank2_2d):
    r"""
    Return the synthetic 3D deviator of a 2D 2-tensor.

    Returns the 3D-synthetic 2D deviator
    :math:`B \in \mathbb{R}^{2 \times 2}`
    of the 2D 2-tensor
    :math:`A \in \mathbb{R}^{2 \times 2}`.

    .. math::

        A &= \begin{pmatrix}
                a_{xx} & a_{xy} \\
                a_{yx} & a_{yy}
            \end{pmatrix}
        \\
        B &= (A)_\mathrm{dev} = \frac{1}{2} (A)_\mathrm{sym}
            - \frac{1}{3} \mathrm{tr}(A) I_{2 \times 2}
    """
    symm = rank2_2d + ufl.transpose(rank2_2d)
    return (
        1/2 * symm
        - (1/3) * ufl.tr(rank2_2d) * ufl.Identity(2)
    )

def sym3d3(rank3_3d):
    r"""
    Return the symmetric part of a 3D 3-tensor.

    Returns the symmetric part :math:`A_{(i j k)}` of a 3D rank-3 tensor
    :math:`A \in \mathbb{R}^{3 \times 3 \times 3}` [STR2005]_ [TOR2018]_.

    .. math::

        A_{(i j k)}=\frac{1}{6}\left(
            A_{i j k}+A_{i k j}+A_{j i k}+A_{j k i}+A_{k i j}+A_{k j i}
        \right)
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

def dev3d3(rank3_3d):
    r"""
    Return the deviator of a 3D 3-tensor.

    Returns the deviator :math:`A_{\langle i j k\rangle}` of a rank-3 tensor
    :math:`A \in \mathbb{R}^{3 \times 3 \times 3}` [TOR2018]_.

    .. math::

        A_{\langle i j k\rangle}=A_{(i j k)}-\frac{1}{5}\left[A_{(i l)}
        \delta_{j k}+A_{(l j l)} \delta_{i k}+A_{(l l k)} \delta_{i j}\right]

    A gradient :math:`\frac{\partial S_{\langle i j}}{\partial x_{k \rangle}}`
    of a symmetric 2-tensor :math:`S \in \mathbb{R}^{3 \times 3}` has the
    deviator [STR2005]_:

    .. math::

        \frac{\partial S_{\langle i j}}{\partial x_{k \rangle}}
        =&
        \frac{1}{3}\left(
             \frac{\partial S_{i j}}{\partial x_{k}}
            +\frac{\partial S_{i k}}{\partial x_{j}}
            +\frac{\partial S_{j k}}{\partial x_{i}}
        \right)
        -\frac{1}{15}
        \biggl[
            \left(
                 \frac{\partial S_{i r}}{\partial x_{r}}
                +\frac{\partial S_{i r}}{\partial x_{r}}
                +\frac{\partial S_{r r}}{\partial x_{i}}
            \right) \delta_{j k}
        \\
        &+
        \left(
             \frac{\partial S_{j r}}{\partial x_{r}}
            +\frac{\partial S_{j r}}{\partial x_{r}}
            +\frac{\partial S_{r r}}{\partial x_{j}}
        \right) \delta_{i k}
        +\left(
             \frac{\partial S_{k r}}{\partial x_{r}}
            +\frac{\partial S_{k r}}{\partial x_{r}}
            +\frac{\partial S_{r r}}{\partial x_{k}}
        \right) \delta_{i j}
        \biggr]
    """
    i, j, k, l = ufl.indices(4)
    delta = df.Identity(3)

    sym_ijk = sym3d3(rank3_3d)[i, j, k]
    traces_ijk = 1/5 * (
        + sym3d3(rank3_3d)[i, l, l] * delta[j, k]
        + sym3d3(rank3_3d)[l, j, l] * delta[i, k]
        + sym3d3(rank3_3d)[l, l, k] * delta[i, j]
    )
    tracefree_ijk = sym_ijk - traces_ijk
    return ufl.as_tensor(tracefree_ijk, (i, j, k))

def gen3dTF2(rank2_2d):
    r"""
    Generate a 3D tracefree 2-tensor from a 2D 2-tensor.

    Returns the synthetic 3D version
    :math:`A \in \mathbb{R}^{3 \times 3}`
    of a 2D rank-2 tensor
    :math:`B \in \mathbb{R}^{2 \times 2}`.

    .. math::

        B &= \begin{pmatrix}
                b_{xx} & b_{xy} \\
                b_{yx} & b_{yy}
            \end{pmatrix}
        \\
        A &= \begin{pmatrix}
                b_{xx} & b_{xy} & 0 \\
                b_{yx} & b_{yy} & 0 \\
                0      & 0      & -(b_{yx}+b_{yy})
            \end{pmatrix}

    Example
    -------

    >>> t2d = df.Expression((("1","2"),("3","4")), degree=1)
    >>> t3d = gen3dTF2(t2d)
    >>> print(
    ...     t3d[0,0]==t2d[0,0] and
    ...     t3d[0,1]==t2d[0,1] and
    ...     t3d[0,2]==0 and
    ...     t3d[1,0]==t2d[1,0] and
    ...     t3d[1,1]==t2d[1,1] and
    ...     t3d[1,2]==0 and
    ...     t3d[2,0]==0 and
    ...     t3d[2,1]==0 and
    ...     t3d[2,2]==-t2d[0,0]-t2d[1,1]
    ... )
    True

    """
    return df.as_tensor([
        [rank2_2d[0, 0], rank2_2d[0, 1], 0],
        [rank2_2d[1, 0], rank2_2d[1, 1], 0],
        [0, 0, -rank2_2d[0, 0]-rank2_2d[1, 1]]
    ])

def gen3d2(rank2_2d):
    r"""
    Generate a 3D 2-tensor from a 2D 2-tensor (add zeros to last dimensions).

    Returns the synthetic 3D version
    :math:`A \in \mathbb{R}^{3 \times 3}`
    of a 2D rank-2 tensor
    :math:`B \in \mathbb{R}^{2 \times 2}`.
    The 3D-components are set to zero.

    .. math::

        B &= \begin{pmatrix}
                b_{xx} & b_{xy} \\
                b_{yx} & b_{yy}
            \end{pmatrix}
        \\
        A &= \begin{pmatrix}
                b_{xx} & b_{xy} & 0 \\
                b_{yx} & b_{yy} & 0 \\
                0      & 0      & 0
            \end{pmatrix}

    Example
    -------

    >>> t2d = df.Expression((("1","2"),("3","4")), degree=1)
    >>> t3d = gen3d2(t2d)
    >>> print(
    ...     t3d[0,0]==t2d[0,0] and
    ...     t3d[0,1]==t2d[0,1] and
    ...     t3d[0,2]==0 and
    ...     t3d[1,0]==t2d[1,0] and
    ...     t3d[1,1]==t2d[1,1] and
    ...     t3d[1,2]==0 and
    ...     t3d[2,0]==0 and
    ...     t3d[2,1]==0 and
    ...     t3d[2,2]==0
    ... )
    True

    """
    return df.as_tensor([
        [rank2_2d[0, 0], rank2_2d[0, 1], 0],
        [rank2_2d[1, 0], rank2_2d[1, 1], 0],
        [0, 0, 0]
    ])

def grad3dOf2(rank2_3d):
    r"""
    Return 3D gradient of 3D 2-tensor.

    Returns the 3D gradient (w.r.t. :math:`x,y,z`)
    :math:`\frac{\partial A_{ij}}{\partial x_{k}}`
    of a 3D :math:`z`-independent rank-2 tensor
    :math:`A(x,y) \in \mathbb{R}^{3 \times 3}`
    where

    .. math::

        A_{ij} = \begin{pmatrix}
                    a_{xx}(x,y) & a_{xy}(x,y) & a_{xz}(x,y) \\
                    a_{yx}(x,y) & a_{yy}(x,y) & a_{yz}(x,y) \\
                    a_{zx}(x,y) & a_{zy}(x,y) & a_{zz}(x,y)
                \end{pmatrix}

    The function then returns the 3-tensor with components
    :math:`\frac{\partial A_{ij}}{\partial x_{k}}`
    where

    .. math::

        \frac{\partial A_{ij}}{\partial x_{1}}
        =
        \frac{\partial A_{ij}}{\partial x}
        &= \begin{pmatrix}
                \partial_x a_{xx}(x,y) &
                \partial_x a_{xy}(x,y) &
                \partial_x a_{xz}(x,y) \\
                \partial_x a_{yx}(x,y) &
                \partial_x a_{yy}(x,y) &
                \partial_x a_{yz}(x,y) \\
                \partial_x a_{zx}(x,y) &
                \partial_x a_{zy}(x,y) &
                \partial_x a_{zz}(x,y)
            \end{pmatrix}
        , \\
        \frac{\partial A_{ij}}{\partial x_{2}}
        =
        \frac{\partial A_{ij}}{\partial y}
        &= \begin{pmatrix}
                \partial_y a_{xx}(x,y) &
                \partial_y a_{xy}(x,y) &
                \partial_y a_{xz}(x,y) \\
                \partial_y a_{yx}(x,y) &
                \partial_y a_{yy}(x,y) &
                \partial_y a_{yz}(x,y) \\
                \partial_y a_{zx}(x,y) &
                \partial_y a_{zy}(x,y) &
                \partial_y a_{zz}(x,y)
            \end{pmatrix}
        , \\
        \frac{\partial A_{ij}}{\partial x_{3}}
        =
        0
        &= \begin{pmatrix}
                0 &
                0 &
                0 \\
                0 &
                0 &
                0 \\
                0 &
                0 &
                0
            \end{pmatrix}
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
