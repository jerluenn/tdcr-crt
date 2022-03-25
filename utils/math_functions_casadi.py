from casadi import *
import numpy as np
import time


def VecTose3(V):
    """
    Returns 4x4 se(3) matrix given a 6x1 vector.
    6x1 vector is in the form of [ v omega ]^T
    """

    V_hat = SX(4, 4)
    V_hat[0:3, 0:3] = skew(V[3:6])
    V_hat[0:3, 3] = V[0:3]

    return V_hat


def se3ToVec(V_hat):
    """
    Returns 6x1 vector given a 4x4 se(3) matrix.
    6x1 vector is in the form of [ v omega ]^T
    """

    V = SX(6, 1)
    V[3:6] = inv_skew(V_hat[0:3, 0:3])
    V[0:3] = V_hat[0:3, 3]

    return V


def SE3ToVecT(SE3):
    """
    Returns 12x1 vector given a 4x4 SE(3) matrix.
    """

    vecT = SX(12, 1)
    vecT[3:12] = reshape(SE3[0:3, 0:3], 9, 1)
    vecT[0:3] = SE3[0:3, 3]

    return vecT


def VecTtoSE3(vecT):
    """
    Returns 4x4 SE(3) matrix given 12x1 matrix.
    """

    SE3 = SX(4, 4)
    SE3[0:3, 0:3] = reshape(vecT[3:12], 3, 3)
    SE3[0:3, 3] = vecT[0:3]
    SE3[3, 3] = 1

    return SE3


def SXtoNumpy_vec(SX_):

    NP = np.zeros(SX_.shape)

    for i in range(SX_.shape[0]):

        NP[i] = float(SX_[i])

    return NP


def SXtoNumpy_mat(SX_):

    s1, s2 = SX_.shape

    SX__ = reshape(SX_, SX_.shape[0] * SX_.shape[1], 1)

    NP = np.zeros((SX_.shape[0] * SX_.shape[1], 1))

    for i in range(SX__.shape[0]):

        NP[i] = float(SX_[i])

    NP = np.transpose(np.reshape(NP, (s2, s1)))

    return NP

def quaternion_multiply(Q0,Q1):
    """
    Multiplies two quaternions.
 
    Input
    :param Q0: A 4 element array containing the first quaternion (q01,q11,q21,q31) 
    :param Q1: A 4 element array containing the second quaternion (q02,q12,q22,q32) 
 
    Output
    :return: A 4 element array containing the final quaternion (q03,q13,q23,q33) 
 
    """
    # Extract the values from Q0
    x0 = Q0[0]
    y0 = Q0[1]
    z0 = Q0[2]
    w0 = Q0[3]
     
    # Extract the values from Q1
    x1 = Q1[0]
    y1 = Q1[1]
    z1 = Q1[2]
    w1 = Q1[3]
     
    # Computer the product of the two quaternions, term by term
    Q0Q1_w = w0 * w1 - x0 * x1 - y0 * y1 - z0 * z1
    Q0Q1_x = w0 * x1 + x0 * w1 + y0 * z1 - z0 * y1
    Q0Q1_y = w0 * y1 - x0 * z1 + y0 * w1 + z0 * x1
    Q0Q1_z = w0 * z1 + x0 * y1 - y0 * x1 + z0 * w1
     
    # Create a 4 element array containing the final quaternion
    final_quaternion = np.array([Q0Q1_x, Q0Q1_y, Q0Q1_z, Q0Q1_w])
     
    # Return a 4 element array containing the final quaternion (q02,q12,q22,q32) 
    return final_quaternion

    # def quaternion_inverse(Q): 

    #     """Q: quaternion as a numpy.array.
    #     Returns: Q^-1"""

    #     w = Q[0]
    #     x = Q[1]
    #     y = Q[2]
    #     z = Q[3]



if __name__ == "__main__":

    # Test functions

    b = SX([1, 2, 3])
    a = SX([[0, 0, 9.42, 1], [0, 0, 0, 2], [-9.4546, 0, 0, 3]])
    print(a)
    # a = SX([[10, 11, 12, 1, 0, 0, 0, 1, 0, 0, 0, 1]])
    # b = (VecTtoSE3(a))
    # print(inv(b)@b)

    print(SXtoNumpy_vec(b))
    a = SXtoNumpy_mat(a)
    print(a)
