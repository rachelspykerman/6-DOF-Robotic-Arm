#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
import numpy as np
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        
        # Define DH param symbols
        # link_offset_i: signed distance from x_(i-1) to x_i along z_i
        # link_length_(i-1): distance from z_(i-1) to z_i along x_i, where x_i is perpendicular to both z_(i-1) and z_i
        # twist_angle_(i-1): angle between z_(i-1) and z_i about x_(i-1) in the right hand sense
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8') # link_offset_i
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7') # link_length_(i-1)
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')# twist_angle_(i-1)

        # theta_i: joint angle between x_(i-1) and x_i about z_i in the right hand sense
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # theta_i

        # Modified DH params
        s = {alpha0:     0,   a0:      0,   d1:  0.75,
         alpha1: -pi/2,   a1:   0.35,   d2:     0,   q2: q2-pi/2,
         alpha2:     0,   a2:   1.25,   d3:     0,
         alpha3: -pi/2,   a3: -0.054,   d4:   1.5,
         alpha4:  pi/2,   a4:      0,   d5:     0,
         alpha5: -pi/2,   a5:      0,   d6:     0,
         alpha6:     0,   a6:      0,   d7: 0.303,   q7:       0}

        # Create individual transformation matrices
        T0_1 = Transformation_Matrix(q1, alpha0, a0, d1)
        T0_1 = T0_1.subs(s)

        T1_2 = Transformation_Matrix(q2, alpha1, a1, d2)
        T1_2 = T1_2.subs(s)

        T2_3 = Transformation_Matrix(q3, alpha2, a2, d3)
        T2_3 = T2_3.subs(s)

        R0_3 = simplify(T0_1 * T1_2 * T2_3)

        # Intrinsic rotation correction for end-effector transformation matrix 
        R_z_intrinsic = Matrix([[    cos(pi),    -sin(pi),               0],
                            [    sin(pi),     cos(pi),               0],
                            [          0,           0,               1]])
        R_y_intrinsic = Matrix([[ cos(-pi/2),           0,      sin(-pi/2)],
                            [          0,           1,               0],
                            [-sin(-pi/2),           0,      cos(-pi/2)]])
        ZY_intrinsic_rot = R_z_intrinsic * R_y_intrinsic
        
        
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()
            
            # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

	    # Calculate transformation matrix from base link to end-effector
	    endeffector_trans = Matrix([[px],
		                        [py],
		                        [pz]])
	    R_z_extrinsic = Matrix([[   cos(yaw),   -sin(yaw),               0],
		                    [   sin(yaw),    cos(yaw),               0],
		                    [          0,           0,               1]])
	    R_y_extrinsic = Matrix([[ cos(pitch),           0,      sin(pitch)],
		                    [          0,           1,               0],
		                    [-sin(pitch),           0,      cos(pitch)]])
	    R_x_extrinsic = Matrix([[          1,           0,               0],
		                    [          0,   cos(roll),      -sin(roll)],
		                    [          0,   sin(roll),      cos(roll)]])
	    XYZ_extrinsic_rot = R_z_extrinsic * R_y_extrinsic * R_x_extrinsic
	    R0_6 = XYZ_extrinsic_rot * ZY_intrinsic_rot

	    # Find Wrist Center Location
	    d_7 = s[d7]
	    wx = px - (d_7 * R0_6[0, 2])
	    wy = py - (d_7 * R0_6[1, 2])
	    wz = pz - (d_7 * R0_6[2, 2])

	    # Finding theta 1-3
	    a_3 = s[a3]
	    d_4 = s[d4]
	    d_1 = s[d1]
	    a_1 = s[a1]
	    a_2 = s[a2]

	    theta1 = (atan2(wy, wx)).evalf()
	    
	    s1 = sqrt(wx**2 + wy**2) - a_1
	    s2 = wz - d_1
	    s3 = sqrt(s2**2 + s1**2)
	    s4 = sqrt(a_3**2 + d_4**2)
	    beta1 = atan2(s2, s1)

	    D2 = (a_2**2 + s3**2 - s4**2) / (2 * a_2 * s3)
	    beta2 = atan2(sqrt(1 - D2**2), D2)
	    
	    D3 = (a_2**2 + s4**2 - s3**2) / (2 * a_2 * s4)
	    beta3 = atan2(sqrt(1 - D3**2), D3)
	    
	    beta4 = atan2(-a_3, d_4)

	    theta2 = ((pi / 2) - beta2 - beta1).evalf()
	    theta3 = ((pi / 2) - beta4 - beta3).evalf()

	    # Finding theta 4-6

	    R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
	    R0_3 = R0_3[0:3, 0:3]
	    R0_3_inv = R0_3.inv("LU")
	    R3_6 = R0_3_inv * R0_6

	    r13 = R3_6[0, 2]
	    r33 = R3_6[2, 2]
	    r23 = R3_6[1, 2]
	    r21 = R3_6[1, 0]
	    r22 = R3_6[1, 1]
	    r12 = R3_6[0, 1]
	    r32 = R3_6[2, 1]
	    
	    theta5 = (atan2(sqrt(r13**2 + r33**2), r23)).evalf()
	    
	    if (sin(theta5) < 0):
	        print("BELOW!!!")
	        theta4 = (atan2(-r33, r13)).evalf()
	        theta6 = (atan2(r22, -r21)).evalf()
	    elif (theta5 == 0):
	        print("EQUAL!!!")
	        theta4 = 0
	        theta6 = (atan2(-r12, -r32)).evalf()
	    else:
	        print("ELSE!!!!!")
	        theta4 = (atan2(r33, -r13)).evalf()
	        theta6 = (atan2(-r22, r21)).evalf()
	        
	    while (theta4 > pi):
	        theta4 = theta4 - 2*pi
	    while (theta4 < -pi):
	        theta4 = 2*pi + theta4
	        
	    while (theta5 > pi):
	        theta5 = theta5 - 2*pi
	    while (theta5 < -pi):
	        theta5 = 2*pi + theta5
	        
	    while (theta6 > pi):
	        theta6 = theta6 - 2*pi
	    while (theta6 < -pi):
	        theta6 = 2*pi + theta6
	        

	    print('wx: ' + str(wx))
	    print('wy: ' + str(wy))
	    print('wz: ' + str(wz))
	    print('theta1: ' + str(theta1))
	    print('theta2: ' + str(theta2))
	    print('theta3: ' + str(theta3))
	    print('theta4: ' + str(theta4))
	    print('theta5: ' + str(theta5))
	    print('theta6: ' + str(theta6))

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)

def Transformation_Matrix(q, alpha, a, d):
    T = Matrix([[           cos(q),           -sin(q),           0,             a],
                [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                [                0,                 0,           0,             1]])
    return(T)

def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
