import turtle
import math

import numpy as np
from sympy import solve
from sympy.abc import x
#print(np.__version__)
import cmath
from tkinter import *  # Python 3

import time
#TODO: look into math.exp() or numpy exp()
# create Readme
def v_infty(r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    '''
    i= 1j
    e = math.e
    theta = angleDeg / 180 * math.pi
    return i/(1 - r*e**(i*theta))

def v_0(r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    '''
    i= 1j
    e = math.e
    theta = angleDeg / 180 * math.pi
    return i * (1+r*e**(-i*theta)) / (1-r**2)

def v_(N, r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    '''
    i= 1j
    e = math.e
    theta = angleDeg / 180 * math.pi
    v_N = i * ( (1+r*e**(-i*theta))/(1-r**2) * (r*e**(i*theta))**N + (1-(r*e**(i*theta))**N)/(1-r*e**(i*theta)) )
    return v_N

def w_(N, r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    '''
    i= 1j
    e = math.e
    theta = angleDeg / 180 * math.pi
    w_N = i * (1-(r*e**(i*theta))**N) / (1-r*e**(i*theta))
    return w_N 

def C(r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    !here m=4 so type-4 for binary fractal tree
    '''
    i= 1j
    return 2 * v_infty(r, angleDeg) - i

def binary_type_3_C(r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    !here m=3 so type-3 for binary fractal tree
    this goes horizontally to the right.
    '''
    return -2*(v_infty(r, angleDeg).real)*math.sqrt(3)

def binary_type_2_C(r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    !here m=2 so type-2 for binary fractal tree
    this goes up
    '''
    return -4*(v_infty(r, angleDeg).real)*1j

def binary_type_2_D(r, angleDeg):
    '''
    important: branchLen is assumed to be 1!
    !here m=2 so type-2 for binary fractal tree
    '''
    e = math.e
    i = 1j
    #translate everything by -v_infty
    v_0_prime = v_0(r, angleDeg) - v_infty(r, angleDeg)
    v_1_prime = v_(1, r, angleDeg) - v_infty(r, angleDeg)

    return v_0_prime * e ** (i*math.pi/2) - v_1_prime * e ** (-i*math.pi/2)
    

def arrowhead(t, small_angle, length):
    og_down = t.isdown()
    turn_angle = 180-small_angle
    
    t.rt(turn_angle)
    t.down()
    t.fd(length)
    t.up()
    t.bk(length)
    t.lt(2 * turn_angle)
    t.down()
    t.fd(length)
    t.up()
    t.bk(length)
    t.rt(turn_angle)
    if og_down:
        t.pendown()

def dashed(tt_turtle_obj, length, dashLen):
    num = length // (2*dashLen)
    remainder = length % (2*dashLen)
    for i in range(int(num)):
        tt_turtle_obj.penup()
        tt_turtle_obj.forward(dashLen)
        tt_turtle_obj.pendown()
        tt_turtle_obj.forward(dashLen)
    tt_turtle_obj.penup()
    tt_turtle_obj.fd(remainder)

def drawCircle(t, r, angleDeg):
    og_isdown = t.isdown()
    t.penup()
    t.fd(r)
    t.lt(90)
    t.pendown()
    t.circle(r, angleDeg, 30)
    t.penup()
    t.circle(r, 360-angleDeg, 30)
    t.rt(90)
    t.bk(r)
    t.pendown()
    if not og_isdown:
        t.penup()
    

def tree(branchLen,t,angleDeg,r,stopLen):
    
    if branchLen > stopLen:
        #t.penup()
        t.forward(branchLen)
        # This draws the right portion of the tree
        t.right(angleDeg)
        # Where we make a recursive call
        tree(branchLen *r,t,angleDeg,r,stopLen)
        # Left portion of the tree
        t.left(2 * angleDeg)
        # Another recursive call
        tree(branchLen *r,t,angleDeg,r,stopLen)
        t.right(angleDeg)
        '''
        if branchLen>10:
            t.color("red")
            arrowhead(t, 20, branchLen/10)'''
        t.up()
        t.backward(branchLen)
        t.down()
        '''
        if branchLen > 10:
            t.color("blue")
            drawCircle(t, branchLen/30, 360)
            t.color("black")'''
    '''else:
        t.pendown()
        drawCircle(t, branchLen/10, 360)
        t.penup()'''

def tree_with_stopK(branchLen,t,angleDeg,r,stopK, annotate):
    og_pos, og_heading, og_width = t.pos(), t.heading(), t.width()
    
    
    if stopK >= 0:
        if annotate and stopK == 0:
            t.pencolor("red")
            #t.pensize(1.5)
        else:
            t.pencolor("black")
            #t.pensize(.25)
        t.forward(branchLen)
        if annotate:
            arrowhead(t, 20, branchLen/6)
        # This draws the right portion of the tree
        t.right(angleDeg)
        # Where we make a recursive call
        tree_with_stopK(branchLen *r,t,angleDeg,r,stopK - 1, annotate)
        # Left portion of the tree
        t.left(angleDeg)
        if annotate and stopK >= 1:
            #t.pensize(1)
            t.pencolor("blue")
            dashed(t,branchLen,.05*branchLen+og_width)
            t.penup()
            t.bk(branchLen)
            t.pendown()
            t.rt(angleDeg)
            drawCircle(t, .2*branchLen, 2*angleDeg)
            t.lt(angleDeg)
            t.pendown()
        t.left(angleDeg)
        # Another recursive call
        tree_with_stopK(branchLen *r,t,angleDeg,r,stopK - 1, annotate)
        t.right(angleDeg)
        t.up()
        t.backward(branchLen)

def ternary_tree(branchLen,t,angleDeg,r,stopK, annotate):
    t.down()
    if stopK >= 0:
        if annotate and stopK == 0:
            t.pencolor("red")
            t.pensize(1.5)
        else:
            t.pencolor("black")
            t.pensize(.25)
        t.forward(branchLen)
        if annotate:
            arrowhead(t, 20, branchLen/6)
        # This draws the center portion of the tree
        ternary_tree(branchLen *r,t,angleDeg,r,stopK - 1, annotate)
        # This draws the right portion of the tree
        t.right(angleDeg)
        # Where we make a recursive call
        ternary_tree(branchLen *r,t,angleDeg,r,stopK - 1, annotate)
        # Left portion of the tree
        t.left(angleDeg)
        if annotate and stopK >= 1:
            t.pensize(1)
            t.pencolor("blue")
            dashed(t,branchLen,4)
            t.penup()
            t.bk(branchLen)
            t.pendown()
            t.rt(angleDeg)
            drawCircle(t, 15, 2*angleDeg)
            t.lt(angleDeg)
            t.pendown()
        t.left(angleDeg)
        # Another recursive call
        ternary_tree(branchLen *r,t,angleDeg,r,stopK - 1, annotate)
        t.right(angleDeg)
        t.up()
        t.backward(branchLen)


def nary_tree(branchLen,t,angleDeg,r,stopK, annotate, n):
    '''
    angleDeg is the maximum turn-angle,
    which represents the greatest extent the left-most and right-most branches will turn
    it will be assumed that everything in between is equally spaced and has the same ratio r.
    n must be at least 2
    '''
    if stopK >= 0:
        if annotate and stopK == 0:
            t.pencolor("red")
            t.pensize(t.width()+1.25)
        else:
            t.pencolor("black")
            t.pensize(t.width())
        t.forward(branchLen)
        if annotate:
            arrowhead(t, 20, branchLen/6)

        # This points in the direction of the right-most portion of the tree
        t.right(angleDeg)
        # math
        tot = angleDeg*2
        spacing = tot/(n-1)
        for kk in range(n-1):
            nary_tree(branchLen *r,t,angleDeg,r,stopK - 1, annotate, n)
            t.left(spacing)
        # Left-most portion of the tree
        nary_tree(branchLen *r,t,angleDeg,r,stopK - 1, annotate, n)

        t.right(angleDeg)
        t.down()
        t.backward(branchLen)
        t.up()

        
def pruned_nary_tree(branchLen,t,angleDeg,r,stopK, annotate, n):
    '''
    angleDeg is the maximum turn-angle in deg,
    which represents the greatest extent the left-most and right-most branches will turn
    it will be assumed that everything in between is equally spaced and has the same ratio r.
    n must be at least 2
    '''
    if stopK >= 0 and n>=2:
        if annotate and stopK == 0:
            t.pencolor("red")
            t.pensize(1.5)
        else:
            t.pencolor("black")
            t.pensize(.25)
        t.forward(branchLen)
        if annotate:
            arrowhead(t, 20, branchLen/6)

        # This points in the direction of the right-most portion of the tree
        t.right(angleDeg)
        # math
        tot = angleDeg*2
        spacing = tot/(n-1)
        for kk in range(n-1):
            pruned_nary_tree(branchLen *r,t,angleDeg,r,stopK - 1, annotate, n-1)
            t.left(spacing)
        # Left-most portion of the tree
        pruned_nary_tree(branchLen *r,t,angleDeg,r,stopK - 1, annotate, n-1)

        t.right(angleDeg)
        t.backward(branchLen)

def crown(branchLen, t, angleDeg, r, stopK):
    #drawCircle(t,branchLen/10,360)
    e = math.e
    i = 1j
    theta = angleDeg * math.pi / 180

    v_infty = branchLen*i/(1 - r*e**(i*theta))
    v_0 = branchLen * i * (1+r*e**(-i*theta)) / (1-r**2)
    #print(v_0)
    
    t.penup()
    t.fd(branchLen)
    t.lt(angleDeg)
    #t.pendown()
    #drawCircle(t,branchLen/10,360)
    for j in range(2):
        t.penup()
        for i in range(stopK):
            b=branchLen * (r ** (i+1))
            t.fd(b)
            t.lt(angleDeg-j*(2*angleDeg))
            #drawCircle(t,b/10,360)
        t.pendown()
        for i in range(stopK):
            t.rt(angleDeg-j*(2*angleDeg))
            t.bk(branchLen * (r ** (stopK - i)))
        t.rt(2*angleDeg)
    t.lt(3*angleDeg)
    t.penup()
    t.bk(branchLen)
    t.fd(v_0.imag)
    t.rt(90)
    topLen = 2 * v_0.real
    for j in range(2):
        t.up()
        t.fd(v_0.real)
        t.rt(angleDeg-j*(2*angleDeg))
        #drawCircle(t,topLen/10, 360)
        for i in range(stopK):
            top=topLen * (r ** (i+1))
            t.fd(top)
            t.rt(angleDeg-j*(2*angleDeg))
            #drawCircle(t,top/10, 360)
        t.down()
        for i in range(stopK):
            t.lt(angleDeg-j*(2*angleDeg))
            t.bk(topLen * (r ** (stopK - i)))
        t.lt(angleDeg-j*(2*angleDeg))
        t.bk(v_0.real)
        t.rt(180)
    t.lt(90)
    t.penup()
    t.bk(v_0.imag)
    t.pendown()

def log_crown_BFT(branchLen,t,angleDeg, r, stopK, annotate, step, cutFirstSeg=True):
    og_pos, og_heading = t.pos(), t.heading()

    log_spiral(branchLen, t, angleDeg, r, stopK, annotate, step, cutFirstSeg)
    log_spiral(branchLen, t, -angleDeg, r, stopK, annotate, step, cutFirstSeg)
    v_0_b = v_0(r, angleDeg) * branchLen
    t.up()
    t.fd(v_0_b.imag)
    t.lt(90)
    t.bk(v_0_b.real)
    t.down()
    topLen = 2 * v_0_b.real
    log_spiral(topLen, t, angleDeg, r, stopK, annotate, step, cutFirstSeg)
    t.up()
    t.fd(topLen)
    t.rt(180)
    t.down()
    log_spiral(topLen, t, -angleDeg, r, stopK, annotate, step, cutFirstSeg)

    t.up()
    t.goto(og_pos)
    t.setheading(og_heading)
    t.down()
    t.pencolor("black")
    t.fd(branchLen)


def log_spiral(branchLen, t, angleDeg, r, stopK, annotate, step, cutFirstSeg):
    og_pos, og_heading = t.pos(), t.heading()
    og_color, og_isdown = t.pencolor(), t.isdown()
    
    e = math.e
    i = 1j
    theta = angleDeg * math.pi / 180
    def v_(N):
        v_N = i * ( (1+r*e**(-i*theta))/(1-r**2) * (r*e**(i*theta))**N + (1-(r*e**(i*theta))**N)/(1-r*e**(i*theta)) )
        return v_N
    def w_(N):
        w_N = i * (1-(r*e**(i*theta))**N) / (1-r*e**(i*theta))
        return w_N
    
    v_infty = branchLen*i/(1 - r*e**(i*theta))

    #this draws the relevant portion of a fractal tree approximating the log spiral
    if annotate:
        t.down()
        t.pencolor("red")
        for kk in range(stopK):
            b=branchLen * (r ** kk)
            t.fd(b)
            t.lt(angleDeg)
        t.up()
        t.goto(og_pos)
        t.setheading(og_heading)
    
    #print(i)
    #t.goto(v_infty.real+og_pos[0],v_infty.imag+og_pos[1])
    #drawCircle(t,branchLen/10,360)
    
    og_heading_rad = og_heading * math.pi/180
    
    t.down()
    width=step
    t.pencolor("blue")
    for jj in np.arange(0 if not cutFirstSeg else 1,stopK+width,width):
        w_N_transformed = branchLen*w_(jj) * e **(i * (og_heading_rad-math.pi/2))
        if cutFirstSeg and (jj==1):
            t.up()
        else:
            t.down()
        t.goto(w_N_transformed.real+og_pos[0],w_N_transformed.imag+og_pos[1])
        #drawCircle(t,branchLen/40,360)
    t.penup()
    t.goto(og_pos)
    t.setheading(og_heading)
    if og_isdown:
        t.pendown()
    t.pencolor(og_color)

    
def binary_type_m_unit_parameters(m,t):
    #returns r and theta needed to generate a binary fractal tree given m and t value
    if not m == 4:
        theta = math.pi*(2/m - 1/2)/t
        r = solve(x**(-t-1) - x**(-t+1)-2*math.sin(theta), x)
        print("for m =" , m, "and t =", t)
        print("theta:", theta)
        print("r:", r)
        print()
        return r, theta


def binary_type_4_unit_radius(angleDeg):
    "returns r given angleDeg in degrees for a type-4 binary fractal tree"
    theta = angleDeg * math.pi / 180
    return .5 * (math.sqrt(2)*math.sqrt(3-math.cos(2*theta))-2*math.sin(theta))

def F(branchLen,t,angleDeg,r,stopK, annotate, n, m, isTree=True, isCrown=True, treeColor="black", crownColor=(55,105,255)):
    turtle.colormode(255)
    og_pos, og_heading = t.pos(), t.heading()
    e = math.e
    i = 1j
    v_inf = v_infty(r, angleDeg) * branchLen
    v_inf_phase = cmath.phase(v_inf)
    v_inf_abs = abs(v_inf)
    t.up()
    if isTree:
        t.pencolor(treeColor)
        for kk in range(m):
            t.setheading(kk*360/m)
            t.lt(math.degrees(v_inf_phase - math.pi/2))
            t.fd(-v_inf_abs)
            t.rt(math.degrees(v_inf_phase - math.pi/2))
            t.down()
            
            nary_tree(branchLen,t,angleDeg,r,stopK, annotate, n)

            
            t.up()
            t.lt(math.degrees(v_inf_phase - math.pi/2))
            t.fd(v_inf_abs)
    if isCrown:
        t.pencolor(crownColor)
        for kk in range(m):
            t.setheading(kk*360/m)
            t.lt(math.degrees(v_inf_phase - math.pi/2))
            t.fd(-v_inf_abs)
            t.rt(math.degrees(v_inf_phase - math.pi/2))
            t.down()
            
            crown(branchLen, t, angleDeg, r, stopK)

            t.up()
            t.lt(math.degrees(v_inf_phase - math.pi/2))
            t.fd(v_inf_abs)
            
    t.setheading(og_heading)
    t.goto(og_pos)
        
def linear_tile(branchLen,t,angleDeg,r,stopK, annotate, n, m, linear_c, linear_num, isTree=True, isCrown=True):
    '''linear_c is a complex number descirbing the displacement between adjacent units in the linear tiling
    linear_c assumes that trunk_len = 1
    linear_num is the number of units in the linear tiling
    '''
    og_pos, og_heading, og_width = t.pos(), t.heading(), t.width()
    fdLen = abs(linear_c)*branchLen
    t.setheading(np.degrees(cmath.phase(linear_c)))

    t.pensize(og_width*2)
    t.pencolor("red")
    t.down()
    
    for kk in range(linear_num-1):
        t.fd(fdLen)
        arrowhead(t, 20, fdLen/5)
    t.up()
    t.goto(og_pos)
    t.width(og_width)
    
    for kk in range(linear_num):
        F(branchLen,t,angleDeg,r,stopK, annotate, n, m, isTree, isCrown)
        t.fd(fdLen)
        
    t.goto(og_pos)
    t.setheading(og_heading)
    t.width(og_width)

def two_d_tile(branchLen,t,angleDeg,r,stopK, annotate, n, m, linear_c, two_d_c, linear_num, two_d_num, isTree=True, isCrown=True):
    '''linear_c is a complex number descirbing the displacement between adjacent units in the linear tiling
    linear_c assumes that trunk_len = 1
    linear_num is the number of units in the linear tiling
    similar thing for two_d_c and two_d_num
    '''
    print("linear_c:",linear_c)
    print("two_d_c:",two_d_c)
    og_pos, og_heading, og_width = t.pos(), t.heading(), t.width()
    fdLen = abs(two_d_c)*branchLen
    
    t.setheading(np.degrees(cmath.phase(two_d_c)))
    t.pensize(og_width*2)
    t.pencolor("lime")
    t.down()
    
    for kk in range(two_d_num-1):
        t.fd(fdLen)
        arrowhead(t, 20, fdLen/5)
    t.up()
    t.goto(og_pos)
    t.width(og_width)
    
    for kk in range(two_d_num):
        linear_tile(branchLen,t,angleDeg,r,stopK, annotate, n, m, linear_c, linear_num, isTree, isCrown)
        t.fd(fdLen)
        
    t.goto(og_pos)
    t.setheading(og_heading)
    t.width(og_width)

def different_tree_comparisons(t, scale):
    #making the most simple type-3
    tree_with_stopK(scale*100,t,60,1/(2 ** 1),3, False)
    t.rt(90)
    t.fd(200)
    t.lt(90)
    t.down()
    nary_tree(scale*100,t,60,1/(2 ** 1),3, False,2)
    t.rt(90)
    t.fd(200)
    t.lt(90)
    t.down()    
    ternary_tree(scale*100,t,60,1/(2 ** 1),3, False)
    t.rt(90)
    t.fd(200)
    t.lt(90)
    t.down()
    nary_tree(scale*100,t,60,1/(2 ** 1),3, False, 3)
    
def different_n_ary_trees(t):
    for i in range(5):
        n=i+2
        t.up()
        t.rt(90)
        t.fd(260)
        t.lt(90)
        t.down()
        nary_tree(135,t,60,1/(2 ** 1),6, False, n)

def ternary_type_m_unit(t):
    n=5
    th = math.pi*(1-4/n)
    ratio = 1 - (2*(1-math.cos(th)))**.5
    nary_tree(270,t,th/math.pi*180,ratio,4, False, 3)

def probe_binary_type_m_unit_parameters():
    r_theta_values = [(("M","T"),("r","theta"))]
    print(r_theta_values)
    for M in range(1,6):
        if M == 4:
            continue
        for T in range(-6,6):
            if T == 0:
                continue
            r, theta = binary_type_m_unit_parameters(m=M,t=T)
            r_theta_values.append(((M,T),(r, theta)))
            
    print(r_theta_values)
    return r_theta_values

def valid_binary_type_m_unit_parameters():
    counter =-1
    valid_r_theta_values = [(("M","T"),("r","theta"))]
    for element in r_theta_values:
        counter += 1
        if counter == 0:
            continue
        #print(element)
        ((M,T),(r, theta)) = element
        #print(((M,T),(r, theta)))
        for i in range(len(r)):
            ratio = r[i]
            print(ratio)
            if (abs(ratio.imag) > 0.0000001):
                print("complex ratio not allowed")
                continue
            elif (abs(ratio.real) >= 1) or (ratio.real <= 0):
                print("ratio is real but out of bounds")
                continue
            else:
                ratio = ratio.real
                print("Good")
                valid_r_theta_values.append(((M,T),(ratio, theta)))
    print(valid_r_theta_values)
    return valid_r_theta_values

def draw_valid_binary_type_m_unit_parameters(t, valid_r_theta_values):
    t.up()
    bL = 20
    t.goto(-925,500)
    t.setheading(0)
    counter = -1
    num_in_row = 5
    for elem in valid_r_theta_values:
        counter += 1
        if counter == 0:
            continue
        t.goto(-700+(counter % num_in_row)*6*bL, 240-(counter // num_in_row)*6.1*bL)
        ((M,T),(r,theta)) = elem
        print(((M,T),(r,theta)))
        angleDeg = theta * 180/math.pi
        ratio = r
        
        info = "m=" +str(M)+" t=" + str(T) + "\nangleDeg="+str(round(angleDeg,2)) + "° r=" + str(round(ratio,4))
        t.write(info, False, align="center")
        t.setheading(90)
        v_inf = v_infty(ratio, angleDeg)
        print(v_inf)
        t.fd(7*bL/(1 if M == 1 else 2))
        t.down()
        trunk_len = bL*5*(1-ratio)/(1 if M == 1 else 2)
        #tip_len = 4
        #k_max = int(math.log(tip_len/trunk_len)/math.log(ratio))
        #error = [12,12,12,12,12]
        #error = [20,20,20,20,20]
        #k_max = int(math.log(error[M-1]*(trunk_len/(1-ratio))**-1)/math.log(r))
        k_max = 8
        #F(trunk_len,t,angleDeg,ratio,k_max, False, 2,M, isTree =True, isCrown=False)
        #F(trunk_len,t,angleDeg,ratio,int(k_max*5), False, 2,M, isTree =False, isCrown=True)
        two_d_tile(
            branchLen=trunk_len,
            t=t,
            angleDeg=angleDeg,
            r=ratio,
            stopK=int(k_max),
            annotate=False,
            n=2,
            m=M,
            linear_c=1-1j,two_d_c = 1+2j,
            linear_num = 1, two_d_num = 1,
            isTree=True,
            isCrown=True
            )
        t.up()
    ts = t.getscreen()
    ts.getcanvas().postscript(file="draw_valid_binary_type_m_unit_parameters.eps")
        
def draw_binary_type_4_tiling_eps(t):
    t.up()
    t.setheading(0)
    t.goto(-600,0)
    t.down()
    #turtle.update()
    for angleDeg in (30,):#(10,20,30,50,90):
        r = binary_type_4_unit_radius(angleDeg)
        bL = 25/abs(v_infty(r, angleDeg).real)
        #F(50/abs(v_infty(r, angleDeg).real),t,angleDeg,r,9, False, 2,4, isTree =True, isCrown=False)
        #F(50/abs(v_infty(r, angleDeg).real),t,angleDeg,r,27, False, 2,4, isTree =False, isCrown=True)
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=7,
            annotate=False,
            n=2,
            m=4,
            linear_c=C(r, angleDeg),
            two_d_c = C(r, angleDeg) * math.e ** (-1j*math.pi/2),
            linear_num=3,
            two_d_num=3,
            isTree=True,
            isCrown=False)
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=27,
            annotate=False,
            n=2,
            m=4,
            linear_c=C(r, angleDeg),
            two_d_c = C(r, angleDeg) * math.e ** (-1j*math.pi/2),
            linear_num=3,
            two_d_num=3,
            isTree=False,
            isCrown=True)
        og_heading = t.heading()

        info_offset = 150
        t.setheading(-90)
        t.fd(info_offset)
        info = "m=4 t=0\nangleDeg="+str(round(angleDeg,2)) + "° r=" + str(round(r,4))
        t.pencolor("black")
        t.write(info, False, align="center")
        t.setheading(90)
        t.fd(info_offset)
        t.setheading(og_heading)
        t.fd(300)
        print("phase",cmath.phase(v_infty(r, angleDeg)))
    ts = t.getscreen()
    ts.getcanvas().postscript(file="binary_type-4_tiling.eps")

def draw_binary_type_4_units_eps(t):
    t.up()
    t.setheading(0)
    #turtle.update()
    counter = -1
    num_in_row = 3
    for angleDeg in (10,20,30,45,60,90):
        counter += 1
        
        t.goto(-600+(counter % num_in_row)*300, 270-(counter // num_in_row)*370)

        r = binary_type_4_unit_radius(angleDeg)
        bL = 25/abs(v_infty(r, angleDeg).real)
        F(60/abs(v_infty(r, angleDeg).real),t,angleDeg,r,9, False, 2,4, isTree =True, isCrown=False)
        F(60/abs(v_infty(r, angleDeg).real),t,angleDeg,r,27, False, 2,4, isTree =False, isCrown=True)

        og_heading = t.heading()
        info_offset = 190
        t.setheading(-90)
        t.up()
        t.fd(info_offset)
        info = "m=4 t=0\nangle="+str(round(angleDeg,2)) + "° r=" + str(round(r,4))
        t.pencolor("black")
        t.write(info, False, align="center")
        t.setheading(90)
        t.up()
        t.fd(info_offset)
        t.setheading(og_heading)
        print("phase",cmath.phase(v_infty(r, angleDeg)))
    ts = t.getscreen()
    ts.getcanvas().postscript(file="binary_type-4_units.eps")

def draw_binary_type_4_tiling_visual_proof_eps(t, sc):
    sc.setup(width=.5,height=0.9, startx=0, starty=0)

    t.up()
    t.goto(0,-150)
    t.setheading(0)
    counter = -1
    num_in_row = 2
    angleDeg = 30
    names = ("two_units", "linear_tiling", "double_linear_tiling", "complete_tiling")
    for (linear_num, two_d_num) in ((2,1),(3,1),(3,2),(3,3)):
        counter += 1
        #t.goto(-550+(counter % num_in_row)*500, 220-(counter // num_in_row)*400)

        r = binary_type_4_unit_radius(angleDeg)
        bL = 45/abs(v_infty(r, angleDeg).real)
        #F(60/abs(v_infty(r, angleDeg).real),t,angleDeg,r,9, False, 2,4, isTree =True, isCrown=False)
        #F(60/abs(v_infty(r, angleDeg).real),t,angleDeg,r,27, False, 2,4, isTree =False, isCrown=True)
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=10,
            annotate=False,
            n=2,
            m=4,
            linear_c=C(r, angleDeg),
            two_d_c = C(r, angleDeg) * math.e ** (-1j*math.pi/2),
            linear_num=linear_num,
            two_d_num=two_d_num,
            isTree=True,
            isCrown=False)
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=27,
            annotate=False,
            n=2,
            m=4,
            linear_c=C(r, angleDeg),
            two_d_c = C(r, angleDeg) * math.e ** (-1j*math.pi/2),
            linear_num=linear_num,
            two_d_num=two_d_num,
            isTree=False,
            isCrown=True)
        print("phase",cmath.phase(v_infty(r, angleDeg)))
        ts = t.getscreen()
        ts.getcanvas().postscript(file="binary_type-4_" + names[counter] +".eps")
        #turtle.update()
        #time.sleep(1)
        t.clear()

def draw_binary_type_2_tiling_eps(t, valid_r_theta_values, stopK=7):
    t.up()
    t.goto(0,500)
    t.setheading(0)
    counter = -1
    num_in_row = 3
    scale=40
    for elem in valid_r_theta_values:
        if isinstance(elem, str):
            continue
        ((M,T),(r,theta)) = elem
        if not M == 2:
            continue
        
        counter += 1
        if counter > 0:
            break
        angleDeg = theta * 180/math.pi
        bL = 15/abs(v_infty(r, angleDeg).real)*scale/40
        
        t.goto(-500+(counter % num_in_row)*11*scale, 140-(counter // num_in_row)*9*scale)
        t.goto(500,-250)
        
        og_heading = t.heading()
        info_offset = 100*scale/40
        t.setheading(-90)
        t.fd(info_offset)
        info = "m=" +str(M)+" t=" + str(T) + "\nangle="+str(round(angleDeg,2)) + "° r=" + str(round(r,4)) 
        t.pencolor("black")
        t.write(info, False, align="center")
        t.setheading(90)
        t.fd(info_offset)
        t.setheading(og_heading)
        
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=stopK,
            annotate=False,
            n=2,
            m=2,
            linear_c=binary_type_2_C(r, angleDeg),
            two_d_c = binary_type_2_D(r, angleDeg),
            linear_num=1,
            two_d_num=1,
            isTree=True,
            isCrown=False)
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=stopK,
            annotate=False,
            n=2,
            m=2,
            linear_c=binary_type_2_C(r, angleDeg),
            two_d_c = binary_type_2_D(r, angleDeg),
            linear_num=1,
            two_d_num=1,
            isTree=False,
            isCrown=True)
        
    og_heading = t.heading()
    ts = t.getscreen()
    ts.getcanvas().postscript(file="binary_type-2_tiling.eps")
''' #interesting edge case where there is parallels between (type2, t=1) and type 3 linearly independent basis vectors. 
    t.up()
    t.goto(0,500)
    t.setheading(0)
    counter = -1
    num_in_row = 3
    scale=50
    for elem in valid_r_theta_values:
        if isinstance(elem, str):
            continue
        ((M,T),(r,theta)) = elem
        if not M == 2:
            continue
        
        counter += 1
        angleDeg = theta * 180/math.pi
        bL = 25/abs(v_infty(r, angleDeg).real)
        
        t.goto(-600+(counter % num_in_row)*9*scale, 180-(counter // num_in_row)*7*scale)
        
        
        og_heading = t.heading()
        info_offset = 100
        t.setheading(-90)
        t.fd(info_offset)
        info = "m=" +str(M)+" t=" + str(T) + "\nangle="+str(round(angleDeg,2)) + "° r=" + str(round(r,4))
        t.pencolor("black")
        t.write(info, False, align="center")
        t.setheading(90)
        t.fd(info_offset)
        t.setheading(og_heading)
        
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=7,
            annotate=False,
            n=2,
            m=2,
            linear_c=binary_type_2_C(r, angleDeg),
            two_d_c = binary_type_2_C(r, angleDeg) * math.e ** (1j*2*math.pi/3),
            linear_num=3,
            two_d_num=2,
            isTree=False,
            isCrown=True)
    og_heading = t.heading()
    ts = t.getscreen()
    ts.getcanvas().postscript(file="binary_type-2_tiling.eps")
'''

def draw_binary_type_3_tiling_eps(t, valid_r_theta_values):
    t.up()
    t.goto(0,500)
    t.setheading(0)
    counter = -1
    num_in_row = 3
    scale=50
    for elem in valid_r_theta_values:
        if isinstance(elem, str):
            continue
        ((M,T),(r,theta)) = elem
        if not M == 3:
            continue
        
        counter += 1
        if counter >0:
            break
        angleDeg = theta * 180/math.pi
        bL = 25/abs(v_infty(r, angleDeg).real)
        
        t.goto(-600+(counter % num_in_row)*9*scale, 180-(counter // num_in_row)*7*scale)
        
        
        og_heading = t.heading()
        info_offset = 100
        t.setheading(-90)
        t.fd(info_offset)
        info = "m=" +str(M)+" t=" + str(T) + "\nangle="+str(round(angleDeg,2)) + "° r=" + str(round(r,4))
        t.pencolor("black")
        t.write(info, False, align="center")
        t.setheading(90)
        t.fd(info_offset)
        t.setheading(og_heading)
        
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=7,
            annotate=False,
            n=2,
            m=3,
            linear_c=binary_type_3_C(r, angleDeg),
            two_d_c = binary_type_3_C(r, angleDeg) * math.e ** (1j*2*math.pi/3),
            linear_num=1,
            two_d_num=1,
            isTree=True,
            isCrown=False)
        two_d_tile(
            branchLen = bL,
            t=t,
            angleDeg=angleDeg,
            r=r,
            stopK=17,
            annotate=False,
            n=2,
            m=3,
            linear_c=binary_type_3_C(r, angleDeg),
            two_d_c = binary_type_3_C(r, angleDeg) * math.e ** (1j*2*math.pi/3),
            linear_num=1,
            two_d_num=1,
            isTree=False,
            isCrown=True)
    og_heading = t.heading()
    ts = t.getscreen()
    ts.getcanvas().postscript(file="binary_type-3_tiling.eps")
    
def draw_binary_type_4_for_3D_printing(t):
    b = 200
    t.up()
    t.goto(0,-300)
    t.pensize(5)
    angleDeg = 20
    r = binary_type_4_unit_radius(angleDeg)
    print(r)
    t.down()
    nary_tree(branchLen=b,t=t,angleDeg=angleDeg,r=r,stopK=12, annotate=False, n=2)
    #crown(branchLen=b, t=t, angleDeg=angleDeg, r=r, stopK=13)
    #ts = t.getscreen()
    #ts.getcanvas().postscript(file="binary_type-4_tree.eps")

def generate_annotated_BFT(t):
    t.pensize(5)
    t.goto(0,-200)
    tree_with_stopK(branchLen=200,t=t,angleDeg=45,r=.75,stopK=3, annotate=True)
    ts = t.getscreen()
    ts.getcanvas().postscript(file="annotated_BFT.eps")
    
def main():
    #setup
    t = turtle.Turtle()
    #branchLen,t,angle,r,stopLen = 100,t,45,1/2*(-math.sqrt(2) + math.sqrt(6)),1
    turtle.tracer(0, 0)
    t.speed(0)
    t.hideturtle()
    t.setheading(90)

    sc = turtle.Screen()

    sc.setup(width=1.0,height=1.0, startx=0, starty=0)
    
    t.color("black")
    #log_spiral(branchLen=400,t=t,angleDeg=60,r=.75,stopK=5)
    #log_spiral(branchLen=200,t=t,angleDeg=30,r=.75,stopK=5)
    #t.pensize(.5) 
    scale=1.35
    #different_n_ary_trees(t)
    #different_tree_comparisons(t, scale)
    #ternary_type_m_unit(t)
    #nary_tree(135,t,45,.75,8, False, 2)
    #generate_annotated_BFT(t)
    #pruned_nary_tree(branchLen=200,t=t,angleDeg=60,r=.5,stopK=100, annotate=False, n=5)

    #probe_binary_type_m_unit_parameters()
    #the result of the above is:
    I=1j
    r_theta_values = [(('M', 'T'), ('r', 'theta')), ((1, -6), ([1.22765711495880, -1.07310941113238 - 0.356320455206363*I, -1.07310941113238 + 0.356320455206363*I, -0.247825927533573 - 0.916290471288863*I, -0.247825927533573 + 0.916290471288863*I, 0.707106781186548 - 0.707106781186548*I, 0.707106781186548 + 0.707106781186548*I], -0.7853981633974483)), ((1, -5), ([-1.27201964951407, 1.27201964951407, -0.587785252292473 - 0.809016994374947*I, -0.587785252292473 + 0.809016994374947*I, 0.587785252292473 - 0.809016994374947*I, 0.587785252292473 + 0.809016994374947*I], -0.9424777960769379)), ((1, -4), ([1.33372826484302, -1.0495475647866 - 0.532784057366145*I, -1.0495475647866 + 0.532784057366145*I, 0.38268343236509 - 0.923879532511286*I, 0.38268343236509 + 0.923879532511286*I], -1.1780972450961724)), ((1, -3), ([-1.41421356237310, 1.41421356237310, -1.0*I, 1.0*I], -1.5707963267948966)), ((1, -2), ([1.41421356237310, -0.707106781186548 - 0.707106781186549*I, -0.707106781186548 + 0.707106781186549*I], -2.356194490192345)), ((1, -1), ([-1.0*I, 1.0*I], -4.71238898038469)), ((1, 1), ([-1.0*I, 1.0*I], 4.71238898038469)), ((1, 2), ([0.707106781186547, -0.707106781186546 - 0.707106781186547*I, -0.707106781186546 + 0.707106781186547*I], 2.356194490192345)), ((1, 3), ([-0.707106781186548, 0.707106781186548, -1.0*I, 1.0*I], 1.5707963267948966)), ((1, 4), ([0.749777916806537, -0.757572390768358 - 0.384568080231984*I, -0.757572390768358 + 0.384568080231984*I, 0.38268343236509 - 0.923879532511287*I, 0.38268343236509 + 0.923879532511287*I], 1.1780972450961724)), ((1, 5), ([-0.786151377757423, 0.786151377757423, -0.587785252292474 - 0.809016994374948*I, -0.587785252292474 + 0.809016994374948*I, 0.587785252292474 - 0.809016994374948*I, 0.587785252292474 + 0.809016994374948*I], 0.9424777960769379)), ((2, -6), ([1.13115789607201, -0.965925826289068 - 0.258819045102521*I, -0.965925826289068 + 0.258819045102521*I, -0.214861133349849 - 0.770337392744067*I, -0.214861133349849 + 0.770337392744067*I, 0.615208011602914 - 0.580526487089546*I, 0.615208011602914 + 0.580526487089546*I], -0.2617993877991494)), ((2, -5), ([-1.15878094977619, 1.15878094977619, -0.503509198308116 - 0.651849796780663*I, -0.503509198308116 + 0.651849796780663*I, 0.503509198308116 - 0.651849796780663*I, 0.503509198308116 + 0.651849796780663*I], -0.3141592653589793)), ((2, -4), ([1.20083289796768, -0.923879532511287 - 0.38268343236509*I, -0.923879532511287 + 0.38268343236509*I, 0.323463083527445 - 0.729886958379422*I, 0.323463083527445 + 0.729886958379422*I], -0.39269908169872414)), ((2, -3), ([-1.27201964951407, 1.27201964951407, -0.786151377757423*I, 0.786151377757423*I], -0.5235987755982988)), ((2, -2), ([1.41421356237310, -0.707106781186548 - 0.707106781186549*I, -0.707106781186548 + 0.707106781186549*I], -0.7853981633974483)), ((2, -1), ([-1.73205080756888, 1.73205080756888], -1.5707963267948966)), ((2, 1), ([-0.577350269189626, 0.577350269189626], 1.5707963267948966)), ((2, 2), ([0.707106781186547, -0.707106781186546 - 0.707106781186547*I, -0.707106781186546 + 0.707106781186547*I], 0.7853981633974483)), ((2, 3), ([-0.786151377757423, 0.786151377757423, -1.27201964951407*I, 1.27201964951407*I], 0.5235987755982988)), ((2, 4), ([0.832755333146204, -0.923879532511287 - 0.38268343236509*I, -0.923879532511287 + 0.38268343236509*I, 0.507501865938185 - 1.1451662095779*I, 0.507501865938185 + 1.1451662095779*I], 0.39269908169872414)), ((2, 5), ([-0.862975871490758, 0.862975871490758, -0.742168599514382 - 0.960821475349482*I, -0.742168599514382 + 0.960821475349482*I, 0.742168599514382 - 0.960821475349482*I, 0.742168599514382 + 0.960821475349482*I], 0.3141592653589793)), ((3, -6), ([1.06243562660615, -0.871217616535367 - 0.120833673289935*I, -0.871217616535367 + 0.120833673289935*I, -0.182850042720857 - 0.634485202675273*I, -0.182850042720857 + 0.634485202675273*I, 0.522849845953147 - 0.461561172786661*I, 0.522849845953147 + 0.461561172786661*I], -0.08726646259971645)), ((3, -5), ([-1.07533558402498, 1.07533558402498, -0.416546397727301 - 0.501681383544138*I, -0.416546397727301 + 0.501681383544138*I, 0.416546397727301 - 0.501681383544138*I, 0.416546397727301 + 0.501681383544138*I], -0.10471975511965974)), ((3, -4), ([1.09492926357674, -0.805274243228795 - 0.169714279536309*I, -0.805274243228795 + 0.169714279536309*I, 0.257809611440426 - 0.534382364977169*I, 0.257809611440426 + 0.534382364977169*I], -0.13089969389957468)), ((3, -3), ([-1.12820632141335, 1.12820632141335, -0.522349981982434*I, 0.522349981982434*I], -0.1745329251994329)), ((3, -2), ([1.19686831751571, -0.598434158757855 - 0.27270923545906*I, -0.598434158757855 + 0.27270923545906*I], -0.26179938779914935)), ((3, -1), ([-1.41421356237310, 1.41421356237310], -0.5235987755982987)), ((3, 1), ([-0.707106781186548, 0.707106781186548], 0.5235987755982987)), ((3, 2), ([0.835513803285944, -1.38368272793204 - 0.630550668490396*I, -1.38368272793204 + 0.630550668490396*I], 0.26179938779914935)), ((3, 3), ([-0.886362698931927, 0.886362698931927, -1.91442525987036*I, 1.91442525987036*I], 0.1745329251994329)), ((3, 4), ([0.913301007896494, -1.18900111628179 - 0.250586020246396*I, -1.18900111628179 + 0.250586020246396*I, 0.732350612333544 - 1.51800101642723*I, 0.732350612333544 + 1.51800101642723*I], 0.13089969389957468)), ((3, 5), ([-0.929942256962243, 0.929942256962243, -0.979659422035065 - 1.1798851146718*I, -0.979659422035065 + 1.1798851146718*I, 0.979659422035065 - 1.1798851146718*I, 0.979659422035065 + 1.1798851146718*I], 0.10471975511965974)), ((5, -6), ([-1.04177784907464, 0.752546723047550, 0.914544313276232, -0.481684981395069 - 0.412752857042638*I, -0.481684981395069 + 0.412752857042638*I, 0.169028387770496 - 0.578131105585955*I, 0.169028387770496 + 0.578131105585955*I], 0.052359877559829876)), ((5, -5), ([-0.89851176154925 - 3.85120355388607e-101*I, -0.708756031838887 - 5.06312004929442e-102*I, -1.80308856246648e-30 + 0.556469674205471*I, 1.80308856246648e-30 - 0.556469674205471*I, 0.708756031838887 - 1.89867001848541e-102*I, 0.89851176154925 + 8.02334073726265e-103*I], 0.06283185307179585)), ((5, -4), ([-1.06326864031738, 0.645566628718154, 0.875149293705979, -0.228723641053378 - 0.457062231990563*I, -0.228723641053378 + 0.457062231990563*I], 0.07853981633974481)), ((5, -3), ([-0.838059632977374, -0.545578639220625, 0.545578639220625, 0.838059632977374], 0.10471975511965975)), ((5, -2), ([-1.12999008363554 + 0.e-23*I, 0.359223789782852 + 0.e-20*I, 0.770766293852684 + 0.e-21*I], 0.15707963267948963)), ((5, -1), ([-0.618033988749895, 0.618033988749895], 0.31415926535897926)), ((5, 1), ([-1.61803398874990, 1.61803398874990], -0.31415926535897926)), ((5, 2), ([-0.884963518248482 + 0.e-23*I, 1.29741013323441 + 0.e-20*I, 2.7837799957639 - 0.e-23*I], -0.15707963267948963)), ((5, 3), ([-1.83291633526659, -1.19323251073113, 1.19323251073113, 1.83291633526659], -0.10471975511965975)), ((5, 4), ([-0.940496091092753, 1.14266218026106, 1.54902678595022, -0.875596437559262 - 1.74971883199615*I, -0.875596437559262 + 1.74971883199615*I], -0.07853981633974481)), ((5, 5), ([-1.41092273656631 - 2.67727902983097e-102*I, -1.11295148577216 - 4.96909966932209e-103*I, -5.82282925320979e-30 - 1.79704312086333*I, 5.82282925320979e-30 + 1.79704312086333*I, 1.11295148577216 + 4.96909966932209e-103*I, 1.41092273656631 + 2.67727902983097e-102*I], -0.06283185307179585))]

    #valid_binary_type_m_unit_parameters()
    #output of above:
    
    valid_r_theta_values = [(("M","T"),("r","theta")),((1, 2), (0.707106781186547, 2.356194490192345)), ((1, 3), (0.707106781186548, 1.5707963267948966)), ((1, 4), (0.749777916806537, 1.1780972450961724)), ((1, 5), (0.786151377757423, 0.9424777960769379)), ((2, 1), (0.577350269189626, 1.5707963267948966)), ((2, 2), (0.707106781186547, 0.7853981633974483)), ((2, 3), (0.786151377757423, 0.5235987755982988)), ((2, 4), (0.832755333146204, 0.39269908169872414)), ((2, 5), (0.862975871490758, 0.3141592653589793)), ((3, 1), (0.707106781186548, 0.5235987755982987)), ((3, 2), (0.835513803285944, 0.26179938779914935)), ((3, 3), (0.886362698931927, 0.1745329251994329)), ((3, 4), (0.913301007896494, 0.13089969389957468)), ((3, 5), (0.929942256962243, 0.10471975511965974)), ((5, -6), (0.75254672304755, 0.052359877559829876)), ((5, -6), (0.914544313276232, 0.052359877559829876)), ((5, -5), (0.708756031838887, 0.06283185307179585)), ((5, -5), (0.89851176154925, 0.06283185307179585)), ((5, -4), (0.645566628718154, 0.07853981633974481)), ((5, -4), (0.875149293705979, 0.07853981633974481)), ((5, -3), (0.545578639220625, 0.10471975511965975)), ((5, -3), (0.838059632977374, 0.10471975511965975)), ((5, -2), (0.359223789782852, 0.15707963267948963)), ((5, -2), (0.770766293852684, 0.15707963267948963)), ((5, -1), (0.618033988749895, 0.31415926535897926))]

    #draw_valid_binary_type_m_unit_parameters(t, valid_r_theta_values)
    #draw_binary_type_4_units_eps(t)
    

    #draw_binary_type_4_tiling_visual_proof_eps(t, sc)
    #draw_binary_type_4_tiling_eps(t)
    #draw_binary_type_2_tiling_eps(t, valid_r_theta_values, 7)
    #draw_binary_type_3_tiling_eps(t, valid_r_theta_values)
    
    #log_crown_BFT(branchLen=50, t=t, angleDeg=60, r=.7, stopK=10, annotate=True, step=.1)
    #draw_binary_type_4_for_3D_printing(t)
    


        
            
    

main()
 
