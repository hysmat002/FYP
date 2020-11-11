# Matthew Hayes
# 31/10/2020
# A bunch of functions for the constraints and parameters of the model with 2 legs, using radau collocation.
# Does not contain the equations of motion.

from numpy import pi
import pickle as pkl

# for scaling of angles
to_degree = 180/pi
to_rad = pi/180

# Other things needed
DOFs = ['x','z','th']
signs = ['ps','ng']
contacts = ['foot1','foot2','wheel']
ground_constraints = ['contact', 'friction', 'slip_ps', 'slip_ng']

#------------------------------------------------------------------------
# get parameters

from pyomo.environ import*

def get_len(m,l):
    """
    Lengths of the leg segments and body.
    """
    if l > 1:
        return 0.2
    else:
        return 1

def get_m(m,l):
    """
    Mass of the leg segments and body.
    """
    if l == 1:
        return 6.35 + 2*m.motor_mass
    if l == 2 or l == 4:  
        return (5.0*5.0*m.len[l]*100 * 2.7)/1000 + m.motor_mass
    if l == 3 or l == 5:
        return (5.0*5.0*m.len[l]*100 * 2.7)/1000

def get_I(m,l):
    """
    Inertia of the leg segments and body.
    """
    if l == 2 or l == 4:
        return (m.m[l]-m.motor_mass)*m.len[l]**2/12 + m.motor_mass*(m.len[l]/2)**2
    if l == 3 or l ==5:
        return m.m[l]*m.len[l]**2/12
    else:
        return 0.514

def get_wheel_lengths(m,dof):
    if dof == 'x':
        return 0.4
    elif dof == 'z':
        return 0.1
    else:
        return Constraint.Skip  

#------------------------------------------------------------------------
# integration rules

def BEuler_p(m,n,l,dof):
    P = len(m.P)
    if n > 1:
        if l == 1:
            if dof == 'x':
                return m.X[n,1,'x',P] == m.X[n-1,1,'x',P] + m.hm*m.h[n] * (m.dX[n,1,'x',P]*cos(m.X[n,1,'th',P]*to_rad) + m.dX[n,1,'z',P]*sin(m.X[n,1,'th',P]*to_rad))
            elif dof == 'z':
                return m.X[n,1,'z',P] == m.X[n-1,1,'z',P] + m.hm*m.h[n] * (-m.dX[n,1,'x',P]*sin(m.X[n,1,'th',P]*to_rad) + m.dX[n,1,'z',P]*cos(m.X[n,1,'th',P]*to_rad))
            else:
                return m.X[n,1,'th',P] == m.X[n-1,1,'th',P] + m.hm*m.h[n] * m.dX[n,1,'th',P]
        else:
            return m.X[n,l,dof,P] == m.X[n-1,l,dof,P] + m.hm*m.h[n]*m.dX[n,l,dof,P]
    else:
        return Constraint.Skip

def BEuler_v(m,n,l,dof):
    P = len(m.P)
    if n > 1:
        return m.dX[n,l,dof,P] == m.dX[n-1,l,dof,P] + m.hm*m.h[n]*m.ddX[n,l,dof,P]
    else:
        return Constraint.Skip

# Radau Stuff

# Integration constraints 
radau_mat_3pt = [[0.19681547722366, -0.06553542585020, 0.02377097434822],
 [0.39442431473909, 0.29207341166523, -0.04154875212600],
 [0.37640306270047, 0.51248582618842, 0.11111111111111]]

radau_mat = [
        [0.416666125000187, -0.083333125000187],
        [0.749999625000187,  0.250000374999812],
    ]

def get_acc(m,n,l,dof,p):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    return m.acc[n,l,dof,p] == sum([radau_mat[p-1][pp-1]*m.ddX[n,l,dof,pp] for pp in range(1,P+1)])


def get_vel(m,n,l,dof,p):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    if n > 1:
        if l == 1:
            if dof == 'x':
                return m.vel[n,l,dof,p] == sum([radau_mat[p-1][pp-1]* (m.dX[n,l,'x',pp]*cos(m.X[n,l,'th',pp]*to_rad) +m.dX[n,l,'z',pp]*sin(m.X[n,l,'th',pp]*to_rad)) for pp in range(1,P+1)] )
            if dof == 'z':
                return m.vel[n,l,dof,p] == sum([radau_mat[p-1][pp-1]*(-m.dX[n,l,'x',pp]*sin(m.X[n,l,'th',pp]*to_rad) +m.dX[n,l,'z',pp]*cos(m.X[n,l,'th',pp]*to_rad)) for pp in range(1,P+1)])
            if dof == 'th':
                return m.vel[n,l,dof,p] == sum([radau_mat[p-1][pp-1]*m.dX[n,l,'th',pp] for pp in range(1,P+1)])
        else:
            return m.vel[n,l,dof,p] == sum([radau_mat[p-1][pp-1]*m.dX[n,l,dof,pp] for pp in range(1,P+1)])


def Radau_p(m,n,l,dof,p):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    if n > 1:
        return m.X[n,l,dof,p] == m.X[n-1,l,dof,P] + m.hm*m.h[n] * m.vel[n,l,dof,p]


def Radau_v(m,n,l,dof,p):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    if n > 1:
        return m.dX[n,l,dof,p] == m.dX[n-1,l,dof,P] + m.hm*m.h[n] * m.acc[n,l,dof,p]

#------------------------------------------------------------------------
# constraints

def def_tail_z(m, n, p):
    P = len(m.P)
    if n == 1 and p<P:
        return Constraint.Skip
    return m.tail_z[n,p] == m.X[n,1,'z',p] + 1*sin(m.X[n,1,'th',p]*to_rad) + 0.125*sin(m.X[n,1,'th',p]*to_rad + m.de[n]*to_rad)

def calc_V(m,n, p):
    P = len(m.P)
    if n == 1 and p<P:
        return Constraint.Skip 
    return m.V[n,p]**2 == m.dX[n,1,'x',p]**2 + m.dX[n,1,'z',p]**2


def calc_a(m,n,p):
    P = len(m.P)
    if n == 1 and p<P:
        return Constraint.Skip
    return tan(m.alpha[n,p]*to_rad)*m.dX[n,1,'x',p] == m.dX[n,1,'z',p]

def top_links(m,n,l,dof,p):
    if n == 1 and p<len(m.P):
        return Constraint.Skip
    if l > 1:
        if dof == 'x':
            return m.top[n,l,dof,p] == m.X[n,l,dof,p] - sin(m.X[n,l,'th',p]*to_rad)*0.5*m.len[l]
        if dof == 'z':
            return m.top[n,l,dof,p] == m.X[n,l,dof,p] - cos(m.X[n,l,'th',p]*to_rad)*0.5*m.len[l]
    return Constraint.Skip

def bottom_links(m,n,l,dof,p):
    if n == 1 and p<len(m.P):
        return Constraint.Skip
    if l > 1:
        if dof == 'x':
            return m.bottom[n,l,'x',p] == m.X[n,l,'x',p] + sin(m.X[n,l,'th',p]*to_rad)*0.5*m.len[l]
        if dof == 'z':
            return m.bottom[n,l,'z',p] == m.X[n,l,'z',p] + cos(m.X[n,l,'th',p]*to_rad)*0.5*m.len[l]
    return Constraint.Skip

def connect(m,n,l,dof,p):
    if (n == 1 and p<len(m.P)) or dof == 'th':
        return Constraint.Skip
    else:
        if l == 1:
            return Constraint.Skip
        if l == 2 or l == 4:
            return m.top[n,l,dof,p] - m.X[n,1,dof,p] == 0 
        if l == 3 or l == 5:
            return m.top[n,l,dof,p] - m.bottom[n,l-1,dof,p] ==  0 

def def_th_joint(m,n,l,p):
    if n == 1 and p<len(m.P):
        return Constraint.Skip
    if l == 4: # 4 is connected to 1 not 3
        return m.th_joint[n,l,p] == m.X[n,l,'th',p] - m.X[n,1,'th',p]
    elif l > 1:
        return m.th_joint[n,l,p] == m.X[n,l,'th',p] - m.X[n,l-1,'th',p]
    else:
        return Constraint.Skip

def def_dth_joint(m,n,l,p):
    P = len(m.P)
    if n == 1 and p<P:
        return Constraint.Skip
    if l == 4:
        return m.dth_joint[n,l,P] == m.dX[n,l,'th',P] - m.dX[n,1,'th',P]
    elif l > 1:
        return m.dth_joint[n,l,P] == m.dX[n,l,'th',P] - m.dX[n,l-1,'th',P]
    else:
        return Constraint.Skip

def min_torque(m,n,l):
    mbody = sum( m.m[l] for l in range(1,len(m.L)+1) ) 
    BW = mbody*m.g
    P = len(m.P)
    if l == 1: 
        return Constraint.Skip
    return m.Tc[n,l]*BW  >= - m.tauMax - m.tauMax/m.omegaMax * m.dth_joint[n,l,P]

def max_torque(m,n,l):
    mbody = sum( m.m[l] for l in range(1,len(m.L)+1) ) 
    BW = mbody*m.g
    P = len(m.P)
    if l == 1: 
        return Constraint.Skip
    return m.Tc[n,l]*BW  <= m.tauMax - m.tauMax/m.omegaMax * m.dth_joint[n,l,P]
        

def def_contact_p(m,n,contact,p):
    if n == 1 and p<len(m.P):
        return Constraint.Skip
    if contact == 'foot1':
        return m.contact_p[n,contact,p] == - m.bottom[n,3,'z',p]
    elif contact == 'foot2':
        return m.contact_p[n,contact,p] == - m.bottom[n,5,'z',p]
    else:
        return m.contact_p[n,contact,p] == -(m.X[n,1,'z',p] + m.wheel_leg['z']*cos(m.X[n,1,'th',p]*to_rad) - m.wheel_leg['x']*sin(m.X[n,1,'th',p]*to_rad))

def def_contact_v(m,n,contact,p):
    if n == 1 and p<len(m.P):
        return Constraint.Skip
    if contact == 'foot1':
        return m.contact_v[n,contact, 'ps',p] -  m.contact_v[n,contact, 'ng',p] == m.dX[n,3,'x',p] + cos(m.X[n,3,'th',p]*to_rad)*m.dX[n,3,'th',p]*to_rad*0.5*m.len[3]
    elif contact == 'foot2':
        return m.contact_v[n,contact, 'ps',p] -  m.contact_v[n,contact, 'ng',p] == m.dX[n,5,'x',p] + cos(m.X[n,5,'th',p]*to_rad)*m.dX[n,5,'th',p]*to_rad*0.5*m.len[5]
    else:
        return m.contact_v[n,contact, 'ps',p] -  m.contact_v[n,contact, 'ng',p] == ( m.dX[n,1,'x',p]*cos(m.X[n,3,'th',p]*to_rad) - m.dX[n,1,'z',p]*sin(m.X[n,3,'th',p]*to_rad) ) - m.wheel_leg['x']*sin(m.X[n,1,'th',p]*to_rad)*m.dX[n,1,'th',p]*to_rad + m.wheel_leg['z']*cos(m.X[n,1,'th',p]*to_rad)*m.dX[n,1,'th',p]*to_rad

def def_friction_cone(m,n, contact,p):
    P = len(m.P)
    if n == 1 and p<P:
        return Constraint.Skip
    return m.friction_cone[n, contact,p]  == m.mu*m.GRF[n,contact, 'z', 'ng',p] - (m.GRF[n,contact, 'x', 'ps',p] + m.GRF[n,contact, 'x', 'ng',p])


def def_ground_varA(m,n, contact, gc):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    if gc == 'contact':
        if n > N:
            return m.ground_varA[n,contact, gc] == sum(m.contact_p[n+1,contact,p] for p in range(1,P+1))
        else:
            return m.ground_varA[n,contact, gc] == sum(m.contact_p[n,contact,p] for p in range(1,P+1))
    if gc == 'friction':
        return m.ground_varA[n,contact, gc] == sum((m.contact_v[n,contact,'ps',p]+m.contact_v[n,contact,'ng',p]) for p in range(1,P+1))
    if gc == 'slip_ps':
        return m.ground_varA[n,contact, gc] == sum(m.contact_v[n,contact,'ps',p] for p in range(1,P+1))
    if gc == 'slip_ng':
        return m.ground_varA[n,contact, gc] == sum(m.contact_v[n,contact,'ng',p] for p in range(1,P+1))



def def_ground_varB(m,n, contact, gc):
    mbody = sum( m.m[l] for l in range(1,len(m.L)+1) ) 
    BW = mbody*m.g

    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    if gc == 'contact':
        return m.ground_varB[n,contact, gc] == sum(BW*m.GRF[n,contact, 'z', 'ng',p] for p in range(1,P+1))
    if gc == 'friction':
        return m.ground_varB[n,contact, gc] == sum((BW*m.friction_cone[n, contact,'ps',p]+BW*m.friction_cone[n, contact,'ng',p]) for p in range(1,P+1))
    if gc == 'slip_ps':
        return m.ground_varB[n,contact, gc] == sum(BW*m.GRF[n,contact,'x','ps',p] for p in range(1,P+1))
    if gc == 'slip_ng':
        return m.ground_varB[n,contact, gc] == sum(BW*m.GRF[n,contact,'x','ng',p] for p in range(1,P+1))


def ground_comp(m,n,contact,gc):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    return m.ground_penalties[n,contact, gc] == m.ground_varA[n,contact,gc]*m.ground_varB[n,contact,gc]

def contact(m,n,contact):
    P = len(m.P)
    if n == 1:
       return Constraint.Skip
    if n < len(m.N):
        A = sum(m.contact_p[n+1,contact,p] for p in range(1,P+1))
        B = sum(m.GRF[n,contact, 'z', 'ng',p] for p in range(1,P+1))
        return A*B <= m.ground_penalties[n,contact,'contact']
    else:
        return Constraint.Skip

def contact_moving(m,n,contact):
    P = len(m.P)
    if n == 1:
       return Constraint.Skip
    if n < len(m.N):
        A = sum(m.phi[n+1,contact,p] for p in range(1,P+1))
        B = sum(m.GRF[n,contact, 'z', 'ng',p] for p in range(1,P+1))
        return A*B <= m.ground_penalties[n,contact,'contact']
    else:
        return Constraint.Skip

def friction(m,n,contact):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    A = sum((m.contact_v[n,contact,'ps',p]+m.contact_v[n,contact,'ng',p]) for p in range(1,P+1))
    B = sum((m.friction_cone[n, contact,p]) for p in range(1,P+1))
    return A*B <= m.ground_penalties[n,contact,'friction']

def slip_ps(m,n,contact):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    A = sum(m.contact_v[n,contact,'ps',p] for p in range(1,P+1))
    B = sum(m.GRF[n,contact,'x','ps',p] for p in range(1,P+1))
    return A*B <= m.ground_penalties[n,contact,'slip_ps']

def slip_ng(m,n,contact):
    P = len(m.P)
    if n == 1:
        return Constraint.Skip
    A = sum(m.contact_v[n,contact,'ng',p] for p in range(1,P+1))
    B = sum(m.GRF[n,contact,'x','ng',p] for p in range(1,P+1))
    return A*B <= m.ground_penalties[n,contact,'slip_ng']

#------------------------------------------------------------------------
# Objective functions

def minPenalty(m):
    P = len(m.P)
    penalty_sum = 0
    for n in range(1,len(m.N)+1):
        for contact in contacts:
            for gc in ground_constraints:
                penalty_sum += m.ground_penalties[n,contact, gc]
    return (10*penalty_sum)

def minTime(m):
    P = len(m.P)
    T = sum(m.h[n] for n in range(1,len(m.N)+1))
    penalty_sum = 0
    for n in range(1,len(m.N)+1):
        for contact in contacts:
            for gc in ground_constraints:
                penalty_sum += m.ground_penalties[n,contact, gc]
    return (T + 1000*penalty_sum )

def LQR(m):
    control_sum = 0
    error = 0
    penalty_sum = 0
    for n in range(1,len(m.N)+1):
        for contact in contacts:
            for gc in ground_constraints:
                penalty_sum += m.ground_penalties[n,contact, gc]

        #control_sum += m.T[n]**2 + m.de[n]**2

        control_sum +=  m.T[n]**2 + m.de[n]**2 + m.df[n]**2 + m.Tc[n,2]**2+ m.Tc[n,3]**2 + m.Tc[n,4]**2+ m.Tc[n,5]**2
        
    error += (m.dX[n,2,'th'])**2 + (m.dX[n,3,'th'])**2 + (m.dX[n,4,'th'])**2 + (m.dX[n,5,'th'])**2
    
    return (control_sum + error + 10000*penalty_sum) 

def minDistance(m):
    P = len(m.P)
    penalty_sum = 0
    for n in range(1,len(m.N)+1):
        for contact in contacts:
            for gc in ground_constraints:
                penalty_sum += m.ground_penalties[n,contact, gc]
    return (penalty_sum*10000 + ((m.X[len(m.N),1,'x',P])**2))

#------------------------------------------------------------------------
# setting things up

def init_opt_old(tol = 1e-6):
    print("2.12.12")
    opt = SolverFactory('ipopt',executable = '/Users/matthayes/anaconda3/envs/EEE4022S/bin/ipopt')
    opt.options["print_level"] = 5 
    opt.options["max_iter"] = 6000000 
    opt.options["max_cpu_time"] = 100000 
    opt.options["Tol"] = tol
    opt.options["halt_on_ampl_error"] = "yes"
    opt.options["OF_acceptable_obj_change_tol"] = 1e-4
    opt.options["OF_ma86_scaling"] = 'none'

    return opt

def init_opt(tol = 1e-6):
    opt = SolverFactory('ipopt',executable = '/Users/matthayes/anaconda3/envs/Optimization/bin/ipopt')
    opt.options["print_level"] = 5 
    opt.options["max_iter"] = 6000000 
    opt.options["max_cpu_time"] = 100000 
    opt.options["Tol"] = tol
    opt.options["halt_on_ampl_error"] = "yes"
    opt.options["OF_acceptable_obj_change_tol"] = 1e-4
    opt.options["OF_ma86_scaling"] = 'none'

    return opt

#------------------------------------------------------------------------
# animation
import numpy as np
def plot_UAV_legs(i,m,ax1, ax2, xmin, xmax, ymin, ymax):
    '''
    Plots a single frame for the animation of Aircraft with two leg Pyomo model.
    i    - the current node (int),
    m    - the Pyomo model
    ax1  - the axis for the "big picture" (Pyplot axis)
    ax2  - a zoomed in on the "robot" (Pyplot axis)
    xmin - lower bound on the x axis (int)
    xmax - upper bound on the x axis (int)
    ymin - lower bound on the y axis (int)
    ymax - upper bound on the y axis (int)
    '''
    P = len(m.P)
    x = list(m.X[:,1,'x',P].value)
    z = list(m.X[:,1,'z',P].value)
    
    L1 = 0.825
    L2 = 0.45
    Le = 0.125

    # body positions
    front_x = m.X[i,1,'x',P].value + L2* np.cos(m.X[i,1,'th',P].value*to_rad)
    front_z = m.X[i,1,'z',P].value - L2* np.sin(m.X[i,1,'th',P].value*to_rad)
    
    back_x = m.X[i,1,'x',P].value - L1* np.cos(m.X[i,1,'th',P].value*to_rad)
    back_z = m.X[i,1,'z',P].value + L1* np.sin(m.X[i,1,'th',P].value*to_rad)
    
    # elevator positions
    e_x = back_x - Le* np.cos(m.de[i].value*to_rad+ m.X[i,1,'th',P].value*to_rad) # m.de is divided by 1000 because it is scaled by 1000 in the model 
    e_z = back_z + Le* np.sin(m.de[i].value*to_rad+ m.X[i,1,'th',P].value*to_rad)
    
    # wing positions
    w_front_x = m.X[i,1,'x',P].value + 0.19 * np.cos(m.X[i,1,'th',P].value*to_rad)
    w_front_z = m.X[i,1,'z',P].value - 0.19 * np.sin(m.X[i,1,'th',P].value*to_rad)
    
    w_back_x = m.X[i,1,'x',P].value - 0.1 * np.cos(m.X[i,1,'th',P].value*to_rad)
    w_back_z = m.X[i,1,'z',P].value + 0.1 * np.sin(m.X[i,1,'th',P].value*to_rad)
    
    # flap
    
    flap_x = w_back_x - 0.09 * np.cos(m.df[i].value*to_rad+ m.X[i,1,'th',P].value*to_rad)
    flap_z = w_back_z + 0.09 * np.sin(m.df[i].value*to_rad+ m.X[i,1,'th',P].value*to_rad)
    
    # wheel position
    wheel_x = m.X[i,1,'x',P].value + m.wheel_leg['x']*np.cos(m.X[i,1,'th',P].value*to_rad) + m.wheel_leg['z']*np.sin(m.X[i,1,'th',P].value*to_rad)
    wheel_z = -m.contact_p[i,'wheel',P].value
    
    # wheel leg
    wheel_top_x = m.X[i,1,'x',P].value + m.wheel_leg['x']*np.cos(m.X[i,1,'th',P].value*to_rad)
    wheel_top_z = m.X[i,1,'z',P].value - m.wheel_leg['x']*np.sin(m.X[i,1,'th',P].value*to_rad)
    
    # ------------------------
    ax1.clear(); ax1.grid()
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    
    # ground
    ax1.plot([-100,200],[0,0], 'C2')
    
    # body
    ax1.plot([front_x,back_x],[front_z,back_z],'grey')
    
    # elevator
    ax1.plot([back_x,e_x],[back_z,e_z],'C0', linewidth = 5)
    ax1.plot(back_x,back_z,'.',color="Grey")
    
    # wing
    ax1.plot([flap_x,w_back_x,w_front_x],[flap_z,w_back_z,w_front_z],'C0', linewidth = 5)
    
    # com
    ax1.plot(m.X[i,1,'x',P].value,m.X[i,1,'z',P].value,'k.')
    
    # legs
    ax1.plot([m.top[i,2,'x',P].value,m.bottom[i,2,'x',P].value],[m.top[i,2,'z',P].value,m.bottom[i,2,'z',P].value], 'k')
    ax1.plot([m.top[i,3,'x',P].value,m.bottom[i,3,'x',P].value],[m.top[i,3,'z',P].value,m.bottom[i,3,'z',P].value], 'k')
    
    ax1.plot([m.top[i,4,'x',P].value,m.bottom[i,4,'x',P].value],[m.top[i,4,'z',P].value,m.bottom[i,4,'z',P].value], color = 'grey')
    ax1.plot([m.top[i,5,'x',P].value,m.bottom[i,5,'x',P].value],[m.top[i,5,'z',P].value,m.bottom[i,5,'z',P].value], color = 'grey')
    
    ax1.plot(m.top[i,3,'x',P].value,m.top[i,3,'z',P].value, '.k')
    ax1.plot(m.top[i,5,'x',P].value,m.top[i,5,'z',P].value, color = 'grey')
    
    # wheel
    ax1.plot([wheel_top_x, wheel_x],[wheel_top_z, wheel_z], color = 'grey')
    ax1.plot(wheel_x, wheel_z,'.k')
    
    # path
    ax1.plot(x[0:i],z[0:i],':', color = "Grey")
    
    # ------------------------
    
    ax2.clear()
    ax2.set_xlim([m.X[i,1,'x',P].value-1,m.X[i,1,'x',P].value+1])
    ax2.set_ylim([m.X[i,1,'z',P].value+1,m.X[i,1,'z',P].value-1])
    
    # ground
    ax2.plot([-100,200],[0,0], 'C2')
    
    # body
    ax2.plot([front_x,back_x],[front_z,back_z],'grey')
    ax2.plot(m.X[i,1,'x',P].value,m.X[i,1,'z',P].value,'k.')
    
    # elevator
    ax2.plot([back_x,e_x],[back_z,e_z],'C0', linewidth = 5)
    ax2.plot(back_x,back_z,'.',color="Grey")
    
    # wing
    ax2.plot([flap_x,w_back_x,w_front_x],[flap_z,w_back_z,w_front_z],'C0', linewidth = 5)
    
    # legs
    ax2.plot([m.top[i,4,'x',P].value,m.bottom[i,4,'x',P].value],[m.top[i,4,'z',P].value,m.bottom[i,4,'z',P].value], color = 'grey',linewidth = 2)
    ax2.plot([m.top[i,5,'x',P].value,m.bottom[i,5,'x',P].value],[m.top[i,5,'z',P].value,m.bottom[i,5,'z',P].value], color = 'grey',linewidth = 2)

    ax2.plot(m.top[i,4,'x',P].value,m.top[i,4,'z',P].value, '.', color = 'grey')
    ax2.plot(m.top[i,5,'x',P].value,m.top[i,5,'z',P].value, '.', color = 'grey')
    
    
    ax2.plot([m.top[i,2,'x',P].value, m.bottom[i,2,'x',P].value],[m.top[i,2,'z',P].value,m.bottom[i,2,'z',P].value],'k',linewidth = 2)
    ax2.plot([m.top[i,3,'x',P].value, m.bottom[i,3,'x',P].value],[m.top[i,3,'z',P].value,m.bottom[i,3,'z',P].value],'k',linewidth = 2)
    
    ax2.plot(m.top[i,2,'x',P].value,m.top[i,2,'z',P].value, '.k')
    ax2.plot(m.top[i,3,'x',P].value,m.top[i,3,'z',P].value, '.k')
    
    
    # com
    ax2.plot(m.X[i,1,'x',P].value,m.X[i,1,'z',P].value,'ko')
    
    # wheel
    ax2.plot([wheel_top_x, wheel_x],[wheel_top_z, wheel_z], color = 'grey')
    if wheel_z >= -1e-7:
        ax2.plot(wheel_x, wheel_z,'ro')
    else:
        ax2.plot(wheel_x, wheel_z,'ko')
        
    if m.bottom[i,3,'z',P].value >= -1e-7:
        ax2.plot(m.bottom[i,3,'x',P].value, m.bottom[i,3,'z',P].value,'.r')
    else:
        ax2.plot(m.bottom[i,3,'x',P].value, m.bottom[i,3,'z',P].value,'.k')
        
    if m.bottom[i,5,'z',P].value >= -1e-7:
        ax2.plot(m.bottom[i,5,'x',P].value, m.bottom[i,5,'z',P].value,'.r')
    else:
        ax2.plot(m.bottom[i,5,'x',P].value, m.bottom[i,5,'z',P].value,'.', color = 'grey')
    
    ax2.plot(x[0:i],z[0:i],':', color = "Grey")
    
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.grid()


def save_data_radau(m,N,hm,h_var,name, legs = True):

    if legs:
        data = {
            'P' : len(m.P),
            'Cost': float(m.Cost.expr()),
            'Cost_type':"Min Penalty",
            'l1':m.len[1],
            'l2':m.len[2],
            'l3':m.len[3],
            'm1':m.m[1],
            'm2':m.m[2],
            'm3':m.m[3],
            'x1': list(m.X[:,1,'x',:].value),
            'z1': list(m.X[:,1,'z',:].value),
            'th1': list(m.X[:,1,'th',:].value),
            'dx1': list(m.dX[:,1,'x',:].value),
            'dz1': list(m.dX[:,1,'z',:].value),
            'dth1': list(m.dX[:,1,'th',:].value),
            'ddx1': list(m.ddX[:,1,'x',:].value),
            'ddz1': list(m.ddX[:,1,'z',:].value),
            'ddth1': list(m.ddX[:,1,'th',:].value),
            'x2': list(m.X[:,2,'x',:].value),
            'z2': list(m.X[:,2,'z',:].value),
            'th2': list(m.X[:,2,'th',:].value),
            'dx2': list(m.dX[:,2,'x',:].value),
            'dz2': list(m.dX[:,2,'z',:].value),
            'dth2': list(m.dX[:,2,'th',:].value),
            'ddx2': list(m.ddX[:,2,'x',:].value),
            'ddz2': list(m.ddX[:,2,'z',:].value),
            'ddth2': list(m.ddX[:,2,'th',:].value),
            'x3': list(m.X[:,3,'x',:].value),
            'z3': list(m.X[:,3,'z',:].value),
            'th3': list(m.X[:,3,'th',:].value),
            'dx3': list(m.dX[:,3,'x',:].value),
            'dz3': list(m.dX[:,3,'z',:].value),
            'dth3': list(m.dX[:,3,'th',:].value),
            'ddx3': list(m.ddX[:,3,'x',:].value),
            'ddz3': list(m.ddX[:,3,'z',:].value),
            'ddth3': list(m.ddX[:,3,'th',:].value),
            'x4': list(m.X[:,4,'x',:].value),
            'z4': list(m.X[:,4,'z',:].value),
            'th4': list(m.X[:,4,'th',:].value),
            'dx4': list(m.dX[:,4,'x',:].value),
            'dz4': list(m.dX[:,4,'z',:].value),
            'dth4': list(m.dX[:,4,'th',:].value),
            'ddx4': list(m.ddX[:,4,'x',:].value),
            'ddz4': list(m.ddX[:,4,'z',:].value),
            'ddth4': list(m.ddX[:,4,'th',:].value),
            'x5': list(m.X[:,5,'x',:].value),
            'z5': list(m.X[:,5,'z',:].value),
            'th5': list(m.X[:,5,'th',:].value),
            'dx5': list(m.dX[:,5,'x',:].value),
            'dz5': list(m.dX[:,5,'z',:].value),
            'dth5': list(m.dX[:,5,'th',:].value),
            'ddx5': list(m.ddX[:,5,'x',:].value),
            'ddz5': list(m.ddX[:,5,'z',:].value),
            'ddth5': list(m.ddX[:,5,'th',:].value),
            'V': list(m.V[:,:].value),
            'a': list(m.alpha[:,:].value),
            'T': list(m.T[:].value),
            'de': list(m.de[:].value),
            'df': list(m.df[:].value),
            #'th_ub':45*np.pi/180,
            #'th_lb':-45*np.pi/180,
            'tau1':list(m.Tc[:,2].value),
            'tau2':list(m.Tc[:,3].value),
            'tau3':list(m.Tc[:,4].value),
            'tau4':list(m.Tc[:,5].value),
            #'alpha_up':19*np.pi/180,
            #'alpha_lb':-19*np.pi/180,
            'N':N,
            'hm':hm,
            'h_var': h_var,
            'wheel penaltys': {
                'contact': list(m.ground_penalties[:,'wheel','contact'].value),
                'friction': list(m.ground_penalties[:,'wheel','friction'].value),
                'slip_ps': list(m.ground_penalties[:,'wheel','slip_ps'].value),
                'slip_ng': list(m.ground_penalties[:,'wheel','slip_ng'].value)
                },
            'foot1 penaltys': {
                'contact': list(m.ground_penalties[:,'foot1','contact'].value),
                'friction': list(m.ground_penalties[:,'foot1','friction'].value),
                'slip_ps': list(m.ground_penalties[:,'foot1','slip_ps'].value),
                'slip_ng': list(m.ground_penalties[:,'foot1','slip_ng'].value)
                },
            'foot2 penaltys': {
                'contact': list(m.ground_penalties[:,'foot2','contact'].value),
                'friction': list(m.ground_penalties[:,'foot2','friction'].value),
                'slip_ps': list(m.ground_penalties[:,'foot2','slip_ps'].value),
                'slip_ng': list(m.ground_penalties[:,'foot2','slip_ng'].value)
                },
            'wheel forces':{
                'normal': list(m.GRF[:,'wheel','z','ng',:].value),
                'x_ps': list(m.GRF[:,'wheel','x','ps',:].value),
                'x_ng': list(m.GRF[:,'wheel','x','ng',:].value)
                },
            'foot1 forces':{
                'normal': list(m.GRF[:,'foot1','z','ng',:].value),
                'x_ps': list(m.GRF[:,'foot1','x','ps',:].value),
                'x_ng': list(m.GRF[:,'foot1','x','ng',:].value)
                },
            'foot2 forces':{
                'normal': list(m.GRF[:,'foot2','z','ng',:].value),
                'x_ps': list(m.GRF[:,'foot2','x','ps',:].value),
                'x_ng': list(m.GRF[:,'foot2','x','ng',:].value)
                    }
            }
    else:
        data = {
            'P': m.P,
            'Cost': float(m.Cost.expr()),
            'x': list(m.X[:,1].value), 
            'z': list(m.X[:,2].value),
            'th': list(m.X[:,3].value),
            'dx': list(m.dX[:,1].value),
            'dz': list(m.dX[:,2].value),
            'dth': list(m.dX[:,3].value),
            'ddx': list(m.ddX[:,1].value),
            'ddz': list(m.ddX[:,2].value),
            'ddth': list(m.ddX[:,3].value),
            'T': list(m.T[:].value),
            'de': list(m.de[:].value),
            'df': list(m.de[:].value),
            'alpha': list(m.alpha[:].value),
            'V': list(m.V[:].value),
            'N':N,
            'hm':hm,
            'h_var': 0.2,
            'Front penaltys': {
                'normal': list(m.ground_penalty[:,'front','normal'].value),
                'friction': list(m.ground_penalty[:,'front','friction'].value),
                'slip_ps': list(m.ground_penalty[:,'front','slip_ps'].value),
                'slip_ng': list(m.ground_penalty[:,'front','slip_ng'].value)
            },
            'Back penaltys': {
                'normal': list(m.ground_penalty[:,'back','normal'].value),
                'friction': list(m.ground_penalty[:,'back','friction'].value),
                'slip_ps': list(m.ground_penalty[:,'back','slip_ps'].value),
                'slip_ng': list(m.ground_penalty[:,'back','slip_ng'].value)
            },
            'Front forces':{
                'normal': list(m.GRF[:,'front',2,'ng'].value),
                'x_ps': list(m.GRF[:,'front',1,'ps'].value),
                'x_ng': list(m.GRF[:,'front',1,'ng'].value)
            },
            'Back forces':{
                'normal': list(m.GRF[:,'back',2,'ng'].value),
                'x_ps': list(m.GRF[:,'back',1,'ps'].value),
                'x_ng': list(m.GRF[:,'back',1,'ng'].value)
            }
        }
            
    outfile = open(name,'wb')
    pkl.dump(data,outfile)
    outfile.close()
    return("Data Saved")