import pickle as pkl
import numpy as np

def save_data(m, N, hm, h_var, name, legs):
    '''
    Saves most of the data for Aircraft with one leg Pyomo model.
    m - Pyomo model
    N - number of nodes
    hm - master time step
    h_var - howmuch the time step is allowed to vary
    name - name of the pickle file created
    legs - 0, 1 or 2 (int) number of legs
    '''
    if legs > 3 or legs < 0:
        #todo : possibly throw error
        return "Nop"

    if legs == 0:
        data = {
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
    
    if legs == 1:
        data = {
            'Cost': float(m.Cost.expr()),
            'l1':m.len[1],
            'l2':m.len[2],
            'l3':m.len[3],
            'm1':m.m[1],
            'm2':m.m[2],
            'm3':m.m[3],
            'x1': list(m.X[:,1,'x'].value),
            'z1': list(m.X[:,1,'z'].value),
            'th1': list(m.X[:,1,'th'].value),
            'dx1': list(m.dX[:,1,'x'].value),
            'dz1': list(m.dX[:,1,'z'].value),
            'dth1': list(m.dX[:,1,'th'].value),
            'ddx1': list(m.ddX[:,1,'x'].value),
            'ddz1': list(m.ddX[:,1,'z'].value),
            'ddth1': list(m.ddX[:,1,'th'].value),
            'x2': list(m.X[:,2,'x'].value),
            'z2': list(m.X[:,2,'z'].value),
            'th2': list(m.X[:,2,'th'].value),
            'dx2': list(m.dX[:,2,'x'].value),
            'dz2': list(m.dX[:,2,'z'].value),
            'dth2': list(m.dX[:,2,'th'].value),
            'ddx2': list(m.ddX[:,2,'x'].value),
            'ddz2': list(m.ddX[:,2,'z'].value),
            'ddth2': list(m.ddX[:,2,'th'].value),
            'x3': list(m.X[:,3,'x'].value),
            'z3': list(m.X[:,3,'z'].value),
            'th3': list(m.X[:,3,'th'].value),
            'dx3': list(m.dX[:,3,'x'].value),
            'dz3': list(m.dX[:,3,'z'].value),
            'dth3': list(m.dX[:,3,'th'].value),
            'ddx3': list(m.ddX[:,3,'x'].value),
            'ddz3': list(m.ddX[:,3,'z'].value),
            'ddth3': list(m.ddX[:,3,'th'].value),
            'V': list(m.V[:].value),
            'a': list(m.alpha[:].value),
            'T': list(m.T[:].value),
            'de': list(m.de[:].value),
            'df': list(m.df[:].value),
            #'th_ub':45*np.pi/180,
            #'th_lb':-45*np.pi/180,
            'tau1':list(m.Tc[:,2].value),
            'tau2':list(m.Tc[:,3].value),
            #'alpha_up':19*np.pi/180,
            #'alpha_lb':-19*np.pi/180,
            'N':N,
            'hm':hm,
            'h_var': h_var,
            'wheel penaltys': {
                'normal': list(m.ground_penalties[:,'wheel','contact'].value),
                'friction': list(m.ground_penalties[:,'wheel','friction'].value),
                'slip_ps': list(m.ground_penalties[:,'wheel','slip_ps'].value),
                'slip_ng': list(m.ground_penalties[:,'wheel','slip_ng'].value)
            },
            'foot penaltys': {
                'normal': list(m.ground_penalties[:,'foot','contact'].value),
                'friction': list(m.ground_penalties[:,'foot','friction'].value),
                'slip_ps': list(m.ground_penalties[:,'foot','slip_ps'].value),
                'slip_ng': list(m.ground_penalties[:,'foot','slip_ng'].value)
            },
            'wheel forces':{
                'normal': list(m.GRF[:,'wheel','z','ng'].value),
                'x_ps': list(m.GRF[:,'wheel','x','ps'].value),
                'x_ng': list(m.GRF[:,'wheel','x','ng'].value)
            },
            'foot forces':{
                'normal': list(m.GRF[:,'foot','z','ng'].value),
                'x_ps': list(m.GRF[:,'foot','x','ps'].value),
                'x_ng': list(m.GRF[:,'foot','x','ng'].value)
            }
        }
    
    if legs == 2:
        data = {
        'Cost': float(m.Cost.expr()),
        'Cost_type':"Min Penalty",
        'l1':m.len[1],
        'l2':m.len[2],
        'l3':m.len[3],
        'm1':m.m[1],
        'm2':m.m[2],
        'm3':m.m[3],
        'x1': list(m.X[:,1,'x'].value),
        'z1': list(m.X[:,1,'z'].value),
        'th1': list(m.X[:,1,'th'].value),
        'dx1': list(m.dX[:,1,'x'].value),
        'dz1': list(m.dX[:,1,'z'].value),
        'dth1': list(m.dX[:,1,'th'].value),
        'ddx1': list(m.ddX[:,1,'x'].value),
        'ddz1': list(m.ddX[:,1,'z'].value),
        'ddth1': list(m.ddX[:,1,'th'].value),
        'x2': list(m.X[:,2,'x'].value),
        'z2': list(m.X[:,2,'z'].value),
        'th2': list(m.X[:,2,'th'].value),
        'dx2': list(m.dX[:,2,'x'].value),
        'dz2': list(m.dX[:,2,'z'].value),
        'dth2': list(m.dX[:,2,'th'].value),
        'ddx2': list(m.ddX[:,2,'x'].value),
        'ddz2': list(m.ddX[:,2,'z'].value),
        'ddth2': list(m.ddX[:,2,'th'].value),
        'x3': list(m.X[:,3,'x'].value),
        'z3': list(m.X[:,3,'z'].value),
        'th3': list(m.X[:,3,'th'].value),
        'dx3': list(m.dX[:,3,'x'].value),
        'dz3': list(m.dX[:,3,'z'].value),
        'dth3': list(m.dX[:,3,'th'].value),
        'ddx3': list(m.ddX[:,3,'x'].value),
        'ddz3': list(m.ddX[:,3,'z'].value),
        'ddth3': list(m.ddX[:,3,'th'].value),
        'x4': list(m.X[:,4,'x'].value),
        'z4': list(m.X[:,4,'z'].value),
        'th4': list(m.X[:,4,'th'].value),
        'dx4': list(m.dX[:,4,'x'].value),
        'dz4': list(m.dX[:,4,'z'].value),
        'dth4': list(m.dX[:,4,'th'].value),
        'ddx4': list(m.ddX[:,4,'x'].value),
        'ddz4': list(m.ddX[:,4,'z'].value),
        'ddth4': list(m.ddX[:,4,'th'].value),
        'x5': list(m.X[:,5,'x'].value),
        'z5': list(m.X[:,5,'z'].value),
        'th5': list(m.X[:,5,'th'].value),
        'dx5': list(m.dX[:,5,'x'].value),
        'dz5': list(m.dX[:,5,'z'].value),
        'dth5': list(m.dX[:,5,'th'].value),
        'ddx5': list(m.ddX[:,5,'x'].value),
        'ddz5': list(m.ddX[:,5,'z'].value),
        'ddth5': list(m.ddX[:,5,'th'].value),
        'V': list(m.V[:].value),
        'a': list(m.alpha[:].value),
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
            'normal': list(m.GRF[:,'wheel','z','ng'].value),
            'x_ps': list(m.GRF[:,'wheel','x','ps'].value),
            'x_ng': list(m.GRF[:,'wheel','x','ng'].value)
            },
        'foot1 forces':{
            'normal': list(m.GRF[:,'foot1','z','ng'].value),
            'x_ps': list(m.GRF[:,'foot1','x','ps'].value),
            'x_ng': list(m.GRF[:,'foot1','x','ng'].value)
            },
        'foot2 forces':{
            'normal': list(m.GRF[:,'foot2','z','ng'].value),
            'x_ps': list(m.GRF[:,'foot2','x','ps'].value),
            'x_ng': list(m.GRF[:,'foot2','x','ng'].value)
                }
        }
            
    outfile = open(name,'wb')
    pkl.dump(data,outfile)
    outfile.close()
    return("Data Saved")

def plot_UAV_leg(i,m,ax1, ax2, xmin, xmax, ymin, ymax, scale = 100):
    '''
    Plots a single frame for the animation of Aircraft with one leg Pyomo model.
    i    - the current node (int),
    m    - the Pyomo model
    ax1  - the axis for the "big picture" (Pyplot axis)
    ax2  - a zoomed in on the "robot" (Pyplot axis)
    xmin - lower bound on the x axis (int)
    xmax - upper bound on the x axis (int)
    ymin - lower bound on the y axis (int)
    ymax - upper bound on the y axis (int)
    scale - value by which delta_e and delta_f are scaled by. (int) default = 100
    '''
    x = list(m.X[:,1,'x'].value)
    z = list(m.X[:,1,'z'].value)
    
    L1 = 0.625
    L2 = 0.5
    Le = 0.125

    # body positions
    front_x = m.X[i,1,'x'].value + L2* np.cos(m.X[i,1,'th'].value)
    front_z = m.X[i,1,'z'].value - L2* np.sin(m.X[i,1,'th'].value)
    
    back_x = m.X[i,1,'x'].value - L1* np.cos(m.X[i,1,'th'].value)
    back_z = m.X[i,1,'z'].value + L1* np.sin(m.X[i,1,'th'].value)
    
    # elevator positions
    e_x = back_x - Le* np.cos(m.de[i].value/scale+ m.X[i,1,'th'].value) # m.de is divided by 1000 because it is scaled by 1000 in the model
    e_z = back_z + Le* np.sin(m.de[i].value/scale+ m.X[i,1,'th'].value)
    
    # wing positions
    w_front_x = m.X[i,1,'x'].value + 0.19 * np.cos(m.X[i,1,'th'].value)
    w_front_z = m.X[i,1,'z'].value - 0.19 * np.sin(m.X[i,1,'th'].value)
    
    w_back_x = m.X[i,1,'x'].value - 0.19 * np.cos(m.X[i,1,'th'].value)
    w_back_z = m.X[i,1,'z'].value + 0.19 * np.sin(m.X[i,1,'th'].value)
    
    # wheel position
    wheel_x = m.X[i,1,'x'].value + m.wheel_leg['x']*np.cos(m.X[i,1,'th'].value) + m.wheel_leg['z']*np.sin(m.X[i,1,'th'].value)
    wheel_z = -m.contact_p[i,'wheel'].value
    
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
    ax1.plot([w_back_x,w_front_x],[w_back_z,w_front_z],'C0', linewidth = 5)
    
    # com
    ax1.plot(m.X[i,1,'x'].value,m.X[i,1,'z'].value,'k.')
    
    # legs
    ax1.plot([m.top[i,2,'x'].value,m.bottom[i,2,'x'].value],[m.top[i,2,'z'].value,m.bottom[i,2,'z'].value], 'k')
    ax1.plot([m.top[i,3,'x'].value,m.bottom[i,3,'x'].value],[m.top[i,3,'z'].value,m.bottom[i,3,'z'].value], 'k')
    
    ax1.plot(m.top[i,3,'x'].value,m.top[i,3,'z'].value, '.k')
    
    # wheel
    ax1.plot(wheel_x, wheel_z,'.k')
    
    # path
    ax1.plot(x[0:i],z[0:i],':', color = "Grey")
    
    # ------------------------
    
    ax2.clear()
    ax2.set_xlim([m.X[i,1,'x'].value-1,m.X[i,1,'x'].value+1])
    ax2.set_ylim([m.X[i,1,'z'].value+1,m.X[i,1,'z'].value-1])
    
    # ground
    ax2.plot([-100,200],[0,0], 'C2')
    
    # body
    ax2.plot([front_x,back_x],[front_z,back_z],'grey')
    ax2.plot(m.X[i,1,'x'].value,m.X[i,1,'z'].value,'k.')
    
    # elevator
    ax2.plot([back_x,e_x],[back_z,e_z],'C0', linewidth = 5)
    ax2.plot(back_x,back_z,'.',color="Grey")
    
    # wing
    ax2.plot([w_back_x,w_front_x],[w_back_z,w_front_z],'C0', linewidth = 5)
    
    # com
    ax2.plot(m.X[i,1,'x'].value,m.X[i,1,'z'].value,'ko')
    
    # legs
    ax2.plot([m.top[i,2,'x'].value, m.bottom[i,2,'x'].value],[m.top[i,2,'z'].value,m.bottom[i,2,'z'].value],'k',linewidth = 2)
    ax2.plot([m.top[i,3,'x'].value, m.bottom[i,3,'x'].value],[m.top[i,3,'z'].value,m.bottom[i,3,'z'].value],'k',linewidth = 2)
    
    ax2.plot(m.top[i,2,'x'].value,m.top[i,2,'z'].value, '.k')
    ax2.plot(m.top[i,3,'x'].value,m.top[i,3,'z'].value, '.k')
    
    # wheel
    if wheel_z >= -1e-7:
        ax2.plot(wheel_x, wheel_z,'ro')
    else:
        ax2.plot(wheel_x, wheel_z,'ko')
        
    if m.bottom[i,3,'z'].value >= -1e-7:
        ax2.plot(m.bottom[i,3,'x'].value, m.bottom[i,3,'z'].value,'.r')
    else:
        ax2.plot(m.bottom[i,3,'x'].value, m.bottom[i,3,'z'].value,'.k')
    
    ax2.plot(x[0:i],z[0:i],':', color = "Grey")
    
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.grid()


def print_results(m, results, start_time = None, finish_time = None):
    '''
        Prints out the result from pyomo opimizer, the solver status, the final value of the cost function
        and the time the solver took in minutes. If the result is infesable then the constraints violated are
        printed as well.
        
        Inputs:
        results - Pyomo
        
        Optional:
        start_time - start time in seconds (I used time.time() from the time library)
        finish_time - finish time in seconds (I used time.time() from the time library)
        '''
    print(results.solver.status)
    print(results.solver.termination_condition)
    print(m.Cost.expr())
    if start_time != None and finish_time != None:
        print('Time:',(finish_time-start_time)/60,'min')
    if results.solver.termination_condition == "infeasible":
        print('___________________________________________________________________')
        from pyomo.util.infeasible import log_infeasible_constraints
        log_infeasible_constraints(m)

def init_from_solution(m, solution_path, legs):
    '''
        Initializes pyomo model from a previous solution.
        
        -m is the pyomo model
        -solution_path is the path to the pickle file with the info of the prevoius solution.
        - legs refers to the number of legs the model has, it which model, (int)(0,1,or 2)
        '''
    infile = open(solution_path, 'rb')
    data = pkl.load(infile)
    infile.close()
    
    N = data['N']
    if N > len(m.N):
        print("The solution has more nodes thant the current model has.")
        print("Adding values for all nodes of the model")
        N = len(m.N)
    
    if N < len(m.N):
        print("The solution has less nodes thant the current model has.")
        print("Adding values up to node ",N)

    if legs < 2:
        print("Sorry Still need to do this one, nothing done")
        return(m)
    if legs == 2:
        for n in range(1,N+1):
            m.alpha[n].value = data['a'][n-1]
            m.de[n].value = data['de'][n-1]
            m.df[n].value = data['df'][n-1]
            m.V[n].value = data['V'][n-1]
            m.T[n].value = data['T'][n-1]
            
            m.X[n,1,'x'].value = data['x1'][n-1]
            m.X[n,2,'x'].value = data['x2'][n-1]
            m.X[n,3,'x'].value = data['x3'][n-1]
            m.X[n,4,'x'].value = data['x4'][n-1]
            m.X[n,5,'x'].value = data['x5'][n-1]
            
            m.X[n,1,'z'].value = data['z1'][n-1]
            m.X[n,2,'z'].value = data['z2'][n-1]
            m.X[n,3,'z'].value = data['z3'][n-1]
            m.X[n,4,'z'].value = data['z4'][n-1]
            m.X[n,5,'z'].value = data['z5'][n-1]
            
            m.X[n,1,'th'].value = data['th1'][n-1]
            m.X[n,2,'th'].value = data['th2'][n-1]
            m.X[n,3,'th'].value = data['th3'][n-1]
            m.X[n,4,'th'].value = data['th4'][n-1]
            m.X[n,5,'th'].value = data['th5'][n-1]
            #---
            m.dX[n,1,'x'].value = data['dx1'][n-1]
            m.dX[n,2,'x'].value = data['dx2'][n-1]
            m.dX[n,3,'x'].value = data['dx3'][n-1]
            m.dX[n,4,'x'].value = data['dx4'][n-1]
            m.dX[n,5,'x'].value = data['dx5'][n-1]
            
            m.dX[n,1,'z'].value = data['dz1'][n-1]
            m.dX[n,2,'z'].value = data['dz2'][n-1]
            m.dX[n,3,'z'].value = data['dz3'][n-1]
            m.dX[n,4,'z'].value = data['dz4'][n-1]
            m.dX[n,5,'z'].value = data['dz5'][n-1]
            
            m.dX[n,1,'th'].value = data['dth1'][n-1]
            m.dX[n,2,'th'].value = data['dth2'][n-1]
            m.dX[n,3,'th'].value = data['dth3'][n-1]
            m.dX[n,4,'th'].value = data['dth4'][n-1]
            m.dX[n,5,'th'].value = data['dth5'][n-1]
            #---
            m.ddX[n,1,'x'].value = data['ddx1'][n-1]
            m.ddX[n,2,'x'].value = data['ddx2'][n-1]
            m.ddX[n,3,'x'].value = data['ddx3'][n-1]
            m.ddX[n,4,'x'].value = data['ddx4'][n-1]
            m.ddX[n,5,'x'].value = data['ddx5'][n-1]
            
            m.ddX[n,1,'z'].value = data['ddz1'][n-1]
            m.ddX[n,2,'z'].value = data['ddz2'][n-1]
            m.ddX[n,3,'z'].value = data['ddz3'][n-1]
            m.ddX[n,4,'z'].value = data['ddz4'][n-1]
            m.ddX[n,5,'z'].value = data['ddz5'][n-1]
            
            m.ddX[n,1,'th'].value = data['ddth1'][n-1]
            m.ddX[n,2,'th'].value = data['ddth2'][n-1]
            m.ddX[n,3,'th'].value = data['ddth3'][n-1]
            m.ddX[n,4,'th'].value = data['ddth4'][n-1]
            m.ddX[n,5,'th'].value = data['ddth5'][n-1]
            #---
            m.Tc[n,2].value = data['tau1'][n-1]
            m.Tc[n,3].value = data['tau2'][n-1]
            m.Tc[n,4].value = data['tau3'][n-1]
            m.Tc[n,5].value = data['tau4'][n-1]
            
            
            for gc in m.ground_constraints.keys():
                m.ground_penalties[n,'wheel',gc].value = data['wheel penaltys'][gc][n-1]
                m.ground_penalties[n,'foot1',gc].value = data['foot1 penaltys'][gc][n-1]
                m.ground_penalties[n,'foot2',gc].value = data['foot2 penaltys'][gc][n-1]
        
            m.GRF[n,'wheel','z','ng'].value = data['wheel forces']['normal'][n-1]
            m.GRF[n,'wheel','x','ng'].value = data['wheel forces']['x_ng'][n-1]
            m.GRF[n,'wheel','x','ps'].value = data['wheel forces']['x_ps'][n-1]
            
            m.GRF[n,'foot1','z','ng'].value = data['foot1 forces']['normal'][n-1]
            m.GRF[n,'foot1','x','ng'].value = data['foot1 forces']['x_ng'][n-1]
            m.GRF[n,'foot1','x','ps'].value = data['foot1 forces']['x_ps'][n-1]
            
            m.GRF[n,'foot2','z','ng'].value = data['foot2 forces']['normal'][n-1]
            m.GRF[n,'foot2','x','ng'].value = data['foot2 forces']['x_ng'][n-1]
            m.GRF[n,'foot2','x','ps'].value = data['foot2 forces']['x_ps'][n-1]

    print("Values updated")
    return(m)

def init_linear(m,x_param,z_param, trim_data, BW):
    '''
        Initialtes a linear path from the start point to the finish point.
        Top part of the leg stays at 0 degrees form o to 3N/2 and them moves from 0 to -90 degrees.
        Bottom part of the leg stays at 45 degrees.
        
        Inputs:
        m - Pyomo model
        x_param - start and finish x values (tuple of floats)
        z_param - start and finish z values (tuple of floats)
        trim_data - dictionary for trim conditions
        BW - body weight of the model (Newtoms, none of this Kg stuff)
        
        Returns:
        m - Pyomo model
        '''
    N = len(m.N)
    start_height = z_param[0]
    finish_height = z_param[1]
    
    start_x = x_param[0]
    finish_x = x_param[1]
    
    alpha_t = trim_data['alpha_t']
    delta_e_t = trim_data['delta_e_t']
    theta_t = alpha_t
    T_t = trim_data['T_t']
    V_t = trim_data['V_t']
    
    guess_x1 = np.linspace(start_x,finish_x, N)
    guess_z1 = np.linspace(start_height,finish_height, N)
    
    guess_th1 = np.ones(N)*theta_t
    guess_de = np.ones(N)*delta_e_t
    guess_df = np.ones(N)* np.pi/18
    guess_V = np.ones(N)*V_t - guess_x1*18/finish_x
    guess_T = np.ones(N)*T_t - guess_x1*T_t/finish_x
    
    guess_th2 = np.zeros(round(N*2/3))
    guess_th2 = np.append(guess_th2, np.linspace(0.0, -np.pi/2, round(N/3)))
    guess_th3 = np.ones(N)*np.pi/4
    
    guess_x2 = guess_x1 + 0.5*m.len[2]*np.sin(guess_th2)
    guess_x3 = guess_x1 + m.len[2]*np.sin(guess_th2) + 0.5 * m.len[3]*np.sin(guess_th3)
    
    guess_z2 = guess_z1 + 0.5*m.len[2]*np.cos(guess_th2)
    guess_z3 = guess_z1 + m.len[2]*np.cos(guess_th2) + 0.5 * m.len[3]*np.cos(guess_th3)
    
    guess_tau2 = - m.m[2]*m.g*0.5*m.len[2]*np.sin(guess_th2)/BW
    guess_tau3 = - m.m[3]*m.g*0.5*m.len[3]*np.sin(guess_th3)/BW
    
    for n in range(1,N+1):
        m.alpha[n].value = alpha_t
        m.de[n].value = guess_de[n-1]
        m.df[n].value = guess_df[n-1]
        m.V[n].value = guess_V[n-1]
        m.T[n].value = guess_T[n-1]
        
        m.X[n,1,'x'].value = guess_x1[n-1]
        m.X[n,2,'x'].value = guess_x2[n-1]
        m.X[n,3,'x'].value = guess_x3[n-1]
        m.X[n,4,'x'].value = guess_x2[n-1]
        m.X[n,5,'x'].value = guess_x3[n-1]
        
        m.X[n,1,'z'].value = guess_z1[n-1]
        m.X[n,2,'z'].value = guess_z2[n-1]
        m.X[n,3,'z'].value = guess_z3[n-1]
        m.X[n,4,'z'].value = guess_z2[n-1]
        m.X[n,5,'z'].value = guess_z3[n-1]
        
        m.X[n,1,'th'].value = guess_th1[n-1]
        m.X[n,2,'th'].value = guess_th2[n-1]
        m.X[n,3,'th'].value = guess_th3[n-1]
        m.X[n,4,'th'].value = guess_th2[n-1]
        m.X[n,5,'th'].value = guess_th3[n-1]
        
        m.Tc[n,2].value = guess_tau2[n-1]
        m.Tc[n,3].value = guess_tau3[n-1]
        m.Tc[n,4].value = guess_tau2[n-1]
        m.Tc[n,5].value = guess_tau3[n-1]
    print("Values Updated :)")
    return m


def plot_UAV_legs(i,m,ax1, ax2, xmin, xmax, ymin, ymax, angle_scale):
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
    angle_scale - value by which angles have been scaled
    '''
    x = list(m.X[:,1,'x'].value)
    z = list(m.X[:,1,'z'].value)
    
    L1 = 0.625
    L2 = 0.5
    Le = 0.125

    # body positions
    front_x = m.X[i,1,'x'].value + L2* np.cos(m.X[i,1,'th'].value/angle_scale)
    front_z = m.X[i,1,'z'].value - L2* np.sin(m.X[i,1,'th'].value/angle_scale)
    
    back_x = m.X[i,1,'x'].value - L1* np.cos(m.X[i,1,'th'].value/angle_scale)
    back_z = m.X[i,1,'z'].value + L1* np.sin(m.X[i,1,'th'].value/angle_scale)
    
    # elevator positions
    e_x = back_x - Le* np.cos(m.de[i].value/angle_scale+ m.X[i,1,'th'].value/angle_scale) # m.de is divided by 1000 because it is scaled by 1000 in the model
    e_z = back_z + Le* np.sin(m.de[i].value/angle_scale+ m.X[i,1,'th'].value/angle_scale)
    
    # wing positions
    w_front_x = m.X[i,1,'x'].value + 0.19 * np.cos(m.X[i,1,'th'].value/angle_scale)
    w_front_z = m.X[i,1,'z'].value - 0.19 * np.sin(m.X[i,1,'th'].value/angle_scale)
    
    w_back_x = m.X[i,1,'x'].value - 0.1 * np.cos(m.X[i,1,'th'].value/angle_scale)
    w_back_z = m.X[i,1,'z'].value + 0.1 * np.sin(m.X[i,1,'th'].value/angle_scale)
    
    # flap
    
    flap_x = w_back_x - 0.09 * np.cos(m.df[i].value/angle_scale+ m.X[i,1,'th'].value/angle_scale)
    flap_z = w_back_z + 0.09 * np.sin(m.df[i].value/angle_scale+ m.X[i,1,'th'].value/angle_scale)
    
    # wheel position
    wheel_x = m.X[i,1,'x'].value + m.wheel_leg['x']*np.cos(m.X[i,1,'th'].value/angle_scale) + m.wheel_leg['z']*np.sin(m.X[i,1,'th'].value/angle_scale)
    wheel_z = -m.contact_p[i,'wheel'].value
    
    # wheel leg
    wheel_top_x = m.X[i,1,'x'].value + m.wheel_leg['x']*np.cos(m.X[i,1,'th'].value/angle_scale)
    wheel_top_z = m.X[i,1,'z'].value - m.wheel_leg['x']*np.sin(m.X[i,1,'th'].value/angle_scale)
    
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
    ax1.plot(m.X[i,1,'x'].value,m.X[i,1,'z'].value,'k.')
    
    # legs
    ax1.plot([m.top[i,2,'x'].value,m.bottom[i,2,'x'].value],[m.top[i,2,'z'].value,m.bottom[i,2,'z'].value], 'k')
    ax1.plot([m.top[i,3,'x'].value,m.bottom[i,3,'x'].value],[m.top[i,3,'z'].value,m.bottom[i,3,'z'].value], 'k')
    
    ax1.plot([m.top[i,4,'x'].value,m.bottom[i,4,'x'].value],[m.top[i,4,'z'].value,m.bottom[i,4,'z'].value], color = 'grey')
    ax1.plot([m.top[i,5,'x'].value,m.bottom[i,5,'x'].value],[m.top[i,5,'z'].value,m.bottom[i,5,'z'].value], color = 'grey')
    
    ax1.plot(m.top[i,3,'x'].value,m.top[i,3,'z'].value, '.k')
    ax1.plot(m.top[i,5,'x'].value,m.top[i,5,'z'].value, color = 'grey')
    
    # wheel
    ax1.plot([wheel_top_x, wheel_x],[wheel_top_z, wheel_z], color = 'grey')
    ax1.plot(wheel_x, wheel_z,'.k')
    
    # path
    ax1.plot(x[0:i],z[0:i],':', color = "Grey")
    
    # ------------------------
    
    ax2.clear()
    ax2.set_xlim([m.X[i,1,'x'].value-1,m.X[i,1,'x'].value+1])
    ax2.set_ylim([m.X[i,1,'z'].value+1,m.X[i,1,'z'].value-1])
    
    # ground
    ax2.plot([-100,200],[0,0], 'C2')
    
    # body
    ax2.plot([front_x,back_x],[front_z,back_z],'grey')
    ax2.plot(m.X[i,1,'x'].value,m.X[i,1,'z'].value,'k.')
    
    # elevator
    ax2.plot([back_x,e_x],[back_z,e_z],'C0', linewidth = 5)
    ax2.plot(back_x,back_z,'.',color="Grey")
    
    # wing
    ax2.plot([flap_x,w_back_x,w_front_x],[flap_z,w_back_z,w_front_z],'C0', linewidth = 5)
    
    # legs
    ax2.plot([m.top[i,4,'x'].value,m.bottom[i,4,'x'].value],[m.top[i,4,'z'].value,m.bottom[i,4,'z'].value], color = 'grey',linewidth = 2)
    ax2.plot([m.top[i,5,'x'].value,m.bottom[i,5,'x'].value],[m.top[i,5,'z'].value,m.bottom[i,5,'z'].value], color = 'grey',linewidth = 2)
    
    ax2.plot(m.top[i,4,'x'].value,m.top[i,4,'z'].value, '.', color = 'grey')
    ax2.plot(m.top[i,5,'x'].value,m.top[i,5,'z'].value, '.', color = 'grey')
    
    
    ax2.plot([m.top[i,2,'x'].value, m.bottom[i,2,'x'].value],[m.top[i,2,'z'].value,m.bottom[i,2,'z'].value],'k',linewidth = 2)
    ax2.plot([m.top[i,3,'x'].value, m.bottom[i,3,'x'].value],[m.top[i,3,'z'].value,m.bottom[i,3,'z'].value],'k',linewidth = 2)
    
    ax2.plot(m.top[i,2,'x'].value,m.top[i,2,'z'].value, '.k')
    ax2.plot(m.top[i,3,'x'].value,m.top[i,3,'z'].value, '.k')
    
    
    # com
    ax2.plot(m.X[i,1,'x'].value,m.X[i,1,'z'].value,'ko')
    
    # wheel
    ax2.plot([wheel_top_x, wheel_x],[wheel_top_z, wheel_z], color = 'grey')
    if wheel_z >= -1e-7:
        ax2.plot(wheel_x, wheel_z,'ro')
    else:
        ax2.plot(wheel_x, wheel_z,'ko')
        
    if m.bottom[i,3,'z'].value >= -1e-7:
        ax2.plot(m.bottom[i,3,'x'].value, m.bottom[i,3,'z'].value,'.r')
    else:
        ax2.plot(m.bottom[i,3,'x'].value, m.bottom[i,3,'z'].value,'.k')
        
    if m.bottom[i,5,'z'].value >= -1e-7:
        ax2.plot(m.bottom[i,5,'x'].value, m.bottom[i,5,'z'].value,'.r')
    else:
        ax2.plot(m.bottom[i,5,'x'].value, m.bottom[i,5,'z'].value,'.', color = 'grey')
    
    ax2.plot(x[0:i],z[0:i],':', color = "Grey")
    
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.grid()
