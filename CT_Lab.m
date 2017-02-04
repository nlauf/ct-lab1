%!/usr/bin/env python
%
% Electrostatic PIC code in a 1D cyclic domain

    %from scipy import arange, concatenate, zeros, linspace, floor, array, pi
    %from scipy import sin, cos, sqrt, random, histogram

    %from scipy import matplotlib.pyplot as plt %# Matplotlib plotting library

try
    %from scipy import matplotlib.gridspec as gridspec  %# For plot layout grid
    %got_gridspec = True
catch
    %got_gridspec = False
end
% Need an FFT routine, either from SciPy or NumPy
try
    %from scipy.fftpack import fft, ifft
catch
    % No SciPy FFT routine. Import NumPy routine instead
    %from numpy.fft import fft, ifft
end

def rk4step(f, y0, dt, args=()):
    %""" Takes a single step using RK4 method """
    k1 = f(y0, args)
    k2 = f(y0 + 0.5*dt*k1, args)
    k3 = f(y0 + 0.5*dt*k2, args)
    k4 = f(y0 + dt*k3, args)

    disp (y0 + (k1 + 2.*k2 + 2.*k3 + k4)*dt / 6.)


def calc_density(position, ncells, L):
    %""" Calculate charge density given particle positions
    
    %Input
     % position  - Array of positions, one for each particle
      %            assumed to be between 0 and L
      %ncells    - Number of cells
      %L         - Length of the domain

    %Output
     % density   - contains 1 if evenly distributed
    %"""
    % This is a crude method and could be made more efficient
    
    density = zeros (ncells)
    nparticles = len(position)
    
    dx = L / ncells       % Uniform cell spacing
    for p = position / dx    % Loop over all the particles, converting position into a cell number
        plower = int(p)        % Cell to the left (rounding down)
        offset = p - plower    % Offset from the left
        density(plower) = 1. - offset
        density(plower) = density(plower) + 1
        density(plower + 1) = mod(density(plower+1),ncells)
        density(plower+1) = offset +1
       % nparticles now distributed amongst ncells
    density = density*(float(ncells) / float(nparticles))  % Make average density equal to 1
    disp (density)
    end
    
def periodic_interp(y, x):
   % """
    %Linear interpolation of a periodic array y at index x
    
    %Input

    %y - Array of values to be interpolated
    %x - Index where result required. Can be an array of values
    
    %Output
    
    %y(x) with non-integer x
    %"""
    ny = len(y)
    if len(x) > 1
        y = array(y) % Make sure it's a NumPy array for array indexing
    xl = floor(x).astype(int) % Left index
    dx = x - xl
    xl = mod(((mod(xl,ny)) + ny),ny)  % Ensures between 0 and ny-1 inclusive
    y(x)= (y(xl)*(1. - dx) + mod((y(xl+1)),ny)*dx)
    end

def fft_integrate(y):
    %""" Integrate a periodic function using FFTs
    %"""
    n = len(y) % Get the length of y
    
    f = fft(y) % Take FFT
    % Result is in standard layout with positive frequencies first then negative
    % n even: [ f(0), f(1), ... f(n/2), f(1-n/2) ... f(-1) ]
    % n odd:  [ f(0), f(1), ... f((n-1)/2), f(-(n-1)/2) ... f(-1) ]
    
    if mod(n,2) == 0 
        % If an even number of points;
        k = [(arange(0, n/2+1)), (arange(1-n/2, 0))]
    else
        k = [(arange(0, (n-1)/2+1)), (arange( -(n-1)/2, 0))]
    k = 2.*pi*k/n
    
    % Modify frequencies by dividing by ik
    f(1,1) = f(1)/(1j * k(1)) 
    f(0) = 0. % Set the arbitrary zero-frequency term to zero
    
    f = ifft(f).real % Reverse Fourier Transform
    end

def pic(f, ncells, L):
    %""" f contains the position and velocity of all particles
    %"""
    nparticles = len(f)/2     % Two values for each particle
    pos = f(1,nparticles) % Position of each particle
    vel = f(nparticles,1)      % Velocity of each particle

    dx = L / float(ncells)    % Cell spacing

    % Ensure that pos is between 0 and L
    pos = mod((mod(pos,L) + L),L)
    
    % Calculate number density, normalised so 1 when uniform
    density = calc_density(pos, ncells, L)
    
    % Subtract ion density to get total charge density
    rho = density - 1.
    
    % Calculate electric field
    E = -fft_integrate(rho)*dx
    
    % Interpolate E field at particle locations
    accel = -periodic_interp(E, pos/dx)

    % Put back into a single array
    accel = [vel, accel]

%

def run(pos, vel, L, ncells=None, out=[], output_times=linspace(0,20,100), cfl=0.5):
    
    if ncells == None
        ncells = int(sqrt(len(pos))) 
        % A sensible default

    dx = L / float(ncells)
    
    f = [pos, vel];   
    % Starting state
    end
    nparticles = len(pos)
    
    time = 0.0
    for tnext = output_times
        % Advance to tnext
        stepping = True;
        while stepping
            % Maximum distance a particle can move is one cell
            dt = cfl * dx / max(abs(vel));
            if time + dt >= tnext
                % Next time will hit or exceed required output time
                stepping = False
                dt = tnext - time
            %print "Time: ", time, dt
            f = rk4step(pic, f, dt,ncells, L)
            time = time + dt
            end
        % Extract position and velocities
        pos = mod((mod((f(0,nparticles)),L) + L),L)
        vel = f(nparticles,1)
        end
    end
        
        % Send to output functions
        for func = out
            func(pos, vel, ncells, L, time)
        
    g = [pos, vel]
        end
%
% 
% Output functions and classes
%

class Plot:
    %"""
    %Displays three plots: phase space, charge density, and velocity distribution
    %"""
    def __init__(self, pos, vel, ncells, L):
        
        d = calc_density(pos, ncells, L)
        vhist, bins  = histogram(vel, int(sqrt(len(vel))))
        vbins = 0.5*(bins(1,0)+bins(0,-1))
        
        % Plot initial positions
        %if got_gridspec:
            self.fig = plt.figure()
            self.gs = gridspec.GridSpec(4, 4)
            %ax = self.fig.add_subplot(self.gs[0:3,0:3])
            %self.phase_plot = ax.plot(pos, vel, '.')[0]
            ax.set_title('Phase space')
            
            %ax = self.fig.add_subplot(self.gs[3,0:3])
           % self.density_plot = ax.plot(linspace(0, L, ncells), d)[0]
            
            %ax = self.fig.add_subplot(self.gs[0:3,3])
            %self.vel_plot = ax.plot(vhist, vbins)[0]
        %else
            self.fig = plt.figure()
            self.phase_plot = plt.plot(pos, vel, '.')%[0]
            
            self.fig = plt.figure()
            %self.density_plot = plt.plot(linspace(0, L, ncells), d)[0]
            
            self.fig = plt.figure()
            %self.vel_plot = plt.plot(vhist, vbins)[0]
        plt.ion()
        plt.show()
        
    def __call__(self, pos, vel, ncells, L, t):
        d = calc_density(pos, ncells, L)
        vhist, bins  = histogram(vel, int(sqrt(len(vel))))
        vbins = 0.5*(bins(1,0)+bins(0,-1))
        
        self.phase_plot.set_data(pos, vel) % Update the plot
        self.density_plot.set_data(linspace(0, L, ncells), d)
        self.vel_plot.set_data(vhist, vbins)
        plt.draw()
        

class Summary:
    def __init__(self):
        self.t = []
        self.firstharmonic = []
        
    def __call__(self, pos, vel, ncells, L, t):
        % Calculate the charge density
        d = calc_density(pos, ncells, L)
        
        % Amplitude of the first harmonic
        fh = 2.*abs(fft(d)) / ncells
        
        print ('Time:', t, 'First:', fh)
        
        self.t.append(t)
        self.firstharmonic.append(fh)

%
% 
% Functions to create the initial conditions


def landau(npart, L, alpha=0.2):
   % """
    %Creates the initial conditions for Landau damping
    
    %"""
    % Start with a uniform distribution of positions
    posa = random.uniform(0., L, npart)
    pos0 = pos.copy()
    k = 2.*pi / L
    for i = range(10) % Adjust distribution using Newton iterations
        posa = posa - ( pos + alpha*sin(k*pos)/k - pos0 ) / ( 1. + alpha*cos(k*pos) )
        
    % Normal velocity distribution
    vel = random.normal(0.0, 1.0, npart)
    
    h= [posa, vel]
    end
def twostream(npart, L, vbeam=2):
    % Start with a uniform distribution of positions
    pos1 = random.uniform(0., L, npart)
    % Normal velocity distribution
    vel = random.normal(0.0, 1.0, npart)
    
    np2 = int(npart / 2)
    vel(0,np2) = vel(0,np2)+vbeam  % Half the particles moving one way
    vel(np2,0) = vel(np2,0)-vbeam  % and half the other
    
    k= [pos1,vel]
%

%if __name__ == "__main__":
    % Generate initial condition
    %
    if False
        % 2-stream instability
        L = 100
        ncells = 20
        pos, vel = twostream(10000, L, 3.)
    else
        % Landau damping
        L = 4.*pi
        ncells = 20
        npart = 1000
        pos, vel = landau(npart, L)
    end
    % Create some output classes
    p = Plot(pos, vel, ncells, L) % This displays an animated figure
    s = Summary()                 % Calculates, stores and prints summary info
    
    % Run the simulation
    pos, vel = run(pos, vel, L, ncells) 
                   out=[p, s],                      % These are called each output
                   output_times=linspace(0.,20,50) % The times to output
    
    % Summary stores an array of the first-harmonic amplitude
    % Make a semilog plot to see exponential damping
    plt.figure()
    plt.plot(s.t, s.firstharmonic)
    plt.xlabel('Time [Normalised]')
    plt.ylabel('First harmonic amplitude [Normalised]')
    plt.yscale('log')
    
    plt.ioff() % This so that the windows stay open
    plt.show()
    
    

