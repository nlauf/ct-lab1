; Electrostatic Particle Cyclic 1D
;
; For studying 2-stream instability and Landau damping
;
; Landau damping
; 
; deltaf, len=4.*!PI, alpha=0.2
;
; For 2-stream
; 
; deltaf, len=100, Vbeam=3

; Function to convert a value into a string
FUNCTION str, var
  RETURN, STRTRIM(STRING(var),2)
END

; Calculate the density given particle positions
FUNCTION calc_density, position, weight, ncells, L
  density = FLTARR(ncells)
  
  dx = FLOAT(L) / ncells

  ; Get the lower cell index
  plower = LONG(position)
  
  FOR i=0,ncells-1 DO BEGIN
    w = WHERE(plower EQ i, count)
    IF count GT 0 THEN BEGIN
      ni = TOTAL( weight[w]*(FLOAT(i) - position[w] + 1.0) )
      density[i]   = density[i] + ni
      
      ip = (i+1) MOD ncells
      ni = TOTAL( weight[w]*(position[w] - FLOAT(i)) )
      density[ip] = density[ip] + ni
    ENDIF
  ENDFOR
  
  ; Need to divide by dx to get a number density
  RETURN, density / dx
END

; Interpolate a periodic set of values v at points x
FUNCTION periodic_interp, v, x
  nv = N_ELEMENTS(v)
  
  xl = LONG(x)
  
  RETURN, v[xl]*(xl - x + 1.0) + v[(xl + 1) MOD nv] * (x - xl)
END

; Integrate a periodic function using FFTs
FUNCTION fft_integrate, y
  on_error,2  ; If an error occurs, return to caller
  
  n = N_ELEMENTS(y)

  F = FFT(y)
  imag = complex(0.0, 1.0)

  F[0] = 0.0

  IF (n MOD 2) EQ 0 THEN BEGIN
      ; even number of points

      FOR i=1l, n/2-1 DO BEGIN
          a = imag*2.0*!PI*FLOAT(i)/FLOAT(n)
          F[i] = F[i] / a         ; positive frequencies
          F[n-i] = - F[n-i] / a   ; negative frequencies
      ENDFOR

      F[n/2] = F[n/2] / (imag*!PI)
  ENDIF ELSE BEGIN
      ; odd number of points

      FOR i=1l, (n-1)/2 DO BEGIN
          a = imag*2.0*!PI*FLOAT(i)/FLOAT(n)
          F[i] = F[i] / a
          F[n-i] = - F[n-i] / a 
      ENDFOR
  ENDELSE

  result = FFT(F, 1)
  
  RETURN, REAL_PART(result)
END

; Calculate starting distribution function and derivatives
; Gaussian distribution in velocity, constant in space
FUNCTION f0, x, v, dfdx=dfdx, dfdv=dfdv
  f = exp(-0.5*v^2) / sqrt(2.*!PI)   ; Gaussian, density and standard deviation of 1
  
  dfdv = -v * f     ; Derivative of a Gaussian
  dfdx = 0.0  ; Zero, as f0 is constant in x
  
  RETURN, f
END

; Evaluates time-derivatives for RK4 method
FUNCTION eval_derivs, t, y, debug=debug
  COMMON epc_com, L, ncells, nparticles, density, E

  dx = L / FLOAT(ncells)
  
  ; Extract position and velocity
  pos = ((y[0:(nparticles-1)] MOD ncells) + ncells) MOD ncells
  vel = y[nparticles:(2*nparticles-1)]
  weight = y[(2*nparticles):(3*nparticles-1)]
  
  ; Gather particles to get density
  density = calc_density(pos, weight, ncells, L)
  
  ; Calculate charge density assuming stationary ions
  n0 = FLOAT(nparticles) / FLOAT(L)
  rho = density / n0

  ; Calculate electric field on cells
  E = -fft_integrate(rho)*dx

  ; Calculate E field at particle locations
  accel = -periodic_interp(E, pos) / dx ; convert into cell
  
  ; get the starting distribution
  f = f0(pos, vel*dx, dfdx=dfdx, dfdv=dfdv)

  dwdt = (weight - 1) * (vel * dfdx  + accel*dx * dfdv) / f
  
  IF KEYWORD_SET(debug) THEN STOP

  RETURN, [vel, accel, dwdt]
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Main program                                                ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO deltaf, len=len, $               ; Length of the box
           ncell=ncell, $           ; Number of cells to use
           npart=npart, $           ; Number of particles
           outputstep=outputstep, $ ; Time between outputs
           noutputs=noutputs, $     ; Number of outputs
           seed=seed, $             ; Initial random seed
           Vbeam=Vbeam, $           ; Velocity of each beam
           alpha=alpha, $           ; Perturbation amplitude
           kwave=kwave, $           ; The k number of the wave (1 = lowest harmonic)
           delay=delay, $           ; Time delay after plotting 
           result=result, $         ; Structure containing the result
           filename=filename        ; Name of the file for the output
  
  ; Common block, used to pass variables through the RK4 routine
  ; to eval_derivs. Not a nice feature of IDL but necessary here
  COMMON epc_com, L, ncells, nparticles, density, E

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Simulation parameters
  
  IF NOT KEYWORD_SET(len) THEN len = 4.*!PI
  L = len
  
  IF NOT KEYWORD_SET(ncell) THEN ncell = 20 ; Default number of cells
  ncells = ncell
  
  IF NOT KEYWORD_SET(npart) THEN npart = 200000L ; Default number of particles
  nparticles = npart
  
  IF NOT KEYWORD_SET(outputstep) THEN outputstep = 1.0 ; Time between outputs
  IF NOT KEYWORD_SET(noutputs) THEN noutputs = 40 ; Number of outputs
  
  IF NOT KEYWORD_SET(seed) THEN seed = 1234 ; Random number seed
  
  IF NOT KEYWORD_SET(Vbeam) THEN Vbeam = 0.  ; Beam velocity
  
  IF NOT KEYWORD_SET(alpha) THEN alpha = 0.  ; Perturbation amplitude
  
  IF NOT KEYWORD_SET(kwave) THEN kwave = 1.  ; Number of periods
  k = kwave * 2*!PI/L

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Print out header
  
  PRINT, "Electrostatic PIC code"
  PRINT, "   Ben Dudson, University of York, Oct 2009"
  PRINT, ""  
  PRINT, "Domain size   [Debye lengths]: "+STRTRIM(STRING(L),2)
  PRINT, "Number of cells    : "+STRTRIM(STRING(ncells),2)
  PRINT, "Number of particles: "+STRTRIM(STRING(nparticles),2)
  PRINT, ""
  PRINT, "      Time [Wp]     K.E.         P.E.      Total Energy"
  PRINT, ""
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Open a file for the result

  OPENW, fid, "epc1d.csv", /GET_LUN
  
  PRINTF, fid, "len= "       + STR(L)
  PRINTF, fid, "ncells= "    + STR(ncells)
  PRINTF, fid, "nparticles= "+ STR(nparticles)
  PRINTF, fid, "alpha= "     + STR(alpha)
  PRINTF, fid, "kwave= "     + STR(kwave)
  PRINTF, fid, "seed= "      + STR(seed)
  
  PRINTF, fid, "Time [Wp], K.E., P.E., Total Energy, n(1), n(2), n(3), n(4)"
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Array of particle positions and velocities.
  ; Both are in units of grid cells
  
  position = LONARR(nparticles)
  velocity = LONARR(nparticles)
  weight   = LONARR(nparticles)

  ; Set initial conditions
  
  dx = FLOAT(L) / FLOAT(ncells)

  ; Position is uniform across domain
  position = RANDOMU(seed, nparticles) * ncells
  ; Gaussian velocity distribution
  velocity = RANDOMN(seed+1, nparticles) / dx
  
  ; Perturb the distribution by setting weights
  ; Uses alpha and k (kwave) to set perturbation.
  weight = alpha * COS(kwave*2.*!PI*position/ncells)
  
  ; Make sure average is zero
  weight = weight - MEAN(weight)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Start time integration loop
  
  max_timestep = 1.0
  time = 0.
  first = 1
  step = 1
  nextoutput = outputstep
  ; evaluate derivatives to get acceleration
  dydt = eval_derivs(time, [position, velocity, weight])
  accel = dydt[nparticles:(2*nparticles-1)]
  WINDOW, xsize=400, ysize=600

  REPEAT BEGIN
    ; Calculate timestep
    tstep = max_timestep
    ; Maximum distance a particle can move is one cell
    vlim = 1. / MAX(ABS(velocity))
    IF vlim LT tstep THEN tstep = vlim
    ; Acceleration limit (using last calculated accel)
    alim = SQRT(1./MAX(ABS(accel)))
    IF alim LT tstep THEN tstep = alim

    outnow = 0
    IF time + tstep GE nextoutput THEN BEGIN
      ; Going to step over next output time
      tstep = nextoutput - time ; Stop at the next time
      outnow = 1 ; Signal to output
    ENDIF

    WRITEU, -1, 13, "Time = "+STRTRIM(STRING(time),2) + $
      " step = "+STRTRIM(STRING(tstep),2)

    ; Now evolve system using Runge-Kutta 4th order method
    
    y = [position, velocity, weight]
    dydt = eval_derivs(time, y)
    y = RK4(y, dydt, time, tstep, 'eval_derivs')
    time = time + tstep

    ; Extract new positions and velocities from y array
    position = ((y[0:(nparticles-1)] MOD ncells) + ncells) MOD ncells
    velocity = y[nparticles:(2*nparticles-1)]
    weight   = y[(2*nparticles):*]
    accel    = dydt[nparticles:(2*nparticles-1)]
    
    ; Plot the current state
    !P.multi=[0,0,3,0,0]
    plot, L*position/FLOAT(ncells), L*velocity/FLOAT(ncells), psym=3, xr=[0,L], xstyle=1, $
      xtitle="Position [Debye lengths]", ytitle="Velocity [Debye length * Wp]", chars=1.5, $
      title="Time = "+STRTRIM(STRING(time),2)+" / Wp"
    plot, L*FINDGEN(ncells)/FLOAT(ncells), L * density / FLOAT(nparticles), yr=[-0.5, 0.5], $
      xtitle="Position [Debye lengths]", ytitle="Density", chars=1.5
    plot, L*FINDGEN(ncells)/FLOAT(ncells), E, $
      xtitle="Position [Debye lengths]", ytitle="Electric field", chars=1.5
    IF KEYWORD_SET(delay) THEN WAIT, delay
    
    IF outnow THEN BEGIN
      ; Output 
      step = step + 1
      nextoutput = FLOAT(step) * outputstep
      
      n0 = FLOAT(nparticles) / FLOAT(L)

      ; Calculate kinetic energy
      ke = 0.5*TOTAL(velocity^2)*dx*dx / n0 ; convert from cell index

      ; Calculate potential
      rho = density / n0 - 1.
      phi = fft_integrate(fft_integrate(rho))*dx*dx
      ; Calculate potential energy of electrons
      pe = -TOTAL(periodic_interp(phi, position)) / n0
      
      WRITEU, -1, 13
      PRINT, time, ke, pe, ke + pe
      
      fd = FFT(density)

      ; Write data to result file
      PRINTF, fid, str(time)+", "+str(ke)+", "+str(pe)+", "+str(ke + pe) + $
        ", " + str(ABS(fd[1]))+", "+str(ABS(fd[2]))+", "+str(ABS(fd[3]))+", "+str(ABS(fd[4]))
    ENDIF
    
    ;file = "epc1d" + STRTRIM(STRING(frame, FORMAT='(I04)'),2)+".bmp"
    ;frame = frame + 1
    ;PRINT, "Writing file: " +file
    ;WRITE_BMP, file, TVRD(TRUE=1)

    ; Calculate amplitude of fundamental harmonic for result
    fd = FFT(density)
    IF first THEN BEGIN
      result_time = [time]
      result_d1 = [ABS(fd[1])] ; First harmonic
      result_d2 = [ABS(fd[2])] ; Second harmonic
      
      first = 0
    ENDIF ELSE BEGIN
      result_time = [result_time, time]
      result_d1 = [result_d1, ABS(fd[1])]
      result_d2 = [result_d2, ABS(fd[2])]
    ENDELSE
    
  ENDREP UNTIL step GT noutputs

  result = {time:result_time, d1:result_d1, d2:result_d2}

  ; Close the output file
  CLOSE, fid
  
  PRINT, "finished"
  
END
