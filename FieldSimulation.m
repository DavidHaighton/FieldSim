#==== Description ====
#{
A basic classical particle simualator for classical electric, magnetic and gravitational fields.
It uses Eulers method to update the simulation.
#}

#==== Constants ====

#Gravitational Constant
global G = 6.67408*10^-11; #m^3*kg^-1^s-2
#Coulombs Constant
global K = 8.9875517923*10^9; # kg*m^3*s*-2*C*-2.
# Permeability of space
global U= 4*pi*10^-7; #T*m*C*s^-1
global MAGNETIC_CONST = U/(4*pi);#used to make calculations faster

#==== Functions ====

function ret = createParticle(name, mass, charge, pos, vel)
  #used to create a particle with given initial conditions
  ret.name = name;
  ret.mass= mass;
  ret.charge = charge;
  ret.pos = pos;
  ret.vel = vel;
  ret.acc = [0,0,0];
endfunction


function force = getForce(particle1, particle2)
  #returns force in direction of particle 2
  global G;
  global K;
  global MAGNETIC_CONST;
  #intermediate calculations
  posDiff = particle2.pos-particle1.pos;
  velDiff = particle2.vel = particle1.vel;
  chargeProd = particle2.charge * particle1.charge;
  massProd = particle2.mass*particle1.mass;
  rrs = posDiff/norm(posDiff)^3; #reciprical r^2 (r^-2)
  #calculate forces
  magneticForce = MAGNETIC_CONST*chargeProd*cross(velDiff, rrs);
  otherForces = (K*chargeProd+G*massProd)*rrs; #faster to do both at once
  #calculate total force
  force = magneticForce+otherForces;
endfunction

function particle = update(particleParam, dt)
    #updates position and velocity
    particleParam.pos += particleParam.vel*dt;
    particleParam.vel += particleParam.acc*dt;
    particle = particleParam;
endfunction

function particles = updateAll(particlesParam, dt)
  #updates all particles
  for i=1:length(particlesParam)
    particlesParam(i) = update(particlesParam(i), dt);
  end
  particles = particlesParam;
endfunction

function particles = simulationStep(particlesParam, dt)
  #runs the simlation for a specified time dt
  particles = particlesParam;
  for i=1:length(particlesParam)
    for j=i+1:length(particlesParam)
        force = getForce(particlesParam(i), particlesParam(j));
        particles(i).acc = force/particles(i).mass;
        particles(j).acc = -force/particles(j).mass;
    end
  end
  #update their velocities and positions
  particles = updateAll(particles, dt);
endfunction



function results = runSimulation(particles, simulationTime, dt)
  #runs and plots the simulation given the particles
    hold on;
  xlabel("x");
  ylabel("y");
  zlabel("z");
  results = zeros(length(particles),simulationTime/dt,3);
  for i=1:length(particles)
    results(i,1,:) = particles.pos;
  end
  for i=2:ceil(simulationTime/dt);
    t= i*dt;
    particles=simulationStep(particles,dt);
    for particleIndex=1:length(particles)
      particle = particles(particleIndex);
      results(particleIndex, i, :) = particle.pos;
      end
  end
endfunction

function plotSimulationResult(result, particles)
  #plot positions of all particles
  hold on;
  
  title("Particles tragectories graphed against eachother")
  for i=1:length(result(:,1,1))
      plot3(result(i,:,1),result(i,:,2),result(i,:,3));
  end
  legend(particles(:).name);
  hold off;
endfunction


#====Example====
mass = 5.972*10^24; #kg
r = 384_400_000; #m
v=sqrt(G*mass/r); #m/s
time = 2_629_800; #s
step = 100; # s
earth = createParticle("earth", mass, 0, [0,0,0], [0,0,0]);
moon = createParticle("moon" , 1, 0, [r,0,0],[0,v,0]);

particles = [earth, moon];

result=runSimulation(particles, time, step);
plotSimulationResult(result, particles);
