README

AUTHORS:
Mark Corrigan	s1304215
Neil McBlane	s1305440

TERMINAL FORMAT:
run as $: java ParticleNBody XBody.txt Parameters.txt TrajX.xyz ChecksX.txt OrbitsX.txt
	where X is Sol for solar system or 3 for three body system

FILE DETAILS:
SolarSys.txt:
Contains the Particle3D input in the following format:
		name1 x1 y1 z1 vx1 vy1 vz1 mass1
		name2 x2 y2 z2 vx2 vy2 vz2 mass2
		and so on...

Where 	name is a string labelling each particle
	x,y,z are doubles representing the initial position, units [AU]
	vx,vy,vz are doubles representing the initial velocity, units [AU day^-1]
	mass is a double representing the particle mass, units [SolarMass]


Parameters.txt
Contains the simulation parameters input in the following format:
		G
		timestep
		number of timesteps

Where	G = 4pi^2/365^2, units [AU^3 day^-2 solarMass^-1]
	timestep = 0.01 a double dicating the separation of each calculation found to
	 be the best reasonable tradeoff between accuracy and efficiency, units [day]
	number of timesteps = 10,000,000 is an integer dictating the total number of
	 steps simulated an is sufficient to cover just over a full Pluto orbit


Traj.xyz
Contains the VMD output file in the default format, units [AU]


Checks.txt
Contains the parameters which verify the accuracy of the simulation
 in the following format:
		Time, totalEnergy, totalMomentum
Where 	Time is the time the simulation has been running for, units [day]
	totalEnergy is the kinetic and pair potentials of all the particles and should
	 remain relatively constant, units arbitrary
	totalMomentum is the total momentum of all the particles and should be zero,
	 units arbitrary

Orbits.txt
Contains the characteristics of each orbit in the following format:
		 Planet, No. Orbits, Aphelion, Perihelion, Eccentricity, KeplerIII
Where 	Planet is the name of each particle
	No. Orbits counts orbits completed to one decimal place
	Aphelion records the particle’s furthest distance from its orbit centre, 
	 units [AU]
	Perihelion records the particle’s shortest distance from its orbit centre,
	 units [AU]
	Eccentricity records the particle’s eccentricity
	KeplerIII compares the period and semimajor axis of the orbit to verify
	 Kepler’s third law, a value of 1.000 represents perfect agreement
