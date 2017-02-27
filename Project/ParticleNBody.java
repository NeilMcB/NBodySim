/**
 * Computer Modelling, N-Body Sumulation Submission
 *
 * @author M. Corrigan
 * @author N. McBlane
 * @version 3/2016
 */

// Packages required for file writing, maths and scanning
import java.io.*;
import java.lang.Math.*;
import java.util.*;

public class ParticleNBody {
    public static void main (String [] argv) throws IOException {
        
        /*
         * Initialize Scanners
         */

        // To track number of particles
        Scanner lnCount = new Scanner(new File(argv[0]));
        int c = 0;
        while (lnCount.hasNextLine()){
            c++;
            lnCount.nextLine();
        }

        // To read in particles
        Scanner scParticle = new Scanner(new File(argv[0]));
        // Array to store c particles
        Particle3D[] ps;
        ps = new Particle3D[c];
        // Loop through input file to create particles
        for(int i=0; i<ps.length; i++){
            ps[i] = new Particle3D(scParticle);
        }
        
        // To read simulation algorithm parameters
        Scanner scProps = new Scanner(new File(argv[1]));
        // Initialize parameters
        double g = scProps.nextDouble();
        double dt = scProps.nextDouble();
        int numstep = scProps.nextInt();
        double t = 0.0;
        
        // Initialise orbit counter
        double[] oCounts;
        oCounts = new double[ps.length];

        /*
         * Initialize Printers
         */
        
        // For trajectory file output
        PrintWriter pwTraj = new PrintWriter(new File(argv[2]));       
        // For energy and total momentum checks file output
        PrintWriter pwCheck = new PrintWriter(new File(argv[3]));        
        // For orbit information file output
        PrintWriter pwOrbit = new PrintWriter(new File(argv[4]));

        /*
         * Calculate initial properties
         */

        // Initial force on each particle
        Vector3D[] currentFors;
        currentFors = Particle3D.manyForce(ps, g);
        
        // Initial total momentum of the system
        Vector3D totP = new Vector3D();
        totP = Particle3D.manyMomentum(ps);
        
        // Adjust initial velocities to correct for centre of
        // mass motion using v=P/M
        ps = Particle3D.correctMomentum(ps);
        
        // Corrected initial total energy and momentum of
        // the system
        double totE = Particle3D.manyEnergy(ps, g);
        totP = Particle3D.manyMomentum(ps);

        // Determine the initial position and distance of
        // particle from origin
        Vector3D[] initialPos;
        initialPos = new Vector3D[ps.length];
        for(int i=0; i<ps.length; i++){
            initialPos[i] = ps[i].getPosition();
        }

        // Initialize matrix to store perihelions
        // and aphelions
        double[][] ap;
        ap = new double[ps.length][2];
        for(int i=0; i<ps.length; i++){
            // Aphelion
	    // Define Moon relative to Earth
            if(i==4){
                ap[i][0]=(Vector3D.subVec(initialPos[i],initialPos[i-1])).mag();
            }
            else{
                ap[i][0] = initialPos[i].mag();
            }   
            // Perihelion         
            ap[i][1] = initialPos[i].mag();
        }

        // Perform Verlet algorithm
        /*
         * This is the start of the velocity Verlet algorithm
         */
       
        // Print the initial position to file
        Particle3D.toVMD(ps, 0, pwTraj);
        
        // Print the time, initial energy, initial position and
        // adjusted total momentum to file
        pwCheck.println("T, E, P");
        pwCheck.printf("%.3f, %.3f, %.3f%n", t, totE, totP.mag());

        /*
         * The main loop
         */
        for (int i=1;i<numstep;i++){
            // Record current position
            Vector3D[] currentPos;
            currentPos = new Vector3D[ps.length];
            for(int j=0; j<ps.length; j++){
		currentPos[j] = ps[j].getPosition();
            }

            // Update the postion using current velocity and force
            Particle3D.manyPosition(ps, dt, g);
            
            // Update the velocity, based on average of current and
            // new forces
            Particle3D.manyVelocity(ps, currentFors, dt, g);
            
            // Calculate the energy and momentum of the system
            totE = Particle3D.manyEnergy(ps, g);
            totP = Particle3D.manyMomentum(ps);
            
            // Check for aphelion/perihelion and orbital position
            for(int j=0; j<ps.length; j++){
                // Define Moon relative to Earth
                if(j==4){
                    currentPos[j]=Vector3D.subVec(currentPos[j], currentPos[j-1]);
                    ps[j].setPosition(Particle3D.separation(ps[j], ps[j-1]));}
                
                // Update ap/peri-helion Matrix
                ap[j] = Particle3D.apPeri(ps[j], ap[j]);
                // Update orbit counter
                oCounts[j] += Particle3D.oDisp(currentPos[j], ps[j]);
                  
                // Redefine Moon relative to origin
                if(j==4){
                    ps[j].setPosition(Vector3D.addVec(ps[j].getPosition(), ps[j-1].getPosition()));
                }
            }           
            
            // Increase the time
            t += dt;
            
            // Print checks and position to file daily
            if (i%100==0){
                Particle3D.toVMD(ps, i, pwTraj);
                pwCheck.printf("%.1f, %.3f, %.3f,%n", t, totE, totP.mag());
            }
        }
        // Calculate the eccentricity of each orbit
        double[] es;
        es = new double[ps.length];
        for(int i=0;i<ps.length;i++){
            es[i] = Particle3D.eccentric(ap[i]);
        }

    	// Verify Kepler's Third Law for each orbit
    	double[] kp;
    	kp = new double[ps.length];
    	for(int i=0;i<ps.length;i++){
            // Define Moon Relative to Earth
            if(i==4){
                kp[i] = Particle3D.keplerThree(ps[3].getMass(), t, oCounts[i], g, ap[i]);
            }else{
                kp[i] = Particle3D.keplerThree(ps[0].getMass(), t, oCounts[i], g, ap[i]);
            }
        }

        //Print orbit characteristics to file
        pwOrbit.printf("Planet, No. Orbits, Aphelion, Perihelion, Eccentricity, KeplerIII%n");
        for(int i=0; i<c; i++){
            pwOrbit.printf("%s, %.1f, %.3f, %.3f, %.3f, %.3f %n", 
			   ps[i].getName(), oCounts[i], ap[i][0], ap[i][1], es[i], kp[i]);
        }
        pwOrbit.printf("NOTE: Moon orbit defined relative to Earth");
        
        //Close output files
        pwTraj.close();
        pwCheck.close();
        pwOrbit.close();
    }
}
