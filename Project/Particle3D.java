s
/**
 * Computer Modelling, Exercise 3: Particle.3D class.
 *
 * @author N. McBlane
 * @author M. Corrigan
 * @version "11/2015"
 *
 */
import java.io.*; 
import java.util.Scanner;
 
public class Particle3D {

    /* ******************************************
     * Properties
     ********************************************/
    private double mass;
    private Vector3D position = new Vector3D();
    private Vector3D velocity = new Vector3D();
    private String name;

    // Setters and Getters
    
    /** Get the position of a particle.
     *
     * @return a Vector3D representing the position.
     */
    public Vector3D getPosition() { return position; }

    /** Get the velocity of a particle.
     *
     * @return a Vector3D representing the velocity.
     */
    public Vector3D getVelocity() { return velocity; }

    /** Get the mass of a particle.
     *
     * @return a double representing the mass.
     */
    public double getMass()     { return mass; }

    /** Get the name of a particle.
     *
     * @return a String representing the name.
     */
    public String getName() {return name;}

    /** Set the position of a particle.
     *
     * @param p a Vector3D representing the position.
     */
    public void setPosition(Vector3D p) 
    { 
	position.setX(p.getX());
	position.setY(p.getY());
	position.setZ(p.getZ()); 
    }
    
    /** Set the velocity of a particle.
     *
     * @param v a Vector3D representing the velocity.
     */
    public void setVelocity(Vector3D v)
    { 
	velocity.setX(v.getX());
	velocity.setY(v.getY());
	velocity.setZ(v.getZ()); 
    }
       
    /** Set the mass of a particle.
     *
     * @param m a double representing the mass.
     */
    public void setMass(double m)     { this.mass = m; }
    
    /** Set the name of a particle.
     *
     *
     * @param n a string representing the name
     */
    public void setName(String n)     {this.name = n; }

    /* ******************************************
     * Constructors
     ********************************************/
    
    /** Default constructor. Sets mass to zero 
     * and leaves name blank. Places the particle
     * at the origin with no velocity.
     */
    public Particle3D() {
        this.setMass(0.);
        this.setPosition(new Vector3D());
        this.setVelocity(new Vector3D());
        this.setName("");
    }

    /** Explicit constructor. Constructs a new Particle3D with
     * explicitly given position, velocity, name and mass.
     *
     * @param m a double that defines the mass.
     * @param p a Vector3D that defines the position.
     * @param v a Vector3D that defines the velocity.
     * @param n a string that defines the name of the particle.
     */
    public Particle3D(double m, Vector3D p, Vector3D v, String n) {
        this.setMass(m);
        this.setPosition(p);
        this.setVelocity(v);
        this.setName(n);
    }

    /** Scanner Constructor. Constructs a new Particle3D with
     * postition, velocity, name and mass taken from a scanner.
     *
     * @param vecIn a scanner that contains information representing a particle
     */

    public Particle3D(Scanner vecIn) {
        this.setName(vecIn.next());
        this.setPosition(new Vector3D(vecIn.nextDouble(),vecIn.nextDouble(),vecIn.nextDouble()));
	this.setVelocity(new Vector3D(vecIn.nextDouble(),vecIn.nextDouble(),vecIn.nextDouble()));
        this.setMass(vecIn.nextDouble());
    }
    /* ******************************************
     * toString Method
     ********************************************/
    
    /** Returns a String representation of Particle3D.
     * Used to print a Particle3D instance using the "%s"
     * format identifier.
     *
     * @return a string representing the particle.
     */
    public String toString() {
        double x = this.getPosition().getX();
        double y = this.getPosition().getY();
        double z = this.getPosition().getZ();
	String n = this.getName();
	return n + " " + x + " " + y + " " + z;
    }



    /* ******************************************
     * Instance Methods
     ********************************************/
    
    /** Returns the kinetic energy of a Particle3D,
     * calculated as 1/2*m*v^2.
     *
     * @return a double that is the kinetic energy.
     */
    public double kineticEnergy() { 
	return 0.5*mass*velocity.mag()*velocity.mag(); 
    }

    /** Time integration support: evolve the velocity
     * according to dv = f/m * dt.
     *
     * @param dt a double that is the timestep.
     * @param force a double that is the current force on the particle.
     */
    public void leapVelocity(double dt, Vector3D force) {
        velocity = Vector3D.addVec(velocity,
				   force.scalarMult(dt/mass));
    }
    
    /** Time integration support: evolve the position
     * according to dx = v * dt.
     *
     * @param dt a double that is the timestep.
     */
    public void leapPosition(double dt) {
        position = Vector3D.addVec(position, velocity.scalarMult(dt));
    }

    /** Time integration support: evolve the position
     * according to dx = v * dt + 0.5 * a * dt**2.
     *
     * @param dt a double that is the timestep.
     * @param force a double that is the current force.
     */
    public void leapPosition(double dt, Vector3D force) {
        position = Vector3D.addVec(position, 
				   Vector3D.addVec(velocity.scalarMult(dt), 
						   force.scalarMult(0.5*dt*dt/mass)));
    }
     
    /* ******************************************
     * Static Methods
     ********************************************/

    /** Returns the vector separating two particles.
     *
     * @param a a particle.
     * @param b another particle.
     *
     * @return a Vector3D that represents the displacement
     * of particle b from a.
     **/
    public static Vector3D separation(Particle3D a, Particle3D b) {
        return Vector3D.subVec(a.getPosition(), b.getPosition());
    }

    /** Returns the gravitational energy between a 
     * pair of Particle3Ds.
     *
     * @param a a Particle3D.
     * @param b a Particle3D.
     * @param g a double representing the gravitational constant.
     *
     * @return a double representing the pair potential of a and b.
     *
     */
    public static double pairEnergy(Particle3D a, Particle3D b, double g){
        double aMass = a.getMass();
        double bMass = b.getMass();
        double magrAB = -g*(Particle3D.separation(a,b)).mag();
        return -aMass*bMass/magrAB;
    }

    /** Returns the gravitational force between a pair
     * pair of Particle3Ds.
     *
     * @param a a Particle3D.
     * @param b a Particle3D.
     * @param g a double representing the gravitational constant.
     *
     * @return a Vector3D representing the force on a due to b.
     *
     */
    public static Vector3D gravForce(Particle3D a, Particle3D b, double g){
        double aMass = a.getMass();
        double bMass = b.getMass();
        Vector3D r = separation(a,b);
        double magrAB = r.mag();
        double gravScalar = -g*(aMass*bMass)/(Math.pow(magrAB, 3.));
        return r.scalarMult(gravScalar);
    }


    /* ******************************************
     * Array Methods
     ********************************************/


    /** Prints VMD trajectory format to a specified
     *  printwriter.
     *
     * @param ps an array of particles.
     * @param i an integer specifying the current step of the simulation.
     * @param pw the trajectory file printwriter.
     *
     **/
    public static void toVMD(Particle3D[] ps, int i, PrintWriter pw ){
        pw.printf("%d \nPoint = %d\n", ps.length, i);
        for(int j=0; j<ps.length; j++){
            pw.println(ps[j].toString());
        }
    }

    /** Calculates the gravitational force between
     *  an array of Particle3Ds.
     *
     * @param ps an array of particles.
     * @param g a double representing the gravitational constant.
     *
     * @return forces a Vector3D array representing the force on each particle.
     *
     **/
    public static Vector3D[] manyForce(Particle3D[] ps, double g){
        Vector3D[] forces;
        forces = new Vector3D[ps.length];
        //initialize array of zero vectors
        for (int i=0; i<ps.length; i++) {
            forces[i] = new Vector3D();
        }
        for (int i=0; i<ps.length; i++){
            //Calculate total force due to all other particles
            for (int j=(i+1); j<ps.length; j++){
                forces[i] = Vector3D.addVec(forces[i], gravForce(ps[i],ps[j], g));
                // By Newton III
                forces[j] = Vector3D.addVec(forces[j], gravForce(ps[i],ps[j], g).scalarMult(-1.));
	    }
        }
        return forces;
    }

    /** Calculates the total energy of an array
     *  of Particle3Ds.
     *
     * @param ps an array of particles.
     * @param g a double representing the gravitational constant.
     *
     * @return energy a double representing the total energy.
     *
     **/    
    public static double manyEnergy(Particle3D[] ps, double g){
        double energy = 0.0;
        for (int i=0; i<ps.length; i++){
            energy += ps[i].kineticEnergy();
            for (int j=i+1; j<ps.length; j++){
                energy += pairEnergy(ps[i],ps[j], g);
            }
        }
        return energy;
    }

    /** Calculates the total momentum of an array
     *  of Particle3Ds.
     *
     * @param ps an array of particles.
     *
     * @return totP a Vector3D representing the total momentum.
     *
     **/
    public static Vector3D manyMomentum(Particle3D[] ps){
        Vector3D totP = new Vector3D();
        for (int i=0;i<ps.length;i++){
            totP = Vector3D.addVec((ps[i].getVelocity()).scalarMult(ps[i].getMass()),totP);
        }
        return totP;
    }

    /** Corrects the total momentum of an array
     *  of Particle3Ds to zero.
     *
     * @param ps an array of particles.
     *
     * @return ps and array of particles with corrected total momentum.
     *
     **/
    public static Particle3D[] correctMomentum(Particle3D[] ps){
        double totM = 0.;
        for (int i=0;i<ps.length;i++){
            totM += ps[i].getMass();
        }
        Vector3D totP = new Vector3D();
        totP = manyMomentum(ps);
        Vector3D comV = new Vector3D();
        // V = P/M
        comV = totP.scalarMult(1/totM);
        // Adjust veloctiy of each particle
        for (int i=0;i<ps.length;i++) {
            ps[i].setVelocity(Vector3D.addVec(ps[i].getVelocity(),comV.scalarMult(-1.)));
        }
        return ps;
    }

    /** Compares the current position to the aphelion and 
     *  perihelion of an orbit.
     *
     * @param p a Particle3D.
     * @param ap a array containting the particle ap- and perihelion.
     *
     * @return ap an array containing the updated ap- and perihelion.
     *
     **/
    public static double[] apPeri(Particle3D p, double[] ap){
        //aphelion
        if ((p.getPosition()).mag() > ap[0]){
	    ap[0] = (p.getPosition()).mag();}
        //perihelion
        if ((p.getPosition()).mag() < ap[1]){
	    ap[1] = (p.getPosition()).mag();}
        return ap;
    }

    /** Calculate the angular displacement between a Particle3D
     *  and its previous position.
     *
     * @param prevPos a Vector3D representing the Particle3D's 
     * position in the previous timestep.
     * @param p a Particle3D.
     *
     * @return oAngle a double representing the fraction of 
     * an orbit completed in the step.
     *
     **/
    public static double oDisp(Vector3D prevPos, Particle3D p){
        double oDot = Vector3D.dotVec(prevPos, p.getPosition());
        double oAngle = Math.acos(oDot/(prevPos.mag()*((p.getPosition()).mag())));
        return oAngle/(2*Math.PI);
    }

    /** Calculate the eccentricity of a Particle3D orbit.
     *
     * @param ap an array containing an orbit ap- and perihelion.
     *
     * @return e a double representing orbit eccentricity.
     *
     **/
    public static double eccentric(double[] ap){
        return (ap[0]-ap[1])/(ap[0]+ap[1]);
    }

    /** Verify Kepler's Third law for a Particle3D orbit.
     *
     * @param m a double representing central body mass.
     * @param t a double representing total simulation time.
     * @param orbits a double representing number of Particle3D orbits.
     * @param g a double representing the gravitational constant.
     * @param ap an array containing orbit ap- and perihelion.
     *
     * @return kp a double representing the agreement with Kepler III; 1 is a perfect agreement.
     *
     **/      
    public static double keplerThree(double m,  double t, double orbits, double g, double[] ap){
	return (4*Math.pow(Math.PI,2)   
                *  Math.pow(((ap[0]+ap[1])/2),3))    
	    /  (g*m*(Math.pow((t/orbits),2)));
    }

    /** Time integration support: evolve the velocity
     * according to dv = f/m * dt for an array.
     *
     * @param ps an array of particles.
     * @param forces an array that is the current forces between the particles.
     * @param dt a double that is the timestep.
     * @param g a double that represents the gravitational constant.
     *
     * @return ps an array of particles with updated velocity.
     *
     */
    public static Particle3D[] manyVelocity(Particle3D[] ps, Vector3D[] forces, double dt, double g){
        // Calculate force due to all other particles
        Vector3D[] new_forces;
        new_forces = manyForce(ps, g);
        for (int i=0; i<ps.length; i++){
            //Update the velocity using the average of the current and new forces
            ps[i].leapVelocity(dt, (Vector3D.addVec(new_forces[i], forces[i])).scalarMult(0.5));
            //Update the current force to the new force
            forces[i] = new_forces[i];
        }
        return ps;
    }   

    /** Time integration support: evolve the position
     * according to dx = v * dt + 0.5 * a * dt**2 for an array.
     *
     * @param ps an array of particles.
     * @param dt a double that is the timestep.
     * @param g a double that represents the gravitational constant.
     *
     * @return ps an array of particles with updated position.
     *
     */
    public static Particle3D[] manyPosition(Particle3D[] ps, double dt, double g){
        Vector3D[] forces;
        forces = manyForce(ps, g);
        for (int i=0; i<ps.length; i++){
            ps[i].leapPosition(dt,forces[i]);
        }
        return ps;
    }
};
