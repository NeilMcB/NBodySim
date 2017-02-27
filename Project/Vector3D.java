 /**
 * A class for 3D vectors. Includes constructors, setters and getters,
 * instance methods to calculate magnitude squared and magnitude, and 
 * static methods to implement simple arithmetic operations.
 *
 * @author M. Corrigan
 * @author N. McBlane
 * @version "10/2015"
 *
 */

public class Vector3D {
   
    /*
     * Vector Elements
     * Protected with "private" keyword
     */
    private double x;
    private double y;
    private double z;

    /*
     * Constructors
     */

    /** Defualt constructor. Constructs a new Vector3D with
     *  x, y and z components set to zero.
     */ 
    public Vector3D()
    {
	this.setX(0.0);
	this.setY(0.0);
	this.setZ(0.0);
    }

    /** Explicit constructor. Constructs a new Vector3D from given
     *  x, y and z components.
     *
     * @param xx a double giving the x component of the new Vector3D 
     * @param yy a double giving the y component of the new Vector3D
     * @param zz a double giving the z component of the new Vector3D
     */

    public Vector3D(double xx, double yy, double zz){
	this.setX(xx);
	this.setY(yy);
	this.setZ(zz);
    }

    /** Copy constructor. Constructs a new Vector3D by copying
     *  all components of another Vector3D.
     *
     * @param original the Vector3D to be copied
     */
    public Vector3D(Vector3D original){
        this.setX(original.getX());
	this.setY(original.getY());
	this.setZ(original.getZ());
    }
    /*
     * Getters and Setters
     */

    // Getters provide access to private internal variables
    
    /** Gets the x component of a Vector3D.
     *
     *  @return a double instance representing the Vector3D x component
     */
    public double getX() {return this.x;}
   
    /** Gets the y component of a Vector3D.
     *
     *  @return a double instance representing the Vector3D y component
     */
    public double getY() {return this.y;}
   
    /** Gets the z component of a Vector3D.
     *
     *  @return a double instance representing the Vector3D z component
     */
    public double getZ() {return this.z;}

    // Setters provide the internal variables

    /** Sets the x component of a Vector3D.
     *
     *  @param xx a double to set the x component
     */
    public void setX(double xx) {this.x = xx;}
    /** Sets the y component of a Vector3D.
     *
     *  @param yy a double to set the y component
     */
    public void setY(double yy) {this.y = yy;}
    /** Sets the z component of a Vector3D.
     *
     *  @param zz a double to set the z component
     */
    public void setZ(double zz) {this.z = zz;}

    /** Sets all three components of a Vector3D.
     *
     *  @param xx a double to set the x component
     *  @param yy a double to set the y component
     *  @param zz a double to set the z component
     */
    public void setVector3D(double xx, double yy, double zz){
	this.x = xx;
	this.y = yy;
	this.z = zz;
    }
    /*
     * Instance Methods
     *
     * These operate on a particular instance of an object.
     */
    
    /** Calculates the square magnitude A.A of the Vector3D A.
     *
     * @return a double representing the Vector3D square magnitude.
     */
    public double magSqrd(){
	return this.getX()*this.getX()+
	    this.getY()*this.getY()+
	    this.getZ()*this.getZ();
    }
    /** Calculates the magnitude |A| of the Vector3D A.
     *
     * @return a double representing the Vector3D magnitude.
     */
    public double mag(){
	return Math.sqrt(this.magSqrd());
    }
    
     /** Returns a String representation of the Vector3D. Methods 
     * called 'toString' are automagically called when an object
     * is printed.<br>
     * The output takes the sign of the component (i.e. positive/negative)
     * in to account.
     *
     * @return a string representation of the Vector3D instance
     */
    public String toString(){
	double a = this.getX();
	double b = this.getY();
	double c = this.getZ();

	return "(" + a + "," + b + "," + c + ")";
    }

    /** Takes a Vector3D and multiplies each of its components by
     *  a given scalar.
     *
     *  @param scal a double that multiplies each of the Vector3D's components
     *
     *  @return a reference to this Vector3D that represents a multiplication by scalar
     */
    public Vector3D scalarMult(double scal){
	Vector3D mult = new Vector3D(this.getX()*scal
				     ,this.getY()*scal
				     ,this.getZ()*scal);

	return mult;
    }

    /** Takes a Vector3D and divide each of its components by
     *  a given scalar.
     *
     *  @param scal a double that divides each of the Vector3D's components
     *
     *  @return a reference to the this Vector3D that represents a division by scalar 
     */
    public Vector3D scalarDiv(double scal){
	Vector3D div = new Vector3D(this.getX()/scal
				     ,this.getY()/scal
				     ,this.getZ()/scal);

	return div;
    }
    
    
    /*
    * Static Methods
    *
    * These allow us to build compound statements.
    */
    
    /** Takes two Vector3Ds and adds their corresponding components
     *  then returns a new Vector3D
     *
     *  @param a a Vector3D to which we add components
     *  @param b a Vector3D from which we take components to add
     *  
     *  @return a Vector3D which has the added componens a+b
     */
    public static Vector3D addVec(Vector3D a, Vector3D b){
	return new Vector3D(a.getX()+b.getX(), a.getY()+b.getY(), a.getZ()+b.getZ());
    }
    
    /** Takes two Vector3Ds and subtracts their corresponding components
     *  then returns a new Vector3D
     *
     *  @param a a Vector3D from which we subtract components
     *  @param b a Vector3D from which we take components to subtract
     *  
     *  @return a Vector3D which has the added componens a+b
     */
    public static Vector3D subVec(Vector3D a, Vector3D b){
	return new Vector3D(a.getX()-b.getX(), a.getY()-b.getY(), a.getZ()-b.getZ());
    }
    
     /** Takes two Vector3Ds and calculates the dot product by
     *  multiplying each Vector3D's corresponding components
     *
     *  @param a a Vector3D on which we perform the dot product
     *  @param b a Vector3D on which we perform the dot product
     *
     *  @return a double which represents the dot product a.b
     */
     public static double dotVec(Vector3D a, Vector3D b){
        return a.getX()*b.getX() +
               a.getY()*b.getY() +
               a.getZ()*b.getZ();
    }

    /** Takes two Vector3Ds and performs the vector cross product
     *  with their extracted components. Returns a new Vector3D
     *  containing this information.
     *
     *  @param a the first Vector3D in the expression
     *  @param b the second Vector3D in the expression
     *
     *  @return a Vector3D which represents the vector cross product aXb
     */ 
    public static Vector3D crossVec(Vector3D a, Vector3D b){
	return new Vector3D((a.getY()*b.getZ())-(a.getZ()*b.getY()),
		     a.getZ()*b.getX()-a.getX()*b.getZ(),
		     a.getX()*b.getY()-a.getY()*b.getX());
    }

}
