ArrayList<PVector> point_array = new ArrayList<PVector>();
float region_width = 400.0;
boolean is_tracking_camera = true;

void setup() {
    size(800,800,P3D);
    //generate_torus(20000);
    generate_torus_sector(20000, -PI/2.0,PI);
    //generate_torus_sector(20000, 0.0, 1.0*PI/2.0);
    /*
    **
    **  Currently not working for phi>PI --> From atan2 implementation
    **
    */
}

void draw() {
    background(40,40,40);

    if(is_tracking_camera) {
        track_camera(); }

    display_point_array();
}

void track_camera() {
    float eye_factor = 800.0;
    float theta = float(mouseY)/float(height)*(1.0*PI);
    float phi   = float(mouseX)/float(width)*(2.0*PI);
    camera(
        eye_factor*cos(phi)*sin(theta), eye_factor*sin(phi)*sin(theta), eye_factor*cos(theta),    //Eye pos
        0.0, 0.0, 0.0,      //Centre pos
        0.0, 0.0, 1.0       //Up direction
    );
}

void display_point_array() {
    /*
    . Displays the point_array particles
    */
    float point_size = 3;
    PVector point_col = new PVector(200,100,100);
    for(int i=0; i<point_array.size(); i++) {
        PVector point = point_array.get(i);
        pushMatrix();
        translate(point.x, point.y, point.z);
        pushStyle();
        noStroke();
        fill(point_col.x, point_col.y, point_col.z);
        sphere(point_size);
        popStyle();
        popMatrix();
    }
}

void generate_torus(int N_points) {
    /*
    . Generates a torus through a random point cloud being sampled
    */
    float centre_R = 100.0;     //Radius to centre of tube
    float tube_R   = 50.0;      //Radius of the tube
    for(int i=0; i<N_points; i++) {
        PVector point = new PVector( 
            region_width*random(-0.5, 0.5), 
            region_width*random(-0.5, 0.5), 
            region_width*random(-0.5, 0.5) 
        );
        boolean withinShape = pow( centre_R-sqrt( pow(point.x,2) + pow(point.y,2) ) ,2) +pow(point.z,2) <= pow(tube_R,2);
        if(withinShape) {
            point_array.add(point);
        }
    }
}
void generate_torus_sector(int N_points, float phi_min, float phi_max) {
    /*
    . Generates a torus sector through a random point cloud being sampled
    */
    float centre_R = 100.0;     //Radius to centre of tube
    float tube_R   = 50.0;      //Radius of the tube

    for(int i=0; i<N_points; i++) {
        PVector point = new PVector( 
            region_width*random(-0.5, 0.5), 
            region_width*random(-0.5, 0.5), 
            region_width*random(-0.5, 0.5) 
        );
        float phi = atan2(point.y, point.x);
        boolean withinShape  = pow( centre_R-sqrt( pow(point.x,2) + pow(point.y,2) ) ,2) +pow(point.z,2) <= pow(tube_R,2);
        boolean withinSector = (phi_min < phi) && (phi < phi_max);
        if(withinShape && withinSector) {
            point_array.add(point);
        }
    }
}

void keyPressed() {
    if(key == '1') {
        is_tracking_camera = !is_tracking_camera;
    }
}
