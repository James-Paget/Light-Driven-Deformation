PVector gridDim = new PVector(60, 60, 60);
PVector gridDim_coords = new PVector(800, 800, 800);
ArrayList<grid> grid_set = new ArrayList<grid>();

float zoom_value = 0.5;
boolean zoomIn = false;
boolean zoomOut = false;

boolean trackCamera = true;

void setup() {
    size(800,800, P3D);

    generateFullGrid();
    generateRandomGrid();
}

void draw() {
    updateZoom();
    calculateCameraTrack(zoom_value, trackCamera);
    drawBackground();

    pushMatrix();
    //This translate moves the axes origin the actual centre of the grid, not for the corner or bounding box -> preference which you prefer
    translate(
        gridDim_coords.x/2.0,
        gridDim_coords.y/2.0,
        gridDim_coords.z/2.0
    );
    drawAxes();
    popMatrix();
    
    drawFullGrid();
    drawFullVectors();
}

class grid {
    boolean active = false;
    boolean activeVector = false;
    PVector gridVector = new PVector(0,0,0);

    grid() {
        generateRandomVector();
    }
    
    void generateRandomVector() {
        gridVector.x = random(-1.0, 1.0);
        gridVector.y = random(-1.0, 1.0);
        gridVector.z = random(-1.0, 1.0);
        float mag = sqrt( pow(gridVector.x, 2) + pow(gridVector.y, 2) + pow(gridVector.z, 2) );
        gridVector.x /= mag;
        gridVector.y /= mag;
        gridVector.z /= mag;
    }
}

void updateZoom() {
    if(zoomIn) {
        zoom_value -= 0.01;}
    if(zoomOut) {
        zoom_value += 0.01;}
}
void calculateCameraTrack(float zoom, boolean trackCamera) {
    if(trackCamera) {
        float phi = 2.0*PI*(float(mouseX)/gridDim_coords.x);   //0->2*PI from 0->width
        float theta = PI*(float(mouseY)/gridDim_coords.y); //0->PI from 0->height
        //println("phi= ",phi,",   theta= ",theta);
        camera(
            gridDim_coords.x/2.0 +zoom*1.5*gridDim_coords.x*cos(phi)*sin(theta), gridDim_coords.y/2.0 +zoom*1.5*gridDim_coords.y*sin(phi)*sin(theta), gridDim_coords.z/2.0 +zoom*1.5*gridDim_coords.z*cos(theta),    //Eye pos
            gridDim_coords.x/2.0, gridDim_coords.y/2.0, gridDim_coords.z/2.0,     //Centre pos
            0.0, 0.0, 1.0   //Up direction
        );
    }
}
void drawAxes() {
    pushStyle();
    noFill();
    strokeWeight(10.0);
    //X axis
    stroke(255,0,0);
    line(
        0.0, 0.0, 0.0,
        1.2*gridDim_coords.x, 0.0, 0.0
    );
    //Y axis
    stroke(0,255,0);
    line(
        0.0, 0.0, 0.0,
        0.0, 1.2*gridDim_coords.y, 0.0
    );
    //Z axis
    stroke(0,0,255);
    line(
        0.0, 0.0, 0.0,
        0.0, 0.0, 1.2*gridDim_coords.z
    );
    popStyle();
}

void drawBackground() {
    background(255);
}
void generateFullGrid() {
    for(int k=0; k<gridDim.z; k++) {
        for(int j=0; j<gridDim.y; j++) {
            for(int i=0; i<gridDim.x; i++) {
                grid_set.add( new grid() );
            }
        }
    }
}
void generateRandomGrid() {
    for(int i=0; i<grid_set.size(); i++) {
        int x = floor(i % gridDim.x);
        int y = floor(floor(i / gridDim.x) % (gridDim.y));
        int z = floor(i / (gridDim.x*gridDim.y));
        //Sphere
        //boolean cond = ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) + pow(z-gridDim.z/2.0, 2) ) < 0.9*gridDim.x;
        //Spheroid
        //float factor_spheroid = 5.0;
        //boolean cond = ( ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) )/factor_spheroid  + pow(z-gridDim.z/2.0, 2)) < 0.9*gridDim.x;
        //Torus
        //float c_factor = 18.0;
        //float rad_factor = 4.0;
        //boolean cond = pow( c_factor - sqrt(pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2)), 2) + pow(z-gridDim.z/2.0, 2) < pow(rad_factor, 2);
        //Spheroid + Torus
        //float factor_spheroid = 5.0;
        //boolean cond_spheroid = ( ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) )/factor_spheroid  + pow(z-gridDim.z/2.0, 2)) < 0.9*gridDim.x;
        //float c_factor = 10.0;
        //float rad_factor = 3.0;
        //boolean cond_torus = pow( c_factor - sqrt(pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2)), 2) + pow(z+gridDim.z/11.0 -gridDim.z/2.0, 2) < pow(rad_factor, 2);
        //boolean cond = cond_spheroid || cond_torus;
        //Multi sphere ring
        //boolean cond = ringOfSpheresCond(new PVector(x,y,z), 10, 0.3*gridDim.x, 0.3*gridDim.x);
        //Sphere & Sphere
        //boolean cond_s1 = ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) + pow(z-gridDim.z/2.0, 2) ) < 1.4*gridDim.x;
        //boolean cond_s2 = ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) + pow(z +gridDim.z/6.0 -gridDim.z/2.0, 2) ) < 0.3*gridDim.x;
        //boolean cond = cond_s1 || cond_s2;
        //Sphere & Anti-Sphere
        //boolean cond_s1 = ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) + pow(z-gridDim.z/2.0, 2) ) < 1.4*gridDim.x;
        //boolean cond_s2 = ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) + pow(z +gridDim.z/6.0 -gridDim.z/2.0, 2) ) < 0.5*gridDim.x;
        //boolean cond = cond_s1 && !cond_s2;
        //Sphere Tri Setup
        boolean cond_ring = ringOfSpheresCond(new PVector(x,y,z), 3, 0.4*gridDim.x, 0.15*gridDim.x);
        boolean cond_other = ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) + pow(z-gridDim.z/2.0 +0.15*gridDim.x, 2) ) < 0.4*gridDim.x;
        boolean cond = cond_ring || cond_other;
        //Random
        //boolean cond = (random(0.0, 1.0) < 0.2);  //## Bug Fixing ##
        if(cond) {
            grid_set.get(i).active = true;
            //grid_set.get(i).activeVector = true;
        }
    }
}

boolean ringOfSpheresCond(PVector xyz, int nSpheres, float radi, float dist) {
    float theta_spacing = 2.0*PI/nSpheres;
    for(float theta=0.0; theta<2.0*PI; theta+=theta_spacing) {
        PVector offset = new PVector(dist*cos(theta), dist*sin(theta), 0.0);
        boolean cond = ( pow(xyz.x +offset.x -gridDim.x/2.0, 2) + pow(xyz.y +offset.y -gridDim.y/2.0, 2) + pow(xyz.z-gridDim.z/2.0, 2) ) < radi;
        if(cond) {
            return true;
        }
    }
    return false;
}

void drawFullGrid() {
    for(int k=0; k<gridDim.z; k++) {
        for(int j=0; j<gridDim.y; j++) {
            for(int i=0; i<gridDim.x; i++) {
                if(grid_set.get(int(k*(gridDim.x*gridDim.y) +j*gridDim.x +i)).active) {
                    drawGrid(new PVector(i*(gridDim_coords.x/gridDim.x), j*(gridDim_coords.y/gridDim.y), k*(gridDim_coords.z/gridDim.z)), grid_set.get(int(k*(gridDim.x*gridDim.y) +j*gridDim.x +i)));
                }
            }
        }
    }
}
void drawGrid(PVector pos, grid cGrid) {
    pushMatrix();
    translate(pos.x, pos.y, pos.z);
    pushStyle();
    //fill(200,120,120);      //Fixed Color
    float reduced_gridPos = pos.z;// -gridDim_coords.z/2.0;
    float reduced_gridHeight = 0.65*gridDim_coords.z;   //Used to tune color gradient
    //fill(255*pow(cos( (PI/2.0)*(reduced_gridPos/reduced_gridHeight) ), 1), 80, 255*(reduced_gridPos/reduced_gridHeight));      //Z Color        //pow(sin( (PI/2.0)*(reduced_gridPos/reduced_gridHeight) ), 2)
    //fill(255*cos( (PI/2.0)*(reduced_gridPos/reduced_gridHeight) ), 120, 255*(reduced_gridPos/reduced_gridHeight));
    fill(255*cos( (PI/2.0)*(pos.z/gridDim_coords.z) ), 100, 255*(pos.z/gridDim_coords.z));
    box(gridDim_coords.x/gridDim.x, gridDim_coords.y/gridDim.y, gridDim_coords.z/gridDim.z);
    popStyle();
    popMatrix();
}
void drawFullVectors() {
    for(int k=0; k<gridDim.z; k++) {
        for(int j=0; j<gridDim.y; j++) {
            for(int i=0; i<gridDim.x; i++) {
                if(grid_set.get(int(j*gridDim.x +i)).activeVector) {
                    drawVector(
                        new PVector(i*(gridDim_coords.x/gridDim.x), j*(gridDim_coords.y/gridDim.y), k*(gridDim_coords.z/gridDim.z)), 
                        grid_set.get(int(k*(gridDim.x*gridDim.y) +j*gridDim.x +i))
                    );
                }
            }
        }
    }
}
void drawVector(PVector pos, grid cGrid) {
    /*
    pushStyle();
    fill(255);
    stroke(255);

    pushMatrix();
    translate(
        pos.x +(gridDim_coords.x/gridDim.x)/2.0,
        pos.y +(gridDim_coords.y/gridDim.y)/2.0,
        pos.z +(gridDim_coords.z/gridDim.z)/2.0
    );
    ellipse(
        (gridDim_coords.x/gridDim.x)*0.2,
        (gridDim_coords.x/gridDim.x)*0.2
    );
    popMatrix();

    pushMatrix();
    translate(
        pos.x +(gridDim_coords.x/gridDim.x)/2.0,
        pos.y +(gridDim_coords.y/gridDim.y)/2.0,
        pos.z +(gridDim_coords.z/gridDim.z)/2.0
    );
    line(pos.x +(gridDim_coords.x/gridDim.x)/2.0, 
        pos.y +(gridDim_coords.y/gridDim.y)/2.0, 
        pos.x +(gridDim_coords.x/gridDim.x)/2.0 +cGrid.gridVector.x*(gridDim_coords.x/gridDim.x)*0.4, 
        pos.y +(gridDim_coords.y/gridDim.y)/2.0 +cGrid.gridVector.y*(gridDim_coords.y/gridDim.y)*0.4
    );
    popMatrix();

    popStyle();
    */
}

void keyPressed() {
    if(key == '1') {
        int index = int(floor(gridDim.x*(mouseX/gridDim_coords.x)) + gridDim.x*floor(gridDim.y*(mouseY/gridDim_coords.y)));
        grid_set.get(index).active = !grid_set.get(index).active;
    }
    if(key == '2') {
        int index = int(floor(gridDim.x*(mouseX/gridDim_coords.x)) + gridDim.x*floor(gridDim.y*(mouseY/gridDim_coords.y)));
        grid_set.get(index).activeVector = !grid_set.get(index).activeVector;
    }
    if(key == '3') {
        for(int i=0; i<grid_set.size(); i++) {
            int x = floor(i % gridDim.x);
            int y = floor(i / gridDim.x);
            int x_mouse = floor((mouseX/gridDim_coords.x)*gridDim.x);
            int y_mouse = floor((mouseY/gridDim_coords.y)*gridDim.y);
            boolean cond = ( pow(x-x_mouse, 2) + pow(y-y_mouse/2.0, 2) ) < 5;
            if(cond) {
                grid_set.get(i).active = true;
                //grid_set.get(i).activeVector = true;
            }
        }
    }
    if(key == '4') {
        trackCamera = !trackCamera;
    }

    if(key == 'w') {
        zoomIn = true;
    }
    if(key == 's') {
        zoomOut = true;
    }
}
void keyReleased() {
    if(key == 'w') {
        zoomIn = false;
    }
    if(key == 's') {
        zoomOut = false;
    }
}