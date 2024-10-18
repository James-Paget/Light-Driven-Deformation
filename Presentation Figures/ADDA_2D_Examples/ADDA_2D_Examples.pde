PVector gridDim = new PVector(20, 20);
PVector gridDim_coords = new PVector(800, 800);
ArrayList<grid> grid_set = new ArrayList<grid>();

void setup() {
    size(800,800);

    generateFullGrid();
    generateRandomGrid();
}

void draw() {
    drawBackground();
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
        gridVector.z = 0.0;
        float mag = sqrt( pow(gridVector.x, 2) + pow(gridVector.y, 2) + pow(gridVector.z, 2) );
        gridVector.x /= mag;
        gridVector.y /= mag;
        gridVector.z /= mag;
    }
}

void drawBackground() {
    background(255,255,255);
}
void generateFullGrid() {
    for(int j=0; j<gridDim.y; j++) {
        for(int i=0; i<gridDim.x; i++) {
            grid_set.add( new grid() );
        }
    }
}
void generateRandomGrid() {
    for(int i=0; i<grid_set.size(); i++) {
        int x = floor(i % gridDim.x);
        int y = floor(i / gridDim.x);
        boolean cond = ( pow(x-gridDim.x/2.0, 2) + pow(y-gridDim.y/2.0, 2) ) < 20;
        if(cond) {
            grid_set.get(i).active = true;
            //grid_set.get(i).activeVector = true;
        }
    }
}
void drawFullGrid() {
    for(int j=0; j<gridDim.y; j++) {
        for(int i=0; i<gridDim.x; i++) {
            if(grid_set.get(int(j*gridDim.x +i)).active) {
                drawGrid(new PVector(i*(gridDim_coords.x/gridDim.x), j*gridDim_coords.y/gridDim.y), grid_set.get(int(j*gridDim.x +i)));
            }
        }
    }
}
void drawGrid(PVector pos, grid cGrid) {
    pushStyle();
    fill(200,100,100);
    rect(pos.x, pos.y, gridDim_coords.x/gridDim.x, gridDim_coords.y/gridDim.y);
    popStyle();
}
void drawFullVectors() {
    for(int j=0; j<gridDim.y; j++) {
        for(int i=0; i<gridDim.x; i++) {
            if(grid_set.get(int(j*gridDim.x +i)).activeVector) {
                drawVector(new PVector(i*(gridDim_coords.x/gridDim.x), j*gridDim_coords.y/gridDim.y), grid_set.get(int(j*gridDim.x +i)));
            }
        }
    }
}
void drawVector(PVector pos, grid cGrid) {
    pushStyle();
    fill(255);
    stroke(255);
    ellipse(
        pos.x +(gridDim_coords.x/gridDim.x)/2.0, 
        pos.y +(gridDim_coords.y/gridDim.y)/2.0,
        (gridDim_coords.x/gridDim.x)*0.2,
        (gridDim_coords.x/gridDim.x)*0.2
    );
    line(pos.x +(gridDim_coords.x/gridDim.x)/2.0, 
        pos.y +(gridDim_coords.y/gridDim.y)/2.0, 
        pos.x +(gridDim_coords.x/gridDim.x)/2.0 +cGrid.gridVector.x*(gridDim_coords.x/gridDim.x)*0.4, 
        pos.y +(gridDim_coords.y/gridDim.y)/2.0 +cGrid.gridVector.y*(gridDim_coords.y/gridDim.y)*0.4
    );
    popStyle();
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
}