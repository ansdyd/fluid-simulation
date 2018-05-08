float dt = 0.01;
int n = 128;
int gridSize = n + 2;
int pixelSize = 4;

int JACOBI_ITERATION = 10;
float FORCE_RADIUS = n / 25.0;

float[][] u = new float[gridSize][gridSize];
float[][] v = new float[gridSize][gridSize];
float[][] uPrev = new float[gridSize][gridSize];
float[][] vPrev = new float[gridSize][gridSize];
float[][] d = new float[gridSize][gridSize];

float[][] currentPressure = new float[gridSize][gridSize];
float[][] prevPressure = new float[gridSize][gridSize];

color[][] prevColor = new color[gridSize][gridSize];
color[][] currentColor = new color[gridSize][gridSize];

int previousMouseX;
int previousMouseY;

void setup() {
  background(0);
  noStroke();
  size(gridSize * pixelSize, gridSize * pixelSize);

  initColor();
  initVelocity();
  initPressure();
}

void draw() {
  // Following the instructions given in the jos stam paper.

  // TODO: Add forces.
  addForces();

  // Advect u.
  advectVelocity();

  // Calculate divergence of advected_u.
  getDivergence();

  // Calculate pressure based on divergence, use jacobi iteration.
  getPressure();

  // Subtract pressure from advected_u to create divergence free vector field.
  subtractPressureGradient();

  // Update color.
  swapColor();
  updateColor();

  // Draw color.
  drawColor();

  swapVelocity();
}

void subtractPressureGradient() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      int up = (i + 1) % gridSize;
      int down = (i - 1 + gridSize) % gridSize;
      int left = (j - 1 + gridSize) % gridSize;
      int right = (j + 1) % gridSize;

      u[i][j] -= (dt / 2.0) * (currentPressure[up][j] - currentPressure[down][j]);
      v[i][j] -= (dt / 2.0) * (currentPressure[i][right] - currentPressure[i][left]);
    }
  }
}

void getPressure() {
  for (int k = 0; k < JACOBI_ITERATION; k++) {
    for (int i = 0; i < gridSize; i++) {
      for (int j = 0; j < gridSize; j++) {
        int up = (i + 2) % gridSize;
        int down = (i - 2 + gridSize) % gridSize;
        int left = (j - 2 + gridSize) % gridSize;
        int right = (j + 2) % gridSize;
        currentPressure[i][j] = (d[i][j] +
        prevPressure[up][j] + prevPressure[down][j] + prevPressure[i][right] +
        prevPressure[i][left]) / 4.0;
      }
    }
    // TODO: swap pressure
    swapPressure();
  }
}

void swapPressure() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      prevPressure[i][j] = currentPressure[i][j];
    }
  }
}

void getDivergence() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      int right = (i + 1) % gridSize;
      int up = (j + 1) % gridSize;
      int left = (i - 1 + gridSize) % gridSize;
      int down = (j - 1 + gridSize) % gridSize;
      d[i][j] = (-2.0 / dt) * (u[right][j] - u[left][j] + v[i][up] - v[i][down]);
    }
  }
}

void swapVelocity() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      uPrev[i][j] = u[i][j];
      vPrev[i][j] = v[i][j];
    }
  }
}

void swapColor() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      prevColor[i][j] = currentColor[i][j];
    }
  }
}

void updateColor() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      float uu = u[i][j];
      float vv = v[i][j];

      float du = uu * dt;
      float dv = vv * dt;

      float x = ((i - du) + gridSize * gridSize) % gridSize;
      float y = ((j - dv) + gridSize * gridSize) % gridSize;

      int x1 = floor(x);
      int y1 = floor(y);

      int x2 = (x1 + 1) % gridSize;
      int y2 = (y1 + 1) % gridSize;

      int x2prime = x1 + 1;
      int y2prime = y1 + 1;

      float rInterpolated = abs((y2prime - y) / (y2prime - y1)) *
      (abs((x2prime - x) / (x2prime - x1)) * red(prevColor[x1][y1]) + abs((x - x1) / (x2prime - x1)) * red(prevColor[x2][y1])) +
      abs((y - y1) / (y2prime - y1)) * (abs((x2prime - x) / (x2prime - x1)) * red(prevColor[x1][y2]) +
      abs((x - x1) / (x2prime - x1)) * red(prevColor[x2][y2]));

      float gInterpolated = abs((y2prime - y) / (y2prime - y1)) *
      (abs((x2prime - x) / (x2prime - x1)) * green(prevColor[x1][y1]) + abs((x - x1) / (x2prime - x1)) * green(prevColor[x2][y1])) +
      abs((y - y1) / (y2prime - y1)) * (abs((x2prime - x) / (x2prime - x1)) * green(prevColor[x1][y2]) +
      abs((x - x1) / (x2prime - x1)) * green(prevColor[x2][y2]));

      float bInterpolated = abs((y2prime - y) / (y2prime - y1)) *
      (abs((x2prime - x) / (x2prime - x1)) * blue(prevColor[x1][y1]) + abs((x - x1) / (x2prime - x1)) * blue(prevColor[x2][y1])) +
      abs((y - y1) / (y2prime - y1)) * (abs((x2prime - x) / (x2prime - x1)) * blue(prevColor[x1][y2]) +
      abs((x - x1) / (x2prime - x1)) * blue(prevColor[x2][y2]));

      currentColor[i][j] = color(clamp(0, 255, rInterpolated),
        clamp(0, 255, gInterpolated),
        clamp(0, 255, bInterpolated));
    }
  }
}

float clamp(float minVal, float maxVal, float val) {
  return max(minVal, min(maxVal, val));
}

// Determine the previous i and j using previous velocity.
void advectVelocity() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      float prevU = uPrev[i][j];
      float prevV = vPrev[i][j];

      // Added some randomness here to make it look a bit more natural.
      float du = prevU * dt + random();
      float dv = prevV * dt + random();

      float x = ((i - du) + gridSize * gridSize) % gridSize;
      float y = ((j - dv) + gridSize * gridSize) % gridSize;

      int x1 = floor(x);
      int y1 = floor(y);

      int x2 = (x1 + 1) % gridSize;
      int y2 = (y1 + 1) % gridSize;

      int x2prime = x1 + 1;
      int y2prime = y1 + 1;

      float uInterpolated = abs((y2prime - y) / (y2prime - y1)) *
      (abs((x2prime - x) / (x2prime - x1)) * uPrev[x1][y1] + abs((x - x1) / (x2prime - x1)) * uPrev[x2][y1]) +
      abs((y - y1) / (y2prime - y1)) * (abs((x2prime - x) / (x2prime - x1)) * uPrev[x1][y2] +
      abs((x - x1) / (x2prime - x1)) * uPrev[x2][y2]);

      float vInterpolated = abs((y2prime - y) / (y2prime - y1)) *
      (abs((x2prime - x) / (x2prime - x1)) * vPrev[x1][y1] + abs((x - x1) / (x2prime - x1)) * vPrev[x2][y1]) +
      abs((y - y1) / (y2prime - y1)) * (abs((x2prime - x) / (x2prime - x1)) * vPrev[x1][y2] +
      abs((x - x1) / (x2prime - x1)) * vPrev[x2][y2]);

      u[i][j] = uInterpolated;
      v[i][j] = vInterpolated;
    }
  }
}

void initColor() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      int x = floor(i / 16.0);
      int y = floor(j / 16.0);
      int mod = (x + y) % 2;
      if (mod == 0) {
        prevColor[i][j] = color(255, 255, 255);
        currentColor[i][j] = color(255, 255, 255);
      } else {
        prevColor[i][j] = color(0, 0, 0);
        currentColor[i][j] = color(0, 0, 0);
      }
    }
  }
}

// This is some arbitrary initial velocity.
void initVelocity() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      //uPrev[i][j] = 91 * sin(2 * PI * (float)(j) / (gridSize / 3.141));
      //vPrev[i][j] = 91 * sin(2 * PI * (float)(i) / (gridSize / 3.141));
      uPrev[i][j] = 0.0;
      vPrev[i][j] = 0.0;
    }
  }
}

void initPressure() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      prevPressure[i][j] = 0.0;
    }
  }
}

void addForces() {
  boolean mouseMoved = (mouseX != pmouseX) || (mouseY != pmouseY);
  if (mousePressed && mouseMoved) {
    int x = (gridSize * mouseY) / height;
    int y = (gridSize * mouseX) / width;

    addForce(uPrev, x, y, x - previousMouseX, FORCE_RADIUS);
    addForce(vPrev, x, y, y - previousMouseY, FORCE_RADIUS);
  }

  previousMouseX = mouseX;
  previousMouseY = mouseY;
}

void addForce(float[][] field, int x, int y, float diff, float radius) {
  int i,j, dx, dy;
  float f;

    for ( i = int(clamp(x-radius, 0, gridSize - 1)); i <= int(clamp(x+radius, 0 , gridSize - 1)); i++ ) {
      dx = x - i;
      for ( j = int(clamp(y-radius, 0, gridSize - 1)); j <= int(clamp(y+radius, 0, gridSize - 1)); j++ ) {
        dy = y - j;
        f = 1 - ( sqrt(dx*dx + dy*dy) / radius );
        field[i][j] += clamp(f,0,1) * diff;
      }
    }
}

void drawColor() {
  loadPixels();
  for (int k = 0; k < pixels.length; k++) {
    int i = k / width;
    int j = k % width;

    // Sample from our color array.
    int iSample = (int) (((float) i / height) * gridSize);
    int jSample = (int) (((float) j / width) * gridSize);
    pixels[k] = currentColor[iSample][jSample];
  }
  updatePixels();
}
