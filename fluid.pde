float dt = 0.01;
int n = 256;
int gridSize = n + 2;
int pixelSize = 2;

int JACOBI_ITERATION = 10;
float FORCE_RADIUS = n / 25.0;

char lastKey;

float[][] u = new float[gridSize][gridSize];
float[][] v = new float[gridSize][gridSize];
float[][] uPrev = new float[gridSize][gridSize];
float[][] vPrev = new float[gridSize][gridSize];
float[][] d = new float[gridSize][gridSize];

float[][] currentPressure = new float[gridSize][gridSize];
float[][] prevPressure = new float[gridSize][gridSize];

color[][] prevColor = new color[gridSize][gridSize];
color[][] currentColor = new color[gridSize][gridSize];

PImage finklesteinImage;

void setup() {
  background(0);
  noStroke();
  size(gridSize * pixelSize, gridSize * pixelSize);

  initColor();
  initVelocity();
  initPressure();

  finklesteinImage = loadImage("finklestein.jpeg");
}

void draw() {
  if (keyPressed) {
    lastKey = key;
  }

  // Following the instructions given in the jos stam paper.
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

void keyReleased() {
  if (lastKey == 'f') {
    addFinklestein();
  } else if (lastKey == 'k') {
    makeCheckerBoard();
  } else if (lastKey == 'c') {
    clear();
  } else if (lastKey == 'r') {
    randomSplat();
  } else if (lastKey == 'v') {
    makeVortex();
  } else if (lastKey == 't') {
    stopMotion();
  }
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
      float du = prevU * dt + random(-1, 1);
      float dv = prevV * dt + random(-1, 1);

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
      prevColor[i][j] = color(0, 0, 0);
      currentColor[i][j] = color(0, 0, 0);
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
    float x = ((float)gridSize * mouseY) / height;
    float y = ((float)gridSize * mouseX) / width;
    float scale = 10.0

    addForce(uPrev, x, y, (mouseY - pmouseY) * scale, FORCE_RADIUS);
    addForce(vPrev, x, y, (mouseX - pmouseX) * scale, FORCE_RADIUS);
  }
}

void addForce(float[][] field, int x, int y, float diff, float radius) {
  int i,j, dx, dy;
  float f;

    for ( i = int(clamp(0, gridSize -1, x-radius)); i <= int(clamp(0 , gridSize - 1, x+radius)); i++ ) {
      dx = x - i;
      for ( j = int(clamp(0, gridSize - 1, y-radius)); j <= int(clamp(0, gridSize - 1, y+radius)); j++ ) {
        dy = y - j;
        f = 1 - ( sqrt(dx*dx + dy*dy) / radius );
        field[i][j] += clamp(0,1,f) * diff;
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

void addFinklestein() {
  int imageWidth = finklesteinImage.width;
  int imageHeight = finklesteinImage.height;

  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      int projectedI = round(((float) i / gridSize) * imageHeight);
      int projectedJ = round(((float) j / gridSize) * imageWidth);

      color c = finklesteinImage.get(projectedJ, projectedI);

      prevColor[i][j] = c;
      currentColor[i][j] = c;
    }
  }
}

void randomSplat() {
  for (int k = 0; k < 7; k++) {
    color randomColor = color((int) 255 * random(), (int) 255 * random(), (int) 255 * random());
    float x = (float) gridSize * random();
    float y = (float) gridSize * random();
    for (int i = (int) clamp(0, gridSize - 1, x - FORCE_RADIUS); i <= (int) clamp(0, gridSize - 1, x + FORCE_RADIUS); i++) {
      for (int j = (int) clamp(0, gridSize - 1, y - FORCE_RADIUS); j <= (int) clamp(0, gridSize - 1, y + FORCE_RADIUS); j++) {
        dx = x - i;
        dy = y - j;
        if (sqrt(dx*dx + dy*dy) < FORCE_RADIUS) {
          currentColor[i][j] = randomColor;
        }
      }
    }
  }
}

void clear() {
  initColor();
  initVelocity();
  initPressure();
}

void makeCheckerBoard() {
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

void makeVortex() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      uPrev[i][j] = (-(j - gridSize / 2.0) - (i - gridSize / 2.0)) * 3.0;
      vPrev[i][j] = ((i - gridSize/2.0) - (j - gridSize/2.0)) * 3.0;
    }
  }
}

void stopMotion() {
  initVelocity();
  initPressure();
}

void makeVortex() {
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      float x = i - gridSize/2.0 + random(-0.001, 0.001);
      float y = j - gridSize/2.0 + random(-0.001, 0.001);
      float theta = atan(y/x);
      if ((x < 0.0 && y > 0.0) || (x < 0.0 && y < 0.0)) {
        theta += PI;
      }
      float velDirection = theta + PI/2.0;
      float r = sqrt(x*x + y*y);
      uPrev[i][j] = 5.0 * r * cos(velDirection);
      vPrev[i][j] = 5.0 * r * sin(velDirection);
    }
  }
}
