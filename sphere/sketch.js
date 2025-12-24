"use strict";

const strokeWeight = 1;

// Initial Radius
const INITIAL_RADIUS = 100;
// Incr determines the number of lines of lat and lon
const INCR = 1.0;
// R is the radius, which is controlled by the slider
let r = INITIAL_RADIUS;
const MIN_RADIUS = 10,
  MAX_RADIUS = 6800;

// Rendering constants
const ELLIPSE_SEGMENTS = 50;
const CAMERA_FOV_DEGREES = 60;

// Laser constants
const INCIDENT_COLOR = [255, 0, 0]; // Red color for incident ray
const REFLECTED_COLOR = [255, 0, 0]; // Red color for reflected ray
const NORMAL_COLOR = [0, 255, 0]; // Green color for surface normal
const LASER_STROKE_WEIGHT = 3;
const NORMAL_LENGTH = 100; // Length of normal vector to draw
// Increase reflection length to reach grid
const REFLECTION_LENGTH = 2000; // Length of reflected beam to draw
const LASER_DISTANCE = 400;
const CAMERA_OFFSET_ANGLE = 15; // Degrees

// Grid Constants
const GRID_DISTANCE = 1200;
const GRID_SIZE = 400;
const GRID_RES = 20; // Number of lines

let slider;
let cameraYSlider;
let cameraXSlider;

const HEIGHT = 800;
const WIDTH = HEIGHT;

// Vector math helper functions
const vec3 = {
  dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
  add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
  scale: (v, s) => [v[0] * s, v[1] * s, v[2] * s],
  length: (v) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]),
  normalize: (v) => {
    const len = vec3.length(v);
    return len > 0 ? vec3.scale(v, 1 / len) : [0, 0, 0];
  },
  cross: (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ]
};

const raySphereIntersection = (rayOrigin, rayDir, sphereRadius, sphereCenter = [0, 0, 0]) => {
  const originToCenter = vec3.sub(rayOrigin, sphereCenter);
  const a = vec3.dot(rayDir, rayDir);
  const b = 2 * vec3.dot(originToCenter, rayDir);
  const c = vec3.dot(originToCenter, originToCenter) - sphereRadius * sphereRadius;
  const discriminant = b * b - 4 * a * c;

  if (discriminant < 0) return null;

  const sqrtDiscriminant = Math.sqrt(discriminant);
  const t1 = (-b - sqrtDiscriminant) / (2 * a);
  const t2 = (-b + sqrtDiscriminant) / (2 * a);

  let t = null;
  if (t1 > 0) t = t1;
  else if (t2 > 0) t = t2;
  else return null;

  const point = vec3.add(rayOrigin, vec3.scale(rayDir, t));
  return { hit: true, point, t };
};

const rayPlaneIntersection = (rayOrigin, rayDir, planePoint, planeNormal) => {
  const denom = vec3.dot(rayDir, planeNormal);
  if (Math.abs(denom) < 1e-6) return null; // Parallel
  const t = vec3.dot(vec3.sub(planePoint, rayOrigin), planeNormal) / denom;
  if (t < 0) return null; // Behind ray
  return { hit: true, point: vec3.add(rayOrigin, vec3.scale(rayDir, t)), t };
};

const reflect = (incident, normal) => {
  const dot = vec3.dot(incident, normal);
  return vec3.sub(incident, vec3.scale(normal, 2 * dot));
};

const getLaserRay = (sphereRadius) => {
  const startX = LASER_DISTANCE;
  const startY = 0;
  const startZ = 0;

  const reflectionAngleRad = (CAMERA_OFFSET_ANGLE / 2) * (Math.PI / 180);

  const targetX = sphereRadius * Math.cos(reflectionAngleRad);
  const targetY = 0;
  const targetZ = sphereRadius * Math.sin(reflectionAngleRad);

  const start = [startX, startY, startZ];
  const target = [targetX, targetY, targetZ];

  const direction = vec3.normalize(vec3.sub(target, start));

  return { origin: start, direction };
};

// Helper to calculate coordinate system for a plane defined by normal
const getPlaneBasis = (normal) => {
  let up = [0, 1, 0];
  if (Math.abs(vec3.dot(normal, up)) > 0.9) {
    up = [0, 0, 1];
  }
  const right = vec3.normalize(vec3.cross(normal, up));
  const newUp = vec3.normalize(vec3.cross(right, normal));
  return { right, up: newUp };
};

// Calculate the Static Grid Parameters (Fixed World State)
// This uses INITIAL_RADIUS and offset 0 to enforce a fixed grid position
const getStaticGridParams = () => {
  const sphereRadius = INITIAL_RADIUS / 2;
  const sphereCenter = [0, 0, 0];

  const laserRay = getLaserRay(sphereRadius);
  const intersection = raySphereIntersection(laserRay.origin, laserRay.direction, sphereRadius, sphereCenter);

  if (intersection && intersection.hit) {
    const normal = vec3.normalize(vec3.sub(intersection.point, sphereCenter));
    const reflectedDir = reflect(laserRay.direction, normal);
    const gridCenter = vec3.add(intersection.point, vec3.scale(reflectedDir, GRID_DISTANCE));
    const basis = getPlaneBasis(reflectedDir);
    return { center: gridCenter, normal: reflectedDir, basis, valid: true };
  }
  return { valid: false };
};

// Shared scene drawing function
const drawScene = (sketch, currentRadius, sphereYOffset) => {
  const doubleRuled = 5;
  const sphereCenter = [0, sphereYOffset, 0];
  const sphereRadius = currentRadius / 2;

  // 1. Get Static Grid (Always draw the grid in the same place)
  const staticGrid = getStaticGridParams();

  // Draw Static Grid
  if (staticGrid.valid) {
    const { center: gridCenter, normal: gridNormal, basis } = staticGrid;

    // Lazy initialize texture for this sketch instance
    if (!sketch.gridTexture) {
      const texSize = 1024;
      const pg = sketch.createGraphics(texSize, texSize);
      const scaleFactor = texSize / GRID_SIZE;

      pg.translate(texSize / 2, texSize / 2);
      pg.scale(scaleFactor);

      // Clear and Background
      pg.clear();
      pg.noStroke();
      pg.fill(200, 200, 255, 100);
      pg.rectMode(sketch.CENTER);
      pg.rect(0, 0, GRID_SIZE, GRID_SIZE);

      const step = GRID_SIZE / GRID_RES; // 20 units
      const halfSize = GRID_SIZE / 2;    // 200 units
      const halfRes = GRID_RES / 2;      // 10 lines

      pg.textSize(24);
      pg.textAlign(sketch.LEFT, sketch.BOTTOM);

      // Loop from -10 to 10
      for (let i = -halfRes; i <= halfRes; i++) {
        const pos = i * step;

        // Draw Lines
        if (i === 0) {
          // Axes
          pg.stroke(0);
          pg.strokeWeight(3);
        } else {
          // Grid lines
          pg.stroke(0, 0, 0, 100);
          pg.strokeWeight(1.5);
        }

        // Vertical Line (constant X = pos)
        pg.line(pos, -halfSize, pos, halfSize);
        // Horizontal Line (constant Y = pos)
        pg.line(-halfSize, pos, halfSize, pos);

        // Labels
        pg.noStroke();
        pg.fill(0);

        // Origin
        if (i === 0) {
          pg.text("(0,0)", 5, -5);
        }
      }
      sketch.gridTexture = pg;
    }

    sketch.push();
    sketch.translate(gridCenter[0], gridCenter[1], gridCenter[2]);

    // Align local plane with world grid basis
    const R = basis.right;
    const U = basis.up;
    const N = gridNormal;

    sketch.applyMatrix(
      R[0], U[0], N[0], 0,
      R[1], U[1], N[1], 0,
      R[2], U[2], N[2], 0,
      0, 0, 0, 1
    );

    // Offset slightly towards source (-Z in local) to avoid z-fighting/overlap
    sketch.translate(0, 0, -1);

    sketch.texture(sketch.gridTexture);
    sketch.noStroke();
    sketch.plane(GRID_SIZE, GRID_SIZE);
    sketch.pop();
  }

  // 2. Calculate Dynamic Laser Physics
  const laserRay = getLaserRay(sphereRadius);
  const intersection = raySphereIntersection(laserRay.origin, laserRay.direction, sphereRadius, sphereCenter);

  if (intersection && intersection.hit) {
    const hitPoint = intersection.point;
    const normal = vec3.normalize(vec3.sub(hitPoint, sphereCenter));
    const reflectedDir = reflect(laserRay.direction, normal);

    // Incident
    sketch.push();
    sketch.stroke(...INCIDENT_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      laserRay.origin[0], laserRay.origin[1], laserRay.origin[2],
      hitPoint[0], hitPoint[1], hitPoint[2]
    );
    sketch.pop();

    // Normal
    const normalEnd = vec3.add(hitPoint, vec3.scale(normal, NORMAL_LENGTH));
    sketch.push();
    sketch.stroke(...NORMAL_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      hitPoint[0], hitPoint[1], hitPoint[2],
      normalEnd[0], normalEnd[1], normalEnd[2]
    );
    sketch.pop();

    // Reflected (Dynamic)
    // Find where it hits static grid
    let reflectionEnd = vec3.add(hitPoint, vec3.scale(reflectedDir, REFLECTION_LENGTH));
    let gridHitPoint = null;

    if (staticGrid.valid) {
      const gridHit = rayPlaneIntersection(hitPoint, reflectedDir, staticGrid.center, staticGrid.normal);
      // We only care if it hits in general direction (t > 0), simple plane intersection
      if (gridHit && gridHit.hit) {
        // Let's extend the beam at least to the grid, but if it goes past, fine.
        // If gridHit.t < REFLECTION_LENGTH, maybe stop there? Or go through?
        // "see the laser as it sweeps across"
        // Drawing the spot is crucial.
        gridHitPoint = gridHit.point;

        // If the hit point is further than standard length, extend line
        if (gridHit.t > REFLECTION_LENGTH) {
          reflectionEnd = gridHitPoint;
        }
      }
    }

    sketch.push();
    sketch.stroke(...REFLECTED_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      hitPoint[0], hitPoint[1], hitPoint[2],
      reflectionEnd[0], reflectionEnd[1], reflectionEnd[2]
    );
    sketch.pop();

    // Draw Dynamic Spot on Static Grid
    if (gridHitPoint) {
      sketch.push();
      sketch.translate(gridHitPoint[0], gridHitPoint[1], gridHitPoint[2]);
      // Move slightly towards source to avoid z-fighting
      const spotOffset = vec3.scale(reflectedDir, -2);
      sketch.translate(spotOffset[0], spotOffset[1], spotOffset[2]);
      sketch.fill(255, 0, 0);
      sketch.noStroke();
      sketch.sphere(8); // Visible dot
      sketch.pop();
    }
  }

  // Laser source
  sketch.push();
  sketch.translate(laserRay.origin[0], laserRay.origin[1], laserRay.origin[2]);
  sketch.fill(255, 255, 0);
  sketch.noStroke();
  sketch.sphere(8);
  sketch.push();
  const targetDir = laserRay.direction;
  sketch.rotateZ(Math.atan2(targetDir[1], targetDir[0]) * (180 / Math.PI));
  sketch.rotateY(-Math.asin(targetDir[2]) * (180 / Math.PI));
  sketch.fill(255, 200, 0);
  sketch.translate(10, 0, 0);
  sketch.rotateZ(90);
  sketch.cone(6, 15);
  sketch.pop();
  sketch.pop();

  // Sphere
  sketch.push();
  sketch.translate(0, sphereYOffset, 0);

  // Longitude
  for (let lon = 0; lon < 180; lon += INCR) {
    sketch.push();
    sketch.rotateY(lon);
    if (lon % doubleRuled == 0) {
      sketch.strokeWeight(2 * strokeWeight);
      sketch.stroke("black");
    } else {
      sketch.strokeWeight(strokeWeight);
      sketch.stroke("darkgrey");
    }
    sketch.ellipse(0, 0, currentRadius, currentRadius, ELLIPSE_SEGMENTS);
    sketch.pop();
  }

  // Latitude
  sketch.rotateX(90);
  for (let lat = -90; lat < 90; lat += INCR / 2) {
    sketch.push();
    if (Math.abs(lat) % (doubleRuled / 2) == 0) {
      sketch.strokeWeight(2 * strokeWeight);
      sketch.stroke("black");
    } else {
      sketch.strokeWeight(strokeWeight);
      sketch.stroke("darkgrey");
    }
    let dz = (currentRadius * Math.cos((Math.PI / 90) * lat)) / 2;
    sketch.translate(0, 0, dz);
    let latRadius = Math.abs(currentRadius * Math.sin((Math.PI / 90) * lat));
    sketch.ellipse(0, 0, latRadius, latRadius, ELLIPSE_SEGMENTS);
    sketch.pop();
  }
  sketch.pop(); // End sphere

  return { laserRay, gridInfo: staticGrid };
};

// --- VIEW 1: Main Reflection View ---
const topCanvas = (sketch) => {
  let cam;
  sketch.setup = () => {
    let canvas = sketch.createCanvas(HEIGHT, WIDTH, sketch.WEBGL);
    canvas.parent('canvas-container');
    sketch.angleMode(sketch.DEGREES);
    sketch.strokeWeight(strokeWeight);
    cam = sketch.createCamera();
  };
  sketch.draw = () => {
    if (slider) r = Math.exp(slider.value());
    sketch.background("white");

    // Stop oscillation by setting offset to 0
    const sphereYOffset = Math.sin(sketch.frameCount * 0.05) * 1;

    drawScene(sketch, r, sphereYOffset);

    // Camera tracks the DYNAMIC laser geometry if we want, or static?
    // Original request was "reflection off sphere shines directly into camera".
    // If r changes, reflection changes.
    // If we want it to shine into camera, camera must move OR laser/sphere must adjust.
    // Current code moves camera to match reflection.
    let camRad = CAMERA_OFFSET_ANGLE * (Math.PI / 180);
    let camX = LASER_DISTANCE * Math.cos(camRad);
    let camY = 0;
    let camZ = LASER_DISTANCE * Math.sin(camRad);

    cam.setPosition(camX, camY, camZ);
    cam.lookAt(0, 0, 0);
    cam.perspective((CAMERA_FOV_DEGREES * INITIAL_RADIUS) / r);
  };
};

// --- VIEW 2: Laser View ---
const laserCanvas = (sketch) => {
  let cam;
  sketch.setup = () => {
    let canvas = sketch.createCanvas(HEIGHT, WIDTH, sketch.WEBGL);
    canvas.parent('canvas-container');
    sketch.angleMode(sketch.DEGREES);
    sketch.strokeWeight(strokeWeight);
    cam = sketch.createCamera();
  };
  sketch.draw = () => {
    sketch.background("white");
    const sphereYOffset = Math.sin(sketch.frameCount * 0.05) * 1;
    const sceneInfo = drawScene(sketch, r, sphereYOffset);

    const laserRay = sceneInfo.laserRay;
    const backUpDist = 50;
    const camPos = vec3.sub(laserRay.origin, vec3.scale(laserRay.direction, backUpDist));
    const lookAtPos = vec3.add(laserRay.origin, vec3.scale(laserRay.direction, 100));

    cam.setPosition(camPos[0], camPos[1], camPos[2]);
    cam.lookAt(lookAtPos[0], lookAtPos[1], lookAtPos[2]);
    cam.perspective(60 * Math.PI / 180);
  };
};

// --- VIEW 3: Grid View ---
const gridCanvas = (sketch) => {
  let cam;
  sketch.setup = () => {
    let canvas = sketch.createCanvas(HEIGHT, WIDTH, sketch.WEBGL);
    canvas.parent('canvas-container');
    sketch.angleMode(sketch.DEGREES);
    sketch.strokeWeight(strokeWeight);
    cam = sketch.createCamera();
  };
  sketch.draw = () => {
    sketch.background("white");
    const sphereYOffset = Math.sin(sketch.frameCount * 0.05) * 1;
    // drawScene returns the static grid info as gridInfo because we changed implementation
    const sceneInfo = drawScene(sketch, r, sphereYOffset);

    if (sceneInfo.gridInfo && sceneInfo.gridInfo.valid) {
      const gridCenter = sceneInfo.gridInfo.center;
      const gridNormal = sceneInfo.gridInfo.normal;

      const distFromGrid = 3000;

      const basis = sceneInfo.gridInfo.basis;
      const offsetScale = 200;

      const camOffset = vec3.add(
        vec3.scale(gridNormal, -distFromGrid),
        vec3.add(vec3.scale(basis.right, offsetScale), vec3.scale(basis.up, offsetScale))
      );

      const camPos = vec3.add(gridCenter, camOffset);

      cam.setPosition(camPos[0], camPos[1], camPos[2]);
      cam.lookAt(gridCenter[0], gridCenter[1], gridCenter[2]);
    }
    cam.perspective(600 * Math.PI / 180);
  };
};

// --- VIEW 4: Overhead View ---
const overheadCanvas = (sketch) => {
  let cam;
  sketch.setup = () => {
    let canvas = sketch.createCanvas(HEIGHT, WIDTH, sketch.WEBGL);
    canvas.parent('canvas-container');
    sketch.angleMode(sketch.DEGREES);
    sketch.strokeWeight(strokeWeight);
    cam = sketch.createCamera();
    // Position: Above the scene
    // Target: Midpoint of scene
    // Up: Z-axis (so North is "up" on screen, or East is "right")
    cam.camera(600, -4000, 0, 600, 0, 0, 0, 0, 1);
  };
  sketch.draw = () => {
    sketch.background("white");
    const sphereYOffset = Math.sin(sketch.frameCount * 0.05) * 1;
    drawScene(sketch, r, sphereYOffset);
    // Camera is static, set in setup
  };
};

let topP5 = new p5(topCanvas);
let laserP5 = new p5(laserCanvas);
let gridP5 = new p5(gridCanvas);
let overheadP5 = new p5(overheadCanvas);

// Controls
let radiusTxt = (radius) => {
  return `radius: ${radius.toFixed(1).padStart(4)}`;
};

let bottomCanvas = (sketch) => {
  const X = WIDTH * 0.05, Y = 50, W = 200, L = 20;

  sketch.setup = () => {
    let canvas = sketch.createCanvas(400, 200);
    canvas.parent('controls-container');

    let sliderStep = 0.01;
    slider = sketch.createSlider(
      Math.log(MIN_RADIUS),
      Math.log(MAX_RADIUS),
      Math.log(INITIAL_RADIUS),
      sliderStep
    );
    slider.parent('controls-container');
    slider.style('width', '200px');
    slider.fValue = () => { return Math.exp(slider.value()); };
  };

  sketch.draw = () => {
    sketch.background("white");
    sketch.textSize(15);
    sketch.textAlign(sketch.LEFT, sketch.CENTER);
    sketch.noStroke();
    sketch.fill(0);

    let radius = slider.fValue();
    sketch.text(radiusTxt(radius), 10, 40);
    sketch.text("Radius Control", 10, 10);
  };
};

let bottomP5 = new p5(bottomCanvas);
