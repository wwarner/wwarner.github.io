'use strict';

const strokeWeight = 1.5;

// Initial Radius
const INITIAL_RADIUS = 50;
// Incr determines the number of lines of lat and lon
const INCR = 10.0;
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

// Glow Constants
const GLOW_RES = 64;
const GLOW_DECAY = 0.94;
let glowMap = new Float32Array(GLOW_RES * GLOW_RES);

// Camera Vantage Points
const GRID_VANTAGE_POINT = {
  distFromGrid: 7500,
  offsetScale: 150,
  fov: 95
};

const REFLECTION_VANTAGE_POINT = {
  camDistance: 150,  // Distance back along incident ray
  upOffset: 30,      // Offset above intersection
  fov: 60
};

// Camera animation timing (in seconds)
const PAUSE_DURATION = 10;  // Pause at each vantage point
const TRANSITION_DURATION = 15;  // Total time to move between vantage points
const PHASE_DURATION = TRANSITION_DURATION / 4; // Each transition has 4 phases (3.75s each)
const CYCLE_DURATION = (PAUSE_DURATION + TRANSITION_DURATION) * 2; // 50 seconds total

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
    a[0] * b[1] - a[1] * b[0],
  ],
  lerp: (a, b, t) => [
    a[0] + (b[0] - a[0]) * t,
    a[1] + (b[1] - a[1]) * t,
    a[2] + (b[2] - a[2]) * t,
  ],
};

// Smooth easing function (ease-in-out cubic)
const easeInOutCubic = (t) => {
  return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
};

// Linear interpolation for scalars
const lerp = (a, b, t) => a + (b - a) * t;

// Spherical linear interpolation (slerp) for smooth rotation between two directions
const slerp = (a, b, t) => {
  const dotProduct = vec3.dot(vec3.normalize(a), vec3.normalize(b));
  const clampedDot = Math.max(-1, Math.min(1, dotProduct));
  const theta = Math.acos(clampedDot) * t;
  
  const relativeVec = vec3.normalize(vec3.sub(b, vec3.scale(a, clampedDot)));
  
  return vec3.add(
    vec3.scale(a, Math.cos(theta)),
    vec3.scale(relativeVec, Math.sin(theta))
  );
};

const getDeterministicRandom = (input) => {
  return 0.05;
  let x = Math.sin(input * 10000) * 10000;
  return (x - Math.floor(x)) * 0.005;
};

const raySphereIntersection = (
  rayOrigin,
  rayDir,
  sphereRadius,
  sphereCenter = [0, 0, 0]
) => {
  const originToCenter = vec3.sub(rayOrigin, sphereCenter);
  const a = vec3.dot(rayDir, rayDir);
  const b = 2 * vec3.dot(originToCenter, rayDir);
  const c =
    vec3.dot(originToCenter, originToCenter) - sphereRadius * sphereRadius;
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
  const intersection = raySphereIntersection(
    laserRay.origin,
    laserRay.direction,
    sphereRadius,
    sphereCenter
  );

  if (intersection && intersection.hit) {
    const normal = vec3.normalize(vec3.sub(intersection.point, sphereCenter));
    const reflectedDir = reflect(laserRay.direction, normal);
    const gridCenter = vec3.add(
      intersection.point,
      vec3.scale(reflectedDir, GRID_DISTANCE)
    );
    const basis = getPlaneBasis(reflectedDir);
    return { center: gridCenter, normal: reflectedDir, basis, valid: true };
  }
  return { valid: false };
};

// Shared scene drawing function
const drawScene = (sketch, currentRadius, sphereYOffset, isMaster = false) => {
  // Update Glow Physics (Decay) if master
  if (isMaster) {
    for (let i = 0; i < glowMap.length; i++) {
      if (glowMap[i] > 0.001) glowMap[i] *= GLOW_DECAY;
      else glowMap[i] = 0;
    }
  }

  const doubleRuled = 5;
  const sphereCenter = [0, sphereYOffset, 0];
  const sphereRadius = currentRadius / 2;

  // 1. Get Static Grid (Always draw the grid in the same place)
  const staticGrid = getStaticGridParams();

  // Draw Static Grid
  if (staticGrid.valid) {
    const { center: gridCenter, normal: gridNormal, basis } = staticGrid;

    sketch.push();
    sketch.translate(gridCenter[0], gridCenter[1], gridCenter[2]);

    // Align local plane with world grid basis
    const R = basis.right;
    const U = basis.up;
    const N = gridNormal;

    sketch.applyMatrix(
      R[0],
      U[0],
      N[0],
      0,
      R[1],
      U[1],
      N[1],
      0,
      R[2],
      U[2],
      N[2],
      0,
      0,
      0,
      0,
      1
    );

    // Draw grid lines using p5 primitives (like the sphere)
    sketch.push();
    sketch.translate(0, 0, -1);
    
    // Set line properties to match sphere
    sketch.strokeWeight(strokeWeight);
    sketch.stroke('#444444');
    
    const halfSize = GRID_SIZE / 2;
    const step = GRID_SIZE / GRID_RES;
    
    // Draw vertical lines
    for (let i = 0; i <= GRID_RES; i++) {
      const x = -halfSize + i * step;
      sketch.line(x, -halfSize, x, halfSize);
    }
    
    // Draw horizontal lines
    for (let i = 0; i <= GRID_RES; i++) {
      const y = -halfSize + i * step;
      sketch.line(-halfSize, y, halfSize, y);
    }
    
    sketch.pop();

    // Draw Glow Overlay
    if (!sketch.glowImg) {
      sketch.glowImg = sketch.createImage(GLOW_RES, GLOW_RES);
    }
    sketch.glowImg.loadPixels();
    for (let i = 0; i < GLOW_RES * GLOW_RES; i++) {
      const val = glowMap[i];
      if (val > 0.01) {
        const idx = i * 4;
        // Red Phosphor color to match laser
        sketch.glowImg.pixels[idx] = 255; // R
        sketch.glowImg.pixels[idx + 1] = 0; // G
        sketch.glowImg.pixels[idx + 2] = 0; // B
        sketch.glowImg.pixels[idx + 3] = Math.min(255, val * 255); // Alpha
      } else {
        const idx = i * 4;
        sketch.glowImg.pixels[idx + 3] = 0;
      }
    }
    sketch.glowImg.updatePixels();

    sketch.push();
    sketch.translate(0, 0, -2); // Slightly above grid
    sketch.texture(sketch.glowImg);
    sketch.noStroke();
    sketch.plane(GRID_SIZE, GRID_SIZE);
    sketch.pop();

    sketch.pop();
  }

  // 2. Calculate Dynamic Laser Physics
  const laserRay = getLaserRay(sphereRadius);
  const intersection = raySphereIntersection(
    laserRay.origin,
    laserRay.direction,
    sphereRadius,
    sphereCenter
  );

  if (intersection && intersection.hit) {
    const hitPoint = intersection.point;
    const normal = vec3.normalize(vec3.sub(hitPoint, sphereCenter));
    const reflectedDir = reflect(laserRay.direction, normal);

    // Incident
    sketch.push();
    sketch.stroke(...INCIDENT_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      laserRay.origin[0],
      laserRay.origin[1],
      laserRay.origin[2],
      hitPoint[0],
      hitPoint[1],
      hitPoint[2]
    );
    sketch.pop();

    // Normal
    const normalEnd = vec3.add(hitPoint, vec3.scale(normal, NORMAL_LENGTH));
    sketch.push();
    sketch.stroke(...NORMAL_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      hitPoint[0],
      hitPoint[1],
      hitPoint[2],
      normalEnd[0],
      normalEnd[1],
      normalEnd[2]
    );
    sketch.pop();

    // Reflected (Dynamic)
    // Find where it hits static grid
    let reflectionEnd = vec3.add(
      hitPoint,
      vec3.scale(reflectedDir, REFLECTION_LENGTH)
    );
    let gridHitPoint = null;

    if (staticGrid.valid) {
      const gridHit = rayPlaneIntersection(
        hitPoint,
        reflectedDir,
        staticGrid.center,
        staticGrid.normal
      );
      if (gridHit && gridHit.hit) {
        gridHitPoint = gridHit.point;

        // If the hit point is further than standard length, extend line
        if (gridHit.t > REFLECTION_LENGTH) {
          reflectionEnd = gridHitPoint;
        }

        // Update Glow Map if master
        if (isMaster) {
          // Project hit point to local grid coordinates
          const diff = vec3.sub(gridHitPoint, staticGrid.center);
          const uLocal = vec3.dot(diff, staticGrid.basis.right);
          const vLocal = vec3.dot(diff, staticGrid.basis.up);

          // Map to texture UV [0..1]
          // Grid ranges from -GRID_SIZE/2 to +GRID_SIZE/2
          const u = (uLocal + GRID_SIZE / 2) / GRID_SIZE;
          const v = (vLocal + GRID_SIZE / 2) / GRID_SIZE;

          if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
            // Map to array index
            const xInd = Math.floor(u * GLOW_RES);
            const yInd = Math.floor(v * GLOW_RES); // Invert y? p5 plane UV is standard.

            // Splat glow
            const rSplat = 2;
            for (let dy = -rSplat; dy <= rSplat; dy++) {
              for (let dx = -rSplat; dx <= rSplat; dx++) {
                const nx = xInd + dx;
                const ny = yInd + dy;
                if (nx >= 0 && nx < GLOW_RES && ny >= 0 && ny < GLOW_RES) {
                  const idx = nx + ny * GLOW_RES;
                  const distSq = dx * dx + dy * dy;
                  const amount = Math.exp(-distSq * 0.5); // Gaussian-ish
                  glowMap[idx] = Math.min(1.0, glowMap[idx] + amount * 0.2); // Add brightness
                }
              }
            }
          }
        }
      }
    }

    sketch.push();
    sketch.stroke(...REFLECTED_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      hitPoint[0],
      hitPoint[1],
      hitPoint[2],
      reflectionEnd[0],
      reflectionEnd[1],
      reflectionEnd[2]
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
  sketch.rotateY(sketch.frameCount * 0.2);

  // Longitude
  for (let lon = 0; lon < 180; lon += INCR) {
    sketch.push();
    sketch.rotateY(lon);
    if (lon % doubleRuled == 0) {
      sketch.strokeWeight(strokeWeight);
      sketch.stroke('#444444');
    }
    sketch.ellipse(0, 0, currentRadius, currentRadius, ELLIPSE_SEGMENTS);
    sketch.pop();
  }

  // Latitude
  sketch.rotateX(90);
  for (let lat = -90; lat < 90; lat += INCR / 2) {
    sketch.push();
    if (Math.abs(lat) % (doubleRuled / 2) == 0) {
      sketch.strokeWeight(strokeWeight);
      sketch.stroke('#444444');
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
    sketch.background('white');

    // Stop oscillation by setting offset to 0
    const freq = getDeterministicRandom(sketch.frameCount);
    const sphereYOffset = Math.sin(sketch.frameCount * freq) * 1;

    // isMaster = true for topCanvas to drive glow physics
    drawScene(sketch, r, sphereYOffset, true);

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
    sketch.background('white');
    const freq = getDeterministicRandom(sketch.frameCount);
    const sphereYOffset = Math.sin(sketch.frameCount * freq) * 1;
    const sceneInfo = drawScene(sketch, r, sphereYOffset, false);

    const laserRay = sceneInfo.laserRay;
    const backUpDist = 50;
    const camPos = vec3.sub(
      laserRay.origin,
      vec3.scale(laserRay.direction, backUpDist)
    );
    const lookAtPos = vec3.add(
      laserRay.origin,
      vec3.scale(laserRay.direction, 100)
    );

    cam.setPosition(camPos[0], camPos[1], camPos[2]);
    cam.lookAt(lookAtPos[0], lookAtPos[1], lookAtPos[2]);
    cam.perspective((60 * Math.PI) / 180);
  };
};

// 4-Phase camera transition calculator
// Phase 1: Point camera at sphere center AND move to R1 (GRID radius)
// Phase 2: Move along great circle at constant R1 (looking at sphere center)
// Phase 3: Change radius from R1 to destination radius
// Phase 4: Point camera at destination target
const calculateTransitionState = (
  startPos, startLookAt, startFov,
  endPos, endLookAt, endFov,
  sphereCenter, transitionTime
) => {
  const phase1End = PHASE_DURATION;
  const phase2End = PHASE_DURATION * 2;
  const phase3End = PHASE_DURATION * 3;
  const phase4End = PHASE_DURATION * 4;

  // R1 is the GRID_VANTAGE_POINT radius (the larger radius for traveling)
  const startDir = vec3.normalize(vec3.sub(startPos, sphereCenter));
  const startRadius = vec3.length(vec3.sub(startPos, sphereCenter));
  const endDir = vec3.normalize(vec3.sub(endPos, sphereCenter));
  const endRadius = vec3.length(vec3.sub(endPos, sphereCenter));
  
  // Always use the larger radius (GRID_VANTAGE_POINT) for great circle travel
  const R1 = Math.max(startRadius, endRadius);

  let camPos, lookAt, fov;

  if (transitionTime < phase1End) {
    // Phase 1: Point camera at sphere center AND move radius to R1
    const t = easeInOutCubic(transitionTime / PHASE_DURATION);
    
    // Interpolate lookAt to sphere center
    lookAt = vec3.lerp(startLookAt, sphereCenter, t);
    
    // Keep direction, interpolate radius to R1
    const currentRadius = lerp(startRadius, R1, t);
    camPos = vec3.add(sphereCenter, vec3.scale(startDir, currentRadius));
    
    fov = startFov;
  } else if (transitionTime < phase2End) {
    // Phase 2: Move along great circle at constant radius R1
    const t = easeInOutCubic((transitionTime - phase1End) / PHASE_DURATION);
    
    // Spherical interpolation for direction, constant radius R1
    const currentDir = slerp(startDir, endDir, t);
    camPos = vec3.add(sphereCenter, vec3.scale(currentDir, R1));
    lookAt = sphereCenter;
    fov = startFov;
  } else if (transitionTime < phase3End) {
    // Phase 3: Change radius from R1 to destination radius
    const t = easeInOutCubic((transitionTime - phase2End) / PHASE_DURATION);
    
    // Keep final direction, interpolate radius from R1 to endRadius
    const currentRadius = lerp(R1, endRadius, t);
    camPos = vec3.add(sphereCenter, vec3.scale(endDir, currentRadius));
    lookAt = sphereCenter;
    fov = lerp(startFov, endFov, t);
  } else {
    // Phase 4: Point camera at destination target
    const t = easeInOutCubic((transitionTime - phase3End) / PHASE_DURATION);
    camPos = endPos;
    lookAt = vec3.lerp(sphereCenter, endLookAt, t);
    fov = endFov;
  }

  return [camPos, lookAt, fov];
};

// --- VIEW 3: Grid View ---
const gridCanvas = (sketch) => {
  let cam;
  sketch.setup = () => {
    let container = document.getElementById('canvas-container');
    // Use container dimensions, falling back to window dimensions
    let w = container ? container.clientWidth : sketch.windowWidth;
    let h = container ? container.clientHeight : sketch.windowHeight;

    // Ensure minimum dimensions for mobile
    if (w === 0 || w < 100) w = Math.max(sketch.windowWidth, 300);
    if (h === 0 || h < 100) h = Math.max(sketch.windowHeight, 300);

    // Enable WebGL antialiasing and quality settings
    sketch.setAttributes('antialias', true);
    sketch.setAttributes('preserveDrawingBuffer', true);
    
    // Use higher pixel density for sharper rendering (especially on retina displays)
    sketch.pixelDensity(Math.min(window.devicePixelRatio || 1, 2));

    let canvas = sketch.createCanvas(w, h, sketch.WEBGL);
    canvas.parent('canvas-container');
    
    // Enable smoothing for better line quality
    sketch.smooth();
    
    sketch.angleMode(sketch.DEGREES);
    sketch.strokeWeight(strokeWeight);
    cam = sketch.createCamera();
  };

  sketch.windowResized = () => {
    let container = document.getElementById('canvas-container');
    if (container) {
      sketch.resizeCanvas(container.clientWidth, container.clientHeight);
    }
  };

  sketch.draw = () => {
    sketch.background('white');
    const freq = getDeterministicRandom(sketch.frameCount);
    const sphereYOffset = Math.sin(sketch.frameCount * freq) * 0.1;
    const sceneInfo = drawScene(sketch, r, sphereYOffset, true);

    // Calculate time in animation cycle
    const elapsedTime = sketch.millis() / 1000;
    const cycleTime = elapsedTime % CYCLE_DURATION;

    // Calculate laser intersection point
    const sphereRadius = r / 2;
    const sphereCenter = [0, sphereYOffset, 0];
    const laserRay = sceneInfo.laserRay;
    const intersection = raySphereIntersection(
      laserRay.origin,
      laserRay.direction,
      sphereRadius,
      sphereCenter
    );

    if (intersection && intersection.hit && sceneInfo.gridInfo && sceneInfo.gridInfo.valid) {
      const hitPoint = intersection.point;
      const gridCenter = sceneInfo.gridInfo.center;
      const gridNormal = sceneInfo.gridInfo.normal;
      const basis = sceneInfo.gridInfo.basis;

      // Calculate both vantage points
      // Grid vantage point
      const gridCamOffset = vec3.add(
        vec3.scale(gridNormal, -GRID_VANTAGE_POINT.distFromGrid),
        vec3.add(
          vec3.scale(basis.right, GRID_VANTAGE_POINT.offsetScale),
          vec3.scale(basis.up, GRID_VANTAGE_POINT.offsetScale)
        )
      );
      const gridCamPos = vec3.add(gridCenter, gridCamOffset);
      const gridLookAt = gridCenter;

      // Reflection vantage point
      const reflectionCamPos = vec3.add(
        hitPoint,
        vec3.add(
          vec3.scale(vec3.scale(laserRay.direction, -1), REFLECTION_VANTAGE_POINT.camDistance),
          [0, REFLECTION_VANTAGE_POINT.upOffset, 0]
        )
      );
      const reflectionLookAt = hitPoint;

      // Determine which transition we're in and calculate camera state
      let camPos, lookAt, fov;
      let isTransitioning = false;
      let isForward = false; // forward = grid to reflection, backward = reflection to grid

      if (cycleTime < PAUSE_DURATION) {
        // At grid vantage point
        camPos = gridCamPos;
        lookAt = gridLookAt;
        fov = GRID_VANTAGE_POINT.fov;
      } else if (cycleTime < PAUSE_DURATION + TRANSITION_DURATION) {
        // Transitioning from grid to reflection
        isTransitioning = true;
        isForward = true;
        const transitionTime = cycleTime - PAUSE_DURATION;
        [camPos, lookAt, fov] = calculateTransitionState(
          gridCamPos, gridLookAt, GRID_VANTAGE_POINT.fov,
          reflectionCamPos, reflectionLookAt, REFLECTION_VANTAGE_POINT.fov,
          sphereCenter, transitionTime
        );
      } else if (cycleTime < PAUSE_DURATION * 2 + TRANSITION_DURATION) {
        // At reflection vantage point
        camPos = reflectionCamPos;
        lookAt = reflectionLookAt;
        fov = REFLECTION_VANTAGE_POINT.fov;
      } else {
        // Transitioning from reflection back to grid
        isTransitioning = true;
        isForward = false;
        const transitionTime = cycleTime - PAUSE_DURATION * 2 - TRANSITION_DURATION;
        [camPos, lookAt, fov] = calculateTransitionState(
          reflectionCamPos, reflectionLookAt, REFLECTION_VANTAGE_POINT.fov,
          gridCamPos, gridLookAt, GRID_VANTAGE_POINT.fov,
          sphereCenter, transitionTime
        );
      }

      cam.setPosition(camPos[0], camPos[1], camPos[2]);
      cam.lookAt(lookAt[0], lookAt[1], lookAt[2]);

      let aspect = sketch.width / sketch.height;
      cam.perspective((fov * Math.PI) / 180, aspect, 1, 10000);
    }
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
    cam.camera(600, -4000, 0, 600, 0, 0, 0, 0, 1);
  };
  sketch.draw = () => {
    sketch.background('white');
    const freq = getDeterministicRandom(sketch.frameCount);
    const sphereYOffset = Math.sin(sketch.frameCount * freq) * 1;
    drawScene(sketch, r, sphereYOffset, false);
    // Camera is static, set in setup
  };
};

// let topP5 = new p5(topCanvas);
// let laserP5 = new p5(laserCanvas);
let gridP5 = new p5(gridCanvas);
// let overheadP5 = new p5(overheadCanvas);

// Controls
let radiusTxt = (radius) => {
  return `radius: ${radius.toFixed(1).padStart(4)}`;
};

let bottomCanvas = (sketch) => {
  const X = WIDTH * 0.05,
    Y = 50,
    W = 200,
    L = 20;

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
    slider.fValue = () => {
      return Math.exp(slider.value());
    };
  };

  sketch.draw = () => {
    sketch.background('white');
    sketch.textSize(15);
    sketch.textAlign(sketch.LEFT, sketch.CENTER);
    sketch.noStroke();
    sketch.fill(0);

    let radius = slider.fValue();
    sketch.text(radiusTxt(radius), 10, 40);
    sketch.text('Radius Control', 10, 10);
  };
};

// let bottomP5 = new p5(bottomCanvas);
