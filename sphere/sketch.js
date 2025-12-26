'use strict';

const GLOBAL_STROKE_WEIGHT = 1.5;

// Initial Radius
const INITIAL_RADIUS = 50;
// Incr determines the number of lines of lat and lon
const INCR = 10.0;

// Rendering constants
const ELLIPSE_SEGMENTS = 50;

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
  fov: 95,
};

const REFLECTION_VANTAGE_POINT = {
  camDistance: 150, // Distance back along incident ray
  upOffset: 30, // Offset above intersection
  fov: 60,
};

// Camera animation timing (in seconds)
const PAUSE_DURATION = 10; // Pause at each vantage point
const TRANSITION_DURATION = 15; // Total time to move between vantage points
const PHASE_DURATION = TRANSITION_DURATION / 4; // Each transition has 4 phases (3.75s each)
const CYCLE_DURATION = (PAUSE_DURATION + TRANSITION_DURATION) * 2; // 50 seconds total

// Smooth easing function (ease-in-out cubic)
const easeInOutCubic = (t) => {
  return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
};

// Linear interpolation for scalars
const lerp = (a, b, t) => a + (b - a) * t;

// Spherical linear interpolation (slerp) for smooth rotation between two directions
const slerp = (a, b, t) => {
  const v1 = a.copy().normalize();
  const v2 = b.copy().normalize();
  const dotProduct = v1.dot(v2);
  const clampedDot = Math.max(-1, Math.min(1, dotProduct));
  const theta = Math.acos(clampedDot);

  if (theta < 0.0001) return a.copy();

  const factor1 = Math.sin((1 - t) * theta) / Math.sin(theta);
  const factor2 = Math.sin(t * theta) / Math.sin(theta);

  return p5.Vector.add(
    p5.Vector.mult(v1, factor1),
    p5.Vector.mult(v2, factor2)
  ).mult(p5.Vector.lerp(a, b, t).mag()); // Scale by lerped magnitude
};

const raySphereIntersection = (
  rayOrigin,
  rayDir,
  sphereRadius,
  sphereCenter = new p5.Vector(0, 0, 0)
) => {
  const originToCenter = p5.Vector.sub(rayOrigin, sphereCenter);
  const a = rayDir.dot(rayDir);
  const b = 2 * originToCenter.dot(rayDir);
  const c = originToCenter.dot(originToCenter) - sphereRadius * sphereRadius;
  const discriminant = b * b - 4 * a * c;

  if (discriminant < 0) return null;

  const sqrtDiscriminant = Math.sqrt(discriminant);
  const t1 = (-b - sqrtDiscriminant) / (2 * a);
  const t2 = (-b + sqrtDiscriminant) / (2 * a);

  let t = null;
  if (t1 > 0) t = t1;
  else if (t2 > 0) t = t2;
  else return null;

  const point = p5.Vector.add(rayOrigin, p5.Vector.mult(rayDir, t));
  return { hit: true, point, t };
};

const rayPlaneIntersection = (rayOrigin, rayDir, planePoint, planeNormal) => {
  const denom = rayDir.dot(planeNormal);
  if (Math.abs(denom) < 1e-6) return null; // Parallel
  const t = p5.Vector.sub(planePoint, rayOrigin).dot(planeNormal) / denom;
  if (t < 0) return null; // Behind ray
  return {
    hit: true,
    point: p5.Vector.add(rayOrigin, p5.Vector.mult(rayDir, t)),
    t,
  };
};

const reflect = (incident, normal) => {
  const dot = incident.dot(normal);
  return p5.Vector.sub(incident, p5.Vector.mult(normal, 2 * dot));
};

const getLaserRay = (sphereRadius) => {
  const startX = LASER_DISTANCE;
  const startY = 0;
  const startZ = 0;

  const reflectionAngleRad = (CAMERA_OFFSET_ANGLE / 2) * (Math.PI / 180);

  const targetX = sphereRadius * Math.cos(reflectionAngleRad);
  const targetY = 0;
  const targetZ = sphereRadius * Math.sin(reflectionAngleRad);

  const start = new p5.Vector(startX, startY, startZ);
  const target = new p5.Vector(targetX, targetY, targetZ);

  const direction = p5.Vector.sub(target, start).normalize();

  return { origin: start, direction };
};

// Helper to calculate coordinate system for a plane defined by normal
const getPlaneBasis = (normal) => {
  let up = new p5.Vector(0, 1, 0);
  if (Math.abs(normal.dot(up)) > 0.9) {
    up = new p5.Vector(0, 0, 1);
  }
  const right = p5.Vector.cross(normal, up).normalize();
  const newUp = p5.Vector.cross(right, normal).normalize();
  return { right, up: newUp };
};

// Calculate the Static Grid Parameters (Fixed World State)
// This uses INITIAL_RADIUS and offset 0 to enforce a fixed grid position
const computeStaticGridParams = () => {
  const sphereRadius = INITIAL_RADIUS / 2;
  const sphereCenter = new p5.Vector(0, 0, 0);

  const laserRay = getLaserRay(sphereRadius);
  const intersection = raySphereIntersection(
    laserRay.origin,
    laserRay.direction,
    sphereRadius,
    sphereCenter
  );

  if (intersection && intersection.hit) {
    const normal = p5.Vector.sub(intersection.point, sphereCenter).normalize();
    const reflectedDir = reflect(laserRay.direction, normal);
    const gridCenter = p5.Vector.add(
      intersection.point,
      p5.Vector.mult(reflectedDir, GRID_DISTANCE)
    );
    const basis = getPlaneBasis(reflectedDir);
    return { center: gridCenter, normal: reflectedDir, basis, valid: true };
  }
  return { valid: false };
};

// Cache the static grid params since they never change (computed lazily on first use)
let staticGridParams = null;

// Shared scene drawing function
const drawScene = (sketch, currentRadius, sphereYOffset, isMaster = false) => {
  // Update Glow Physics (Decay) if master
  if (isMaster) {
    for (let i = 0; i < glowMap.length; i++) {
      if (glowMap[i] > 0.001) glowMap[i] *= GLOW_DECAY;
      else glowMap[i] = 0;
    }
  }

  const sphereCenter = new p5.Vector(0, sphereYOffset, 0);
  const sphereRadius = currentRadius / 2;

  // 1. Get Static Grid (Always draw the grid in the same place)
  // Lazy initialization: compute on first use to avoid p5.Vector issues at module load
  if (!staticGridParams) {
    staticGridParams = computeStaticGridParams();
  }
  const staticGrid = staticGridParams;

  // Draw Static Grid
  if (staticGrid.valid) {
    const { center: gridCenter, normal: gridNormal, basis } = staticGrid;

    sketch.push();
    sketch.translate(gridCenter.x, gridCenter.y, gridCenter.z);

    // Align local plane with world grid basis
    const R = basis.right;
    const U = basis.up;
    const N = gridNormal;

    sketch.applyMatrix(
      R.x,
      U.x,
      N.x,
      0,
      R.y,
      U.y,
      N.y,
      0,
      R.z,
      U.z,
      N.z,
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
    sketch.strokeWeight(GLOBAL_STROKE_WEIGHT);
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
    const normal = p5.Vector.sub(hitPoint, sphereCenter).normalize();
    const reflectedDir = reflect(laserRay.direction, normal);

    // Incident
    sketch.push();
    sketch.stroke(...INCIDENT_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      laserRay.origin.x,
      laserRay.origin.y,
      laserRay.origin.z,
      hitPoint.x,
      hitPoint.y,
      hitPoint.z
    );
    sketch.pop();

    // Normal
    const normalEnd = p5.Vector.add(
      hitPoint,
      p5.Vector.mult(normal, NORMAL_LENGTH)
    );
    sketch.push();
    sketch.stroke(...NORMAL_COLOR);
    sketch.strokeWeight(LASER_STROKE_WEIGHT);
    sketch.line(
      hitPoint.x,
      hitPoint.y,
      hitPoint.z,
      normalEnd.x,
      normalEnd.y,
      normalEnd.z
    );
    sketch.pop();

    // Reflected (Dynamic)
    // Find where it hits static grid
    let reflectionEnd = p5.Vector.add(
      hitPoint,
      p5.Vector.mult(reflectedDir, REFLECTION_LENGTH)
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
          const diff = p5.Vector.sub(gridHitPoint, staticGrid.center);
          const uLocal = diff.dot(staticGrid.basis.right);
          const vLocal = diff.dot(staticGrid.basis.up);

          // Map to texture UV [0..1]
          // Grid ranges from -GRID_SIZE/2 to +GRID_SIZE/2
          const u = (uLocal + GRID_SIZE / 2) / GRID_SIZE;
          const v = (vLocal + GRID_SIZE / 2) / GRID_SIZE;

          if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
            // Map to array index
            const xInd = Math.floor(u * GLOW_RES);
            const yInd = Math.floor(v * GLOW_RES);

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
      hitPoint.x,
      hitPoint.y,
      hitPoint.z,
      reflectionEnd.x,
      reflectionEnd.y,
      reflectionEnd.z
    );
    sketch.pop();

    // Draw Dynamic Spot on Static Grid
    if (gridHitPoint) {
      sketch.push();
      sketch.translate(gridHitPoint.x, gridHitPoint.y, gridHitPoint.z);
      // Move slightly towards source to avoid z-fighting
      const spotOffset = p5.Vector.mult(reflectedDir, -2);
      sketch.translate(spotOffset.x, spotOffset.y, spotOffset.z);
      sketch.fill(255, 0, 0);
      sketch.noStroke();
      sketch.sphere(8); // Visible dot
      sketch.pop();
    }
  }

  // Laser source
  sketch.push();
  sketch.translate(laserRay.origin.x, laserRay.origin.y, laserRay.origin.z);
  sketch.fill(255, 255, 0);
  sketch.noStroke();
  sketch.sphere(8);
  sketch.push();
  const targetDir = laserRay.direction;
  sketch.rotateZ(Math.atan2(targetDir.y, targetDir.x) * (180 / Math.PI));
  sketch.rotateY(-Math.asin(targetDir.z) * (180 / Math.PI));
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
    sketch.strokeWeight(GLOBAL_STROKE_WEIGHT);
    sketch.stroke('#444444');
    sketch.ellipse(0, 0, currentRadius, currentRadius, ELLIPSE_SEGMENTS);
    sketch.pop();
  }

  // Latitude
  sketch.rotateX(90);
  for (let lat = -90; lat < 90; lat += INCR / 2) {
    sketch.push();
    sketch.strokeWeight(GLOBAL_STROKE_WEIGHT);
    sketch.stroke('#444444');
    let dz = (currentRadius * Math.cos((Math.PI / 90) * lat)) / 2;
    sketch.translate(0, 0, dz);
    let latRadius = Math.abs(currentRadius * Math.sin((Math.PI / 90) * lat));
    sketch.ellipse(0, 0, latRadius, latRadius, ELLIPSE_SEGMENTS);
    sketch.pop();
  }
  sketch.pop(); // End sphere

  return { laserRay, gridInfo: staticGrid };
};

// 4-Phase camera transition calculator
// Phase 1: Point camera at sphere center AND move to R1 (GRID radius)
// Phase 2: Move along great circle at constant R1 (looking at sphere center)
// Phase 3: Change radius from R1 to destination radius
// Phase 4: Point camera at destination target
const calculateTransitionState = (
  startPos,
  startLookAt,
  startFov,
  endPos,
  endLookAt,
  endFov,
  sphereCenter,
  transitionTime
) => {
  const phase1End = PHASE_DURATION;
  const phase2End = PHASE_DURATION * 2;
  const phase3End = PHASE_DURATION * 3;
  const phase4End = PHASE_DURATION * 4;

  // R1 is the GRID_VANTAGE_POINT radius (the larger radius for traveling)
  const startDir = p5.Vector.sub(startPos, sphereCenter).normalize();
  const startRadius = p5.Vector.dist(startPos, sphereCenter);
  const endDir = p5.Vector.sub(endPos, sphereCenter).normalize();
  const endRadius = p5.Vector.dist(endPos, sphereCenter);

  // Always use the larger radius (GRID_VANTAGE_POINT) for great circle travel
  const R1 = Math.max(startRadius, endRadius);

  let camPos, lookAt, fov;

  if (transitionTime < phase1End) {
    // Phase 1: Point camera at sphere center AND move radius to R1
    const t = easeInOutCubic(transitionTime / PHASE_DURATION);

    // Interpolate lookAt to sphere center
    lookAt = p5.Vector.lerp(startLookAt, sphereCenter, t);

    // Keep direction, interpolate radius to R1
    const currentRadius = lerp(startRadius, R1, t);
    camPos = p5.Vector.add(
      sphereCenter,
      p5.Vector.mult(startDir, currentRadius)
    );

    fov = startFov;
  } else if (transitionTime < phase2End) {
    // Phase 2: Move along great circle at constant radius R1
    const t = easeInOutCubic((transitionTime - phase1End) / PHASE_DURATION);

    // Spherical interpolation for direction, constant radius R1
    const currentDir = slerp(startDir, endDir, t);
    camPos = p5.Vector.add(sphereCenter, p5.Vector.mult(currentDir, R1));
    lookAt = sphereCenter.copy();
    fov = startFov;
  } else if (transitionTime < phase3End) {
    // Phase 3: Change radius from R1 to destination radius
    const t = easeInOutCubic((transitionTime - phase2End) / PHASE_DURATION);

    // Keep final direction, interpolate radius from R1 to endRadius
    const currentRadius = lerp(R1, endRadius, t);
    camPos = p5.Vector.add(sphereCenter, p5.Vector.mult(endDir, currentRadius));
    lookAt = sphereCenter.copy();
    fov = lerp(startFov, endFov, t);
  } else {
    // Phase 4: Point camera at destination target
    const t = easeInOutCubic((transitionTime - phase3End) / PHASE_DURATION);
    camPos = endPos.copy();
    lookAt = p5.Vector.lerp(sphereCenter, endLookAt, t);
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
    sketch.strokeWeight(GLOBAL_STROKE_WEIGHT);
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
    const sphereYOffset = Math.sin(sketch.frameCount * 0.05) * 0.1;
    const sceneInfo = drawScene(sketch, INITIAL_RADIUS, sphereYOffset, true);

    // Calculate time in animation cycle
    const elapsedTime = sketch.millis() / 1000;
    const cycleTime = elapsedTime % CYCLE_DURATION;

    // Calculate laser intersection point
    const sphereRadius = INITIAL_RADIUS / 2;
    const sphereCenter = new p5.Vector(0, sphereYOffset, 0);
    const laserRay = sceneInfo.laserRay;
    const intersection = raySphereIntersection(
      laserRay.origin,
      laserRay.direction,
      sphereRadius,
      sphereCenter
    );

    if (
      intersection &&
      intersection.hit &&
      sceneInfo.gridInfo &&
      sceneInfo.gridInfo.valid
    ) {
      const hitPoint = intersection.point;
      const gridCenter = sceneInfo.gridInfo.center;
      const gridNormal = sceneInfo.gridInfo.normal;
      const basis = sceneInfo.gridInfo.basis;

      // Calculate both vantage points
      // Grid vantage point
      const gridCamOffset = p5.Vector.add(
        p5.Vector.mult(gridNormal, -GRID_VANTAGE_POINT.distFromGrid),
        p5.Vector.add(
          p5.Vector.mult(basis.right, GRID_VANTAGE_POINT.offsetScale),
          p5.Vector.mult(basis.up, GRID_VANTAGE_POINT.offsetScale)
        )
      );
      const gridCamPos = p5.Vector.add(gridCenter, gridCamOffset);
      const gridLookAt = gridCenter.copy();

      // Reflection vantage point
      const reflectionCamPos = p5.Vector.add(
        hitPoint,
        p5.Vector.add(
          p5.Vector.mult(
            p5.Vector.mult(laserRay.direction, -1),
            REFLECTION_VANTAGE_POINT.camDistance
          ),
          new p5.Vector(0, REFLECTION_VANTAGE_POINT.upOffset, 0)
        )
      );
      const reflectionLookAt = hitPoint.copy();

      // Determine which transition we're in and calculate camera state
      let camPos, lookAt, fov;

      if (cycleTime < PAUSE_DURATION) {
        // At grid vantage point
        camPos = gridCamPos;
        lookAt = gridLookAt;
        fov = GRID_VANTAGE_POINT.fov;
      } else if (cycleTime < PAUSE_DURATION + TRANSITION_DURATION) {
        // Transitioning from grid to reflection
        const transitionTime = cycleTime - PAUSE_DURATION;
        [camPos, lookAt, fov] = calculateTransitionState(
          gridCamPos,
          gridLookAt,
          GRID_VANTAGE_POINT.fov,
          reflectionCamPos,
          reflectionLookAt,
          REFLECTION_VANTAGE_POINT.fov,
          sphereCenter,
          transitionTime
        );
      } else if (cycleTime < PAUSE_DURATION * 2 + TRANSITION_DURATION) {
        // At reflection vantage point
        camPos = reflectionCamPos;
        lookAt = reflectionLookAt;
        fov = REFLECTION_VANTAGE_POINT.fov;
      } else {
        // Transitioning from reflection back to grid
        const transitionTime =
          cycleTime - PAUSE_DURATION * 2 - TRANSITION_DURATION;
        [camPos, lookAt, fov] = calculateTransitionState(
          reflectionCamPos,
          reflectionLookAt,
          REFLECTION_VANTAGE_POINT.fov,
          gridCamPos,
          gridLookAt,
          GRID_VANTAGE_POINT.fov,
          sphereCenter,
          transitionTime
        );
      }

      cam.setPosition(camPos.x, camPos.y, camPos.z);
      cam.lookAt(lookAt.x, lookAt.y, lookAt.z);

      let aspect = sketch.width / sketch.height;
      cam.perspective((fov * Math.PI) / 180, aspect, 1, 10000);
    }
  };
};

let gridP5 = new p5(gridCanvas);
