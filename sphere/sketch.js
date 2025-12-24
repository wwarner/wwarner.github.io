'use strict';

const strokeWeight = 1;

// Initial Radius
const INITIAL_RADIUS = 50;
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

// Glow Constants
const GLOW_RES = 64;
const GLOW_DECAY = 0.94;
let glowMap = new Float32Array(GLOW_RES * GLOW_RES);

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

    // Lazy initialize texture using createImage (avoids Chrome/iOS iframe issues with createGraphics)
    if (!sketch.gridTexture) {
      // Use smaller texture for iOS compatibility
      const texSize = 256;
      const img = sketch.createImage(texSize, texSize);
      img.loadPixels();

      const step = texSize / GRID_RES;
      const lineHalfWidth = 1;
      const axisHalfWidth = 2;

      // Fill with background color (light blue, fully opaque for iOS)
      for (let i = 0; i < img.pixels.length; i += 4) {
        img.pixels[i] = 220;     // R
        img.pixels[i + 1] = 220; // G
        img.pixels[i + 2] = 255; // B
        img.pixels[i + 3] = 255; // A (fully opaque for iOS)
      }

      // Draw non-axis grid lines first (use gray color for iOS - no alpha)
      for (let g = 0; g <= GRID_RES; g++) {
        if (g === GRID_RES / 2) continue; // Skip axes for now
        const pos = Math.floor(g * step);
        const hw = lineHalfWidth;

        // Horizontal and vertical lines
        for (let offset = -hw; offset <= hw; offset++) {
          // Vertical line at x = pos
          for (let y = 0; y < texSize; y++) {
            const x = pos + offset;
            if (x >= 0 && x < texSize) {
              const idx = (y * texSize + x) * 4;
              img.pixels[idx] = 150;     // Gray instead of black with alpha
              img.pixels[idx + 1] = 150;
              img.pixels[idx + 2] = 150;
              img.pixels[idx + 3] = 255; // Fully opaque
            }
          }
          // Horizontal line at y = pos
          for (let x = 0; x < texSize; x++) {
            const y = pos + offset;
            if (y >= 0 && y < texSize) {
              const idx = (y * texSize + x) * 4;
              img.pixels[idx] = 150;     // Gray instead of black with alpha
              img.pixels[idx + 1] = 150;
              img.pixels[idx + 2] = 150;
              img.pixels[idx + 3] = 255; // Fully opaque
            }
          }
        }
      }

      // Draw axes on top (solid black, not interrupted)
      const axisPos = Math.floor((GRID_RES / 2) * step);
      for (let offset = -axisHalfWidth; offset <= axisHalfWidth; offset++) {
        // Y axis (vertical line at x = axisPos)
        for (let y = 0; y < texSize; y++) {
          const x = axisPos + offset;
          if (x >= 0 && x < texSize) {
            const idx = (y * texSize + x) * 4;
            img.pixels[idx] = 0;
            img.pixels[idx + 1] = 0;
            img.pixels[idx + 2] = 0;
            img.pixels[idx + 3] = 255;
          }
        }
        // X axis (horizontal line at y = axisPos)
        for (let x = 0; x < texSize; x++) {
          const y = axisPos + offset;
          if (y >= 0 && y < texSize) {
            const idx = (y * texSize + x) * 4;
            img.pixels[idx] = 0;
            img.pixels[idx + 1] = 0;
            img.pixels[idx + 2] = 0;
            img.pixels[idx + 3] = 255;
          }
        }
      }

      img.updatePixels();
      sketch.gridTexture = img;
    }

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

    // Static Grid Plane
    sketch.push();
    sketch.translate(0, 0, -1);
    sketch.texture(sketch.gridTexture);
    sketch.noStroke();
    sketch.plane(GRID_SIZE, GRID_SIZE);
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
      sketch.strokeWeight(2 * strokeWeight);
      sketch.stroke('black');
    } else {
      sketch.strokeWeight(strokeWeight);
      sketch.stroke('darkgrey');
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
      sketch.stroke('black');
    } else {
      sketch.strokeWeight(strokeWeight);
      sketch.stroke('darkgrey');
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

    // Fix for Chrome WebGL issues in iframes with CSS transforms
    sketch.setAttributes('preserveDrawingBuffer', true);

    let canvas = sketch.createCanvas(w, h, sketch.WEBGL);
    canvas.parent('canvas-container');
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
    const sphereYOffset = Math.sin(sketch.frameCount * freq) * 1;
    // drawScene returns the static grid info as gridInfo because we changed implementation
    const sceneInfo = drawScene(sketch, r, sphereYOffset, true);

    if (sceneInfo.gridInfo && sceneInfo.gridInfo.valid) {
      const gridCenter = sceneInfo.gridInfo.center;
      const gridNormal = sceneInfo.gridInfo.normal;

      const distFromGrid = 3000;

      const basis = sceneInfo.gridInfo.basis;
      const offsetScale = 200;

      const camOffset = vec3.add(
        vec3.scale(gridNormal, -distFromGrid),
        vec3.add(
          vec3.scale(basis.right, offsetScale),
          vec3.scale(basis.up, offsetScale)
        )
      );

      const camPos = vec3.add(gridCenter, camOffset);

      cam.setPosition(camPos[0], camPos[1], camPos[2]);
      cam.lookAt(gridCenter[0], gridCenter[1], gridCenter[2]);
    }
    cam.perspective((600 * Math.PI) / 180);
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
