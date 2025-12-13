"use strict";

const strokeWeight = 1;

// Initial Radius
const INITIAL_RADIUS = 200;
// Incr determines the number of lines of lat and lon
const INCR = 1.0;
// R is the radius, which is controlled by the slider
let r = INITIAL_RADIUS;
const MIN_RADIUS = 100,
  MAX_RADIUS = 6800;

// Cam is the camera
let cam;

// Rendering constants
const ELLIPSE_SEGMENTS = 50;
const CAMERA_FOV_DEGREES = 60;

// Laser constants
const INCIDENT_COLOR = [255, 0, 0]; // Red color for incident ray
const REFLECTED_COLOR = [255, 0, 0]; // Red color for reflected ray
const NORMAL_COLOR = [0, 255, 0]; // Green color for surface normal
const LASER_STROKE_WEIGHT = 3;
const NORMAL_LENGTH = 100; // Length of normal vector to draw
const REFLECTION_LENGTH = 1000; // Length of reflected beam to draw

let slider;
let cameraYSlider;
let cameraXSlider;

const HEIGHT = 800;
const WIDTH = HEIGHT;

const topCanvas = (sketch) => {
  const doubleRuled = 5;
  
  // Vector math helper functions
  const vec3 = {
    // Dot product
    dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
    // Subtract vectors
    sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
    // Add vectors
    add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
    // Multiply vector by scalar
    scale: (v, s) => [v[0] * s, v[1] * s, v[2] * s],
    // Vector length
    length: (v) => Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]),
    // Normalize vector
    normalize: (v) => {
      const len = vec3.length(v);
      return len > 0 ? vec3.scale(v, 1 / len) : [0, 0, 0];
    }
  };
  
  // Calculate ray-sphere intersection using quadratic formula
  // Returns { hit: boolean, point: [x, y, z], t: number } or null
  const raySphereIntersection = (rayOrigin, rayDir, sphereRadius, sphereCenter = [0, 0, 0]) => {
    // Ray equation: P(t) = O + tD
    // Sphere equation: |P - C|² = r² (where C is sphere center)
    // Substituting: |(O - C) + tD|² = r²
    // Results in quadratic: at² + bt + c = 0
    // where a = D·D, b = 2((O-C)·D), c = (O-C)·(O-C) - r²
    
    const originToCenter = vec3.sub(rayOrigin, sphereCenter);
    
    const a = vec3.dot(rayDir, rayDir);
    const b = 2 * vec3.dot(originToCenter, rayDir);
    const c = vec3.dot(originToCenter, originToCenter) - sphereRadius * sphereRadius;
    
    const discriminant = b * b - 4 * a * c;
    
    if (discriminant < 0) {
      return null; // No intersection
    }
    
    const sqrtDiscriminant = Math.sqrt(discriminant);
    const t1 = (-b - sqrtDiscriminant) / (2 * a);
    const t2 = (-b + sqrtDiscriminant) / (2 * a);
    
    // Use the closest positive intersection (first hit point)
    let t = null;
    if (t1 > 0) {
      t = t1;
    } else if (t2 > 0) {
      t = t2;
    } else {
      return null; // Both intersections are behind the ray origin
    }
    
    // Calculate intersection point
    const point = vec3.add(rayOrigin, vec3.scale(rayDir, t));
    return { hit: true, point, t };
  };
  
  // Calculate reflection vector using: R = I - 2(I·N)N
  // where I is incident direction, N is surface normal
  const reflect = (incident, normal) => {
    const dot = vec3.dot(incident, normal);
    return vec3.sub(incident, vec3.scale(normal, 2 * dot));
  };
  
  // Get laser start position and direction
  const getLaserRay = () => {
    // Laser is stationary in world space, positioned toward the right
    // This creates a fixed laser that the rotating camera view shows from different angles
    const startX = 400;  // Right side
    const startY = 0;    // Middle height
    const startZ = 200;  // Slightly forward
    
    // Aim at a point offset from center to create an angle of incidence
    const targetOffset = [-20, 10, -15]; // Offset from sphere center
    const start = [startX, startY, startZ];
    const direction = vec3.normalize(vec3.sub(targetOffset, start));
    
    return { origin: start, direction };
  };
  
  sketch.setup = () => {
    sketch.createCanvas(HEIGHT, WIDTH, sketch.WEBGL);
    sketch.angleMode(sketch.DEGREES);
    sketch.strokeWeight(strokeWeight);
    cam = sketch.createCamera();
  };
  sketch.draw = () => {
    if (slider) {
      r = Math.exp(slider.value());
    }
    sketch.background("white");
    
    // Animate sphere moving up and down
    const sphereYOffset = Math.sin(sketch.frameCount * 0.02) * 50; // Oscillates between -50 and 50
    const sphereCenter = [0, sphereYOffset, 0];
    
    // Draw laser beam with reflection (before any transformations)
    const laserRay = getLaserRay();
    // Note: ellipse() uses width/height, so actual sphere radius is r/2
    const sphereRadius = r / 2;
    const intersection = raySphereIntersection(laserRay.origin, laserRay.direction, sphereRadius, sphereCenter);
    
    if (intersection && intersection.hit) {
      const hitPoint = intersection.point;
      
      // Calculate surface normal at intersection point
      // For a sphere, normal = (point - center) / |point - center|
      const normal = vec3.normalize(vec3.sub(hitPoint, sphereCenter));
      
      // Calculate reflected direction
      const reflectedDir = reflect(laserRay.direction, normal);
      
      // Draw incident laser beam (RED - from start to intersection point)
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
      
      // Draw surface normal (GREEN - perpendicular to sphere surface)
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
      
      // Draw reflected laser beam (RED - bouncing off surface)
      const reflectionEnd = vec3.add(hitPoint, vec3.scale(reflectedDir, REFLECTION_LENGTH));
      
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
    }
    
    // Draw laser source indicator
    sketch.push();
    sketch.translate(laserRay.origin[0], laserRay.origin[1], laserRay.origin[2]);
    sketch.fill(255, 255, 0); // Yellow
    sketch.noStroke();
    sketch.sphere(8); // Slightly larger than before
    // Draw a small cone to show direction
    sketch.push();
    // Calculate rotation to point cone toward target
    const targetDir = laserRay.direction;
    sketch.rotateZ(Math.atan2(targetDir[1], targetDir[0]) * (180 / Math.PI));
    sketch.rotateY(-Math.asin(targetDir[2]) * (180 / Math.PI));
    sketch.fill(255, 200, 0); // Orange
    sketch.translate(10, 0, 0);
    sketch.rotateZ(90);
    sketch.cone(6, 15);
    sketch.pop();
    sketch.pop();
    
    // Draw longitude lines (with sphere offset)
    sketch.push();
    sketch.translate(0, sphereYOffset, 0);
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
      sketch.ellipse(0, 0, r, r, ELLIPSE_SEGMENTS);
      sketch.pop();
    }
    sketch.pop();
    
    // Draw latitude lines (with sphere offset)
    sketch.push();
    sketch.translate(0, sphereYOffset, 0);
    sketch.rotateX(90);
    for (let lat = -90; lat < 90; lat += INCR / 2) {
      sketch.push();
      // Use Math.abs to handle negative latitudes correctly
      if (Math.abs(lat) % (doubleRuled / 2) == 0) {
        sketch.strokeWeight(2 * strokeWeight);
        sketch.stroke("black");
      } else {
        sketch.strokeWeight(strokeWeight);
        sketch.stroke("darkgrey");
      }
      // Move pen down z-axis, r*cos(lat_radians)/2 units
      let dz = (r * Math.cos((Math.PI / 90) * lat)) / 2;
      sketch.translate(0, 0, dz);
      // Radius at this latitude is r*sin(lat_radians), use abs to ensure positive radius
      let latRadius = Math.abs(r * Math.sin((Math.PI / 90) * lat));
      sketch.ellipse(0, 0, latRadius, latRadius, ELLIPSE_SEGMENTS);
      sketch.pop();
    }
    sketch.pop();
    
    // Fixed camera position (no rotation)
    let camX = cameraXSlider ? cameraXSlider.value() : r * 1.5,  // Controlled by slider
      camY = cameraYSlider ? cameraYSlider.value() : r * 0.3,    // Controlled by slider
      camZ = r * 0.8;    // Forward
    // camera points to sphere's center
    cam.setPosition(camX, camY, camZ);
    cam.lookAt(0, 0, 0);
    cam.perspective((CAMERA_FOV_DEGREES * INITIAL_RADIUS) / r);
  };
};
let topP5 = new p5(topCanvas);

let radiusTxt = (radius) => {
  // Format radius with 1 decimal place, padding to 4 characters
  return `radius: ${radius.toFixed(1).padStart(4)}`;
};

let bottomCanvas = (sketch) => {
  const X = WIDTH * 0.05,
    Y = 50,
    W = 200,
    L = 20;
  const SLIDER_XPOS = WIDTH * 0.05,
    SLIDER_YPOS = HEIGHT + 20;
  const CAMERA_Y_SLIDER_YPOS = HEIGHT + 60;
  const CAMERA_X_SLIDER_YPOS = HEIGHT + 100;
  const SLIDER_WIDTH_FACTOR = 0.9;
  sketch.setup = () => {
    sketch.createCanvas(WIDTH, HEIGHT / 5);
    let sliderStep = 0.01;
    
    // Radius slider
    slider = sketch.createSlider(
      Math.log(MIN_RADIUS),
      Math.log(MAX_RADIUS),
      Math.log(INITIAL_RADIUS),
      sliderStep
    );
    slider.position(SLIDER_XPOS, SLIDER_YPOS);
    slider.style("width", WIDTH * SLIDER_WIDTH_FACTOR + "px");
    slider.fValue = () => {
      return Math.exp(slider.value());
    };
    
    // Camera Y slider
    cameraYSlider = sketch.createSlider(-400, 400, 60, 1);
    cameraYSlider.position(SLIDER_XPOS, CAMERA_Y_SLIDER_YPOS);
    cameraYSlider.style("width", WIDTH * SLIDER_WIDTH_FACTOR + "px");
    
    // Camera X slider
    cameraXSlider = sketch.createSlider(-800, 800, 300, 1);
    cameraXSlider.position(SLIDER_XPOS, CAMERA_X_SLIDER_YPOS);
    cameraXSlider.style("width", WIDTH * SLIDER_WIDTH_FACTOR + "px");
  };
  sketch.draw = () => {
    sketch.background("white");
    sketch.textSize(15);
    sketch.textAlign(sketch.LEFT, sketch.CENTER);
    let radius = slider.fValue();
    sketch.text(radiusTxt(radius), X, Y, W, L);
    
    let camYValue = cameraYSlider.value();
    sketch.text(`camera Y: ${camYValue}`, X, Y + 40, W, L);
    
    let camXValue = cameraXSlider.value();
    sketch.text(`camera X: ${camXValue}`, X, Y + 80, W, L);
  };
};
let bottomP5 = new p5(bottomCanvas);
