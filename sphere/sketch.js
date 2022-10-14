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

// rotation rate of the sphere
const OMEGA = 1 / 400;

let slider;

const HEIGHT = 400;
const WIDTH = HEIGHT;

const topCanvas = (sketch) => {
  const doubleRuled = 5;
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
    for (var lon = 0; lon < 180; lon += INCR) {
      sketch.push();
      sketch.rotateY(lon);
      if (lon % doubleRuled == 0) {
        sketch.strokeWeight(2 * strokeWeight);
        sketch.stroke("black");
      } else {
        sketch.strokeWeight(strokeWeight);
        sketch.stroke("darkgrey");
      }
      sketch.ellipse(0, 0, r, r, 50);
      sketch.pop();
    }
    sketch.rotateX(90);
    for (var lat = -90; lat < 90; lat += INCR / 2) {
      sketch.push();
      if (lat % (doubleRuled / 2) == 0) {
        sketch.strokeWeight(2 * strokeWeight);
        sketch.stroke("black");
      } else {
        sketch.strokeWeight(strokeWeight);
        sketch.stroke("darkgrey");
      }
      // Move pen down z-axis, r*cos(lat)/2 units
      let dz = (r * Math.cos((Math.PI / 90) * lat)) / 2;
      sketch.translate(0, 0, dz);
      // Radius at this latitude is r*sin(lat)
      let latRadius = r * Math.sin((Math.PI / 90) * lat);
      sketch.ellipse(0, 0, latRadius, latRadius, 50);
      sketch.pop();
    }
    let camX = r * Math.cos(sketch.frameCount * OMEGA),
      camY = 0,
      camZ = r * Math.sin(sketch.frameCount * OMEGA);
    // camera points to sphere's center
    cam.setPosition(camX, camY, camZ);
    cam.lookAt(0, 0, 0);
    cam.perspective((60 * INITIAL_RADIUS) / r);
  };
};
let topP5 = new p5(topCanvas);

let radiusTxt = (radius) => {
  return sprintf("radius: %4.1f", radius);
};

let bottomCanvas = (sketch) => {
  const X = WIDTH * 0.05,
    Y = 50,
    W = 100,
    L = 20;
  const SLIDER_XPOS = WIDTH * 0.05,
    SLIDER_YPOS = 420;
  sketch.setup = () => {
    sketch.createCanvas(WIDTH, HEIGHT / 5);
    let sliderStep = 0.01;
    slider = sketch.createSlider(
      Math.log(MIN_RADIUS),
      Math.log(MAX_RADIUS),
      Math.log(INITIAL_RADIUS),
      sliderStep
    );
    slider.position(SLIDER_XPOS, SLIDER_YPOS);
    slider.style("width", WIDTH * 0.9 + "px");
    slider.fValue = () => {
      return Math.exp(slider.value());
    };
  };
  sketch.draw = () => {
    sketch.background("white");
    sketch.textSize(15);
    sketch.textAlign(sketch.LEFT, sketch.CENTER);
    let radius = slider.fValue();
    sketch.text(radiusTxt(radius), X, Y, W, L);
  };
};
let bottomP5 = new p5(bottomCanvas);
