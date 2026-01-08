// NPZD Model Parameters (defaults)
let defaultParams = {
    vm: 1.5,      // Maximum uptake rate
    kn: 0.5,      // Half-saturation constant
    i0: 1.0,      // Light intensity
    rm: 1.0,      // Maximum grazing rate
    lambda: 1.0,  // Ivlev grazing constant
    alpha: 0.3,   // Fraction to nutrients
    beta: 0.3,    // Fraction to zooplankton
    epsilon: 0.05,// Phytoplankton mortality
    g: 0.1,       // Zooplankton mortality
    r: 0.05,      // Phytoplankton sinking
    phi: 0.1      // Detritus remineralization
};

// Current parameters
let params = {...defaultParams};

// State variables
let N = 10.0;  // Nutrients
let P = 1.0;   // Phytoplankton
let Z = 0.5;   // Zooplankton
let D = 0.5;   // Detritus

// Initial conditions (for reset)
let N0 = 10.0;
let P0 = 1.0;
let Z0 = 0.5;
let D0 = 0.5;

// Simulation parameters
let dt = 0.01;  // Time step
let time = 0;
let running = false;
let vizMode = 'timeseries';  // 'timeseries' or 'flow'
let speedMultiplier = 1;  // Speed multiplier (1x, 2x, or 4x)

// History for plotting
let history = {
    t: [],
    N: [],
    P: [],
    Z: [],
    D: []
};

let maxHistoryLength = 1000;

// Canvas dimensions
let canvasWidth = 1200;
let canvasHeight = 800;
let margin = { top: 40, right: 40, bottom: 60, left: 60 };
let mainPlotHeight = 450;  // Height of main plot
let functionalPlotHeight = 250;  // Height of functional form plots
let functionalPlotY = mainPlotHeight + 50;  // Y position where functional plots start

function setup() {
    let canvas = createCanvas(canvasWidth, canvasHeight);
    canvas.parent('canvas-container');

    // Initialize history
    history.t.push(time);
    history.N.push(N);
    history.P.push(P);
    history.Z.push(Z);
    history.D.push(D);

    // Setup parameter controls
    setupControls();

    frameRate(30);
}

function draw() {
    background(255);

    if (running) {
        // Run multiple steps per frame for faster simulation
        let stepsPerFrame = 5 * speedMultiplier;
        for (let i = 0; i < stepsPerFrame; i++) {
            rk4Step();
            time += dt;

            // Store history
            history.t.push(time);
            history.N.push(N);
            history.P.push(P);
            history.Z.push(Z);
            history.D.push(D);

            // Limit history length
            if (history.t.length > maxHistoryLength) {
                history.t.shift();
                history.N.shift();
                history.P.shift();
                history.Z.shift();
                history.D.shift();
            }
        }
    }

    // Draw the appropriate visualization
    if (vizMode === 'timeseries') {
        drawTimeSeries();
        displayCurrentValues();
    } else {
        drawFlowDiagram();
    }

    // Draw functional form plots below
    drawFunctionalForms();
}

// Light limitation function (simple linear for now)
function f_light(I) {
    return I;
}

// Differential equations
function dNdt(N, P, Z, D) {
    let uptake = params.vm * (N / (params.kn + N)) * f_light(params.i0) * P;
    let grazing = params.rm * (1 - Math.exp(-params.lambda * P)) * Z;

    return -uptake + params.alpha * grazing + params.epsilon * P + params.g * Z + params.phi * D;
}

function dPdt(N, P, Z, D) {
    let uptake = params.vm * (N / (params.kn + N)) * f_light(params.i0) * P;
    let grazing = params.rm * (1 - Math.exp(-params.lambda * P)) * Z;

    return uptake - grazing - params.epsilon * P - params.r * P;
}

function dZdt(N, P, Z, D) {
    let grazing = params.rm * (1 - Math.exp(-params.lambda * P)) * Z;

    return params.beta * grazing - params.g * Z;
}

function dDdt(N, P, Z, D) {
    let grazing = params.rm * (1 - Math.exp(-params.lambda * P)) * Z;

    return params.r * P + (1 - params.alpha - params.beta) * grazing - params.phi * D;
}

// Runge-Kutta 4th order integration
function rk4Step() {
    // k1
    let k1_N = dNdt(N, P, Z, D);
    let k1_P = dPdt(N, P, Z, D);
    let k1_Z = dZdt(N, P, Z, D);
    let k1_D = dDdt(N, P, Z, D);

    // k2
    let k2_N = dNdt(N + 0.5 * dt * k1_N, P + 0.5 * dt * k1_P, Z + 0.5 * dt * k1_Z, D + 0.5 * dt * k1_D);
    let k2_P = dPdt(N + 0.5 * dt * k1_N, P + 0.5 * dt * k1_P, Z + 0.5 * dt * k1_Z, D + 0.5 * dt * k1_D);
    let k2_Z = dZdt(N + 0.5 * dt * k1_N, P + 0.5 * dt * k1_P, Z + 0.5 * dt * k1_Z, D + 0.5 * dt * k1_D);
    let k2_D = dDdt(N + 0.5 * dt * k1_N, P + 0.5 * dt * k1_P, Z + 0.5 * dt * k1_Z, D + 0.5 * dt * k1_D);

    // k3
    let k3_N = dNdt(N + 0.5 * dt * k2_N, P + 0.5 * dt * k2_P, Z + 0.5 * dt * k2_Z, D + 0.5 * dt * k2_D);
    let k3_P = dPdt(N + 0.5 * dt * k2_N, P + 0.5 * dt * k2_P, Z + 0.5 * dt * k2_Z, D + 0.5 * dt * k2_D);
    let k3_Z = dZdt(N + 0.5 * dt * k2_N, P + 0.5 * dt * k2_P, Z + 0.5 * dt * k2_Z, D + 0.5 * dt * k2_D);
    let k3_D = dDdt(N + 0.5 * dt * k2_N, P + 0.5 * dt * k2_P, Z + 0.5 * dt * k2_Z, D + 0.5 * dt * k2_D);

    // k4
    let k4_N = dNdt(N + dt * k3_N, P + dt * k3_P, Z + dt * k3_Z, D + dt * k3_D);
    let k4_P = dPdt(N + dt * k3_N, P + dt * k3_P, Z + dt * k3_Z, D + dt * k3_D);
    let k4_Z = dZdt(N + dt * k3_N, P + dt * k3_P, Z + dt * k3_Z, D + dt * k3_D);
    let k4_D = dDdt(N + dt * k3_N, P + dt * k3_P, Z + dt * k3_Z, D + dt * k3_D);

    // Update
    N += (dt / 6.0) * (k1_N + 2 * k2_N + 2 * k3_N + k4_N);
    P += (dt / 6.0) * (k1_P + 2 * k2_P + 2 * k3_P + k4_P);
    Z += (dt / 6.0) * (k1_Z + 2 * k2_Z + 2 * k3_Z + k4_Z);
    D += (dt / 6.0) * (k1_D + 2 * k2_D + 2 * k3_D + k4_D);

    // Ensure non-negative values
    N = Math.max(0, N);
    P = Math.max(0, P);
    Z = Math.max(0, Z);
    D = Math.max(0, D);
}

function drawTimeSeries() {
    if (history.t.length < 2) return;

    // Calculate plot dimensions
    let plotWidth = canvasWidth - margin.left - margin.right;
    let plotHeight = mainPlotHeight - margin.top - margin.bottom;

    // Find min/max for scaling
    let allValues = [...history.N, ...history.P, ...history.Z, ...history.D];
    let minVal = 0;
    let maxVal = Math.max(...allValues) * 1.1;

    let minTime = Math.min(...history.t);
    let maxTime = Math.max(...history.t);

    // Draw axes
    stroke(200);
    strokeWeight(1);
    line(margin.left, margin.top, margin.left, mainPlotHeight - margin.bottom);
    line(margin.left, mainPlotHeight - margin.bottom, canvasWidth - margin.right, mainPlotHeight - margin.bottom);

    // Draw grid
    stroke(240);
    for (let i = 0; i <= 5; i++) {
        let y = map(i, 0, 5, mainPlotHeight - margin.bottom, margin.top);
        line(margin.left, y, canvasWidth - margin.right, y);
    }

    // Draw time series
    drawLine(history.t, history.N, minTime, maxTime, minVal, maxVal, color(0, 100, 200), 'N');
    drawLine(history.t, history.P, minTime, maxTime, minVal, maxVal, color(0, 200, 0), 'P');
    drawLine(history.t, history.Z, minTime, maxTime, minVal, maxVal, color(200, 0, 0), 'Z');
    drawLine(history.t, history.D, minTime, maxTime, minVal, maxVal, color(150, 75, 0), 'D');

    // Draw axis labels
    fill(0);
    noStroke();
    textAlign(CENTER);
    textSize(14);
    text('Time', canvasWidth / 2, mainPlotHeight - 20);

    push();
    translate(15, mainPlotHeight / 2);
    rotate(-PI / 2);
    text('Concentration', 0, 0);
    pop();

    // Draw title
    textSize(16);
    textAlign(CENTER);
    text('NPZD Model Time Series', canvasWidth / 2, 20);

    // Draw y-axis labels
    textSize(11);
    textAlign(RIGHT);
    for (let i = 0; i <= 5; i++) {
        let y = map(i, 0, 5, mainPlotHeight - margin.bottom, margin.top);
        let val = map(i, 0, 5, minVal, maxVal);
        text(val.toFixed(2), margin.left - 10, y + 4);
    }

    // Draw x-axis labels
    textAlign(CENTER);
    for (let i = 0; i <= 5; i++) {
        let x = map(i, 0, 5, margin.left, canvasWidth - margin.right);
        let val = map(i, 0, 5, minTime, maxTime);
        text(val.toFixed(1), x, mainPlotHeight - margin.bottom + 20);
    }

    // Draw legend
    drawLegend();
}

function drawLine(timeData, valueData, minTime, maxTime, minVal, maxVal, col, label) {
    stroke(col);
    strokeWeight(2);
    noFill();

    beginShape();
    for (let i = 0; i < timeData.length; i++) {
        let x = map(timeData[i], minTime, maxTime, margin.left, canvasWidth - margin.right);
        let y = map(valueData[i], minVal, maxVal, mainPlotHeight - margin.bottom, margin.top);
        vertex(x, y);
    }
    endShape();
}

function drawLegend() {
    let legendX = canvasWidth - margin.right - 100;
    let legendY = margin.top + 20;
    let lineLength = 30;
    let spacing = 25;

    textAlign(LEFT);
    textSize(12);

    // N
    stroke(0, 100, 200);
    strokeWeight(2);
    line(legendX, legendY, legendX + lineLength, legendY);
    fill(0);
    noStroke();
    text('Nutrients', legendX + lineLength + 5, legendY + 4);

    // P
    legendY += spacing;
    stroke(0, 200, 0);
    strokeWeight(2);
    line(legendX, legendY, legendX + lineLength, legendY);
    fill(0);
    noStroke();
    text('Phytoplankton', legendX + lineLength + 5, legendY + 4);

    // Z
    legendY += spacing;
    stroke(200, 0, 0);
    strokeWeight(2);
    line(legendX, legendY, legendX + lineLength, legendY);
    fill(0);
    noStroke();
    text('Zooplankton', legendX + lineLength + 5, legendY + 4);

    // D
    legendY += spacing;
    stroke(150, 75, 0);
    strokeWeight(2);
    line(legendX, legendY, legendX + lineLength, legendY);
    fill(0);
    noStroke();
    text('Detritus', legendX + lineLength + 5, legendY + 4);
}

function displayCurrentValues() {
    let x = margin.left + 20;
    let y = mainPlotHeight - margin.bottom - 80;

    fill(255, 255, 255, 200);
    stroke(200);
    strokeWeight(1);
    rect(x - 10, y - 20, 200, 90, 5);

    fill(0);
    noStroke();
    textAlign(LEFT);
    textSize(12);
    text(`Time: ${time.toFixed(2)}`, x, y);
    text(`N: ${N.toFixed(4)}`, x, y + 20);
    text(`P: ${P.toFixed(4)}`, x, y + 40);
    text(`Z: ${Z.toFixed(4)}`, x, y + 60);
    text(`D: ${D.toFixed(4)}`, x + 100, y + 20);
    text(`Total: ${(N + P + Z + D).toFixed(4)}`, x + 100, y + 40);
}

function drawFlowDiagram() {
    // Calculate all flow rates
    let uptake = params.vm * (N / (params.kn + N)) * f_light(params.i0) * P;  // N → P
    let grazing = params.rm * (1 - Math.exp(-params.lambda * P)) * Z;  // Total grazing
    let P_mortality = params.epsilon * P;  // P → N
    let P_sinking = params.r * P;  // P → D
    let Z_excretion = params.alpha * grazing;  // Z → N (from grazing)
    let Z_mortality = params.g * Z;  // Z → N (mortality)
    let Z_assimilation = params.beta * grazing;  // grazing → Z
    let grazing_to_D = (1 - params.alpha - params.beta) * grazing;  // grazing → D
    let remineralization = params.phi * D;  // D → N

    // Node positions (in a roughly square layout)
    let centerX = canvasWidth / 2;
    let centerY = mainPlotHeight / 2;
    let spacing = 150;

    let positions = {
        N: { x: centerX - spacing, y: centerY - spacing },  // Top left
        P: { x: centerX + spacing, y: centerY - spacing },  // Top right
        Z: { x: centerX + spacing, y: centerY + spacing },  // Bottom right
        D: { x: centerX - spacing, y: centerY + spacing }   // Bottom left
    };

    // Colors matching time series
    let colors = {
        N: color(0, 100, 200),
        P: color(0, 200, 0),
        Z: color(200, 0, 0),
        D: color(150, 75, 0)
    };

    // Find max concentration for scaling circles
    let maxConc = Math.max(N, P, Z, D);
    let minRadius = 30;
    let maxRadius = 100;

    // Find max flow rate for scaling arrows
    let allFlows = [uptake, P_mortality, P_sinking, Z_excretion + Z_mortality,
                    Z_assimilation, grazing_to_D, remineralization];
    let maxFlow = Math.max(...allFlows, 0.001);  // Avoid division by zero

    // Draw title
    fill(0);
    noStroke();
    textAlign(CENTER);
    textSize(16);
    text('NPZD Flow Diagram', canvasWidth / 2, 25);

    // Draw arrows (behind circles)
    drawFlowArrow(positions.N, positions.P, uptake, maxFlow, 'N→P (uptake)');
    drawFlowArrow(positions.P, positions.N, P_mortality, maxFlow, 'P→N (mort)');
    drawFlowArrow(positions.P, positions.D, P_sinking, maxFlow, 'P→D (sink)');
    drawFlowArrow(positions.P, positions.Z, Z_assimilation, maxFlow, 'P→Z (graz)');
    drawFlowArrow(positions.Z, positions.N, Z_excretion + Z_mortality, maxFlow, 'Z→N');
    drawFlowArrow(positions.Z, positions.D, grazing_to_D, maxFlow, 'Z→D');
    drawFlowArrow(positions.D, positions.N, remineralization, maxFlow, 'D→N (remin)');

    // Draw circles for each compartment
    drawCompartment(positions.N.x, positions.N.y, N, maxConc, minRadius, maxRadius, colors.N, 'N');
    drawCompartment(positions.P.x, positions.P.y, P, maxConc, minRadius, maxRadius, colors.P, 'P');
    drawCompartment(positions.Z.x, positions.Z.y, Z, maxConc, minRadius, maxRadius, colors.Z, 'Z');
    drawCompartment(positions.D.x, positions.D.y, D, maxConc, minRadius, maxRadius, colors.D, 'D');

    // Display current values and time
    displayFlowInfo();
}

function drawCompartment(x, y, concentration, maxConc, minRadius, maxRadius, col, label) {
    // Calculate radius based on concentration
    let radius = map(concentration, 0, maxConc, minRadius, maxRadius);

    // Draw circle
    fill(col);
    stroke(0);
    strokeWeight(2);
    circle(x, y, radius * 2);

    // Draw label
    fill(255);
    noStroke();
    textAlign(CENTER, CENTER);
    textSize(18);
    text(label, x, y - 10);
    textSize(14);
    text(concentration.toFixed(3), x, y + 10);
}

function drawFlowArrow(from, to, flowRate, maxFlow, label) {
    if (flowRate <= 0) return;  // Don't draw zero or negative flows

    // Calculate arrow thickness based on flow rate
    let minThickness = 1;
    let maxThickness = 15;
    let thickness = map(flowRate, 0, maxFlow, minThickness, maxThickness);

    // Arrow color - semi-transparent black
    stroke(0, 0, 0, 150);
    strokeWeight(thickness);

    // Calculate arrow positions (offset from center of circles)
    let angle = Math.atan2(to.y - from.y, to.x - from.x);
    let offset = 50;  // Approximate radius offset to start/end outside circles

    let startX = from.x + offset * Math.cos(angle);
    let startY = from.y + offset * Math.sin(angle);
    let endX = to.x - offset * Math.cos(angle);
    let endY = to.y - offset * Math.sin(angle);

    // Draw line
    line(startX, startY, endX, endY);

    // Draw arrowhead
    let arrowSize = thickness + 5;
    let arrowAngle = Math.PI / 6;  // 30 degrees

    push();
    translate(endX, endY);
    rotate(angle);
    fill(0, 0, 0, 150);
    noStroke();
    triangle(0, 0, -arrowSize, -arrowSize/2, -arrowSize, arrowSize/2);
    pop();

    // Draw flow rate label (optional, for major flows)
    if (flowRate > maxFlow * 0.1) {  // Only label significant flows
        let midX = (startX + endX) / 2;
        let midY = (startY + endY) / 2;

        fill(255, 255, 255, 220);
        noStroke();
        let boxWidth = 60;
        let boxHeight = 18;
        rect(midX - boxWidth/2, midY - boxHeight/2, boxWidth, boxHeight, 3);

        fill(0);
        textAlign(CENTER, CENTER);
        textSize(10);
        text(flowRate.toFixed(3), midX, midY);
    }
}

function displayFlowInfo() {
    let x = 20;
    let y = mainPlotHeight - 120;

    fill(255, 255, 255, 230);
    stroke(200);
    strokeWeight(1);
    rect(x - 10, y - 20, 280, 100, 5);

    fill(0);
    noStroke();
    textAlign(LEFT);
    textSize(12);
    text(`Time: ${time.toFixed(2)}`, x, y);
    text(`Total mass: ${(N + P + Z + D).toFixed(4)}`, x, y + 20);
    text(`Arrow thickness = flow rate`, x, y + 40);
    text(`Circle diameter = concentration`, x, y + 60);
}

function drawFunctionalForms() {
    // Plot dimensions
    let plotWidth = (canvasWidth - 3 * margin.left) / 2;
    let plotHeight = functionalPlotHeight - margin.top - margin.bottom;

    // Left plot: Michaelis-Menten
    let mmX = margin.left;
    let mmY = functionalPlotY;
    drawMichaelisMenten(mmX, mmY, plotWidth, plotHeight);

    // Right plot: Ivlev grazing
    let ivlevX = canvasWidth / 2 + margin.left / 2;
    let ivlevY = functionalPlotY;
    drawIvlevGrazing(ivlevX, ivlevY, plotWidth, plotHeight);
}

function drawMichaelisMenten(x0, y0, w, h) {
    // Title
    fill(0);
    noStroke();
    textAlign(CENTER);
    textSize(14);
    text('Michaelis-Menten Nutrient Limitation', x0 + w/2, y0 - 5);

    // Axes
    stroke(200);
    strokeWeight(1);
    line(x0, y0, x0, y0 + h);
    line(x0, y0 + h, x0 + w, y0 + h);

    // Grid
    stroke(240);
    for (let i = 0; i <= 5; i++) {
        let y = map(i, 0, 5, y0 + h, y0);
        line(x0, y, x0 + w, y);
    }

    // Plot the function N/(K_N + N)
    let numPoints = 100;
    let maxN = 5.0;  // Maximum N to plot

    stroke(0, 100, 200);
    strokeWeight(2);
    noFill();
    beginShape();
    for (let i = 0; i <= numPoints; i++) {
        let n = map(i, 0, numPoints, 0, maxN);
        let uptakeFactor = n / (params.kn + n);
        let px = map(n, 0, maxN, x0, x0 + w);
        let py = map(uptakeFactor, 0, 1, y0 + h, y0);
        vertex(px, py);
    }
    endShape();

    // Mark current K_N on the plot
    let kn_x = map(params.kn, 0, maxN, x0, x0 + w);
    let kn_y = map(0.5, 0, 1, y0 + h, y0);  // At K_N, uptake = 0.5
    fill(200, 0, 0);
    noStroke();
    circle(kn_x, kn_y, 8);

    // Mark current N on the plot
    let currentUptake = N / (params.kn + N);
    let n_x = map(N, 0, maxN, x0, x0 + w);
    let n_y = map(currentUptake, 0, 1, y0 + h, y0);
    fill(0, 200, 0);
    circle(n_x, n_y, 8);

    // Axes labels
    fill(0);
    textSize(12);
    textAlign(CENTER);
    text('N (Nutrient)', x0 + w/2, y0 + h + 25);

    push();
    translate(x0 - 35, y0 + h/2);
    rotate(-PI / 2);
    text('N/(K_N + N)', 0, 0);
    pop();

    // Y-axis values
    textSize(10);
    textAlign(RIGHT);
    for (let i = 0; i <= 5; i++) {
        let val = i / 5.0;
        let y = map(val, 0, 1, y0 + h, y0);
        text(val.toFixed(1), x0 - 5, y + 3);
    }

    // X-axis values
    textAlign(CENTER);
    for (let i = 0; i <= 5; i++) {
        let val = map(i, 0, 5, 0, maxN);
        let x = map(i, 0, 5, x0, x0 + w);
        text(val.toFixed(1), x, y0 + h + 15);
    }

    // Legend
    textSize(10);
    textAlign(LEFT);
    fill(200, 0, 0);
    circle(x0 + 10, y0 + 15, 8);
    fill(0);
    text(`K_N = ${params.kn.toFixed(2)}`, x0 + 20, y0 + 18);

    fill(0, 200, 0);
    circle(x0 + 10, y0 + 30, 8);
    fill(0);
    text(`N = ${N.toFixed(2)}`, x0 + 20, y0 + 33);
}

function drawIvlevGrazing(x0, y0, w, h) {
    // Title
    fill(0);
    noStroke();
    textAlign(CENTER);
    textSize(14);
    text('Ivlev Grazing Function', x0 + w/2, y0 - 5);

    // Axes
    stroke(200);
    strokeWeight(1);
    line(x0, y0, x0, y0 + h);
    line(x0, y0 + h, x0 + w, y0 + h);

    // Grid
    stroke(240);
    for (let i = 0; i <= 5; i++) {
        let y = map(i, 0, 5, y0 + h, y0);
        line(x0, y, x0 + w, y);
    }

    // Plot the function 1 - e^(-λP)
    let numPoints = 100;
    let maxP = 5.0;  // Maximum P to plot

    stroke(200, 0, 0);
    strokeWeight(2);
    noFill();
    beginShape();
    for (let i = 0; i <= numPoints; i++) {
        let p = map(i, 0, numPoints, 0, maxP);
        let grazingFactor = 1 - Math.exp(-params.lambda * p);
        let px = map(p, 0, maxP, x0, x0 + w);
        let py = map(grazingFactor, 0, 1, y0 + h, y0);
        vertex(px, py);
    }
    endShape();

    // Mark current P on the plot
    let currentGrazing = 1 - Math.exp(-params.lambda * P);
    let p_x = map(P, 0, maxP, x0, x0 + w);
    let p_y = map(currentGrazing, 0, 1, y0 + h, y0);
    fill(0, 200, 0);
    noStroke();
    circle(p_x, p_y, 8);

    // Axes labels
    fill(0);
    textSize(12);
    textAlign(CENTER);
    text('P (Phytoplankton)', x0 + w/2, y0 + h + 25);

    push();
    translate(x0 - 35, y0 + h/2);
    rotate(-PI / 2);
    text('1 - e^(-λP)', 0, 0);
    pop();

    // Y-axis values
    textSize(10);
    textAlign(RIGHT);
    for (let i = 0; i <= 5; i++) {
        let val = i / 5.0;
        let y = map(val, 0, 1, y0 + h, y0);
        text(val.toFixed(1), x0 - 5, y + 3);
    }

    // X-axis values
    textAlign(CENTER);
    for (let i = 0; i <= 5; i++) {
        let val = map(i, 0, 5, 0, maxP);
        let x = map(i, 0, 5, x0, x0 + w);
        text(val.toFixed(1), x, y0 + h + 15);
    }

    // Legend
    textSize(10);
    textAlign(LEFT);
    fill(0, 200, 0);
    circle(x0 + 10, y0 + 15, 8);
    fill(0);
    text(`P = ${P.toFixed(2)}`, x0 + 20, y0 + 18);

    fill(0);
    text(`λ = ${params.lambda.toFixed(2)}`, x0 + 10, y0 + 33);
}

function setupControls() {
    // Update parameter displays and values
    const controls = ['vm', 'kn', 'i0', 'rm', 'lambda', 'alpha', 'beta', 'epsilon', 'g', 'r', 'phi'];

    controls.forEach(param => {
        let slider = document.getElementById(param);
        let display = document.getElementById(param + '-val');

        slider.addEventListener('input', function() {
            params[param] = parseFloat(this.value);
            display.textContent = parseFloat(this.value).toFixed(2);
        });
    });
}

function toggleSimulation() {
    running = !running;
}

function toggleVisualization() {
    vizMode = (vizMode === 'timeseries') ? 'flow' : 'timeseries';

    // Update button text
    let btn = document.getElementById('viz-toggle');
    if (btn) {
        btn.textContent = vizMode === 'timeseries' ? 'Show Flow Diagram' : 'Show Time Series';
    }
}

function resetSimulation() {
    N = N0;
    P = P0;
    Z = Z0;
    D = D0;
    time = 0;

    history.t = [time];
    history.N = [N];
    history.P = [P];
    history.Z = [Z];
    history.D = [D];

    running = false;

    // Reset parameters to defaults
    params = {...defaultParams};

    // Reset all sliders and their displays
    const controls = ['vm', 'kn', 'i0', 'rm', 'lambda', 'alpha', 'beta', 'epsilon', 'g', 'r', 'phi'];

    controls.forEach(param => {
        let slider = document.getElementById(param);
        let display = document.getElementById(param + '-val');

        if (slider && display) {
            slider.value = defaultParams[param];
            display.textContent = defaultParams[param].toFixed(2);
        }
    });
}

function setSpeed(multiplier) {
    speedMultiplier = multiplier;

    // Update button states
    let buttons = document.querySelectorAll('.speed-btn');
    buttons.forEach(btn => {
        btn.classList.remove('active');
    });

    // Find and activate the clicked button
    buttons.forEach(btn => {
        if (btn.textContent === `${multiplier}x`) {
            btn.classList.add('active');
        }
    });
}
