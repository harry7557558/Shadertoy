// Power spiral, generalize Archimedean spiral to power k
// Might be more desireable than logarithmic spiral because it has a finite number of layers


#define PI 3.1415926

// k: the power of the spiral
// size: maximum distance of the spiral from the origin
// n: number of spiral layers
float distPowerSpiral(vec2 p, float k, float size, float n) {
    // r = a*θ^k
    float a = size / pow(2.0*PI*n, k);
    // polar coordinates
    float r = length(p), theta = atan(p.y, p.x);
    if (theta > 0.0) theta -= 2.0*PI;
    // index of layer
    float i = (pow(r/a, 1.0/k)-theta)/(2.0*PI);
    i = min(i, n);
    float i0 = floor(i), i1 = ceil(i);
    // two bounds along ray from origin
    float theta0 = max(2.0*PI*i0+theta,0.);
    float theta1 = max(2.0*PI*i1+theta,0.);
    float r0 = a*pow(theta0, k), r1 = a*pow(theta1, k);
    // distance
    float d = r-r0;
    if (theta1 < 2.0*PI*n) d = min(d, r1-r);
    // clipping at the end
    float re = a*pow(2.0*PI*n,k);
    vec2 pe = re*vec2(cos(2.0*PI*n), sin(2.0*PI*n));
    d = min(abs(d), length(p-pe));
    return d;
}


#define assertGrad(val,grad) return length(vec2(dFdx(val),dFdy(val))/(5.0/min(iResolution.x,iResolution.y))-(grad))

float distPowerSpiral1(vec2 p, float k, float size, float n) {
    // r = a*θ^k
    float a = size / pow(2.0*PI*n, k);
    // polar coordinates
    float r = length(p), theta = atan(p.y, p.x);
    if (theta > 0.0) theta -= 2.0*PI;
    // index of layer
    float i = (pow(r/a, 1.0/k)-theta)/(2.0*PI);
    i = min(i, n);
    float i0 = floor(i), i1 = ceil(i);
    // two bounds along ray from origin
    float theta0 = max(2.0*PI*i0+theta,-0.);
    float theta1 = max(2.0*PI*i1+theta,-0.);
    float r0 = a*pow(theta0, k), r1 = a*pow(theta1, k);
    // distance
    float d = r-r0;
    float adj = 0.5;
    if (theta0!=0.) d /= max(length(normalize(p)-adj*r0*k/theta0*vec2(-p.y,p.x)/dot(p,p)), 1.0);
    if (theta1 < 2.0*PI*n) {
        float d1 = r1-r;
        if (theta1!=0.) d1 /= max(length(normalize(p)-adj*r1*k/theta1*vec2(-p.y,p.x)/dot(p,p)), 1.0);
        d = min(d, d1);
    }
    // clipping at the end
    float re = a*pow(2.0*PI*n,k);
    vec2 pe = re*vec2(cos(2.0*PI*n), sin(2.0*PI*n));
    d = min(abs(d), length(p-pe));
    return d;
}


float sdf(in vec2 p) {
    float k = mix(0.5,4.0,0.5+0.5*sin(iTime));
    float size = 2.5;
    float n = mix(1.5, 2.5, 0.5+0.5*sin(0.6*iTime));
    return distPowerSpiral(p, k, size, n);
}




// 2D SDF Visualizer Template
// yellow-blue: positive/negative
// red-green: high/low gradient
// white: zero isoline

#define BOX_RADIUS 2.5
#define ISOLINE_WIDTH 0.1
#define SHOW_GRAD 1
#define HIGH_GRAD_HIGHLIGHT 1.0
#define LOW_GRAD_HIGHLIGHT 1.0

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    float e = (2.0*BOX_RADIUS)/min(iResolution.x,iResolution.y);
    vec2 p = (fragCoord-0.5*iResolution.xy)*e;

    float d = sdf(p);

#if SHOW_GRAD
    float h = 0.5*e;
    vec2 grad = vec2(
        sdf(p+vec2(h,0))-sdf(p-vec2(h,0)),
        sdf(p+vec2(0,h))-sdf(p-vec2(0,h))) / (2.0*h);
    float red = 1.0-exp(-HIGH_GRAD_HIGHLIGHT*max(length(grad)-1.0,0.));
    float green = 1.0-exp(-LOW_GRAD_HIGHLIGHT*max(1.0-length(grad),0.));
#else
    float red = 0.0;
    float green = 0.0;
#endif

    vec3 col = d>0.0 ? vec3(1.0,0.8,0.2) : vec3(0.3,0.6,1.0);  // sign
    col *= 0.7+0.3*cos(2.0*PI*d/ISOLINE_WIDTH);  // isolines
    col *= 1.0-0.8*exp(-0.5*abs(d)/ISOLINE_WIDTH);  // fade near zero
    col = mix(col, vec3(1.0), clamp(2.0-6.0*abs(d/ISOLINE_WIDTH),0.,1.));  // zero isoline
    col *= mix(vec3(1.0), vec3(2,0,0), red);  // high grad
    col *= mix(vec3(1.0), vec3(0,2,0), green);  // low grad
    fragColor = vec4(col, 1.0);
}
