#define PI 3.1415926

// RYB/RGB
#define RYB 1

// smoothed color gradient or not
#define SMOOTH 0


// Color space conversion

float hue2rgb(float p, float q, float t) {
    if (t < 0.) t += 1.;
    if (t > 1.) t -= 1.;
    if (t < 1./6.) return mix(p, q, 6. * t);
    if (t < 1./2.) return q;
    if (t < 2./3.) return mix(p, q, (2./3. - t) * 6.);
    return p;
}
vec3 hsl2rgb(float h, float s, float l) {
    float q = l < 0.5 ? l * (1.0 + s) : l + s - l * s;
    float p = 2. * l - q;
    return vec3(
        hue2rgb(p, q, h + 1./3.),
        hue2rgb(p, q, h),
        hue2rgb(p, q, h - 1./3.)
    );
}

vec3 ryb2rgb(vec3 ryb) {
    const vec3 ryb000 = vec3(1, 1, 1);
    const vec3 ryb001 = vec3(0.163, 0.373, 0.6);
    const vec3 ryb010 = vec3(1, 1, 0);
    const vec3 ryb100 = vec3(1, 0, 0);
    const vec3 ryb011 = vec3(0, 0.66, 0.2);
    const vec3 ryb101 = vec3(0.5, 0, 0.5);
    const vec3 ryb110 = vec3(1, 0.5, 0);
    const vec3 ryb111 = vec3(0, 0, 0);
    return mix(mix(
        mix(ryb000, ryb001, ryb.z),
        mix(ryb010, ryb011, ryb.z),
        ryb.y), mix(
        mix(ryb100, ryb101, ryb.z),
        mix(ryb110, ryb111, ryb.z),
        ryb.y), ryb.x);
}


// Mask layer - please ignore my terrible tracing skill

float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) - k * h * (1.0 - h);
}
float smax(float a, float b, float k) {
    return -smin(-a, -b, k);
}

float mask(float x, float y) {
    float head = smax(x+0.25*y-0.65, -0.27*x+y-0.94-0.01*cos(3.0*x), 0.6); // head back
    float face = -x-0.08*y-0.6 + 0.03*cos(6.0*y-1.0);  // face
    face -= -0.2+0.2*tanh(18.0*(y+0.4)) + 0.1*exp(-3.0*y-2.2); // chin/neck
    face -= 0.04*exp(-100.0*(y-0.1)*y)*(0.8-x+10.0*y);  // nose
    face -= 0.03*sin(40.0*(y-0.03))*exp(-100.0*y*y);  // mouth
    float d = smax(head, face, 0.2);
    d = smax(d, 0.3*x-y-0.87 - 0.1*x*x, 0.05);  // bottom
    float back = x+0.7*y+0.05 - 0.1*exp(-3.0*y-2.2);  // bottom-right block
    back = smax(-back, -0.6*x+y+0.72-0.05*sin(8.0*x), 0.02);
    d = smax(d, -back, 0.02);
    return d;
}


// Noise

vec2 hash22(vec2 p) {
    // from David Hoskins's https://www.shadertoy.com/view/4djSRW
    vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

float SimplexNoise(vec2 xy) {
    const float K1 = 0.3660254038;  // (sqrt(3)-1)/2
    const float K2 = 0.2113248654;  // (-sqrt(3)+3)/6
    vec2 p = xy + (xy.x + xy.y)*K1;
    vec2 i = floor(p);
    vec2 f1 = xy - (i - (i.x + i.y)*K2);
    vec2 s = f1.x < f1.y ? vec2(0.0, 1.0) : vec2(1.0, 0.0);
    vec2 f2 = f1 - s + K2;
    vec2 f3 = f1 - 1.0 + 2.0*K2;
    vec2 n1 = 2.0 * hash22(i) - 1.0;
    vec2 n2 = 2.0 * hash22(i + s) - 1.0;
    vec2 n3 = 2.0 * hash22(i + 1.0) - 1.0;
    vec3 v = vec3(dot(f1, n1), dot(f2, n2), dot(f3, n3));
    vec3 w = max(-vec3(dot(f1, f1), dot(f2, f2), dot(f3, f3)) + 0.5, vec3(0.0));
    return dot((w*w*w*w) * v, vec3(32.0));
}


// Main

void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // Cartesian coordinate
    vec2 pos = 2.0 * fragCoord/iResolution.xy - 1.0;
    pos *= iResolution.xy / min(iResolution.x, iResolution.y);

    // polar coordinate
    vec2 p = pos - vec2(-0.0, 0.15);
    float r = length(p);
    float a = atan(p.x, p.y) / (2.0*PI);
    if (a < 0.0) a += 1.0;

    // color
    float noise = SimplexNoise(vec2(40.0*a, 8.0*sqrt(r)));
#if !SMOOTH
    a = round(12.0*a)/12.0;
    r = round(3.0*r)/3.0;
#endif
    vec3 col = hsl2rgb(a, 1.0, (1.0-0.1*noise)*(1.0-0.8*exp(-0.8*r)));
#if RYB
    col = ryb2rgb(col);
#endif

    // apply mask
    float m = mask(pos.x, pos.y);
    if (m > 0.0) col = vec3(1.0);

    fragColor = vec4(col,1.0);
}
