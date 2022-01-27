// Visualize the SDF (approximate) and AABB (exact) of an oriented L^n-ellipsoid
// Equation: (|x|/rx)^n + (|y|/ry)^n + (|z|/rz)^n = 1

// SDF inspired by iq's ellipsoid SDF approximation
// Doesn't work well inside the shape and for n < 1.5

// Exact AABB derived using Lagrange multiplier
// I'm not sure when I will use this, I derived it as a math practice and for fun

// Forked from a 3D SDF visualizer - https://www.shadertoy.com/view/ssKGWR (v1),
// https://github.com/harry7557558/Shadertoy/blob/master/sdf_visualizer_2.glsl (v2)
//  - Visualize the SDF using isolines (volume rendering);
//  - Red highlight high directional gradient, which causes issues in raymarching;
//  - Magenta highlight high surface gradient, green highlight low surface gradient;

// There are people who reported that my previous shaders look very dim
// while they look appropriate on my screen.
// If the shader appears too dark/bright, adjust the following value:
#define GAMMA 0.8


// START OF ELLIPSOID STUFF ================================

// Approximation of the SDF of an L^n-norm ellipsoid
// inspired by iq's ellipsoid SDF approximation
float sdLnNormEllipsoid(vec3 p, vec3 r, float n) {
    p = abs(p);
    float k1 = pow(dot(pow(p/r,vec3(n)),vec3(1)),1.0/n);
    float k2 = pow(dot(pow(p/(r*r),vec3(n)),vec3(1)),1.0/n);
    return k1*(k1-1.0)/k2;
}

// Axes-aligned bounding box of a rotated L^n-norm ellipsoid
// Derived using Lagrange multiplier, return the radius of the box
vec3 aabbLnNormEllipsoid(vec3 r, float n, mat3 rotation) {
    float res[3];
    for (int i=0; i<3; i++) {  // iterate through each dimension
        vec3 d = abs(rotation[i]);
        float n_lambda = pow(dot(pow(d*r,vec3(n/(n-1.))),vec3(1)),(n-1.)/n);
        res[i] = dot(d, pow(d*pow(r,vec3(n))/n_lambda,vec3(1./(n-1.))));
    }
    return vec3(res[0], res[1], res[2]);
}

// END OF ELLIPSOID STUFF ================================


// Get the radius, power of norm, and inverse rotation matrix of the shape, for demonstration purpose
void getParameters(out vec3 r, out float n, out mat3 rotation) {
    // radius among axes
    r = vec3(1.0+0.5*cos(1.6*iTime), 1.0+0.5*cos(1.5*iTime), 1.0+0.5*cos(1.4*iTime));
    // power of norm
    n = mix(1.5, 6.0, 0.5+0.5*cos(1.0*iTime));
    // rotation matrix (axis-angle)
    vec3 u = normalize(vec3(sin(1.3*iTime), sin(1.2*iTime), sin(1.1*iTime)));
    float a = sin(0.9*iTime), s = sin(a), c = cos(a), oc = 1.0-c;
    rotation =  mat3(
        oc*u.x*u.x+c, oc*u.x*u.y-u.z*s, oc*u.z*u.x+u.y*s,
        oc*u.x*u.y+u.z*s, oc*u.y*u.y+c, oc*u.y*u.z-u.x*s,
        oc*u.z*u.x-u.y*s, oc*u.y*u.z+u.x*s, oc*u.z*u.z+c
    );
}

// SDF and numerical gradient
float sdf(vec3 p) {
    vec3 r; float n; mat3 rotation;
    getParameters(r, n, rotation);
    return sdLnNormEllipsoid(rotation*p, r, n);
}
vec3 sdfGrad(vec3 p, float e) {
	float a = sdf(p+vec3(e,e,e));
	float b = sdf(p+vec3(e,-e,-e));
	float c = sdf(p+vec3(-e,e,-e));
	float d = sdf(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}

// 3D SDF visualizer template
// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

// constants
#define PI 3.1415926
#define ZERO min(iTime,0.)

// raymarching parameters
#define STEP 0.1
#define MIN_STEP 0.002

// rendering parameters
#define FIELD_EMISSION (insideAABB(ro+rd*t)?0.4:0.05)
#define ISOSURFACE_FREQUENCY 6.0
#define DISCONTINUITY_OPACITY 0.2
#define SURFACE_GRADIENT 10.0

// projection parameters
#define PERSPECTIVE 10.0  /* larger: less perspective effect */
#define SCALE 6.0  /* image appears smaller when this is set larger */

// light direction as a global variable
vec3 light;

// colormaps - https://www.shadertoy.com/view/NsSSRK
vec3 colorSdf(float t) {
    float r = .385+.619*t+.238*cos(4.903*t-2.61);
    float g = -5.491+.959*t+6.089*cos(.968*t-.329);
    float b = 1.107-.734*t+.172*cos(6.07*t-2.741);
    return clamp(vec3(r, g, b), 0.0, 1.0);
}
vec3 colorNormal(float t) {
    float r = .529-.054*t+.55*cos(5.498*t+2.779);
    float g = .21+.512*t+.622*cos(4.817*t-1.552);
    float b = .602-.212*t+.569*cos(5.266*t+2.861);
    return clamp(vec3(r, g, b), 0.0, 1.0);
}

// test if inside or outside AABB
bool insideAABB(vec3 p) {
    vec3 r; float n_; mat3 rotation;
    getParameters(r, n_, rotation);
    vec3 box_radius = aabbLnNormEllipsoid(r, n_, rotation);
    return abs(p.x)<box_radius.x && abs(p.y)<box_radius.y && abs(p.z)<box_radius.z;
}

// raymarhing, return RGB
vec3 render(in vec3 ro, in vec3 rd, float t0, float t1) {
    float t = t0;
    vec3 totcol = vec3(0.0);
    float totabs = 1.0;
    float v_old = sdf(ro+rd*t), v;
    float dt = min(STEP, abs(v_old));
    for (float i=ZERO; i<200.; i++) {
        t += dt;
        if (t > t1) return totcol;
        v = sdf(ro+rd*t);
        if (v*v_old<0.) break;
        vec3 col = colorSdf(0.5+0.5*sin(ISOSURFACE_FREQUENCY*PI*0.5*(v_old+v)));
        float grad = abs(v-v_old)/dt;
        float grad_abs = (1.0-grad)/dt;
        col = mix(vec3(1,0,0), col, clamp(exp(grad_abs),0.0,1.0));
        float absorb = FIELD_EMISSION+DISCONTINUITY_OPACITY*max(-grad_abs,0.0);
        totabs *= exp(-absorb*dt);
        totcol += col*absorb*totabs*dt;
        dt = clamp(abs(v_old=v), MIN_STEP, STEP);
    }
    vec3 grad = sdfGrad(ro+rd*t, 1e-3);
    float grad_col = SURFACE_GRADIENT*(0.5*length(grad)-0.5);
    vec3 col = colorNormal(1.0-1.0/(1.0+exp(2.0*grad_col)));  // 0.5+0.5*tanh(grad_col)
    col = 0.2+0.05*grad.y+col*max(dot(normalize(grad), light),0.0);
    return totcol + col * totabs;
}

// ray intersection with a box
bool boxIntersection(vec3 ro, vec3 rd, out float tn, out float tf) {
    vec3 inv_rd = 1.0 / rd;
    vec3 n = inv_rd*(ro);
    vec3 k = abs(inv_rd)*vec3(2.5);
    vec3 t1 = -n - k, t2 = -n + k;
    tn = max(max(t1.x, t1.y), t1.z);
    tf = min(min(t2.x, t2.y), t2.z);
    if (tn > tf) return false;
    return true;
}

// main
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // set camera
    float rx = iMouse.z!=0. ? 3.14*(iMouse.y/iResolution.y)-1.57 : 0.3;
    float rz = iMouse.z!=0. ? -iMouse.x/iResolution.x*4.0*3.14 : -1.0;
    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);
    vec3 ro = SCALE*PERSPECTIVE*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, PERSPECTIVE*length(iResolution.xy));
    rd = normalize(rd);
    // calculate pixel color
    light = normalize(w+0.5*u+0.1*v);
    float t0, t1;
    if (!boxIntersection(ro, rd, t0, t1)) {
        fragColor = vec4(vec3(0.0), 1.0);
        return;
    }
    vec3 col = render(ro, rd, t0, t1);
    col = pow(col, vec3(GAMMA));
    fragColor = vec4(col, 1.0);
}
