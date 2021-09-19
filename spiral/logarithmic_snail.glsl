// *INCOMPLETE*

// Fork https://www.shadertoy.com/view/ssKGWR (3D SDF Visualizer)
// Should also mention https://www.shadertoy.com/view/sdVGWh (Nautilus Shell)

// Updates on SDF renderer:
// Uses dt=min(step,dist) instead of dt=step
// Allows a usually larger average marching step without missing thin surface
// (expected to be faster than the original one)

// Maximum near-surface gradient that maintains stability:
// at least (STEP+thickness)/STEP

// About the snail SDF:

// Inspired by logarithmic spiral
// 


#define PI 3.1415926
#define ZERO min(iTime,0.)


float mapShell(in vec3 p) {
    p.x-=0.5;

    float a = 0.8;
    float b = 0.15;
    vec2 param_p = vec2(0.0, 1.6);
    vec2 param_d = vec2(0.2, -1.5) + vec2(1.,0.);
    a=1.0, b=0.17, param_p = vec2(0.0), param_d = vec2(1.,0.);

    // convert to cylindrical coordinate
    float r = length(p.xy);
    float t = atan(p.y, p.x);
    vec2 q = vec2(r, p.z);

    // shell opening
    float rad = a*exp(b*PI);
    vec2 ro = param_p+param_d*rad*vec2(cos(PI),1.);
    float d = length(vec2(length(p.xz-ro)-rad,p.y));
    float u = t, dx = r-ro.x, dy = p.z-ro.y;  // can be used to apply texture
    //return d-0.2;

    // spiral
#if 0
    for (float n = -4.0; n <= 0.0; n += 1.0) {
        float rad = a*exp(b*(2.*PI*n+t));
        vec2 c = param_p + rad*param_d;
        float d1 = abs(length(q-c)-rad);
        if (d1 < d) d = d1, u = 2.*PI*n+t, dx = q.x-c.x, dy = q.y-c.y;
    }
#else
    // r(n) = a*exp(b*(2.*PI*n+t)), (q-(param_p+param_d*r(n)))^2=r(n)^2, solve for n
    float n;
    {
        float a = dot(param_d, param_d)-1.;
        float b = dot(param_d, q-param_p);
        float c = dot(q-param_p, q-param_p);
        float r = (sqrt(b*b-a*c)-b)/a;
        if (!isnan(r)) return 0.99;
        n = (log(r/a)/b-t)/(2.0*PI);
    }
    n = min(n, 0.0);
    float n0 = floor(n), n1 = ceil(n);
    float r0 = a*exp(b*(2.*PI*n0+t)), r1 = a*exp(b*(2.*PI*n1+t));
    vec2 c0 = param_p+r0*param_d, c1 = param_p+r1*param_d;
    float d0 = abs(length(q-c0)-r0);
    float d1 = abs(length(q-c1)-r1);
    if (d0 < d) d = d0, u = 2.*PI*n0+t, dx = q.x-c0.x, dy = q.y-c0.y;
    if (d1 < d) d = d1, u = 2.*PI*n1+t, dx = q.x-c1.x, dy = q.y-c1.y;
#endif
    d = abs(d);

    d += r*0.001*cos(64.0*u) + r*0.001*sin(64.0*atan(dy,dx));

    d = d-0.01;
    d = max(d, p.x);
    //d = max(d, p.z);
    return d;
}

float sdf(in vec3 p) {
    return mapShell(p);
}

vec3 sdfGrad(in vec3 p, in float e) {
	float a = sdf(p+vec3(e,e,e));
	float b = sdf(p+vec3(e,-e,-e));
	float c = sdf(p+vec3(-e,e,-e));
	float d = sdf(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}



// raymarching parameters
#define BOX_RADIUS vec3(3.0, 3.0, 2.0)
#define STEP 0.1
#define MAX_STEP 200.

// rendering parameters
#define FIELD_EMISSION 0.2
#define ISOSURFACE_FREQUENCY 8.0
#define DISCONTINUITY_OPACITY 0.5
#define SURFACE_GRADIENT 10.0

// projection parameters
#define PERSPECTIVE 10.0  /* larger: less perspective effect */
#define SCALE 6.0  /* image appears smaller when this is set larger */


// light direction as global variable
vec3 light = normalize(vec3(0.5,0.5,1.0));

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


// modified from a volume rendering demo
// https://github.com/harry7557558/Graphics/blob/master/raytracing/webgl_volume/fs-source.glsl
vec3 render(in vec3 ro, in vec3 rd, float t0, float t1) {
    float t = t0;
    vec3 totcol = vec3(0.0);
    float totabs = 1.0;
    float v_old = sdf(ro+rd*t), v;
    float dt = min(STEP, abs(v_old));
    for (float i=ZERO; i<MAX_STEP; i++) {
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
        if ((dt=min(STEP,abs(v_old=v))) < 1e-3) break;
    }
    if (v*v_old<0.) {
        for (int s = int(ZERO); s < 4; s += 1) {
            v_old = v, dt *= -0.5;
            for (int i = int(ZERO); i < 2; i++) {
                t += dt, v = sdf(ro+rd*t);
                if (v*v_old<0.0) break;
            }
        }
    }
    vec3 grad = sdfGrad(ro+rd*t, 1e-3);
    vec3 col = colorNormal(0.5+0.5*tanh(SURFACE_GRADIENT*(0.5*length(grad)-0.5)));
    col = 0.2+0.05*grad.y+col*max(dot(normalize(grad), light),0.0);
    return totcol + col * totabs;
}


// ray intersection with a box
bool boxIntersection(vec3 ro, vec3 rd, out float tn, out float tf) {
    vec3 inv_rd = 1.0 / rd;
    vec3 n = inv_rd*(ro);
    vec3 k = abs(inv_rd)*BOX_RADIUS;
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
    float rz = iMouse.z!=0. ? -iMouse.x/iResolution.x*4.0*3.14 : -0.6;

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
    vec3 col = render(ro, rd, t0, t1);;
    fragColor = vec4(col, 1.0);
}
