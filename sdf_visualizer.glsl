// Debugging an SDF, visualize its discontinuity

// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

#define PI 3.1415926


// a debugging SDF, heavily inspired by iq's Snail
float mapShell(in vec3 p0) {
    p0 -= vec3(1.3, 0, 0);
    vec3 p = p0;

    const float b = 0.17;

    float r = length(p.xy);
    float a = mix(0.0, 0.8, smoothstep(0.0, 1.0, 0.5*(r-0.6)));
    p.xy = mat2(cos(a),-sin(a),sin(a),cos(a))*p.xy;
    float t = atan(p.y, p.x);
 
    float n = (log(r)/b-t)/(2.0*PI);
    n = min(n, 0.0);

    float n0 = floor(n), n1 = ceil(n);
    float x0 = exp(b*(t+2.0*PI*n0));
    float x1 = exp(b*(t+2.0*PI*n1));
    float r0 = 1.0*x0;
    float r1 = 1.0*x1;

    float h0 = p.z + 0.4*(x0-1.0);
    float h1 = p.z + 0.4*(x1-1.0);
    float d0 = length(vec2(x0-r,h0)) - r0;
    float d1 = length(vec2(x1-r,h1)) - r1;

    float d, dx, dy;
    if (d0 < 0.0) d = d0, dx = x0-r, dy = h0;
    else if (d1 < 0.0 && d1<-d0) d = -d0, dx = x0-r, dy = h0;
    else if (d1 < 0.0) d = d1, dx = x1-r, dy = h1;
    else if (d0 < d1) d = d0, dx = x0-r, dy = h0;
    else d = d1, dx = x1-r, dy = h1;

    d += 0.002*r*sin(40.*t);
    d += 0.002*r*sin(40.*atan(dy,dx));

    d = abs(d)-0.1*r;
    d = max(d, p0.x);
    return d;
}

// test SDF
float mapTest(vec3 p) {
    vec3 r = vec3(1.5,1.0,0.6);
    float k1 = length(p/r);
    float k2 = length(p/(r*r));
    return k1*(k1-1.0)/k2+0.1*sign(p.x+p.y+p.z);
}

float sdf(in vec3 p) {
    //return mapTest(p);
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
const vec3 BoxRadius = vec3(2.0, 2.0, 2.0);
#define STEP 0.01
#define MAX_STEP 120.


// rendering parameters
#define FIELD_EMISSION 0.3
#define DISCONTINUITY_OPACITY 0.1
#define SURFACE_GRADIENT 10.0

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
    float step_count = min(ceil((t1-t0)/STEP), MAX_STEP);
    float t = t0, dt = (t1-t0) / step_count;
    vec3 totcol = vec3(0.0);
    float totabs = 1.0;
    float v_old, v;
    for (t = t0; t < t1; t += dt) {
        v = sdf(ro+rd*t);
        vec3 col = colorSdf(0.5+0.5*sin(8.0*PI*v));
        //float grad = length(sdfGrad(ro+rd*t,dt));
        float grad = t==t0 ? 0.0 : abs(v-v_old)/dt;
        float grad_abs = (1.0-grad)/dt;
        col = mix(vec3(1,0,0), col, clamp(exp(grad_abs),0.0,1.0));
        float absorb = FIELD_EMISSION+DISCONTINUITY_OPACITY*max(-grad_abs,0.0);
        totabs *= exp(-absorb*dt);
        totcol += col*absorb*totabs*dt;
        if (v < 0.0) break;
        v_old = v;
    }
    if (v > 0.0) return totcol;
    for (int s = 0; s < 4; s += 1) {
        v_old = v;
        dt *= -0.5;
        for (int i = 0; i < 2; i++) {
            t += dt;
            v = sdf(ro+rd*t);
            if (v*v_old<0.0) break;
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
    vec3 k = abs(inv_rd)*BoxRadius;
    vec3 t1 = -n - k, t2 = -n + k;
    tn = max(max(t1.x, t1.y), t1.z);
    tf = min(min(t2.x, t2.y), t2.z);
    if (tn > tf) return false;
    return true;
}

// main
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // set camera
    float rx = iMouse.z>0. ? 3.14*(iMouse.y/iResolution.y)-1.57 : 0.3;
    float rz = iMouse.z>0. ? -iMouse.x/iResolution.x*4.0*3.14 : -0.6;

    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);

    vec3 ro = 10.0*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 2.0*length(iResolution.xy));
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
