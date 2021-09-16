// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

#define PI 3.1415926


// \left[\cos \left(u\right)\cdot e^{b\cdot u}\cdot \left(1+w\cdot \cos \left(v\right)\right),\sin \left(u\right)\cdot e^{b\cdot u}\cdot \left(1+w\cdot \cos \left(v\right)\right),e^{b\cdot u}\cdot w\cdot \sin \left(v\right)\right]
float mapShell(in vec3 p0) {
    p0 -= vec3(0.7, 0, 0);
    vec3 p = p0;

    const float b = 0.17;

    float r = length(p.xy);
    float a = mix(0.0, 0.45, smoothstep(0.0, 1.0, 0.5*(r-0.6)));
    p.xy = mat2(cos(a),-sin(a),sin(a),cos(a))*p.xy;
    float t = atan(p.y, p.x);
    
    // shell opening
    float ro = exp(b*PI);
    float d = length(vec2(length(p.xz-vec2(-ro,0))-ro,p.y));
    float u = t, dx = r-ro, dy = p.z;

    // spiral
#if 0
    for (float n = -4.0; n <= 0.0; n += 1.0) {
        float rad = exp(b*(2.*PI*n+t));
        float d1 = abs(length(vec2(r-rad,p.z))-rad);
        if (d1 < d) d = d1, u = 2.*PI*n+t, dx = r-rad, dy = p.z;
    }
#else
    // r(n) = exp(b*(2.*PI*n+t)), (x-r)^2+y^2=r^2, solve for n
    float n = (log((r*r+p.z*p.z)/(2.*r))/b-t)/(2.0*PI);
    n = min(n, 0.0);
    float n0 = floor(n), n1 = ceil(n);
    float r0 = exp(b*(2.*PI*n0+t)), r1 = exp(b*(2.*PI*n1+t));
    float d0 = abs(length(vec2(r-r0,p.z))-r0);
    float d1 = abs(length(vec2(r-r1,p.z))-r1);
    if (d0 < d) d = d0, u = 2.*PI*n0+t, dx = r-r0, dy = p.z;
    if (d1 < d) d = d1, u = 2.*PI*n1+t, dx = r-r1, dy = p.z;
#endif

    // cells, possibly cause discontinuities
    const float f = 2.4;
    float s0 = t + 2.0*PI*(n0+0.5);
    float v = fract(n);
    float s = f*s0 + 1.0*pow(0.25-(v-0.5)*(v-0.5), 0.5)+0.5*v;
    s += pow(min(1.0/(40.0*length(vec2(v-0.5,p.z))+1.0), 0.5), 2.0);
    float sf = fract(s);
    sf = s0>-1.8 ? abs(s+3.25) : min(sf, 1.0-sf);
    float w = sf/f*exp(b*(s0+PI));
    if (length(p*vec3(1,1,1.5))<3.0) d = min(d, 0.5*w+0.012);

    // texture
    d += 0.002*r*sin(40.*u);
    //d += 0.002*r*sin(40.*atan(dy,dx));

    d = abs(d)-max(0.02*pow(r,0.4),0.02);
    //d = max(d, p0.x);
    d = max(d, p0.z);
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
const vec3 BoxRadius = vec3(3.0, 3.0, 2.0);
#define STEP 0.008
#define MAX_STEP 160.


// rendering parameters
#define FIELD_EMISSION 0.2
#define DISCONTINUITY_OPACITY 0.5
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
    float rx = iMouse.z!=0. ? 3.14*(iMouse.y/iResolution.y)-1.57 : 0.3;
    float rz = iMouse.z!=0. ? -iMouse.x/iResolution.x*4.0*3.14 : -0.6;

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
