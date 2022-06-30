#define _uv float u, float v
#define _xyz float x, float y, float z

#define PI 3.1415926
const float p = PI; // pi
const float a_o = 0.16*p; // half of opening angle
const float b = 0.6; // r=e^bt
float s_min(float a, float b, float k) { return -1./k*log(exp(-k*a)+exp(-k*b)); } // smoothed minimum

// Cross section
float C_m(_uv) { return 1.-(1.-0.01*exp(sin(12.*p*(u+2.*v))))*exp(-5.*v*5.*v); } // mid rod
float C_s(_uv) {
    float _x = u-exp(-16.*v);
    float _y = v*(1.-0.2*exp(-4.*sqrt(u*u+.1*.1)))-0.5+0.5*exp(-v)*sin(4.*u)+.2*cos(2.*u)*exp(-v);
    return (sqrt(_x*_x+_y*_y)-0.55)*tanh(5.*sqrt(2.*u*u+(v-1.2)*(v-1.2)))+.01*sin(40.*u)*sin(40.*v)*exp(-(u*u+v*v));
}
float C_0(_uv) { return abs(C_s(u,v))*C_m(u,v); } // single layer
float n_1(_uv) { return log(sqrt(u*u+v*v))/b+2.; } // index of layer
float a_1(_uv) { return atan(v,u)/a_o; } // opening angle, 0-1
float d_1(_uv, float s_d) {
    float n = n_1(u,v);
    return 0.5*sqrt(u*u+v*v)*C_0(n>0.?n-s_d:fract(n)-s_d,a_1(u,v));
}
float C(_uv) { return min(d_1(u,v,0.5),d_1(u,v,1.5)); } // cross section

// Spiral
float l_p(float x, float y) { return exp(b*atan(y,x)/(2.*p)); } // a multiplying factor
float U(_xyz) { return exp(log(-z)+b*atan(y,x)/(2.*p)); } // xyz to cross section u
float V(_xyz) { return sqrt(x*x+y*y)*l_p(x,y); } // xyz to cross section v
float S_s(_xyz) { return C(U(x,y,z),V(x,y,z))/l_p(x,y); } // body
float S_o(_xyz) { return sqrt(pow(C(exp(log(-z)-b/2.),-x*exp(-b/2.))*exp(b/2.),2.)+y*y); } // opening
float S_t(_xyz) { return d_1(-z,sqrt(x*x+y*y),0.5); } // tip
float S_a(_xyz) { return -z>0.?min(S_s(x,y,z),S_o(x,y,z)):S_t(x,y,z); } // body+tip
float S_0(_xyz) { return S_a(x,y,z)-0.01-0.01*pow(x*x+y*y+z*z,0.4)
    -0.02*sqrt(x*x+y*y)*exp(cos(8.*atan(y,x)))
    -0.007*(0.5-0.5*tanh(10.*(z+1.+8.*sqrt(3.*x*x+y*y)))); } // subtract thickness
float S_r(_xyz) { return -s_min(-S_0(x,y,z),z+1.7,10.); } // clip bottom
float r_a(_xyz) { return -0.1*sin(3.*z)*tanh(2.*(x*x+y*y-z-1.5)); }
float S(_xyz) { return S_r(x-r_a(x,y,z)*y,y+r_a(x,y,z)*x,z-0.8); }


#define ZERO min(iTime,0.)

// rotation matrices
mat2 rot2(float a) { return mat2(cos(a), sin(a), -sin(a), cos(a)); }
mat3 rotx(float a) { return mat3(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a)); }
mat3 rotz(float a) { return mat3(cos(a), sin(a), 0, -sin(a), cos(a), 0, 0, 0, 1); }
mat3 roty(float a) { return mat3(cos(a), 0, -sin(a), 0, 1, 0, sin(a), 0, cos(a)); }

// calculate signed distance only
float sdf(vec3 p) {
    p.z += 0.5;
    vec3 q = rotz(0.125*PI)*rotx(0.395*PI)*(0.7*p-vec3(0,0,0.275));  // orientation/position
    float d = p.z;
    float bound = length(vec3(vec2(1.2,1.4)*exp(q.z*q.z),1.)*q)/exp(q.z*q.z)-1.0;
    bound = max(bound, length(vec3(1.2,1.4,1)*(q+vec3(0,0.1,0)))-1.);
    float boundw = 0.2;  // bounding box
    if (bound > 0.0) d = min(d, bound+boundw);
    else {
        float v = S(q.x,q.y,q.z);
        float k = 1.0-0.9/length(vec3(4.*q.xy,1.0*abs(q.z+0.7)+1.));  // reduce gradient
        k = 0.7*mix(k, 1.0, clamp(10.*max(-q.x,q.z-.7*q.x+0.5), 0., 1.));
        v = k*v/0.7;  // dividing by 0.7 is due to scaling
        v = mix(v, bound+boundw, smoothstep(0.,1.,(bound+boundw)/boundw));  // continuous transition
        d = min(d, v);
    }
    return d;
}

vec3 sdfGrad(in vec3 p, in float e) {
	float a = sdf(p+vec3(e,e,e));
	float b = sdf(p+vec3(e,-e,-e));
	float c = sdf(p+vec3(-e,e,-e));
	float d = sdf(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}

vec3 sdfColor(vec3 p, vec3 g) {
    float t = 0.5-0.5*cos(2.0*log(0.6*length(g)));
    t += 0.05*sin(40.*p.x)*sin(40.*p.y)*sin(20.*p.z);
    vec3 col = mix(vec3(0.9,0.9,0.85), vec3(0.75,0.55,0.3), t);
    if (p.z+0.5<1e-2) col = vec3(0.0);
    return col;
}

#define COLOR 1


// raymarching parameters
const vec3 BoxRadius = vec3(1.5,1.5,1.5);
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
    t -= dt * v/(v-v_old);
    vec3 grad = sdfGrad(ro+rd*t, 1e-3);
    vec3 col = colorNormal(0.5+0.5*tanh(SURFACE_GRADIENT*(0.5*length(grad)-0.5)));
#if COLOR
    col = sdfColor(ro+rd*t, grad);
#endif
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
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 3.0*length(iResolution.xy));
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
