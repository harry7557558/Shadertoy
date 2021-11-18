#iChannel0 "self"


uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }

#define PI 3.1415926
#define ZERO min(iTime,0.)
#define EPSILON 1e-4

const int MAT_BACKGROUND = -1;
const int MAT_NONE = 0;
const int MAT_PLANE = 1;
const int MAT_GLASS = 2;
const int MAT_CONTENT = 3;


// Model object

#ifndef __CPLUSPLUS
float sdLnNormEllipsoid(vec3 p, vec3 r, vec3 n) {
    float d = pow(dot(pow(abs(p)/r,n),vec3(1)),1.0/max(max(n.x,n.y),n.z))-1.0;
    float m = min(min(r.x,r.y),r.z);
    return d * m;
}
vec4 mapContent(vec3 p, bool col_required) {  // demo
    float d = sdLnNormEllipsoid(p-vec3(0,0,-0.0), vec3(1.2,1.2,2.4), vec3(1.5));
    return vec4(1.0,0.5,0.0, d);
}
vec3 gradContent(in vec3 p) {
    const float e = 0.001;
	float a = mapContent(p+vec3(e,e,e), false).w;
	float b = mapContent(p+vec3(e,-e,-e), false).w;
	float c = mapContent(p+vec3(-e,e,-e), false).w;
	float d = mapContent(p+vec3(-e,-e,e), false).w;
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}
#endif


// Model glass

#ifndef __CPLUSPLUS
float sqr(float x) { return x*x; }
float mapGlass(vec3 p) {
    vec3 r = 0.85*vec3(vec2(2.2,1.8)/(1.0+sqr(0.1*(p.z-2.0))),3.6);
    p.z -= r.z+0.01;
    return sdLnNormEllipsoid(p, vec3(r), vec3(18.0,18.0,24.0));
}
vec3 gradGlass(in vec3 p) {
    const float e = 0.001;
	float a = mapGlass(p+vec3(e,e,e));
	float b = mapGlass(p+vec3(e,-e,-e));
	float c = mapGlass(p+vec3(-e,e,-e));
	float d = mapGlass(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}
// to be compatible with mesh generator
const vec3 BOX_RADIUS = vec3(8.0, 8.0, 8.0);
const float STEP = 1.0;
const float MIN_STEP = 0.01;
vec4 map(vec3 p, bool col_required) { return vec4(0,0,0, mapGlass(p)); }
#endif


// Ray sampling and scattering

vec3 sampleCosWeighted(vec3 n) {
    vec3 u = normalize(cross(n, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, n);
    float rn = rand01();
    float an = 2.0*PI*rand01();
    vec2 rh = sqrt(rn) * vec2(cos(an), sin(an));
    float rz = sqrt(1. - rn);
    return rh.x * u + rh.y * v + rz * n;
}

vec3 sampleFresnelDielectric(vec3 rd, vec3 n, float n1, float n2) {
    float eta = n1 / n2;
    float ci = -dot(n, rd);
    if (ci < 0.0) ci = -ci, n = -n;
    float ct = 1.0 - eta * eta * (1.0 - ci * ci);
    if (ct < 0.0) return rd + 2.0*ci*n;
    ct = sqrt(ct);
    float Rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct);
    float Rp = (n1 * ct - n2 * ci) / (n1 * ct + n2 * ci);
    float R = 0.5 * (Rs * Rs + Rp * Rp);
    return rand01() > R ?
        rd * eta + n * (eta * ci - ct)  // refraction
        : rd + 2.0*ci*n;  // reflection
}

vec3 sampleCookTorrance(
    vec3 wi, vec3 n,
    float alpha,  // roughness
    float f0,  // ratio of reflection along the normal
    float lambertian,  // ratio of lambertian coefficient
    vec3 lambert_col,  // lambertian color
    vec3 microfacet_col,  // microfacet color
    inout vec3 m_col
    ) {

    // importance sample Lambertian
    if (rand01() < lambertian) {
        vec3 wo = sampleCosWeighted(n);
        m_col *= lambert_col;
        return wo;
    }

    // transform
    vec3 u = normalize(cross(n, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, n);
    wi = vec3(dot(wi, u), dot(wi, v), dot(wi, n));
    vec3 wo, m;  // out and half vector

    // GGX importance sampling
    float su = 2.0*PI*rand01();
    float sv = rand01();
    sv = atan(alpha*sqrt(sv/(1.0-sv)));
    m = vec3(sin(sv)*vec2(cos(su),sin(su)), cos(sv));
    wo = -(wi-2.0*dot(wi,m)*m);
    float D = wo.z<0. ? 0. : 4.0*dot(wi, m);

    // Geometry
    float tan2_theta_i = (1.0-wi.z*wi.z)/(wi.z*wi.z);
    float tan2_theta_o = (1.0-wo.z*wo.z)/(wo.z*wo.z);
    float lambda_i = 0.5*(sqrt(1.0+alpha*alpha*tan2_theta_i)-1.0);
    float lambda_o = 0.5*(sqrt(1.0+alpha*alpha*tan2_theta_o)-1.0);
    float G = 1.0/(1.0+lambda_i+lambda_o);

    // Fresnel
    float F = f0 + (1.0-f0)*pow(1.0-dot(wi, m), 5.0);

    // Put all together
    float Fr = D*G*F / (4.0*wi.z*wo.z+EPSILON);
    float Fr_cos = Fr * wo.z;  // wo is the direction of light in path tracing
    m_col *= Fr_cos * microfacet_col;
    return wo.x * u + wo.y * v + wo.z * n;
}


// Intersection functions

#ifndef __CPLUSPLUS

bool intersectGlass(vec3 ro, vec3 rd, inout float t, in float t1, out vec3 n) {
    float v_old = mapGlass(ro+rd*t), v;
    float dt = min(STEP, abs(v_old));
    for (int i=int(ZERO); i<128; i++) {
        t += dt;
        if (t > t1) return false;
        v = mapGlass(ro+rd*t);
        if (v*v_old<0.) break;
        dt = clamp(abs(v_old=v), MIN_STEP, STEP);
    }
    if (v*v_old<0.) {
        float t0 = t-dt, t1 = t;
        float v0 = v_old, v1 = v;
        for (int s = int(ZERO); s < 8; s++) {
            float t = 0.5*(t0+t1);
            float v = mapGlass(ro+rd*t);
            if (v*v0 < 0.0) t1 = t, v1 = v;
            else t0 = t, v0 = v;
        }
        t = t0;
        n = gradGlass(ro+rd*t);
        return true;
    }
    return false;
}

bool intersectContent(vec3 ro, vec3 rd, inout float t, in float t1, out vec3 n, out vec3 col) {
    const float STEP = 1.0, MIN_STEP = 0.005;
    float v_old = mapContent(ro+rd*t, false).w, v;
    float dt = min(STEP, abs(v_old));
    for (int i=int(ZERO); i<128; i++) {
        t += dt;
        if (t > t1) return false;
        v = mapContent(ro+rd*t, false).w;
        if (v*v_old<0.) break;
        dt = clamp(abs(v_old=v), MIN_STEP, STEP);
    }
    if (v*v_old<0.) {
        float t0 = t-dt, t1 = t;
        float v0 = v_old, v1 = v;
        for (int s = int(ZERO); s < 8; s++) {
            float t = 0.5*(t0+t1);
            float v = mapContent(ro+rd*t, false).w;
            if (v*v0 < 0.0) t1 = t, v1 = v;
            else t0 = t, v0 = v;
        }
        t = t0;
        n = gradContent(ro+rd*t);
        col = mapContent(ro+rd*t, true).xyz;
        return true;
    }
    return false;
}

#endif


// Rendering

vec3 light(vec3 rd) {
    const vec3 sunpos = normalize(vec3(0.2, -0.5, 0.2));
    //vec3 col = sin(20.0*rd.x)*sin(20.0*rd.y)*sin(20.0*rd.z)>0.0 ? vec3(1.0) : vec3(0.8);
    vec3 col = vec3(0.9+0.2*sin(2.0*rd.x+rd.y));
    vec3 amb = vec3(1.0) + vec3(2.0) * pow(max(dot(rd, sunpos), 0.), 4.);
    vec3 sun = (dot(rd,sunpos)>0.9 ? 1.0 : 0.0) * vec3(8.0);
    return col * 0.5*amb + sun;
}

vec3 mainRender(vec3 ro, vec3 rd) {

    vec3 m_col = vec3(1.0), t_col = vec3(0.0), col;
    bool inside_glass = false, inside_object = false;

    for (int iter = 0; iter < 64; iter++) {
        ro += EPSILON*rd;
        float t, min_t = 1e12;
        vec3 n, min_n;
        vec3 min_ro = ro, min_rd = rd;
        vec3 min_emi = vec3(0.0);
        int material = MAT_BACKGROUND;

        // plane
        t = -ro.z / rd.z;
        if (t > 0.0) {
            min_t = t, min_n = vec3(0, 0, 1);
            min_ro = ro + rd * t, min_rd = rd;
            col = vec3(0.5);
            material = MAT_PLANE;
        }

        // glass
        t = 0.0;
        if (intersectGlass(ro, rd, t, min_t, min_n)) {
            min_t = t;
            min_ro = ro + rd * t, min_rd = rd;
            min_n = normalize(min_n);
            col = vec3(1.0);
            material = MAT_GLASS;
        }

        // content
        t = 0.0;
        if (inside_glass) {
            if (intersectContent(ro-vec3(0,0,2.8), rd, t, min_t, min_n, col)) {
                min_t = t;
                min_ro = ro + rd * t, min_rd = rd;
                min_n = normalize(min_n);
                material = MAT_CONTENT;
            }
        }

        // update ray
        if (material == MAT_BACKGROUND) {
            if (iter == 0) return vec3(0.0);
            col = light(rd);
            return m_col * col + t_col;
        }
        if (inside_object);
        else if (inside_glass) m_col *= exp(-0.1*vec3(0.0,0.2,0.4)*min_t);
        min_n = dot(rd, min_n) < 0. ? min_n : -min_n;  // ray hits into the surface
        ro = min_ro, rd = min_rd;
        if (material == MAT_PLANE) {
            rd = sampleCookTorrance(-rd, min_n, 0.01, 0.1, 0.01, col, col, m_col);
        }
        else if (material == MAT_GLASS) {
            vec2 eta = inside_glass ? vec2(1.5, 1.0) : vec2(1.0, 1.5);
            rd = sampleFresnelDielectric(rd, min_n, eta.x, eta.y);
            if (dot(rd, min_n) < 0.0) inside_glass = !inside_glass;
        }
        else if (material == MAT_CONTENT) {
            rd = sampleCookTorrance(-rd, min_n, 0.5, 0.8, 0.4, vec3(1.0), vec3(1.0), m_col);
            m_col *= 1.6*pow(col, vec3(1.0));
            if (dot(rd, min_n) < 0.0) inside_object = !inside_object;
        }
        if (m_col == vec3(0.0)) return t_col;
        if (inside_object) return 1e12f*vec3(1,-1,-1);  // red warning
    }
    return m_col + t_col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // random number seed
    seed = uint(fragCoord.x)*uint(fragCoord.y)*uint(iFrame+1);
    seed = randu() + 161u*uint(fragCoord.y);
    seed = randu() + 239u*uint(fragCoord.x);
    seed = randu() + 197u*uint(iFrame+1);

    // camera
    float rx = iMouse.y==0.0 ? 0.33 : 1.8*(iMouse.y/iResolution.y)-0.3;
    float rz = iMouse.x==0.0 ? -7.6 : -iMouse.x/iResolution.x*4.0*3.14;
    rx = 0.33, rz = -7.6;
    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);
    vec3 ro = 20.0*w + vec3(0, 0, 3.0);
    vec2 uv = 5.0*(2.0*(fragCoord.xy+vec2(rand01(),rand01())-0.5)/iResolution.xy - vec2(1.0));
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 8.0*length(iResolution.xy));
    rd = normalize(rd);

    // calculate pixel color
    vec3 col = mainRender(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(int(fragCoord.x), int(fragCoord.y)), 0);
    if (iMouse.z>0.) rgbn.w = 0.0;
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}
