#iChannel0 "self"

#iChannel1 "../cubemaps/shadertoy_uffizi_gallery/{}.jpg"
#iChannel1::Type "CubeMap"


#define PI 3.1415926
#define ZERO min(iTime, 0.)


uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }



vec3 light(vec3 rd) {
    //return vec3(0.0);
    vec3 col = texture(iChannel1, rd).xyz;
    vec3 bri = vec3(1.0) + vec3(2.0) * pow(max(dot(rd, normalize(vec3(-0.2, -0.5, 0.5))), 0.), 4.);
    return col * bri;
}


// sphere intersection function
bool intersectSphere(vec3 o, float r, vec3 ro, vec3 rd,
        inout float t, inout vec3 n) {
    ro -= o;
    float b = -dot(ro, rd), c = dot(ro, ro) - r * r;
    float delta = b * b - c;
    if (delta < 0.0) return false;
    delta = sqrt(delta);
    float t1 = b - delta, t2 = b + delta;
    if (t1 > t2) t = t1, t1 = t2, t2 = t;
    if (t1 > t || t2 < 0.) return false;
    t = t1 > 0. ? t1 : t2;
    n = normalize(ro + rd * t);
    return true;
}


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

vec3 sampleUniformSphere() {
    float u = 2.0*PI*rand01();
    float v = 2.0*rand01()-1.0;
    return vec3(vec2(cos(u), sin(u))*sqrt(1.0-v*v), v);
}

vec3 sampleHenyeyGreenstein(vec3 wi, float g) {
    if (g == 0.0) return sampleUniformSphere();
    if (g >= 1.0) return wi;
    if (g <= -1.0) return -wi;
    float us = rand01();
    float vs = 2.0*PI*rand01();
    float z = (1.0+g*g-pow((1.0-g*g)/(2.0*g*(us+(1.0-g)/(2.0*g))),2.0))/(2.0*g);
    vec2 xy = vec2(cos(vs), sin(vs)) * sqrt(1.0-z*z);
    vec3 u = normalize(cross(wi, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, wi);
    vec3 wo = normalize(xy.x*u + xy.y*v + z*wi);
    return wo;
}



const int mat_none = -1;
const int mat_background = 0;
const int mat_lambertian = 1;
const int mat_refractive = 2;


// https://www.shadertoy.com/view/XljGzV
vec3 hsl2rgb(float h, float s, float l) {
    vec3 rgb = clamp(abs(mod(h*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0);
    return l + s * (rgb-0.5)*(1.0-abs(2.0*l-1.0));
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

void calcAbsorb(in vec3 p, out vec3 emi, out vec3 tabs, out vec3 sabs, out float k, out float g) {
    p -= vec3(0, 0, 1);  // center of sphere
    vec3 q = p; q.xz = mat2(0.9,0.2,-0.2,0.9)*q.xz;
    emi = (abs(q.z+0.01)<0.02 && max(length(q.xy)-0.6, sin(7.0*PI*length(q.xy)))<0.) || length(abs(q)-vec3(0,0,0.25))<0.1 ?
         8.0*ryb2rgb(hsl2rgb(atan(q.y,q.x)/(2.0*PI)+0.5, 1.0, 0.5)) : vec3(0.0);
    tabs = vec3(0.0);
    sabs = vec3(0.9, 0.95, 0.99);
    k = clamp(-100.0*p.z+0.5, 0., 20.);
    g = 0.5;
}

void calcScatter(in vec3 ro, inout vec3 rd,
        inout float mt, inout vec3 m_col, out vec3 m_emi, out int material, inout vec3 min_n) {
    float p = 1.0, hit_p = rand01();
    float dt = 0.01;
    m_emi = vec3(0.0);
    for (float t = 1e-4; t < mt-dt; t += dt) {
        vec3 emi, tabs, sabs; float k, g;
        calcAbsorb(ro + rd * (t+0.5*dt), emi, tabs, sabs, k, g);
        float dp = exp(-k*dt);
        if (p * dp < hit_p) {
            dt *= log(p/hit_p)/k;
            mt = t + dt;
            rd = sampleHenyeyGreenstein(rd, g);
            m_col *= sabs * exp(-tabs*dt);
            material = mat_none;
            min_n = vec3(0.0);
            return;
        }
        p *= dp;
        m_col *= exp(-tabs*dt);
        m_emi += m_col * emi * dt;
    }
    material = mat_refractive;
}


vec3 mainRender(vec3 ro, vec3 rd) {

    vec3 m_col = vec3(1.0), t_col = vec3(0.0), col;
    bool is_inside = false;

    for (int iter = int(ZERO); iter < 128; iter++) {
        ro += 1e-4f*rd;
        if (is_inside != (length(ro-vec3(0,0,1))<1.0)) return vec3(1, 0, 0);

        vec3 n, min_n;
        float t, min_t = 1e12;
        vec3 min_ro = ro, min_rd = rd;
        vec3 min_emi = vec3(0.0);
        int material = mat_background;

        // plane
        t = -ro.z / rd.z;
        if (t > 0.0) {
            min_t = t, min_n = vec3(0, 0, 1);
            min_ro = ro + rd * t, min_rd = rd;
            col = vec3(0.9, 0.95, 0.98);
            material = mat_lambertian;
        }

        // object
        t = min_t;
        if (intersectSphere(vec3(0.0, 0.0, 1.0), 1.0, ro, rd, t, n)) {
            min_t = t, min_n = n;
            if (is_inside) {
                col = vec3(1.0);
                min_rd = rd;
                calcScatter(ro, min_rd, min_t, col, min_emi, material, min_n);
                min_ro = ro + rd * min_t;
            }
            else {
                min_ro = ro + rd * t, min_rd = rd;
                col = vec3(1.0);
                material = mat_refractive;
            }
        }

        // update ray
        if (material == mat_background) {
            col = light(rd);
            return m_col * col + t_col;
        }
        ro = min_ro, rd = min_rd;
        min_n = dot(rd, min_n) < 0. ? min_n : -min_n;  // ray hits into the surface
        if (material == mat_lambertian) {  // diffuse
            rd = sampleCosWeighted(min_n);
        }
        else if (material == mat_refractive) {  // steel ball
            vec2 eta = is_inside ? vec2(1.5, 1.0) : vec2(1.0, 1.5);
            rd = sampleFresnelDielectric(rd, min_n, eta.x, eta.y);
        }
        m_col = m_col * col;
        t_col += min_emi;
        if (dot(rd, min_n) < 0.0) {
            is_inside = !is_inside;
        }
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
    float rx = 2.0*(iMouse.y/iResolution.y)-0.5;
    //float rx = 3.14*(iMouse.y/iResolution.y)-1.57;
    float rz = -iMouse.x/iResolution.x*4.0*3.14;
    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);
    vec3 ro = 10.0*w + vec3(0, 0, 0.7);
    vec2 uv = 2.0*(fragCoord.xy+vec2(rand01(),rand01())-0.5)/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 3.0*length(iResolution.xy));
    rd = normalize(rd);

    // calculate pixel color
    vec3 col = mainRender(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(fragCoord), 0);
    if (iMouse.z>0.) rgbn.w = 0.0;
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}
