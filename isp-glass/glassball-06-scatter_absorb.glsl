#iChannel0 "self"

#iChannel1 "../cubemaps/shadertoy_uffizi_gallery/{}.jpg"
#iChannel1::Type "CubeMap"


#define PI 3.1415926
#define ZERO min(iTime, 0.)


uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }



vec3 light(vec3 rd) {
    vec3 col = texture(iChannel1, rd).xyz;
    //col = vec3(0.5, 0.5, 0.6);
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

vec3 sampleUniformSphere() {
    float u = 2.0*PI*rand01();
    float v = 2.0*rand01()-1.0;
    return vec3(vec2(cos(u), sin(u))*sqrt(1.0-v*v), v);
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



const int mat_none = -1;
const int mat_background = 0;
const int mat_lambertian = 1;
const int mat_refractive = 2;


void calcAbsorb(in vec3 p, out vec3 tabs, out vec3 sabs, out float k) {
    p -= vec3(0, 0, 1);  // center of sphere
    tabs = vec3(0.9, 0.3, 0.2);
    float c = 0.5+0.3*sin(2.0*p.x)*cos(3.0*p.y)*sin(2.0*p.z)+dot(p, vec3(0.1));
    c = clamp(c, 0., 1.);
    sabs = pow(vec3(0.8*pow(c, 0.5), 0.4*pow(c, 0.6), c), vec3(0.5));
    k = 0.1+10.0*pow(c, 5.0);
    k = max(k, 1e-3);
}

void calcScatter(in vec3 ro, inout vec3 rd,
        inout float mt, inout vec3 m_col, out int material, inout vec3 min_n) {
    float p = 1.0, hit_p = rand01();
    float dt = 0.02;
    for (float t = 1e-4; t < mt-dt; t += dt) {
        vec3 tabs, sabs; float k;
        calcAbsorb(ro + rd * (t+0.5*dt), tabs, sabs, k);
        float dp = exp(-k*dt);
        if (p * dp < hit_p) {
            dt *= log(p/hit_p)/k;
            mt = t + dt;
            rd = sampleUniformSphere();
            m_col *= sabs * exp(-tabs*dt);
            material = mat_none;
            min_n = vec3(0.0);
            return;
        }
        p *= dp;
        m_col *= exp(-tabs*dt);
    }
    material = mat_refractive;
}


vec3 mainRender(vec3 ro, vec3 rd) {

    vec3 m_col = vec3(1.0), col;
    bool is_inside = false;

    for (int iter = int(ZERO); iter < 128; iter++) {
        ro += 1e-4f*rd;
        if (is_inside != (length(ro-vec3(0,0,1))<1.0)) return vec3(1, 0, 0);

        vec3 n, min_n;
        float t, min_t = 1e12;
        vec3 min_ro = ro, min_rd = rd;
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
                calcScatter(ro, min_rd, min_t, col, material, min_n);
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
            //if (iter == 0) return vec3(0.f);
            col = light(rd);
            return m_col * col;
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
        m_col *= col;
        if (dot(rd, min_n) < 0.0) {
            is_inside = !is_inside;
        }
    }
    return m_col;
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
