#iChannel0 "self"

#iChannel1 "../cubemaps/shadertoy_uffizi_gallery/{}.jpg"
#iChannel1::Type "CubeMap"



// Left: without importance sampling
// Right: with importance sampling for GGX normal
// The two versions look slightly different. Should be the same.

// if there is no bug with GGX importance sampling,
// left and right should be the same when this is set to 1
#define DEBUG_GGX 0

// References:
// http://www.codinglabs.net/article_physically_based_rendering_cook_torrance.aspx
// https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models
// https://agraphicsguy.wordpress.com/2015/11/01/sampling-microfacet-brdf/
// https://www.shadertoy.com/view/3slSzn
// https://computergraphics.stackexchange.com/questions/4394/path-tracing-the-cook-torrance-brdf

// A Desmos implementation of the BRDF:
// https://www.desmos.com/calculator/zccocbuygk


#define PI 3.1415926
#define ZERO min(iTime,0.)


// random number generator
uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }


// cubemap lighting
vec3 light(vec3 rd) {
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


// sample hemisphere distributions
vec3 sampleCosWeighted(vec3 n) {
    vec3 u = normalize(cross(n, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, n);
    float rn = rand01();
    float an = 2.0*PI*rand01();
    vec2 rh = sqrt(rn) * vec2(cos(an), sin(an));
    float rz = sqrt(1. - rn);
    return rh.x * u + rh.y * v + rz * n;
}
vec3 sampleUniformHemisphere(vec3 n) {
    vec3 u = normalize(cross(n, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, n);
    float rz = rand01();
    float an = 2.0*PI*rand01();
    vec2 rh = sqrt(1.0-rz*rz) * vec2(cos(an), sin(an));
    return rh.x * u + rh.y * v + rz * n;
}


// sample GGX, return weight (Fr divided by PDF)
// At least one of the two functions have a bug
float sampleGgxSimple(in vec3 wi, in float alpha, out vec3 wo) {
    wo = sampleUniformHemisphere(vec3(0, 0, 1));
    vec3 m = normalize(wi+wo);
    float denom = (alpha*alpha-1.)*m.z*m.z+1.;
    float Fr = alpha*alpha / (PI * denom*denom);
    return Fr / (1.0/(2.0*PI));
}
float sampleGgxImportance(in vec3 wi, in float alpha, out vec3 wo) {
    float su = 2.0*PI*rand01();
    float sv = rand01();
    //sv = acos(sqrt((1.0-sv)/((alpha*alpha-1.)*sv+1.)));
    sv = atan(alpha*sqrt(sv/(1.0-sv)));
    vec3 h = vec3(sin(sv)*vec2(cos(su),sin(su)), cos(sv));
    wo = -(wi-2.0*dot(wi,h)*h);
    return wo.z<0. ? 0. : 4.0*dot(wi, h);
}


// as global variable
bool ImportanceSampling;


// return random output ray direction, multiply m_col by weight
vec3 sampleCookTorrance(
    vec3 wi, vec3 n,
    float alpha,  // roughness
    float f0,  // ratio of reflection along the normal
    float lambertian,  // ratio of lambertian coefficient
    vec3 lambert_col,  // lambertian color
    vec3 microfacet_col,  // microfacet color
    inout vec3 m_col
    ) {

#if !DEBUG_GGX
    if (ImportanceSampling) {
        if (rand01() < lambertian) {
            vec3 wo = sampleCosWeighted(n);
            m_col *= lambert_col;
            return wo;
        }
    }
#endif

    vec3 u = normalize(cross(n, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, n);
    wi = vec3(dot(wi, u), dot(wi, v), dot(wi, n));
    vec3 wo, m;  // out and half vector

    // GGX divided by PDF
    float D;
    if (ImportanceSampling) D = sampleGgxImportance(wi, alpha, wo);
    else D = sampleGgxSimple(wi, alpha, wo);
    // uncomment to see if lambert sampling is ok
    //D = sampleGgxSimple(wi, alpha, wo);

    // debug GGX importance sampling
#if DEBUG_GGX
    m_col *= 0.25*D*wo.z; return wo.x*u+wo.y*v+wo.z*n;
#endif
    m = normalize(wi+wo);

    // Geometry
    float tan2_theta_i = (1.0-wi.z*wi.z)/(wi.z*wi.z);
    float tan2_theta_o = (1.0-wo.z*wo.z)/(wo.z*wo.z);
    float lambda_i = 0.5*(sqrt(1.0+alpha*alpha*tan2_theta_i)-1.0);
    float lambda_o = 0.5*(sqrt(1.0+alpha*alpha*tan2_theta_o)-1.0);
    float G = 1.0/(1.0+lambda_i+lambda_o);

    // Fresnel
    float F = f0 + (1.0-f0)*pow(1.0-dot(wi, m), 5.0);

    // Put all together
    float Fr = D*G*F / (4.0*wi.z*wo.z+1e-4);
    float Fr_cos = Fr * wo.z;  // wo is the direction of light in path tracing
    if (ImportanceSampling) {
        m_col *= Fr_cos * microfacet_col;
    }
    else {
        vec3 col = lambertian * lambert_col / (1.0/2.0) + (1.0-lambertian) * microfacet_col * Fr;
        m_col *= col * wo.z;
    }
    return wo.x * u + wo.y * v + wo.z * n;
}


// path tracing
vec3 mainRender(vec3 ro, vec3 rd) {

    const int background = -1;
    const int lambertian = 0;
    const int mat_floor = 1;

    vec3 m_col = vec3(1.0), col;
    bool is_inside = false;

    for (int step = int(ZERO); step < 64; step++) {
        ro += 1e-4f*rd;
        vec3 n, min_n;
        float t, min_t = 1e12;
        int material = background;

        // plane
        t = -ro.z / rd.z;
        if (t > 0.0) {
            min_t = t, min_n = vec3(0, 0, 1);
            col = vec3(0.9, 0.95, 0.98);
            //material = lambertian;
            material = mat_floor;
        }

        // objects
        for (float i = 0.; i < 6.; i++) {
            t = min_t;
            vec3 pos = vec3(2.2*vec2(cos(2.*PI*i/6.), sin(2.*PI*i/6.)), 1.0+1e-4);
            if (intersectSphere(pos, 1.0, ro, rd, t, n)) {
                min_t = t, min_n = n;
                col = vec3(1.0);
                material = int(i)+2;
            }
        }

        // update ray
        if (material == background) {
            col = light(rd);
            return m_col * col;
        }
        ro = ro + rd * min_t;
        min_n = dot(rd, min_n) < 0. ? min_n : -min_n;
        if (material == lambertian) rd = sampleCosWeighted(min_n), m_col *= col;
        if (material == mat_floor) rd = sampleCookTorrance(-rd, min_n, 0.3, 0.4, 0.2, col, col, m_col);
        if (material == 2) rd = sampleCookTorrance(-rd, min_n, 0.8, 0.4, 0.0, vec3(1.0), col, m_col);
        if (material == 3) rd = sampleCookTorrance(-rd, min_n, 0.4, 0.4, 0.0, vec3(1.0), col, m_col);
        if (material == 4) rd = sampleCookTorrance(-rd, min_n, 0.1, 0.4, 0.0, vec3(1.0), col, m_col);
        if (material == 5) rd = sampleCookTorrance(-rd, min_n, .02, 0.4, 0.0, vec3(1.0), col, m_col);
        if (material == 6) rd = sampleCookTorrance(-rd, min_n, 0.1, 0.4, 0.5, vec3(0.8,0.4,0.1), vec3(1.0), m_col);
        if (material == 7) rd = sampleCookTorrance(-rd, min_n, .01, 0.4, 0.2, vec3(0.8,1.0,0.8), vec3(0.6,0.6,1.0), m_col);
        if (m_col == vec3(0.0)) break;
        if (dot(rd, min_n) < 0.0) is_inside = !is_inside;
        if (is_inside) return vec3(1, 0, 0);  // should not happen
    }
    return m_col;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    ImportanceSampling = fragCoord.x > 0.5*iResolution.x;

    // https://www.shadertoy.com/view/4djSRW by David Hoskins, MIT licence
    vec3 p3 = fract(vec3(fragCoord, iFrame) * .1031);
    p3 += dot(p3, p3.zyx + 31.32);
    float h = fract((p3.x + p3.y) * p3.z);
    seed = uint(16777216.*h);

    // camera
    vec2 mouse_uv = iMouse.z==0. ? vec2(0.6,0.4) : iMouse.xy/iResolution.xy;
    float rx = 1.8*mouse_uv.y-0.3;
    float rz = -mouse_uv.x*4.0*3.14;
    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);
    vec3 ro = 10.0*w + vec3(0, 0, 0.7);
    vec2 uv = 2.0*(fragCoord.xy+vec2(rand01(),rand01())-0.5)/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 1.8*length(iResolution.xy));
    rd = normalize(rd);

    // calculate pixel color
    vec3 col = mainRender(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(fragCoord), 0);
    if (iMouse.z>0.) rgbn.w = 0.0;
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}
