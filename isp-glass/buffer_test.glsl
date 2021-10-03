// test buffer with vscode Shader-toy plugin

#iChannel0 "self"

#define PI 3.14159265358979

// set SAMPLE to a smaller value if you have a slow machine
#define SAMPLE 1
#define EPSILON 1e-4
#define MAX_STEP 64
#define MAX_DIST 10.0
#define MAX_RECU 50

vec2 CENTER = vec2(0.0, 0.0);
float SCALE = 100.0;
float BULB = 2.0;
float INDEX = 1.5;


float hash(vec2 x){
    return fract(sin(dot(x,vec2(12.9898,78.233)))*43758.5453);
}

vec2 Refract(vec2 I, vec2 N, float n1, float n2, out float R) {
    float eta = n1 / n2;
    float ci = -dot(N, I);
    if (ci < 0.) ci = -ci, N = -N;
    float ct = sqrt(1.0 - eta * eta * (1.0 - ci * ci));
    vec2 r = I * eta + N * (eta * ci - ct);
    float Rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct);
    float Rp = (n1 * ct - n2 * ci) / (n1 * ct + n2 * ci);
    R = 0.5 * (Rs * Rs + Rp * Rp);
    return r;
}


float sdBulb(vec2 p) {
    //return length(p - vec2(3.0)) - 1.0;
    //return length(vec2(abs(p.x) - 3.0, abs(p.y) - 3.0)) - 0.5;
    return length(vec2(p.x, p.y - 3.0)) - 1.0;
}

float sdObj(vec2 p) {
    //return length(p)-1.0;    // circle
    //return (abs(p.x)>0.8?length(vec2(abs(p.x)-0.8,p.y)):abs(p.y))-0.8;    // capsule
    //return max(abs(p.x) - 1.2, abs(p.y) - 0.75);    // rectangle
    //return max(length(vec2(p.x,p.y-0.4))-1.0,p.y-0.5);    // semi-circle
    return min(max(abs(p.x),abs(p.y)-1.2), max(abs(p.x)-0.7,abs(p.y-0.5)))-0.2;    // cross
    //p=abs(p); return min((p.y>1.0?length(p-vec2(0.65,1.0)):abs(p.x-0.65))-0.2, max(p.x-0.65,p.y)-0.2);    // letter H
}

vec2 gradient(vec2 p) {
    float k = 0.001;
    float u = sdObj(vec2(p.x + k, p.y)) - sdObj(vec2(p.x - k, p.y));
    float v = sdObj(vec2(p.x, p.y + k)) - sdObj(vec2(p.x, p.y - k));
    return vec2(u, v) * (0.5 / k);
}

float traceRay(vec2 p, vec2 d) {
    int N = 0;
    while (N++ < MAX_RECU) {
        float t = 10.0*EPSILON, dt, sdb, sdo, ot, it, R;
        vec2 q, n, r;
        int i; for (i = 0; i < MAX_STEP; i++) {
            q = p + d * t;
            sdb = sdBulb(q);
            if (sdb <= EPSILON) return BULB;
            sdo = sdObj(q);
            dt = sdb > sdo ? sdo : sdb;
            if (abs(dt) <= EPSILON) {
                n = normalize(gradient(q)), r;
                if (dt >= 0.0) r = Refract(d, n, 1.0, INDEX, R);
                else r = Refract(d, n, INDEX, 1.0, R);
                if (isnan(R)) R = 1.0;  // bug fixed: 0.0*R!=0.0 got optimized
                break;
            }
            t += abs(dt);
            if (t > MAX_DIST) return 0.0;
        }
        if (i == MAX_STEP) return 0.0;

        t = hash(q + d * R + iTime);
        if (t < R) p = q, d = reflect(d, n);
        else p = q, d = r;
    }
    return 0.0;
}

float Sample(vec2 p) {
    float c = 0.0;
    float s = 1.0 / SCALE, h = -0.5 / SCALE;
    for (int i = 0; i < SAMPLE + min(iFrame, 0); i++) {
        float a = 2.0 * PI * (float(i) + hash(p + vec2(i) + iTime)) / float(SAMPLE);
        vec2 d = vec2(cos(a), sin(a));
        c += traceRay(p + vec2(hash(p + iTime - float(i))) * s, d);
    }
    return c / float(SAMPLE);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    SCALE = 0.2*sqrt(iResolution.x*iResolution.y);
    vec2 p = (fragCoord - iResolution.xy * 0.5) * (1.0 / SCALE) + CENTER * 0.5;
    float c = Sample(p);
    vec3 col = vec3(c);
    vec4 rgbn = texelFetch(iChannel0, ivec2(fragCoord), 0);
    if (iMouse.z>0.) rgbn.w = 0.0;
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}
