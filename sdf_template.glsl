float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) - k * h * (1.0 - h);
}
float smax(float a, float b, float k) {
    return -smin(-a, -b, k);
}

float map(in vec3 p) {
    float r = 1.0 + 0.001*sin(64.0*atan(length(p.xy), p.z));
    float d1 = length(p) - r;
    float b = 0.7 + 0.05*sin(10.0*p.x)*sin(10.0*p.y)*sin(10.0*p.z);
    float d2 = max(max(abs(p.x + 2.0*sin(2.0*iTime)), abs(p.y)), abs(p.z)) - b;
    return smin(d1, d2, 0.5);
}

vec3 mapNormal(vec3 p) {
    mat3 k = mat3(p,p,p) - mat3(0.001);
    return normalize(map(p) - vec3(map(k[0]),map(k[1]),map(k[2])));
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    float rx = iMouse.z!=0.?3.14*(iMouse.y/iResolution.y)-1.57:0.3;
    float rz = iMouse.z!=0.?-iMouse.x/iResolution.x*4.0*3.14:0.5*iTime-2.0;

    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);

    vec3 ro = 10.0*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 2.0*length(iResolution.xy));
    rd = normalize(rd);

    float t0 = 0.01;
    float t1 = 2.0*length(ro);
    float t = t0;
    for (int i=min(int(iTime),0); i<100; i++) {
        float dt = map(ro+rd*t);
        t += dt;
        if (abs(dt) < 1e-3) break;
        if (t > t1) {
            fragColor = vec4(0, 0, 0, 1);
            return;
        }
    }
    vec3 light = normalize(vec3(0.5, 0.5, 1.0));
    //vec3 light = normalize(w+u+v);
    vec3 n = mapNormal(ro+rd*t);
    float col = 0.2+0.1*n.y+0.7*max(dot(n, light), 0.0);
    fragColor = vec4(vec3(col), 1.0);
    return;
}
