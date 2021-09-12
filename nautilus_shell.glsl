#define PI 3.1415926

float map( in vec3 p0, out vec3 col ) {
    p0 -= vec3(0.7, 0, 0);
    vec3 p = p0;

    const float b = 0.17;

    float r = length(p.xy);
    float a = mix(0.0, 0.5, smoothstep(0.0, 1.0, 0.5*(r-0.6)));
    p.xy = mat2(cos(a),-sin(a),sin(a),cos(a))*p.xy;
    float t = atan(p.y, p.x);
    
    // shell opening
    float ro = exp(b*PI);
    float d = length(vec2(length(p.xz-vec2(-ro,0))-ro,p.y));
    float u = t, dx = r-ro, dy = p.z;

    // spiral
    // r(n) = exp(b*(2.*PI*n+t)), (x-r)^2+y^2=r^2, solve for n
    float n = (log((r*r+p.z*p.z)/(2.*r))/b-t)/(2.0*PI);
    n = min(n, 0.0);
    float n0 = floor(n), n1 = ceil(n);
    float r0 = exp(b*(2.*PI*n0+t)), r1 = exp(b*(2.*PI*n1+t));
    float d0 = abs(length(vec2(r-r0,p.z))-r0);
    float d1 = abs(length(vec2(r-r1,p.z))-r1);
    if (d0 < d) d = d0, u = 2.*PI*n0+t, dx = r-r0, dy = p.z;
    if (d1 < d) d = d1, u = 2.*PI*n1+t, dx = r-r1, dy = p.z;

    // cells, possibly cause discontinuities
    const float f = 2.4;
    float s0 = t + 2.0*PI*(n0+0.5);
    float v = fract(n);
    float s = f*s0 + 1.0*pow(0.25-(v-0.5)*(v-0.5), 0.5)+0.5*v;
    s += 0.5/(40.0*length(vec2(v-0.5,p.z))+1.0);
    float sf = fract(s);
    sf = s0>-1.8 ? abs(s+3.25) : min(sf, 1.0-sf);
    float w = sf/f*exp(b*(s0+PI));
    if (length(p*vec3(1,1,1.5))<3.0) d = min(d, 0.5*w+0.012);

    // texture
    d += 0.002*r*sin(40.*u);
    //d += 0.002*r*sin(40.*atan(dy,dx));

    d = abs(d)-0.8*max(0.02*pow(r,0.4),0.02);
    //d = max(d, p0.x);
    //d = max(d, p0.z);
    d = max(d, abs(p0.z)-0.1);

    col = vec3(0.8);

    return 0.8*d;
}


vec3 mapNormal(vec3 p) {
    mat3 k = mat3(p,p,p) - mat3(0.001);
    vec3 c;
    return normalize(map(p,c) - vec3(map(k[0],c),map(k[1],c),map(k[2],c)));
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    float rx = 3.14*(iMouse.y/iResolution.y)-1.57;
    float rz = -iMouse.x/iResolution.x*4.0*3.14;

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
    vec3 col;
    for (int i=min(int(iTime),0); i<100; i++) {
        float dt = map(ro+rd*t, col);
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
    col = 0.2+0.1*n.y+col*max(dot(n, light), 0.0);
    fragColor = vec4(vec3(col), 1.0);
    return;
}
