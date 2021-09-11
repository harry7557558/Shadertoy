#define PI 3.1415926

float map( in vec3 p, out vec3 col ) {

    const float b = 0.17;
    
    float r = length(p.xy);
    float a = mix(0.0, 0.8, smoothstep(0.0, 1.0, 0.5*(r-0.6)));
    //p.xy = mat2(cos(a),-sin(a),sin(a),cos(a))*p.xy;
    float t = atan(p.y, p.x);
 
    float n = (log(r)/b-t)/(2.0*PI);
    n = min(n, 0.0);
    
    float n0 = floor(n), n1 = ceil(n);
    
    float x0 = exp(b * (t + 2.0*PI*n0));
    float x1 = exp(b * (t + 2.0*PI*n1));
    float r0 = 1.0*x0;
    float r1 = 1.0*x1;

    float h0 = p.z + 0.0*(x0-1.0);
    float h1 = p.z + 0.0*(x1-1.0);
    float d0 = length(vec2(x0-r,h0)) - r0;
    float d1 = length(vec2(x1-r,h1)) - r1;

    float d, dx, dy;
    if (d0 < 0.0) d = d0, dx = x0-r, dy = h0, col=vec3(0,0,1);
    else if (d1 < 0.0 && d1<-d0) d = -d0, dx = x0-r, dy = h0, col=vec3(0,1,1);
    else if (d1 < 0.0) d = d1, dx = x1-r, dy = h1, col=vec3(1,0,1);
    else if (d0 < d1) d = d0, dx = x0-r, dy = h0, col=vec3(1,0,0);
    else d = d1, dx = x1-r, dy = h1, col=vec3(1,1,0);

    //col = vec3(dx,0,0);
    //col = mix(col, vec3(0.5+0.5*sin(40.0*atan(dy,dx))), 0.5);

    //d += 0.002*r*sin(40.*t);
    d += 0.002*r*sin(40.*atan(dy,dx));

    d = abs(d)-0.1*r;
    //d = max(d, p.z);
    d = max(d, p.x);

    return 0.5*d;
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
