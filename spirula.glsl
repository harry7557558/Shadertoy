// https://www.google.com/search?q=spirula&tbm=isch
// Not exactly the same...

#define PI 3.1415926

float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) - k * h * (1.0 - h);
}
float smax(float a, float b, float k) {
    return -smin(-a, -b, k);
}

float map(in vec3 p) {
    const float k = 0.15;  // r = e^kθ
    // https://swiftcoder.wordpress.com/2010/06/21/logarithmic-spiral-distance-field/
    float r = length(p.xy);  // Cartesian to polar
    float a = atan(-p.y, p.x);  // clockwise
    float n = (log(r)/k-a)/(2.*PI);  // approximate index of revolution
    n = min(n, 1.0);  // θ < 1.0*2π
    float r0 = exp(k*(2.*PI*floor(n)+a));  // distance to inner curve
    float r1 = exp(k*(2.*PI*ceil(n)+a));  // distance to outer curve
    float d = min(abs(r1-r), abs(r0-r));  // choose the minimum distance
    vec2 op = exp(k*(3.*PI))*vec2(-1.,0.);  // opening of the shell
    d = min(d,length(p.xy-op));  // fix discontinuity
    d = length(vec2(d, p.z+0.0*r));  // 2d to 3d
    float w = cos(20.0*a);  // "wrinkles" of the shell
    r = 0.4*pow(r,0.7);  // offset radius
    r *= 1.0+0.02*(1.0-w);  // apply "wrinkles"
    d = d - max(r, 0.0);  // offset
    float h = d + 0.5*r*(w+1.0);  // "holes" inside
    h = smin(h, length(p-vec3(op-vec2(0.15,-0.5),0.))-0.4*exp(2.0*PI*k), 0.2);  // opening
    d = smax(d, -0.1-h, 0.1);  // holes
    //d = max(d, p.z);  // cut it through to see its inside
    return 0.7*d;  // distance for raymarching
}

vec3 mapNormal(vec3 p) {
    mat3 k = mat3(p,p,p) - mat3(0.001);
    return normalize(map(p) - vec3(map(k[0]),map(k[1]),map(k[2])));
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // set camera

    float rx = iMouse.z>0.?3.14*(iMouse.y/iResolution.y)-1.57:0.5;
    float rz = iMouse.z>0.?-iMouse.x/iResolution.x*4.0*3.14:-4.0;

    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);

    vec3 ro = 20.0*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 2.0*length(iResolution.xy));
    rd = normalize(rd);

    // raymarching
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

    // shading
    vec3 light = normalize(vec3(0.5, 0.5, 1.0));
    vec3 n = mapNormal(ro+rd*t);
    float fresnel = 0.2+1.4*pow(1.0+dot(rd,n),2.0);
    float specular = pow(max(dot(rd-2.0*dot(rd,n)*n,light),0.0),100.);
    float diffuse = max(dot(n,light),0.0);
    //vec3 col = vec3(0.2)*fresnel + vec3(0.15)*specular + vec3(0.8,0.85,0.6)*diffuse;
    vec3 col = vec3(0.2)+vec3(0.2,0.15,0.1)*n.y+vec3(0.6,0.65,0.7)*diffuse;
    fragColor = vec4(col, 1.0);
}
