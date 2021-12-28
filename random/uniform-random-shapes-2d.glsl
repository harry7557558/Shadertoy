#iChannel0 "self"

// Showcase generating uniformly random points inside 2D shapes
// Each random function takes exactly two random values between 0 and 1, no rejection sampling

// Functions are derived based on the following principles:
// For circles and other curves, solve for a constant Jacobian (equals the reciprocal of shape area)
// For polygons/stars, divide into multiple triangles
// For results of boolean operations of circles, parameterize and apply inverse CDF sampling


#define PI 3.1415926


vec3 intersectRay(vec3 ro, vec3 rd) {
    // inspired by https://www.shadertoy.com/view/MdfBRX by BigWIngs
    float t = -(ro.z+2.0)/rd.z;
    vec3 p = ro+rd*t;
    vec3 col = rd.z > 0.0 ?
         mix(vec3(0.0), 0.3*vec3(0.2,0.5,0.8), pow(rd.z,1.0)) :  // sky
         abs(p.x)<12.0 ? 0.3*vec3(0.05,0.06,0.06) : 0.3*vec3(0.04,0.05,0.02);  // road
    // all distance are in meters
    for (float d=10.0; d<=120.0; d+=10.0) {  // street light
        t = -(ro.y+d)/rd.y;
        p = ro+rd*t;
        p.x = abs(p.x);
        if (length(p.xz-vec2(8.0,5.0))<0.15)
            col += 20.0*vec3(1.0,0.8,0.5);
    }
    for (float d=4.0; d<=54.0; d+=10.0) {  // head light
        t = -(ro.y+d)/rd.y;
        p = ro+rd*t;
        p.x = abs(p.x+4.0);
        if (length(p.xz-vec2(2.0,-0.5))<0.1)
            col += 15.0*vec3(0.8,0.8,1.0);
    }
    for (float d=6.0; d<=56.0; d+=10.0) {  // tail light
        t = -(ro.y+d)/rd.y;
        p = ro+rd*t;
        p.x = abs(p.x-4.0);
        if (length(p.xz-vec2(2.0,-0.5))<0.08)
            col += 10.0*vec3(1.0,0.1,0.0);
    }
    return col;
}


// circle with radius 1
vec2 randCircle(float rand1, float rand2) {
    float u = 2.0*PI*rand1;  // θ
    float v = sqrt(rand2);  // r
    return v*vec2(cos(u), sin(u));
}

// sector with abs(θ)<angle
vec2 randSector(float angle, float rand1, float rand2) {
    float u = angle*(2.0*rand1-1.0);  // θ
    float v = sqrt(rand2);  // r
    return v*vec2(cos(u), sin(u));
}

// x²-|x|y+y² < 1
vec2 randHeart(float rand1, float rand2) {
    float u = 2.0*PI*rand1;  // θ
    float v = sqrt(rand2);  // r
    vec2 c = v*vec2(cos(u), sin(u));  // unit circle
    c = mat2(1.0,1.0,-0.577,0.577)*c;  // ellipse
    if (c.x<0.0) c.y=-c.y;  // mirror
    return c;
}

// ellipse with major and minor radius
vec2 randEllipse(float rx, float ry, float rand1, float rand2) {
    float u = 2.0*PI*rand1;  // θ
    float v = sqrt(rand2);  // r
    vec2 circ = v * vec2(cos(u), sin(u));  // unit circle
    return vec2(rx, ry) * circ;  // linear transform
}

// ring formed by two concentric circles with radius r0 and r1
vec2 randConcentric(float r0, float r1, float rand1, float rand2) {
    float u = 2.0*PI*rand1;  // θ
    float v = sqrt(mix(r0*r0, r1*r1, rand2));  // r
    return v * vec2(cos(u), sin(u));  // polar to Cartesian
}

// intersection of two circles with centers (0,±c) and radius r
vec2 randIntersection(float c, float r, float rand1, float rand2) {
    // https://www.desmos.com/calculator/sctxdxh1td
    float u = rand1;
    float v = 2.0*rand2-1.0;
    float x1 = sqrt(r*r-c*c);
    float i1 = 0.5*r*r*asin(x1/r)+0.5*x1*sqrt(r*r-x1*x1)-c*x1;
    u = 2.0*i1*u - i1;
    float x = 0.0;
    for (int iter=0; iter<6; iter++) {  // Newton-Raphson
        float cdf = 0.5*(r*r*asin(x/r)+x*sqrt(r*r-x*x))-c*x;
        float pdf = sqrt(r*r-x*x)-c;
        x -= (cdf-u)/pdf;
    }
    float y = (sqrt(r*r-x*x)-c) * v;
    return vec2(x, y);
}

// union of two circles with centers (±c,0) and radius r
vec2 randUnion(float c, float r, float rand1, float rand2) {
    // https://www.desmos.com/calculator/ddyum0rkgw
    float u = 2.0*rand1-1.0;
    float v = 2.0*rand2-1.0;
    float s = sign(u);
    float x1 = r+c;
    float h0 = sqrt(r*r-c*c);
    float i0 = -0.5*(c*h0+r*r*atan(c/h0));
    float i1 = 0.25*PI*r*r;
    u = mix(i0, i1, abs(u));
    float x = c;  // start at a point of inflection
    for (int iter=0; iter<6; iter++) {  // Newton-Raphson
        float pdf = sqrt(r*r-(x-c)*(x-c));
        float cdf = 0.5*((x-c)*pdf+r*r*atan((x-c)/pdf));
        x -= (cdf-u)/pdf;
    }
    float y = sqrt(r*r-(x-c)*(x-c)) * v;
    return vec2(s*x, y);
}

// subtraction of a circle with center (c,0) from a circle with (0,0)
vec2 randSubtraction(float c, float r, float rand1, float rand2) {
    // https://www.desmos.com/calculator/jk1okdxoks
    float u = rand1;
    float v = 2.0*rand2-1.0;
    float x1 = 0.5*c;
    float y1 = sqrt(r*r-0.25*c*c);
    float a1 = y1*c;
    float i1 = 0.5*(x1*sqrt(r*r-x1*x1)+r*r*asin(x1/r))-x1*y1;
    float a2 = 2.0*i1;
    float th = a1/(a1+a2);
    if (u < th) {
        float y = y1*v;
        float x = c*(u/th) - sqrt(r*r-y*y);
        return vec2(x, y);
    }
    u = i1 * (2.*((u-th)/(1.-th)) - 1.);
    float x = 0.0;
    for (int iter=0; iter<6; iter++) {  // Newton-Raphson
        float cdf = 0.5*(x*sqrt(r*r-x*x)+r*r*asin(x/r))-y1*x;
        float pdf = sqrt(r*r-x*x)-y1;
        x -= (cdf-u)/pdf;
    }
    float y = sign(v) * (y1+(sqrt(r*r-x*x)-y1)*abs(v));
    return vec2(x, y);
}

// -1 <= x,y < 1
vec2 randSquare(float rand1, float rand2) {
    float u = 2.0*rand1-1.0;
    float v = 2.0*rand2-1.0;
    return vec2(u, v);
}

// convex quadrilateral defined by ccw vertices
vec2 randQuad(vec2 p0, vec2 p1, vec2 p2, vec2 p3, float rand1, float rand2) {
    float u = rand1;
    float v = rand2;
    p1 -= p0, p2 -= p0, p3 -= p0;  // one point and three vectors
    float area1 = 0.5*(p1.x*p2.y-p1.y*p2.x);  // area of the first triangle
    float area2 = 0.5*(p2.x*p3.y-p2.y*p3.x);  // area of the second triangle
    float th = area1 / (area1+area2);  // threshold to decide which triangle
    if (u > th) p1 = p2, p2 = p3, u = (u-th)/(1.0-th);  // use the second triangle
    else u /= th;  // use the first triangle
    return p0 + mix(p1, p2, u) * sqrt(v);  // sample inside triangle
    if (u+v>1.) u=1.-u, v=1.-v; return p0 + p1*u + p2*v;  // avoid square root, may or may not be faster
}

// regular n-gon with radius 1
vec2 randPolygon(float n, float rand1, float rand2) {
    float u = n*rand1;
    float v = rand2;
    float ui = floor(u);  // index of triangle
    float uf = fract(u);  // interpolating in triangle
    vec2 v0 = vec2(cos(2.*PI*ui/n), sin(2.*PI*ui/n));  // triangle edge #1
    vec2 v1 = vec2(cos(2.*PI*(ui+1.)/n), sin(2.*PI*(ui+1.)/n));  // triangle edge #2
    return sqrt(v) * mix(v0, v1, uf);  // sample inside triangle
}

// regular n-star with normalized size
vec2 randStar(float n, float rand1, float rand2) {
    float u = n*rand1;
    float v = rand2;
    float ui = floor(u);  // index of triangle
    float uf = fract(u);  // interpolating in rhombus
    vec2 v0 = vec2(cos(2.*PI*ui/n), sin(2.*PI*ui/n));  // rhombus edge #1
    vec2 v1 = vec2(cos(2.*PI*(ui+1.)/n), sin(2.*PI*(ui+1.)/n));  // rhombus edge #2
    vec2 p = v0 * v + v1 * uf;  // sample rhombus
    return p / (n*sin(2.*PI/n)/PI);  // normalize size
}


// random number generator
float vanDerCorput(float n, float b) {
    float x = 0.0;
    float e = 1.0 / b;
    while (n > 0.0) {
        float d = mod(n, b);
        x += d * e;
        e /= b;
        n = floor(n / b);
    }
    return x;
}

// main
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // random number seed
    vec3 p3 = fract(fragCoord/iResolution.xy*.1031).xyx;
    p3 += dot(p3, p3.zyx + 31.32);
    float h = fract((p3.x + p3.y) * p3.z);
    float seed = floor(65536.*h) + float(iFrame);

    // camera parameters
    const vec3 POS = vec3(3.0,3.0,0.0);
    const float SCALE = 1.5;  // larger = smaller (more view field)
    const float DIST = 0.2;  // larger = smaller
    const float VIEW_FIELD = 0.5;  // larger = larger + more perspective
    const float APERTURE = 0.01;  // larger = blurred

    // sample aperture shape
    float rand1 = vanDerCorput(seed, 2.);
    float rand2 = vanDerCorput(seed, 3.);
    vec2 rnd;
    //rnd = randCircle(rand1, rand2);
    //rnd = randSector(0.8*PI, rand1, rand2);
    rnd = randHeart(rand1, rand2);
    rnd = randEllipse(1.2, 0.8, rand1, rand2);
    //rnd = randConcentric(0.6, 1.1, rand1, rand2);
    rnd = randIntersection(0.9, 1.6, rand1, rand2);
    rnd = randUnion(0.6, 0.8, rand1, rand2);
    rnd = randSubtraction(1.0, 1.0, rand1, rand2);
    //rnd = randSquare(rand1, rand2);
    //rnd = randQuad(vec2(-1.0,-1.0), vec2(0.8,-0.9), vec2(1.0,1.2), vec2(-0.4,0.8), rand1, rand2);
    //rnd = randPolygon(5., rand1, rand2);
    //rnd = randStar(5., rand1, rand2);

    // camera
    vec3 ro = POS+vec3(0,DIST,0);
    vec2 randuv = vec2(vanDerCorput(seed,5.), vanDerCorput(seed, 7.));
    vec2 uv = SCALE*(2.0*(fragCoord.xy+randuv-0.5)/iResolution.xy-1.0);
    vec2 sc = iResolution.xy/length(iResolution.xy);
    vec2 offset = APERTURE*rnd;
    ro.xz += offset;
    vec3 rd = vec3(VIEW_FIELD*uv*sc+vec2(-0.2,0.0)-offset/DIST, -1.0).xzy;

    // calculate pixel color
    vec3 col = intersectRay(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(int(fragCoord.x), int(fragCoord.y)), 0);
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}
